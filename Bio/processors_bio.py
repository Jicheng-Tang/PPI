import multiprocessing as mp
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import glob
import os
import requests
import re


#1、过滤去重（基于验证实验类型）
def filter_bio():
    """
    基于验证实验类型(Experimental System Type)过滤并去重 BioGRID PPI 数据。
    该函数递归处理制表符分隔的文本文件，过滤高通量、低置信实验类型，
    使用指定的低通量实验系统，提取交互物 ID 和名称，
    规范化配对以避免顺序重复，并将结果保存到 'high_confidence_Bio_ppi.csv'。
    """
    # 过滤原则 (基于DEBUG, 用Experimental System具体方法)
    available_systems = [
        'Affinity Capture-Western', 'Affinity Capture - Western',  # co-IP
        'Two hybrid', 'Two-hybrid', 'Two hybrid array', 'Two-hybrid array',  # 双杂交
        'FRET',  # 荧光
        'Biochemical Fractionation', 'Biochemical Activity', 'Biochemical Fractionation with biochemical activity'  # 纯化/生化
    ]
    exclude_systems = [
        'Affinity Capture-MS', 'Affinity Capture - MS',  # AP-MS
        'Affinity Capture-Luminescence'  # 高通量
    ]
    physical_type = 'physical'  # Type广义过滤
    low_throughput = 'Low Throughput'

    # 步骤1: 递归搜txt
    txt_files = glob.glob('**/*.txt', recursive=True)
    print(f"Found {len(txt_files)} txt files:")
    for f in txt_files[:5]:
        print(f"  - {f}")
    if len(txt_files) > 5:
        print("  ...")

    all_data = []
    for file_path in txt_files:
        print(f"\nProcessing {os.path.basename(file_path)}...")
        try:
            df_temp = pd.read_csv(file_path, sep='\t', low_memory=False)
        
            # 过滤: Type=='physical' + System in available + 排除 + Low Throughput
            df_filtered = df_temp[
                (df_temp['Experimental System Type'] == physical_type) &
                (df_temp['Experimental System'].isin(available_systems)) &
                (~df_temp['Experimental System'].isin(exclude_systems)) &
                (df_temp['Throughput'] == low_throughput)
            ]
        
            # 提取 A/B ID + Name
            pairs = df_filtered[['BioGRID ID Interactor A', 
                             'BioGRID ID Interactor B', 'Experimental System']].drop_duplicates()
            pairs.columns = ['Bio_A', 'Bio_B', 'Experimental System']
            all_data.append(pairs)
            print(f"  Filtered pairs: {len(pairs)}")
        except Exception as e:
            print(f"  Error: {e}")
    
    #合并id
    df_all = pd.concat(all_data, ignore_index=True)
    df_all['Bio_id'] = df_all.apply(lambda row: str(row['Bio_A']) + '_' + str(row['Bio_B']), axis=1)
    
    # 消除顺序影响（A_B=B_A）
    df_all['normalized_key'] = df_all['Bio_id'].apply(
    lambda x: '_'.join(sorted(str(x).split('_')))
    )

    #去重
    df_deduplicated = df_all.drop_duplicates(subset=['normalized_key']).drop(columns=['normalized_key'])

    #保存结果
    df_deduplicated.to_csv('high_confidence_Bio_ppi.csv', sep='\t', index=False)
    return df_deduplicated

#2、序列映射
def get_id_seq(fasta_file='uniprot_sprot.fasta', id_file="bio_uniprot.csv"):
    """
    从 FASTA 文件将 UniProt ID 映射到蛋白质序列。
    该函数解析 FASTA 文件创建 UniProt ID 到序列的字典，
    将其与映射文件中的 BioGRID ID 匹配，将结果保存到 'bio_uniprot_seq.csv',
    并将未匹配的 ID 记录到 'undo_ids_bio.txt'。
    """
    # 处理fasta文件
    seqdct = SeqIO.to_dict(SeqIO.parse(fasta_file, format='fasta'))
    filter_seqdct = {}
    for k, v in seqdct.items():
        uniport_id = k.split('|')[1]
        filter_seqdct[uniport_id] = str(v.seq)

    # 处理csv文件
    bio2uni = pd.read_csv(id_file)
    res = {}
    undo_ids = set()  # 用于记录无法匹配的ID
    
    for row in tqdm(bio2uni.itertuples(), desc='extract seq...'):
        uniprot_id = row.Uniprot_id
        bio_id = row.Bio_id
        
        try:
            # 尝试获取序列
            seq = filter_seqdct[uniprot_id]
            res[(int(bio_id), str(uniprot_id))] = seq
        except KeyError:
            # 记录无法找到的ID
            undo_ids.add(f"{bio_id}, {uniprot_id}")

    # 保存匹配结果
    df = pd.DataFrame(res, index=['seq']).T
    df.to_csv('bio_uniprot_seq.csv', sep='\t')

    # 保存未匹配的ID
    with open('undo_ids_bio.txt', 'w') as f:
        f.write('bio_id\t uniprot_id\n')
        for id_info in undo_ids:
            f.write(f'{id_info}\n')
    
    print(f"处理完成，成功匹配 {len(res)} 条记录，{len(undo_ids)} 条记录未匹配")
    return res, undo_ids

#3、PPI_SEQ整合
def PPI_SEQ(ppi_file='high_confidence_Bio_ppi.csv',id_seq_file='bio_uniprot_seq.csv'):
    """
    将 PPI 配对与其对应的 UniProt ID 和序列整合。
    该函数读取 PPI 和序列映射文件，为每个交互物配对匹配序列，
    将成功映射保存到 'success_ppi_bio.csv'（包括长度），并将失败记录到 'failed_ppi_bio.txt'。
    """
    #读取文件
    ppi_df=pd.read_csv(ppi_file,sep='\t')
    pep_df=pd.read_csv(id_seq_file,sep='\t',skiprows=1,names=['Bio_id','Uniprot_id','sequence'])

    #提取id与seq
    ID_A_list=ppi_df.iloc[:,0].astype(str).tolist()
    ID_B_list=ppi_df.iloc[:,1].astype(str).tolist()
    lookup_dict=dict(zip(pep_df['Bio_id'].astype(str),zip(pep_df['Uniprot_id'].astype(str),pep_df['sequence'])))

    #映射整合
    with open('success_ppi_bio.csv','w') as success_file, \
        open('failed_ppi_bio.txt','w') as failed_file:
        success_file.write('Bio_id,Uniprot_id,len_A,len_B,seq_A,seq_B\n')

        for i in range(len(ID_A_list)):
            ID_A=ID_A_list[i]
            ID_B=ID_B_list[i]
            Bio_id=f'{ID_A}_{ID_B}'

            if ID_A in lookup_dict and ID_B in lookup_dict:
                uniprot_a,seq_a=lookup_dict[ID_A]
                uniprot_b,seq_b=lookup_dict[ID_B]
                Uniprot_id=f'{uniprot_a}_{uniprot_b}'
                len_A=len(seq_a)
                len_B=len(seq_b)
                line=f'{Bio_id},{Uniprot_id},{len_A},{len_B},{seq_a},{seq_b}'
                success_file.write(line + '\n')
            else:
                failed_file.write(Bio_id + '\n')
    return success_file,failed_file

#4、长度（1024）过滤
def filter_len(ppi_seq_path='success_ppi_bio.csv'):
    """
    按长度阈值过滤 PPI 序列配对。
    该函数读取集成的 PPI 序列文件，并过滤掉任一序列超过 1024 个氨基酸的配对，
    将过滤结果保存到 'ppi_seq_bio.csv'。
    """
    #读取文件
    success_ppi_df=pd.read_csv(ppi_seq_path,sep=',')

    #过滤
    filtered_success_df=success_ppi_df[(success_ppi_df['len_A']<=1024)&(success_ppi_df['len_B']<=1024)]

    #保存结果
    filtered_success_df.to_csv('ppi_seq_bio.csv',sep=",",index=False)
    return filtered_success_df