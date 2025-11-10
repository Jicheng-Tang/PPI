import pandas as pd
from Bio import SeqIO


#1、依据验证实验类型（高置信、低通量）筛选
def experiment_method(raw_data_path):
    """
    依据验证实验类型(高置信、低通量)筛选。
    该函数读取原始数据文件，根据高置信实验方法的 PSI-MI 编号和名称进行筛选，
    保留关键列(如蛋白 ID、实验方法),并保存结果到 'high_confidence_ppi.csv'。
    """
    #定义列名
    columns = [
    "protein1_uniprot", "protein2_uniprot", "protein1_db_ids", "protein2_db_ids",
    "protein1_annotation", "protein2_annotation", "experiment_method", "publication",
    "pub_ids", "taxid_protein1", "taxid_protein2", "interaction_type",
    "database_source", "interaction_ids", "confidence_score"
    ]
    df=pd.read_csv(
    raw_data_path,
    sep='\t',
    header=None,
    names=columns
    )

    # 定义高置信实验方法的PSI-MI编号和名称
    high_confidence_methods = [
    "MI:0107", "Surface plasmon resonance",
    "MI:0065", "Isothermal titration calorimetry",
    "MI:0077", "Nuclear magnetic resonance",
    "MI:0114", "X-ray crystallography",
    "MI:0872", "Atomic force microscopy",
    "MI:0019", "Immunoprecipitation",
    "MI:0096", "Pull down",
    "MI:0676", "Tandem affinity purification",
    "MI:0030", "Cross-linking", "MI:0031", "Cross-linking",
    "MI:0411", "ELISA",
    "MI:0809", "Bimolecular fluorescence complementation",
    "MI:1203", "Split luciferase complementation",
    "MI:0012", "Bioluminescence resonance energy transfer",
    "MI:0055", "Fluorescence resonance energy transfer",
    "MI:0016", "Circular dichroism",
    "MI:0038", "Dynamic light scattering",
    "MI:1311", "Differential scanning calorimetry",
    "MI:0966", "UV-visible spectroscopy"
    ]

    # 构建筛选条件：包含高置信方法的行
    include_high_confidence = df['experiment_method'].str.contains(
    '|'.join(high_confidence_methods),
    case=False,
    na=False
    )

    # 执行筛选，生成高置信实验的DataFrame
    high_conf_df = df[include_high_confidence].copy()

    # 保留关键列（如蛋白ID、实验方法、置信度等）
    key_columns = ["protein1_uniprot", "protein2_uniprot", "experiment_method"]
    high_conf_key_df = high_conf_df[key_columns]

    # 保存结果
    high_conf_key_df.to_csv('high_confidence_ppi.csv', sep='\t',index=False)
    return high_conf_key_df


#2、提取Uniprot_id
def Uniprot_id(input_path='high_confidence_ppi.csv'):
    """
    提取 UniProt ID。
    该函数读取高置信 PPI 文件，筛选前两列均包含 'uniprotkb:' 的行，
    提取并替换为纯 ID,保存结果到 'Extracted_id.csv'。
    """
    # 读取文件
    df = pd.read_csv(input_path, sep='\t')
    
    # 提取前两列均包含'uniprotkb:'的行
    filtered_df = df[df['protein1_uniprot'].str.contains('uniprotkb:', na=False) & 
                 df['protein2_uniprot'].str.contains('uniprotkb:', na=False)]
    
    #提取id
    filtered_df["A_id"] = filtered_df["protein1_uniprot"].str.replace("uniprotkb:", "")
    filtered_df["B_id"] = filtered_df["protein2_uniprot"].str.replace("uniprotkb:", "")
    
    #保存结果
    filtered_df.to_csv('Extracted_id.csv', sep='\t',index=False)
    return filtered_df


#3、去重
def rm_repetition(input_path='Extracted_id.csv'):
    """
    去除重复项。
    该函数读取提取的 ID 文件，合并 A_id 和 B_id 为 Uniprot_id,
    规范化顺序以消除 A_B = B_A 的影响，去重后保存到 'uniprot_id.csv'。
    """
    # 读取文件
    df = pd.read_csv(input_path, sep='\t')

    #合并id
    df['Uniprot_id'] = df.apply(lambda row: row['A_id'] + '_' + row['B_id'], axis=1)

    # 消除顺序影响（A_B=B_A）
    df['normalized_key'] = df['Uniprot_id'].apply(
    lambda x: '_'.join(sorted(str(x).split('_')))
    )

    #去重
    df_deduplicated = df.drop_duplicates(subset=['normalized_key']).drop(columns=['normalized_key'])

    #保存结果
    df_deduplicated.to_csv('uniprot_id.csv', sep='\t', index=False)
    return df_deduplicated


#4、提取ppi序列
def ppi_seq(seq_path, ppi_path='uniprot_id.csv'):
    """
    提取 PPI 序列。
    该函数读取 PPI 文件和 FASTA 序列文件，匹配 UniProt ID 与序列，
    保存成功匹配的 PPI 序列（包括长度）到 'success_ppi.csv',
    失败项保存到 'failed_ppi.csv'。
    """
    #读取ppi文件
    ppi_df=pd.read_csv(ppi_path,sep='\t')

    #提取id序列
    ID_A_list=ppi_df.iloc[:,3].astype(str).tolist()
    ID_B_list=ppi_df.iloc[:,4].astype(str).tolist()

    #从seq_path提取seq字典
    seqdct = SeqIO.to_dict(SeqIO.parse(seq_path, format='fasta'))
    lookup_dict = {}
    for k, v in seqdct.items():
        uniport_id = k.split('|')[1]
        lookup_dict[uniport_id] = str(v.seq)
    
    #提取ppi_seq
    with open('success_ppi.csv','w') as success_file, \
         open('failed_ppi.csv','w') as failed_file:
        success_file.write('Uniprot_id,len_A,len_B,seq_A,seq_B\n')
        failed_file.write('Uniprot_id\n')

        for i in range(len(ID_A_list)):
            ID_A=ID_A_list[i]
            ID_B=ID_B_list[i]
            Uniprot_id=f'{ID_A}_{ID_B}'
            
            if ID_A in lookup_dict and ID_B in lookup_dict:
                seq_a=lookup_dict[ID_A]
                seq_b=lookup_dict[ID_B]
                len_A=len(seq_a)
                len_B=len(seq_b)
                line=f'{Uniprot_id},{len_A},{len_B},{seq_a},{seq_b}'
                success_file.write(line + '\n')
            else:
                failed_file.write(Uniprot_id + '\n')
    return success_file,failed_file

#5、长度（1024）过滤
def filter_len(ppi_seq_path='success_ppi.csv'):
    """
    按长度(1024)过滤。
    该函数读取成功 PPI 序列文件，过滤 len_A 和 len_B 均 <= 1024 的配对，
    保存结果到 'ppi_seq_mint.csv'。
    """
    #读取文件
    success_ppi_df=pd.read_csv(ppi_seq_path,sep=',')

    #过滤
    filtered_success_df=success_ppi_df[(success_ppi_df['len_A']<=1024)&(success_ppi_df['len_B']<=1024)]

    #保存结果
    filtered_success_df.to_csv('ppi_seq_mint.csv',sep=",",index=False)
    return filtered_success_df