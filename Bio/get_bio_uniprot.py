import multiprocessing as mp
import pandas as pd
from tqdm import tqdm
import os
import requests
import re

#将Bio_id映射到Uniprot_id
#由于涉及网络下载，该代码生成文件bio_uniprot.csv，failed_bio_ids.txt已提供在文件夹里

#1、提取Bio_id集合
def unique_bio_id(input='high_confidence_Bio_ppi.csv'):
    """
    从过滤后的 PPI 数据中提取唯一的 BioGRID ID,减少下载量,
    该函数读取PPI CSV文件,将Bio_A和Bio_B列合并去重,
    并将唯一 ID 保存到 'unique_bio_id.csv'。
    """
    df=pd.read_csv(input,sep='\t')

    #提取ID_A和ID_B列并合并
    combined_ids=pd.concat([df['Bio_A'],df['Bio_B']])

    #去重
    unique_ids=combined_ids.unique()

    #保存结果
    unique_ids_df=pd.DataFrame(unique_ids,columns=['Bio_id'])
    unique_ids_df.to_csv('unique_bio_id.csv', sep='\t', index=False)
    return unique_ids_df

#2、Uniprot_id映射
def get_Uniprot_id(BIO_ID_FILE='unique_bio_id.csv',CPU_COUNT=1):
    """
    使用多进程将 BioGRID ID 映射到 UniProt ID。
    该函数为每个 BioGRID ID 从 BioGRID 网页获取 UniProt 访问号，
    具有处理中断运行的恢复功能，以指定的 CPU 数量并行处理，
    并将结果保存到 'bio_uniprot.csv'，同时将失败记录到 'failed_bio_ids.txt'。
    """
    # 配置参数
    #BIO_ID_FILE 存储Bio_id的文件
    OUTPUT_FILE = "bio_uniprot.csv"  # 存储结果的文件，格式：Bio_id,uniprot_id
    FAILED_FILE = "failed_bio_uniprot.txt"  # 存储下载失败的Bio_id
    RESUME_FILE = "resume_point.txt"  # 存储断点续传信息
    # CPU_COUNT使用的CPU核心数

    def get_uniprot_id(bio_id):
        """根据Bio_id获取UniProt ID"""
        try:
            # 步骤1：请求BioGRID页面，提取UniProt accession ID
            biogrid_url = f"https://thebiogrid.org/{bio_id}/summary/anas-platyrhynchos/ddx58.html"
            headers = {
                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
            }
            response = requests.get(biogrid_url, headers=headers, timeout=10)
            if response.status_code != 200:
                return bio_id, None
        
            # 提取UniProt ID
            uniprot_match = re.search(r"uniprot.org/uniprot/(\w+)", response.text)
            if not uniprot_match:
                return bio_id, None
            uniprot_id = uniprot_match.group(1)
        
            return bio_id, uniprot_id
        except Exception as e:
            print(f"处理Bio_id {bio_id} 时出错: {str(e)}")
            return bio_id, None
        
    def main():
        # 检查断点续传文件
        resume_point = 0
        if os.path.exists(RESUME_FILE):
            with open(RESUME_FILE, "r") as f:
                resume_point = int(f.read().strip())
    
        # 读取Bio_id列表
        with open(BIO_ID_FILE, "r") as f:
            bio_ids = [line.strip() for line in f.readlines()][resume_point:]
    
        # 初始化多进程池
        pool = mp.Pool(CPU_COUNT)
        results = []
        failed_bio_ids = []
    
        # 带进度条的多进程处理
        for bio_id, uniprot_id in tqdm(pool.imap_unordered(get_uniprot_id, bio_ids), 
                                        total=len(bio_ids), initial=resume_point, 
                                        desc="获取UniProt ID进度"):
            if uniprot_id:
                results.append(f"{bio_id},{uniprot_id}")
            else:
                failed_bio_ids.append(bio_id)
        
            # 每处理100个条目，保存一次断点和结果
            if len(results) % 100 == 0 or len(bio_ids) - (len(results) + len(failed_bio_ids)) == 0:
                # 保存结果
                with open(OUTPUT_FILE, "a") as f:
                    f.write("\n".join(results) + "\n")
                results = []
            
                # 保存断点
                with open(RESUME_FILE, "w") as f:
                    f.write(str(resume_point + len(results) + len(failed_bio_ids)))
    
        # 处理剩余结果
        if results:
            with open(OUTPUT_FILE, "a") as f:
                f.write("\n".join(results) + "\n")
    
            # 保存失败的Bio_id
            with open(FAILED_FILE, "w") as f:
                f.write("\n".join(failed_bio_ids))
    
        pool.close()
        pool.join()
        print("批量获取完成！")
        print(f"成功获取 {len(bio_ids) - len(failed_bio_ids)} 个UniProt ID，{len(failed_bio_ids)} 个失败。")

    if __name__ == "__main__":
        main()
    return
