import pandas as pd

#将BIOGRID和MINT获得的ppi_seq合并去重
def combine_bio_mint(bio='Bio/ppi_seq_bio.csv',mint='mint/ppi_seq_mint.csv'):
    """
    将BIOGRID和MINT获得的ppi_seq合并去重,结果保存在combined_ppi_seq.csv中。
    """
    # 读取两个CSV文件
    bio_df = pd.read_csv(bio,sep=',')
    mint_df = pd.read_csv(mint,sep=',')

    # 定义要提取的列
    columns_to_extract = ['Uniprot_id', 'len_A', 'len_B', 'seq_A', 'seq_B']

    #从每个DataFrame中提取指定的列
    bio_subset = bio_df[columns_to_extract]
    mint_subset = mint_df[columns_to_extract]   

    # 将子集合并成一个DataFrame
    combined_df = pd.concat([bio_subset, mint_subset], ignore_index=True)

    # 消除顺序影响（A_B=B_A）
    combined_df['normalized_key'] =combined_df['Uniprot_id'].apply(
        lambda x: '_'.join(sorted(str(x).split('_')))
        )
    
    #去重
    df_deduplicated =combined_df.drop_duplicates(subset=['normalized_key']).drop(columns=['normalized_key'])

    #保存结果
    df_deduplicated.to_csv('combined_ppi_seq.csv',sep=',',index=False)
    return df_deduplicated

if __name__ == "__main__":
    combine_bio_mint()