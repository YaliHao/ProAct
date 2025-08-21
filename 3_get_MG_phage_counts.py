#!/usr/bin/env python3
"""
Script: 3_get_MG_phage_counts.py 
Description:
  对单个样本执行基因与噬菌体区段深度统计。

Usage:
  python 3_get_MG_phage_counts.py  \
    phage_info.txt \
    gene_annotation.tsv \
    sample.depth \
    ./result_dir
"""
import os
import sys
import pandas as pd

if len(sys.argv) != 5:
    print(__doc__)
    sys.exit(1)

phage_info_path, gene_path, depth_path, output_dir = sys.argv[1:]

# 读取 phage_info
phage_df = pd.read_csv(
    phage_info_path,
    sep="\t", header=None,
    names=["srr","host_genome","chromosome","positions","phage_ids"]
)

# 确保输出目录存在
os.makedirs(output_dir, exist_ok=True)

def load_depth(path):
    df = pd.read_csv(path, sep="\t", header=None,
                     names=['chromosome','position','depth'], dtype={'chromosome':str})
    df = df.dropna(subset=['position','depth'])
    df['position'] = pd.to_numeric(df['position'], errors='coerce')
    df['depth']    = pd.to_numeric(df['depth'],    errors='coerce')
    return df.dropna()

# 读取 depth 与 gene annotation
depth_df = load_depth(depth_path)
gene_df  = pd.read_csv(gene_path, sep="\t")

gene_stats = []
phage_stats = []

# 基因统计
for _, gr in gene_df.iterrows():
    gid    = gr['Gene Id']
    chrom  = gid.split('_')[0]
    start  = int(gr['Start Position'])
    stop   = int(gr['Stop Position'])
    seg    = depth_df[
                (depth_df['chromosome']==str(chrom)) &
                (depth_df['position']>=start) &
                (depth_df['position']<=stop)
             ]
    tot    = seg['depth'].sum()
    length = stop - start + 1
    per    = tot/length if length>0 else 0
    med    = seg['depth'].median() if not seg.empty else 0
    gene_stats.append({
        'Gene Id':       gid,
        'Total_Counts':  tot,
        'Per_Counts':    per,
        'Median_Depth':  med,
        'Region_Length': length
    })

# 噬菌体统计（假设 phage_info.txt 只有一行或多行，均可）
for _, pr in phage_df.iterrows():
    chrom0   = pr['chromosome']
    pos_list = str(pr['positions']).split(',')
    id_list  = str(pr['phage_ids']).split(',')
    for pid, pos in zip(id_list, pos_list):
        try:
            st, sp = map(int, pos.split('-'))
        except ValueError:
            continue
        seg    = depth_df[
                    (depth_df['chromosome']==str(chrom0)) &
                    (depth_df['position']>=st) &
                    (depth_df['position']<=sp)
                 ]
        tot    = seg['depth'].sum()
        length = sp - st + 1
        per    = tot/length if length>0 else 0
        med    = seg['depth'].median() if not seg.empty else 0
        phage_stats.append({
            'Phage_Id':      pid,
            'Chromosome':    chrom0,
            'Start':         st,
            'Stop':          sp,
            'Total_Counts':  tot,
            'Per_Counts':    per,
            'Median_Depth':  med,
            'Region_Length': length
        })

# 写出结果
gene_out  = os.path.join(output_dir, f"marker_gene_counts.tsv")
phage_out = os.path.join(output_dir, f"phage_counts.tsv")

pd.DataFrame(gene_stats).to_csv(gene_out, sep="\t", index=False)
pd.DataFrame(phage_stats).to_csv(phage_out, sep="\t", index=False)

# —— 新增板块：计算 host 基因深度中位数 ——
# 从 gene_stats 中提取第二列 Total_Counts 的中位数
gene_df_out = pd.DataFrame(gene_stats)
median_total = gene_df_out['Per_Counts'].median()

# 构造 Sample_ID: host_genome--srr
host  = phage_df.iloc[0]['host_genome']
srr   = phage_df.iloc[0]['srr']
sample_id = f"{host}--{srr}"

# 写出 host counts
host_out = os.path.join(output_dir, f"host_counts.tsv")
host_counts_df = pd.DataFrame([{
    'Sample_ID':      sample_id,
    'Median_of_MG':   median_total
}])
host_counts_df.to_csv(host_out, sep="\t", index=False)

print(f"Done. Gene counts -> {gene_out}\n     Phage counts -> {phage_out}\n      Host median -> {host_out}")

