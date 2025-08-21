#!/usr/bin/env python3
"""
Script:4_calculate_PtoH.py

读取：
  1) phage_counts.tsv（噬菌体区段深度平均值）
  2) host_counts.tsv（host 基因深度中位数）

计算：
  每条 phage 记录的 PtoH = Per_Counts / Median_of_host，
  Quality = 若 Per_Counts<10 或 Median_of_MG<10 则 "low"，否则 "high"。

输出：
  PtoH.tsv，包含列：
    Sample_ID, Phage_Id, Chromosome, Start, Stop,
    Total_Counts, Per_Counts, Median_of_MG, PtoH, Quality
"""

import sys
import os
import pandas as pd

def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <phage_counts.tsv> <host_counts.tsv> <output.tsv>")
        sys.exit(1)

    phage_file, host_file, output_file = sys.argv[1:]

    # 读取输入
    phage_df = pd.read_csv(phage_file, sep="\t")
    host_df = pd.read_csv(host_file, sep="\t")

    # 仅取第一条 host 记录
    sample_id   = host_df.loc[0, "Sample_ID"]
    median_mg   = host_df.loc[0, "Median_of_MG"]

    # 扩展到每条 phage 记录
    phage_df["Sample_ID"]     = sample_id
    phage_df["Median_of_MG"]  = median_mg
    phage_df["PtoH"]          = phage_df["Per_Counts"] / phage_df["Median_of_MG"]

    # 计算 Quality
    phage_df["Quality"] = phage_df.apply(
        lambda r: "low" if (r["Per_Counts"] < 10 or r["Median_of_MG"] < 10) else "high",
        axis=1
    )
    # 计算 Activity
    phage_df["Activity"] = phage_df["PtoH"].apply(
        lambda x: "active" if x >= 1.5 else ("inactive" if x < 1 else "low")
    )

    # 重新排序列
    cols = [
        "Sample_ID", "Phage_Id", "Chromosome", "Start", "Stop",
        "Total_Counts", "Per_Counts", "Median_of_MG", "PtoH", "Quality"
    ]
    out_df = phage_df[cols]

    # 确保输出目录存在
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    out_df.to_csv(output_file, sep="\t", index=False)
    print(f"Written PtoH table to {output_file}")

if __name__ == "__main__":
    main()

