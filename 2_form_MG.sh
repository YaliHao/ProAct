#!/usr/bin/env bash
# ================================================================================
# Script:2_form_MG.sh 
# Description: 使用 GTDB-Tk 标识 marker genes 并整合 Pfam/TIGRFAM 注释及蛋白位点
# Usage:
#   bash 2_form_MG.sh \
#     -g host_genome.fna \
#     -o output_dir \
#     [-c cpus]
# Requirements:
#   - gtdbtk >= 2.x
#   - python3 环境（包含 pandas、biopython）
# ================================================================================

set -euo pipefail

# 默认参数
CPUS=4

# 帮助信息
usage() {
  cat <<EOF
Usage: bash $0 \
  -g host_genome.fna \
  -o output_dir \
  [-c cpus]

Options:
  -g  输入宿主基因组 FASTA 文件（.fna）路径
  -o  GTDB-Tk 输出目录
  -c  使用的 CPU 数，默认: $CPUS
  -h  查看帮助
EOF
  exit 1
}

# 解析参数
tmp=$(getopt -o g:o:c:h --long genome:,outdir:,cpus:,help -n "$0" -- "$@")
if [ $? != 0 ]; then usage; fi
eval set -- "$tmp"
while true; do
  case "$1" in
    -g|--genome)   HOST_GENOME="$2"; shift 2;;
    -o|--outdir)   OUTDIR="$2";     shift 2;;
    -c|--cpus)     CPUS="$2";       shift 2;;
    -h|--help)     usage;            shift;;
    --) shift; break;;
    *) echo "Unknown option $1"; usage;;
  esac
done

# 参数检查
if [ -z "${HOST_GENOME:-}" ] || [ -z "${OUTDIR:-}" ]; then
  echo "Error: -g and -o are required." >&2
  usage
fi

# 准备目录
gtdb_in=$(dirname "${HOST_GENOME}")
mkdir -p "${OUTDIR}"

echo "[GTDB-Tk] identify: genomes in $gtdb_in -> $OUTDIR"
gtdbtk identify --genome_dir "$gtdb_in" --out_dir "$OUTDIR" -x fna --cpus "$CPUS"

echo "[GTDB-Tk] removing .sha256 files"
find "$OUTDIR/identify/intermediate_results/marker_genes/" -name "*.sha256" -delete

echo "[GTDB-Tk] merging Pfam/TIGRFAM annotation and protein coordinates"
for dir in "$OUTDIR"/identify/intermediate_results/marker_genes/*/; do
  file_id=$(basename "$dir")
  pfam_tsv="$dir/${file_id}_pfam_tophit.tsv"
  tigr_tsv="$dir/${file_id}_tigrfam_tophit.tsv"
  prot_fna="$dir/${file_id}_protein.fna"
  out_tsv="$OUTDIR/${file_id}.tsv"

  echo "  - processing $file_id"
  python3 <<EOF
import os
import pandas as pd
from Bio import SeqIO

def load_tophit(path):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=['Gene Id','Family Id'])
    df = pd.read_csv(path, sep="\t", comment='#', header=0)
    if df.shape[1] < 2:
        return pd.DataFrame(columns=['Gene Id','Family Id'])
    gene_col, hit_col = df.columns[0], df.columns[1]
    df = df[[gene_col, hit_col]].rename(columns={gene_col: 'Gene Id'})
    df['Family Id'] = df[hit_col].astype(str).apply(lambda s: s.split(',')[0] if isinstance(s, str) else '')
    return df[['Gene Id','Family Id']]

# 读取注释表
df1 = load_tophit(r"${pfam_tsv}")
df2 = load_tophit(r"${tigr_tsv}")
# 合并注释
df = pd.concat([df1, df2], ignore_index=True)

# 读取蛋白 FASTA 并获取位点
gene_info = []
if os.path.isfile(r"${prot_fna}"):
    for rec in SeqIO.parse(r"${prot_fna}", 'fasta'):
        parts = rec.description.split(' # ')
        if len(parts) >= 3:
            try:
                start = int(parts[1]); stop = int(parts[2])
                gene_info.append({'Gene Id': rec.id,
                                  'Start Position': start,
                                  'Stop Position': stop})
            except ValueError:
                pass
info_df = pd.DataFrame(gene_info)

# 合并并去重
if not df.empty and not info_df.empty:
    final = pd.merge(df, info_df, on='Gene Id').drop_duplicates('Gene Id')
elif not df.empty:
    final = df.copy()
    final['Start Position'] = None; final['Stop Position'] = None
elif not info_df.empty:
    final = info_df.copy()
    final['Family Id'] = None
else:
    final = pd.DataFrame(columns=['Gene Id','Family Id','Start Position','Stop Position'])

# 保存结果
os.makedirs(os.path.dirname(r"${out_tsv}"), exist_ok=True)
final.to_csv(r"${out_tsv}", sep="\t", index=False)
EOF

done

echo "[Done] annotations in $OUTDIR"

