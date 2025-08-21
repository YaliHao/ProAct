#!/usr/bin/env bash
# ================================================================================
# Script: run_pipeline.sh
# Description:
#   一条命令执行四步：生成 depth、注释 marker genes、统计 MG/Phage 深度、计算 PtoH。
#   所有输出都放在 -o 指定的目录
#
# Usage:
#   bash run_pipeline.sh \
#     -g host_genome.fna \
#     -1 reads_R1.fastq \
#     -2 reads_R2.fastq \
#     -r contig,coords  \# 可多次指定，如 -r CP024621.1,292609-333351 \  \
#     [-e conda_env] \# 可选：指定 conda 环境名 \
#     [-t threads] \# 可选：map threads，默认8 \
#     [-c cpus]    \# 可选：GTDB-Tk CPU，默认4 \
#     -o outdir    \# 输出目录（必选）
#
# 输出目录结构：
#   outdir/
#     sample.depth
#     sample_phage_info.txt
#     MG/*.tsv
#     counts/*_gene_counts.tsv, *_phage_counts.tsv, *_host_counts.tsv
#     PtoH.tsv
#
# Dependencies:
#   bash>=4.0, bwa, samtools, bbmap(bbduk.sh), igvtools, gtdbtk,
#   python3 (pandas, biopython, pysam)
# ================================================================================
set -euo pipefail

# 默认值
CONDA_ENV=""
THREADS=8
CPUS=4
OUTPUT_DIR=""
declare -a REGIONS=()

usage() {
  grep '^#' "$0" | sed 's/^#//'
  exit 1
}

# 解析参数
while [[ $# -gt 0 ]]; do
  case "$1" in
    -g) HOST="$2"; shift 2;;
    -1) R1="$2"; shift 2;;
    -2) R2="$2"; shift 2;;
    -r) REGIONS+=("$2"); shift 2;;
    -e) CONDA_ENV="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    -c) CPUS="$2"; shift 2;;
    -o) OUTPUT_DIR="$2"; shift 2;;
    -h) usage;;
    *) echo "Unknown option: $1"; usage;;
  esac
done

# 校验必填
if [[ -z "${HOST:-}" || -z "${R1:-}" || -z "${R2:-}" || -z "${OUTPUT_DIR:-}" || ${#REGIONS[@]} -eq 0 ]]; then
  echo "Error: -g, -1, -2, -r (>=1次) and -o are required." >&2
  usage
fi

# 激活 conda 环境（可选）
eval "$(conda shell.bash hook)"
conda activate "$(conda info --envs | grep '*' | cut -d ' ' -f 1)"

# 推断 sample 名称
sample=$(basename "$R1" .fastq)
sample=${sample%_R1}

# 创建输出目录
if [[ ! -d "$OUTPUT_DIR" ]]; then
  mkdir -p "$OUTPUT_DIR"
fi

# 生成 phage_info.txt
PHAGE_INFO="$OUTPUT_DIR/${sample}_phage_info.txt"
> "$PHAGE_INFO"
declare -A COORDS MAP_IDS
count=1
for entry in "${REGIONS[@]}"; do
  contig=${entry%%,*}
  coords=${entry#*,}
  id="P${count}"
  ((count++))
  if [[ -z "${COORDS[$contig]:-}" ]]; then
    COORDS[$contig]="$coords"
    MAP_IDS[$contig]="$id"
  else
    COORDS[$contig]+=",$coords"
    MAP_IDS[$contig]+=",$id"
  fi
done
for contig in "${!COORDS[@]}"; do
  echo -e "${sample}\t${HOST}\t${contig}\t${COORDS[$contig]}\t${MAP_IDS[$contig]}" >> "$PHAGE_INFO"
done

echo "[+] Generated phage info: $PHAGE_INFO"

# 步骤 1：生成 depth
echo "[1] Generating depth"
bash 1_form_depth.sh -g "$HOST" -1 "$R1" -2 "$R2" -t "$THREADS" -o "$sample"
mv "${sample}.depth" "$OUTPUT_DIR/"
DEPTH_FILE="$OUTPUT_DIR/${sample}.depth"
echo "    -> Depth file: $DEPTH_FILE"

# 步骤 2：GTDB-Tk 注释 marker genes
echo "[2] Annotating marker genes"
MG_DIR="$OUTPUT_DIR/MG"
mkdir -p "$MG_DIR"
bash 2_form_MG.sh -g "$HOST" -o "$MG_DIR" -c "$CPUS"
GENE_TSV="$MG_DIR/$(basename "$HOST" .fna).tsv"
echo "    -> Gene annotation TSV: $GENE_TSV"

# 步骤 3：统计 MG & Phage 深度
echo "[3] Counting MG & Phage depths"
COUNT_DIR="$OUTPUT_DIR/counts"
mkdir -p "$COUNT_DIR"
# 调用统计脚本
BASE=$(basename "$PHAGE_INFO" .txt)
python3 3_get_MG_phage_counts.py \
  "$PHAGE_INFO" \
  "$GENE_TSV" \
  "$DEPTH_FILE" \
  "$COUNT_DIR"
echo "    -> Counts in: $COUNT_DIR/*"

# 步骤 4：计算 PtoH
echo "[4] Calculate PtoH"
PTOH_FILE="$OUTPUT_DIR/PtoH.tsv"
python3 4_calculate_PtoH.py \
  "$COUNT_DIR/phage_counts.tsv" \
  "$COUNT_DIR/host_counts.tsv" \
  "$PTOH_FILE"
echo "    -> PtoH result: $PTOH_FILE"

# 完成
echo "=== ProAct completed successfully ==="
echo "Output directory: $OUTPUT_DIR"
