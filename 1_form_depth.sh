#!/usr/bin/env bash
# ================================================================================
# Script: 1_form_depth.sh 
# Description: 基因组质控、比对、过滤、排序流程
# Requirements:
#   - bash >=4.0
#   - bwa, samtools, igvtools, bbduk.sh
#   - conda 环境 vip-seq
# Usage:
#   bash 1_form_depth.sh \
#     -g <host_genome.fna> \
#     -1 <raw_read1.fastq> \
#     -2 <raw_read2.fastq> \
#     [-p <phage.fna>] \
#     [-e <conda_env_path>] \
#     [-t <threads>] \
#     [-o <output_prefix>]
# ================================================================================

set -euo pipefail
# 设置 bbduk.sh 的路径
export PATH=~/anaconda3/envs/ProAct/bin:$PATH

# 默认参数
CONDA_ENV=""
THREADS=4
OUTPUT_PREFIX="pipeline"
PHAGE_FASTA=""

# 打印帮助信息
usage() {
  cat <<EOF
Usage: bash $0 \
  -g host_genome.fna \
  -1 raw_R1.fastq \
  -2 raw_R2.fastq  \
  [-e conda_env_path] \
  [-t threads] \
  [-o output_prefix]

Options:
  -g  参考基因组 FASTA 文件路径
  -1  原始测序 Read1 FASTQ
  -2  原始测序 Read2 FASTQ
  -e  conda 环境路径，默认: $CONDA_ENV
  -t  并行线程数，默认: $THREADS
  -o  所有输出文件前缀，默认: $OUTPUT_PREFIX
  -h  显示本帮助信息
EOF
  exit 1
}

# 参数解析
tmp_opts=$(getopt -o g:1:2:p:e:t:o:h --long genome:,read1:,read2:,phage:,env:,threads:,out:,help -n "$0" -- "$@")
if [ $? != 0 ]; then usage; fi
eval set -- "$tmp_opts"
while true; do
  case "$1" in
    -g|--genome)    HOST_GENOME="$2"; shift 2;;
    -1|--read1)     RAW_R1="$2";      shift 2;;
    -2|--read2)     RAW_R2="$2";      shift 2;;
    -e|--env)       CONDA_ENV="$2";    shift 2;;
    -t|--threads)   THREADS="$2";      shift 2;;
    -o|--out)       OUTPUT_PREFIX="$2"; shift 2;;
    -h|--help)      usage;              shift;;
    --) shift; break;;
    *) echo "Unknown option: $1"; usage;;
  esac
done

# 检查必填参数
if [ -z "${HOST_GENOME:-}" ] || [ -z "${RAW_R1:-}" ] || [ -z "${RAW_R2:-}" ]; then
  echo "Error: -g, -1 and -2 are required." >&2
  usage
fi

# 日志文件
LOG_OUT="${OUTPUT_PREFIX}.out"
LOG_ERR="${OUTPUT_PREFIX}.err"

# 生成样本清单（包含读段路径）
LIST_TXT="${OUTPUT_PREFIX}_input.txt"
printf "#sample_prefix\tref_genome\tread1\tread2\n" > "$LIST_TXT"
SAMPLE_PREFIX=$(basename "$RAW_R1" | sed -E 's/_R1.*//')
printf "%s\t%s\t%s\t%s\n" "$SAMPLE_PREFIX" "$HOST_GENOME" "$RAW_R1" "$RAW_R2" >> "$LIST_TXT"

# 激活 conda 环境
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

START_TIME=$(date +%s)
echo "开始时间: $(date)" 

# 主循环：读取四列（样本、基因组、R1、R2）
grep -v '^#' "$LIST_TXT" | while IFS=$'\t' read -r sample REF IN1 IN2; do
  echo "--- Processing sample: $sample ---"

  # 阶段 1：质控
  echo "[QC] Running bbduk for $sample"
  ~/anaconda3/envs/ProAct/bin/bbduk.sh in1="$IN1" in2="$IN2" \
           out1="${sample}_1_clean.fastq.gz" \
           out2="${sample}_2_clean.fastq.gz" \
           ref=adapters,phix trimq=14 maq=20 maxns=1 minlen=45

  # 阶段 2：比对
  echo "[ALIGN] Index & align to $REF"
  [[ -f "${REF}.bwt" ]] || bwa index "$REF"
  bwa mem -a -t "$THREADS" "$REF" \
      "${sample}_1_clean.fastq.gz" \
      "${sample}_2_clean.fastq.gz" | \
    samtools view -b -F 4 -o "${sample}_temp.bam"

  # 阶段 3：过滤
  echo "[FILTER] Filtering alignments"
  TEMP_BAM="${sample}_temp.bam"
  FILTERED_BAM="${sample}_filtered.bam"
  python3 <<EOF
import pysam
min_al, min_id, min_cov = 45, 0.97, 0.80
in_bam = pysam.AlignmentFile("${TEMP_BAM}", "rb")
out_bam = pysam.AlignmentFile("${FILTERED_BAM}", "wb", template=in_bam)
for read in in_bam:
    if read.is_unmapped:
        continue
    stats = read.get_cigar_stats()[0]
    al_len = sum(stats[:3])
    q_len = read.infer_read_length() or 0
    if al_len < min_al or q_len <= 0:
        continue
    cov = (stats[0] + stats[1]) / q_len
    nm = read.get_tag("NM") if read.has_tag("NM") else al_len
    if (al_len - nm) / al_len >= min_id and cov >= min_cov:
        out_bam.write(read)
in_bam.close()
out_bam.close()
EOF

  # 阶段 4：排序
  echo "[SORT] Sorting & generating TDF"
  samtools sort -@ "$THREADS" "${sample}_filtered.bam" -o "${sample}_sorted.bam"
  samtools depth -aa "${sample}_sorted.bam" > "./${sample}.depth"
  #igvtools count "${sample}_sorted.bam" "${sample}_sorted.tdf" "$REF"

  # 阶段 5：清理
  echo "[CLEAN] Removing intermediate files"
  rm -f "${sample}_temp.bam" "${sample}_filtered.bam" "${sample}_sorted.bam" \
        "${sample}_1_clean.fastq.gz" "${sample}_2_clean.fastq.gz" \
        ${REF}.{amb,ann,bwt,pac,sa} || true

done
