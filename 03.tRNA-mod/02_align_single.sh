#!/usr/bin/env bash
# ============================================================
# Step 2 子脚本: 对单条 tRNA 做 MAFFT 多序列比对
# 用法: 02_align_single.sh <tRNA_ID>
# 依赖: seqkit, MAFFT (linsi)
# ============================================================
set -euo pipefail

tRNA_ID="$1"

# 提取目标 tRNA 序列，拼接参考序列 ref.fa
seqkit grep -p "${tRNA_ID}" tRNA.fa | cat - ref.fa > "${tRNA_ID}.fasta"

# MAFFT 本地比对
linsi --quiet "${tRNA_ID}.fasta" > "${tRNA_ID}.afa"

# 清理中间文件
rm -f "${tRNA_ID}.fasta"