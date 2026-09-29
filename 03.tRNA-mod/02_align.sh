#!/usr/bin/env bash
# ============================================================
# Step 2: 对某氨基酸下所有 tRNA 并行做 MAFFT 比对
# 用法: 02_align.sh <amino_acid>
# 依赖: seqkit, MAFFT (linsi)
# ============================================================
set -euo pipefail

source "$(dirname "$0")/config.sh"

aa="$1"
ALIGN_SCRIPT="${WORK_DIR}/02_align_single.sh"

echo "=========================================="
echo "Step 2: ${aa} tRNA 比对"
echo "=========================================="

cd "${aa}"

# 提取所有 tRNA ID，跳过 ref.fa
nseq=$(grep -c '^>' tRNA.fa || true)
echo "  ${aa}: 共 ${nseq} 条 tRNA，${NPROC_ALIGN} 线程并行比对..."

seqkit fx2tab -n tRNA.fa | \
    xargs -I{} -P "${NPROC_ALIGN}" bash "${ALIGN_SCRIPT}" {}

nafa=$(ls -1 *.afa 2>/dev/null | wc -l)
echo "  ${aa}: 生成 ${nafa} 个比对文件"

cd "${WORK_DIR}"
echo ""