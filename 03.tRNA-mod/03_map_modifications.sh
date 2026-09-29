#!/usr/bin/env bash
# ============================================================
# Step 3: 修饰位点映射与统计
# 用法: 03_map_modifications.sh <amino_acid>
# 依赖: Python 3, biopython, tqdm, numpy
# ============================================================
set -euo pipefail

source "$(dirname "$0")/config.sh"

aa="$1"

echo "=========================================="
echo "Step 3: ${aa} 修饰位点映射"
echo "=========================================="

mkdir -p "${RESULTS_DIR}"

python "${WORK_DIR}/map_modifications.py" \
    "${MOD_SITES_TSV}" \
    "${aa}" \
    "${RESULTS_DIR}/" \
    "${NPROC_MOD_MAP}"

echo "  ${aa} 修饰映射完成。"
echo ""