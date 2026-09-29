#!/usr/bin/env bash
# ============================================================
# Step 1: 从总表提取各氨基酸类型的 tRNA 序列
# 输入: all.RNA_features.1211.tsv
# 输出: 各氨基酸目录下的 tRNA.fa (FASTA 格式)
# ============================================================
set -euo pipefail

source "$(dirname "$0")/config.sh"

echo "=========================================="
echo "Step 1: 提取 tRNA 序列"
echo "=========================================="

# 检查输入文件
if [[ ! -f "${FEATURES_TSV}" ]]; then
    echo "错误: 输入文件不存在: ${FEATURES_TSV}" >&2
    exit 1
fi

# TSV 列: $1=ID, $2=RNA_Subtype(氨基酸), $3=RNA_Feature(密码子), $4=Species, $5=Locations, $6=Sequence
while read -r aa; do
    [[ -z "${aa}" ]] && continue
    mkdir -p "${aa}"
    echo "  提取 ${aa} tRNA 序列..."
    awk -F'\t' -v aa="${aa}" \
        '$2 == aa { printf ">%s\n%s\n", $1, $6 }' \
        "${FEATURES_TSV}" > "${aa}/tRNA.fa"
    nseq=$(grep -c '^>' "${aa}/tRNA.fa" || true)
    echo "    ${aa}: ${nseq} 条序列"
done < "${AA_LIST}"

echo ""
echo "Step 1 完成。"