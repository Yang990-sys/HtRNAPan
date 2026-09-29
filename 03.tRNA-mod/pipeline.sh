#!/usr/bin/env bash
# ============================================================
# tRNA 修饰位点映射 Pipeline - 主入口
#
# 流程:
#   Step 1: 从总表提取各氨基酸 tRNA 序列
#   Step 2: 对每条 tRNA 与参考序列做 MAFFT 比对
#   Step 3: 将已知修饰位点映射到比对结果，统计支持率
#
# 用法:
#   bash pipeline.sh              # 运行全部步骤
#   bash pipeline.sh --from 2     # 从 Step 2 开始（跳过提取）
#   bash pipeline.sh --only 3     # 只运行 Step 3（比对已完成）
#   bash pipeline.sh --aa Gly     # 只处理指定氨基酸
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# ---- 解析参数 ----
FROM_STEP=1
ONLY_STEP=0
TARGET_AA=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --from)  FROM_STEP="$2"; shift 2 ;;
        --only)  ONLY_STEP="$2"; shift 2 ;;
        --aa)    TARGET_AA="$2";  shift 2 ;;
        *)       echo "未知参数: $1" >&2; exit 1 ;;
    esac
done

# ---- 确定要处理的氨基酸列表 ----
if [[ -n "${TARGET_AA}" ]]; then
    AA_LIST_PROC=(${TARGET_AA})
else
    mapfile -t AA_LIST_PROC < "${AA_LIST}"
fi

echo "######################################################"
echo "#  tRNA 修饰位点映射 Pipeline"
echo "######################################################"
echo "  工作目录: ${WORK_DIR}"
echo "  氨基酸:   ${AA_LIST_PROC[*]}"
echo "  起始步骤: ${FROM_STEP}"
echo ""

# ---- Step 1: 提取 tRNA 序列 ----
if [[ "${ONLY_STEP}" -eq 0 || "${ONLY_STEP}" -eq 1 ]] && [[ "${FROM_STEP}" -le 1 ]]; then
    bash "${SCRIPT_DIR}/01_extract_tRNA.sh"
fi

# ---- Step 2: MAFFT 比对 ----
if [[ "${ONLY_STEP}" -eq 0 || "${ONLY_STEP}" -eq 2 ]] && [[ "${FROM_STEP}" -le 2 ]]; then
    for aa in "${AA_LIST_PROC[@]}"; do
        [[ -z "${aa}" ]] && continue
        bash "${SCRIPT_DIR}/02_align.sh" "${aa}"
    done
fi

# ---- Step 3: 修饰位点映射 ----
if [[ "${ONLY_STEP}" -eq 0 || "${ONLY_STEP}" -eq 3 ]] && [[ "${FROM_STEP}" -le 3 ]]; then
    for aa in "${AA_LIST_PROC[@]}"; do
        [[ -z "${aa}" ]] && continue
        bash "${SCRIPT_DIR}/03_map_modifications.sh" "${aa}"
    done
fi

echo "######################################################"
echo "#  Pipeline 全部完成！"
echo "#  结果目录: ${RESULTS_DIR}"
echo "######################################################"