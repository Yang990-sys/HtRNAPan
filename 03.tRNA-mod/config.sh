#!/usr/bin/env bash
# ============================================================
# tRNA 修饰位点映射 Pipeline - 配置文件
# ============================================================

# ---- 输入文件 ----
# 待预测的 tRNA 序列总表（含 tRNA_ID, amino_acid, sequence 等列）
FEATURES_TSV="../../new_data/all.RNA_features.1211.tsv"
# 参考序列（无修饰的标准碱基 FASTA，923 条）
REF_FA="unmodified.fa"
# 已知修饰位点表（sequence_id, position, modification, base）
MOD_SITES_TSV="modification_sites.tsv"

# ---- 工作目录 ----
WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${WORK_DIR}/results"

# ---- 比对参数 ----
ALIGNER="linsi"           # MAFFT 算法：linsi (本地比对，适合 tRNA)
NPROC_ALIGN=30            # 单氨基酸并行比对线程数

# ---- 修饰映射参数 ----
NPROC_MOD_MAP=32          # map_modifications.py 并行进程数
BATCH_SIZE=5000            # 每批处理文件数

# ---- 氨基酸列表 ----
AA_LIST="${WORK_DIR}/aalist"