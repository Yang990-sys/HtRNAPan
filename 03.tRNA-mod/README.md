# tRNA Modification Site Mapping Pipeline

将已知 tRNA 修饰位点通过多序列比对映射到待预测的 tRNA 序列上，统计每个位点的修饰支持率。

## 流程概览

```
all.RNA_features.tsv (待预测序列)
        │
        ▼ Step 1: 按氨基酸类型提取
  各氨基酸目录/tRNA.fa
        │
        ▼ Step 2: 与参考序列做 MAFFT 比对
  各氨基酸目录/tRNA_xxx.afa
        │
        ▼ Step 3: 映射修饰位点 + 统计
  results/tRNA_xxx_site_stats.tsv
```

## 输入文件

| 文件 | 说明 | 格式 |
|------|------|------|
| `all.RNA_features.tsv` | 待预测 tRNA 序列总表 | TSV: ID, RNA_Subtype, RNA_Feature, Species, Locations, Sequence, 2D |
| `unmodified.fa` | 参考序列（无修饰，标准碱基） | FASTA |
| `modification_sites.tsv` | 已知修饰位点 | TSV: sequence_id, position, modification, base |
| `{aa}/ref.fa` | 各氨基酸参考序列 | FASTA（从 unmodified.fa 按氨基酸拆分） |

## 依赖

- **seqkit** — FASTA 操作
- **MAFFT (linsi)** — 多序列比对
- **Python 3** + `biopython`, `tqdm`, `numpy`

```bash
conda install -c bioconda seqkit mafft
pip install biopython tqdm numpy
```

## 用法

### 1. 配置

编辑 `config.sh`，设置输入文件路径：

```bash
FEATURES_TSV="path/to/all.RNA_features.tsv"
REF_FA="unmodified.fa"
MOD_SITES_TSV="modification_sites.tsv"
```

### 2. 准备参考目录

将 `unmodified.fa` 按氨基酸类型拆分到各目录的 `ref.fa`：

```bash
while read aa; do
    mkdir -p $aa
    # 根据 sequence ID 或外部注释拆分
    seqkit grep -p "$aa" unmodified.fa > $aa/ref.fa
done < aalist
```

### 3. 运行 Pipeline

```bash
# 完整流程
bash pipeline.sh

# 只处理某个氨基酸
bash pipeline.sh --aa Gly

# 从 Step 2 开始（Step 1 已完成）
bash pipeline.sh --from 2

# 只运行 Step 3（比对已完成）
bash pipeline.sh --only 3
```

### 4. 输出

每个 tRNA 比对文件生成一个 `results/tRNA_xxx_site_stats.tsv`：

| 列 | 说明 |
|----|------|
| query_id | 查询序列 ID |
| query_pos | 查询序列上的位置 |
| query_base | 碱基 |
| aln_column | 比对列号 |
| mod_type | 修饰类型（如 pGm, pD, xU） |
| match | 匹配该修饰的参考序列数 |
| mismatch | 不匹配的参考序列数 |
| valid | 有效比对的参考序列数 |
| support_rate | 支持率 = match / valid |

## 目录结构

```
├── config.sh                 # 配置文件
├── aalist                    # 23 个氨基酸列表
├── pipeline.sh               # 主入口
├── 01_extract_tRNA.sh        # Step 1: 提取序列
├── 02_align.sh               # Step 2: 并行比对
├── 02_align_single.sh        # Step 2: 单条比对
├── 03_map_modifications.sh   # Step 3: 修饰映射
├── map_modifications.py      # 核心 Python 脚本
└── example/                  # 示例数据
    ├── all.RNA_features.example.tsv
    ├── modification_sites.example.tsv
    ├── unmodified.example.fa
    └── ref.example.fa
```

## 示例数据

`example/` 目录下提供小规模示例文件，可用于验证流程格式正确性。