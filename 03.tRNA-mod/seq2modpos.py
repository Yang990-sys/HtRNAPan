import pandas as pd
import sys

def read_modification_table(file_path):
    """读取修饰对照表文件（TSV格式）"""
    try:
        df = pd.read_csv(file_path, sep='\t')
        return dict(zip(df['mod'], df['base']))
    except Exception as e:
        print(f"读取修饰对照表失败: {e}")
        sys.exit(1)

def find_modifications(modified_seq, unmodified_seq, modification_table, seq_id):
    """查找序列中的修饰位点，增强错误提示"""
    modifications = []
    mod_ptr = 0
    unmod_ptr = 0
    
    while mod_ptr < len(modified_seq):
        found_mod = False
        max_match_len = 0
        best_mod = None
        
        # 找出所有可能的匹配，选择最长的
        for mod in modification_table:
            if modified_seq[mod_ptr:].startswith(mod):
                if len(mod) > max_match_len:
                    max_match_len = len(mod)
                    best_mod = mod
        
        if best_mod:
            # 记录修饰
            position = unmod_ptr + 1
            base = modification_table[best_mod]
            modifications.append((position, best_mod, base))
            
            # 更新指针
            mod_ptr += max_match_len
            unmod_ptr += 1
            found_mod = True
        else:
            # 处理普通碱基
            if mod_ptr < len(modified_seq) and unmod_ptr < len(unmodified_seq):
                if modified_seq[mod_ptr] == unmodified_seq[unmod_ptr]:
                    mod_ptr += 1
                    unmod_ptr += 1
                else:
                    # 输出上下文信息
                    context_start = max(0, unmod_ptr - 5)
                    context_end = min(len(modified_seq), unmod_ptr + 5)
                    
                    mod_context = modified_seq[context_start:context_end]
                    unmod_context = unmodified_seq[context_start:context_end]
                    context_pos = unmod_ptr - context_start + 1
                    
                    print(f"\n序列比对错误: 序列 '{seq_id}' 位置 {unmod_ptr+1} 处不匹配")
                    print(f"修饰序列片段:   {mod_context}")
                    print(f"未修饰序列片段: {unmod_context}")
                    print(f"{' ' * (context_pos * 2 + 17)}↑ 错误位置")
                    
                    mod_ptr += 1
                    unmod_ptr += 1
            else:
                break  # 序列结束
    
    return modifications

def read_fasta(file_path):
    """手动读取FASTA文件，保留所有字符"""
    sequences = {}
    current_id = None
    current_seq = []
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue  # 跳过空行
                if line.startswith('>'):
                    # 新序列的ID
                    if current_id is not None:
                        sequences[current_id] = ''.join(current_seq)
                    current_id = line[1:]  # 去掉>符号
                    current_seq = []
                else:
                    # 序列内容
                    current_seq.append(line)
        
        # 添加最后一个序列
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
        
        return sequences
    
    except Exception as e:
        print(f"读取FASTA文件失败: {e}")
        sys.exit(1)

def main():
    # 参数设置
    MOD_TABLE_FILE = 'mod_base.tsv'
    MODIFIED_FASTA = 'modified.fa'
    UNMODIFIED_FASTA = 'unmodified.fa'
    OUTPUT_FILE = 'modification_sites.tsv'
    
    # 读取修饰对照表
    mod_table = read_modification_table(MOD_TABLE_FILE)
    
    # 读取FASTA文件（不使用Bio库）
    modified_records = read_fasta(MODIFIED_FASTA)
    unmodified_records = read_fasta(UNMODIFIED_FASTA)
    
    # 检查序列ID是否匹配
    common_ids = set(modified_records.keys()) & set(unmodified_records.keys())
    if not common_ids:
        print("错误: 修饰和未修饰序列文件中没有共同的序列标识")
        sys.exit(1)
    
    # 处理每个序列
    all_results = []
    for seq_id in common_ids:
        modified_seq = modified_records[seq_id]
        unmodified_seq = unmodified_records[seq_id]
        
        
        # 查找修饰位点
        modifications = find_modifications(modified_seq, unmodified_seq, mod_table, seq_id)
        
        # 添加到结果列表
        for pos, mod, base in modifications:
            all_results.append([seq_id, pos, mod, base])
    
    # 输出到TSV文件
    if all_results:
        result_df = pd.DataFrame(all_results, columns=['sequence_id', 'position', 'modification', 'base'])
        result_df.to_csv(OUTPUT_FILE, sep='\t', na_rep='nan', index=False)
        print(f"成功输出结果到 {OUTPUT_FILE}")
        print(f"共发现 {len(all_results)} 个修饰位点")
    else:
        print("未发现任何修饰位点")

if __name__ == "__main__":
    main()
