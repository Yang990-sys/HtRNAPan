import re
import argparse
from collections import defaultdict

def read_sequences(file_path):
    """从文件读取序列数据，返回字典"""
    sequences = {}
    with open(file_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    current_id = None
    current_seq = []
    
    for line in lines:
        if line.startswith('>'):
            if current_id:
                sequences[current_id] = ''.join(current_seq)
            current_id = line[1:]
            current_seq = []
        else:
            current_seq.append(re.sub(r'\s+', '', line))
    
    if current_id:
        sequences[current_id] = ''.join(current_seq)
    
    return sequences

def read_modifications(file_path):
    """从文件读取修饰数据，格式：tRNA_id 位置 修饰类型"""
    modifications = []
    with open(file_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):  # 跳过注释和空行
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            try:
                tRNA_id, pos_str, mod_type = parts
                pos = int(pos_str)
                modifications.append((tRNA_id, pos, mod_type))
            except ValueError:
                print(f"警告: 修饰文件第{line_num}行位置不是整数，已跳过")
    return modifications

def get_all_aligned_positions(sequences):
    """获取所有序列的对齐位置及其对应的实际位置和碱基"""
    aligned_data = {}
    for seq_id, seq in sequences.items():
        positions = []  # 存储实际位置
        bases = []      # 存储对应位置的碱基
        aligned_pos = 0
        for char in seq:
            if char != '-':  # 忽略gap
                aligned_pos += 1
            positions.append(aligned_pos if char != '-' else None)
            bases.append(char if char != '-' else None)
        aligned_data[seq_id] = {'positions': positions, 'bases': bases}
    return aligned_data

def get_valid_positions(aligned_data, exclude_target='target'):
    """获取所有序列中存在碱基（非gap）的对齐位置"""
    valid_positions = set()
    for seq_id, data in aligned_data.items():
        if seq_id == exclude_target:  # 跳过目标序列
            continue
        for align_index, pos in enumerate(data['positions']):
            if pos is not None:  # 只考虑有实际位置的对齐点
                valid_positions.add(align_index)
    return sorted(valid_positions)

def predict_modifications(sequences, modifications, target_id='target', min_samples=2, base_match_threshold=0.5):
    """
    预测目标序列的修饰比例，包含碱基匹配校验
    
    参数:
        sequences: 所有序列的字典
        modifications: 已知修饰的列表，每个元素格式为 (tRNA_id, 位置, 修饰类型)
        target_id: 目标序列的ID
        min_samples: 最小样本量
        base_match_threshold: 碱基匹配阈值，低于此值的位置不进行修饰预测
    
    返回:
        预测的目标序列修饰比例字典
    """
    if target_id not in sequences:
        raise ValueError(f"目标序列 {target_id} 不在输入序列中")
    
    # 获取所有序列的对齐数据
    aligned_data = get_all_aligned_positions(sequences)
    
    # 收集所有tRNA ID（排除目标序列）
    tRNA_ids = [seq_id for seq_id in sequences.keys() if seq_id != target_id]
    
    # 获取所有有效对齐位置
    valid_align_positions = get_valid_positions(aligned_data, target_id)
    
    # 构建修饰字典：{对齐位置: {tRNA_id: 修饰类型}}
    mod_dict = defaultdict(dict)
    for tRNA_id, pos, mod_type in modifications:
        if tRNA_id not in sequences or tRNA_id == target_id:
            continue
            
        # 找到该位置对应的对齐索引
        positions = aligned_data[tRNA_id]['positions']
        align_index = None
        for i, p in enumerate(positions):
            if p == pos:
                align_index = i
                break
                
        if align_index is not None:
            mod_dict[align_index][tRNA_id] = mod_type
    
    # 计算每个对齐位置的修饰比例
    target_predictions = {}
    target_data = aligned_data[target_id]
    target_positions = target_data['positions']
    target_bases = target_data['bases']
    target_seq = sequences[target_id]
    
    for align_index in valid_align_positions:
        # 检查目标序列在该位置是否有碱基
        if align_index >= len(target_seq) or target_seq[align_index] == '-':
            continue
            
        target_pos = target_positions[align_index]
        target_base = target_bases[align_index]
        if target_pos is None or target_base is None:
            continue
        
        # 先检查该位置的碱基匹配情况
        base_matches = 0
        total_bases = 0
        base_counts = defaultdict(int)
        
        for tRNA_id in tRNA_ids:
            tRNA_data = aligned_data[tRNA_id]
            if align_index >= len(tRNA_data['positions']):
                continue
                
            tRNA_pos = tRNA_data['positions'][align_index]
            tRNA_base = tRNA_data['bases'][align_index]
            
            if tRNA_pos is None or tRNA_base is None:
                continue  # 该tRNA在这个对齐位置是gap，不纳入统计
            
            total_bases += 1
            base_counts[tRNA_base] += 1
            if tRNA_base == target_base:
                base_matches += 1
        
        # 计算碱基匹配率
        if total_bases == 0:
            continue
            
        base_match_rate = base_matches /total_bases
        if base_match_rate < base_match_threshold:
            print(f"信息: 位置 {target_pos} 碱基匹配率低 ({base_match_rate:.2%})，已跳过修饰预测")
            continue
        
        # 统计该位置的所有tRNA（包括未修饰的），但只考虑碱基匹配的
        total = 0
        mod_counts = defaultdict(int)
        
        for tRNA_id in tRNA_ids:
            tRNA_data = aligned_data[tRNA_id]
            if align_index >= len(tRNA_data['positions']):
                continue
                
            tRNA_pos = tRNA_data['positions'][align_index]
            tRNA_base = tRNA_data['bases'][align_index]
            
            if tRNA_pos is None or tRNA_base is None:
                continue  # 该tRNA在这个对齐位置是gap，不纳入统计
            
            # 只考虑碱基匹配的tRNA
            if tRNA_base != target_base:
                continue
                
            total += 1
            # 检查是否有修饰
            if tRNA_id in mod_dict.get(align_index, {}):
                mod_type = mod_dict[align_index][tRNA_id]
                mod_counts[mod_type] += 1
            else:
                mod_counts['未修饰'] += 1  # 记录未修饰的情况
        
        # 过滤样本量不足的位置
        if total < min_samples:
            print(f"信息: 位置 {target_pos} 碱基匹配的样本量不足 ({total} < {min_samples})，已跳过")
            continue
        
        # 过滤全部未修饰的位置
        if len(mod_counts) == 1 and '未修饰' in mod_counts:
            continue  # 全未修饰，不输出
        
        # 计算比例
        mods_with_ratio = {}
        for mod_type, count in mod_counts.items():
            if mod_type!="未修饰":
                mods_with_ratio[mod_type] = count / total
 
        # 添加碱基信息
        target_predictions[target_pos] = {
            'modifications': mods_with_ratio,
            'target_base': target_base,
            'base_match_rate': base_match_rate,
            'base_counts': dict(base_counts)
        }
    
    return target_predictions

def write_output(predictions, output_file):
    """将预测结果写入输出文件"""
    with open(output_file, 'w') as f:
        f.write("Pos\tBase\tBaseMatch_Rate\tModification:Rate\n")
        for pos in sorted(predictions.keys()):
            data = predictions[pos]
            mod_strings = [f"{mod_type}:{ratio:.2%}" for mod_type, ratio in data['modifications'].items()]
            f.write(f"{pos}\t{data['target_base']}\t{data['base_match_rate']:.2%}\t{', '.join(mod_strings)}\n")

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='根据已知tRNA修饰预测目标序列的修饰比例（带碱基校验）')
    parser.add_argument('-s', '--sequences', required=True, help='包含目标序列和参考tRNA序列的FASTA文件')
    parser.add_argument('-m', '--modifications', required=True, help='包含修饰信息的文件（格式：tRNA_id 位置 修饰类型）')
    parser.add_argument('-o', '--output', required=True, help='输出结果文件路径')
    parser.add_argument('-t', '--target', default='target', help='目标序列的ID（默认：target）')
    parser.add_argument('-min', '--min_samples', type=int, default=2, help='最小样本量（默认：2）')
    parser.add_argument('-minP', '--base_threshold', type=float, default=0.1, 
                      help='碱基匹配阈值，低于此值不进行预测（默认：0.5）')
    
    args = parser.parse_args()
    
    # 读取输入数据
    print("读取序列数据...")
    sequences = read_sequences(args.sequences)
    
    print("读取修饰数据...")
    modifications = read_modifications(args.modifications)
    
    # 预测修饰
    print("计算修饰比例（带碱基校验）...")
    predictions = predict_modifications(
        sequences, 
        modifications, 
        target_id=args.target, 
        min_samples=args.min_samples,
        base_match_threshold=args.base_threshold
    )
    
    # 写入输出
    print(f"写入结果到 {args.output}...")
    write_output(predictions, args.output)
    
    print("完成!")

if __name__ == "__main__":
    main()
    
