#!/usr/bin/env python3
"""
大规模处理版：优化内存和性能，适合20万+文件处理
用法: python map_site_mods_bulk.py <修饰文件> <FASTA目录> [输出目录] [进程数]
"""

import os
import sys
import glob
import gzip
import pickle
import time
import multiprocessing as mp
from collections import defaultdict
from pathlib import Path
from tqdm import tqdm
from Bio import SeqIO
import numpy as np

class ModCache:
    """修饰信息缓存，避免重复加载"""
    _instance = None
    _cache_file = ".mod_cache.pkl"
    
    def __new__(cls, mod_file):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance.load(mod_file)
        return cls._instance
    
    def load(self, mod_file):
        """加载或从缓存读取修饰信息"""
        cache_path = Path(self._cache_file)
        
        # 检查缓存是否有效（修饰文件比缓存新）
        if cache_path.exists():
            cache_mtime = cache_path.stat().st_mtime
            mod_mtime = Path(mod_file).stat().st_mtime
            
            if cache_mtime > mod_mtime:
                print(f"从缓存加载修饰信息...")
                with open(cache_path, 'rb') as f:
                    self.mod_dict, self.ref_mods = pickle.load(f)
                print(f"  缓存加载完成")
                return
        
        print("加载修饰信息...")
        mod_dict = defaultdict(list)
        ref_mods = defaultdict(dict)
        
        # 使用numpy加速读取大文件
        with open(mod_file, 'r') as f:
            header = f.readline()
            # 预分配内存
            lines = f.readlines()
            
            for line in tqdm(lines, desc="解析修饰数据"):
                if not line.strip():
                    continue
                parts = line.strip().split()
                if len(parts) >= 4:
                    ref_id = parts[0]
                    position = int(parts[1])
                    mod_type = parts[2]
                    base = parts[3]
                    
                    key = (ref_id, position, base)
                    mod_dict[key].append(mod_type)
                    
                    if position not in ref_mods[ref_id]:
                        ref_mods[ref_id][position] = {'base': base, 'mods': []}
                    ref_mods[ref_id][position]['mods'].append(mod_type)
        
        self.mod_dict = dict(mod_dict)
        self.ref_mods = dict(ref_mods)
        
        # 保存缓存
        with open(cache_path, 'wb') as f:
            pickle.dump((self.mod_dict, self.ref_mods), f)
        
        total_mods = len(self.mod_dict)
        multi_mod_count = sum(1 for mods in self.mod_dict.values() if len(mods) > 1)
        print(f"  总修饰记录: {total_mods:,}")
        if multi_mod_count > 0:
            print(f"  其中 {multi_mod_count:,} 个位置有多个修饰类型")

def process_single_file(args):
    """处理单个文件的工作函数"""
    fasta_file, output_dir, mod_dict, ref_mods = args
    
    try:
        fa_path = Path(fasta_file)
        base_name = fa_path.stem
        stats_file = Path(output_dir) / f"{base_name}_site_stats.tsv"
        
        # 如果输出文件已存在且比输入文件新，跳过
        if stats_file.exists():
            if stats_file.stat().st_mtime > fa_path.stat().st_mtime:
                return f"跳过（已处理）: {fa_path.name}", True
        
        # 读取序列（使用生成器节省内存）
        records = list(SeqIO.parse(fasta_file, "fasta"))
        if len(records) < 2:
            return f"跳过（序列不足）: {fa_path.name}", True
        
        query_record = records[0]
        query_id = query_record.id
        query_seq = str(query_record.seq).upper()
        ref_records = records[1:]
        
        # 快速计算查询序列长度
        query_length = query_seq.count('A') + query_seq.count('C') + query_seq.count('G') + query_seq.count('T') + query_seq.count('U')
        
        # 预计算查询序列的非gap位置
        query_positions = []
        query_chars = 0
        for i, char in enumerate(query_seq):
            if char != '-':
                query_chars += 1
                query_positions.append((query_chars, i, char))  # (原始位置, 比对位置, 碱基)
        
        # 收集所有结果
        results = []
        
        for query_pos, query_aln_pos, query_base in query_positions:
            # 跳过gap
            if query_base == '-':
                continue
            
            # 统计变量
            total_refs = len(ref_records)
            valid_count = 0
            invalid_count = 0
            gap_ref_count = 0
            
            # 按修饰类型计数
            mod_type_counts = defaultdict(int)
            
            for ref_record in ref_records:
                ref_id = ref_record.id
                ref_seq = str(ref_record.seq).upper()
                
                if query_aln_pos >= len(ref_seq):
                    ref_base = '-'
                else:
                    ref_base = ref_seq[query_aln_pos]
                
                if ref_base == '-':
                    gap_ref_count += 1
                elif query_base == ref_base:
                    # 有效比对，计算参考位置
                    ref_chars = 0
                    ref_pos = None
                    for i, char in enumerate(ref_seq):
                        if char != '-':
                            ref_chars += 1
                            if i == query_aln_pos:
                                ref_pos = ref_chars
                                break
                    
                    if ref_pos is not None:
                        valid_count += 1
                        
                        # 检查修饰
                        if ref_id in ref_mods and ref_pos in ref_mods[ref_id]:
                            ref_data = ref_mods[ref_id][ref_pos]
                            if ref_data['base'] == ref_base:
                                for mod_type in ref_data['mods']:
                                    mod_type_counts[mod_type] += 1
                else:
                    invalid_count += 1
            
            # 生成统计行
            for mod_type, match_count in mod_type_counts.items():
                mismatch_count = valid_count - match_count
                support_rate = match_count / valid_count if valid_count > 0 else 0
                
                results.append('\t'.join([
                    query_id,
                    str(query_pos),
                    query_base,
                    str(query_aln_pos + 1),
                    mod_type,
                    str(match_count),
                    str(mismatch_count),
                    str(valid_count),
                    str(invalid_count),
                    str(gap_ref_count),
                    str(total_refs),
                    f"{support_rate:.2f}"
                ]))
        
        # 写入文件
        if results:
            with open(stats_file, 'w') as f:
                f.write('query_id\tquery_pos\tquery_base\taln_column\tmod_type\tmatch\tmismatch\tvalid\tinvalid\tgap_ref\ttotal_refs\tsupport_rate\n')
                f.write('\n'.join(results))
            
            # 压缩输出文件以节省空间
            # with gzip.open(f"{stats_file}.gz", 'wt') as f:
            #     f.write('\n'.join([header] + results))
            # stats_file.unlink()  # 删除未压缩文件
        
        return f"完成: {fa_path.name} ({len(results)}个位点)", True
        
    except Exception as e:
        return f"错误: {fa_path.name} - {str(e)}", False

def process_files_parallel(fasta_files, output_dir, mod_cache, num_processes):
    """并行处理文件"""
    print(f"\n开始并行处理 {len(fasta_files):,} 个文件...")
    print(f"使用 {num_processes} 个进程")
    
    # 准备参数
    args_list = [(f, output_dir, mod_cache.mod_dict, mod_cache.ref_mods) 
                 for f in fasta_files]
    
    # 进度条
    completed = 0
    failed = 0
    
    with mp.Pool(processes=num_processes) as pool:
        with tqdm(total=len(fasta_files), desc="处理进度") as pbar:
            for result, success in pool.imap_unordered(process_single_file, args_list):
                completed += 1
                if not success:
                    failed += 1
                
                # 每1000个文件或最后显示一次状态
                if completed % 1000 == 0 or completed == len(fasta_files):
                    pbar.update(completed - pbar.n)
                    pbar.set_postfix({
                        '完成': f"{completed:,}",
                        '失败': f"{failed:,}",
                        '成功率': f"{(completed-failed)/completed*100:.1f}%"
                    })
                
                # 记录失败日志
                if not success and failed <= 100:  # 只记录前100个错误
                    with open(Path(output_dir) / "error_log.txt", 'a') as f:
                        f.write(f"{result}\n")
    
    return completed, failed

def main():
    if len(sys.argv) < 3:
        print("大规模tRNA修饰统计工具")
        print("用法: python map_site_mods_bulk.py <修饰文件> <FASTA目录> [输出目录] [进程数]")
        print("示例: python map_site_mods_bulk.py modification_sites.tsv ./afa_files/ ./results/ 32")
        sys.exit(1)
    
    mod_file = sys.argv[1]
    fasta_dir = sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "./site_mod_stats"
    num_processes = int(sys.argv[4]) if len(sys.argv) > 4 else mp.cpu_count()
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 清空错误日志
    error_log = Path(output_dir) / "error_log.txt"
    if error_log.exists():
        error_log.unlink()
    
    # 加载修饰信息（单例缓存）
    start_time = time.time()
    mod_cache = ModCache(mod_file)
    
    # 查找文件（分批处理，避免内存爆炸）
    print(f"\n扫描目录: {fasta_dir}")
    fasta_files = []
    
    # 分批查找，避免一次加载太多文件路径
    for batch in tqdm(range(100), desc="扫描文件批次"):
        pattern = os.path.join(fasta_dir, f"tRNA_*_sp{batch:03d}*.afa")
        batch_files = glob.glob(pattern)
        fasta_files.extend(batch_files)
    
    # 如果没有找到特定模式，尝试通用查找
    if not fasta_files:
        print("未找到tRNA_*_sp*.afa模式，尝试通用查找...")
        fasta_files = list(Path(fasta_dir).glob("tRNA*.afa"))
        fasta_files = [str(f) for f in fasta_files]
    
    print(f"找到 {len(fasta_files):,} 个tRNA比对文件")
    
    if not fasta_files:
        print("错误: 没有找到任何tRNA*.afa文件")
        sys.exit(1)
    
    # 分批处理，每批5000个文件
    batch_size = 5000
    total_completed = 0
    total_failed = 0
    
    for i in range(0, len(fasta_files), batch_size):
        batch = fasta_files[i:i + batch_size]
        print(f"\n处理批次 {i//batch_size + 1}/{(len(fasta_files)+batch_size-1)//batch_size}")
        print(f"本批文件: {len(batch):,} 个")
        
        completed, failed = process_files_parallel(
            batch, output_dir, mod_cache, min(num_processes, len(batch))
        )
        
        total_completed += completed
        total_failed += failed
        
        # 每批完成后保存进度
        with open(Path(output_dir) / "progress.txt", 'w') as f:
            f.write(f"已处理: {total_completed:,}/{len(fasta_files):,}\n")
            f.write(f"失败: {total_failed:,}\n")
            f.write(f"成功率: {(total_completed-total_failed)/total_completed*100:.1f}%\n")
    
    # 最终统计
    total_time = time.time() - start_time
    avg_time_per_file = total_time / total_completed if total_completed > 0 else 0
    
    print(f"\n{'='*50}")
    print(f"处理完成！")
    print(f"总文件数: {len(fasta_files):,}")
    print(f"成功处理: {total_completed - total_failed:,}")
    print(f"失败: {total_failed:,}")
    print(f"成功率: {(total_completed-total_failed)/total_completed*100:.1f}%")
    print(f"总耗时: {total_time:.1f}秒 ({total_time/3600:.2f}小时)")
    print(f"平均每个文件: {avg_time_per_file:.3f}秒")
    print(f"结果目录: {output_dir}")
    
    if total_failed > 0:
        print(f"错误日志: {error_log}")
        print(f"（只显示前100个错误，更多错误请查看日志文件）")

if __name__ == "__main__":
    # 设置递归限制，避免大文件出错
    sys.setrecursionlimit(1000000)
    
    # 设置进程启动方法（Linux下推荐'forkserver'）
    if sys.platform != "win32":
        mp.set_start_method('forkserver', force=True)
    
    main()
