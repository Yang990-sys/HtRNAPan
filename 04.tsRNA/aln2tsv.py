import argparse

def main():
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='处理tsRNA的aln文件，整理格式并分类')
    parser.add_argument('-i', '--input', required=True, help='输入的aln文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出的结果文件路径')
    
    args = parser.parse_args()
    
    with open(args.input, "r") as infile, open(args.output, "w") as outfile:
        current_tRNA = None
        for line in infile:
            line = line.strip()
            if not line:
                continue
            # 识别tRNA标识行
            if line.startswith(">"):
                current_tRNA = line[1:].strip()
            # 处理tsRNA行
            elif current_tRNA:
                parts = line.split()
                if len(parts) >= 2:
                    # 去除序列中的点
                    tsrna_seq = parts[0].replace('.', '')
                    # 将tsrna_info中的空格替换为制表符
                    tsrna_info = '\t'.join(parts[1:])
                    info_parts = parts[1:]
                    if len(info_parts) >= 5:
                        try:
                            ninth_value = int(info_parts[4])
                            if ninth_value >= 35:
                                category = "5tRF"
                            else:
                                category = "3tRF"
                        except ValueError:
                            category = "unknown"  # 无法转换为数字时的默认值
                    else:
                        category = "unknown"  # 列数不足时的默认值
                    
                    # 组合所有字段并写入，字段间用制表符分隔
                    outfile.write(f"{current_tRNA}\t{tsrna_seq}\t{tsrna_info}\t{category}\n")

if __name__ == "__main__":
    main()

