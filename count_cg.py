import argparse

# 创建命令行参数解析器
parser = argparse.ArgumentParser(description='Count the number of C and G in the FASTA format genome file')
parser.add_argument('filename', metavar='FILENAME', type=str,
                    help='genome file path')

# 解析命令行参数
args = parser.parse_args()

# 打开基因组文件
with open(args.filename, 'r') as file:
    # 跳过文件头
    file.readline()
    # 初始化计数器
    count_c = 0
    count_g = 0
    # 遍历文件中的每一行
    for line in file:
        # 忽略注释行
        if line.startswith('>'):
            continue
        # 统计C和G
        count_c += line.upper().count('C')
        count_g += line.upper().count('G')
    count_cg = count_c + count_g

# 打印结果
print(f'{count_c} C, {count_g} G and {count_cg} CG in genome')

