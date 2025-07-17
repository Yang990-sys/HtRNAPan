from bs4 import BeautifulSoup
import re
import argparse
import json
import os

parser = argparse.ArgumentParser(description='从 HTML 文件中提取分子信息')
parser.add_argument('-i', '--input', required=True, help='输入 HTML 文件路径')
parser.add_argument('-o', '--output', required=True, help='输出结果文件路径')
args = parser.parse_args()

# 文件名
file_name = args.input

try:
    # 读取 HTML 文件内容
    with open(file_name, 'r', encoding='utf-8') as file:
        html_content = file.read()
except FileNotFoundError:
    print(f"未找到文件: {file_path}")

class CustomEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, datetime):
            return obj.isoformat()  # 将日期转换为字符串
        return super().default(obj)

# 解析 HTML 内容
soup = BeautifulSoup(html_content, 'html.parser')

# 提取 ID Card: 表格
summary_table = soup.find('h3', string='ID Card:').find_next('table')
# 提取Publications:
publications=soup.find('h3', string='Publications:').find_next('table')

# 定义要提取的信息
keys = [
    'ID',
    'Sequence',
    'Full name:',
    'UniProt:',
    'Structures:',
    'Enzyme type:',
    'Position of modification - modification:',
    'Publications:'
]

extracted_info={}

for key in keys:
        extracted_info[key] = ''

# 提取Protein sequence:
sequence=soup.find('h3', string='Protein sequence:').find_next('pre')
if sequence:
    extracted_info['Sequence']=soup.find('h3', string='Protein sequence:').find_next('pre').text
extracted_info['ID'] = os.path.basename(file_name).split(".")[0]



# 遍历 Summary 表格提取信息
for row in summary_table.find_all('tr'):
    key = row.find('th').text.strip()
    value =row.find('td').text.strip()
    for target_key in keys:
        if key == target_key:
            extracted_info[f'{target_key}'] = value


#遍历 Publications: 表格提取信息
headers = [th.text.strip() for th in publications.find("tr").find_all("th")]

data_all=[]
for tr in publications.find_all("tr")[1:]:  # 跳过表头行
    data=[]
    td_tags = tr.find_all("td")  # 获取该行所有单元格
    for i, td in enumerate(td_tags, 1):  # 按列索引处理（假设列名自定义）
        a_tag = td.find("a")  # 每个单元格可能只有一个链接
        if a_tag:
            result = {}
            result[a_tag.text.strip()]=a_tag.get("href")
            data.append(result)
        else:
            data.append(td.text.strip())
    data_all.append(dict(zip(headers, data)))  # 合并表头和数据为字典

extracted_info['Publications:']=data_all
# 提取结构式图片链接


# 输出文件名
print(f"文件名: {file_name}")

# 输出提取的信息
#for key, value in extracted_info.items():
#    print(f"{key} {value}")

with open(args.output, "w") as f:
    json.dump(extracted_info, f, cls=CustomEncoder, indent=2)
