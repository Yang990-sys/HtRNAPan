from bs4 import BeautifulSoup
import re
import argparse
import json

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

# 提取 Summary 表格
summary_table = soup.find('h3', string='Summary').find_next('table')
# 提取 Chemical information 表格
chemical_info_table = soup.find('h3', string='Chemical information').find_next('table')
# 提取Download Structures
DS_table=soup.find('h3', string='Download Structures').find_next('table')
# 提取LC-MS Information
LCMS_table=soup.find('h3', string='LC-MS Information').find_next('table')
# 提取
publications=soup.find('h3', string='LC-MS Publications')
# 定义要提取的信息
keys = [
    ('Summary', 'Full name'),
    ('Summary', 'Short name'),
    ('Summary', 'MODOMICS code'),
    ('Summary', 'IUPAC name'),
    ('Summary', 'Synonyms'),
    ('Summary', 'Enzymes'),
    ('Chemical information', 'Sum formula'),
    ('Chemical information', 'SMILES'),
    ('Chemical information', 'Search the molecule in external databases'),
    ('Download Structures','2D'),
    ('Download Structures','3D'),
    ('LC-MS Information','Monoisotopic mass'),
    ('LC-MS Information','Average mass'),
    ('LC-MS Information','[M+H]+'),
    ('LC-MS Information','Product ions') 
    ]

extracted_info = {}

# 遍历 Summary 表格提取信息
for row in summary_table.find_all('tr'):
    columns = row.find_all('td')
    if len(columns) == 2:
        key = columns[0].text.strip().replace('</b>', ',').replace('<b>', ',')
        value = columns[1].text.replace('\n', ' ')
        value=re.sub(r' +', ',', value)
        value=re.sub(r'^,', '', value)
        for section, target_key in keys:
            if section == 'Summary' and key == target_key:
                extracted_info[f'{section} - {target_key}'] = value

# 遍历 Chemical information 表格提取信息
for row in chemical_info_table.find_all('tr'):
    columns = row.find_all('td')
    if len(columns) == 2:
        key = columns[0].text.strip().replace('</b>', '').replace('<b>', '')
        value = columns[1].text.strip()
        for section, target_key in keys:
            if section == 'Chemical information' and key == target_key:
                extracted_info[f'{section} - {target_key}'] = value
            if section == 'Chemical information' and target_key == 'Search the molecule in external databases':
                extracted_info[f'{section} - {target_key}'] = ",".join([a['href'].strip() for a in  columns[1].find_all('a')])
# 遍历 Structure 表格提取信息
for row in DS_table.find_all('tr'):
    columns = row.find_all('td')
    if len(columns) == 2:
        key = columns[0].text.strip().replace('</b>', '').replace('<b>', '')
        value = ",".join(["https://iimcb.genesilico.pl/"+a['href'] for a in  columns[1].find_all('a')])
        for section, target_key in keys:
            if section == 'Download Structures' and key == target_key:
                extracted_info[f'{section} - {target_key}'] = value

# 遍历 LCMS 表格提取信息
for row in LCMS_table.find_all('tr'):
    columns = row.find_all('td')
    if len(columns) == 2:
        key = columns[0].text.strip().replace('</b>', '').replace('<b>', '')
        value = columns[1].text.strip()
        for section, target_key in keys:
            if section == 'LC-MS Information' and key == target_key:
                extracted_info[f'{section} - {target_key}'] = value

# 遍历 LCMS Publications 表格提取信息
if publications :
    publications=publications.find_next('table')
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
else:
    extracted_info['Publications:']=[]

# 提取结构式图片链接
image = soup.find('img', id='image')
image_url = "https://iimcb.genesilico.pl//modomics/"+image['src'] if image else ''
extracted_info['Structural formula image URL'] = image_url

# 确保所有需要的键都存在，如果不存在则设为空字符串
for section, key in keys:
    full_key = f'{section} - {key}'
    if full_key not in extracted_info:
        extracted_info[full_key] = ''

# 输出文件名
print(f"文件名: {file_name}")

# 输出提取的信息
#for key, value in extracted_info.items():
#    print(f"{key}: {value}")

with open(args.output, "w") as f:
    json.dump(extracted_info, f, cls=CustomEncoder, indent=2)
