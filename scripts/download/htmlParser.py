from bs4 import BeautifulSoup
import csv
import glob
def getPuropose(infile:str,purpose:str):
    with open(infile, encoding="utf-8") as f:
        html = f.read()

    soup = BeautifulSoup(html, "lxml")

    table = soup.find_all(purpose)
    return table

def table2csv(table,outfile:str):
    data = []
    for row in table.find_all("tr"):
        cells = row.find_all(["td", "th"])
        values = []
        for cell in cells:
            a_tag = cell.find("a")
            if a_tag:
                text = a_tag.get_text(strip=True)
                href = a_tag.get("href", "")
                if text:
                    values.append(text)  
                elif href:
                    values.append(href)  
                else:
                    values.append("")
            else:
                # 没有 <a>，提取纯文本
                values.append(cell.get_text(strip=True))
        if values:
            data.append(values)

    for row in data:
        print(row)
    with open(outfile, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerows(data)

def getTrInofromation(
        html_path:str,
        text:str
):
    """
    提取GSM html页面信息,为"Characteristics"设计
    """
    with open(html_path, encoding="utf-8") as f:
        html = f.read()
    soup = BeautifulSoup(html, "html.parser")
    result = {}
    for tr in soup.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) < 2:
            continue
        key = tds[0].get_text(strip=True)
        value_td = tds[1]
        if key == text:
            information = value_td.get_text(separator = "\n",strip=True)
            for line in information.split("\n"):
                if ":" in line:
                    sub_key, sub_value = line.split(":", 1)
                    result[sub_key.strip()] = sub_value.strip()
    return result
    



if __name__ == '__main__':
    # infile = "/home/luosg/Data/genomeStability/data/GSM/GSM2786643.html"
    infiles = glob.glob("/disk5/luosg/GCN2_20251224/data/GCN2pub/html/*.html")
    outfile = "/disk5/luosg/GCN2_20251224/data/GCN2pub/Characteristics.csv"

    all_samples = {} 

    for infile in infiles:
        print(f"Processing {infile}")
        GSM = infile.split("/")[-1].replace(".html","")
        dictInformation = getTrInofromation(infile,"Characteristics")
        all_samples[GSM] = dictInformation
    all_keys = set()
    for info in all_samples.values():
        if info is None:
            continue
        all_keys |= info.keys()
    all_keys = sorted(all_keys)
    with open(outfile, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["GSM"] + all_keys)
        writer.writeheader()

        for gsm, info in all_samples.items():
            row = {"GSM": gsm}
            if info is None:
                writer.writerow(row)
                continue
            row.update(info)  
            writer.writerow(row)

        





