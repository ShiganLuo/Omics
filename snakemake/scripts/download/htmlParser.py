from bs4 import BeautifulSoup
import csv

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

if __name__ == '__main__':
    infile = "/ChIP_seq_2/StemCells/ChIPseq_Sox2/data/sample.html"
    tables = getPuropose(infile,'table')

    table2csv(tables[2],"/ChIP_seq_2/StemCells/ChIPseq_Sox2/data/sample.csv")



