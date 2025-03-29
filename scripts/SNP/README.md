# Description

for TE: temporary only in favor of TEcount,not TElocal

annovar.sh: This script is originally used to perform annovar analysis
vcfIntersectBed.sh: This script is originally used to perform intersect between vcf.gz and bed

20250327:

    - commonExpression.py 增加命令行参数支持，以整合annovar分析流程进snakemake

20250318:

    - getBed.py 增加TE和基因bed文件分开写入功能，原来的命令无需修改也能正确执行。目的：适应StringTie流程
    - getBed.py fix:增加双引号处理机制，防止基因bed输出为空，匹配不上
        之前忽略了TEtranscripts在输出计数时，基因行id会加上双引号，而转座子不会。之前用python合并样本表达值时没出现，
        那是因为pandas读取文件时会默认去除字符外面的双引号。
