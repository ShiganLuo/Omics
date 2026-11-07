# 转座元件 (TE) 统计报告

输入 GTF：3 个；层级链 `class_id -> family_id -> gene_id (subfamily)`。

## 1. 总览

| 物种/组装 | 记录总数 | Class 数 | Family 数 | Subfamily 数 | TE 总长 (Mb) |
|---|---:|---:|---:|---:|---:|
| Human_GRCh38 | 4,832,405 | 13 | 125 | 1278 | 1584.9 |
| Mouse_GRCm39 | 3,799,629 | 10 | 102 | 1226 | 1121.5 |
| Rhesus_rheMac10 | 4,536,790 | 12 | 97 | 1274 | 1492.6 |

## 2. 各物种 Class 水平构成

### Human_GRCh38

| Class | 拷贝数 | 总长 (Mb) | Family 数 |
|---|---:|---:|---:|
| SINE | 1,901,998 | 422.00 | 3 |
| LINE | 1,614,481 | 681.31 | 10 |
| LTR | 770,551 | 286.15 | 6 |
| DNA | 512,404 | 107.92 | 20 |
| Satellite | 9,133 | 79.84 | 5 |
| Retroposon | 5,974 | 4.64 | 1 |
| LTR? | 5,777 | 1.21 | 1 |
| Unknown | 5,732 | 0.78 | 77 |
| DNA? | 3,364 | 0.46 | 3 |
| RC | 1,820 | 0.38 | 1 |
| RNA | 721 | 0.12 | 1 |
| RC? | 417 | 0.06 | 1 |
| SINE? | 33 | 0.00 | 1 |

### Mouse_GRCm39

| Class | 拷贝数 | 总长 (Mb) | Family 数 |
|---|---:|---:|---:|
| SINE | 1,492,624 | 202.27 | 6 |
| LTR | 1,105,064 | 331.86 | 6 |
| LINE | 976,368 | 542.90 | 10 |
| DNA | 164,783 | 28.75 | 21 |
| Satellite | 37,228 | 5.14 | 3 |
| Retroposon | 17,253 | 7.59 | 1 |
| Unknown | 5,786 | 2.88 | 53 |
| RC | 416 | 0.07 | 1 |
| Unspecified | 98 | 0.01 | 1 |
| DNA? | 9 | 0.00 | 2 |

### Rhesus_rheMac10

| Class | 拷贝数 | 总长 (Mb) | Family 数 |
|---|---:|---:|---:|
| SINE | 1,819,940 | 405.73 | 3 |
| LINE | 1,477,445 | 619.99 | 9 |
| LTR | 711,148 | 259.93 | 6 |
| DNA | 498,518 | 104.31 | 21 |
| Satellite | 11,341 | 99.25 | 4 |
| LTR? | 7,669 | 1.50 | 1 |
| Unknown | 4,427 | 0.78 | 51 |
| DNA? | 3,906 | 0.60 | 4 |
| RC | 1,726 | 0.36 | 1 |
| RC? | 363 | 0.06 | 1 |
| Retroposon | 162 | 0.05 | 1 |
| Unspecified | 145 | 0.02 | 1 |

## 3. 跨物种 Class 拷贝数对比

| Class | Human_GRCh38 | Mouse_GRCm39 | Rhesus_rheMac10 | 合计 |
|---|---|---|---|---|
| SINE | 1,901,998 | 1,492,624 | 1,819,940 | 5,214,562 |
| LINE | 1,614,481 | 976,368 | 1,477,445 | 4,068,294 |
| LTR | 770,551 | 1,105,064 | 711,148 | 2,586,763 |
| DNA | 512,404 | 164,783 | 498,518 | 1,175,705 |
| Satellite | 9,133 | 37,228 | 11,341 | 57,702 |
| Retroposon | 5,974 | 17,253 | 162 | 23,389 |
| Unknown | 5,732 | 5,786 | 4,427 | 15,945 |
| LTR? | 5,777 | 0 | 7,669 | 13,446 |
| DNA? | 3,364 | 9 | 3,906 | 7,279 |
| RC | 1,820 | 416 | 1,726 | 3,962 |
| RC? | 417 | 0 | 363 | 780 |
| RNA | 721 | 0 | 0 | 721 |
| Unspecified | 0 | 98 | 145 | 243 |
| SINE? | 33 | 0 | 0 | 33 |

## 4. 各物种主要 Family (Top 15)

### Human_GRCh38

| Family | Class | 拷贝数 | 总长 (Mb) | Subfamily 数 |
|---|---|---:|---:|---:|
| Alu | SINE | 1,282,821 | 333.25 | 54 |
| L1 | LINE | 1,031,524 | 554.71 | 132 |
| MIR | SINE | 616,589 | 88.48 | 5 |
| L2 | LINE | 486,431 | 108.92 | 17 |
| ERVL-MaLR | LTR | 367,470 | 117.61 | 86 |
| hAT-Charlie | DNA | 270,942 | 48.38 | 83 |
| ERV1 | LTR | 188,489 | 91.38 | 306 |
| ERVL | LTR | 174,276 | 61.14 | 128 |
| TcMar-Tigger | DNA | 123,237 | 38.42 | 81 |
| CR1 | LINE | 69,603 | 12.46 | 26 |
| hAT-Tip100 | DNA | 49,880 | 9.23 | 53 |
| Gypsy | LTR | 24,975 | 5.13 | 24 |
| hAT-Blackjack | DNA | 20,062 | 3.48 | 9 |
| TcMar-Mariner | DNA | 16,871 | 2.93 | 6 |
| RTE-X | LINE | 15,944 | 3.65 | 4 |

### Mouse_GRCm39

| Family | Class | 拷贝数 | 总长 (Mb) | Subfamily 数 |
|---|---|---:|---:|---:|
| L1 | LINE | 877,139 | 528.17 | 149 |
| Alu | SINE | 558,808 | 66.19 | 15 |
| ERVL-MaLR | LTR | 534,111 | 128.71 | 65 |
| B4 | SINE | 380,688 | 57.41 | 4 |
| B2 | SINE | 370,294 | 59.68 | 5 |
| ERVK | LTR | 370,031 | 138.11 | 285 |
| MIR | SINE | 120,815 | 14.62 | 5 |
| ERVL | LTR | 113,067 | 31.93 | 116 |
| hAT-Charlie | DNA | 105,223 | 18.35 | 80 |
| L2 | LINE | 75,565 | 11.85 | 14 |
| ERV1 | LTR | 69,131 | 20.58 | 121 |
| ID | SINE | 60,640 | 4.24 | 5 |
| Satellite | Satellite | 32,960 | 4.56 | 6 |
| TcMar-Tigger | DNA | 27,962 | 5.26 | 57 |
| L1-dep | Retroposon | 17,253 | 7.59 | 3 |

### Rhesus_rheMac10

| Family | Class | 拷贝数 | 总长 (Mb) | Subfamily 数 |
|---|---|---:|---:|---:|
| Alu | SINE | 1,248,211 | 323.88 | 44 |
| L1 | LINE | 902,425 | 497.98 | 136 |
| MIR | SINE | 569,270 | 81.58 | 5 |
| L2 | LINE | 485,746 | 105.72 | 14 |
| ERVL-MaLR | LTR | 336,312 | 106.71 | 86 |
| hAT-Charlie | DNA | 252,393 | 45.23 | 82 |
| ERV1 | LTR | 170,291 | 79.44 | 321 |
| ERVL | LTR | 159,361 | 55.43 | 124 |
| TcMar-Tigger | DNA | 121,578 | 36.57 | 81 |
| CR1 | LINE | 64,201 | 11.48 | 24 |
| hAT-Tip100 | DNA | 62,571 | 11.72 | 54 |
| Gypsy | LTR | 23,025 | 4.68 | 25 |
| hAT-Blackjack | DNA | 16,772 | 2.47 | 10 |
| TcMar-Mariner | DNA | 15,670 | 2.79 | 6 |
| ERVK | LTR | 14,950 | 12.19 | 56 |

## 5. 各物种主要 Subfamily (Top 15)

### Human_GRCh38

| Subfamily | Family | Class | 拷贝数 | 总长 (Mb) |
|---|---|---|---:|---:|
| MIRb | MIR | SINE | 230,248 | 34.88 |
| MIR | MIR | SINE | 180,671 | 28.59 |
| L2a | L2 | LINE | 177,912 | 44.32 |
| L2c | L2 | LINE | 143,542 | 28.75 |
| AluJb | Alu | SINE | 136,483 | 33.74 |
| AluSx | Alu | SINE | 133,671 | 35.98 |
| AluSx1 | Alu | SINE | 126,342 | 36.03 |
| AluY | Alu | SINE | 114,541 | 33.44 |
| AluSz | Alu | SINE | 110,813 | 31.15 |
| MIRc | MIR | SINE | 104,038 | 13.49 |
| L2b | L2 | LINE | 100,555 | 18.93 |
| AluJr | Alu | SINE | 91,919 | 21.37 |
| MIR3 | MIR | SINE | 89,435 | 10.30 |
| AluJo | Alu | SINE | 82,993 | 18.42 |
| AluSq2 | Alu | SINE | 65,121 | 18.66 |

### Mouse_GRCm39

| Subfamily | Family | Class | 拷贝数 | 总长 (Mb) |
|---|---|---|---:|---:|
| B3 | B2 | SINE | 147,452 | 24.82 |
| RSINE1 | B4 | SINE | 109,830 | 13.78 |
| ID_B1 | B4 | SINE | 108,580 | 14.39 |
| B4A | B4 | SINE | 99,738 | 18.56 |
| B1_Mus1 | Alu | SINE | 96,512 | 12.96 |
| B3A | B2 | SINE | 91,151 | 13.76 |
| B2_Mm2 | B2 | SINE | 87,324 | 13.60 |
| B1_Mus2 | Alu | SINE | 70,143 | 9.08 |
| ORR1B1-int | ERVL-MaLR | LTR | 65,146 | 4.88 |
| B4 | B4 | SINE | 62,540 | 10.67 |
| PB1D10 | Alu | SINE | 61,994 | 5.98 |
| RLTR17B_Mm | ERVK | LTR | 60,405 | 8.99 |
| Lx8 | L1 | LINE | 57,160 | 18.87 |
| MTD | ERVL-MaLR | LTR | 54,578 | 14.68 |
| Lx9 | L1 | LINE | 46,584 | 15.09 |

### Rhesus_rheMac10

| Subfamily | Family | Class | 拷贝数 | 总长 (Mb) |
|---|---|---|---:|---:|
| MIRb | MIR | SINE | 212,246 | 32.11 |
| MIR | MIR | SINE | 167,944 | 26.53 |
| L2a | L2 | LINE | 157,201 | 40.17 |
| AluSz | Alu | SINE | 123,532 | 33.60 |
| L2c | L2 | LINE | 114,748 | 23.54 |
| AluJb | Alu | SINE | 112,282 | 27.46 |
| AluJr | Alu | SINE | 110,188 | 26.65 |
| AluY | Alu | SINE | 108,544 | 31.53 |
| AluSx | Alu | SINE | 105,228 | 29.38 |
| AluSx1 | Alu | SINE | 98,804 | 28.47 |
| MIRc | MIR | SINE | 95,619 | 12.38 |
| L2b | L2 | LINE | 84,628 | 16.23 |
| MIR3 | MIR | SINE | 82,269 | 9.44 |
| AluYRa1 | Alu | SINE | 71,566 | 21.24 |
| AluSp | Alu | SINE | 61,345 | 17.24 |

