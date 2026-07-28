# 加载必要的包
library(Biobase)
library(Mfuzz)
library(tidyverse)
# 1. 读取数据
group_data <- read_csv("group.csv")
expression_data <- read_csv("testdata.csv")
# 2. 计算组内均值
group_mean <- expression_data %>%
  pivot_longer(-gene, names_to = "Sample", values_to = "Expression") %>%
  left_join(group_data, by = "Sample") %>%
  group_by(gene, Group) %>%
  summarise(Mean = mean(Expression)) %>%
  pivot_wider(names_from = Group, values_from = Mean) %>%
  column_to_rownames("gene")
# 3. 创建ExpressionSet对象
eset <- ExpressionSet(assayData = as.matrix(group_mean))
# 4. 过滤缺失值（超过25%组别缺失的基因）
eset <- filter.NA(eset, thres = 0.25)
# 5. 用均值填充缺失值
eset <- fill.NA(eset, mode = "mean")
# 根据最小标准差过滤基因；
eset <- filter.std(eset,min.std=0.1)
# 6. 标准化数据（Mfuzz要求标准化）
eset <- standardise(eset)
# 7. 设置聚类参数
#使用mestimate函数估计m值；
m1 <- mestimate(eset)
set.seed(123)
cl <- mfuzz(eset, c = 9, m = m1)
# 8. 可视化聚类结果（3x3布局）
mfuzz.plot(
  eset,
  cl,
  mfrow = c(3, 3),
new.window = FALSE,
  time.labels = colnames(eset))