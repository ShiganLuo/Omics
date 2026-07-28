# 加载必要的包
library(ggpubr)
library(reshape2)
library(dplyr)

# 读取数据
data  <- read.csv("test.csv")
# 将数据从宽格式转换为长格式
data_long  <- melt(data, variable.name =  "Group", value.name =  "Value")
# 设置组的顺序（对照组在前）
data_long$Group <- factor(data_long$Group, levels = c("IgG",  "aITGA11",  "aCHI3L1",  "aITGA11_aCHI3L1"))
# 进行组间比较的统计检验
stat_test  <- compare_means(Value ~ Group, data = data_long,  
                                         ref.group =  "IgG", method =  "t.test")
mycolor  <- c("#0077c1",  "#00a99e",  "#6bc72b",  "#ff5a20")
# 绘制箱形图
p  <- ggboxplot(data_long, x =  "Group", y =  "Value",  
                       fill  =  "Group",
                       color  =  "Group",
                       width  = 0.45,
                       palette  = c("#0077c1",  "#00a99e",  "#6bc72b",  "#ff5a20"),
                       add  =  "jitter",
                       add.params = list(size = 1.5, alpha = 0.7),
                       xlab  =  "",
                       ylab  =  "LYVE-1 positive vessel/tumor") +
scale_colour_manual(values=alpha(mycolor,1))+
scale_fill_manual(values=alpha(mycolor,0.1))+
theme(legend.position =  "none",
             axis.text.x = element_text(angle = 45, hjust = 1))

# 添加显著性标记(星号)
p  + stat_pvalue_manual(stat_test,  
                                   label  =  "p.signif",
                                   y.position = max(data_long$Value) * 1.05,
                                   step.increase = 0.1)

# 添加显著性标记(pvalue)
p  + stat_pvalue_manual(stat_test,  
                                   label  =  "p.adj",
                                   y.position = max(data_long$Value) * 1.05,
                                   step.increase = 0.1)