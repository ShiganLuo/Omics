library(sf)
library(ggplot2)
###---物种和地图数据---###
###地图数据
Jiujiang = st_read("./datum/map/Jiujiang.json")
country=read.csv("./datum/species_distribution/Jiujiang_tude.csv",sep="\t",header=TRUE)#获取区县坐标信息，便于添加文本
###物种数据
species = read.csv("./datum/species_distribution/species.csv",sep=" ",header=TRUE)
# head(species)
# species[,c("scientific_name")]
####读取每个物种学名，以好循环
species_name = read.table('./datum/species_distribution/scientific_name.txt', sep = '\t', stringsAsFactors = FALSE)
chinese_name = read.table('./datum/species_distribution/chinese_name.txt',sep = '\t', stringsAsFactors = FALSE)
species_name = species_name$V1#species_names必须是个向量
chinese_name = chinese_name$V1
i=1
###---绘制图形---###
for(a in species_name){
    #从species中获取某个物种的数据
    specie = species[species$scientific_name == a,]
    filename = paste0("./output/map/",chinese_name[i])
    filename = paste0(filename,".jpeg")
    titlename = paste0(chinese_name[i],"分布地点")
    jpeg(width = 320, height = 200, units = "mm", res = 200, file = filename)
    p = ggplot(data = Jiujiang) +
            geom_sf(color = "black", fill = "lightgrey")+
            geom_point(data=specie, aes(x=decimalLongitude, y=decimalLatitude),
                        shape=21,color="black", fill="red",alpha=0.6)+
            xlab("Longitude") + ylab("Latitude")+
            ggtitle(titlename)+
            theme_bw()+
            geom_text(data = country, aes(x=decimalLongitude, y=decimalLatitude, label = country_names), size = 3, nudge_y = 0.1)  # 添加县级区域名称
    print(p)
    dev.off()
    i = i+1

}
