suppressPackageStartupMessages({
  library(argparse)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
})
plot_sv_circos_from_files <- function(file_list, genome = "mm39", cytoband_file = NULL, show_legend = TRUE, outImg = "sv_circos.png") {
  
  # 1. 图像输出设置
  png(filename = outImg, width = 2400, height = 2400, res = 300)
  circos.clear()

  # 2. 初始化 Ideogram
  if (!is.null(cytoband_file)) {
    cyto_data <- read.cytoband(cytoband_file)$df
    circos.initializeWithIdeogram(cyto_data)
  } else {
    tryCatch({
      circos.initializeWithIdeogram(species = genome)
    }, error = function(e) {
      stop("Genome '", genome, "' not found. Please provide a cytoband_file.")
    })
  }
  
  # 内部辅助函数：确保染色体以 "chr" 开头
  fix_chr <- function(df) {
    if (nrow(df) == 0) return(df)
    df[[1]] <- ifelse(grepl("^chr", df[[1]]), df[[1]], paste0("chr", df[[1]]))
    if (ncol(df) >= 3 && is.character(df[[3]])) { # 针对 TRA 的第二坐标列
       df[[3]] <- ifelse(grepl("^chr", df[[3]]), df[[3]], paste0("chr", df[[3]]))
    }
    return(df)
  }

  # 3. BLOCKS (DEL, DUP, INV)
  block_configs <- list(
    list(type = "DEL", col = "#FF413680"),
    list(type = "DUP", col = "#0074D980"),
    list(type = "INV", col = "#B10DC980")
  )

  for (cfg in block_configs) {
    if (!is.null(file_list[[cfg$type]]) && file.exists(file_list[[cfg$type]])) {
      data <- fix_chr(read.table(file_list[[cfg$type]], header = FALSE))
      circos.genomicTrack(data, ylim = c(0, 1), panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, ybottom = 0, ytop = 1, col = cfg$col, border = NA)
      }, track.height = 0.06)
    }
  }
  
  # 4. POINTS (INS)
    if (!is.null(file_list$INS) && file.exists(file_list$INS)) {
        ins_data <- fix_chr(read.table(file_list$INS, header = FALSE))
        
        # 强制设为 1bp 范围，防止因 VCF 的 END 导致横向占满
        ins_clean <- data.frame(
        chr = ins_data[[1]],
        start = ins_data[[2]],
        end = ins_data[[2]]
        )
        
        # 使用透明度颜色，重叠越多颜色越深，视觉效果更好
        ins_col <- "#2ECC4044" 

        circos.genomicTrack(ins_clean, ylim = c(0, 1), panel.fun = function(region, value, ...) {
        # 为当前扇区（Sector）的点生成随机 Y 坐标实现抖动
        # nrow(region) 获取当前轨道在该扇区的点数
        y_jitter = runif(nrow(region), min = 0.1, max = 0.9)
        
        # 绘制抖动后的点
        circos.genomicPoints(region, y_jitter, col = ins_col, pch = 16, cex = 0.1)
        }, track.height = 0.08) # 稍微加高轨道，给抖动留出空间
    }

    
    # 5. LINKS (TRA)
  tra_col <- "#FF000080" # 改为红色半透明，更容易看见
  if (!is.null(file_list$TRA) && file.exists(file_list$TRA)) {
    tra_raw <- fix_chr(read.table(file_list$TRA, header = FALSE))
    
    # 构建 region 对象
    region1 <- data.frame(chr = tra_raw[,1], start = tra_raw[,2], end = tra_raw[,2])
    region2 <- data.frame(chr = tra_raw[,3], start = tra_raw[,4], end = tra_raw[,4])
    
    # 过滤掉不在 Ideogram 中的染色体防止报错
    sectors <- get.all.sector.index()
    keep <- region1$chr %in% sectors & region2$chr %in% sectors
    
    if (any(keep)) {
      circos.genomicLink(region1[keep,], region2[keep,], col = tra_col, lwd = 0.8)
    } else {
      warning("No TRA links were plotted. Check chromosome naming (e.g., 'chr1' vs '1').")
    }
  }

  # 6. 添加图例
  if (show_legend) {
    lgd_list = Legend(
      labels = c("Deletion", "Duplication", "Inversion", "Insertion", "Translocation"),
      type = c("box", "box", "box", "points", "lines"),
      legend_gp = gpar(fill = c("#FF413680", "#0074D980", "#B10DC980", NA, NA), 
                       col = c(NA, NA, NA, "#2ECC40", tra_col)),
      pch = c(NA, NA, NA, 16, NA),
      title = "SV Types",
      nrow = 5 # 垂直排列适合放在侧面
    )
    # 绘制在右上方
    draw(lgd_list, x = unit(0.95, "npc"), y = unit(0.95, "npc"), just = c("right", "top"))
  }

  dev.off()
}

parser <- ArgumentParser(description='High-performance Circos Plotter for SV visualization')

# 设置参数组
group_input <- parser$add_argument_group("Input/Output Arguments")
group_input$add_argument("-i", "--input_dir", required=TRUE, help="Directory containing .bed files (DEL, TRA, INS, etc.)")
group_input$add_argument("-o", "--output", default="sv_circos.png", help="Output image path [default %(default)s]")

group_genome <- parser$add_argument_group("Genome Configuration")
group_genome$add_argument("-g", "--genome", default="mm39", help="Genome version (mm39, hg38, etc.) [default %(default)s]")
group_genome$add_argument("-c", "--cytoband", help="Path to local cytoband file (optional)")

group_visual <- parser$add_argument_group("Visual Settings")
group_visual$add_argument("--no_legend", action="store_true", help="Disable legend rendering")
group_visual$add_argument("-s", "--size", type="integer", default=2400, help="Image size in pixels [default %(default)s]")
group_visual$add_argument("-r", "--res", type="integer", default=300, help="Resolution DPI [default %(default)s]")

# 解析参数
args <- parser$parse_args()

main <- function() {
  # 构建文件路径列表
  file_list <- list(
    DEL = file.path(args$input_dir, "blocks_del.bed"),
    TRA = file.path(args$input_dir, "links_tra.bed"),
    INS = file.path(args$input_dir, "points_ins.bed"),
    DUP = file.path(args$input_dir, "points_dup.bed"),
    INV = file.path(args$input_dir, "points_inv.bed")
  )
  
  # 调用绘图函数
  plot_sv_circos_from_files(
    file_list = file_list,
    genome = args$genome,
    cytoband_file = args$cytoband,
    show_legend = !args$no_legend,
    outImg = args$output
  )
}
main()

# file_list <- list(
#   DEL = "/data/pub/zhousha/Totipotent20251031/PacBio/circos/DMSO06/blocks_del.bed",
#   TRA = "/data/pub/zhousha/Totipotent20251031/PacBio/circos/DMSO06/links_tra.bed",
#   INS = "/data/pub/zhousha/Totipotent20251031/PacBio/circos/DMSO06/points_ins.bed",
#   DUP = "/data/pub/zhousha/Totipotent20251031/PacBio/circos/DMSO06/points_dup.bed",
#   INV = "/data/pub/zhousha/Totipotent20251031/PacBio/circos/DMSO06/points_inv.bed"
# )
# plot_sv_circos_from_files(file_list, genome = "mm39")