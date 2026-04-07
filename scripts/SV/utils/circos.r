suppressPackageStartupMessages({
  library(argparse)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
})
log_msg <- function(level = c("INFO","WARN","ERROR"), ..., quit = FALSE) {
  level <- match.arg(level)
  msg <- paste(...)
  prefix <- paste0(
    "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
    "[", level, "] "
  )
  message(prefix, msg)
  if (quit) stop(msg, call. = FALSE)
}
plot_sv_circos_from_files <- function(
  file_list, 
  genome = "mm39", 
  cytoband_file = NULL, 
  show_legend = TRUE, 
  outImg = "sv_circos.png",
  ins_plot_type = "points" , # "points" or "bar"
  bin_size = 1e6 # 1Mb bins for bar plot
  ) {
  
  # 1. 图像输出设置
  png(filename = outImg, width = 2400, height = 2400, res = 300)
  circos.clear()
  # cytoband = read.cytoband(species = genome)$df
  # write.table(cytoband, file = "debug_cytoband.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  # 2. 初始化 Ideogram
  if (!is.null(cytoband_file)) {
    cyto_data <- read.cytoband(cytoband_file)$df
    circos.initializeWithIdeogram(cyto_data)
  } else {
    tryCatch({
      circos.initializeWithIdeogram(species = genome)
    }, error = function(e) {
      log_msg("ERROR", paste("Failed to load cytoband data for genome:", genome), quit = TRUE)
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
    if (!is.null(file_list[[cfg$type]]) && file.exists(file_list[[cfg$type]]) && file.info(file_list[[cfg$type]])$size > 0) {
      data <- fix_chr(read.table(file_list[[cfg$type]], header = FALSE))
      circos.genomicTrack(data, ylim = c(0, 1), panel.fun = function(region, value, ...) {
        circos.genomicPoints(region, runif(nrow(region)), ybottom = 0, ytop = 1, col = cfg$col, border = NA,pch=16,cex=0.5)
      }, track.height = 0.06)
    } else {
      log_msg("INFO", paste("No", cfg$type, "data found or file is empty. Skipping this track."))
    }
  }
  
  # 4. POINTS (INS)
    if (!is.null(file_list$INS) && file.exists(file_list$INS) && file.info(file_list$INS)$size > 0) {
        ins_data <- fix_chr(read.table(file_list$INS, header = FALSE))
        
        # 强制设为 1bp 范围，防止因 VCF 的 END 导致横向占满
        ins_clean <- data.frame(
        chr = ins_data[[1]],
        start = ins_data[[2]],
        end = ins_data[[2]]
        )
        if (ins_plot_type == "points") {
          # 使用透明度颜色，重叠越多颜色越深，视觉效果更好
          ins_col <- "#2ECC4044" 

          circos.genomicTrack(ins_clean, ylim = c(0, 1), panel.fun = function(region, value, ...) {
          # 为当前扇区（Sector）的点生成随机 Y 坐标实现抖动
          # nrow(region) 获取当前轨道在该扇区的点数
          y_jitter = runif(nrow(region), min = 0.1, max = 0.9)
          
          # 绘制抖动后的点
          circos.genomicPoints(region, y_jitter, col = ins_col, pch = 16, cex = 0.1)
          }, track.height = 0.08) # 稍微加高轨道，给抖动留出空间
        } else if (ins_plot_type == "bar") {
          if (bin_size <= 0 || !is.numeric(bin_size)) {
                  log_msg("ERROR", "bin_size must be positive number for bar plot.", quit = TRUE)  
              }
              
              # 1. 计算 Bins，并显式确保 start < end
              ins_bins <- ins_clean |>
                  dplyr::group_by(
                      chr,
                      bin = floor(start / bin_size)
                  ) |>
                  dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
                  dplyr::mutate(
                      start = as.numeric(bin * bin_size),
                      end   = as.numeric((bin + 1) * bin_size) # 使用 bin+1 确保区间跨度
                  ) |>
                  # circlize 默认前三列为 genomic regions
                  dplyr::select(chr, start, end, count) |>
                  # 确保 start 永远小于 end (预防万一)
                  dplyr::filter(end > start)

              # 4. 调试输出：检查是否有负值或异常
              # write.table(ins_bins, file = "debug_ins_bins.txt", ...)
              max_val <- if(nrow(ins_bins) > 0) max(ins_bins$count) else 1
              # 5. 绘图
              circos.genomicTrack(
                ins_bins, 
                ylim = c(0, max_val * 1.2), 
                track.height = 0.12,
                panel.fun = function(region, value, ...) {
                  # type = "h" 是关键，它会从每个 bin 的位置向底部画垂线. genomicReact 画矩形在这里不合适
                  circos.genomicLines(
                    region, 
                    value, 
                    type = "h", 
                    col = "#2ECC40", 
                    lwd = 1   # 如果柱子太稀疏，调大 lwd；太挤则调小
                  )
                }
              )
        } else {
          log_msg("WARN", paste("Unknown ins_plot_type:", ins_plot_type, ". Skipping INS points."))
        }

    } else {
        log_msg("INFO", "No INS data found or file is empty. Skipping INS points.")
    }

    
    # 5. LINKS (TRA)
  tra_col <- "#FF000080" # 改为红色半透明，更容易看见
  if (!is.null(file_list$TRA) && file.exists(file_list$TRA) && file.info(file_list$TRA)$size > 0) {
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
  } else {
    log_msg("INFO", "No TRA data found or file is empty. Skipping TRA links.")
  }

  # 6. 添加图例
  if (show_legend) {
    lgd_list = Legend(
      labels = c("Deletion", "Duplication", "Inversion", "Insertion", "Translocation"),
      type   = "points",
      pch    = rep(15, 5),            # 实心方块
      size   = unit(4, "mm"),         # 控制色块大小
      legend_gp = gpar(
        col = c(
          "#FF4136",  # Deletion
          "#0074D9",  # Duplication
          "#B10DC9",  # Inversion
          "#2ECC40",  # Insertion
          tra_col     # Translocation
        )
      ),
      title = "SV Types",
      nrow  = 5,
      background = "white"            # 明确白底
    )

    draw(
      lgd_list,
      x = unit(0.95, "npc"),
      y = unit(0.95, "npc"),
      just = c("right", "top")
    )
  }


  dev.off()
}

parse_bin_size <- function(x) {
  x <- tolower(x)
  if (grepl("k$", x)) return(as.integer(as.numeric(sub("k$", "", x)) * 1e3))
  if (grepl("m$", x)) return(as.integer(as.numeric(sub("m$", "", x)) * 1e6))
  if (grepl("g$", x)) return(as.integer(as.numeric(sub("g$", "", x)) * 1e9))
  as.integer(as.numeric(x))  # 支持 1e6
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
group_ins <- parser$add_argument_group("Insertion (INS) Plot Settings")
group_ins$add_argument("-n","--ins_plot_type", choices=c("points","bar"), default="points", help="INS plot type: points or bar [default %(default)s]")
group_ins$add_argument("-b","--ins_bin_size", type="character", default=1e6, help="Bin size for INS bar plot (e.g. 500k, 1e6, 2M)")

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
  args$ins_bin_size <- parse_bin_size(args$ins_bin_size)
  log_msg("INFO", paste("INS bar plot bin size set to:", args$ins_bin_size))
  # 调用绘图函数
  plot_sv_circos_from_files(
    file_list = file_list,
    genome = args$genome,
    cytoband_file = args$cytoband,
    show_legend = !args$no_legend,
    outImg = args$output,
    ins_plot_type = args$ins_plot_type,
    bin_size = args$ins_bin_size
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