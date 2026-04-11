suppressPackageStartupMessages(library(XenofilteR))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(BiocParallel))
# 获取命令行参数
parser <- ArgumentParser(description='XenofilteR: 去除比对到人基因组的bam文件中小鼠基因组的污染')
parser$add_argument('-i', '--inputFile', type='character', required=TRUE,
                    help='input bam sample csv.')
parser$add_argument('-o', '--outputDir', type='character', required=TRUE,
                    help='output dir.')
parser$add_argument('-r', '--renameSamples', nargs="+", type='character', required=FALSE,
                    help='file to rename output bam')
parser$add_argument('-m', '--MM', type='integer', default=4,
                    help='MM threshold,suggested by XenofilteR(150bp : 8,default == 4)')
parser$add_argument('-w', '--workers', type='integer', default=1,
                    help='The number of workers represents the number of CPUs used for the analysis and thereby the number of samples analyses simultaneously.')

# 解析命令行参数
args <- parser$parse_args()
csv_file <- args$inputFile  # 第一个参数：CSV 文件路径
output_folder <- args$outputDir # 第二个参数：输出文件夹路径

MM <- as.numeric(args$MM)  # 第四个参数：MM,根据read长度调整,150默认为8
workers <- as.numeric(args$workers)  # 第五个参数：workers,默认为1

if (!dir.exists(output_folder)) {
    cat("Output directory does not exist. Creating it: ", output_folder, "\n")
    dir.create(output_folder, recursive = TRUE)
}
# 创建并设置并行处理的参数
bp.param <- SnowParam(workers = workers, type = "SOCK")

# 执行 XenofilteR 函数
if(is.null(args$renameSamples)){
    print("---1----")
    sample <- read.table(csv_file, sep = ",", header = FALSE,stringsAsFactors = FALSE)
    XenofilteR(sample.list = sample, destination.folder = output_folder, bp.param = bp.param, MM_threshold = MM)
    
} else {
    print("---2----")
    sample <- read.table(csv_file, sep = ",", header = FALSE, stringsAsFactors = FALSE)
    output <- args$renameSamples
    XenofilteR(sample.list = sample, destination.folder = output_folder, bp.param = bp.param, output.names = output, MM_threshold = MM)
    
}




