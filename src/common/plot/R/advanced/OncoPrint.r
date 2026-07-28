library(ComplexHeatmap)
library(grid)
sv_plot = function(
    oncoprint_matrix_file,
    outpng = "oncoprint_plot.png",
    alteration_colors = c(
        "INS" = "#e41a1c",
        "DEL" = "#1f77b4",
        "DUP" = "#ff7f0e",
        "INV" = "#2ca02c",
        "TRA" = "#ffcc00",
        "OTHER" = "#9467bd"
    )
) {


    mat = read.csv(oncoprint_matrix_file, row.names = 1, check.names = FALSE)

    # 提取实际存在的SV类型
    sv_types = unique(unlist(strsplit(as.vector(as.matrix(mat)), ",")))
    sv_types = sv_types[sv_types != "" & !is.na(sv_types)]

    # 仅保留存在的颜色
    col = alteration_colors[names(alteration_colors) %in% sv_types]

    # 构建 alter_fun
    alter_fun = list(
        background = function(x, y, w, h) {
            grid.rect(x, y, w, h, gp = gpar(fill = "white", col = NA))
        }
    )

    n = length(sv_types)

    for (i in seq_along(sv_types)) {

        sv = sv_types[i]

        alter_fun[[sv]] = local({
            svtype = sv
            idx = i

            function(x, y, w, h) {

                h_unit = h / n
                y_pos = y - h/2 + h_unit * (idx - 0.5)

                grid.rect(
                    x,
                    y_pos,
                    w * 0.9,
                    h_unit * 0.9,
                    gp = gpar(fill = col[svtype], col = NA)
                )
            }
        })
    }

    png(outpng, width = 2000, height = 1500, res = 300, bg = "white")
    ht = oncoPrint(
        mat,
        get_type = function(x) strsplit(x, ",")[[1]],
        alter_fun = alter_fun,
        col = col,
        show_column_names = TRUE,
        remove_empty_columns = TRUE,
        remove_empty_rows = TRUE,
        column_names_rot = 0,
        column_names_centered = TRUE,
        column_names_gp = gpar(
            fontsize = 12      
        ),
        column_order = colnames(mat),
        row_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(
            title = "Alterations",
            at = names(col), 
            labels = names(col),
            # 2. 修复颜色掩盖问题：将内部变量名重命名为 type_nm 避免与坐标轴 x 冲突
            graphics = lapply(names(col), function(type_nm) {
                local({
                    t_color = col[type_nm]
                    function(x, y, w, h) {
                        grid.rect(x, y, w * 0.8, h * 0.8, gp = gpar(fill = t_color, col = NA))
                    }
                })
            }),
            grid_height = unit(5, "mm"), 
            grid_width = unit(5, "mm"),   
            gap = unit(2, "mm"),          
            labels_gp = gpar(fontsize = 10) 
        )
    )

    draw(ht, heatmap_legend_side = "right")

    dev.off()
}
sv_plot(
    oncoprint_matrix_file = "/data/pub/zhousha/Totipotent20251031/PacBio/OncoPrint/oncoprint_matrix.csv",
    outpng = "/data/pub/zhousha/Totipotent20251031/PacBio/OncoPrint/plot/oncoprint_plot.png"
)

