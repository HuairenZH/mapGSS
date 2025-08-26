#!/usr/bin/env Rscript

# ========= 命令行参数：--gene =========
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default=NULL){
  i <- which(args == flag)
  if (length(i) == 1 && i < length(args)) return(args[i+1])
  default
}
gene <- get_arg("--gene")
if (is.null(gene) || gene == "") stop("请用 --gene 指定基因，例如：--gene SLC32A1")

# ========= 解析脚本所在目录 & 项目根目录 =========
get_script_dir <- function(){
  cmd <- commandArgs(trailingOnly = FALSE)
  i <- grep("^--file=", cmd)
  if (length(i) == 1) return(dirname(normalizePath(sub("^--file=", "", cmd[i]))))
  if ("rstudioapi" %in% rownames(installed.packages())) {
    if (rstudioapi::isAvailable()) {
      p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error=function(e) "")
      if (nzchar(p)) return(dirname(normalizePath(p)))
    }
  }
  getwd()
}
SCRIPT_DIR   <- normalizePath(get_script_dir(), winslash = "/", mustWork = FALSE)
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), winslash = "/", mustWork = FALSE)

# ========= 便捷工具 =========
gcol  <- function(suffix) paste0(gene, "_", suffix)      # gcol("GSS") -> "<gene>_GSS"
glab  <- function(suffix) paste(gene, suffix)            # glab("GSS") -> "<gene> GSS"
gpath <- function(prefix="mouse_E16.5_E1S1", ext="csv",
                  base=file.path(PROJECT_ROOT, "data", "Genes")) {
  file.path(base, gene, paste0(prefix, "_", gene, ".", ext))
}

# ========= 依赖 =========
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr)
  library(ggpubr);  library(ggridges); library(RColorBrewer)
  library(forcats); library(stringr);  library(rlang); library(scales)
  library(patchwork); library(grid)
})

# ========= 字体 =========
library(extrafont)
loadfonts() # 加载字体以便使用
theme(text = element_text(family = "Times New Roman"))

# ========= 全局缩放 =========
# 你的原脚本里有 scale_factor（=2），很多字号已乘它。
# 这里新增 ui_scale=0.25，用于将“所有视觉元素”等比例缩小为原来的 1/4。
scale_factor <- 2
ui_scale <- 0.25

# ========= 全局主题 & 配色 =========
theme_set(theme_bw())
theme_update(
  plot.title  = element_text(family="Times New Roman", size=25*scale_factor*ui_scale, face="bold", colour="black"),
  axis.title.x= element_text(family="Times New Roman", size=25*scale_factor*ui_scale, colour="black"),
  axis.title.y= element_text(family="Times New Roman", size=25*scale_factor*ui_scale, colour="black"),
  axis.text.x = element_text(family="Times New Roman", size=25*scale_factor*ui_scale, colour="black"),
  axis.text.y = element_text(family="Times New Roman", size=25*scale_factor*ui_scale, colour="black"),
  legend.text = element_text(family="Times New Roman", size=25*scale_factor*ui_scale),
  legend.title= element_text(family="Times New Roman", size=25*scale_factor*ui_scale),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  legend.key = element_blank(),
  strip.background = element_blank(),
  strip.text = element_text(family="Times New Roman", color="black", size=25*scale_factor*ui_scale),
  panel.border = element_blank(),
  axis.line = element_line(size = 0.5*ui_scale),
  legend.position = 'right'
)
theme_paper <- function(base=18){
  theme_classic(base_size = base*ui_scale) +
    theme(plot.title = element_text(family="Times New Roman", face="bold", hjust=0.5),
          axis.text  = element_text(family="Times New Roman", color="black"),
          axis.title.x = element_text(family="Times New Roman", margin = margin(t=8)),
          plot.margin  = margin(10,14,10,14))
}

color_pval <- rev(c("#a50026","#d73027","#f46d43","#fdae61",
                    "#fee090","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695"))
Colormap <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(32)
blue <- "#17BECF"; red <- "#E74C3C"

# ========= 读数据 =========
embryo_csv <- gpath(prefix = "mouse_E16.5_E1S1", ext = "csv")
adult_csv  <- gpath(prefix = "mouse_Adult_Mouse_brain_cell_bin", ext = "csv")
if (!file.exists(embryo_csv)) stop("未找到文件：", embryo_csv, "\n请先运行 extract_gene_GSS 生成数据。")
if (!file.exists(adult_csv))  stop("未找到文件：", adult_csv,  "\n请先运行 extract_gene_GSS 生成数据。")

show_embryo <- suppressMessages(read_csv(embryo_csv, show_col_types = FALSE))
show_adult  <- suppressMessages(read_csv(adult_csv,  show_col_types = FALSE))

# ========= 列检查 =========
req_cols_embryo <- c("x","y","annotation","cell_name", gcol("GSS"))
req_cols_adult  <- c("x","y","annotation","region","cell_name", gcol("GSS"))
for (nm in req_cols_embryo) if (!nm %in% names(show_embryo)) stop("embryo 缺少列：", nm)
for (nm in req_cols_adult)  if (!nm %in% names(show_adult))  stop("adult  缺少列：", nm)

# ========= 1) 胚胎面板（左）=========
color_all <- c(
  "Adipose tissue"="#6d32e6ff","Adrenal gland"="#bf024fff","AGM"="#d147a3ff",
  "Blood vessel"="#b3a726ff","Bone"="#103a14ff","Brain"="#ef833aff",
  "Branchial arch"="#b38b5cff","Cartilage"="#35586dff","Cartilage primordium"="#3cb44bff",
  "Cavity"="#dfdce0ff","Choroid plexus"="#bd3addff","Connective tissue"="#0bd3b1ff",
  "Dermomyotome"="#ff4374ff","Dorsal root ganglion"="#b74c11ff","Epidermis"="#036df4ff",
  "Facial nerve"="#dd7936ff","GI tract"="#5c5ca6ff","Head mesenchyme"="#be9b72ff",
  "Heart"="#d3245aff","Inner ear"="#03fff4ff","Jaw and tooth"="#f062f9ff",
  "Kidney"="#62cfe8ff","Liver"="#c923b1ff","Lung"="#7ec136ff","Lung primordium"="#7ec136ff",
  "Meninges"="#dfca43ff","Mesenchyme"="#86733aff","Mesentery"="#e71d36ff",
  "Mesothelium"="#ff8383ff","Mucosal epithelium"="#2f7dd1ff","Muscle"="#af1041ff",
  "Neural crest"="#d4e727ff","Notochord"="#f07f92ff","Olfactory epithelium"="#60a9eaff",
  "Ovary"="#55afd9ff","Pancreas"="#739b1eff","Primitive gut tube"="#93a8edff",
  "Sclerotome"="#2d739eff","Smooth muscle"="#fc5151ff","Spinal cord"="#f9d5baff",
  "Submandibular gland"="#ab32e6ff","Surface ectoderm"="#428fe2ff",
  "Sympathetic nerve"="#cc5a0dff","Thymus"="#be50ffff","Urogenital ridge"="#9e7bffff"
)

obs_e <- show_embryo
obs_e[[gcol("GSS")]] <- suppressWarnings(as.numeric(obs_e[[gcol("GSS")]]))
obs_e <- obs_e |> filter(!is.na(.data[[gcol("GSS")]]))

figa <- ggplot(obs_e, aes(x=x, y=y))+
  geom_point(aes(col=annotation), size=0.01*ui_scale)+
  scale_color_manual(values=color_all)+
  ggtitle("E16.5")+
  theme(
    plot.title = element_text(size=15*scale_factor*ui_scale),
    axis.title = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
    legend.key.height = unit(1.5*ui_scale,"cm"), legend.position = "right", legend.title = element_blank(),
    panel.background = element_rect(fill = "black"),
    plot.background  = element_rect(fill = "white")
  )

fig_1 <- ggplot(obs_e, aes(x = x, y = y)) +
  ggtitle(gene) +
  geom_point(aes(col = .data[[gcol("GSS")]]), size = 2*ui_scale, shape = 16) +
  scale_color_gradientn(
    colours = color_pval,
    limits  = c(0, max(obs_e[[gcol("GSS")]], na.rm = TRUE)),
    labels  = scales::number_format(accuracy = 1)
  ) +
  guides(color = guide_colorbar(title.position = "top")) +
  labs(color = glab("GSS")) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    legend.title.align = 0.5,
    legend.position = c(0.045, 0.01),
    legend.direction = "horizontal",
    legend.justification = c(0.1, 0),
    legend.key.width  = unit(3*ui_scale, "cm"),
    legend.key.height = unit(2*ui_scale, "cm"),
    legend.text  = element_text(color = "white", size = 20*scale_factor*ui_scale),
    legend.title = element_text(color = "white", size = 20*scale_factor*ui_scale),
    panel.background = element_rect(fill = "black"),
    plot.background  = element_rect(fill = "white"),
    legend.background     = element_rect(fill = "transparent", colour = NA),
    legend.key            = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    plot.title = element_text(color = "white", size = 25*scale_factor*ui_scale, hjust = 0.5)
  )

# annotation 的密度哑铃 + 比值柱（胚胎）
dat_annot <- obs_e |> filter(!is.na(annotation))
gss_sym   <- rlang::sym(gcol("GSS"))

sum_annot <- dat_annot |>
  group_by(annotation) |>
  summarise(
    n_lt = sum((!!gss_sym) < 1, na.rm = TRUE),
    n_gt = sum((!!gss_sym) > 1, na.rm = TRUE),
    n_eq = sum((!!gss_sym) == 1, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    label = annotation,
    ratio = case_when(
      n_lt == 0 & n_gt > 0 ~ Inf,
      n_lt == 0 & n_gt == 0 ~ NA_real_,
      TRUE ~ n_gt / n_lt
    )
  ) |>
  arrange(desc(ratio)) |>
  mutate(label = factor(label, levels = rev(label)))

xmax_annot <- max(sum_annot$ratio[is.finite(sum_annot$ratio)], na.rm = TRUE)
if (!is.finite(xmax_annot)) xmax_annot <- 1
sum_annot2 <- sum_annot |>
  mutate(
    ratio_plot = ifelse(is.infinite(ratio), xmax_annot * 1.05, ratio),
    ratio_lab  = ifelse(is.infinite(ratio), "\u221E", scales::number(ratio, accuracy = 0.01))
  )

lab_lt <- paste(glab("GSS"), "< 1"); lab_gt <- paste(glab("GSS"), "> 1")

p_annot_db <- ggplot(sum_annot, aes(y = label)) +
  geom_segment(aes(x = n_lt, xend = n_gt, yend = label), color = "grey70", linewidth = 5*ui_scale) +
  geom_point(aes(x = n_lt, color = lab_lt), size = 15*ui_scale) +
  geom_point(aes(x = n_gt, color = lab_gt), size = 15*ui_scale) +
  scale_color_manual(name = "Number of spots",
                     values = setNames(c(blue, red), c(lab_lt, lab_gt))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = scales::label_number(big.mark=",")) +
  labs(title = "", x = "Count", y = NULL) +
  theme_paper(18) +
  theme(
    axis.text.y  = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale),
    axis.text.x  = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(size = 1.5*ui_scale),
    plot.title   = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale, face = "bold"),
    axis.title.x = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale, face = "bold"),
    axis.title.y = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale, face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale),
    legend.text  = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(color = "black", size = 2.5*ui_scale)
  ) +
  guides(color = guide_legend(override.aes = list(size = 12*ui_scale)))

# 右侧柱图
xmax_lim_annot <- xmax_annot * 1.35
p_annot_bar <- ggplot(sum_annot2, aes(y = label, x = ratio_plot, fill = ratio_plot)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ratio_lab),
            hjust = -0.05, nudge_x = 0.02 * xmax_lim_annot,
            size = 20*ui_scale, family = "Times New Roman") +
  scale_fill_gradient(low = "#04DAF2", high = "#E74C3C", name = "Spots ratio") +
  scale_x_continuous(limits = c(0, xmax_lim_annot),
                     expand = expansion(mult = c(0, 0.14))) +
  labs(x = paste0("ratio"), y = NULL) +
  theme_paper(18) +
  theme(
    axis.text.y  = element_blank(),
    axis.text.x  = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale),
    axis.ticks.x = element_line(size = 1.5*ui_scale),
    axis.title.x = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale, face = "bold"),
    axis.ticks.y = element_blank(),
    plot.margin  = margin(10, 14, 10, 4),
    legend.position = "none"
  ) +
  coord_cartesian(clip = "off")

fig4 <- p_annot_db + p_annot_bar + patchwork::plot_layout(widths = c(3, 1))

# 左侧总拼
blank_gap <- grid::nullGrob()
fig_1 <- fig_1 + coord_cartesian(clip = "off") + theme(plot.margin = margin(0, 0, 0, 100))
fig4  <- fig4  + coord_cartesian(clip = "off") + theme(plot.margin = margin(0, 50, 0, 0))
combo_plot <- ggarrange(fig_1, blank_gap, fig4, ncol=3, nrow=1, widths=c(1, 0.05, 1.75))

combo_name <- paste0("mouse_E16.5_E1S1_", gene)
assign(combo_name, combo_plot)

# ========= 2) 鼠脑面板（右）=========
show_a <- show_adult
show_a$annotation[show_a$annotation == 'Glu Neu']  <- 'Glu-neuron'
show_a$annotation[show_a$annotation == 'DA Neu']   <- 'Da-neuron'
show_a$annotation[show_a$annotation == 'GABA Neu'] <- 'Gaba-neuron'
show_a[[gcol("GSS")]] <- suppressWarnings(as.numeric(show_a[[gcol("GSS")]]))

colors_ct = c("Astr"='#377EB8',"Da-neuron"='#4DAF4A',"End"='#984EA3',"Epe"='#FF7F00',
              "Ery"='#d53e4f',"Gaba-neuron"='#A65628',"Glu-neuron"='#999999',"GN DG"='#66C2A5',
              "Meninge"='#FC8D62',"Mic"='#8DA0CB',"Olig"='#E78AC3',"OPC"='#A6D854',
              "SM"='#FFD92F',"Unknown"='#E5C494')
colors_rg = c("CA1"='#377EB8',"CA3"='#984EA3',"CAA"='#d53e4f',"Cortex"='#A65628',
              "DG"='#F781BF',"Ependymal"='#999999',"FT"='#66C2A5',"Meninge"='#FC8D62',
              "Midb"='#8DA0CB',"PAN"='#E78AC3',"SL/R CA1"='#A6D854',"SN/VTA"='#FFD92F',
              "Subi"='#E5C494',"Thal"='#31A354')

# Brain region 图
figa_rg <- ggplot(show_a, aes(x = x, y = y)) +
  ggtitle("Brain region") +
  geom_point(aes(col = region), size = 1*ui_scale) +
  scale_color_manual(values = colors_rg) +
  guides(col = guide_legend(override.aes = list(size = 15*ui_scale))) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    legend.key.height = unit(1.5*ui_scale, "cm"),
    legend.position = "left",
    legend.title = element_blank(),
    legend.background     = element_rect(fill = "transparent", colour = NA),
    legend.key            = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "black"),
    plot.background  = element_rect(fill = "white"),
    plot.title = element_text(color = "white", size = 25*scale_factor*ui_scale, hjust = 0.5)
  )

# Cell-type 图
figa_ct <- ggplot(show_a, aes(x = x, y = y)) +
  ggtitle("Cell-type") +
  geom_point(aes(col = annotation), size = 1*ui_scale) +
  scale_color_manual(values = colors_ct) +
  guides(col = guide_legend(override.aes = list(size = 15*ui_scale))) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    legend.key.height = unit(1.5*ui_scale, "cm"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.background     = element_rect(fill = "transparent", colour = NA),
    legend.key            = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "black"),
    plot.background  = element_rect(fill = "white"),
    plot.title = element_text(color = "white", size = 25*scale_factor*ui_scale, hjust = 0.5)
  )

fige_2 <- ggplot(show_a, aes(x = x, y = y)) +
  ggtitle(gene) +
  geom_point(aes(col = .data[[gcol("GSS")]]), size = 2*ui_scale, shape = 16) +
  scale_color_gradientn(
    colours = color_pval,
    limits  = c(0, max(show_a[[gcol("GSS")]], na.rm = TRUE)),
    labels  = scales::number_format(accuracy = 1)
  ) +
  guides(color = guide_colorbar(title.position = "top")) +
  labs(color = glab("GSS")) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    legend.title.align = 0.5,
    legend.position = c(0.045, 0.01),
    legend.direction = "horizontal",
    legend.justification = c(0.1, 0),
    legend.key.width  = unit(3*ui_scale, "cm"),
    legend.key.height = unit(2*ui_scale, "cm"),
    legend.text  = element_text(color = "white", size = 20*scale_factor*ui_scale),
    legend.title = element_text(color = "white", size = 20*scale_factor*ui_scale),
    panel.background = element_rect(fill = "black"),
    plot.background  = element_rect(fill = "white"),
    legend.background     = element_rect(fill = "transparent", colour = NA),
    legend.key            = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    plot.title = element_text(color = "white", size = 25*scale_factor*ui_scale, hjust = 0.5)
  )

figall1 <- ggarrange(figa_rg, figa_ct, fige_2, ncol = 3, widths = c(1.25,1.25,1))

# region 的密度山脊图（成人）
group_data <- show_a
gss <- rlang::sym(gcol("GSS"))
fig1 <- ggplot(group_data, aes(x = !!gss, y = region)) +
  stat_density_ridges(aes(height = after_stat(density), fill = after_stat(density)),
                      geom = "density_ridges_gradient", scale = 1, rel_min_height = 0,
                      panel_scaling = FALSE, adjust = 1, linewidth = 0, color = NA) +
  stat_density_ridges(aes(height = after_stat(density)),
                      geom = "density_ridges", scale = 1, rel_min_height = 0,
                      panel_scaling = FALSE, adjust = 1, fill = NA, color = "black", linewidth = 3*ui_scale) +
  scale_fill_gradientn(name = "Density", colours = Colormap) +
  scale_y_discrete(expand = c(0.01, 0.9)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 2*ui_scale) +
  labs(title = "", x = glab("GSS"), y = "Region") +
  theme_paper(18)

# 选一个 annotation 的组成
anno_of_interest <- "Gaba-neuron"
if (!(anno_of_interest %in% show_a$annotation)) {
  anno_of_interest <- stats::na.omit(show_a$annotation)[1]
}
comp_df <- show_a |>
  filter(annotation == anno_of_interest) |>
  count(region) |>
  mutate(prop = n / sum(n)) |>
  arrange(desc(prop))
p_comp <- ggplot(comp_df, aes(x = fct_reorder(region, prop), y = prop)) +
  geom_col(width = 0.7, fill = blue, color = "black", linewidth = 1.5*ui_scale) +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),
            hjust = -0.05, size = 12*ui_scale, fontface = "bold") +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(title = paste0("Region Composition within Annotation = ", anno_of_interest),
       x = NULL, y = "Proportion") +
  theme_paper(18)

fig1_2 <- ggarrange(fig1, grid::nullGrob(), p_comp, ncol = 1, nrow = 3, heights = c(1, 0.05, 1))
figX   <- ggarrange(fig1, grid::nullGrob(), fig1_2, ncol = 3, nrow = 1, widths = c(1, 0.05, 1))

# region 比值（成人）
sum_region <- show_a |>
  filter(!is.na(region), !is.na(.data[[gcol("GSS")]])) |>
  group_by(region) |>
  summarise(
    n_lt = sum((.data[[gcol("GSS")]]) < 1, na.rm = TRUE),
    n_gt = sum((.data[[gcol("GSS")]]) > 1, na.rm = TRUE),
    n_eq = sum((.data[[gcol("GSS")]]) == 1, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    label = region,
    ratio = case_when(
      n_lt == 0 & n_gt > 0 ~ Inf,
      n_lt == 0 & n_gt == 0 ~ NA_real_,
      TRUE ~ n_gt / n_lt
    )
  ) |>
  arrange(desc(ratio)) |>
  mutate(label = factor(label, levels = rev(label)))

xmax <- max(sum_region$ratio[is.finite(sum_region$ratio)], na.rm = TRUE)
sum_region2 <- sum_region |>
  mutate(ratio_plot = ifelse(is.infinite(ratio), xmax * 1.05, ratio),
         ratio_lab  = ifelse(is.infinite(ratio), "∞", scales::number(ratio, accuracy = 0.01)))

p_region_db <- ggplot(sum_region, aes(y = label)) +
  geom_segment(aes(x = n_lt, xend = n_gt, yend = label), color = "grey70", linewidth = 5*ui_scale) +
  geom_point(aes(x = n_lt, color = paste(glab("GSS"), "< 1")), size = 15*ui_scale) +
  geom_point(aes(x = n_gt, color = paste(glab("GSS"), "> 1")), size = 15*ui_scale) +
  scale_color_manual(name = "Number of spots",
                     values = setNames(c(blue, red), c(paste(glab("GSS"), "< 1"), paste(glab("GSS"), "> 1")))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = scales::label_number(big.mark=",")) +
  labs(title = "", x = "Count", y = NULL) +
  theme_paper(18) +
  theme(
    axis.text.y  = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale),
    axis.text.x  = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(size = 1.5*ui_scale),
    plot.title   = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale, face = "bold"),
    axis.title.x = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale, face = "bold"),
    axis.title.y = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale, face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale),
    legend.text  = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(color = "black", size = 2.5*ui_scale)
  ) +
  guides(color = guide_legend(override.aes = list(size = 12*ui_scale)))

# 右侧柱图
xmax_lim_region <- xmax * 1.35
p_region_bar <- ggplot(sum_region2, aes(y = label, x = ratio_plot, fill = ratio_plot)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ratio_lab),
            hjust = -0.05, nudge_x = 0.02 * xmax_lim_region,
            size = 20*ui_scale, family = "Times New Roman") +
  scale_fill_gradient(low = "#04daf2ff", high = "#E74C3C", name = "Spots ratio") +
  scale_x_continuous(limits = c(0, xmax_lim_region),
                     expand = expansion(mult = c(0, 0.14))) +
  labs(x = paste0("ratio"), y = NULL) +
  theme_paper(18) +
  theme(
    axis.text.y  = element_blank(),
    axis.text.x  = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale),
    axis.ticks.x = element_line(size = 1.5*ui_scale),
    axis.title.x = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale, face = "bold"),
    axis.ticks.y = element_blank(),
    plot.margin  = margin(10, 14, 10, 4),
    legend.position = "none"
  ) +
  coord_cartesian(clip = "off")

fig3 <- p_region_db + p_region_bar + plot_layout(widths = c(3, 1))

# ==== 自动识别 gene 列名（保留原逻辑）====
if (!exists("gene")) {
  src <- if (exists("show_a")) show_a else if (exists("show")) show else NULL
  if (is.null(src)) stop("请先定义 show_a 或 show 数据框，或手动设置 gene。")
  gss_cols <- grep("_GSS$", names(src), value = TRUE)
  if (length(gss_cols) == 0) stop("未找到 *_GSS 列，请手动设置 gene。")
  gene <- sub("_GSS$", "", gss_cols[1])
}
gcol <- function(suffix) paste0(gene, "_", suffix)
glab <- function(suffix) paste(gene, suffix)
gss <- rlang::sym(gcol("GSS"))

# ==== 自适应画布高度 ====
n_cat        <- dplyr::n_distinct(show_a$annotation)
inch_per_cat <- 0.45
min_height   <- 8
plot_width   <- 12
plot_height  <- max(min_height, n_cat * inch_per_cat)
options(repr.plot.width = plot_width, repr.plot.height = plot_height)

wrap_lab <- function(x) stringr::str_wrap(x, width = 20)

# ==== 再做一轮 annotation 统计（与上面一致，保留）====
sum_annot <- show_a %>%
  group_by(annotation) %>%
  summarise(
    n_lt = sum((!!gss) < 1, na.rm = TRUE),
    n_gt = sum((!!gss) > 1, na.rm = TRUE),
    n_eq = sum((!!gss) == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(annotation),
    ratio = dplyr::case_when(
      n_lt == 0 & n_gt > 0 ~ Inf,
      n_lt == 0 & n_gt == 0 ~ NA_real_,
      TRUE ~ n_gt / n_lt
    )
  ) %>%
  arrange(desc(ratio)) %>%
  mutate(label = factor(label, levels = rev(label)))

xmax_annot <- max(sum_annot$ratio[is.finite(sum_annot$ratio)], na.rm = TRUE)
if (!is.finite(xmax_annot)) xmax_annot <- 1
sum_annot2 <- sum_annot %>%
  mutate(
    ratio_plot = ifelse(is.infinite(ratio), xmax_annot * 1.05, ratio),
    ratio_lab  = ifelse(is.infinite(ratio), "\u221E", scales::number(ratio, accuracy = 0.01))
  )
sum_annot2$ratio_lab <- enc2utf8(sum_annot2$ratio_lab)
sum_annot2$ratio_lab <- gsub('"', "", sum_annot2$ratio_lab, fixed = TRUE)

blue <- "#17BECF"; red <- "#E74C3C"
lab_lt <- paste(glab("GSS"), "< 1")
lab_gt <- paste(glab("GSS"), "> 1")

p_annot_db <- ggplot(sum_annot, aes(y = label)) +
  geom_segment(aes(x = n_lt, xend = n_gt, yend = label),
               color = "grey70", linewidth = 5*ui_scale) +
  geom_point(aes(x = n_lt, color = lab_lt), size = 15*ui_scale) +
  geom_point(aes(x = n_gt, color = lab_gt), size = 15*ui_scale) +
  scale_color_manual(
    name = "Number of spots",
    values = setNames(c(blue, red), c(lab_lt, lab_gt))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.15)),
    labels = scales::label_number(big.mark = ",")
  ) +
  scale_y_discrete(labels = wrap_lab) +
  labs(title = "", x = "Count", y = NULL) +
  theme_paper(18) +
  theme(
    axis.text.y  = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale),
    axis.text.x  = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(size = 1.5*ui_scale),
    plot.title   = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale, face = "bold"),
    axis.title.x = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale, face = "bold"),
    axis.title.y = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale, face = "bold"),
    legend.title = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale),
    legend.text  = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(color = "black", size = 2.5*ui_scale)
  ) +
  guides(color = guide_legend(override.aes = list(size = 12*ui_scale)))

p_annot_bar <- ggplot(sum_annot2, aes(y = label, x = ratio_plot, fill = ratio_plot)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ratio_lab), hjust = -0.1, size = 20*ui_scale, family = "Times New Roman") +
  scale_fill_gradient(low = "#04DAF2", high = "#E74C3C", name = "Spots ratio") +
  scale_x_continuous(
    limits = c(0, xmax_annot * 1.15),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(x = paste0("ratio"), y = NULL) +
  theme_paper(18) +
  theme(
    axis.text.y  = element_blank(),
    axis.text.x  = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale),
    axis.ticks.x = element_line(size = 1.5*ui_scale),
    axis.title.x = element_text(family = "Times New Roman", size = 25*scale_factor*ui_scale, face = "bold"),
    axis.ticks.y = element_blank(),
    plot.margin  = margin(10, 14, 10, 4),
    legend.position = "none"
  )

fig4 <- p_annot_db + p_annot_bar + patchwork::plot_layout(widths = c(3, 1))
fig4 <- fig4 + coord_cartesian(clip = "off") + theme(plot.margin = margin(10, 50, 10, 10))

fig3_4 <- ggarrange(fig3, grid::nullGrob(), fig4, ncol = 3, nrow = 1, widths = c(1, 0.05, 1))
figall  <- ggarrange(figall1, grid::nullGrob(), fig3_4, ncol = 1, nrow = 3, heights = c(1, 0.05, 1))

# ========= 最终组合 =========
if (!exists(combo_name, inherits = TRUE)) {
  stop(sprintf("未找到对象：%s。", combo_name))
}
figa_all <- ggarrange(
  get(combo_name),
  grid::nullGrob(),
  figall,
  ncol = 1,
  nrow = 3,
  heights = c(1, 0.05, 2)
)

# ========= 导出 =========
out_dir <- file.path(PROJECT_ROOT, "data", "Genes", gene)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_pdf <- file.path(out_dir, paste0("map_gene_", gene, ".pdf"))
out_png <- file.path(out_dir, paste0("map_gene_", gene, ".png"))

ggplot2::ggsave(filename = out_pdf, plot = figa_all,
                width = 16, height = 18, units = "in", limitsize = FALSE)
ggplot2::ggsave(filename = out_png, plot = figa_all,
                width = 16, height = 18, units = "in", dpi = 300, limitsize = FALSE)

cat("\n",
    "PDF: ", out_pdf, "\nPNG: ", out_png, "\n", sep = "")