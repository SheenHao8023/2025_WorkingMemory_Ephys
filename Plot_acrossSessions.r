# -------- 1. 加载必要包 --------
pkg_list <- c("ggplot2", "readxl", "tidyr", "dplyr", "ggprism", "patchwork")
sapply(pkg_list, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
})

# -------- 2. 读取 Excel 数据 --------
path <- "C:/Users/XinHao/Desktop/"
file_path <- file.path(path, 'data_acrossSessions.xlsx')
df <- readxl::read_excel(file_path, sheet = 1)

# -------- 3. 自定义主题与配色 --------
theme_custom <- theme_prism() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "top"
  )

mycolors <- c("fitting" = "#80d6ff", "cdproj" = "#f47c7c")

# -------- 4. func: 2×2 bar plot --------
plot_2x2 <- function(data, cols, ylab) {
  df_long <- data %>%
    select(all_of(cols)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    mutate(
      ROI = ifelse(grepl("SC", variable), "SC", "SNr"),
      method = ifelse(grepl("cdproj", variable), "cdproj", "fitting")
    )
  
  summary_df <- df_long %>%
    group_by(ROI, method) %>%
    summarise(mean = mean(value, na.rm = TRUE),
              sd = sd(value, na.rm = TRUE),
              .groups = "drop")
  
  ggplot(summary_df, aes(x = ROI, y = mean, fill = method)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.5, color = "black") +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2,
                  position = position_dodge(0.8)) +
    geom_jitter(data = df_long, aes(x = ROI, y = value), color = 'black',
                size = 2, shape = 16, stroke = 0.5, show.legend = FALSE,
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7)) +
    scale_fill_manual(values = mycolors) +
    labs(x = NULL, y = ylab) +
    theme_custom
}

# -------- 5. func: Angle bar plot --------
plot_singlebar <- function(data, col, xlab, ylab) {
  df <- data %>% select(all_of(col)) %>% rename(value = !!col)
  mean_val <- mean(df$value, na.rm = TRUE)
  sd_val <- sd(df$value, na.rm = TRUE)
  
  ggplot(df, aes(x = "SC-SNr", y = value)) +
    geom_bar(stat = "summary", fun = mean, width = 0.3, fill = NA, color = "black") +
    geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val), width = 0.1) +
    geom_jitter(width = 0.05, size = 2, alpha = 1, color = 'black') +
    labs(x = NULL, x = xlab, y = ylab) +
    theme_custom +
    theme(axis.text.x = element_blank())
}

# -------- 6. func: scatter plot --------
plot_scatter <- function(data, x_expr, y, xlab, ylab = "Behavior") {
  ggplot(data, aes(x = {{ x_expr }}, y = .data[[ylab]])) +
    geom_point(color = "black", size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    labs(x = xlab, y = ylab) +
    theme_custom
}

# -------- 7. draw --------
p1 <- plot_2x2(df, c("r_SC_fitting", "r_SC_cdproj", "r_SNr_fitting", "r_SNr_cdproj"), ylab = "Pearson r")
p2 <- plot_2x2(df, c("R2_SC_fitting", "R2_SC_cdproj", "R2_SNr_fitting", "R2_SNr_cdproj"), ylab = expression(R^2))
p3 <- plot_singlebar(df, "Angle", xlab = 'SC-SNr', ylab = "Angle(°)")

p4 <- plot_2x2(df, c("r_SC_fitting_noise", "r_SC_cdproj_noise", "r_SNr_fitting_noise", "r_SNr_cdproj_noise"), ylab = "Pearson r (Noise)")
p5 <- plot_2x2(df, c("R2_SC_fitting_noise", "R2_SC_cdproj_noise", "R2_SNr_fitting_noise", "R2_SNr_cdproj_noise"), ylab = expression(R^2 * " (Noise)"))
p6 <- plot_singlebar(df, "Angle_noise", xlab = 'SC-SNr', ylab = "Angle_Noise(°)")

p7 <- plot_scatter(df, r_SC_fitting - r_SC_cdproj, Behavior, "diff Pearson r SC")
p8 <- plot_scatter(df, r_SNr_fitting - r_SNr_cdproj, Behavior, "diff Pearson r SNr")
p9 <- plot_scatter(df, r_SC_fitting_noise - r_SC_cdproj_noise, Behavior, "diff Pearson r SC_noise")
p10 <- plot_scatter(df, r_SNr_fitting_noise - r_SNr_cdproj_noise, Behavior, "diff Pearson r SNr_noise")

# -------- 8. 拼图组合 --------
final_bar_plot <- (p1 | p2 | p3) / (p4 | p5 | p6)
final_scatter_plot <- (p7 | p8) / (p9 | p10)

# -------- 9. 保存图像 --------
ggsave(file.path(path, 'plot_bar_acrossSessions.png'), final_bar_plot, width = 14, height = 8, dpi = 300)
ggsave(file.path(path, "plot_scatter_acrossSessions.png"), final_scatter_plot, width = 12, height = 10, dpi = 300)
