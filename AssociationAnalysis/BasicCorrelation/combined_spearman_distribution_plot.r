pa_dia_path <- "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/BasicCorrelation/combined_corr_cv_pear_04_thres/PA_diastolic/PA_diastolic_spearman_correlations.tsv"
pa_sys_path <- "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/BasicCorrelation/combined_corr_cv_pear_04_thres/PA_systolic/PA_systolic_spearman_correlations.tsv"
df_dia <- read.table(pa_dia_path, header=TRUE, sep=";")
df_sys <- read.table(pa_sys_path, header=TRUE, sep=";")

df_dia$Type <- "diastolic"
df_sys$Type <- "systolic"
df_combined <- rbind(df_dia, df_sys)

print(head(df_combined))

#density plot for the Spearman correlation values
library(ggplot2)
ggplot(df_combined, aes(x=spearman_corr, fill=Type)) +
    geom_density(alpha=0.6) +
    labs(title="A",
       x="Spearman correlation value", y="Density", fill="Pulmonary artery pressure") +
    theme_bw() +
    coord_cartesian(xlim=c(-0.8, 0.8)) +
    theme(legend.position="bottom",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) +
    scale_fill_manual(values=c("diastolic"="#00aeba", "systolic"="#e7b800"))

codefile_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

ggsave(file.path(codefile_dir, "combined_spearman_distribution_plot.png"), width=6, height=6)



df_combined_filtered <- df_combined %>% filter(p_value < 0.05)

ggplot(df_combined_filtered, aes(x=spearman_corr, fill=Type)) +
    geom_density(alpha=0.6) +
    labs(title="B",
       x="Spearman correlation value \n (p-value < 0.05)", y="Density", fill="Pulmonary artery pressure") +
    theme_bw() +
    theme(legend.position="bottom",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
    ) +
    coord_cartesian(xlim=c(-0.8, 0.8)) +
    scale_fill_manual(values=c("diastolic"="#00aeba", "systolic"="#e7b800"))

ggsave(file.path(codefile_dir, "combined_spearman_distribution_plot_filtered.png"), width=6, height=6.5)
