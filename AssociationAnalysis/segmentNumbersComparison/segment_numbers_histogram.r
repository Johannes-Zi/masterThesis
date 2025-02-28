# Import csv file with segment numbers per clinical parameter before filtering
input_path_1 <- "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/plotRaw04Corrs/runs/11_01_25/num_genes_above_threshold_per_clinical_parameter.csv"
df_segment_numbers_before_filtering <- read.csv(input_path_1, header = TRUE)

print(head(df_segment_numbers_before_filtering))

# Import csv file with segment numbers per clinical parameter after filtering
input_path_2 <- "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/CorrelationVisualizationFiltered/visulizations/v1/number_of_entries_per_param.tsv"
df_segment_numbers_after_filtering <- read.delim(input_path_2, header = TRUE, sep = "\t")

print(head(df_segment_numbers_after_filtering))

# Rename columns to prepare alignment
colnames(df_segment_numbers_before_filtering) <- c("clinical_parameter", "segment_number_before_filtering")
colnames(df_segment_numbers_after_filtering) <- c("clinical_parameter", "segment_number_after_filtering")

#print(head(df_segment_numbers_before_filtering))
#print(head(df_segment_numbers_after_filtering))

df_combined <- merge(df_segment_numbers_before_filtering, df_segment_numbers_after_filtering, by = "clinical_parameter", all = TRUE)

df_combined <- df_combined[!df_combined$clinical_parameter == "SMW", ]

rownames(df_combined) <- 1:nrow(df_combined)

#print(head(df_combined))


library(ggplot2)

#  dataframe with the data for the plot
df_plot <- data.frame(clinical_parameter = df_combined$clinical_parameter,
                      segment_number_before_filtering = df_combined$segment_number_before_filtering,
                      segment_number_after_filtering = df_combined$segment_number_after_filtering)

# df_plot$clinical_parameter <- factor(df_plot$clinical_parameter, 
#     levels = c("age",
#                "AT",
#                "AT_ET", 
#                "AVDO2", 
#                "BNP", 
#                "Cardiac.Index", 
#                "CVP", 
#                "DLCO", 
#                "eGFR", 
#                "FEV1", 
#                "heart.rate", 
#                "height", 
#                "HZV_Fick", 
#                "HZV_Thermodil", 
#                "Kreatinin", 
#                "mPAP", 
#                "NYHA", 
#                "PA_diastolic", 
#                "PA_systolic", 
#                "Paradoxe_Septumbewegung", 
#                "PAWP", 
#                "Perikarderguss", 
#                "PVR", 
#                "RA_area", 
#                "Rrsys", 
#                "RVEDD", 
#                "S", 
#                "sPAP.excl..ZVD", 
#                "TAPSE", 
#                "Tiffeneau.Index.FEV1.VC", 
#                "VCI_diameter", 
#                "ven_SO2", 
#                "weight", 
#                "ZVD"))

#print("### before rename")
#print(head(df_plot))

# Rename all clinical parameters to be more readable
df_plot$clinical_parameter <- gsub("age", "age", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("\\bAT\\b", "pulmonary acceleration time", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("AT_ET", "(pulmonary acceleration time) / (right ventricular ejection time)", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("AVDO2", "arteriovenous oxygen difference", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("BNP", "BNP heart failure marker", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("Cardiac.Index", "caridac performance index", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("CVP", "central venous pressure", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("DLCO", "lung CO diffusion capacity", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("eGFR", "glomerular filtration rate", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("heart.rate", "heart rate", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("height", "height", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("HZV_Fick", "cardiac output - Fick principle", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("HZV_Thermodil", "cardiac output - thermodilution", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("Kreatinin", "kreatinin", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("mPAP", "mean pulmonary artery pressure", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("NYHA", "NYHA heart failure classification", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("PA_diastolic", "diastolic pulmonary artery pressure", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("PA_systolic", "systolic pulmonary artery pressure", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("Paradoxe_Septumbewegung", "ventricular septum movement", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("PAWP", "pulmonary capillary occlusion pressure", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("Perikarderguss", "pericardium fluid accumulation", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("PVR", "pulmonary vascular resistance", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("RA_area", "right atrium size", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("Rrsys", "systolic blood pressure at rest", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("RVEDD", "right ventricle size at diastole", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("\\bS\\b", "TAPSV impaired right ventricular contraction indicator", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("sPAP.excl..ZVD", "estimated systolic pulmonary artery pressure", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("TAPSE", "TAPSE impaired right ventricular contraction marker", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("Tiffeneau.Index.FEV1.VC", "airflow obstruction in lung - Tiffeneau", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("ZVD", "estimated central venous pressure", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("VCI_diameter", "inferior vena cava diameter", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("ven_SO2", "venous oxygen saturation", df_plot$clinical_parameter)
df_plot$clinical_parameter <- gsub("\\bFEV1\\b", "one second forced expiratory volume", df_plot$clinical_parameter)

#print(head(df_plot))
# Add colum to df wich holds the partion of after filtering to before filteing
df_plot$partion <- df_plot$segment_number_after_filtering / df_plot$segment_number_before_filtering
df_plot$partion_percent <- df_plot$partion * 100

#print(head(df_plot))

df_plot_sorted <- df_plot[order(df_plot$partion, decreasing = FALSE), ]

rownames(df_plot_sorted) <- 1:nrow(df_plot_sorted)

# Update the factor levels of clinical_parameter based on the sorted order
df_plot_sorted$clinical_parameter <- factor(df_plot_sorted$clinical_parameter, levels = df_plot_sorted$clinical_parameter)

#print(head(df_plot_sorted))

ggplot(df_plot_sorted, aes(y = clinical_parameter)) +
  geom_col(aes(x = segment_number_before_filtering, fill = "Before filtering"), position = "dodge") +
  geom_col(aes(x = segment_number_after_filtering * 10, fill = "After filtering"), position = "dodge") +
  geom_text(aes(x = segment_number_after_filtering * 10, label = paste0(round(partion_percent, 2), "%")), 
            position = position_dodge(width = 0.9), hjust = -0.1, size = 5, color = "#000000") +
  scale_fill_manual(values = c("Before filtering" = "#e7b800", "After filtering" = "#009e73")) +
  scale_x_continuous(
    name = "segments before filtering",
    sec.axis = sec_axis(~ . / 10, name = "segments after filtering")
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(hjust = 1),
        legend.position = "bottom",
        legend.box = "vertical",
        axis.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 16, colour = "#000000"),
        legend.title = element_text(size = 16, colour = "#000000"),
        axis.text = element_text(size = 16, colour = "black")) +
  labs(title = "",
       y = "", fill = "") +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  geom_point(aes(x = Inf, y = Inf, color = "% Remaining segment set size after filter"), size = 0, show.legend = TRUE) +
  scale_color_manual(name = "", values = c("% Remaining segment set size after filter" = "black"))


output_dir <- "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/segmentNumbersComparison/"

output_path <- paste0(output_dir, "segment_numbers_histogram.png")
ggsave(output_path, width = 12, height = 12, dpi = 300)