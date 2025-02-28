library(ggplot2)

segment_name = "chr14.24162155.24162244"
gene_name = "ENSG00000213928"
clinical_parameter = "PA_systolic"

atac_vector = c(0.945924 , 0 , 2.08935 , 0.925572 , 0.666292 , 0.55749 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0.991903 , 1.31306 , 0 , 1.26098 , 1.27345 , 0.767462 , 2.22566 , 0.919898 , 0.928498 , 1.66963 , 0 , 0)

clinical_vector = c(22 , 85 , 27 , 43 , 32 , 62 , 81 , 107 , 109 , 73 , NA , 93 , 76 , 57 , 52 , 37 , 69 , 46 , 41 , 75 , 40 , 33 , 51 , 66 , 104 , 81)


ggplot(data = data.frame(x = atac_vector, y = clinical_vector), aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = paste("Correlation between ATAC-seq signal and", clinical_parameter, "for", gene_name, "in", segment_name),
       x = "ATAC-seq signal",
       y = clinical_parameter) +
  theme_minimal()





segment_name = "chr9.78235502.78235541"
gene_name = "ENSG00000148019"
clinical_parameter = "PA_systolic"

atac_vector = c(3.70145 , 0 , 1.3929 , 0 , 2.72574 , 1.14032 , 0 , 0 , 0 , 0.538901 , 0 , 0.63404 , 0 , 1.70598 , 0 , 0.868933 , 2.0905 , 3.3877 , 2.49153 , 0 , 2.68271 , 3.5996 , 1.86529 , 1.12139 , 0 , 0)

clinical_vector = c(22 , 85 , 27 , 43 , 32 , 62 , 81 , 107 , 109 , 73 , NA , 93 , 76 , 57 , 52 , 37 , 69 , 46 , 41 , 75 , 40 , 33 , 51 , 66 , 104 , 81)

ggplot(data = data.frame(x = atac_vector, y = clinical_vector), aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = paste("Correlation between ATAC-seq signal and", clinical_parameter, "for", gene_name, "in", segment_name),
       x = "ATAC-seq signal",
       y = clinical_parameter) +
  theme_minimal()





segment_name = "chr14.24189949.24190798"
gene_name = "ENSG00000196497"
clinical_parameter = "PA_systolic"

atac_vector = c(0.870927058823529 , 0 , 0.614515070588235 , 0.300538470588235 , 1.04540455764706 , 0.673459505882353 , 0.249601058823529 , 0 , 0 , 0.634000411764706 , 0.432818882352941 , 0.229746176470588 , 0.413867617647059 , 0.505772294117647 , 0.456630964705882 , 0.815774117647059 , 0.491881011764706 , 0.350726234117647 , 1.16662211764706 , 0.451448191764706 , 0.631227091764706 , 0.842729411764706 , 0.962049214941176 , 0.5303518 , 0.522606223529412 , 0)

clinical_vector = c(22 , 85 , 27 , 43 , 32 , 62 , 81 , 107 , 109 , 73 , NA , 93 , 76 , 57 , 52 , 37 , 69 , 46 , 41 , 75 , 40 , 33 , 51 , 66 , 104 , 81)

ggplot(data = data.frame(x = atac_vector, y = clinical_vector), aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = paste("Correlation between ATAC-seq signal and", clinical_parameter, "for", gene_name, "in", segment_name),
       x = "ATAC-seq signal",
       y = clinical_parameter) +
  theme_minimal()





segment_name = "chr11.64758404.64758493"
gene_name = "ENSG00000068976"
clinical_parameter = "PA_systolic"

atac_vector = c(0.870927058823529 , 0 , 0.614515070588235 , 0.300538470588235 , 1.04540455764706 , 0.673459505882353 , 0.249601058823529 , 0 , 0 , 0.634000411764706 , 0.432818882352941 , 0.229746176470588 , 0.413867617647059 , 0.505772294117647 , 0.456630964705882 , 0.815774117647059 , 0.491881011764706 , 0.350726234117647 , 1.16662211764706 , 0.451448191764706 , 0.631227091764706 , 0.842729411764706 , 0.962049214941176 , 0.5303518 , 0.522606223529412 , 0)

clinical_vector = c(22 , 85 , 27 , 43 , 32 , 62 , 81 , 107 , 109 , 73 , NA , 93 , 76 , 57 , 52 , 37 , 69 , 46 , 41 , 75 , 40 , 33 , 51 , 66 , 104 , 81)

ggplot(data = data.frame(x = atac_vector, y = clinical_vector), aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = paste("Correlation between ATAC-seq signal and", clinical_parameter, "for", gene_name, "in", segment_name),
       x = "ATAC-seq signal",
       y = clinical_parameter) +
  theme_minimal()