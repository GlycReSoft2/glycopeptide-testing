library('ProjectTemplate')
load.project()


scale_valid_plot <- ggplot(ion_meta_data, aes(x = MS2_score, y = MS2_Score)) + geom_point() + 
  labs(y = "Scaled MS2 Score", x = "Unscaled MS2 Score", title = "No distortion by 0-1 Scaling") 
print(scale_valid_plot)
ggsave("graphs/no_scale_distortion_plot.png")


correl_plot <- 
  ggplot() + 
  geom_point(data = ion_score_matrix, 
             mapping = aes(x = ms1_score, 
                           y = ms2_score, size = volume, 
                           color = sequence_with_mod, shape = cluster)) + 
  guides(colour = guide_legend(ncol = 2,
                               byrow = TRUE, 
                               label.theme = element_text(size=8, angle=0)))
print(correl_plot)
ggsave("graphs/correl_plot.png", width=23, height=15)
