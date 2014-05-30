require(RJSONIO)
require(ProjectTemplate)

clean_data <- function(){
  score_data <- fromJSON("data/MS1-matching-output 20131219_005_theoretical_ions_match_frags.csv_score.json")  
  whitelist <- sapply(score_data[[1]], function(col){class(col) != "list"})
  types <- sapply(score_data[[1]], function(col){class(col)})

  pruned_data <- lapply(score_data, function(row){row[whitelist]})
  ion_meta_data <- as.data.frame(do.call(rbind, pruned_data))

  ion_meta_data$MS1_Score <- as.numeric(ion_meta_data$MS1_Score)
  ion_meta_data$MS2_Score <- as.numeric(ion_meta_data$MS2_score)

  # Normalize between 0 and 1
  ion_meta_data$MS2_Score <- (ion_meta_data$MS2_Score - min(ion_meta_data$MS2_Score)) /
                            (max(ion_meta_data$MS2_Score) - min(ion_meta_data$MS2_Score))

  ion_meta_data$MS2_score <- unlist(ion_meta_data$MS2_score)
  
  ion_meta_data$vol <- as.numeric(unlist(ion_meta_data$vol))
  ion_meta_data$Seq_with_mod <- unlist(ion_meta_data$Seq_with_mod)

  ion_score_matrix <- data.frame(ion_meta_data$MS1_Score, ion_meta_data$MS2_Score, ion_meta_data$vol, ion_meta_data$Seq_with_mod)
  colnames(ion_score_matrix) <- c("ms1_score", "ms2_score", "volume", "sequence_with_mod")
  
  ion_score_matrix$calc_mass = as.numeric(unlist(ion_meta_data$Calc_mass))
  ion_score_matrix$obs_mass = as.numeric(unlist(ion_meta_data$Obs_Mass))
  
  set.seed(123)
  # K-Means has a random component. This may not yield the same results on different architectures
  clusters <- (kmeans(x=ion_score_matrix[c("ms1_score", "ms2_score")],centers=6))
  
  ion_score_matrix$cluster = as.factor(clusters$cluster)
  clusters$centers2 <- as.data.frame(clusters$centers)
  clusters$centers2 <- name_rows(clusters$centers2)
  clusters$centers2$volume <- 10e4
  
  save(ion_score_matrix, file="data/ion_score_matrix.RData")
  save(ion_meta_data, file="data/ion_meta_data.RData")
  cat(toJSON(ion_score_matrix, byrow = T, colNames = T, pretty = T), 
      file="graphs/ion_scores.json")
}

clean_data()