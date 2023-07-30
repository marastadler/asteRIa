#' scatterplot_interactions
#'

library(scales)
library(ggrepel)

scatterplot_interactions <- function(coef_matrix, index_int, q = 12){
  coef_int <- coef_matrix[, index_int]
  coef_mat3 <- matrix(nrow = 3, ncol = sum(coef_matrix[(q + 1):nrow(coef_matrix),] != 0))
  rownames(coef_mat3) <- c("Mod. 1", "Mod. 2", "Int") 
  export_csv <- matrix(nrow = length(index_int) + 1
                       , ncol = 6)
  n <- 0
  col_names <- c()
  for(prot in index_int){
    for(x in 1:(q - 1)){
      for(y in (x + 1):q){
        mod1 <- rownames(coef_matrix)[x]
        mod2 <- rownames(coef_matrix)[y]
        if(paste0(mod1, ":", mod2) %in% rownames(coef_matrix)){
          mod1mod2 <- paste0(mod1, ":", mod2)
        }
        else{mod1mod2 <- paste0(mod2, ":", mod1)}
        
        
        if(coef_matrix[mod1mod2, prot] != 0){
          n <- n + 1
          
          coef_mat3[, n] <- c(coef_matrix[mod1, prot], coef_matrix[mod2, prot], 
                              coef_matrix[mod1mod2, prot])
          col_names[n] <- colnames(coef_matrix)[prot]
          
          export_csv[n, 1:6] <-  c(mod1, round(coef_matrix[mod1, prot], 2), mod2, 
                                   round(coef_matrix[mod2, prot], 2), mod1mod2,
                                   round(coef_matrix[mod1mod2, prot], 2))
        }
      }
      
    }
    
    
  }
  colnames(coef_mat3) <- col_names
  rownames(export_csv) <- col_names
  
  df <- as.data.frame(t(coef_mat3))
  df$label <- col_names
  
  plot_int <- ggplot(df, aes(x = df$`Mod. 1`, y = df$`Mod. 2`,
                             colour = df$Int#,
                             #label = rownames(df)
  )) +
    
    theme_minimal() +
    theme(axis.title=element_text(size=15)
    ) +
    theme(
      legend.title=element_text(size=15)) +
    geom_vline(xintercept = 0, color = "grey") +
    geom_hline(yintercept = 0, color = "grey") +
    geom_point(size = 2.5) +
    scale_colour_gradient2(
      low = muted("blue"),
      mid = "white",
      high = muted("red"),
      midpoint = 0
    ) +
    
    geom_text_repel(data = df,
                    aes(x = df$`Mod. 1`, y = df$`Mod. 2`,
                        colour = "black",
                        label = label
                    ), size = 2.5, color = "black", max.overlaps = 1000)
  return(plot_int)
}
