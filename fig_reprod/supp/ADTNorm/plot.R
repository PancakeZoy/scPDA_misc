require(ggplot2)
require(grid)
require(gridExtra)
require(scales)
require(dplyr)
library(tidyr)
require(tibble)

plot_f <- function(df, title, fill){
  min_val <- min(0, min(df$Expression))
  box = annotation_custom(grob = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 2)),
                          xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
  
  is_grouped <- is.vector(fill) && length(fill) == nrow(df)
  
  if (is_grouped) {
    df$fill_var <- fill
    p <- ggplot(df, aes(x = protein, y = Expression, fill = fill_var)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7),
               width = 0.7, color = 'grey40') +
      scale_fill_manual(values = c("ADTnorm" = "sienna1")) +
      guides(fill = guide_legend(title = "Method")) + theme_minimal() +
      theme(legend.position = c(0.02, 0.98),
            legend.justification = c("left", "top"))
  } else {
    p <- ggplot(df, aes(x = protein, y = Expression)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7),
               width = 0.7, fill = fill, color = 'grey40') +
      theme(legend.position = "none") + theme_minimal()
  }
  
  p <- p +
    labs(title = title, x = NULL, y = "NRC") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank()) +
    scale_y_continuous(labels = label_percent(scale = 100)) +
    coord_cartesian(ylim = c(-0.5, 1)) +
    box
  
  return(p)
}

unstain.adtnorm = readRDS('data/dsb_unstain_nrc.rds') %>% enframe(name = "protein", value = "Expression")
nonexpr.adtnorm = readRDS('data/dsb_nonexpr_nrc.rds')
nonexpr = data.frame(ADTnorm = nonexpr.adtnorm,
                     protein = names(nonexpr.adtnorm)) %>% 
  pivot_longer(cols = c(ADTnorm), 
               names_to = "Method", 
               values_to = "Expression")

p1 = plot_f(unstain.adtnorm, title='NRC derived from ADTnorm on unstained cells', fill='sienna1')
p2 = plot_f(nonexpr.adtnorm %>% enframe(name = "protein", value = "Expression"), 
            title='NRC derived from ADTnorm on non-expressive proteins', fill='sienna1')

label_plot <- function(p, label) {
  arrangeGrob(p, top = textGrob(label, x = unit(0, "npc"), y = unit(1, "npc"),
                                just = c("left", "top"),
                                gp = gpar(fontsize = 20, fontface = "bold")))
}

combined_plot <- arrangeGrob(
  label_plot(p1, "a"),
  label_plot(p2, "b"),
  ncol = 1
)

ggsave("ADTnorm.png", combined_plot, width = 12, height = 6, dpi = 200)
