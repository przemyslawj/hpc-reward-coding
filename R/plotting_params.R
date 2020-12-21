library(ggplot2)
library(cowplot)

px2pt = 1/(ggplot2::.pt*72.27/96)
gtheme = theme_minimal() + theme_cowplot() + 
  theme(text = element_text(size=8, family = 'Arial'), 
        panel.background = element_rect(fill = 'white'),
        line = element_line(size = px2pt),
        axis.line.x = element_line(size = 0.8 * px2pt, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.8 * px2pt, linetype = "solid", colour = "black"),
        axis.text.x=element_text(size=8, colour = "black"),
        axis.text.y=element_text(size=8, colour = "black"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        plot.title = element_text(size=8))
