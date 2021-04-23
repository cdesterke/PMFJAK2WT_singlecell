
library(ggplot2)

#####  
df <- data.frame(x = data$mutations, z = data$clusters)
df <- as.data.frame(with(df, prop.table(table(x, z), margin = NULL)))



plot <- ggplot(data = df, aes(x = x, y = Freq, fill = z)) + 
  geom_bar(width = 0.9, position = "fill", stat = "identity") + 
  scale_fill_brewer(palette = "Paired") + 
  scale_y_continuous(expand = c(0.01, 0), labels = scales::percent_format()) + 
  xlab("mutations") + 
  ylab("Percent") + 
  labs(fill = "treatment") + 
  theme_classic(base_size = 12, base_family = "sans") + 
  theme(legend.position = "right") + coord_flip()
plot
