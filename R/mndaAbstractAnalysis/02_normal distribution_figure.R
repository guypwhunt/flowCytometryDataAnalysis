edema <- round(rnorm(10000, mean=50 , sd=20))
neutrophils <- round(rnorm(10000, mean=100 , sd=25))
monocytes <- round(rnorm(10000, mean=150 , sd=30))

df <- data.frame(name = rep("edema", length(edema)), value = edema)
df <- rbind(df, data.frame(name = rep("neutrophils", length(neutrophils)), value = neutrophils))
df <- rbind(df, data.frame(name = rep("monocytes", length(monocytes)), value = monocytes))

ggplot(df, aes(x=value, group=name, color=name, fill=name, label=name)) +
  geom_density(alpha=0.4, show.legend = FALSE) +
  theme(axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        #axis.ticks.y=element_blank()
        ) +
  labs(y= "", x = "")

