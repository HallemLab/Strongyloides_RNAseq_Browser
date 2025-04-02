theme_Publication <- function(base_size=13, base_family="sans") {
        
        (theme_foundation(base_size=base_size, base_family=base_family)
                + theme(plot.title = element_text(face = "bold",
                                                  size = 15, hjust = 0),
                        plot.title.position = "plot",
                        plot.caption.position = "plot",
                        text = element_text(),
                        panel.background = element_rect(colour = NA),
                        plot.background = element_rect(colour = NA),
                        panel.border = element_rect(colour = NA),
                        axis.title = element_text(face = "bold",size = rel(1)),
                        axis.title.y = element_text(angle=90,vjust =2),
                        axis.title.x = element_text(vjust = -0.2),
                        axis.text = element_text(size = rel(1)), 
                        axis.line = element_line(colour="black", size = 0.5, lineend = "square"),
                        axis.ticks = element_line(size = 0.5),
                        axis.ticks.length = unit(1, "mm"),
                        legend.key = element_rect(colour = NA),
                        legend.position = "right",
                        legend.direction = "vertical",
                        legend.key.size= unit(1, "cm"),
                        legend.margin = margin(0,0,0,0, "cm"),
                        legend.spacing.x = unit(0.2,'cm'),
                        legend.title = element_text(face="bold"),
                        panel.grid = element_line(colour = NA),
                        plot.margin=margin(10,5,5,5,"mm"),
                        strip.background=element_rect(colour="black",
                                                      fill="grey90",
                                                      size = .5),
                        strip.text = element_text(face="bold")
                ))
        
}