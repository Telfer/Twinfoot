## Figures for paper

# =============================================================================

## dynamic trials
# import dynamic SPM results
dy_results <- read.csv("Z:/Projects/049 TwinFoot/DATA/Emed/dy_results.csv", 
                       header = FALSE)
mz_results <- read.csv("Z:/Projects/049 TwinFoot/DATA/Emed/mz_results.csv", 
                       header = FALSE)

# overall
dy_means <- colMeans(abs(dy_results))
mz_means <- colMeans(abs(mz_results))
dy_max <- apply(abs(dy_results), 2, max)
mz_max <- apply(abs(mz_results), 2, max)

# Dimensions
area_length <- 54
area_width <- 21

# plot images
dy_mean_image <- plot_spm_pp(dy_means, active_sensors = active_sensors, 
                             area_length = area_length, 
                             area_width = area_width)
mz_mean_image <- plot_spm_pp(mz_means, active_sensors = active_sensors, 
                             area_length = area_length, 
                             area_width = area_width)
dy_max_image <- plot_spm_pp(dy_max, active_sensors = active_sensors, 
                            area_length = area_length, 
                            area_width = area_width)
mz_max_image <- plot_spm_pp(mz_max, active_sensors = active_sensors, 
                            area_length = area_length, 
                            area_width = area_width)

ggarrange(dy_mean_image, mz_mean_image, ncol = 2)
ggarrange(dy_max_image, mz_max_image, ncol = 2)


# split to every 20% of stance
# dy results
dy_20 <- colMeans(abs(dy_results[1:20, ]))
dy_40 <- colMeans(abs(dy_results[21:40, ]))
dy_60 <- colMeans(abs(dy_results[41:60, ]))
dy_80 <- colMeans(abs(dy_results[61:80, ]))
dy_100 <- colMeans(abs(dy_results[81:101, ]))

# mz results
mz_20 <- colMeans(abs(mz_results[1:20, ]))
mz_40 <- colMeans(abs(mz_results[21:40, ]))
mz_60 <- colMeans(abs(mz_results[41:60, ]))
mz_80 <- colMeans(abs(mz_results[61:80, ]))
mz_100 <- colMeans(abs(mz_results[81:101, ]))

# names
time_steps_tw <- c("dy_20", "dy_40", "dy_60", "dy_80", "dy_100", 
                   "mz_20", "mz_40", "mz_60", "mz_80", "mz_100") 

# produce images
dynamic_images <- list()
for (i in seq_along(time_steps_tw)) {
  dynamic_images[[i]] <- plot_spm_pp(get(time_steps_tw[i]), active_sensors, 
                                     area_length = area_length, 
                                     area_width = area_width)
}

ggsave("dynamic_images_dy20.png", dynamic_images[[1]], dpi = 600)
ggsave("dynamic_images_dy40.png", dynamic_images[[2]], dpi = 600)
ggsave("dynamic_images_dy60.png", dynamic_images[[3]], dpi = 600)
ggsave("dynamic_images_dy80.png", dynamic_images[[4]], dpi = 600)
ggsave("dynamic_images_dy100.png", dynamic_images[[5]], dpi = 600)
ggsave("dynamic_images_mz20.png", dynamic_images[[6]], dpi = 600)
ggsave("dynamic_images_mz40.png", dynamic_images[[7]], dpi = 600)
ggsave("dynamic_images_mz60.png", dynamic_images[[8]], dpi = 600)
ggsave("dynamic_images_mz80.png", dynamic_images[[9]], dpi = 600)
ggsave("dynamic_images_mz100.png", dynamic_images[[10]], dpi = 600)


# =============================================================================

# figure 2 - average z-stat throughout stance
summary_df <- data.frame(Mz = rowMeans(abs(mz_results)), 
                         Dz = rowMeans(abs(dy_results)))
summary_df <- gather(summary_df, "Twin_type", "z_stat", 1:2)
summary_df <- cbind(summary_df, Time = rep(1:101, 2))
g <- ggplot(summary_df, aes(x = Time, y = z_stat, colour = Twin_type))
g <- g + geom_line()
g <- g + xlab("% Stance") + ylab("z-statistic")
g <- g + guides(colour = (guide_legend(title = "Twin type")))
g <- g + theme_minimal()
g

ggsave("Figure 2.png", g, dpi = 600, height = 4, width = 8)


# =============================================================================

## make animation of SPM results
for (i in 1:101) {
  filename <- paste0("C:/Users/telfe/Dropbox/My_Projects/Twinfoot/paper", 
                     "/animation", "/twin_ani", i, ".png")
  g1 <- plot_spm_pp(abs(dy_results[i, ]), active_sensors, area_length, 
                    area_width)
  g1 <- g1 + ggtitle("DZ")
  g2 <- plot_spm_pp(abs(mz_results[i, ]), active_sensors, area_length, 
                    area_width)
  g2 <- g2 + ggtitle("MZ")
  g <- ggarrange(g1, g2, ncol = 2)
  ggsave(filename, g, dpi = 300, width = 3, height = 5, bg = "transparent")
}


# =============================================================================

## static trials
# import static SPM results
dy_static_results <- read.csv("Z:/Projects/049 TwinFoot/DATA/Emed/dy_static_results.csv", 
                       header = FALSE)
mz_static_results <- read.csv("Z:/Projects/049 TwinFoot/DATA/Emed/mz_static_results.csv", 
                       header = FALSE)

# overall
dy_static_means <- colMeans(abs(dy_static_results))
mz_static_means <- colMeans(abs(mz_static_results))

# plot results
# Dimensions
area_length <- dim(static_footprints_registered[[1]])[1]
area_width <- dim(static_footprints_registered[[1]])[2]
dy_static_mean_image <- plot_spm_pp(dy_static_means, 
                                    active_sensors = active_sensors, 
                                    area_length = area_length, 
                                    area_width = area_width)
mz_static_mean_image <- plot_spm_pp(mz_static_means, 
                                    active_sensors = active_sensors, 
                                    area_length = area_length, 
                                    area_width = area_width)
dy_static_mean_image
mz_static_mean_image



