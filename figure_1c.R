library(jsonlite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(R.matlab)
library(scales)
library(grid)
library(cowplot)

use <- "agcl"  
lead_to_plot <- 3

if (use == "hg") {
    data_sensor <- readMat("hg_p5.mat")
} else {
    data_sensor <- readMat("agcl_p5.mat")
}

artifacts <- list(
  list(start = 10, end = 14, label = "Motion Artifact 1", position = "left"),
  list(start = 17, end = 21, label = "Motion Artifact 2", position = "right")
)


if (use == "hg"){
    ecg_ch <- data_sensor$hg[,lead_to_plot] 
    ecg_filtered_ch <- data_sensor$hg.all.filters[,lead_to_plot]  
} else {
    ecg_ch <- data_sensor$agcl[,lead_to_plot] 
    ecg_filtered_ch <- data_sensor$agcl.all.filters[,lead_to_plot]  
}

if (use == "hg") {
  sensor_raw_label <- "PPHG Raw"
  sensor_processed_label <- "PPHG Processed"
  sensor_color_values <- c("#000000", "#0066CC")  
} else {
  sensor_raw_label <- "AgCl Raw"
  sensor_processed_label <- "AgCl Processed"
  sensor_color_values <- c("#000000", "#d9020d")  
}

time <- data_sensor$t.vec

ecg_data <- data.frame(
  Time = rep(time, 2),
  Amplitude = c(ecg_ch, ecg_filtered_ch),
  Sensor = factor(rep(c(sensor_raw_label, sensor_processed_label), each = length(time)))
)

y_min <- min(ecg_data$Amplitude)
y_max <- max(ecg_data$Amplitude)
y_range <- y_max - y_min
plot_top <- y_max + 0.01 * y_range 

color_values <- c(sensor_color_values)
names(color_values) <- c(sensor_raw_label, sensor_processed_label)


main_plot <- ggplot(ecg_data, aes(x = Time, y = Amplitude, color = Sensor)) +
  geom_line(size = 0.7) +
  scale_color_manual(values = color_values) +
  coord_cartesian(ylim = c(y_min - 0.12 * y_range, plot_top + 0.2 * y_range)) +
  labs(
    title = "",
    x = "Time (s)",
    y = expression("Amplitude ("*mu*"V)"),
    color = ""
  ) +
  theme_minimal_grid(font_size = 10) + 
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 30),
    legend.key.width = unit(2, "cm"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(5, 10, 5, 10) 
  ) 

for (artifact in artifacts) {
  main_plot <- main_plot + 
    annotate("rect", xmin = artifact$start, xmax = artifact$end, 
             ymin = -Inf, ymax = Inf, fill = "#0d84ac", alpha = 0.2)
}

# Function to create bubble plot
draw_bubble_panel <- function(plot) {
  ggdraw() +
    draw_grob(
      circleGrob(
        x = 0.5, y = 0.5, r = 1,
        gp = gpar(fill = "white", col = "gray50", lwd = 1.8)
      )
    ) +
    draw_plot(plot, x = 0.1, y = 0.1, width = 0.8, height = 0.8)
}


bubble_plots <- list()
x_mains <- numeric(length(artifacts))
y_mains <- numeric(length(artifacts))


for (i in 1:length(artifacts)) {
  artifact <- artifacts[[i]]
  

  zoom_plot <- ggplot(
    subset(ecg_data, Time >= artifact$start & Time <= artifact$end), 
    aes(x = Time, y = Amplitude, color = Sensor)
  ) +
    geom_line(size = 1.2) +
    scale_color_manual(values = color_values) +
    theme_void() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  bubble_plots[[i]] <- draw_bubble_panel(zoom_plot)
  
  x_mains[i] <- mean(c(artifact$start, artifact$end))
  y_mains[i] <- mean(subset(ecg_data, 
                          Time >= artifact$start & 
                          Time <= artifact$end)$Amplitude)
}

bubble_y <- plot_top + 0.1 * y_range
bubble_spacing <- 8  

final_plot <- main_plot

for (i in 1:length(artifacts)) {
  artifact <- artifacts[[i]]
  
  # Add bubble
  final_plot <- final_plot + 
    annotation_custom(
      grob = ggplotGrob(bubble_plots[[i]]),
      xmin = x_mains[i] - 3.5, xmax = x_mains[i] + 3.5,
      ymin = bubble_y - 200, ymax = bubble_y + 200
    )
  
  if (artifact$position == "left") {
    text_x <- x_mains[i] - 3
    hjust_val <- 1
  } else {
    text_x <- x_mains[i] + 3
    hjust_val <- 0
  }
  
  final_plot <- final_plot +
    annotate(
      "text", 
      x = text_x,
      y = bubble_y, 
      label = artifact$label, 
      hjust = hjust_val,
      vjust = 3,
      size = 6, 
      fontface = "bold"
    )
  
  final_plot <- final_plot +
    geom_segment(
      data = data.frame(
        x = x_mains[i],
        y = y_mains[i],
        xend = x_mains[i],
        yend = bubble_y - 200
      ),
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      arrow = arrow(type = "closed", length = unit(0.25, "cm")),
      color = "black", linetype = "dashed", size = 0.8
    )
}

print(final_plot)

ggsave(paste0("R_figures/figure_1c/",use,"_lead", lead_to_plot, "_with_enlarged_bubbles.png"), 
       final_plot, width = 12, height = 8, dpi = 300, bg = "white")