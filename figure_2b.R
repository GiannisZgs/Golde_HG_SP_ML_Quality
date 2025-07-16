library(jsonlite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(grid)
library(cowplot)
library(ggrepel)
library(RColorBrewer)

data_path <- "/home/giannis/Documents/ECG HG paper/results_data/participant_id_results.json"
data <- fromJSON(txt = data_path, simplifyVector = FALSE)
fs <- 200
sensors <- c("AgCl", "HG")

# Create output directory
output_dir <- paste0("R_figures/figure_2b/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

selected_participants_profiles <- c("p1","p5", "p10", "p22", "p39") # Pick 5 participants
channels <- c("ch1", "ch2", "ch3")

# Selected sensor-channel for profile visualization
profile_sensor <- "AgCl"  # Will be displayed as PPHG
profile_channel <- "ch3"  # Lead 3

# Update the bubble panel function to create MUCH larger bubbles with clearer waveforms
draw_bubble_panel <- function(plot) {
  ggdraw() +
    draw_grob(
      circleGrob(
        x = 0.5, y = 0.5, r = 0.5,
        gp = gpar(fill = "white", col = "gray50", lwd = 2.5)  # Even thicker border
      )
    ) +
    draw_plot(plot, x = 0.1, y = 0.1, width = 0.8, height = 0.8)  # Larger content area
}

create_tsne_plots <- function(data, sensor, channel, save_individual = TRUE) {
  # Extract t-SNE data
  tsne_data <- data$tsne_vis[[sensor]][[channel]]
  
  # Extract ground truth labels and predicted clusters
  ground_truth <- unlist(tsne_data$ground_truth)
  predictions <- unlist(tsne_data$kmeans_clusters)
  
  # Extract t-SNE projections
  projection_data <- tsne_data$projection
  x_coords <- sapply(projection_data, function(point) point[[1]])
  y_coords <- sapply(projection_data, function(point) point[[2]])
  
  tsne_projections <- data.frame(
    x = x_coords,
    y = y_coords,
    ground_truth = as.factor(ground_truth),
    predictions = as.factor(predictions)
  )
  
  # Get color palettes - ensure consistent colors between ground truth and predictions
  n_classes <- length(unique(c(ground_truth, predictions)))
  
  if(n_classes <= 8) {
    # For 8 or fewer classes, use the Dark2 palette
    colors <- brewer.pal(max(3, n_classes), "Dark2")
  } else {
    # For more classes, interpolate additional colors
    colors <- colorRampPalette(brewer.pal(8, "Dark2"))(n_classes)
  }
  
  # Create ground truth t-SNE plot with better styling
  gt_plot <- ggplot(tsne_projections, aes(x = x, y = y, color = ground_truth)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = colors) +
    labs(
      title = "Ground Truth",
      x = "t-SNE 1",
      y = "t-SNE 2",
      color = "Class"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # Create predictions t-SNE plot with better styling
  pred_plot <- ggplot(tsne_projections, aes(x = x, y = y, color = predictions)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = colors) +
    labs(
      title = "K-means Clusters",
      x = "t-SNE 1",
      y = "t-SNE 2",
      color = "Cluster"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # Extract metrics
  accuracy <- tsne_data$acc
  macro_f1 <- tsne_data$macro_f1
  micro_f1 <- tsne_data$micro_f1
  
  # Create text annotation with metrics
  metrics_text <- paste0(
    "Accuracy: ", round(accuracy, 3), "\n",
    "Macro F1: ", round(macro_f1, 3), "\n",
    "Micro F1: ", round(micro_f1, 3)
  )
  
  # Combine plots with shared title and metrics
  combined_plot <- gt_plot + pred_plot +
    plot_layout(guides = "collect") +  # Collect legends
    plot_annotation(
      title = paste0(ifelse(sensor == "HG", "PPHG", sensor), " ", gsub("ch", "Lead ", channel)),
      subtitle = metrics_text,
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )
  
  # Save individual plot if requested
  if(save_individual) {
    plot_file <- paste0(output_dir, "tsne_", sensor, "_", channel, ".png")
    ggsave(plot_file, combined_plot, width = 12, height = 8, dpi = 300, bg = "white")
    cat("Saved individual t-SNE plot for", sensor, channel, "to:", plot_file, "\n")
  }
  
  return(combined_plot)
}


create_profile_cluster_visualization <- function(data, sensor, channel, participants) {
  # Extract t-SNE data for the selected sensor-channel
  tsne_data <- data$tsne_vis[[sensor]][[channel]]
  
  # Extract ground truth labels and predicted clusters
  predictions <- unlist(tsne_data$kmeans_clusters)
  
  # Extract t-SNE projections
  projection_data <- tsne_data$projection
  x_coords <- sapply(projection_data, function(point) point[[1]])
  y_coords <- sapply(projection_data, function(point) point[[2]])
  
  # Create t-SNE projection dataframe
  tsne_df <- data.frame(
    x = x_coords,
    y = y_coords,
    cluster = as.factor(predictions)
  )
  
  # Create a list to store participant position information and waveform bubbles
  participant_info <- list()
  waveform_bubbles <- list()
  
  # For each participant, extract their position in t-SNE space and create waveform bubbles
  for (i in seq_along(participants)) {
    participant <- participants[i]
    
    # Check if participant exists in the dataset
    if (!is.null(data$centroids[[sensor]][[channel]][[participant]])) {
      # Get the centroid waveform
      centroid_waveform <- data$centroids[[sensor]][[channel]][[participant]]$kmeans_cluster_mapped_waveform
      
      if (!is.null(centroid_waveform)) {
        # Create data frame for plotting the waveform
        waveform_data <- data.frame(
          Time = (1:length(centroid_waveform) - 1) / fs,
          Amplitude = unlist(centroid_waveform)
        )
        
        # Normalize amplitude for better visualization in bubble
        waveform_data$Amplitude <- scale(waveform_data$Amplitude)
        
        # Create a MUCH larger waveform plot with better styling for bubble
        waveform_plot <- ggplot(waveform_data, aes(x = Time, y = Amplitude)) +
          geom_line(color = ifelse(sensor == "AgCl", "#d9020d", "#0066CC"), size = 3.0) +  # Extra thick line
          labs(title = paste0(participant)) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 20, hjust = 0.5),  # Larger title
            panel.background = element_rect(fill = "transparent", color = NA),
            plot.background = element_rect(fill = "transparent", color = NA),
            axis.title = element_blank(),
            axis.text = element_blank(),
            panel.grid = element_blank(),
            plot.margin = margin(10, 10, 10, 10)  # More internal margin
          )
        
        # Create bubble plot
        bubble_plot <- draw_bubble_panel(waveform_plot)
        waveform_bubbles[[participant]] <- bubble_plot
        
        # Store information about this participant
        participant_info[[participant]] <- list(
          participant = participant,
          label = paste0("P", i),
          bubble = bubble_plot
        )
      }
    }
  }
  
  n_clusters <- length(unique(tsne_df$cluster))

  # Use a reliable color palette approach that works for many clusters
  if(n_clusters <= 12) {
    # For 12 or fewer clusters, use the Paired palette directly
    colors <- brewer.pal(max(3, n_clusters), "Paired")
  } else {
    # For more than 12 clusters, interpolate additional colors
    colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_clusters)
  }

  # Create the base t-SNE scatter plot - ensure all clusters are visible with a good color palette
  tsne_base_plot <- ggplot(tsne_df, aes(x = x, y = y, color = cluster)) +
    geom_point(alpha = 0.7, size = 2.0) +  # Larger, more opaque points
    scale_color_manual(values = colors) +  # Force all clusters to have colors
    labs(
      title = paste0(ifelse(sensor == "HG", "PPHG", sensor), " ", gsub("ch", "Lead ", channel), " - Clusters"),
      x = "t-SNE 1",
      y = "t-SNE 2",
      color = "Cluster"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 40, face = "bold"),
      axis.title = element_text(size = 40),
      axis.text = element_text(size = 40),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.8, "lines"),  # Smaller legend keys
      legend.spacing.y = unit(0.5, "lines")  # Less space between legend items
    ) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))  # Force legend to show all clusters, use 2 columns
  
  # Get the bounds of the t-SNE plot
  x_range <- range(tsne_df$x, na.rm = TRUE)
  y_range <- range(tsne_df$y, na.rm = TRUE)
  
  # Ensure ranges are finite
  if(!all(is.finite(c(x_range, y_range)))) {
    warning("Non-finite values found in t-SNE coordinates. Using default ranges.")
    x_range <- c(0, 100)
    y_range <- c(0, 100)
  }
  
  # Ensure range differences are not zero
  if(diff(x_range) <= 0) x_range <- c(x_range[1], x_range[1] + 100)
  if(diff(y_range) <= 0) y_range <- c(y_range[1], y_range[1] + 100)
  
  # Calculate padding based on range
  x_padding <- diff(x_range) * 0.8
  
  # Position bubbles with fixed spacing to avoid overlap
  x_bubble_pos <- rep(x_range[2] + (x_padding * 0.4), length(participants))
  
  # Use fixed spacing between bubbles based on bubble size
  bubble_size <- 15
  bubble_spacing <- bubble_size * 2.2  # Leave 20% extra space between bubbles
  
  # Calculate total height needed for all bubbles
  total_height_needed <- bubble_spacing * length(participants)
  
  # Calculate vertical positions with fixed spacing
  if(length(participants) == 1) {
    # For a single participant, center the bubble
    y_bubble_pos <- mean(y_range)
    } else {
    # For multiple participants, we need to ensure bubbles stay inside the plot

    # Calculate the necessary margins at top and bottom (as a ratio of bubble size)
    margin_size <- bubble_size * 1.1  # Slightly larger than bubble size

    # Calculate the available plot height for placing bubbles
    available_height <- diff(y_range) - 2 * margin_size

    # Calculate spacing between bubbles
    if(length(participants) > 1) {
        spacing <- available_height / (length(participants) - 1)
    } else {
        spacing <- available_height
    }

    # If spacing is too small, we need to extend the plot height
    min_spacing <- bubble_size * 2.2  # Minimum acceptable spacing

    if(spacing < min_spacing && length(participants) > 1) {
        # Extend the y-range to accommodate all bubbles with proper spacing
        total_needed_height <- margin_size * 2 + (length(participants) - 1) * min_spacing
        
        # Calculate how much we need to extend the range
        extension_needed <- total_needed_height - diff(y_range)
        if(extension_needed > 0) {
        # Extend the range symmetrically
        y_range[1] <- y_range[1] - extension_needed/2
        y_range[2] <- y_range[2] + extension_needed/2
        }
        
        # Recalculate available height and spacing
        available_height <- diff(y_range) - 2 * margin_size
        spacing <- available_height / (length(participants) - 1)
    }

    # Calculate bubble positions with proper margins
    y_start <- y_range[1] + margin_size
    y_bubble_pos <- y_start + ((0:(length(participants)-1)) * spacing)
  }

  # Ensure all positions are finite
  if(!all(is.finite(y_bubble_pos))) {
    warning("Non-finite bubble positions calculated. Using default positions.")
    y_bubble_pos <- seq(y_range[1] + diff(y_range)*0.1, 
                        y_range[2] - diff(y_range)*0.1, 
                        length.out = length(participants))
  }

  # Later in the function, update the coord_cartesian section:
  coord_cartesian(
    xlim = c(min(x_range[1], min(tsne_df$x, na.rm = TRUE)) - diff(x_range)*0.05, 
            max(x_range[2], max(tsne_df$x, na.rm = TRUE)) + x_padding),
    ylim = y_range,  # Use the potentially extended y-range
    expand = FALSE
  )
  
  # Sample specific points from the t-SNE plot to represent participants - try to pick points from different clusters
  set.seed(42) # For reproducibility
  
  # Get unique clusters
  unique_clusters <- unique(as.numeric(as.character(tsne_df$cluster)))
  
  # Try to select points from different clusters if possible
  sample_indices <- numeric(length(participants))
  
  if(length(unique_clusters) >= length(participants)) {
    # If we have enough clusters, pick one point from each cluster
    selected_clusters <- sample(unique_clusters, length(participants))
    for(i in 1:length(participants)) {
      cluster_points <- which(as.numeric(as.character(tsne_df$cluster)) == selected_clusters[i])
      if(length(cluster_points) > 0) {
        sample_indices[i] <- sample(cluster_points, 1)
      } else {
        sample_indices[i] <- sample(1:nrow(tsne_df), 1)
      }
    }
  } else {
    # Otherwise just sample randomly
    sample_indices <- sample(1:nrow(tsne_df), length(participants))
  }
  
  # Create the final plot with bubbles
  final_plot <- tsne_base_plot
  
  # Add bubbles and connecting arrows
  for (i in seq_along(participants)) {
    participant <- participants[i]
    if (!is.null(participant_info[[participant]])) {
      # Get a point in t-SNE space for this participant
      point_index <- sample_indices[i]
      x_point <- tsne_df$x[point_index]
      y_point <- tsne_df$y[point_index]
      
      # Skip if points are not finite
      if(!all(is.finite(c(x_point, y_point, x_bubble_pos[i], y_bubble_pos[i])))) {
        warning(paste("Skipping non-finite coordinates for participant", participant))
        next
      }
      
      # Add bubble annotation - with explicit bounds check
      final_plot <- final_plot +
        annotation_custom(
          grob = ggplotGrob(waveform_bubbles[[participant]]),
          xmin = x_bubble_pos[i] - bubble_size,
          xmax = x_bubble_pos[i] + bubble_size,
          ymin = y_bubble_pos[i] - bubble_size,
          ymax = y_bubble_pos[i] + bubble_size
        )
      
      # Add connecting arrow with bounds check
      final_plot <- final_plot +
        geom_segment(
          data = data.frame(x = x_point, y = y_point, 
                           xend = x_bubble_pos[i] - bubble_size, 
                           yend = y_bubble_pos[i]),
          aes(x = x, y = y, xend = xend, yend = yend),
          inherit.aes = FALSE,
          arrow = arrow(type = "closed", length = unit(0.3, "cm")),
          color = "gray30",
          linetype = "dashed",
          size = 1.0
        )
      
      # [Rest of the point labeling code]
    }
  }
  
  # Update the plot area bounds with explicit limits
  final_plot <- final_plot +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    ) +
    # Use explicit limits to avoid non-finite values
    coord_cartesian(
      xlim = c(min(x_range[1], min(tsne_df$x, na.rm = TRUE)) - diff(x_range)*0.05, 
              max(x_range[2], max(tsne_df$x, na.rm = TRUE)) + x_padding),
      expand = FALSE
    )
  
  return(final_plot)
}

# Create clustering performance visualizations
clustering_plots <- list()

# Replace the clustering performance visualization section with individual plot saving
for (sensor in sensors) {
  for (channel in channels) {
    # Create and save individual plots
    create_tsne_plots(data, sensor, channel, save_individual = TRUE)
  }
}

cat("Individual t-SNE plots saved to:", output_dir, "\n")


profile_vis <- create_profile_cluster_visualization(data, profile_sensor, profile_channel, selected_participants_profiles)
profile_file <- paste0(output_dir, "participant_profiles_", ifelse(profile_sensor == "HG", "PPHG", profile_sensor), "_Lead_", substr(profile_channel, 3, 3), ".png")

# Set a reasonable height based on number of participants
plot_height <- 20 + (length(selected_participants_profiles) - 5) * 3  # Base height 16, +2 for each participant over 5
plot_height <- max(16, min(30, plot_height))  # Constrain between 16 and 30 inches

ggsave(profile_file, profile_vis, width = 30, height = plot_height, dpi = 300, bg = "white", limitsize = FALSE)
cat("Profile visualization saved to:", profile_file, "\n")