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

data_path <- "/home/giannis/Documents/ECG HG paper/results_data/channel_id_results.json"
data <- fromJSON(txt = data_path, simplifyVector = FALSE)
fs <- 200
sensors <- c("AgCl", "HG")

# Create output directory
output_dir <- paste0("R_figures/figure_2c/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

channels <- c("ch1", "ch2", "ch3")

# Print the structure to verify our understanding
cat("Data structure:\n")
cat("tsne_vis keys:", paste(names(data$tsne_vis), collapse=", "), "\n")
for (sensor in sensors) {
  if (!is.null(data$tsne_vis[[sensor]])) {
    cat("Keys in tsne_vis for", sensor, ":", paste(names(data$tsne_vis[[sensor]]), collapse=", "), "\n")
  }
}


draw_bubble_panel <- function(plot) {
  ggdraw() +
    draw_grob(
      circleGrob(
        x = 0.5, y = 0.5, r = 1,
        gp = gpar(fill = "white", col = "gray50", lwd = 2.5)  # Even thicker border
      )
    ) +
    draw_plot(plot, x = 0.1, y = 0.1, width = 0.8, height = 0.8)  # Larger content area
}

# Modified function to create t-SNE plots with correct class-cluster correspondence
create_tsne_plots <- function(data, sensor, save_individual = TRUE) {
  # Extract t-SNE data - directly from sensor level
  tsne_data <- data$tsne_vis[[sensor]]
  
  # Check if we have valid tsne_data
  if (is.null(tsne_data)) {
    warning(paste("No t-SNE data available for sensor:", sensor))
    return(NULL)
  }
  
  # Extract ground truth labels and predicted clusters
  ground_truth <- unlist(tsne_data$ground_truth)
  predictions <- unlist(tsne_data$kmeans_clusters)
  
  # Check if we have valid ground truth and predictions
  if (is.null(ground_truth) || length(ground_truth) == 0) {
    warning(paste("No ground truth labels for", sensor))
    return(NULL)
  }
  
  if (is.null(predictions) || length(predictions) == 0) {
    warning(paste("No predictions for", sensor))
    return(NULL)
  }
  
  # Extract t-SNE projections
  projection_data <- tsne_data$projection
  
  # Check if we have valid projection data
  if (is.null(projection_data) || length(projection_data) == 0) {
    warning(paste("No projection data for", sensor))
    return(NULL)
  }
  
  # Extract coordinates with robust error handling
  x_coords <- numeric(0)
  y_coords <- numeric(0)
  
  for (i in seq_along(projection_data)) {
    point <- projection_data[[i]]
    if (is.list(point) && length(point) >= 2) {
      x_coords <- c(x_coords, point[[1]])
      y_coords <- c(y_coords, point[[2]])
    } else if (is.numeric(point) && length(point) >= 2) {
      x_coords <- c(x_coords, point[1])
      y_coords <- c(y_coords, point[2])
    }
  }
  
  # Check if we have valid coordinates
  if (length(x_coords) == 0 || length(y_coords) == 0) {
    warning(paste("Failed to extract valid coordinates for", sensor))
    return(NULL)
  }
  
  # Make sure dimensions match
  min_length <- min(length(x_coords), length(y_coords), length(ground_truth), length(predictions))
  x_coords <- x_coords[1:min_length]
  y_coords <- y_coords[1:min_length]
  ground_truth <- ground_truth[1:min_length]
  predictions <- predictions[1:min_length]
  
  tsne_projections <- data.frame(
    x = x_coords,
    y = y_coords,
    ground_truth = as.factor(ground_truth),
    predictions = as.factor(predictions)
  )
  
  # Define the color mapping based on the correspondence provided
  # We'll use the same colors for matching classes/clusters
  if (sensor == "AgCl") {
    # AgCl correspondence: 
    # class 1 → cluster 2, class 2 → cluster 1, class 3 → cluster 3
    class_colors <- c("1" = "#377EB8", "2" = "#4DAF4A", "3" = "#E41A1C")
    cluster_colors <- c("1" = "#4DAF4A", "2" = "#377EB8", "3" = "#E41A1C")
  } else if (sensor == "HG") {
    # PPHG correspondence:
    # class 1 → cluster 2, class 2 → cluster 3, class 3 → cluster 1
    class_colors <- c("1" = "#377EB8", "2" = "#4DAF4A", "3" = "#E41A1C")
    cluster_colors <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A")
  } else {
    # Default colors for any other sensors
    class_colors <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A")
    cluster_colors <- class_colors
  }
  
  # Create ground truth t-SNE plot with correct colors
  gt_plot <- ggplot(tsne_projections, aes(x = x, y = y, color = ground_truth)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = class_colors) +
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
  
  # Create predictions t-SNE plot with correct correspondence colors
  pred_plot <- ggplot(tsne_projections, aes(x = x, y = y, color = predictions)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = cluster_colors) +
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
  
  # Add correspondence information
  if (sensor == "AgCl") {
    correspondence_text <- "Class-cluster correspondence: Class 1 → Cluster 2, Class 2 → Cluster 1, Class 3 → Cluster 3"
  } else if (sensor == "HG") {
    correspondence_text <- "Class-cluster correspondence: Class 1 → Cluster 2, Class 2 → Cluster 3, Class 3 → Cluster 1"
  } else {
    correspondence_text <- ""
  }
  
  # Combine plots with shared title, metrics and correspondence information
  combined_plot <- gt_plot + pred_plot +
    plot_layout(guides = "collect") +  # Collect legends
    plot_annotation(
      title = paste0(ifelse(sensor == "HG", "PPHG", sensor), " Channel Clustering"),
      subtitle = metrics_text,
      caption = correspondence_text,
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(size = 10, hjust = 0.5)
      )
    )
  
  # Save individual plot if requested
  if(save_individual) {
    plot_file <- paste0(output_dir, "tsne_", sensor, ".png")
    ggsave(plot_file, combined_plot, width = 12, height = 8, dpi = 300, bg = "white")
    cat("Saved t-SNE plot for", sensor, "to:", plot_file, "\n")
  }
  
  return(combined_plot)
}

create_sensor_cluster_visualization <- function(data, sensor) {
  # Extract t-SNE data - directly from sensor level
  tsne_data <- data$tsne_vis[[sensor]]
  
  # Check if we have valid tsne_data
  if (is.null(tsne_data)) {
    stop(paste("No t-SNE data available for sensor:", sensor))
  }
  
  # Extract predicted clusters
  predictions <- unlist(tsne_data$kmeans_clusters)
  
  # Check if we have valid predictions
  if (is.null(predictions) || length(predictions) == 0) {
    stop(paste("No predictions for", sensor))
  }
  
  # Get unique clusters
  unique_clusters <- sort(unique(predictions))
  
  # Create a list to store cluster waveform bubbles
  cluster_waveform_bubbles <- list()
  
  # Extract t-SNE projections
  projection_data <- tsne_data$projection
  
  # Check if we have valid projection data
  if (is.null(projection_data) || length(projection_data) == 0) {
    stop(paste("No projection data for", sensor))
  }
  
  # Extract coordinates with robust error handling
  x_coords <- numeric(0)
  y_coords <- numeric(0)
  
  for (i in seq_along(projection_data)) {
    point <- projection_data[[i]]
    if (is.list(point) && length(point) >= 2) {
      x_coords <- c(x_coords, point[[1]])
      y_coords <- c(y_coords, point[[2]])
    } else if (is.numeric(point) && length(point) >= 2) {
      x_coords <- c(x_coords, point[1])
      y_coords <- c(y_coords, point[2])
    }
  }
  
  # Check if we have valid coordinates
  if (length(x_coords) == 0 || length(y_coords) == 0) {
    stop(paste("Failed to extract valid coordinates for", sensor))
  }
  
  # Make sure dimensions match
  min_length <- min(length(x_coords), length(y_coords), length(predictions))
  x_coords <- x_coords[1:min_length]
  y_coords <- y_coords[1:min_length]
  predictions <- predictions[1:min_length]
  
  # Create t-SNE projection dataframe
  tsne_df <- data.frame(
    x = x_coords,
    y = y_coords,
    cluster = as.factor(predictions)
  )
  
  # Print data structure to help debug
  cat("Data structure for", sensor, ":\n")
  if (!is.null(data$centroids) && !is.null(data$centroids[[sensor]])) {
    cat("  Found centroids for", sensor, "with channels:", 
        paste(names(data$centroids[[sensor]]), collapse=", "), "\n")
    
    for (channel in channels) {
      if (!is.null(data$centroids[[sensor]][[channel]])) {
        cat("    Channel", channel, "has centroids data available\n")
        if (!is.null(data$centroids[[sensor]][[channel]]$kmeans_cluster_mapped_waveform)) {
          cat("      With kmeans_cluster_mapped_waveform available\n")
        }
      }
    }
  } else {
    cat("  No centroids found for", sensor, "\n")
  }
  
  # Assign specific channel to each cluster
  # We'll distribute the channels evenly among the clusters
  cluster_channel_mapping <- list()
  for (i in seq_along(unique_clusters)) {
    cluster_id <- unique_clusters[i]
    # Assign channel based on cluster index modulo number of channels
    channel_idx <- ((i - 1) %% length(channels)) + 1
    cluster_channel_mapping[[as.character(cluster_id)]] <- channels[channel_idx]
  }
  
  # For each cluster, get the assigned channel's waveform
  waveforms_found <- FALSE
  
  for (cluster_id in unique_clusters) {
    cluster_str <- as.character(cluster_id)
    assigned_channel <- cluster_channel_mapping[[cluster_str]]
    
    # Check if centroids are available for this channel
    if (!is.null(data$centroids) && 
        !is.null(data$centroids[[sensor]]) && 
        !is.null(data$centroids[[sensor]][[assigned_channel]]) &&
        !is.null(data$centroids[[sensor]][[assigned_channel]]$kmeans_cluster_mapped_waveform)) {
      
      # Get centroid data for this channel
      centroid_waveform <- data$centroids[[sensor]][[assigned_channel]]$kmeans_cluster_mapped_waveform
      
      # Check if we have a valid waveform
      if (length(centroid_waveform) > 0) {
        # Create data frame for plotting the waveform
        waveform_data <- data.frame(
          Time = (1:length(centroid_waveform) - 1) / fs,
          Amplitude = unlist(centroid_waveform)
        )
        
        # Normalize amplitude for better visualization
        waveform_data$Amplitude <- scale(waveform_data$Amplitude)
        
        # Create a waveform plot with styling for bubble
        waveform_plot <- ggplot(waveform_data, aes(x = Time, y = Amplitude)) +
          geom_line(color = ifelse(sensor == "AgCl", "#d9020d", "#0066CC"), size = 3.0) +
          labs(title = paste0(gsub("ch", "Lead ", assigned_channel))) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 25, hjust = 0.5),
            panel.background = element_rect(fill = "transparent", color = NA),
            plot.background = element_rect(fill = "transparent", color = NA),
            axis.title = element_blank(),
            axis.text = element_blank(),
            panel.grid = element_blank(),
            plot.margin = margin(10, 10, 10, 10)
          )
        
        # Create bubble plot
        bubble_plot <- draw_bubble_panel(waveform_plot)
        
        # Store in our list - only one entry per cluster now
        cluster_waveform_bubbles[[cluster_str]] <- bubble_plot
        waveforms_found <- TRUE
        
        cat("Using waveform for", sensor, "cluster", cluster_id, "from channel", assigned_channel, "\n")
      }
    }
  }
  
  # Check if we found any waveforms - if not, stop with an error
  if (!waveforms_found) {
    stop(paste("No valid waveform data found for any cluster in", sensor))
  }
  
  # Use a good color palette
  if(length(unique_clusters) <= 12) {
    colors <- brewer.pal(max(3, length(unique_clusters)), "Paired")
  } else {
    colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique_clusters))
  }
  
  # Create the base t-SNE scatter plot
  tsne_plot <- ggplot(tsne_df, aes(x = x, y = y, color = cluster)) +
    geom_point(alpha = 0.7, size = 2.0) +
    scale_color_manual(values = colors) +
    labs(
      title = "",
      x = "t-SNE 1",
      y = "t-SNE 2",
      color = "Cluster"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 30, face = "bold"),
      axis.title = element_text(size = 30),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  # Get the bounds of the t-SNE plot
  x_range <- range(tsne_df$x, na.rm = TRUE)
  y_range <- range(tsne_df$y, na.rm = TRUE)
  
  # Ensure ranges are finite
  if(!all(is.finite(c(x_range, y_range)))) {
    stop("Non-finite values found in t-SNE coordinates.")
  }
  
  # Ensure range differences are not zero
  if(diff(x_range) <= 0) x_range <- c(x_range[1], x_range[1] + 100)
  if(diff(y_range) <= 0) y_range <- c(y_range[1], y_range[1] + 100)
  
  # Calculate padding and positioning for bubbles
  bubble_size <- 15
  x_padding <- diff(x_range) * 0.8
  
  # Calculate cluster centroids for connecting arrows
  cluster_centroids <- data.frame(
    cluster = numeric(0),
    x = numeric(0),
    y = numeric(0)
  )
  
  for (cluster_id in unique_clusters) {
    cluster_points <- which(as.numeric(as.character(tsne_df$cluster)) == cluster_id)
    if (length(cluster_points) > 0) {
      centroid_x <- mean(tsne_df$x[cluster_points], na.rm = TRUE)
      centroid_y <- mean(tsne_df$y[cluster_points], na.rm = TRUE)
      cluster_centroids <- rbind(cluster_centroids, 
                                data.frame(cluster = cluster_id, 
                                          x = centroid_x, 
                                          y = centroid_y))
    }
  }
  
  # Start with base t-SNE plot
  final_plot <- tsne_plot
  
  # Calculate the arrangement of bubbles
  # We'll position bubbles around the scatter plot based on number of clusters
  clusters_with_waveforms <- names(cluster_waveform_bubbles)
  
  if (length(clusters_with_waveforms) == 0) {
    stop("No clusters have valid waveforms")
  }
  
  # For small number of bubbles, place them on the right side
  n_bubbles <- length(clusters_with_waveforms)
  
  if (n_bubbles <= 4) {
    x_bubble_pos <- rep(x_range[2] + (x_padding * 0.4), n_bubbles)
    y_bubble_pos <- seq(y_range[1], y_range[2], length.out = n_bubbles)
  } 
  # For more bubbles, distribute them on right, top, and left sides
  else {
    # Right side
    right_count <- min(4, n_bubbles)
    x_right <- rep(x_range[2] + (x_padding * 0.4), right_count)
    y_right <- seq(y_range[1], y_range[2], length.out = right_count)
    
    # If we have more, add to the top
    remaining <- n_bubbles - right_count
    if (remaining > 0) {
      top_count <- min(3, remaining)
      x_top <- seq(x_range[1] + diff(x_range)*0.25, 
                  x_range[1] + diff(x_range)*0.75, 
                  length.out = top_count)
      y_top <- rep(y_range[2] + (diff(y_range) * 0.4), top_count)
      
      # If still more, add to the left
      remaining <- remaining - top_count
      if (remaining > 0) {
        left_count <- min(4, remaining)
        x_left <- rep(x_range[1] - (x_padding * 0.4), left_count)
        y_left <- seq(y_range[2], y_range[1], length.out = left_count)
        
        x_bubble_pos <- c(x_right, x_top, x_left)
        y_bubble_pos <- c(y_right, y_top, y_left)
      } else {
        x_bubble_pos <- c(x_right, x_top)
        y_bubble_pos <- c(y_right, y_top)
      }
    } else {
      x_bubble_pos <- x_right
      y_bubble_pos <- y_right
    }
  }
  
  # For each cluster with waveforms, add bubble and connecting arrow
  for (i in seq_along(clusters_with_waveforms)) {
    cluster_str <- clusters_with_waveforms[i]
    cluster_id <- as.numeric(cluster_str)
    bubble_plot <- cluster_waveform_bubbles[[cluster_str]]
    
    # Find this cluster's centroid
    centroid_idx <- which(cluster_centroids$cluster == cluster_id)
    if (length(centroid_idx) == 0) next
    
    # Get centroid coordinates
    x_centroid <- cluster_centroids$x[centroid_idx]
    y_centroid <- cluster_centroids$y[centroid_idx]
    
    # Add cluster label near the centroid
    final_plot <- final_plot +
      geom_label(
        data = data.frame(x = x_centroid, y = y_centroid, label = paste0("C", cluster_id)),
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        fill = "white",
        color = "black",
        label.size = 0.5,
        size = 4,
        fontface = "bold"
      )
    
    # Get bubble position
    x_bubble <- x_bubble_pos[i]
    y_bubble <- y_bubble_pos[i]
    
    # Draw arrow from centroid to bubble
    arrow_end_x <- x_bubble
    arrow_end_y <- y_bubble
    
    # Adjust arrow end to meet bubble edge
    if (x_bubble > x_centroid) {
      arrow_end_x <- x_bubble - bubble_size
    } else if (x_bubble < x_centroid) {
      arrow_end_x <- x_bubble + bubble_size
    }
    
    if (y_bubble > y_centroid) {
      arrow_end_y <- y_bubble - bubble_size
    } else if (y_bubble < y_centroid) {
      arrow_end_y <- y_bubble + bubble_size
    }
    
    final_plot <- final_plot +
      geom_segment(
        data = data.frame(x = x_centroid, y = y_centroid, 
                          xend = arrow_end_x, yend = arrow_end_y),
        aes(x = x, y = y, xend = xend, yend = yend),
        inherit.aes = FALSE,
        arrow = arrow(type = "closed", length = unit(0.3, "cm")),
        color = "gray30",
        linetype = "dashed",
        size = 1.5
      )
    
    # Add bubble
    final_plot <- final_plot +
      annotation_custom(
        grob = ggplotGrob(bubble_plot),
        xmin = x_bubble - bubble_size,
        xmax = x_bubble + bubble_size,
        ymin = y_bubble - bubble_size,
        ymax = y_bubble + bubble_size
      )
  }
  

  
  # Calculate expanded plot limits to include bubbles
  x_min <- min(c(x_bubble_pos) - bubble_size) - diff(x_range)*0.1
  x_max <- max(c(x_bubble_pos) + bubble_size) + diff(x_range)*0.1
  y_min <- min(c(y_bubble_pos) - bubble_size) - diff(y_range)*0.1
  y_max <- max(c(y_bubble_pos) + bubble_size) + diff(y_range)*0.1
  
  # Ensure limits include the scatter plot area
  x_min <- min(x_min, x_range[1] - diff(x_range)*0.1)
  x_max <- max(x_max, x_range[2] + diff(x_range)*0.1)
  y_min <- min(y_min, y_range[1] - diff(y_range)*0.1)
  y_max <- max(y_max, y_range[2] + diff(y_range)*0.1)
  
  # Set plot limits to include bubbles
  final_plot <- final_plot +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    coord_cartesian(
      xlim = c(x_min, x_max),
      ylim = c(y_min, y_max),
      expand = FALSE
    )
  
  return(final_plot)
}

# Process each sensor to create t-SNE plots and cluster visualizations
sensor_plots <- list()
tsne_plots <- list()

for (sensor in sensors) {
  cat("\nProcessing sensor:", sensor, "\n")
  
  # Create the t-SNE plots
  cat("Creating t-SNE plots...\n")
  tsne_plot <- tryCatch({
    create_tsne_plots(data, sensor)
  }, error = function(e) {
    cat("Error creating t-SNE plot for", sensor, ":", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(tsne_plot)) {
    tsne_plots[[sensor]] <- tsne_plot
  }
  
  # Create the cluster visualization - with proper error handling but no synthetic data
  cat("Creating cluster visualization...\n")
  cluster_plot <- tryCatch({
    create_sensor_cluster_visualization(data, sensor)
  }, error = function(e) {
    cat("Error creating cluster visualization for", sensor, ":", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(cluster_plot)) {
    sensor_plots[[sensor]] <- cluster_plot
    
    # Save the plot
    plot_file <- paste0(output_dir, "channel_clusters_", sensor, ".png")
    ggsave(plot_file, cluster_plot, width = 24, height = 16, dpi = 300, bg = "white", limitsize = FALSE)
    cat("Saved cluster visualization for", sensor, "to:", plot_file, "\n")
  }
}

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


# Improved function to create a complete visualization for each sensor with side bubbles
create_improved_visualization <- function(data, sensor) {
  # Extract t-SNE data
  tsne_data <- data$tsne_vis[[sensor]]
  
  # Get predictions and ground truth
  predictions <- unlist(tsne_data$kmeans_clusters)
  ground_truth <- unlist(tsne_data$ground_truth)
  
  # Get projections
  projection_data <- tsne_data$projection
  
  # Extract coordinates
  x_coords <- numeric(0)
  y_coords <- numeric(0)
  
  for (i in seq_along(projection_data)) {
    point <- projection_data[[i]]
    if (is.list(point) && length(point) >= 2) {
      x_coords <- c(x_coords, point[[1]])
      y_coords <- c(y_coords, point[[2]])
    } else if (is.numeric(point) && length(point) >= 2) {
      x_coords <- c(x_coords, point[1])
      y_coords <- c(y_coords, point[2])
    }
  }
  
  # Make sure dimensions match
  min_length <- min(length(x_coords), length(y_coords), length(predictions), length(ground_truth))
  x_coords <- x_coords[1:min_length]
  y_coords <- y_coords[1:min_length]
  predictions <- predictions[1:min_length]
  ground_truth <- ground_truth[1:min_length]
  
  # Create dataframe
  tsne_df <- data.frame(
    x = x_coords,
    y = y_coords,
    cluster = as.factor(predictions),
    ground_truth = as.factor(ground_truth)
  )
  
  # Define colors based on correct correspondence
  if (sensor == "AgCl") {
    # AgCl correspondence: class 1 → cluster 2, class 2 → cluster 1, class 3 → cluster 3
    class_colors <- c("1" = "#377EB8", "2" = "#4DAF4A", "3" = "#E41A1C")
    cluster_colors <- c("1" = "#4DAF4A", "2" = "#377EB8", "3" = "#E41A1C")
  } else if (sensor == "HG") {
    # PPHG correspondence: class 1 → cluster 2, class 2 → cluster 3, class 3 → cluster 1
    class_colors <- c("1" = "#377EB8", "2" = "#4DAF4A", "3" = "#E41A1C")
    cluster_colors <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A")
  } else {
    # Default colors
    class_colors <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A")
    cluster_colors <- class_colors
  }
  
  # Create ground truth plot
  gt_plot <- ggplot(tsne_df, aes(x = x, y = y, color = ground_truth)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = class_colors) +
    labs(
      title = "Ground Truth",
      x = "t-SNE 1",
      y = "t-SNE 2",
      color = "Class"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 30, face = "bold"),
      axis.title = element_text(size = 30),
      axis.text = element_blank(),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # Create predictions plot
  pred_plot <- ggplot(tsne_df, aes(x = x, y = y, color = cluster)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = cluster_colors) +
    labs(
      title = "K-means Clusters",
      x = "t-SNE 1",
      y = "t-SNE 2",
      color = "Cluster"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 30, face = "bold"),
      axis.title = element_text(size = 30),
      axis.text = element_blank(),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # Map each class to its corresponding lead (class 1 -> lead 1, etc.)
  unique_classes <- sort(unique(ground_truth))
  class_lead_mapping <- list()
  
  for (i in seq_along(unique_classes)) {
    class_id <- unique_classes[i]
    class_lead_mapping[[as.character(class_id)]] <- channels[i]
  }
  
  # Add cluster labels to the prediction plot
  # Calculate cluster centroids
  cluster_centroids <- data.frame(
    cluster = numeric(0),
    x = numeric(0),
    y = numeric(0)
  )
  
  for (cluster_id in sort(unique(predictions))) {
    cluster_points <- which(as.numeric(as.character(tsne_df$cluster)) == cluster_id)
    if (length(cluster_points) > 0) {
      centroid_x <- mean(tsne_df$x[cluster_points], na.rm = TRUE)
      centroid_y <- mean(tsne_df$y[cluster_points], na.rm = TRUE)
      cluster_centroids <- rbind(cluster_centroids, 
                                data.frame(cluster = cluster_id, 
                                          x = centroid_x, 
                                          y = centroid_y))
    }
  }
  
  # Add cluster labels to the prediction plot
  for (i in 1:nrow(cluster_centroids)) {
    pred_plot <- pred_plot +
      geom_label(
        data = data.frame(
          x = cluster_centroids$x[i], 
          y = cluster_centroids$y[i], 
          label = paste0("C", cluster_centroids$cluster[i])
        ),
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        fill = "white",
        color = "black",
        label.size = 0.5,
        size = 3.5,
        fontface = "bold"
      )
  }
  
  # Create waveform bubbles
  class_waveform_bubbles <- list()
  
  for (class_id in unique_classes) {
    class_str <- as.character(class_id)
    assigned_channel <- class_lead_mapping[[class_str]]
    
    # Check if centroids are available
    if (!is.null(data$centroids) && 
        !is.null(data$centroids[[sensor]]) && 
        !is.null(data$centroids[[sensor]][[assigned_channel]]) &&
        !is.null(data$centroids[[sensor]][[assigned_channel]]$kmeans_cluster_mapped_waveform)) {
      
      # Get waveform data
      centroid_waveform <- data$centroids[[sensor]][[assigned_channel]]$kmeans_cluster_mapped_waveform
      
      if (length(centroid_waveform) > 0) {
        waveform_data <- data.frame(
          Time = (1:length(centroid_waveform) - 1) / fs,
          Amplitude = scale(unlist(centroid_waveform))
        )
        
        # Create waveform plot with class color
        waveform_plot <- ggplot(waveform_data, aes(x = Time, y = Amplitude)) +
          geom_line(color = class_colors[class_str], size = 2.5) +
          labs(title = paste0(gsub("ch", "Lead ", assigned_channel))) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 18, hjust = 0.5),
            panel.background = element_rect(fill = "transparent", color = NA),
            plot.background = element_rect(fill = "transparent", color = NA),
            axis.title = element_blank(),
            axis.text = element_blank(),
            panel.grid = element_blank(),
            plot.margin = margin(5, 5, 5, 5)
          )
        
        # Create bubble
        bubble_plot <- draw_bubble_panel(waveform_plot)
        
        # Store bubble
        class_waveform_bubbles[[class_str]] <- list(
          plot = bubble_plot,
          channel = assigned_channel,
          class_id = class_id
        )
      }
    }
  }

  
  # Create the main panel with the two t-SNE plots stacked vertically
  main_panel <- plot_grid(
    gt_plot + theme(aspect.ratio = 1.4), 
    pred_plot + theme(aspect.ratio = 1.4), 
    ncol = 2,
    align = "h"
  )
  
  # Create small individual bubble panels for each class
  bubble_panels <- list()
  for (class_str in names(class_waveform_bubbles)) {
    bubble_info <- class_waveform_bubbles[[class_str]]
    
    # Create a labeled bubble panel
    bubble_panel <- ggdraw() + 
      draw_grob(ggplotGrob(bubble_info$plot)) 
    
    bubble_panels[[class_str]] <- bubble_panel
  }
  
  # Arrange the bubbles in a vertical grid
  if (length(bubble_panels) > 0) {
    bubbles_grid <- plot_grid(
      plotlist = bubble_panels,
      ncol = 1,
      align = "v"
    )
    
    # Combine the main panel and bubbles grid
    final_panel <- plot_grid(
      main_panel, bubbles_grid,
      ncol = 2,
      rel_widths = c(2, 0.8)  # Main panel gets 2/3 of width, bubbles get 1/3
    )
  } else {
    final_panel <- main_panel
  }
  
  # Add title, metrics and correspondence as top and bottom annotations
  title <- ggdraw() 
  
  
  # Put everything together
  complete_plot <- plot_grid(
    final_panel,
    ncol = 1,
    rel_heights = c(0.99)
  )
  
  return(complete_plot)
}

# Create individual plots for each sensor
cat("\nGenerating individual improved visualizations for each sensor...\n")

for (sensor in sensors) {
  cat("Creating improved visualization for", sensor, "...\n")
  
  improved_plot <- tryCatch({
    create_improved_visualization(data, sensor)
  }, error = function(e) {
    cat("Error creating improved visualization for", sensor, ":", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(improved_plot)) {
    # Save the plot
    plot_file <- paste0(output_dir, "improved_visualization_", sensor, ".png")
    # Use a wider aspect ratio to ensure bubbles fit properly
    ggsave(plot_file, improved_plot, width = 14, height = 10, dpi = 300, bg = "white")
    cat("Saved improved visualization for", sensor, "to:", plot_file, "\n")
  }
}

cat("All visualizations complete.\n")

