# Twinfoot pre-processing code for pressure data
# uses functions from generic_pressure_functions.R

# =============================================================================

# libraries required
library(RNiftyReg)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(ggpubr)


# =============================================================================

#' Register and transform pressure image
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param source Array. A 3D array covering each timepoint of the measurement. z
#'   dimension represents time. This will be registered to the template
#' @param template Array. A 3D array covering each timepoint of the measurement.
#'   z dimension represents time. The source data will be registered to
#' @param match_type String. "max' or "min" pressure image
#' @return Array. A 3D array covering each timepoint of the measurement
#'   (normalized to 101 points)

register_emed <- function(source, template, match_type = "max") {
  # template file footprint
  temp_fprint <- footprint(template, value = match_type)
  temp_fprint[temp_fprint > 0.1] <- 1
  
  # source file footprint
  source_fprint <- footprint(source, value = match_type)
  source_fprint[source_fprint > 0.1] <- 1
  
  # generate transformation matrix
  reg <- niftyreg.linear(source_fprint, temp_fprint, scope = "affine")
  nifty_transform <- forward(reg)
  
  # normalize to 101 points for analysis
  source_101 <- pressure_interp(source, 101)
  template_101 <- pressure_interp(template, 101)
  
  # transform entire trial
  source_reg <- template_101
  for (i in 1:101) {
    y <- RNiftyReg::applyTransform(nifty_transform, source_101[, , i], 
                                   interpolation = 1)
    source_reg[, , i] <- matrix(as.numeric(y), nrow = dim(source_reg)[1], 
                                ncol = dim(source_reg)[2])
  }
  
  # return registered footprint
  return(source_reg)
}


# =============================================================================

#' Find average pressure footprint within subject
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param template_path String. Path to file to be used as template
#' @param list_of_source_paths Vector of strings. Paths to all files that are to
#'   be aligned to templated
#' @param flip Logical. If file is to be mirrored around its vertical axis

mean_subject <- function(template_path, list_of_source_paths, 
                         flip = FALSE) {
  # import template
  template <- load_emed(template_path, flip = flip)

  # import files to be aligned with template
  list_of_sources <- list()
  for (i in seq_along(list_of_source_paths)) {
    list_of_sources[[i]] <- load_emed(list_of_source_paths[i], flip = flip)
  }
  
  # Register
  sources_reg <- list()
  for (i in seq_along(list_of_sources)) {
    sources_reg[[i]] <- register_emed(list_of_sources[[i]], template)
  }
    
  # Average
  dims <- dim(template)
  mean_fprint <- array(NA, c(dims[1], dims[2], 101))
  for (i in 1:101) {
    arr <- array(NA, c(nrow(sources_reg[[1]]), ncol(sources_reg[[1]]),
                       length(sources_reg)))
    for (j in seq_along(sources_reg)) {
      source_array <- sources_reg[[j]]
      arr[, , j] <- source_array[, , i]
    } 
    mean_fprint[, , i] <- footprint(arr, value = "mean")
  }
  
  # return mean pressure trial
  return(mean_fprint)
}


# =============================================================================

#' Export pressure data frames
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_dfs List. List of pressure data frames
#' @param active_sensors Sensors which are to be analyzed
#' @param pair_set String. For output filename. The group
#' @param add_noise Logical. Adds a small amount of noise to zero values.
#'   Defaults to FALSE
export_df <- function(pressure_dfs, active_sensors, pair_set, 
                      add_noise = FALSE) {
  for (i in seq_along(active_sensors)) {
    df <- data.frame(matrix(NA, nrow = length(pressure_dfs), ncol = 101))
    for (j in seq_along(pressure_dfs)) {
      x <- pressure_dfs[[j]]
      df[j, ] <- x[, active_sensors[i]]
    }
    
    # add noise
    if (add_noise == TRUE) {
      noise <- runif(sum(df == 0), min = 0, max = 0.1)
      df[df == 0] <- noise
    }
  
    # output_filename
    output_filename <- paste0("C:/Users/telfe/Dropbox/My_Projects", 
                              "/Twinfoot/data/sensor_dfs/", pair_set, 
                              "sensor_", i , ".csv")
    
    # write output
    write.table(df, sep = ",", output_filename, row.names = FALSE, col.names = FALSE)
  }
}

#' Export pressure data frames
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param pressure_dfs1 List. List of pressure data frames for twin set 1
#' @param pressure_dfs2 List. List of pressure data frames for twin set 2
#' @param active_sensors Sensors which are to be analyzed
#' @param pair_set String. For output filename. The group
#' @param add_noise Logical. Adds a small amount of noise to those. Defaults to
#'   FALSE
export_df2 <- function(pressure_dfs1, pressure_dfs2, active_sensors, pair_set, 
                      add_noise = FALSE) {
  # set up empty list to store sensor dfs for each twin set, and which columns are all zeros
  sens_df_list1 <- list()
  sens_df_list2 <- list()
  zero_col_list1 <- list()
  zero_col_list2 <- list()
  
  # produce dfs for each active sensor
  for (i in seq_along(active_sensors)) {
    df <- data.frame(matrix(NA, nrow = length(pressure_dfs1), ncol = 101))
    for (j in seq_along(pressure_dfs1)) {
      x <- pressure_dfs1[[j]]
      df[j, ] <- x[, active_sensors[i]]
    }
    
    # check if columns are all zeros
    zeros <- which(colSums(df) == 0)
    if (length(zeros) > 0) {zero_col_list1[[i]] <- zeros}
    
    # add noise
    if (add_noise == TRUE) {
      noise <- runif(sum(df == 0), min = 0, max = 0.1)
      df[df == 0] <- noise
    }
    
    # store
    sens_df_list1[[i]] <- df
  }
  
  for (i in seq_along(active_sensors)) {
    df <- data.frame(matrix(NA, nrow = length(pressure_dfs2), ncol = 101))
    for (j in seq_along(pressure_dfs2)) {
      x <- pressure_dfs2[[j]]
      df[j, ] <- x[, active_sensors[i]]
    }
    
    # check if columns are all zeros
    zeros <- which(colSums(df) == 0)
    if (length(zeros) > 0) {zero_col_list2[[i]] <- zeros}
    
    # add noise
    if (add_noise == TRUE) {
      noise <- runif(sum(df == 0), min = 0, max = 0.1)
      df[df == 0] <- noise
    }
    
    # store
    sens_df_list2[[i]] <- df
  }
  
  
  df <- df[, -zero_cols[[i]]]
  
  # output_filename
  output_filename <- paste0("C:/Users/telfe/Dropbox/My_Projects", 
                            "/Twinfoot/data/sensor_dfs/", pair_set, 
                            "sensor_", i , ".csv")
  
  # write output
  write.table(df, sep = ",", output_filename, row.names = FALSE, col.names = FALSE)
  
  # return zero col list
  return(zero_col_list)
}


# =============================================================================

#' Identify active sensors
#' @param pressure_frames Array. A 3D array covering each timepoint of the measurement.
#' @return Vector. Numbers of all sensors that are active during trial 

act_sens <- function(pressure_frames) {
  sens_fprint <- footprint(pressure_frames, "max")
  active_sensors <- which(sens_fprint > 0)
  return(active_sensors)
}


# =============================================================================

#' SPM results plot function
#' @param spm_data Vector. Results from spm analysis
#' @param active_sensors Vector describing which sensors have non-zero data
#' @param area_length Integer. Length (number of rows) in collection area
#' @param area_width Integer. Width (number of columns) in collection area
#' @param plot_colouring String. "Discrete" or "Cont"
#' @param plot_image Logical. Whether to plot image
#' @return ggplot object
plot_spm_pp <- function(spm_data, active_sensors, area_length, area_width, 
                        plot_colouring = "discrete", plot_image = FALSE) {
  # format
  zstat <- rep(0, times = (area_length * area_width))
  zstat[active_sensors] <- spm_data
  if (is.list(zstat) == TRUE) {zstat <- unlist(zstat)}
  
  x_cor <- seq(from = 0.0025, by = 0.005, length.out = area_width)
  x_cor <- rep(x_cor, each = area_length)
  y_cor <- seq(from = 0.0025 + ((area_length - 1) * 0.005), 
               by = -0.005, length.out = area_length)
  y_cor <- rep(y_cor, times = area_width)
  cor <- cbind(x_cor, y_cor)
  cor <- cbind(x_cor, y_cor, zstat)
  cor <- as.data.frame(cor)
  colnames(cor) <- c("x", "y", "z")
  
  # define colours
  plot_cs <- list(cs_breaks = c(-0.1, 0.0001, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 10.0),
                  cs_cols = c("white", "darkblue", "blue", "light blue", "green",
                              "yellow", "red", "deeppink"))
  
  cols <- c()
  for (i in 1:length(plot_cs[[2]])) {cols[i] = plot_cs[[2]][i]}
  name_cs <- as.character(1:length(plot_cs[[2]]))
  names(cols) <- name_cs
  
  colour <- c()
  for (i in 1:(area_width * area_length)) {
    for (j in 1:length(plot_cs[[2]])) {
      if (cor$z[i] >= plot_cs[[1]][j] & cor$z[i] < plot_cs[[1]][j + 1]) {
        colour = append(colour, j)
      }
    }
  }
  cor <- cbind(cor, colour)
  
  # plot 
  if (plot_colouring == "discrete") {
    g <- ggplot()
    g <- g + geom_raster(data = cor, aes(x = x, y = y, 
                                         fill = as.factor(colour)))
    g <- g + scale_fill_manual(values = cols)
    g <- g + scale_x_continuous(expand = c(0, 0))
    g <- g + scale_y_continuous(expand = c(0, 0))
    g <- g + coord_fixed()
    g <- g + theme(axis.line = element_blank(), axis.text.x = element_blank(),
                   axis.text.y = element_blank(), axis.ticks = element_blank(),
                   axis.title.x = element_blank(), 
                   axis.title.y = element_blank(),
                   panel.background = element_blank(),
                   panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = "none")
  }
  
  # continous colour distribution
  if (plot_colouring == "cont") {
    g <- ggplot(data = cor, aes(x = x, y = y))
    g <- g + geom_raster(aes(fill = z))
    g <- g + coord_fixed()
    g <- g + scale_fill_gradient2(low = "blue", mid = "white", high = "red", breaks = c(0.5, 3, 7))
    g <- g + theme(axis.line = element_blank(), axis.text.x = element_blank(),
                   axis.text.y = element_blank(), axis.ticks = element_blank(),
                   axis.title.x = element_blank(), 
                   axis.title.y = element_blank(),
                   panel.background = element_blank(),
                   panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
  }
  
  # plot if required
  if (plot_image == TRUE) {print(g)}
  
  # return
  return(g)
}


# =============================================================================
# =============================================================================
## End of functions
# =============================================================================
# =============================================================================

## Twinfoot study analysis
# subject ids and twin status
subjects <- c("TWINFOOT001A", "TWINFOOT001B", "TWINFOOT003A", "TWINFOOT003B", 
              "TWINFOOT004A", "TWINFOOT004B", "TWINFOOT005A", "TWINFOOT005B",
              "TWINFOOT006A", "TWINFOOT006B", "TWINFOOT007A", "TWINFOOT007B", 
              "TWINFOOT008A", "TWINFOOT008B", "TWINFOOT009A", "TWINFOOT009B", 
              "TWINFOOT010A", "TWINFOOT010B", "TWINFOOT012A", "TWINFOOT012B")
subjects_dy1 <- c(3, 4, 5, 6, 13, 14, 17, 18, 19, 20)
subjects_mono <- c(1, 2, 7, 8, 9, 10, 11, 12, 15, 16)


# =============================================================================

## Dynamic pressure processing
# generate mean footprints within subject
mean_footprints <- list()
for (i in seq_along(subjects)) {
  subject_folder <- paste0("C:/Users/telfe/Dropbox/My_Projects/Twinfoot/data/", 
                           subjects[i])
  
  ## right foot
  foot_folder <- paste0(subject_folder, "/right")
  list_of_source_paths <- list.files(foot_folder, full.names = TRUE)
  template_path <- list_of_source_paths[1]
  
  # register
  mean_fprint <- mean_subject(template_path, list_of_source_paths, 
                              flip = FALSE)
  plot_footprint(mean_fprint)
  
  # store
  mean_footprints[[((i * 2) - 1)]] <- mean_fprint
  
  ## left foot
  foot_folder <- paste0(subject_folder, "/left")
  list_of_source_paths <- list.files(foot_folder, full.names = TRUE)
  template_path <- list_of_source_paths[1]
  
  # register
  mean_fprint <- mean_subject(template_path, list_of_source_paths, 
                              flip = TRUE)
  plot_footprint(mean_fprint)
  
  # store
  mean_footprints[[(i * 2)]] <- mean_fprint
}

# generate registered footprints between subjects
mean_footprints_registered <- list()
template <- mean_footprints[[1]]
for (i in seq_along(mean_footprints)) {
  mean_footprints_registered[[i]] <- register_emed(mean_footprints[[i]], 
                                                   template)
}

# convert matrix lists to data frames
pressure_dfs <- list()
for (i in seq_along(mean_footprints_registered)) {
  x <- mean_footprints_registered[[i]]
  dims <- dim(x)
  pressure_df <- data.frame(matrix(NA, nrow = 101, ncol = (dims[1] * dims[2])))
  counter <- 1
  for (j in 1:dims[2]) {
    for (k in 1:dims[1]) {
      pressure_df[, counter] <- as.vector(x[k, j, ])
      counter <- counter + 1
    }
  }
  pressure_dfs[[i]] <- pressure_df
}

# by group
dy1 <- pressure_dfs[c(5, 6, 9, 10, 25, 26, 33, 34, 37, 38)]
dy2 <- pressure_dfs[c(7, 8, 11, 12, 27, 28, 35, 36, 39, 40)]
mz1 <- pressure_dfs[c(1, 2, 13, 14, 17, 18, 21, 22, 29, 30)]
mz2 <- pressure_dfs[c(3, 4, 15, 16, 19, 20, 23, 24, 31, 32)]


## export for analysis using SPM
# identify active sensors
active_sensors <- act_sens(mean_footprints_registered[[1]])

# dy pair 1
export_df(dy1, active_sensors, "dy1", add_noise = TRUE)

# dy pair 2
export_df(dy2, active_sensors, "dy2", add_noise = TRUE)

# mz pair 1
export_df(mz1, active_sensors, "mz1", add_noise = TRUE)

# mz pair 2
export_df(mz2, active_sensors, "mz2", add_noise = TRUE)


# =============================================================================

# static pressure processing
# import data
static_files <- list()
for (subject in seq_along(subjects)) {
  subject_folder <- paste0("C:/Users/telfe/Dropbox/My_Projects/Twinfoot/data/", 
                           subjects[subject])
  static_files[[subject]] <- load_emed(paste0(subject_folder, "/static.lst"))
}

# rotate left pointing steps to ensure all face the same way
tbr <- c(3, 4, 5, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18)
for (i in seq_along(tbr)) {
  zz <- static_files[[tbr[i]]]
  zz <- zz[, dim(zz)[2]:1, ]
  zz <- zz[dim(zz)[1]:1, , ]
  static_files[[tbr[i]]] <- zz
}

# select footprints
static_footprints <- list()
for (subject in seq_along(subjects)) {
  # right foot
  #right <- select_region(static_files[[subject]], image_value = "max")
  right <- static_files[[subject]]
  right <- right[1:round(nrow(right) / 2),,]
  
  # left foot
  #left <- select_region(static_files[[subject]], image_value = "max")
  left <- static_files[[subject]]
  left <- left[(round(nrow(left) / 2):nrow(left)),,]
  
  # store
  static_footprints[[(subject * 2) -1]] <- right
  static_footprints[[(subject * 2)]] <- left
}

# select second half of trial footprint (also flip left)
static_footprints2 <- list()
for (i in seq_along(static_footprints)) {
  x <- static_footprints[[i]]
  dims <- dim(x)
  x <- x[,, 10:(dims[3] - 10)]
  if (i %% 2 == 0) {
    static_footprints2[[i]] <- x
  } else {
    static_footprints2[[i]] <- x[c(nrow(x):1),,]
  }
}

# generate registered footprints between subjects
static_footprints_registered <- list()
template <- static_footprints2[[1]]
for (i in seq_along(static_footprints)) {
  static_footprints_registered[[i]] <- register_emed(static_footprints2[[i]], 
                                                     template, 
                                                     match_type = "mean")
}

# convert array data to data frames, with columns representing sensors
static_pressure_dfs <- list()
for (i in seq_along(static_footprints_registered)) {
  x <- static_footprints_registered[[i]]
  dims <- dim(x)
  pressure_df <- rep(NA, times = 101)
  for (j in 1:dims[2]) {
    for (k in 1:dims[1]) {
      pressure_df <- cbind(pressure_df, x[k, j, ])
    }
  }
  static_pressure_dfs[[i]] <- pressure_df[, 2:ncol(pressure_df)]
}

# by group
dy1 <- static_pressure_dfs[c(5, 6, 9, 10, 25, 26, 33, 34, 37, 38)]
dy2 <- static_pressure_dfs[c(7, 8, 11, 12, 27, 28, 35, 36, 39, 40)]
mz1 <- static_pressure_dfs[c(1, 2, 13, 14, 17, 18, 21, 22, 29, 30)]
mz2 <- static_pressure_dfs[c(3, 4, 15, 16, 19, 20, 23, 24, 31, 32)]


## export for analysis using SPM
# identify active sensors
active_sensors <- act_sens(static_footprints_registered[[1]])

# dy pair 1
dy1_zero <- export_df(dy1, active_sensors, "dy1_static", add_noise = FALSE)

# dy pair 2
dy2_zero <- export_df(dy2, active_sensors, "dy2_static", add_noise = FALSE)

# mz pair 1
mz1_zero <- export_df(mz1, active_sensors, "mz1_static", add_noise = FALSE)

# mz pair 2
mz2_zero <- export_df(mz2, active_sensors, "mz2_static", add_noise = FALSE)



# =============================================================================

## Analysis of discrete variables (dynamic analysis)
# load pre-processed data
discrete_data_max <- read_csv("C:/Users/telfe/Dropbox/My_Projects/Twinfoot/max_PP_results.csv")
discrete_data_mean <- read_csv("C:/Users/telfe/Dropbox/My_Projects/Twinfoot/mean_PP_results.csv")

# produce bland altman plots for all regions
regions <- c("Hindfoot", "Midfoot", "MH1", "MH2", "MH3", "MH4", "MH5", 
             "Hallux", "Toe2", "Toe3-5")

# data frame for results
discrete_results <- data.frame(Region = regions, Mean_dy = NA, Mean_mz = NA, 
                               Max_dy = NA, Max_mz = NA)

# Mean pressure dataset
for (i in seq_along(regions)) {
  # format data frames
  df_mean <- discrete_data_mean %>% select(Subject, regions[i])
  dy1_mean <- df_mean[c(5, 6, 9, 10, 25, 26, 33, 34, 37, 38), 2]
  dy2_mean <- df_mean[c(7, 8, 11, 12, 27, 28, 35, 36, 39, 40), 2]
  mz1_mean <- df_mean[c(1, 2, 13, 14, 17, 18, 21, 22, 29, 30), 2]
  mz2_mean <- df_mean[c(3, 4, 15, 16, 19, 20, 23, 24, 31, 32), 2]
  
  df_max <- discrete_data_max %>% select(Subject, regions[i])
  dy1_max <- df_max[c(5, 6, 9, 10, 25, 26, 33, 34, 37, 38), 2]
  dy2_max <- df_max[c(7, 8, 11, 12, 27, 28, 35, 36, 39, 40), 2]
  mz1_max <- df_max[c(1, 2, 13, 14, 17, 18, 21, 22, 29, 30), 2]
  mz2_max <- df_max[c(3, 4, 15, 16, 19, 20, 23, 24, 31, 32), 2]
  
  # calculate correlations and store
  discrete_results$Mean_dy[i] <- as.vector(cor(dy1_mean[, 1], dy2_mean[, 1]))
  discrete_results$Mean_mz[i] <- as.vector(cor(mz1_mean[, 1], mz2_mean[, 1]))
  discrete_results$Max_dy[i] <- as.vector(cor(dy1_max[, 1], dy2_max[, 1]))
  discrete_results$Max_mz[i] <- as.vector(cor(mz1_max[, 1], mz2_max[, 1]))
}

## Center of Pressure Excurrsion Index (CPEI)
cpei_results <- data.frame(Subjects = rep(subjects, each = 2), 
                           Foot = rep(c("right", "left"), 
                                      times = length(subjects)), 
                           CPEI = NA, SD = NA)
for (i in seq_along(subjects)) {
  subject_folder <- paste0("C:/Users/telfe/Dropbox/My_Projects/Twinfoot/data/", subjects[i])
  
  ## right foot
  foot_folder <- paste0(subject_folder, "/right")
  list_of_source_paths <- list.files(foot_folder, full.names = TRUE)
  
  # calculate cpei
  cpei_subject_results <- c()
  for (j in seq_along(list_of_source_paths)) {
    x <- load_emed(list_of_source_paths[[j]])
    cpei_subject_results <- c(cpei_subject_results, 
                              cpei(x, side = "right", plot_result = TRUE))
  }
  
  # store results
  cpei_results[i * 2 - 1, 3] <- mean(cpei_subject_results)
  cpei_results[i * 2 - 1, 4] <- sd(cpei_subject_results)
  
  ## left foot
  foot_folder <- paste0(subject_folder, "/left")
  list_of_source_paths <- list.files(foot_folder, full.names = TRUE)
  
  # calculate cpei
  cpei_subject_results <- c()
  for (j in seq_along(list_of_source_paths)) {
    x <- load_emed(list_of_source_paths[[j]])
    cpei_subject_results <- c(cpei_subject_results, 
                              cpei(x, side = "left", plot_result = TRUE))
  }
  
  # store results
  cpei_results[i * 2, 3] <- mean(cpei_subject_results)
  cpei_results[i * 2, 4] <- sd(cpei_subject_results)
}

# correlation results
dy1_cpei <- cpei_results[c(5, 6, 9, 10, 25, 26, 33, 34, 37, 38), 3]
dy2_cpei <- cpei_results[c(7, 8, 11, 12, 27, 28, 35, 36, 39, 40), 3]
mz1_cpei <- cpei_results[c(1, 2, 13, 14, 17, 18, 21, 22, 29, 30), 3]
mz2_cpei <- cpei_results[c(3, 4, 15, 16, 19, 20, 23, 24, 31, 32), 3]

cor_dy <- as.vector(cor(dy1_cpei, dy2_cpei))
cor_mz <- as.vector(cor(mz1_cpei, mz2_cpei))

cbind(dy1_cpei, dy2_cpei)
cbind(mz1_cpei, mz2_cpei)

