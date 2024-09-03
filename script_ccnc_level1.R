##### SCRIPT: CCNC Level 1 #####
#' @author Christopher Rapp
#' @description
#'
#'
# ---------------------------------------------------------------------------- #
##### SECTION: Import, Export, Functions #####
#'

{
  #' Clear global environment
  rm(list = ls())
  gc(verbose = F)

  # Libraries necessary
  library(data.table)
  library(lubridate)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(viridis)

  # Set working directory
  setwd("~/Library/CloudStorage/Box-Box/TAMU Chamber Experiments")
  work.dir <- getwd()

  tz.local = "US/Central"

  #' @import
  #' Specify import directory
  import.ccnc = paste0(work.dir, "/CCNC/export/level0")

  #' @export
  #' Specify where to send plots and text files
  export.data = paste0(work.dir, "/CCNC/export/level1/")
  export.plot = paste0(work.dir, "/CCNC/export/plots/")
}

# ---------------------------------------------------------------------------- #
##### SECTION: File Selection #####
#'
#'

{
  {
    # List all .csv files
    ccn.files <- list.files(path = import.ccnc,
                            recursive = TRUE,
                            full.names = TRUE,
                            pattern = '.csv')

    ccn.dates = str_extract(ccn.files, pattern = "\\d{4}-\\d{2}-\\d{2}")
    ccn.dates = as.POSIXct(ccn.dates)
  }
}

# ---------------------------------------------------------------------------- #
##### SECTION: Level 1 #####
#'

{
  for (n in 1:length(ccn.files)){

    # ------------------------------------------------------------------------ #
    ##### SUBSECTION: CCNC DATA #####
    #'

    {
      {
        # Change timezone of the local string back to local time
        rawCCN.df <- fread(ccn.files[n]) %>%
          mutate(`Time Local` = lubridate::with_tz(`Time Local`, tzone = tz.local))

        # Remove any R generated filler columns
        remove.c <- str_which(colnames(rawCCN.df), "V\\d{1,2}")
        if (length(remove.c) > 0){
          rawCCN.df <- rawCCN.df %>%
            select(!all_of(remove.c))
        }

        # Remove old time column
        # Remove duplicate time values which for whatever reason appear
        rawCCN.df <- rawCCN.df %>%
          select(!`Time`) %>%
          distinct(`Time Local`, .keep_all = T)

        # Unique date string to use later
        date.c <- unique(lubridate::as_date(rawCCN.df$`Time Local`, tz = tz.local))

        print(date.c)
      }

      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: Column A Processing #####
      #'
      #'

      {
        # Subset Column A data
        data <- rawCCN.df %>%
          select(`Time Local`, ends_with("A"))

        # Create a flow ratio variable and remove data that corresponds to data beyond the 95% range of standard 10:1
        data <- data %>%
          mutate(`Flow Ratio A` = `Sheath Flow A`/`Sample Flow A`) %>%
          filter(`Flow Ratio A` >= 10*0.95 & `Flow Ratio A` <= 10*1.05)

        # Optical parameters
        data <- data %>%
          filter(`1st Stage Mon A` <= 0.8)

        # Change temperature flag to a factor
        data$`Temps Stabilized A` <- factor(ifelse(data$`Temps Stabilized A` == 1, T, F))

        # Remove duplicate time stamps
        data <- data %>%
          distinct(., `Time Local`, .keep_all = T)

        # Run length encoding to group by super saturation
        y = as.numeric(unlist((rle(data$`Current SS A`)[1])))
        x = cumsum(y)

        # This is where SS values shift
        switch.indices = append(0, x[1:(length(x)-1)])

        # Units on this percent change between samples
        # i.e. 0.05% growth or decay between samples for temperature readings
        stable.nm = 0.05

        A.ls <- list()
        for (k in 1:length(switch.indices)){

          if (k == length(switch.indices)){
            end = nrow(data)
          } else {
            end = switch.indices[k+1]
          }

          # Select data within range of indices to select SS ratio
          tmp.df = data[(switch.indices[k] + 1):end, ]

          if (nrow(tmp.df) < 5*60){
            next
          }

          # Calculate the rate of change for the thermocouple values
          # Default is in percentage change from the previous value
          tmp.df$T1.change = abs(collapse::fgrowth(tmp.df$`T1 Read A`))
          tmp.df$T2.change = abs(collapse::fgrowth(tmp.df$`T2 Read A`))
          tmp.df$T3.change = abs(collapse::fgrowth(tmp.df$`T3 Read A`))

          # Smooth the values
          tmp1 = frollmean(tmp.df$T1.change, 15)
          tmp2 = frollmean(tmp.df$T2.change, 15)
          tmp3 = frollmean(tmp.df$T3.change, 15)

          # Three thermocouples are used in determining the column temperature
          # Find the first instance at which the variability meets the threshold
          keep.c <- first(c(first(which(tmp1 < stable.nm)), first(which(tmp2 < stable.nm)), first(which(tmp3 < stable.nm)))):nrow(tmp.df)

          # Filter out data
          tmp.df <- tmp.df[keep.c, ]

          # Add a sequence ID
          tmp.df <- tmp.df %>%
            mutate(`ID A` = k, .after = "Time Local") %>%
            select(!contains("change"))

          # Export data out
          A.ls[[k]] <- tmp.df
        }
      }

      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: Column B Processing #####
      #'
      #'

      {
        # Subset Column B data
        data <- rawCCN.df %>%
          select(`Time Local`, ends_with("B"))

        # Create a flow ratio variable and remove data that corresponds to data beyond the 95% range of standard 10:1
        data <- data %>%
          mutate(`Flow Ratio B` = `Sheath Flow B`/`Sample Flow B`) %>%
          filter(`Flow Ratio B` >= 10*0.95 & `Flow Ratio B` <= 10*1.05)

        # Optical parameters
        data <- data %>%
          filter(`1st Stage Mon. B` <= 0.8)

        # Change temperature flag to a factor
        data$`Temps Stabilized B` <- factor(ifelse(data$`Temps Stabilized B` == 1, T, F))

        # Remove duplicate time stamps
        data <- data %>%
          distinct(., `Time Local`, .keep_all = T)

        # Run length encoding to group by super saturation
        y = as.numeric(unlist((rle(data$`Current SS B`)[1])))
        x = cumsum(y)

        # This is where SS values shift
        switch.indices = append(0, x[1:(length(x)-1)])

        # Units on this percent change between samples
        # i.e. 0.05% growth or decay between samples for temperature readings
        stable.nm = 0.05

        B.ls <- list()
        for (k in 1:length(switch.indices)){

          if (k == length(switch.indices)){
            end = nrow(data)
          } else {
            end = switch.indices[k+1]
          }

          # Select data within range of indices to select SS ratio
          tmp.df = data[(switch.indices[k] + 1):end, ]

          if (nrow(tmp.df) < 5*60){
            next
          }

          # Calculate the rate of change for the thermocouple values
          # Default is in percentage change from the previous value
          tmp.df$T1.change = abs(collapse::fgrowth(tmp.df$`T1 Read B`))
          tmp.df$T2.change = abs(collapse::fgrowth(tmp.df$`T2 Read B`))
          tmp.df$T3.change = abs(collapse::fgrowth(tmp.df$`T3 Read B`))

          # Smooth the values
          tmp1 = frollmean(tmp.df$T1.change, 15)
          tmp2 = frollmean(tmp.df$T2.change, 15)
          tmp3 = frollmean(tmp.df$T3.change, 15)

          # Three thermocouples are used in determining the column temperature
          # Find the first instance at which the variability meets the threshold
          keep.c <- first(c(first(which(tmp1 < stable.nm)), first(which(tmp2 < stable.nm)), first(which(tmp3 < stable.nm)))):nrow(tmp.df)

          # Filter out data
          tmp.df <- tmp.df[keep.c, ]

          # Add a sequence ID
          tmp.df <- tmp.df %>%
            mutate(`ID B` = k, .after = "Time Local") %>%
            select(!contains("change"))

          # Export data out
          B.ls[[k]] <- tmp.df
        }
      }

      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: Merge Data #####
      #'
      #'

      {
        # Retrieve data from the original dataframe without the column specific variables
        tmp.df <- rawCCN.df %>%
          select(!ends_with("B")) %>%
          select(!ends_with("A"))

        # Remove empty lists (SS settings that were short)
        A.ls <- A.ls[lengths(A.ls) != 0]
        B.ls <- B.ls[lengths(B.ls) != 0]

        if (length(A.ls) != 0){

          # Apply temperature stabilization DMT uses in addition to the filtering used above
          A.df <- rbindlist(A.ls) %>%
            filter(`Temps Stabilized A` == T) %>%
            filter(`CCN Number Conc A` >= 0) %>%
            filter(!is.na(`Current SS A`))

        } else {
          A.df = data.frame("Time Local" = NA, check.names = F)
        }

        if (length(B.ls) != 0){

          # Apply temperature stabilization DMT uses in addition to the filtering used above
          B.df <- rbindlist(B.ls) %>%
            filter(`Temps Stabilized B` == T) %>%
            filter(`CCN Number Conc B` >= 0) %>%
            filter(!is.na(`Current SS B`))

        } else {
          B.df = data.frame("Time Local" = NA, check.names = F)
        }

        # Merge data
        dataCCN.df <- left_join(tmp.df, A.df, by = "Time Local")
        dataCCN.df <- left_join(dataCCN.df, B.df, by = "Time Local")
      }

      # ---------------------------------------------------------------------- #
      ##### SECTION: Export Data #####
      #'
      {
        export.filename = paste0(export.data, date.c, "_level1.csv")

        print(export.filename)

        data.table::fwrite(dataCCN.df, file = export.filename, showProgress = T)

        gc(verbose = F)
      }

      # ----------------------------------------------------------------------- #
      ##### SECTION: Plotting #####
      #'

      if (nrow(A.df) > 1 & nrow(B.df) > 1){

        # Create factors instead which are better in plotting
        dataCCN.df <- dataCCN.df %>%
          mutate(`Current SS A` = factor(`Current SS A`)) %>%
          mutate(`Current SS B` = factor(`Current SS B`))

        SS.c = sort(unique(c(unique(dataCCN.df$`Current SS A`), unique(dataCCN.df$`Current SS B`))))

        palette_blues <- colorRampPalette(colors = c("gray", "#004b88"))(length(SS.c))

        colors.df <- data.frame(values = palette_blues, labels = SS.c)

        # Plot daily data
        A.gg <- ggplot(dataCCN.df, aes(`Time Local`, `CCN Number Conc A`, color = `Current SS A`)) +
          geom_point() +
          ggtitle(paste0("CCNC Data Level 1 - ", date.c)) +
          ylab(expression("Number Concentration (n cm"^"-3"*")")) +
          scale_x_datetime(date_breaks = "1 hour") +
          scale_y_continuous(n.breaks = 10) +
          theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(colour = "grey80", linewidth = 0.25),
                panel.grid.minor = element_line(colour = "grey80", linewidth = 0.25),
                panel.border = element_rect(colour = "black", fill = NA),
                axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                legend.title = element_blank()) +
          scale_color_manual(values = colors.df$values, breaks = colors.df$labels)

        # Plot daily data
        B.gg <- ggplot(dataCCN.df, aes(`Time Local`, `CCN Number Conc B`, color = `Current SS B`)) +
          geom_point() +
          xlab("\nTime Local") +
          ylab(expression("Number Concentration (n cm"^"-3"*")")) +
          scale_x_datetime(date_breaks = "1 hour", date_labels = "%H:%M") +
          scale_y_continuous(n.breaks = 10) +
          theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(colour = "grey80", size = 0.25),
                panel.grid.minor = element_line(colour = "grey80", size = 0.25),
                panel.border = element_rect(colour = "black", fill = NA),
                axis.text.x = element_text(angle = 90, vjust = 0.5),
                legend.position = "none") +
          scale_color_manual(values = colors.df$values, breaks = colors.df$labels) +
          guides(color = "none")

        data.gg <- A.gg/B.gg & theme(legend.position = "right")
        data.gg <- data.gg + plot_layout(axis_titles = "collect", guides = 'collect')

        # DECLARE EXPORT.DATA AT BEGINNING OF CODE OTHERWISE THIS WONT WORK
        plot.filename = paste0(export.plot, 'CCNC_level1_', date.c, '.png')

        print(plot.filename)

        # This is a save function for the plot
        ggsave(plot.filename, data.gg, width = 10, height = 6, dpi = 300, units = "in")
      }
    }
  }
}
