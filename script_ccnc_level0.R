##### SCRIPT: CCNC Level 0 #####
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

  # Libraries necessary
  library(data.table)
  library(lubridate)
  library(dplyr)
  library(stringr)

  # Set working directory
  setwd("~/Library/CloudStorage/Box-Box/TAMU Chamber Experiments")
  work.dir <- getwd()

  #' @import
  #' Specify import directory
  import.ccnc = paste0(work.dir, "/CCNC/import/CCN/")

  #' @export
  #' Specify where to send plots and text files
  export.data = paste0(work.dir, "/CCNC/export/level0/")
  export.plot = paste0(work.dir, "/CCNC/export/plots/")
}

# ---------------------------------------------------------------------------- #
##### SECTION: File Selection #####
#' Loop through each year and read in available data
#' Merge data by year into a dataframe to apply operations
#'

{
  # List all .csv files
  ccn.files <- list.files(path = import.ccnc,
                          recursive = TRUE,
                          full.names = TRUE,
                          pattern = '.csv')

  # Identify any erroneous files as a position vector
  remove.c = str_which(ccn.files, "duplicate")

  # Only remove files if the vector is non-zero
  if (length(remove.c) != 0){

    ccn.files <- ccn.files[-remove.c]
  }

  # Use regex to extract the dated portion of the filename
  ccn.filename = str_extract(paste0(ccn.files), pattern = '\\d*(?=.csv)')

  # Extract the dates and times from the file pattern
  ccn.starttime = as.POSIXct(ccn.filename, format = "%y%m%d%H%M%S")

  # Change to location specific timezone
  # Changes automatically based on daylight savings time
  ccn.starttime = lubridate::force_tz(ccn.starttime, tzone = "US/Central")

  # Identify unique dates
  dates.c = unique(format(ccn.starttime, '%y%m%d'))
}

# ---------------------------------------------------------------------------- #
##### SECTION: Level 0 #####
#'

{

  # -------------------------------------------------------------------------- #
  ##### SECTION: Data Wrangling #####
  #'

  {
    data.all.ls <- lapply(dates.c, function(x){

      tmp.files <- ccn.files[str_which(ccn.files, x)]

      print(x)

      # Loop through directory and create a data.frame for each file
      # This is necessary if you are looking at multiple text files of data
      # Skip goes to the keyword in the text file and starts reading there
      data.ls <- lapply(tmp.files, function(x) {

        # Read in data, starting on "Current SS"
        data <- fread(paste0(x), header = TRUE, stringsAsFactors = FALSE,
                      fill = TRUE, check.names = FALSE, sep = ",", skip = "Current SS", strip.white = T)

        # Set empty values to NA
        data[data == ""] <- NA

        # Only keep non-empty rows
        data <- data[rowSums(is.na(data)) == 0, ]

        # Extract filename pattern
        filename = str_extract(paste0(x), pattern = '\\d*(?=.csv)')

        # Extract the dates and times from the file pattern
        tmp.date = lubridate::as_date(as.POSIXct(filename, format = "%y%m%d%H%M%S", tz = "US/Central"), tz = "US/Central")

        # Add date column
        data <- data %>%
          mutate(`Time Local` = as.POSIXct(paste0(tmp.date, " ", data$Time), tz = "US/Central"), .before = everything()) %>%
          mutate(`Time UTC` = with_tz(`Time Local`, tz = "UTC")) %>%
          select(`Time Local`, `Time UTC`, !Time, everything())

        return(data)
      })

      # Merge the list of dataframes
      rawCCN.df <- data.table::rbindlist(data.ls, use.names = T, fill = T)

      # ---------------------------------------------------------------------- #
      ##### SUBSECTION: Bin Data #####
      #' Binned aerosol data
      {
        # Rename columns to include sizing information
        bins.nm = colnames(rawCCN.df)[str_which(colnames(rawCCN.df), "Bin \\d+")]

        # Bin information as decribed by https://dropletmeasure.wpenginepowered.com/wp-content/uploads/2020/02/DOC-0128-Rev-F-Manual-Operator-Dual-CCN.pdf
        bins.c = c(0.75, 0.875, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75)

        # Specify what column
        binsA.c = paste(bins.c, "A")
        binsB.c = paste(bins.c, "B")

        # Merge two strings
        bins.c = c(binsA.c, binsB.c)

        # Rename
        setnames(rawCCN.df, old = bins.nm, new = bins.c)
      }

      # ------------------------------------------------------------------------ #
      ##### SUBSECTION: Total Concentration #####
      #'

      {
        columnA <- rawCCN.df %>%
          select(all_of(binsA.c)) %>%
          rowSums()

        columnB <- rawCCN.df %>%
          select(all_of(binsB.c)) %>%
          rowSums()

        # See DMT manual for specifics but this is how they perform their total calculation
        # Despite this there is still an offset for this vs theirs despite being the same...
        columnA <- ((columnA + rawCCN.df$`Overflow A`)/rawCCN.df$`Sample Flow A`)*60
        columnB <- ((columnB + rawCCN.df$`Overflow B`)/rawCCN.df$`Sample Flow B`)*60

        # Replace columns with these calculated ones
        rawCCN.df <- rawCCN.df %>%
          mutate(`CCN Number Conc A` = columnA) %>%
          mutate(`CCN Number Conc B` = columnB)
      }

      # ------------------------------------------------------------------------ #
      ##### SUBSECTION: Export Data #####
      #'

      {
        export.filename = paste0(export.data, unique(as_date(rawCCN.df$`Time Local`, tz = "US/Central")), "_level0.csv")

        print(export.filename)

        data.table::fwrite(rawCCN.df, file = export.filename, showProgress = T)
      }

      return(rawCCN.df)
    })

    gc(verbose = F)
  }
}
