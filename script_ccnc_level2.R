##### SCRIPT: CCNC Level 2 #####
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
  library(latex2exp)

  # Set working directory
  setwd("~/Library/CloudStorage/Box-Box/TAMU Chamber Experiments")
  work.dir <- getwd()

  tz.local = "US/Central"

  #' @import
  #' Specify import directory
  import.ccnc = paste0(work.dir, "/CCNC/export/level1")
  import.cpc = paste0(work.dir, "/CCNC/import/CPC")
  import.dma = paste0(work.dir, "/CCNC/import/DMA")

  #' @export
  #' Specify where to send plots and text files
  export.data = paste0(work.dir, "/CCNC/export/level2/")
  export.plot = paste0(work.dir, "/CCNC/export/plots/")

  #' Function scripts
  source("~/Library/CloudStorage/Box-Box/TAMU Chamber Experiments/scripts/Functions/functions_aerosols.R")
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

  {
    # List all .csv files
    cpc.files <- list.files(path = import.cpc,
                            recursive = TRUE,
                            full.names = TRUE,
                            pattern = '.txt')

    cpc.dates = str_extract(cpc.files, pattern = "(?<=MCPC_)\\d{6}")
    cpc.dates = as.POSIXct(cpc.dates, format = "%y%m%d")
  }

  {
    dma.files <- list.files(path = import.dma,
                            recursive = TRUE,
                            full.names = TRUE)

    dma.dates = str_extract(dma.files, pattern = "(?<=DMA/)\\d{6}")
    dma.dates = as.POSIXct(dma.dates, format = "%y%m%d")
  }

  # Find the intersection of SEMS, PCU, and SPIN data
  dates.c <- Reduce(intersect, list(as.character(ccn.dates), as.character(cpc.dates)))
}

# ---------------------------------------------------------------------------- #
##### SECTION: Level 2 #####
#'

for (n in 1:length(dates.c)){

  date.c = dates.c[n]

  # ------------------------------------------------------------------------ #
  ##### SUBSECTION: CCNC Data #####
  #'

  {
    {
      # Match file to the dates
      tmp <- ccn.files[which(as.character(ccn.dates) %in% date.c)]

      # Change timezone of the local string back to local time
      rawCCN.df <- fread(tmp) %>%
        mutate(`Time Local` = lubridate::with_tz(`Time Local`, tzone = tz.local))
    }
  }

  # ------------------------------------------------------------------------ #
  ##### SUBSECTION: CPC Data #####
  #'

  {
    # Match file to the dates
    tmp <- cpc.files[which(as.character(cpc.dates) %in% date.c)]

    cpc.ls <- lapply(tmp, function(x){

      # Read in data
      tmp <- fread(x, skip = "#YY/MM/DD")

      # Create time object
      tmp <- tmp %>%
        mutate(`Time Local` = as.POSIXct(paste0(tmp$`#YY/MM/DD`, " ", tmp$`HR:MN:SC`), format = "%y/%m/%d %H:%M:%S")) %>%
        mutate(`Time UTC` = with_tz(`Time Local`, tzone = "UTC")) %>%
        select(!`#YY/MM/DD`, !`HR:MN:SC`) %>%
        select(`Time Local`, `Time UTC`, everything())

      return(tmp)
    })

    rawCPC.df <- rbindlist(cpc.ls) %>%
      select(`Time Local`, concent) %>%
      rename(`CPC Concentration` = concent) %>%
      distinct(`Time Local`, .keep_all = T)
  }

  # ------------------------------------------------------------------------ #
  ##### SUBSECTION: DMA Data #####
  #'

  {
    # Match file to the dates
    tmp <- dma.files[which(as.character(dma.dates) %in% date.c)]

    # Some days did not have a DMA measurement
    # They default to 100 nm
    if (length(tmp) != 0){

      dma.ls <- lapply(tmp, function(x){

        # Read in data
        tmp <- fread(x, skip = "#YY/MM/DD", fill = T)

        if (ncol(tmp) == 38){

          # Create time object
          tmp <- tmp %>%
            mutate(`Time Local` = as.POSIXct(paste0(tmp$`#YY/MM/DD`, " ", tmp$`HR:MN:SC`), format = "%y/%m/%d %H:%M:%S")) %>%
            mutate(`Time UTC` = with_tz(`Time Local`, tzone = "UTC")) %>%
            select(!c(`#YY/MM/DD`, `HR:MN:SC`)) %>%
            select(`Time Local`, `Time UTC`, everything())

          return(tmp)
        }
      })

      rawDMA.df <- rbindlist(dma.ls) %>%
        select(`Time Local`, mono_size) %>%
        rename(`Diameter` = mono_size) %>%
        distinct(`Time Local`, .keep_all = T)

      data.check = T

    } else {

      data.check = F

      Diameter = 100 # nm
    }
  }

  # ------------------------------------------------------------------------ #
  ##### SUBSECTION: Merge Data #####
  #'

  {
    tmp.df <- left_join(rawCCN.df, rawCPC.df, by = "Time Local")

    # Merge CPC and DMA data
    if (data.check) {

      tmp.df <- left_join(tmp.df, rawDMA.df, by = "Time Local")

    } else {

      tmp.df <- tmp.df %>%
        mutate(`Diameter` = Diameter)
    }

    # Apply filtering
    dataCCN.df <- tmp.df %>%
      filter(!is.na(`CPC Concentration`)) %>%
      filter(!is.na(`Current SS A`) & !is.na(`Current SS B`)) %>%
      filter(!is.na(`CCN Number Conc A`) & !is.na(`CCN Number Conc B`)) %>%
      filter(!is.na(`Diameter`)) %>%
      select(contains("Time"), `CPC Concentration`, `Diameter`,
             contains("CCN"),
             contains("Current SS"),
             everything())
  }

  # ------------------------------------------------------------------------ #
  ##### SECTION: Statistical Analysis #####
  #'

  if (nrow(dataCCN.df != 0)){

    # Calculate activated fraction
    dataCCN.df <- dataCCN.df %>%
      mutate(`Activated Fraction A` = `CCN Number Conc A`/`CPC Concentration`) %>%
      mutate(`Activated Fraction B` = `CCN Number Conc B`/`CPC Concentration`) %>%
      mutate(`Activated Fraction A` = if_else(`Activated Fraction A` >= 1, NA, `Activated Fraction A`)) %>%
      mutate(`Activated Fraction B` = if_else(`Activated Fraction B` >= 1, NA, `Activated Fraction B`)) %>%
      mutate(`Current SS A` = as.factor(`Current SS A`)) %>%
      mutate(`Current SS B` = as.factor(`Current SS B`))

    {
      dataA.df <- dataCCN.df %>%
        select(`Time Local`, `Diameter`, `CPC Concentration`, ends_with("A")) %>%
        group_by(`ID A`) %>%
        summarize(`Column` = "A",
                  `CPC Mean` = mean(`CPC Concentration`, na.rm = T),
                  `CCN Mean` = mean(`CCN Number Conc A`, na.rm = T),
                  `Activated_Fraction` = `CCN Mean`/`CPC Mean`,
                  `SS` = unique(`Current SS A`),
                  `Samples` = n(),
                  `SD` = sd(`Activated Fraction A`, na.rm = T),
                  `Temperature` = mean(`T2 Read A`, na.rm = T),
                  `Length Diameter` = length(unique(`Diameter`)),
                  `Diameter` = mean(`Diameter`, na.rm = T)) %>%
        select(!`ID A`) %>%
        filter(`Length Diameter` == 1)

      dataB.df <- dataCCN.df %>%
        select(`Time Local`, `Diameter`, `CPC Concentration`, ends_with("B")) %>%
        group_by(`ID B`) %>%
        summarize(`Column` = "B",
                  `CPC Mean` = mean(`CPC Concentration`, na.rm = T),
                  `CCN Mean` = mean(`CCN Number Conc B`, na.rm = T),
                  `Activated_Fraction` = `CCN Mean`/`CPC Mean`,
                  `SS` = unique(`Current SS B`),
                  `Samples` = n(),
                  `SD` = sd(`Activated Fraction B`, na.rm = T),
                  `Temperature` = mean(`T2 Read B`, na.rm = T),
                  `Length Diameter` = length(unique(`Diameter`)),
                  `Diameter` = mean(`Diameter`, na.rm = T)) %>%
        select(!`ID B`) %>%
        filter(`Length Diameter` == 1)

      dataALL.df <- rbind(dataA.df, dataB.df) %>%
        mutate(`SS` = as.numeric(as.character(SS))) %>%
        filter(!is.na(`Activated_Fraction`))

      dataALL.ls <- dataALL.df %>%
        group_split(`Diameter`)

      for (i in 1:length(dataALL.ls)){

        #tmp.df$Activated_Fraction <- scales::rescale(tmp.df$Activated_Fraction, to = c(0, 1))

        Dp = unique(tmp.df$Diameter)

        # Skip loops with no ramp (indicating experiment never started)
        {
          fit.break <<- FALSE

          tryCatch(
            nls(Activated_Fraction ~ SSlogis(SS, Asym, xmid, scal), data = tmp.df),
            error = function(e) {fit.break <<- TRUE})

          # Stop code from continuing if there are no ramps detected which means the dataset will be empty and continue to throw errors
          if (fit.break) {break}
        }

        fit = nls(Activated_Fraction ~ SSlogis(SS, Asym, xmid, scal), data = tmp.df, control = list(maxiter = 100, tol = 1e-7))
        newdata = data.frame(SS = seq(0, 1, length.out = 1000))
        predict <- predict(fit, newdata, interval = c("prediction"))

        newdata$Fit = as.numeric(predict)

        sc = newdata$SS[first(which(newdata$Fit >= fit$m$getAllPars()[1]/2))]
        Sc = sc/100 + 1

        kappa = round(calc.kappa(Sc, Dp/1000, mean(tmp.df$Temperature) + 273.15), 3)

        title = date.c
        subtitle = paste0("$\\kappa$", " = ", kappa, ", ", "$\\D_p$", " = ", Dp , " nm")

        kappa.gg = ggplot(tmp.df, aes(x = `SS`, y = `Activated_Fraction`, col = Column)) +
          geom_point() +
          geom_errorbar(aes(ymin = `Activated_Fraction` - SD, ymax = `Activated_Fraction` + SD)) +
          geom_line(data = newdata, aes(x = `SS`, y = `Fit`), inherit.aes = F) +
          scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
          scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
          labs(title = title, subtitle = TeX(subtitle)) +
          ylab("Activated Fraction") +
          xlab("Supersaturation (%)") +
          theme_minimal() +
          annotate("pointrange", x = fit$m$getAllPars()[2], y = fit$m$getAllPars()[1]/2, xmin = sc - 0.05, xmax = sc + 0.05)

        ggsave(filename = paste0(export.plot, date.c, " ", Dp, " Hygroscopicity.png"),
               plot = kappa.gg, width = 6, height = 6, bg = "white")

        print(paste0("Plotting ", date.c))
      }
    }
  }
}

calc.Sc(0.53, 100/1000, mean(tmp.df$Temperature) + 273.15, T)
