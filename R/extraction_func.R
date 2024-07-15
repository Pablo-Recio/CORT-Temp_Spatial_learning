# func.R - a script for custom functions

#-------------------
# bookkeeping
#-------------------
pacman::p_load(tidyverse, here, brms, rstan, flowCore)
#-------------------
# statistics
#-------------------
#----
#' @title geo.mean
#' @param x a vector
#' @description calculates the geometric mean a robust metric of the central tendency of a population. applicable to values on a logicle scale. is robust to outliers. will only ever be at most equal to the arithmetic mean
geo.mean <- function(x, na.rm = T){
  10^base::mean(log.sc(x), na.rm = na.rm)
}
#----
#' @title se
#' @param x a numeric vector
#' @description calculates the standard error
se <- function(x, na.rm = T){
  stats::sd(x, na.rm = na.rm) / base::sqrt(base::length(x))
}
#----
#' @title r.sd
#' @param x a numeric vector
#' @description calculates the robust standard deviation; based on median absolute deviation.
r.sd <- function(x, na.rm = T) {
  1.4826 * stats::mad(x, constant = 1.4826, na.rm = na.rm)
}
#----
#' @title r.se
#' @param x a numeric vector
#' @description calculates the robust standard error
r.se <- function(x, na.rm = T) {
  r.se = r.sd(x, na.rm = na.rm) / sqrt(n())
}
#----
#' @title log.sc
#' @param x a numeric vector
#' @description scales a vector to logicle. allows for negative or zero values
log.sc <- function(x) {
  base::suppressWarnings(base::ifelse(x <= 0,
                                      -base::log10(base::abs(x) + 1),
                                      base::log10(x + 1)))
}
#----
#' @title pmcmc
#' @param x The vector for the posterior distribution. Note that this will test the null hypothesis that the parameter of interest is significantly different from 0.
#' @param null A numeric value decsribing what the null hypothesis should be
#' @param twotail Whether to conduct a one-tailed hypothesis or a two-tailed hypotheses. Default = true indicating a two-tailed test will be done.
#' @description calculates the pmcmc value for a given posterior distribution (or any combination of posterior distributions) as a metric of significance.
pmcmc <- function(x, null = 0, twotail = TRUE){
  if(twotail){
    2*(1 - max(table(x<=null) / length(x)))
  } else{
    (1 - max(table(x<=null) / length(x)))
  }
}
#----

#-------------------
# data processing
#-------------------
#----
#' @title read.it
#' @param path a string object directing to a folder
#' @param pattern a string with a file extension. defaults to .csv
#' @description reads in all files within a folder to a list. id's are assigned based on file name. currently supports .csv, .exp, and .fcs file extensions.
read.it <- function(path, pattern = ".csv") {
  data <- list()
  files <- list.files(path = path, pattern = pattern, full.names = TRUE)
  names <- list.files(path = path, pattern = pattern, full.names = FALSE)
  for(i in 1:length(files)) {
    if(pattern == ".csv"){
      x <- readr::read_csv(files[i])
    }
    if(pattern == ".exp"){
      x <- metabR::read.exp(files[i])  
    }
    if(pattern == ".fcs"){
      x <- extract.fcs(files[i])
    }
    data[[i]] <- x
    names(data)[[i]] <- names[i]
  } 
  if(pattern != ".fcs"){
    data <- plyr::ldply(data)
  }
  return(data)
}
#----
#' @title re.code
#' @param x a character or factor vector
#' @param mapping a vector of matched channel names to label names
#' @description function for quickly re-labelling columns with new names
re.code <- function(x, mapping) {
  labels <- base::unname(mapping)
  base::names(labels) <- base::names(mapping)
  out <- labels[base::match(x, base::names(mapping))]
  out[base::is.na(out)] <- x[base::is.na(out)]
  return(out)
}
#----
#' @title read.fcs
#' @description simple wrapper from flowCore that reads a .fcs file
#' @param path a character string that indicates a path to a .fcs file
read.fcs <- function(path){
  fcs <- flowCore::read.FCS(path, transformation = F)
  return(fcs)
}
#----
#' @title data.fcs
#' @description a tidyr pipeline for standardizing .fcs data
#' @param fcs a .fcs file read in using read.fcs
#' @return a dataframe with standardized column names and a row for each observation
data.fcs <- function(fcs){
  data <- fcs@exprs %>%
    data.frame() %>%
    # apply standardized naming convention
    dplyr::rename_with(.cols = tidyselect::everything(),
                       .fn = stringr::str_to_lower) %>%
    # remove the time column (if it exists)
    ## [temporary - could be useful later, but right now causes issues]
    dplyr::select(-time) %>%
    # give each row (cell) an identifier number for pivoting
    dplyr::mutate(event = dplyr::row_number()) %>%
    # pivot into long form
    tidyr::pivot_longer(-event,
                        names_to = "ch",
                        values_to = "val") %>%
    # logicle transform the value column
    dplyr::mutate(val = log.sc(val)) %>%
    # separate the channel and parameter columns
    dplyr::mutate(param = stringr::str_extract(ch, ".[a|h|w]$"),
                  param = stringr::str_remove(param, "\\."),
                  ch = stringr::str_remove(ch, ".[a|h|w]$")) %>%
    # arrange columns into standard order
    dplyr::select(event, ch, param, val) %>%
    # arrange rows into standard order
    dplyr::arrange(event, ch, param, val)
  return(data)
}
#----
#' @title process.fcs
#' @description a tidyr pipeline for setting a standard processing protocol on flow cytometry data. log transforms all values, filters offscale events, then gates to most "average" cells in the population - within 1 s.d. of the median in all channels simultaneously.
#' @param data a dataframe
#' @return a dataframe
process.fcs <- function(data){
  data %>% 
    # pivot to wide form
    tidyr::pivot_wider(names_from = "ch",
                       values_from = "val",
                       id_cols = c(event, param)) %>%
    # gate out any off-scale events
    dplyr::filter(dplyr::if_all(.cols = -c(event, param),
                                .fns = ~.x >= 0 & .x <= 5)) %>%
    # gate to within 1 robust s.d. of the median on all channels
    group_by(param) %>%
    dplyr::filter(dplyr::if_all(.cols = -c(event),
                                .fns = ~.x <= stats::median(.x) + r.sd(.x) &
                                  .x >= stats::median(.x) - r.sd(.x))) %>%
    ungroup() %>%
    # pivot back to long form
    tidyr::pivot_longer(cols = -c(event, param),
                        names_to = "ch",
                        values_to = "val")
}
#----
#' @title summary.fcs
#' @param data a dataframe of fcs data read in using read.fcs + data.fcs
#' @param na.rm a logical of whether or not to remove missing values from calculations. defaults to TRUE.
#' @description a tidyverse wrapper for calculating common summary statistics for flow cytometry data. includes number of cells (n), min, median, max, geometric mean (geo.mean), arithmetic mean (mean), standard deviation (sd), standard error (se), robust standard deviation (r.sd), and robust standard error (r.se).
summary.fcs <- function(data, na.rm = TRUE){
  data %>%
    dplyr::group_by(ch, param) %>%
    dplyr::reframe(n = dplyr::n(),
                   min = base::min(val, na.rm = na.rm),
                   median = stats::median(val, na.rm = na.rm),
                   max = base::max(val, na.rm = na.rm),
                   # if data is already log transformed, the mean is the geometric mean
                   geo.mean = base::mean(val, na.rm = na.rm),
                   arith.mean = log10(mean(10^val, na.rm = na.rm)),
                   sd = stats::sd(val, na.rm = na.rm),
                   se = se(val, na.rm = na.rm),
                   r.sd = r.sd(val, na.rm = na.rm),
                   r.se = r.se(val, na.rm = na.rm))
}
#----
#' @title meta.fcs
#' @param fcs an fcs file that has been read in with read.fcs
#' @description extracts metadata on a .fcs file into a user-friendly format
meta.fcs <- function(fcs){
  meta <- fcs@description %>%
    as_tibble() %>%
    slice(1) %>%
    select(`FILENAME` , `FCSversion`, `$CYT`, 
           `CREATOR`, `EXPERIMENT NAME`, `$SRC`, 
           `TUBE NAME`, `$FIL`, `$TOT`, 
           `$DATE`, `$BTIM`, `$ETIM`) %>%
    rename(file = `FILENAME`,
           fcs.version = `FCSversion`,
           cytometer = `$CYT`,
           source = `CREATOR`,
           experiment = `EXPERIMENT NAME`,
           specimen = `$SRC`,
           tube = `TUBE NAME`,
           export.id = `$FIL`,
           events = `$TOT`,
           date = `$DATE`,
           start = `$BTIM`,
           end = `$ETIM`) %>%
    mutate(date = as_date(dmy(date)),
           start = hms(start),
           end = hms(end),
           time = end - start) %>%
    mutate(across(.cols = everything(),
                  .fns = as.character)) %>%
    pivot_longer(everything(),
                 names_to = "info",
                 values_to = "description")
  return(meta)
}
#----
#' @title voltage.fcs
#' @param fcs an fcs file that has been read in with read.fcs
#' @description extracts fluorescent channel and voltage information from a .fcs file
voltage.fcs <- function(fcs){
  voltage <- fcs@description %>%
    as_tibble() %>%
    slice(1) %>%
    select(grep("^\\$P[0-9]{1,2}N$", colnames(.)),
           grep("^\\$P[0-9]{1,2}V$", colnames(.))) %>%
    pivot_longer(cols = everything()) %>%
    mutate(param = str_extract(name, "[N|S|V]$"),
           name = str_remove(name, "[N|S|V]$")) %>%
    pivot_wider(names_from = "param",
                values_from = "value") %>%
    mutate(N = str_remove(N, "\\-[A|H|W]")) %>%
    select(N, V) %>%
    rename(ch = N,
           voltage = V) %>%
    group_by(ch) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(ch = str_to_lower(ch),
           ch = str_replace_all(ch, "(\\-| )", "\\.")) %>%
    dplyr::filter(ch != "time")
  return(voltage)
}
#----
#' @title lasers.fcs
#' @param fcs an fcs file that has been read in with read.fcs
#' @description extracts information on laser delays from a .fcs file into a user-friendly format
lasers.fcs <- function(fcs){
  lasers <- fcs@description %>%
    as_tibble() %>%
    slice(1) %>%
    select(`$CYT`, `FSC ASF`, grep("^LASER[0-9]{1,2}", colnames(.))) %>%
    pivot_longer(cols = -c(`$CYT`)) %>%
    mutate(param = str_extract(name, "(NAME|DELAY|ASF)$"),
           laser = str_extract(name, "([0-9]{1,2}|FSC)")) %>%
    select(-name) %>%
    pivot_wider(names_from = "param",
                values_from = "value") %>%
    rename(cytometer = `$CYT`,
           laser = laser,
           spectra = NAME,
           delay = DELAY,
           asf = ASF) %>%
    mutate(spectra = ifelse(laser == "FSC", "FSC", spectra),
           spectra = str_to_lower(spectra)) %>%
    select(cytometer, laser, spectra, delay, asf)
  return(lasers)
}
#----
#' @title extract.fcs
#' @param path a character string denoting a path to .fcs file
#' @description a wholistic pipeline for extracting flow cytometry data from a .fcs file into a user-friendly list in r. returns metadata on the file, voltage, and laser settings, the raw, gated, and summarised fluorescence data, and some basic plots of the raw and gated data. if you only care about extracting the final summary data (i.e., for analysis), see `read.it()`
extract.fcs <- function(path){
  # read the fcs file
  fcs <- read.fcs(path)
  # extract metadata
  meta <- list(info = meta.fcs(fcs),
               voltage = voltage.fcs(fcs),
               lasers = lasers.fcs(fcs))
  # extract raw data
  raw <- data.fcs(fcs)
  # apply basic gating process to data
  gated <- process.fcs(raw)
  # summarise data
  summary <- summary.fcs(gated)
  # compile data into a list
  data <- list(raw = raw,
               gated = gated,
               summary = summary)
  # create some basic plots
  plots <- list(
    # plot raw data (fsc v. ssc)
    density.raw = plot.density(raw),
    # plot raw data (each channel)
    wave.raw = plot.wave(raw),
    # plot gated data (fsc v. ssc)
    density.gated = plot.density(gated),
    # plot gated data (each channel)
    wave.gated = plot.wave(gated)
  )
  # compile all extracted data into a list
  extracted <- list(meta = meta,
                    data = data,
                    plots = plots)
  return(extracted)
}
#----
#' @title simple.fcs
#' @param path a string that is a path to either a single .fcs file or a folder containing multiple .fcs files. function checks for presence of the .fcs extension to decide whether this is a file or a folder.
#' @description quickly extracts only the most relevant data needed for all downstream analyses. this is a truncated version of the file information and all the summarised data. 
simple.fcs <- function(path){
  # make an empty list for storing data
  x <- list()
  y <- list()
  if(str_detect(path, "\\.fcs$") == T){
    fcs <- extract.fcs(path = path)
    # for each element of the list
    # extract the most relevant metadata
    x[[1]] <- fcs[["meta"]][["info"]] %>%
      pivot_wider(names_from = "info",
                  values_from = "description") %>%
      select(experiment, specimen, tube, events, time)
    # extract the summarised data
    y[[1]] <- fcs[["data"]][["summary"]]
    # apply the file names (for joining)
    names(x)[[1]] <- path
    names(y)[[1]] <- path
  }
  else {
    fcs <- read.it(path = path, pattern = ".fcs")
    # for each element of the list
    for(i in 1:length(fcs)){
      # extract the most relevant metadata
      x[[i]] <- fcs[[i]][["meta"]][["info"]] %>%
        pivot_wider(names_from = "info",
                    values_from = "description") %>%
        select(experiment, specimen, tube, events, time)
      # extract the summarised data
      y[[i]] <- fcs[[i]][["data"]][["summary"]]
      # apply the file names (for joining)
      names(x)[[i]] <- names(fcs)[[i]]
      names(y)[[i]] <- names(fcs)[[i]]
    }
  }
  # unlist into a new dataframe
  x <- plyr::ldply(x, .id = "file")
  y <- plyr::ldply(y, .id = "file")
  # merge the metadata to the data
  z <- left_join(x = x, y = y,
                 by = "file",
                 multiple = "all")
  return(z)
}
#----
#' @name comp.fcs
#' @param data a dataframe of extracted, summarised, and processed flow cytometry data containing information on the experiments, plates, channels, and parameters being measured
#' @description a looping function that builds upon the autospill algorithm with a simple intercept model using the bayesian framework to "compensate" and "normalise" flow cytometry data. this function loops through your dataframe compensating each parameter of each channel and normalising each plate from each experiment. the models (at the current time) default to identifying sample number 1 (which in most cases should be the negative control) as the intercept. following each model, the posterior distributions of each sample are summarised to the mean value and returned as a dataframe.
#' @return a dataframe of posterior distributions summarised to the mean value for each tube or well. re-incorporates the experiment, plate, channel, and parameter information back into the dataframe for later joining to the original data or 
comp.fcs <- function(data){
  # make an empty tibble to hold the compensated data
  comp <- tibble()
  # make a vector of each unique parameter
  params <- unique(
    data$param
  )
  # for each parameter...
  for(h in 1:length(params)){
    # make a vector of each channel
    channels <- unique(
      data$ch[data$param == unique(
        data$param
      )[h]]
    )
    # for each channel...
    for(i in 1:length(channels)){
      # make a vector of the plates
      plates <- unique(
        data$plate[data$ch == unique(
          data$ch[data$param == unique(
            data$param
          )[h]]
        )[i]]
      )
      # for each plate...
      for(j in 1:length(plates)){
        # make a list of the experiments
        experiments <- unique(
          data$experiment[data$plate == unique(
            data$plate[data$ch == unique(
              data$ch[data$param == unique(
                data$param
              )[h]]
            )[i]]
          )[j]]
        )
        # for each experiment...
        for(k in 1:length(experiments)){
          # isolate the data you're going to model
          .data <- data[data$experiment == experiments[k] &
                          data$plate == plates[j] &
                          data$ch == channels[i] &
                          data$param == params[h],]
          # make an intercept model for compensation and normalising
          ## note: 1 - the negative control - is set as the default level
          .bf <- bf(geo.mean | se(r.sd) ~ 1 + .no) + gaussian()
          # run the formula
          .mod <- brm(formula = .bf,
                      data = .data,
                      iter = 10000, warmup = 5000,
                      chains = 4, cores = 4)
          # extract the posterior distributions
          .post <- .mod %>%
            posterior_samples() %>%
            select(grep("^b_", colnames(.))) %>%
            pivot_longer(everything(),
                         names_to = ".no") %>%
            mutate(.no = str_replace(.no, "Intercept", ".no1"),
                   .no = str_remove(.no, "^b_"),
                   .no = str_remove(.no, "^\\.no")) %>%
            group_by(.no) %>%
            reframe(comp = mean(value)) %>%
            mutate(experiment = experiments[k],
                   plate = plates[j],
                   ch = channels[i],
                   param = params[h])
          # add the compensated values to the dataframe
          comp <- rbind(comp, .post)
        }
      }
    }
  }
  return(comp)
}

#----
#' @title extract.comp
#' @param model a brms object of a flow cytometry data set that has been modelled using a compensation formula
#' @description extracts compensated channel values for each sample after running the data through a compensation model
extract.comp <- function(model){
  model %>%
    posterior_samples() %>%
    # select only the "main effects"
    select(grep("^b_", colnames(.))) %>%
    rename_with(.cols = everything(),
                .fn = ~str_remove(.x, "^b_")) %>%
    # summarise each sample by channel
    reframe(across(.cols = everything(),
                   .fns = mean)) %>%
    # pivot to long form
    pivot_longer(cols = everything(),
                 names_to = "ch:sample",
                 values_to = "value") %>%
    # apply the correct channel and sample names
    mutate(`ch:sample` = ifelse(str_detect(`ch:sample`,"sample") == F,
                                str_c(`ch:sample`, "samplescc", sep = ":"),
                                ifelse(str_detect(`ch:sample`, "^ch") == F,
                                       str_c("chaf",`ch:sample`, sep = ":"),
                                       `ch:sample`))) %>%
    # separate the channel and sample columns
    separate(col = `ch:sample`,
             into = c("ch", "sample"),
             sep = ":") %>%
    # clean up the names
    mutate(ch = str_remove(ch, "^ch"),
           sample = str_remove(sample, "^sample")) %>%
    # pivot to wide form for downstream analyses
    pivot_wider(names_from = "ch",
                values_from = "value") 
}

#----
#-------------------
# figures
#-------------------
#----
#' @name theme_custom
#' @description customised theme for making publication-ready plots. simply add to a ggplot object. credit to rory telemeco for original version of code.
# start with the classic theme
theme_custom <- theme_classic(base_size = 24) +
  # adjust axis title position
  theme(axis.title.y = element_text(vjust=1.5), 
        axis.title.x = element_text(vjust=0.2)) + 
  # adjust plot margins and line element size
  theme(plot.margin = unit(c(.3,.3,.6,.6), "cm"), 
        line = element_line(size = 1.25)) + 
  # draw x and y axes
  theme(axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_line(colour = "black")) + 
  # put margins around axis labels so that nothing overlaps
  theme(axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) + 
  # move tickmarks inside the axes and paint black
  theme(axis.ticks.length = unit(-0.3, "cm")) + 
  # spread out facets
  theme(panel.spacing = unit(2, units = "lines")) + 
  # make tick marks black
  theme(axis.ticks = element_line(color = "black")) + 
  # remove border from facet labels
  theme(strip.background = element_blank()) 
#----
#' @title plot.density
#' @description a wrapper for making a ggplot object of the density of cells on two axes. creates a density column that estimates distance of each observation from the center of the channel. the product of density on both channels defines the density of cells. 
#' @param data a dataframe of flow cytometry data
#' @param x an object; the channel to plot on the x axis; defaults to fsc (forward scatter)
#' @param y an object; the channel to plot on the y axis; defaults to ssc (side scatter)
#' @param param a character string specifying the parameter to plot. defaults to area ("a"). other options are height ("h") and width ("w").
#' @return a ggplot object
plot.density <- function(data, x = fsc, y = ssc, param = "a"){
  data %>%
    dplyr::filter(param == param) %>%
    tidyr::pivot_wider(names_from = "ch",
                       values_from = "val") %>%
    dplyr::mutate(.d.x = base::scale({{x}}, scale = F),
                  .d.y = base::scale({{y}}, scale = F),
                  .d.z = base::sqrt(base::abs(.d.x * .d.y))) %>%
    dplyr::mutate(dplyr::across(.cols = .d.x:.d.z,
                                .fns = ~base::abs(.x))) %>%
    ggplot2::ggplot(aes(x = {{x}},
                        y = {{y}},
                        color = .d.z,
                        alpha = 0.1)) +
    ggplot2::geom_point(shape = 21) +
    ggplot2::guides(color = "none", alpha = "none") +
    ggplot2::scale_color_gradient2(low = "darkred",
                                   mid = "gold",
                                   high = "darkblue",
                                   midpoint = 0.5) +
    ggplot2::coord_cartesian(xlim = c(0, 5),
                             ylim = c(0, 5)) +
    theme_custom
}
#----
#' @title plot.wave
#' @description wrapper for making a ggplot object that plots the fluorescent intensities in each channel of a dataframe as density plots (waves). 
#' @param data a dataframe of flow cytometry data
#' @param param a character string specifying the parameter to plot. defaults to area ("a"). other options are height ("h") and width ("w").
plot.wave <- function(data, param = "a"){
  data %>%
    dplyr::filter(param == param) %>%
    ggplot2::ggplot(aes(x = val, fill = ch, alpha = 0.1)) +
    ggplot2::geom_density() +
    ggplot2::facet_grid(rows = vars(ch)) +
    ggplot2::guides(alpha = "none", fill = "none") +
    ggplot2::coord_cartesian(xlim = c(0,5)) +
    theme_custom
}
#----
