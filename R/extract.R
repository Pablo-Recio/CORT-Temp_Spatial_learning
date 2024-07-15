# extraction
# extracting data

# clear workspace
rm(list = ls())

# load libraries
library(tidyverse)
library(here)
library(brms)
library(rstan)

# set working directory to top level
setwd(here())

# source custom functions
source("./r/func.R")

# make a path to your data
path <- "./data/"

# make an empty list to hold the data
data <- list()

# make paths to each file
paths <- list.files(path = path,
                    pattern = "\\.fcs",
                    full.names = T,
                    recursive = T)

# for each file, read in a simple fcs
for(i in 1:length(paths)){
  data[[i]] <- paths[i] %>%
    read.fcs() %>%
    data.fcs() %>%
    process.fcs() %>%
    summary.fcs()
  names(data)[[i]] <- paths[i]
}

# prepare data for compensation
df <- plyr::ldply(data) %>%
  # tidy up .id names
  mutate(.id = str_remove(.id, path),
         .id = str_remove(.id, "/"),
         .id = str_remove(.id, "\\.fcs$")) %>%
  separate(col = .id,
           into = c("experiment", "plate", "sample"),
           sep = "/") %>%
  separate(col = sample,
           into = c("specimen", "spec.no", "tube", "tube.no", ".no"),
           sep = "_") %>%
  # tidy up the acquisition number
  mutate(.no = as.numeric(.no),
         .no = factor(.no))

# r crashes if you try to do too much...
# i'll break this into chunks for now
# starting with only the area values from the medial cortex experiment
df.mc <- df %>%
  dplyr::filter(experiment == "medial cortex" &
                  param == "a")
comp.mc <- comp.fcs(df.mc)

# now olfactory bulbs
df.ob <- df %>%
  dplyr::filter(experiment == "olfactory bulb" &
                  param == "a")
comp.ob <- comp.fcs(df.ob)

# now optical tecta
df.ot <- df %>%
  dplyr::filter(experiment == "optical tecta" &
                  param == "a")
comp.ot <- comp.fcs(df.ot)

# bind and tidy everything up
comp <- bind_rows(comp.mc, comp.ob, comp.ot) %>%
  right_join(df) %>%
  unite(specimen, c(specimen, spec.no), sep = "_") %>%
  select(experiment, plate, specimen, 
         tube, tube.no, .no,
         n:r.se, comp)

# write the extracted, compensated data to a csv
write_csv(comp, "./data/data.csv")
