# Setup ---- 
library(here)
library(tidyverse)
library(plyr)
library(scales)
library(zoo)
library(Cairo)
library(ggforce)

# doing that on post processing files 
# Functions ---- 
## for merging in wide format ----

file_to_table2 <- function(file, indir, pattern, type) {
  
  df <- readr::read_tsv(here::here(indir, file), col_names = c("contig", "pos", "depth")) 
  id <- stringr::str_remove(file, pattern)
  
  new_name <- paste0(id, "_depth_", type) 
  df <- df %>%
    dplyr::rename({{new_name}} := depth )
  
}


VI07096 <- file_to_table2("VI07096_postproc.depth_all_pos", 
                               here::here("R", "data", "raw"), 
                               pattern = "_postproc.depth_all_pos", 
                               type = "postproc")
VI02736 <- file_to_table2("VI02736_postproc.depth_all_pos", 
                               here::here("R", "data", "raw"), 
                               pattern = "_postproc.depth_all_pos", 
                               type = "postproc")

VI02391 <- file_to_table2("VI02391_postproc.depth_all_pos", 
                               here::here("R", "data", "raw"), 
                               pattern = "_postproc.depth_all_pos", 
                               type = "postproc")
VI02761 <- file_to_table2("VI02761_postproc.depth_all_pos", 
                               here::here("R", "data", "raw"), 
                               pattern = "_postproc.depth_all_pos", 
                               type = "postproc")
VI02740 <- file_to_table2("VI02740_postproc.depth_all_pos", 
                               here::here("R", "data", "raw"), 
                               pattern = "_postproc.depth_all_pos", 
                               type = "postproc")
VI02395 <- file_to_table2("VI02395_postproc.depth_all_pos", 
                               here::here("R", "data", "raw"), 
                               pattern = "_postproc.depth_all_pos", 
                               type = "postproc")
test_df <- VI07096 %>%
  full_join(VI02736) %>%
  full_join(VI02391) %>%
  full_join(VI02761) %>%
  full_join(VI02740) %>%
  full_join(VI02395)


rm(list = c("VI07096", "VI02736", "VI02391", "VI02761", "VI02740", "VI02395"))


saveRDS(test_df, here::here("R", "data", "5_post_depths.rds"))
rm(test_df)
test_df <- readRDS(here::here("R", "data", "5_post_depths.rds"))


sliding_window <- 300

test_analysis <- 
  test_df %>% 
  dplyr::filter(contig %in% c("bctg00000002")) %>%
  #dplyr::filter( pos < 5000 & pos > 1000) %>%
  # mean depth sliding Window over all the depth
  dplyr::mutate(across(ends_with("depth_postproc"), ~ zoo::rollmean(., sliding_window, fill = 0))) %>%
  # Z normalisation to be able to compare the depth 
  dplyr::mutate(across(ends_with("depth_postproc"), ~ scales::rescale(.))) %>%
  # transform to long format
  pivot_longer(cols = ends_with("postproc"), names_to = "Sample", values_to = "depth") %>%
  mutate(Sample = stringr::str_remove(Sample, "_depth_postproc")) %>%
  mutate(Pathogenicity = case_when(
    Sample == "VI07096" ~ "high",
    Sample == "VI02736" ~ "89%",
    Sample == "VI02391" ~ "31%",
    Sample == "VI02761" ~ "18%",
    Sample == "VI02740" ~ "0%",
    Sample == "VI02395" ~ "0%")) %>%
  mutate(Pathogenicity = factor(Pathogenicity, levels = c("high", "89%", "31%", "18%", "0%"), ordered = T)) %>%
  arrange(Pathogenicity, pos) 
  
    
head(test_analysis)
# Damned needs to be normalized 

globalplot <- 
  ggplot(test_analysis) +
  geom_line(aes(x = pos, y = depth, color = paste(Sample, Pathogenicity, sep = " ")), alpha = .5, linewidth = .5) +
  labs(title = "Depth for samples of different pathogenecity", 
       x = "Position", 
       y = " Z normalized Depth",
       color = "Sample & Pathogenicity") +
  facet_wrap(~Pathogenicity, nrow = 6, scales = "fixed") +
  #facet_grid(~Pathogenicity, scales = "free") +
  theme_minimal() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 

ggsave(globalplot,
       here::here("R", "figures", "globalplot.png"), 
       width = 10,
       height = 20)


# I think I need a plot split by coordinates ----
test_analysis2 <- 
  test_df %>% 
  dplyr::filter(contig %in% c("bctg00000002")) %>%
  #dplyr::filter( pos < 5000 & pos > 1000) %>%
  # mean depth sliding Window over all the depth
  dplyr::mutate(across(ends_with("depth_postproc"), ~ zoo::rollmean(., sliding_window, fill = 0))) %>%
  # Z normalisation to be able to compare the depth 
  dplyr::mutate(across(ends_with("depth_postproc"), ~ scales::rescale(.))) 
  
#https://forum.posit.co/t/ggplot-best-way-to-split-x-axis-in-intervals/64931
length(test_analysis2$pos) # position stats at 1
length(test_analysis2$pos) / 10000 # 113 groups 
group_vector <- rep(x = 1:188, each = 10000)[1:length(test_analysis2$pos)] # 188 groups of 10000
length(group_vector) == length(test_analysis2$pos) 
# we need to do that when still position sorted 
test_analysis2$group_vector <- group_vector


test_analysis2 <- 
  test_analysis2 %>%
  # transform to long format
  pivot_longer(cols = ends_with("postproc"), names_to = "Sample", values_to = "depth") %>%
  mutate(Sample = stringr::str_remove(Sample, "_depth_postproc")) %>%
  mutate(Pathogenicity = case_when(
    Sample == "VI07096" ~ "high",
    Sample == "VI02736" ~ "89%",
    Sample == "VI02391" ~ "31%",
    Sample == "VI02761" ~ "18%",
    Sample == "VI02740" ~ "0%",
    Sample == "VI02395" ~ "0%")) %>%
  mutate(Pathogenicity = factor(Pathogenicity, levels = c("high", "89%", "31%", "18%", "0%"), ordered = T)) %>%
  arrange(Pathogenicity, pos) 

head(test_analysis2)

# split the data set by group_vector
list_test_analysis2 <- split(test_analysis2, test_analysis2$group_vector)
str(list_test_analysis2)
names(list_test_analysis2[1])

facet_wrap(vars(y_i), ncol = 1, scales = "free_x") 

plot_draw_save <- function(list_element_data){
  
  save_name <- names(list_element_data)
  
  p <- ggplot(list_element_data[[1]]) +
    geom_line(aes(x = pos, y = depth, color = paste(Sample, Pathogenicity, sep = " ")), alpha = .5, linewidth = .5) +
    labs(title = "Normalized mapping depth to compare samples of different pathogenecity", 
         x = "Position", 
         y = " Z normalized Depth",
         color = "Sample & Pathogenicity") +
    facet_wrap(~Pathogenicity, nrow = 6, scales = "fixed") +
    #facet_wrap_paginate(c("Pathogenicity", "group_vector"), nrow = 6, scales = "fixed", page = 1) +
    theme_minimal() +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          legend.position = "bottom") 
  
  svg(sprintf("%s/%s_pathog_depth.svg", here::here("R", "figures"), save_name))
  print(p)
  dev.off()
  
}
  

#plot_draw_save(list_test_analysis2[1])

for (i in 1:length(list_test_analysis2)){
  plot_draw_save(list_test_analysis2[i])
}
