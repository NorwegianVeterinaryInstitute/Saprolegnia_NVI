# Setup  ----
library(here)
library(tidyverse)
library(plyr)
library(scales)
#library(ggridges) # cool plots to do
library(zoo)
library(Cairo)
# Analyses ----

## getting list of files to merge ----

# list_preproc <- 
#   list.files(here::here("R", "data", "raw"), 
#                            pattern = "_preproc.depth_all_pos", 
#                            full.names = F)
# list_preproc
# 
# list_postproc <- 
#   list.files(here::here("R", "data", "raw"), 
#               pattern = "_postproc.depth_all_pos", 
#               full.names = F)
# list_postproc


## Functions ----- 

## for merging in long format ---- 

# file_to_table <- function(file, indir, pattern, type) {
#   
#   df <- readr::read_tsv(here::here(indir, file), col_names = c("contig", "pos", "depth")) 
#   id <- stringr::str_remove(file, pattern)
#   
#   df <- df %>%
#     dplyr::mutate(id = id, 
#            type = type)
#   
# }

## for merging in wide format ----

file_to_table2 <- function(file, indir, pattern, type) {
  
  df <- readr::read_tsv(here::here(indir, file), col_names = c("contig", "pos", "depth")) 
  id <- stringr::str_remove(file, pattern)
  
  new_name <- paste0(id, "_depth_", type) 
  df <- df %>%
    dplyr::rename({{new_name}} := depth )
  
}

# test on one 
# file_to_table(list_preproc[1], 
#               here::here("R", "data", "raw"), 
#               pattern = "_preproc.depth_all_pos", 
#               type = "preproc") %>% head()
#   
  
# now getting that for all pre and all post 
# all_preproc <- 
#   plyr::ldply(
#     list_preproc, 
#     file_to_table, 
#     indir = here::here("R", "data"), 
#     pattern = "_preproc.depth_all_pos", 
#     type = "preproc"
#     ) 
# 
# all_preproc %>% head()

# # saving in case of crash 
# saveRDS() 
# 
# 
# # post processing 
# all_postproc <- 
#   plyr::ldply(
#     list_postproc, 
#     file_to_table, 
#     indir = here::here("R", "data"), 
#     pattern = "_postproc.depth_all_pos", 
#     type = "preproc"
#   ) 

# This is too heavy trying to do plots pre and post for one sample at the time then
VI02736_pre <- file_to_table2("VI02736_preproc.depth_all_pos", 
                             here::here("R", "data", "raw"), 
                             pattern = "_preproc.depth_all_pos", 
                             type = "preproc")

VI02736_post <- file_to_table2("VI02736_postproc.depth_all_pos", 
                              here::here("R", "data", "raw"), 
                              pattern = "_postproc.depth_all_pos", 
                              type = "postproc")
VI0236 <- 
  VI02736_pre %>% 
  full_join(VI02736_post)

rm(list = c("VI02736_pre" , "VI02736_post"))
head(VI0236)

saveRDS(VI0236, here::here("R", "data", "VI0236_depths.rds"))
VI0236 <- readRDS(here::here("R", "data", "VI0236_depths.rds"))

# VI0236 %>% 
#   filter(contig %in% c("bctg00000000", "bctg00000001", "bctg00000002", "bctg00000003")) %>%
#   pull(contig) %>% unique()

test_df <- 
  VI0236 %>% 
  filter(contig %in% c("bctg00000000")) #%>% 
  filter( pos < 1500000 & pos > 1000000)

head(test_df)
# zoo mean sliding window 80 
# https://stats.stackexchange.com/questions/3051/mean-of-a-sliding-window-in-r
sliding_window <- 300
test_df <- 
  test_df %>% 
  dplyr::mutate(preprocessing = zoo::rollmean(VI02736_depth_preproc, sliding_window, fill = 0)) %>%
  dplyr::mutate(postprocessing = zoo::rollmean(VI02736_depth_postproc, sliding_window, fill = 0)) %>%
  # transform to long format
  pivot_longer(cols = c(preprocessing, postprocessing), names_to = "Stage", values_to = "depth")
  
# Example of effect of masking and filtering reads. 
mean_depth_bfafter_VI02736  <- 
  ggplot(test_df) +
  geom_line(aes(x = pos, y = depth, color = Stage), linewidth = .5) +
  labs(title = "Example: Mean depth over a sliding window of 300pb, for VI02736 - Pathogenecity 89%", 
       subtitle = "Before (black) and after (blue) filtering mapped reads.\nLow complexity regions (masked) have a depth of 0.",  
       x = "Position", 
       y = expression(sqrt(Depth))
       ) +
  scale_y_sqrt() +
  #scale_color_manual(values = c(alpha(c("blue", "black"), alpha = 0.5)))+
  scale_color_manual(values = c(alpha(c("blue", "black"))))+
  facet_wrap(~Stage, nrow = 2) +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 

mean_depth_bfafter_VI02736  
ggsave(mean_depth_bfafter_VI02736,
       here::here("R", "figures", "mean_depth_bfafter_VI02736.svg"), 
       dpi = 300, 
       width = 30, hight = 15, units = "cm")
#alpha(c("black", "blue"), alpha = 0.5)
# when in wider format 
# ggplot(test_df) +
#   geom_line(aes(x = pos, y = mean_pre_150), color = "black") +
#   geom_line(aes(x = pos, y = mean_post_150), color = "blue") +
#   labs(title = "Mean depth over a sliding window of 150, for VI02736 - Pathogenecity 89%", 
#        subtitle = "Before (black) and after (blue) filtering mapped reads",  
#        x = "Position", 
#        y = "Depth"
#        ) +
#   scale_y_sqrt() +
#   theme_minimal()

## Now I want to see different samples with different mappling same contigs


VI07096
Vi02736
VI02391
VI02761
VI02740

# Relics TESTS ---- 
#test_plot <- 
  ggplot(VI0236 %>% filter(contig %in% c("bctg00000000")) %>% 
           filter( pos < 1500000 & pos > 1000000))   +
  #ggplot(VI0236 %>% filter(contig %in% c("bctg00000000", "bctg00000001", "bctg00000002", "bctg00000003")))   +
  geom_smooth(aes(x = pos, y = VI02736_depth_preproc), color = "black", method = "loess") +
  #geom_smooth(aes(x = pos, y = VI02736_depth_postproc), color = "blue", method = "loess") +
  #geom_line(aes(x = pos, y = VI02736_depth_preproc), color = "black") +
  #geom_line(aes(x = pos, y = VI02736_depth_postproc), color = "blue") +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)),
  #               na.value=-Inf) +
  scale_y_sqrt() +
  facet_wrap(~contig, nrow = 4) +
  theme_minimal() +
  labs(title = "Depth plot for VI02736 - Pathogenecity 89%", 
    subtitle = "Before (black) and after (blue) filtering mapped reads, for four contigs",  
       x = "Position", 
       y = "sqrt(Depth)"
       ) 

test_plot

# https://r-charts.com/distribution/ggridges/
