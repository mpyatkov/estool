#!/usr/bin/env Rscript

## suppress warnings
options(warn=-1)

library <- function (...) {
   packages <- as.character(match.call(expand.dots = FALSE)[[2]])
   suppressWarnings(suppressMessages(lapply(packages, base::library, character.only = TRUE)))
   return(invisible())
}

## uncomment this if you are not going to use singularity image
## it should install the packages locally

# cran_packages <- c("writexl", "gtools", "openxlsx2", "formattable")
# for (pkg in cran_packages) {
#   if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
#     install.packages(pkg)
#     library(pkg)
#   }
# }

# if (!require("plyranges", character.only = TRUE, quietly = TRUE)) {
#     if (!require("BiocManager", quietly = TRUE))

#         install.packages("BiocManager")

#     BiocManager::install("plyranges", update = FALSE)
#     library(plyranges)
#   }

cat("Init libraries...\n")
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readxl)
library(broom)
library(ggpubr)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(cowplot)
library(gtools) ## for mixedsort
library(plyranges)
library(formattable)
library(openxlsx2)
library(RColorBrewer)

OUTPUTDIR <- "result"
NEED_OUTPUT_BOOL <- T

## create output directory
create_directory <- function(dirname){
  if (dir.exists(dirname)){
    print(str_glue("Cleanup previous results..."))
    unlink(dirname, recursive=TRUE)  
  }
  dir.create(dirname, showWarnings = TRUE, recursive = FALSE, mode = "0777")
  return(dirname)
}

outd <- create_directory(OUTPUTDIR)
print(str_glue("The output '{outd}' directory has been created"))

get_fisher <- function(x1,y1,x2,y2){
  mx <- matrix(c(x1,y1,x2,y2), nr = 2)
  fisher.test(mx) %>% tidy() %>% select(ES=estimate, p.value)
}

bio_group_dict <- list.files(path = "./bio_regions/", pattern = "\\.bed$", full.names = T, recursive = T) %>% 
  map_dfr(\(f){
    tibble(bio_group = f %>% dirname() %>% basename(),
           filename = f %>% basename())
  })

fg_list <- list.files(path = "./foreground", pattern = "bed", full.names = T) %>% 
  map(\(f){
    filename <- basename(f)
    read_tsv(file = f, col_names = F, show_col_types = FALSE) %>% 
      select(seqnames = 1, start = 2, end = 3) %>% 
      distinct(seqnames, start, end) %>% 
      filter(end - start > 1) %>% 
      mutate(fgbg_name = filename,
             total = n(),
             type = "foreground")
  })


bg_list <- list.files(path = "./background/", pattern = "bed", full.names = T) %>% 
  map(\(f){
    filename <- basename(f)
    read_tsv(file = f, col_names = F, show_col_types = FALSE) %>% 
      select(seqnames = 1, start = 2, end = 3) %>% 
      distinct(seqnames, start, end) %>% 
      filter(end - start > 1) %>% 
      mutate(fgbg_name = filename,
             total = n(),
             type = "background")
  })

calc_overlap_v3 <- function(fg,bg,bioreg.df,prepare_output_overlaps) {

  fname <- unique(bioreg.df$bioreg_name)
  sorted_bioreg <- bioreg.df %>% arrange(seqnames, start, end) %>% as_granges()
  
  fg_overlap <- find_overlaps(fg %>% as_granges() %>% arrange(seqnames, start, end), 
                              sorted_bioreg, 
                              minoverlap = 1) %>% 
    as_tibble() %>% 
    distinct()
  
  bg_overlap <- find_overlaps(bg %>% arrange(seqnames, start, end) %>% as_granges(), 
                              sorted_bioreg, 
                              minoverlap = 1) %>% 
    as_tibble() %>% 
    distinct()
  
  summary <- tibble(
    bioreg = fname,
    bioreg.total = nrow(bioreg.df),
    fg.name = unique(fg$fgbg_name),
    bg.name = unique(bg$fgbg_name),
    fg.total = nrow(fg),
    bg.total = nrow(bg),
    fg.overlap = nrow(fg_overlap),
    bg.overlap = nrow(bg_overlap),
    fg.nonoverlap = fg.total - fg.overlap,
    bg.nonoverlap = bg.total - bg.overlap)

  if (prepare_output_overlaps)
    list(dump = bind_rows(fg_overlap, bg_overlap) %>% 
           select(-width,-strand,-total) %>% 
           dplyr::rename("input_type" = "type", "input_name" = "fgbg_name"), 
         summary = summary)
  else{
    list(summary = summary)  
  }
}

#categories <- "01_Male_Chromatin_States"
categories <- sort(bio_group_dict$bio_group %>% unique())

#system.time(zz1 <- map_dfr(sort(bio_group_dict$bio_group %>% unique()), \(d){
system.time(
  
  ## loop over all bioreg categories
  es_scores <- map_dfr(categories, \(category){
    
    print(str_glue("Start processing '{category}' bio regions"))
    
    ## read bio-regions
    bioreg_list <- list.files(path = str_glue("./bio_regions/{category}/"), pattern = "bed", full.names = T) %>% 
      mixedsort(.) %>% 
      map(\(f){
        filename <- basename(f)
        read_tsv(file = f, col_names = F, show_col_types = FALSE) %>% 
          select(seqnames = 1, start = 2, end = 3) %>% 
          distinct(seqnames, start, end) %>% 
          filter(end - start > 1) %>% 
          mutate(bioreg_name = filename)
      })
    
    
    ## create all possible pairs of foregrounds, backgrounds and bioregs 
    ## NEED_OUTPUT_BOOL = T adds dump with overlaps to output list
    input <- tidyr::expand_grid(fg = fg_list,
                                bg = bg_list,
                                bioreg = bioreg_list) %>% 
      mutate(prepare_output_overlaps = NEED_OUTPUT_BOOL)
    
    
    tmp <- pmap(input, calc_overlap_v3) ## list of lists(summary, dump)
    
    ## SIDE EFFECT: write down all overlaps to individual files for each category
    if (NEED_OUTPUT_BOOL){
      ## aggregate all dump tables from one category to one
      agg_dump <- map_dfr(tmp, \(t){t$dump}) %>% distinct_all()
      if (nrow(agg_dump) > 990000) {
        print(str_glue("Too many rows for category table '{category}'. Skipping creating overlap files for this category"))
      } else {
        # write_csv(file = str_glue("{OUTPUTDIR}/overlaps_{category}.csv.gz"), col_names = T)
        openxlsx2::write_xlsx(agg_dump, file = str_glue("{OUTPUTDIR}/overlaps_{category}.xlsx"), col_names = TRUE)
        
      }
    }
    
    ## for each row of data.frame in summary calculate fisher test
    ## aggregate to data.frame
    map_dfr(tmp, \(l)l$summary) %>% mutate(category = category) %>% 
      rowwise() %>% 
      mutate(ft = list(get_fisher(fg.overlap, fg.nonoverlap, bg.overlap, bg.nonoverlap))) %>% 
      unnest(ft) %>% 
      ungroup() %>% 
      mutate(fg.name.fig = str_glue("{fg.name} ({fg.total})"),
             bioreg.fig = str_glue("{bioreg} ({bioreg.total})"),
             bg.name.fig = str_glue("{bg.name} ({bg.total})"))
  }))



## loop over 
es_score_list <- es_scores %>% 
  group_by(category, bg.name) %>% 
  group_split()

## create vector of non consecutive numbers to make neghbour colors on
## plot different (works only for then 3 colors)
# generate_non_consecutive_vector <- function(n, min_val, max_val) {
#   if (max_val - min_val + 1 < n) {
#     stop("Range is too small for the given n without consecutive numbers.")
#   }
#   repeat {
#     vec <- sample(min_val:max_val, n)
#     if (all(abs(diff(vec)) != 1)) return(vec)
#   }
# }

make_barplot <- function(df){

  order_df <- unique(df$bioreg.fig) %>% mixedsort() %>% tibble(bioreg.fig = .) %>% mutate(ix = row_number())
  
  category <- unique(df$category)
  
  bg_name <- df$bg.name.fig %>% unique()

  ## Max ES
  maxes <- max(df$ES[is.finite(df$ES)])
  
  test_df <- df %>%
    left_join(., order_df, join_by(bioreg.fig)) %>%
    mutate(star = case_when(p.value <=1e-100 ~ "***",
                            p.value <= 1e-10 ~ "**",
                            p.value <= 0.001 ~ "*",
                            TRUE ~ "")) %>% 
    mutate(fill_color= ifelse(is.infinite(ES),"Enrichment score = Inf", as.character(fg.name.fig))) %>% 
    mutate(ES = case_when(is.infinite(ES) & p.value<=0.001 ~ maxes,
                          is.infinite(ES) ~ 0,
                          .default = ES),
           fg.pct.overlap = round(100.*fg.overlap/fg.total,2))
  
  ## colors
  num_colors <- length(unique(test_df$fg.name.fig))
  #brewer_colors <- brewer.pal(n = max(3, num_colors), name = "Set2")[1:num_colors]
  
  ## get random colors order if we have a lot of colors
  # brewer_colors <- if (num_colors > 6) {
  #   colors_order <- generate_non_consecutive_vector(num_colors, 1, num_colors)
  #   colorRampPalette(brewer.pal(6,"Set2"))(num_colors)[colors_order]
  # } else {
  #   brewer.pal(n = max(3, num_colors), name = "Set2")[1:num_colors]
  # }
  
  brewer_colors <- c("#66C2A5", "#FFC94A", "#8DA0CB", "#F55FB9", "#A6D854", 
                     "#FFFF00", "firebrick2", "cadetblue2", "darkorchid1", "darkgoldenrod1",
                     "dodgerblue2", "#9BAF8D", "darkorange", "blue", "lightgoldenrod", 
                     "#66C2A5", "#FFC94A", "#8DA0CB", "#F55FB9", "#A6D854", 
                     "#FFFF00", "firebrick2", "cadetblue2", "darkorchid1", "darkgoldenrod1",
                     "dodgerblue2", "#9BAF8D", "darkorange", "blue", "lightgoldenrod")[1:num_colors]
  
  
  custom_colors <- c("gray", brewer_colors)
  names(custom_colors) <- c("Enrichment score = Inf", unique(test_df$fg.name.fig))
  
  # max ylim for ES plot
  legend_ylim <- maxes+maxes*0.15
 
  title <- str_glue("Foreground set enrichments for BioRegions \n",
                    "BioRegion set: {category}\n",
                    "{bg_name} background\n",
                    "p.value <= 1e-100 = ***; 1e-10 = **; 0.001 = *")
  
  es.barplot <- ggplot(data=test_df, aes(x=factor(ix), y=ES, fill=fill_color)) +
    geom_bar(stat="identity", position=position_dodge2(), colour="black") +
    geom_hline(yintercept=1, colour="Red") +
    scale_y_continuous(breaks = scales::pretty_breaks(20), limits = c(0, legend_ylim))+
    scale_fill_manual(
      name = "Foreground",
      values = custom_colors) +
    ggtitle(title) +
    geom_text(aes(label=star,y=ES),position = position_dodge2(width=0.9),size=7,hjust=-.3, vjust=0.75, angle=90)+
    ylab("Enrichment Score") + xlab("Biological Sites") +
    theme(legend.position = c(1.9, 1.05),
          legend.justification = c(1, 1),
          plot.title = element_text(size = 12))
  
  ## same barplot as before but y-axis is fg.pct.overlap not ES
  legend_ylim_pct <- max(test_df$fg.pct.overlap[is.finite(test_df$fg.pct.overlap)])
  legend_ylim_pct <- legend_ylim_pct + legend_ylim_pct*0.15
  
  pct.overlap.barplot <- ggplot(data=test_df, aes(x=factor(ix), y=fg.pct.overlap, fill=fill_color)) +
    geom_bar(stat="identity", position=position_dodge2(), colour="black") +
    geom_hline(yintercept=1, colour="Red") +
    scale_y_continuous(breaks = scales::pretty_breaks(20), limits = c(0, legend_ylim_pct))+
    scale_fill_manual(
      name = "Foreground",
      values = custom_colors) +
    ggtitle(title) +
    geom_text(aes(label=star,y=fg.pct.overlap),position = position_dodge2(width=0.9),size=7,hjust=-.3, vjust=0.75, angle=90)+
    ylab("Foreground Overlap as Percentage of Foreground Total") + xlab("Biological Sites") +
    theme(legend.position = c(1.9, 1.05),
          legend.justification = c(1, 1),
          plot.title = element_text(size = 12))
  
  names_table <- ggtexttable(test_df %>% select(index = ix, Bio.region = bioreg.fig) %>% distinct(), rows = NULL)
  
  list(cowplot::plot_grid(es.barplot+names_table+plot_layout(ncol = 2)), 
       cowplot::plot_grid(pct.overlap.barplot+names_table+plot_layout(ncol = 2)))
}

es_score_list %>% 
  map(make_barplot) %>%
  flatten() %>% 
  marrangeGrob(nrow = 1, ncol = 1) %>%
  ggsave(filename = str_glue("{OUTPUTDIR}/enrichment_barplots.pdf"), width = 15, height = 11.29)
  
### OUTPUT FINAL TABLE
options(scipen=0)
zz1 <- es_scores %>% 
  mutate(fg.pct.overlap = round(100.*fg.overlap/fg.total,2),
         bg.pct.overlap = round(100.*bg.overlap/bg.total,2),
         fg_to_bio.pct.overlap = round(100.*fg.overlap/bioreg.total,2),
         bg_to_bio.pct.overlap = round(100.*bg.overlap/bioreg.total,2),
         
        #p.value = format(p.value, digits = 3, scientific = T, trim = T),
        p.value = formattable::comma(p.value, format = "e", width = 2, digits = 2), 
        #ES = ifelse(is.infinite(ES), 1e300, round(ES,2))
        ES = case_when(is.infinite(ES) & p.value<=0.001 ~ 1e300,
                       is.infinite(ES) ~ 0,
                       .default = round(ES, 2))
        ) %>%   
  select(Bio.region_category = category, 
         Bio.region = bioreg, 
         Bio.region_total_sites = bioreg.total,
         `Enrichment_Score (ES)` = ES, 
         Fishers_Exact_test_pvalue = p.value, 
         Foreground_filename = fg.name, 
         Background_filename = bg.name, 
         Foreground_total_sites = fg.total, 
         Background_total_sites = bg.total, 
         Foreground_overlap_sites = fg.overlap,
         Background_overlap_sites = bg.overlap,
         Foreground_overlap_sites_percent = fg.pct.overlap, 
         Background_overlap_sites_percent = bg.pct.overlap,
         Foreground_overlap_to_Bio.region_total_percent = fg_to_bio.pct.overlap,
         Background_overlap_to_Bio.region_total_percent = bg_to_bio.pct.overlap) 

openxlsx2::write_xlsx(zz1, file = str_glue("{OUTPUTDIR}/overlap_summary.xlsx"), col_names = TRUE)
print("Please check the 'result' directory")

## TODO: stoped here because thought how to aggregate data if we have multiple backgrounds and foregrounds
#pivot_wider(names_from = c("fg.name"), values_from = c("ES","p.value","fg.pct.overlap")) %>% 
# zz2 <- zz1 %>% 
#   group_by(fg.name, bg.name) %>% 
#   group_split() %>% 
#   map_dfr(\(df){
#     df %>% pivot_wider(names_from = c("fg.name"), values_from = c("fg.pct.overlap")) %>% 
#       pivot_wider(names_from = c("bg.name"), values_from = c("bg.pct.overlap"))
#   })
  
# zz3 <- zz2[[1]]  %>% pivot_wider(names_from = c("fg.name"), values_from = c("fg.pct.overlap")) %>% 
#   pivot_wider(names_from = c("bg.name"), values_from = c("bg.pct.overlap"))
# zz3.1 <- zz2[[2]]  %>% pivot_wider(names_from = c("fg.name"), values_from = c("fg.pct.overlap")) %>% 
#   pivot_wider(names_from = c("bg.name"), values_from = c("bg.pct.overlap"))
# zz4 <- left_join(zz3, zz3.1, by = c("category", "bioreg"))

