#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readxl)
library(writexl)
library(broom)

remotes::install_cran("plyranges", upgrade = "never")
library(plyranges)

remotes::install_cran("MAnorm2", upgrade = "never")
library(MAnorm2)

remotes::install_cran("gtools", upgrade = "never")
library(gtools)

remotes::install_cran("openxlsx2", upgrade = "never")
library(openxlsx2)

library(ggpubr)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(cowplot)


OUTPUTDIR <- "result"
NEED_OUTPUT_BOOL <- T

## create output directory
create_directory <- function(dirname){
  if (dir.exists(dirname)){
    unlink(dirname, recursive=TRUE)  
  }
  dir.create(dirname, showWarnings = TRUE, recursive = FALSE, mode = "0777")
  return(dirname)
}


create_directory(OUTPUTDIR)

get_fisher <- function(x1,y1,x2,y2){
  mx <- matrix(c(x1,y1,x2,y2), nr = 2)
  fisher.test(mx) %>% tidy() %>% select(ES=estimate, p.value)
}

bio_group_dict <- list.files(path = "./bio_regions/", pattern = "bed", full.names = T, recursive = T) %>% 
  map_dfr(\(f){
    tibble(bio_group = f %>% dirname() %>% basename(),
           filename = f %>% basename())
  })

fg_list <- list.files(path = "./diffexp", pattern = "bed", full.names = T) %>% 
  map(\(f){
    filename <- basename(f)
    read_tsv(file = f, col_names = F, show_col_types = FALSE) %>% 
      select(seqnames = 1, start = 2, end = 3) %>% 
      mutate(fgbg_name = filename,
             total = n(),
             type = "foreground")
  })


bg_list <- list.files(path = "./background/", pattern = "bed", full.names = T) %>% 
  map(\(f){
    filename <- basename(f)
    read_tsv(file = f, col_names = F, show_col_types = FALSE) %>% 
      select(seqnames = 1, start = 2, end = 3) %>% 
      mutate(fgbg_name = filename,
             total = n(),
             type = "background")
  })

calc_overlap_v3 <- function(fg,bg,bioreg.df,prepare_output_overlaps) {

  fname <- unique(bioreg.df$bioreg_name)
  sorted_bioreg <- bioreg.df %>% arrange(seqnames, start, end) %>% as_granges()
  
  fg_overlap <- find_overlaps(fg %>% as_granges() %>% arrange(seqnames, start, end), 
                              sorted_bioreg, 
                              minoverlap = 2) %>% 
    as_tibble() %>% 
    distinct()
  
  bg_overlap <- find_overlaps(bg %>% arrange(seqnames, start, end) %>% as_granges(), 
                              sorted_bioreg, 
                              minoverlap = 2) %>% 
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
    
    print(category)
    
    ## read bio-regions
    bioreg_list <- list.files(path = str_glue("./bio_regions/{category}/"), pattern = "bed", full.names = T) %>% 
      mixedsort(.) %>% 
      map(\(f){
        filename <- basename(f)
        read_tsv(file = f, col_names = F, show_col_types = FALSE) %>% 
          select(seqnames = 1, start = 2, end = 3) %>% 
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
      map_dfr(tmp, \(t){t$dump}) %>% 
        # write_csv(file = str_glue("{OUTPUTDIR}/overlaps_{category}.csv.gz"), col_names = T)
        openxlsx2::write_xlsx(file = str_glue("{OUTPUTDIR}/overlaps_{category}.xlsx"), col_names = TRUE)
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

make_barplot <- function(df){

  order_df <- unique(df$bioreg.fig) %>% mixedsort() %>% tibble(bioreg.fig = .) %>% mutate(ix = row_number())
  
  category <- unique(df$category)
  
  bg_name <- df$bg.name.fig %>% unique()

  test_df <- df %>% left_join(., order_df) %>% 
    mutate(
      p.value = case_when(is.infinite(ES) ~ 0,
                          TRUE ~ p.value),
      ES = case_when(is.infinite(ES) ~ 0,
                     TRUE ~ ES),
      star = case_when(p.value <=1e-100 & ES > 0 ~ "***",
                       p.value <= 1e-10 & ES > 0 ~ "**",
                       p.value <= 0.001 & ES > 0 ~ "*",
                       TRUE ~ ""))
  
  
  legend_ylim <- max(test_df$ES) + 5
  

  title <- str_glue("delta-DHS Enrichment for \n",
                    "{category}\n",
                    "{bg_name} background\n",
                    "p.value <= 1e-100 = ***; 1e-10 = **; 0.001 = *")
  
  barplot <- ggplot(data=test_df, aes(x=factor(ix), y=ES, fill=fg.name.fig)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_hline(yintercept=1, colour="Red") +
    scale_y_continuous(limits=c(0, legend_ylim), breaks=seq(0, legend_ylim, 1.0)) +
    scale_fill_discrete(name = "foreground")+
    ggtitle(title) +
    geom_text(aes(label=star,y=ES),position = position_dodge(width=0.9),size=7,hjust=-.3, vjust=0.75, angle=90)+
    ylab("Enrichment Score") + xlab("Biological Sites") + 
    theme(legend.position = c(1, 1), 
          legend.justification = c(1, 1),
          plot.title = element_text(size = 12))

  names_table <- ggtexttable(test_df %>% select(index = ix, Bio.region = bioreg.fig) %>% distinct(), rows = NULL)
  
  cowplot::plot_grid(barplot+names_table+plot_layout(ncol = 2))
}

es_score_list %>% 
  map(make_barplot) %>%
  marrangeGrob(nrow = 1, ncol = 1) %>%
  ggsave(filename = str_glue("{OUTPUTDIR}/enrichment_barplots.pdf"), width = 15, height = 11.29)
  

### OUTPUT FINAL TABLE
options(scipen=0)
zz1 <- es_scores %>% 
  mutate(fg.pct.overlap = round(100.*fg.overlap/fg.total,2),
         bg.pct.overlap = round(100.*bg.overlap/bg.total,2),
        #fg.pct.overlap = str_glue("{fg.overlap} ({fg.pct.overlap}%)"),
        #bg.pct.overlap = str_glue("{bg.overlap} ({bg.pct.overlap}%)"),
        #bg.name = str_glue("{bg.name} background"),
        #fg.name = str_glue("{fg.name} foreground"),
        #p.value = format(p.value, digits = 3, scientific = T, trim = T),
        p.value = formattable::comma(p.value, format = "e", width = 2, digits = 2),
        ES = ifelse(is.infinite(ES), 0, round(ES,2))) %>%   
  select(Bio.region_category = category, 
         Bio.region = bioreg, 
         `Enrichment_Score (ES)` = ES, 
         Fishers_Exact_test_pvalue = p.value, 
         Foreground_filename = fg.name, 
         Background_filename = bg.name, 
         Foreground_total_sites = fg.total, 
         Background_total_sites = bg.total, 
         Foreground_overlap_sites = fg.overlap,
         Background_overlap_sites = bg.overlap,
         Foreground_overlap_sites_percent = fg.pct.overlap, 
         Background_overlap_sites_percent = bg.pct.overlap) 

openxlsx2::write_xlsx(zz1, file = str_glue("{OUTPUTDIR}/overlap_summary.xlsx"), col_names = TRUE)


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





