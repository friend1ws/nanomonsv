library(tidyverse)
library(ggrepel)

format_bp <- function(x) {
  x2 <- unlist(strsplit(x, ','))
  paste(x2[1], ":", prettyNum(x2[2], big.mark=","), " (", x2[3], ")", sep = "")
}

args <- commandArgs(trailingOnly = TRUE)

nanomonsv_prefix <- args[1]
genome_file <- paste(nanomonsv_prefix, ".bwa.txt", sep = "")
repeat_file <- paste(nanomonsv_prefix, ".rmsk.txt", sep = "")
output_dir <- args[2]

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

D1 <- read_tsv(genome_file) %>% 
  mutate(Target_Align_Label = paste(Target_Align_Chromosome, Target_Align_Start, Target_Align_End, Query_Align_Strand, Mapping_Quality, sep = "," )) %>%
  select(Contig_ID, Contig_Len, Query_Align_Start, Query_Align_End, Target_Align_Label, Mapping_Quality) %>%
  arrange(Contig_ID, Query_Align_Start)

# contig_list <- D1 %>% pull(Contig_ID) %>% unique()
tcontig_list <- D1 %>% select(Contig_ID, Contig_Len) %>% unique()
contig_id2len <- tcontig_list$Contig_Len
names(contig_id2len) <- tcontig_list$Contig_ID

genome_seg <- data.frame()

temp_contig = ""
tempy = 0
is_odd = 0
for (i in 1:nrow(D1)) {
  if (D1$Contig_ID[i] != temp_contig) {
    tempy <- tempy + 1
    temp_contig <- D1$Contig_ID[i]

  }
  
  genome_seg <- rbind(genome_seg,
                      data.frame(contig = temp_contig,
                                 xstart = D1$Query_Align_Start[i],
                                 xend = D1$Query_Align_End[i],
                                 is_odd = is_odd,
                                 label = D1$Target_Align_Label[i],
                                 mapq = D1$Mapping_Quality[i]))
  if (is_odd == 0) {
    is_odd <- 1
  } else {
    is_odd <- 0
  }
  
}


D2 <- read_tsv(repeat_file) %>% 
  mutate(Repeat_Align_Label = paste(Repeat_NAME, Repeat_Align_Start, Repeat_Align_End, Query_Align_Strand, sep = "," )) %>%
  select(Contig_ID, Contig_Len, Query_Align_Start, Query_Align_End, Repeat_Align_Label, Repeat_Class) %>% 
  arrange(Contig_ID, Query_Align_Start)


repeat_seg <- data.frame()

temp_contig = ""
tempy = 0
is_odd = 0
for (i in 1:nrow(D2)) {
  if (D2$Contig_ID[i] != temp_contig) {
    tempy <- tempy + 1
    temp_contig <- D2$Contig_ID[i]
  }
  
  repeat_seg <- rbind(repeat_seg,
                      data.frame(contig = temp_contig,
                                 xstart = D2$Query_Align_Start[i],
                                 xend = D2$Query_Align_End[i],
                                 is_odd = is_odd,
                                 label = D2$Repeat_Align_Label[i],
                                 class = D2$Repeat_Class[i]))
  if (is_odd == 0) {
    is_odd <- 1
  } else {
    is_odd <- 0
  }
  
}


alpha_mapq_trans <- function() {
  scales::trans_new("alpha_mapq", 
                    trans = function(x) {0.9 * x / 60 + 0.1},
                    inverse = function(x) {60 * (x - 0.1) / 0.9})
}

repeat_class_list <- unique(repeat_seg$class)
cvalues <- rep("#fdc086", length(repeat_class_list))
names(cvalues) <- repeat_class_list
cvalues[str_detect(names(cvalues), "Satellite")] <- "#386cb0"
cvalues[str_detect(names(cvalues), "LINE")] <- "#f0027f"
cvalues[str_detect(names(cvalues), "SINE")] <- "#7fc97f"
cvalues[str_detect(names(cvalues), "Simple_repeat")] <- "#bf5b17"


for (i in 1:length(contig_id2len)) {
  tcontig <- names(contig_id2len)[i]
  tcontig_len <- as.numeric(contig_id2len[i])
  
  tcontig_seg <- data.frame(xstart = 0, xend = tcontig_len)
  
  tgenome_seg <- genome_seg %>% filter(contig == tcontig)
  trepeat_seg <- repeat_seg %>% filter(contig == tcontig)
  
  # dummy
  tgenome_seg <- rbind(tgenome_seg,
                       data.frame(contig = tcontig, xstart = 0, xend = 0, is_odd = 0, label = "", mapq = 0),
                       data.frame(contig = tcontig, xstart = 0, xend = 0, is_odd = 0, label = "", mapq = 60)
                       )
  
  ggplot() + 
    geom_segment(data = tcontig_seg, aes(x = xstart, xend = xend, y = 3, yend = 3), arrow = arrow(length = unit(0.2, "inches"))) +
    geom_segment(data = tgenome_seg, aes(x = xstart, xend = xend, y = 4 + is_odd, yend = 4 + is_odd, alpha = mapq), colour = "#666666", size = 3) +
    geom_text_repel(data = tgenome_seg, aes(x = 0.5 * (xstart + xend), y = 4 + is_odd, label = label), size = 2) +
    geom_segment(data = trepeat_seg, aes(x = xstart, xend = xend, y = 2 - is_odd, yend = 2 - is_odd, colour = class), size = 3) +
    geom_text_repel(data = trepeat_seg, aes(x = 0.5 * (xstart + xend), y = 2 - is_odd, label = label), size = 2) +
    geom_point(data = data.frame(x = c(0), y = c(0)), aes(x = x, y = y), alpha = 0) + 
    geom_point(data = data.frame(x = c(0), y = c(6)), aes(x = x, y = y), alpha = 0) +
    ggtitle(format_bp(tcontig)) +
    labs(x = "Contig coordinate") +
    scale_alpha_continuous(range = c(0.2, 1) ) +
    guides(alpha = FALSE, colour = FALSE) +
    scale_colour_manual(values = cvalues) +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
    ggsave(paste(output_dir, "/", tcontig, ".pdf", sep = ""), width = 20, height = 10, units = "cm")
  
}

warnings()

