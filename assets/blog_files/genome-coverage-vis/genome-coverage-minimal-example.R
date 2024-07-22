# Required packages
library(tidyverse)
library(ggcoverage)

# Packages for added customization
library(scales)
library(gggenes)
library(RColorBrewer)

# this is a folder with .bam files (different samples mapped to same reference)
track.folder = 'path/to/bams/' 

# Create a metadata file from .bam names
# sample.meta = data.frame(SampleName = gsub('.bam', '' , grep('bam$', list.files(track.folder), value = T)))
# sample.meta is just a 1 column df with names of bam files at this points (without the extension) equivalent to:

sample.meta <- data.frame(
  stringsAsFactors = FALSE,
  SampleName = c("CV0857_viridis_North_M.dedup.filtered.unique", "CV0985_concolor_Other_F.dedup.filtered.unique",
                 "CV1081_viridis_Mid_M.dedup.filtered.unique", "CV1082_viridis_South_M.dedup.filtered.unique", 
                 "CV1086_viridis_South_M.dedup.filtered.unique", "CV1087_viridis_North_F.dedup.filtered.unique",
                 "CV1089_viridis_South_M.dedup.filtered.unique", "CV1090_cerberus_Other_M.dedup.filtered.unique", 
                 "CV1095_viridis_North_M.dedup.filtered.unique", "CV1096_viridis_North_F.dedup.filtered.unique")
)

# Add more attributes to the metadata
sample.meta <- sample.meta %>% 
  mutate(SampleID = str_extract(SampleName, "^[^_]+")) %>% 
  mutate(Species = str_extract(SampleName, 'viridis|concolor|cerberus')) %>% 
  arrange(match(Species, c('viridis', 'concolor', 'cerberus')))

# load regions of interest from bam files
chrom = "scaffold-mi2"
start_pos = floor(8923875 / 100 ) * 100 # Rounds the start coordinate to the nearest 100
end_pos = ceiling(8924426 / 500 ) * 500 # Rounds the end coordinate to the nearest 500
track.df = LoadTrackFile(track.folder = track.folder, format = "bam", # change to bam/bw
                         bamcoverage.path = '/Users/sidgopalan/miniconda3/bin/bamCoverage',
                         meta.info = sample.meta,
                         single.nuc = T, # change to T for BAM, F for BigWig
                         single.nuc.region = paste0(chrom, ':', start_pos, '-', end_pos) # short hand for getting a region, change "region" to single.nuc.region for BAM, "region" for
)

#### Basic plot
ggplot() +
  theme_minimal() +
  geom_col(data = track.df, aes(x = start, y = score, fill = Species)) +
  facet_wrap(~factor(SampleID, levels = sample.meta$SampleID), ncol = 1, strip.position = "right") + 
  guides(fill = guide_legend(nrow = 1, title = NULL)) +
  labs(x = "Position on scaffold-mi2", y = "ATAC-seq read density") +
  ggtitle("Regulatory region")


#### Customized plot
ggplot() +
  theme_minimal() +
  geom_col(data = track.df, aes(x = start, y = score, fill = Species)) +
  facet_wrap(~factor(SampleID, levels = sample.meta$SampleID), ncol = 1, strip.position = "right") + # can use scales = "free_y" if you want to have each track to have independently determined y axes
  scale_x_continuous(limits = c(start_pos,end_pos), expand = c(0,0), labels = label_number(scale = 1e-6), breaks = pretty_breaks(n=3)) +
  scale_y_continuous(breaks = pretty_breaks(n=2)) +
  theme(legend.position = "bottom", # the following 3 lines customize the legend
        legend.direction = "horizontal",
        plot.margin = margin(10,5,5,5),
        ) +
  guides(fill = guide_legend(nrow = 1, title = NULL)) + 
  scale_fill_manual(values = c("viridis" = "#2E8B58", # the following 7 lines customize the fill legend
                               "concolor" = "#D9A528",
                               "cerberus" = "black"),
                    labels = c(expression(italic("C. viridis")),
                               expression(italic("C. o. concolor")),
                               expression(italic("C. o. cerberus"))),
                    breaks = c("viridis", "concolor", "cerberus")) +
  labs(x = paste0("pos. on ", chrom, ' (Mb)'), y = "ATAC-seq read density") +
  ggtitle("Regulatory region")

#### Customize plot, added annotations

# ATAC-seq peak
peaks <-
  data.frame(
    start = c(8923901L),
    end = c(8924400L)
  )

# Gene locus
arrow <- data.frame(
  stringsAsFactors = FALSE,
             start = c(8924000L),
               end = c(8924410L),
          SampleID = c("CV0857") # We want SampleID here so we plot this on the top-most sample
)


# Random transcription factor binding sites (TFBSs) data
TFBSs <- data.frame(
  pos = sample(8923901:8924400, 30, replace = TRUE),
  TF = sample(LETTERS[1:5], 30, replace = TRUE),
  SampleID = sample(sample.meta$SampleID, 30, replace = TRUE)
)

# Hi-C contact
HiC <- data.frame(x1 = 8924000, x2 = 8924410, y1 = 0, y2 = 0, SampleID = "CV1090")

ggplot() +
  theme_minimal() +
  geom_col(data = track.df, aes(x = start, y = score, fill = Species)) +
  annotate(geom = "rect",
           xmin = peaks$start,
           xmax = peaks$end,
           ymin = -Inf,
           ymax = +Inf,
           alpha = 0.2) +
  geom_curve(data = HiC, aes(x = x1, y = y1, xend = x2, yend = y2), curvature = -0.2, color = "darkorchid3") + # make a curve showing a Hi-C contact
  geom_point(data = TFBSs, aes(x = pos, y = 100, color = TF)) + # Add coloured dots representing TFBSs
  geom_gene_arrow(data = arrow, aes(xmin = start, xmax = end, y = 180, fill = 'blue'), arrowhead_height = grid::unit(2, "mm"), arrow_body_height = grid::unit(1, "mm")) + # customize the gene arrow
  facet_wrap(~factor(SampleID, levels = sample.meta$SampleID), ncol = 1, strip.position = "right") + # can use scales = "free_y" if you want to have each track to have independently determined y axes
  scale_x_continuous(limits = c(start_pos,end_pos), expand = c(0,0), labels = label_number(scale = 1e-6), breaks = pretty_breaks(n=3)) + # use pretty_breaks to customize breakpoints in the axes
  scale_y_continuous(breaks = pretty_breaks(n=2)) +
  theme(legend.position = "bottom", # the following 4 lines customize the legend
        legend.direction = "horizontal",
        legend.box = "vertical", # stack the 2 legends vertically
        plot.margin = margin(10,5,5,5),
  ) +
  guides(fill = guide_legend(nrow = 1, title = NULL), # remove legend titles
         color = guide_legend(nrow = 1, title = NULL)) + 
  scale_fill_manual(values = c("viridis" = "#2E8B58", # the following 7 lines customize the fill legend
                               "concolor" = "#D9A528",
                               "cerberus" = "black"),
                    labels = c(expression(italic("C. viridis")),
                               expression(italic("C. o. concolor")),
                               expression(italic("C. o. cerberus"))),
                    breaks = c("viridis", "concolor", "cerberus")) +
  scale_color_manual(values = brewer.pal(5, "Dark2")) + # add custom colours for the TFBSs
  labs(x = paste0("pos. on ", chrom, ' (Mb)'), y = "ATAC-seq read density") +
  ggtitle("Regulatory region")


  