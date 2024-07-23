---
title: "Visualizing read alignment data with ggplot"
sub_title: "Methods to make publication-ready genomic coverage figures in R"
categories:
  - Bioinformatics
elements:
  - ggcoverage
  - ggplot
  - dplyr
  - bam
  - bw
  - R
  - ATAC-seq
  - ChIP-seq
  - WGS
  - mapping
  - read alignment
---
 
 Read coverage plots are a readily interpretable way to visualize genomic or epigenomic profiles (RNA-seq, ChIP-seq, ATAC-seq, WGS, etc.) across many samples mapped to the same reference. See some examples [here](https://www.nature.com/articles/s41588-021-00914-y/figures/4 "Bowfin genome ATAC"), [here](https://genome.cshlp.org/content/32/6/1058/F3.expansion.html "Rattlesnake ChIP-ATAC"), and [here](https://academic.oup.com/view-large/figure/474290094/evae110f5.jpg "Rattlesnake ATAC"). A common tool to visualize genomic data in this manner is [IGV](https://igv.org/), which while versatile, often can be challenging to customize for publication-ready figures. Here is a tutorial on using ggplot and R to have much more artistic control over genomic coverage figures. The structure of intermediate objects will be shown so they can be easily replicated with custom data.

*Disclaimer*: this is not a `ggplot` or `dplyr` tutorial, but rather an example of a `ggplot` application for visualizating quantitative genomic mapping data. 

The R script used in this tutorial can be found here: [genome-coverage-minimal-example.R][1].

[1]:{{ site.url }}{{ site.baseurl }}/assets/blog_files/genome-coverage-vis/genome-coverage-minimal-example.R

<br>

## R packages
There are just two required packages for this tutorial: [`ggcoverage`](https://github.com/showteeth/ggcoverage), [`tidyverse`](https://www.tidyverse.org/).

If you want to do a bit more figure customization: [`scales`](https://scales.r-lib.org/), [`gggenes`](https://github.com/wilkox/gggenes/tree/master) and [`RColorBrewer`](https://cran.r-project.org/web/packages/RColorBrewer/index.html) are used in the accompanying script.

<br>

## Making coverage plots
`ggcoverage` is a package that has its own genomic coverage visualization implementation, but you'll get more flexibility with `ggplot`. It does however have a nice function, used in this guided example, to import many [commonly used alignment data formats](https://rdrr.io/cran/ggcoverage/man/LoadTrackFile.html) (bam, wig, bw, bedgraph, txt) efficiently, without overloading your RAM. 

I'm using [data](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1061517/ "Data") from a [recent publication](https://doi.org/10.1093/gbe/evae110 "GBE 2024") representing chromatin accessibility (ATAC-seq) libraries mapped to a common reference genome, producing .bam alignment files. 

<br>

#### 1. Set up metadata, and load in coverage files 

The first step is to configure a metadata object to let `ggcoverage` know where to look for data. You can include any other useful data attributes that you want to represent in the visualization for grouping or colouring.

```R
# this is a folder with .bam files (multi-sample ATAC-seq mapped to common reference)
track.folder = 'path/to/bams/' 

# Create metadata from .bam names
# sample.meta = data.frame(SampleName = gsub('.bam', '' , grep('bam$', list.files(track.folder), value = T)))
# sample.meta is just a 1 column df with basenames of bam files (no .bam extension), equivalent to:

sample.meta <- data.frame(
  stringsAsFactors = FALSE,
        SampleName = c("CV0857_viridis_North_M", "CV0985_concolor_Other_F",
                       "CV1081_viridis_Mid_M", "CV1082_viridis_South_M", 
                       "CV1086_viridis_South_M", "CV1087_viridis_North_F",
                       "CV1089_viridis_South_M", "CV1090_cerberus_Other_M", 
                       "CV1095_viridis_North_M", "CV1096_viridis_North_F")
)

# Add more attributes to the metadata
sample.meta <- sample.meta %>% 
  mutate(SampleID = str_extract(SampleName, "^[^_]+")) %>% 
  mutate(Species = str_extract(SampleName, 'viridis|concolor|cerberus')) %>% 
  arrange(match(Species, c('viridis', 'concolor', 'cerberus')))
```

The basic format for the metadata, `sample.meta`, should be like shown below. The first column needs to be the file name without the extension, then you can add more attributes as additional columns


| SampleName                                | SampleID | Species |
|:----------------------------------------- |:--------:|:-------:|
| CV0857_viridis_North_M | CV0857   | viridis |
| CV1081_viridis_Mid_M   | CV1081   | viridis |
| CV1082_viridis_South_M | CV1082   | viridis |
| CV1086_viridis_South_M | CV1086   | viridis |
| CV1087_viridis_North_F | CV1087   | viridis |
| CV1089_viridis_South_M | CV1089   | viridis |
| CV1095_viridis_North_M | CV1095   | viridis |
| CV0985_concolor_Other_F | CV0985   | concolor|
| CV1090_cerberus_Other_M | CV1090   | cerberus|


<br>
#### 2. Load regions of interest from alignment files
`ggcoverage` requires minor tweaks to arguments in the `LoadTrackFile` function depending on the input format. I'd use `?LoadTrackFile` to understand how to read your in data correctly. I've shown an example with comments specifying differences for .bam and .bw filetypes. 

In this example, we're going to be focusing on a regulatory region of a gene of interest.

```R
# load regions of interest from bam files (change this for custom data)
chrom = "scaffold-mi2"
start_pos = floor(8923875 / 100 ) * 100 # Rounds the start coordinate to the nearest 100
end_pos = ceiling(8924426 / 500 ) * 500 # Rounds the end coordinate to the nearest 500

track.df = LoadTrackFile(track.folder = track.folder, format = "bam", # change between bam/bw
                         bamcoverage.path = '/Users/sidgopalan/miniconda3/bin/bamCoverage', 
                         meta.info = sample.meta,
                         single.nuc = T, # change to T for BAM, F for BigWig
                         single.nuc.region = paste0(chrom, ':', start_pos, '-', end_pos) # bam, if single.nuc = T
                         #region = paste0(chrom, ':', start_pos, '-', end_pos) # bw
)
```
The resulting `track.df` object will provide scores (RPKM normalized by default, this can be changed) representing levels of coverage across genomic coordinates, either nucleotide-by-nucleotide or by larger bins depending on the options used in `LoadTrackFile`. Note that the attributes we added to `sample.meta`, *SampleID* and *Species*, are also included in this object. This will allow us to group and colour our plot. You may notice that this dataframe is already in [`tidy`](https://cran.r-project.org/web/packages/tidyr/vignettes/tidy-data.html) format, ideal for input to `ggplot`. See below for a snapshot of `track.df`: 

| seqnames     |   start   |   end     | score |  SampleID  | Species  |
|:-------------|:---------:|:---------:|:-----:|:------:|:------:|
| scaffold-mi2 |  8924084  |  8924085  |   93  | CV0857 | viridis|
| scaffold-mi2 |  8924085  |  8924086  |   95  | CV0857 | viridis|
| scaffold-mi2 |  8924086  |  8924087  |   92  | CV0857 | viridis|
| scaffold-mi2 |  8924087  |  8924088  |   93  | CV0857 | viridis|
| scaffold-mi2 |  8924088  |  8924089  |   91  | CV0857 | viridis|
| scaffold-mi2 |  8924089  |  8924090  |   90  | CV0857 | viridis|
| scaffold-mi2 |  8924090  |  8924091  |   87  | CV0857 | viridis|
| scaffold-mi2 |  8924091  |  8924092  |   77  | CV0857 | viridis|
| scaffold-mi2 |  8924092  |  8924093  |   77  | CV0857 | viridis|
| scaffold-mi2 |  8924093  |  8924094  |   77  | CV0857 | viridis|


<br>
#### 3. Making a basic coverage plot

First is a demonstration of the most basic plot, then we can increase its complexity by building off the starter code.

The data used here represents ATAC-seq data from ten individuals of three species (*viridis*, *concolor* and *cerberus*). So for these plots, we want to split up coverage tracks to show each individual, and colour by species. This process involves ["facet wrapping"](https://ggplot2.tidyverse.org/reference/facet_wrap.html). Here, we facet wrap by SampleID, then set `fill=Species`. We use `geom_col` because [we want heights to represent data values](https://ggplot2.tidyverse.org/reference/geom_bar.html). 

```R
# Basic plot
ggplot() +
  theme_minimal() +
  geom_col(data = track.df, aes(x = start, y = score, fill = Species)) +
  facet_wrap(~factor(SampleID, levels = sample.meta$SampleID)) +
  labs(x = "Position on scaffold-mi2", y = "ATAC-seq read density") +
  ggtitle("Regulatory region")
```
![center-aligned-image]({{ site.url }}{{ site.baseurl }}/assets/blog_files/genome-coverage-vis/basic_plot.png){: .align-center}

Not quite publication quality, so we can make a few more customizations.

<br>

#### 4. More cutomization with `ggplot`
`ggplot` allows a great deal of customization and complexity. For example, we can add the following:
- Single column facet
- Changed colours, axes labels and legend (requires `scales`)
<br>

```R
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
```

<div style="text-align: center;">
  <img src="{{ site.url }}{{ site.baseurl }}/assets/blog_files/genome-coverage-vis/intermediate_plot@2x.png" class="align-center" alt="center-aligned-image">
  <figcaption>Slightly more customized plot.</figcaption>
</div>

<br>

#### 5. Adding annotations to ggplot

You may want to add data annotation to your plot, representing other biological/genomic information. Adding annotations to ggplot is easy; it's just a matter of creating new data frames with the information you want to plot, then setting the appropriate x, y, and potentially fill and color values in the aes() calls to the various geoms.

Lets add the following (*this data is made up for the sake of illustration*):
- A translucent gray region representing a peak region (such as those called by peak-calling algorithms)
- An arrow representing the gene locus (uses `gggenes`)
- Transcription factors binding sites in each individual (coloured with `RColorBrewer`)
- A Hi-C contact curve

```R
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
```

<div style="text-align: center;">
  <img src="{{ site.url }}{{ site.baseurl }}/assets/blog_files/genome-coverage-vis/complex_plot@2x.png" class="align-center" alt="center-aligned-image">
  <figcaption>Even more customized plot, with four types of annotated data.</figcaption>
</div>

<br>

## Additional steps
Ultimately `ggplot` offers far more customization potential than IGV, making it more useful for making publication-ready figures. 

Because the `track.df` object imported by `ggcoverage::LoadTrackFile` is already `tidy` formatted, this makes it incredibly easy to augment and manipulate with `dplyr`. For instance, you can calculate coverage statistics across some or all samples, and plot that in addition to coverage by creating an additional variable to facet by. You can also add other plots using [`cowplot`](https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html), such as expression profiles, right next to coverage plots. See a demonstration of these figures [here](https://academic.oup.com/view-large/figure/474290098/evae110f6.jpg).

There are limits to `ggplot`, and there may be the need to touch up figures in other vector graphics tools such as Adobe Illustrator or InkScape.