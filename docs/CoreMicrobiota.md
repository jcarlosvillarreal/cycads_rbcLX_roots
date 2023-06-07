Core microbiota of Zamia species
================

### Zamia Coralloid root Core Cyanobacteria

### Contributors

Juan Carlos Villarreal A.

Adriel M. Sierra

We use the two dataset generated with two amplicon markers to identify
the core microbiota.

First we use a subset dataset with only ASVs whose phylum corresponds to
Cyanobacteria in the 16S dataset.

``` r
knitr::opts_knit$set(root.dir = "/Users/adrielsierra/Documents/GitHub/cycads_rbcLX_roots")
root16ps = readRDS("./Data/Zamia_root_bacteria16s_ps_mai23.rds")
ps_cyano = subset_taxa(root16ps, Phylum=="Cyanobacteria")
```

Second, we use the cyanobacterial community identified with the amplicon
marker rbcL-X.

``` r
knitr::opts_knit$set(root.dir = "/Users/adrielsierra/Documents/GitHub/cycads_rbcLX_roots")
rbclx.ps = readRDS("./Data/RBCLX_ps_rootZamia.rds")
rbclx.ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 175 taxa and 25 samples ]
    ## sample_data() Sample Data:       [ 25 samples by 3 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 175 taxa by 6 taxonomic ranks ]

**16S**

*Zamia pseudoparasitica*

``` r
cyano.zpseudo = subset_samples(ps_cyano, Species %in% c("Zamia pseudoparasitica"))
rank.pseudo16 = rankabundance(otu_table(cyano.zpseudo))
rankabunplot(rank.pseudo16, labels=" ",scale="logabun",  lty=2, pch= 20)
```

    ## Warning in xy.coords(x, y, xlabel, ylabel, log): 20 y values <= 0 omitted from
    ## logarithmic plot

![](CoreMicrobiota_files/figure-gfm/rank16S_zp-1.png)<!-- -->

*Zamia nana*

``` r
cyano.znana = subset_samples(ps_cyano, Species %in% c("Zamia nana"))
rank.nana16 = rankabundance(otu_table(cyano.znana))
rankabunplot(rank.nana16, labels=" ",scale="logabun",  lty=2, pch= 20)
```

    ## Warning in xy.coords(x, y, xlabel, ylabel, log): 76 y values <= 0 omitted from
    ## logarithmic plot

![](CoreMicrobiota_files/figure-gfm/rank16S_zn-1.png)<!-- -->

**rbcL-X**

*Zamia pseudoparasitica*

``` r
ps_pseudo = subset_samples(rbclx.ps, Species_name %in% c("Zamia pseudoparasitica"))
rank.pseudo = rankabundance(t(otu_table(ps_pseudo)))
rankabunplot(rank.pseudo, labels=" ",scale="logabun",  lty=2, pch= 20)
```

    ## Warning in xy.coords(x, y, xlabel, ylabel, log): 34 y values <= 0 omitted from
    ## logarithmic plot

![](CoreMicrobiota_files/figure-gfm/rankrbclxzp-1.png)<!-- -->

*Zamia nana*

``` r
ps_nana = subset_samples(rbclx.ps, Species_name %in% c("Zamia nana"))
rank.nana = rankabundance(t(otu_table(ps_nana)))
rankabunplot(rank.nana, labels=" ",scale="logabun",  lty=2, pch= 20)
```

    ## Warning in xy.coords(x, y, xlabel, ylabel, log): 104 y values <= 0 omitted from
    ## logarithmic plot

![](CoreMicrobiota_files/figure-gfm/rankrbclx_zn-1.png)<!-- -->

### Core microbiome.

Indicator species analysis, the labdsv package (version 2.0-1) was used
to identify which ASVs were relatively more abundant and predominantly
found to one of the specific Zamia species, in contrast to the other
host species.

##### The script is adapted from <https://bocasbiome.github.io/wf1.html>

The Indicator species analysis was done for both marker. A core
microbiota was only observed when we analyzed the dataset generated from
the rbcL-X marker. Therefore, here is only presented for the rbcL-X.

``` r
zamia.host = data.frame(t(otu_table(rbclx.ps)))
cyano.ASV = tibble::remove_rownames(zamia.host)
host.group = data.frame(sample_data(rbclx.ps)) %>%
  select(Species_name)

host.group$Status = host.group$Species
host.group = host.group[, c(2,1)]
host.group$Status = str_replace(host.group$Status, "Zamia nana", "1")
host.group$Status = str_replace(host.group$Status, "Zamia pseudoparasitica", "2")
host.group = tibble::rownames_to_column(host.group, "Label")
host.group$Status =  as.integer(host.group$Status)
host.group$Species =  as.character(host.group$Species)
```

Here we calculate the indicator values.

``` r
set.seed(1280)
iva = indval(cyano.ASV, host.group$Status)

gr = iva$maxcls[iva$pval <= 0.01]
iv = iva$indcls[iva$pval <= 0.01]
pv = iva$pval[iva$pval <= 0.01]
fr = apply(cyano.ASV > 0, 2, sum)[iva$pval <= 0.01]
indval.out = data.frame(group = gr, indval = iv, pval = pv, freq = fr)
indval.out = indval.out[order(indval.out$group, -indval.out$indval),]
indval.out
```

    ##       group  indval  pval freq
    ## ASV_2     2 0.74343 0.004   21

``` r
indval.out$prob.corrected = p.adjust(indval.out$pval, "bonferroni")

indval.out
```

    ##       group  indval  pval freq prob.corrected
    ## ASV_2     2 0.74343 0.004   21          0.004

### core microbiome using threshold frequency

For the threshold frequency base method, we used the microbiome R
package (Lahti and Shetty 2012). Due to the lack of consensus about a
fixed threshold (Risely 2020). We tested for the presence of a ‘core’
microbiota as defined by any taxon with a prevalence higher than 50%,
75% or 90% (Jorge et al. 2020; Neu et al. 2021). We pruned out the
phyloseq object to retain those samples with more than 50 reads. We
detected the core microbiota with the function (core_members) and
estimated the total abundance of specific members of the core microbiota
in each sample (sample_sums).

To allow for comparisons of the core microbiomes, we considered samples
that we have amplicon sequencing data for the two markers and from the
same locality.

``` r
ps_cyano_sub = subset_samples(ps_cyano, Locality !="El Cope")
ps_cyano_prun = prune_taxa(taxa_sums(ps_cyano_sub)!=0, ps_cyano_sub)
```

First, we pruned samples with less than 50 reads, and transform
incidence of ASVs in the samples to compositional relative abundance.

``` r
which(!rowSums(otu_table(rbclx.ps)) > 50)
```

    ## ASV_145 ASV_146 ASV_147 ASV_148 ASV_149 ASV_150 ASV_151 ASV_152 ASV_153 ASV_154 
    ##     145     146     147     148     149     150     151     152     153     154 
    ## ASV_155 ASV_156 ASV_157 ASV_158 ASV_159 ASV_160 ASV_161 ASV_162 ASV_163 ASV_164 
    ##     155     156     157     158     159     160     161     162     163     164 
    ## ASV_166 ASV_167 ASV_168 ASV_170 ASV_171 ASV_172 ASV_173 ASV_176 ASV_177 ASV_180 
    ##     165     166     167     168     169     170     171     172     173     174 
    ## ASV_183 
    ##     175

``` r
total = median(sample_sums(rbclx.ps))
standf = function(x, t=total) round(t * (x / sum(x)))
rbclx.ps = transform_sample_counts(rbclx.ps, standf)

ps2prun.ra = microbiome::transform(rbclx.ps, "compositional")

ps2prun.ra2 = prune_taxa(taxa_sums(ps2prun.ra) > 0, ps2prun.ra)
taxonomy = as.data.frame(tax_table(ps2prun.ra2))
```

We first extract the taxonomy table. Subset this taxonomy table to
include only core ASVs with the different threshold. Total core
abundance in each sample (sum of abundances of the core members).

``` r
core50_name = core_members(ps2prun.ra2, detection = 0.0001, prevalence = 50/100, include.lowest = FALSE)

core50_seq = core_members(ps2prun.ra2, detection = 0.0001, prevalence = 50/100, include.lowest = FALSE)

core.abundance.50 = sample_sums(core(ps2prun.ra2, detection = 0.0001, prevalence = 50/100))
DT::datatable(as.data.frame(core.abundance.50))
```

<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-6d9344b01b77ad3b2ecf" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-6d9344b01b77ad3b2ecf">{"x":{"filter":"none","vertical":false,"data":[["Zn-1R","Zn-25R","Zn-26R","Zn-27R","Zn-28R","Zn-3R","Zn-4R","Zn-5R","Zn-6R","Zn-7R","Zp-10R","Zp-11R","Zp-12R","Zp-13R","Zp-14R","Zp-15R","Zp-16R","Zp-17R","Zp-20R","Zp-21R","Zp-22R","Zp-23R","Zp-24R","Zp-8R","Zp-9R"],[1,0.644400537603529,0.944619307999816,0.98595092588338,0.0413316178835638,0.998460672479351,0.999471586275028,0.996151725405505,0.0541400542204659,0.355048074161718,0.964595467025077,0.927398881115668,0.992004686907675,0.989397507322957,0.994888171573639,0.999138445986307,0.997300433079459,0.967640402972902,0.95318728532206,0.933808916612101,0.997197013210798,0.986306877577512,0.990407921792972,0.939266185728071,0.973992854927458]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>core.abundance.50<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":1},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

``` r
core_taxa_id50 = subset(taxonomy, rownames(taxonomy) %in% core50_name)
DT::datatable(core_taxa_id50)
```

<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-692e7cb02b354cadf993" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-692e7cb02b354cadf993">{"x":{"filter":"none","vertical":false,"data":[["ASV_1","ASV_2","ASV_3","ASV_6","ASV_7","ASV_10","ASV_11","ASV_12","ASV_14","ASV_15","ASV_17","ASV_18"],["Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria"],["Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria"],["Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria","Cyanobacteria"],["Nostocales","Nostocales","Nostocales","Nostocales","Nostocales","Nostocales","Nostocales","Nostocales","Nostocales","Nostocales","Nostocales","Nostocales"],["Nostocaceae","Nostocaceae","Nostocaceae","Nostocaceae","Nostocaceae","Nostocaceae","Nostocaceae","Nostocaceae","Nostocaceae","Nostocaceae","Nostocaceae","Nostocaceae"],["ASV_1","ASV_2 _KX92289_P430","ASV_3","ASV_6","ASV_7 _MN86557","ASV_10","ASV_11","ASV_12","ASV_14","ASV_15","ASV_17","ASV_18"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>Family<\/th>\n      <th>Genus<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

``` r
core75_name = core_members(ps2prun.ra2, detection = 0.0001, prevalence = 75/100, include.lowest = FALSE)

core75_seq = core_members(ps2prun.ra2, detection = 0.0001, prevalence = 75/100, include.lowest = FALSE)

core.abundance.75 = sample_sums(core(ps2prun.ra2, detection = 0.0001, prevalence = 75/100))
DT::datatable(as.data.frame(core.abundance.75))
```

<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-ebdd61aec8c7a3f94c10" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-ebdd61aec8c7a3f94c10">{"x":{"filter":"none","vertical":false,"data":[["Zn-1R","Zn-25R","Zn-26R","Zn-27R","Zn-28R","Zn-3R","Zn-4R","Zn-5R","Zn-6R","Zn-7R","Zp-10R","Zp-11R","Zp-12R","Zp-13R","Zp-14R","Zp-15R","Zp-16R","Zp-17R","Zp-20R","Zp-21R","Zp-22R","Zp-23R","Zp-24R","Zp-8R","Zp-9R"],[0.634953820704866,0.636290535650695,0.892340210448927,0.968559022193631,0.0404815512567201,0.998460672479351,0.999471586275028,0.996151725405505,0,6.89235293441926e-05,0.941792742185615,0.916497225764207,0.975692410196322,0.977818620412383,0.984342871583978,0.999138445986307,0.991349898335459,0.945286204955602,0.93079759675585,0.932545289542912,0.989040781160253,0.968650561165294,0.975416709744862,0.824403804622524,0.967514043169104]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>core.abundance.75<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":1},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

``` r
core_taxa_id75 = subset(taxonomy, rownames(taxonomy) %in% core75_name)
DT::datatable(core_taxa_id75)
```

<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-273bf8fd9d1a96e8d711" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-273bf8fd9d1a96e8d711">{"x":{"filter":"none","vertical":false,"data":[["ASV_1","ASV_2","ASV_3"],["Bacteria","Bacteria","Bacteria"],["Cyanobacteria","Cyanobacteria","Cyanobacteria"],["Cyanobacteria","Cyanobacteria","Cyanobacteria"],["Nostocales","Nostocales","Nostocales"],["Nostocaceae","Nostocaceae","Nostocaceae"],["ASV_1","ASV_2 _KX92289_P430","ASV_3"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>Family<\/th>\n      <th>Genus<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

``` r
core90_name = core_members(ps2prun.ra2, detection = 0.0001, prevalence = 90/100, include.lowest = FALSE)

core90_seq = core_members(ps2prun.ra2, detection = 0.0001, prevalence = 90/100, include.lowest = FALSE)

core.abundance.90 = sample_sums(core(ps2prun.ra2, detection = 0.0001, prevalence = 90/100))
DT::datatable(as.data.frame(core.abundance.90))
```

<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-7f40d5d3eeb2d5b5fe0b" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-7f40d5d3eeb2d5b5fe0b">{"x":{"filter":"none","vertical":false,"data":[["Zn-1R","Zn-25R","Zn-26R","Zn-27R","Zn-28R","Zn-3R","Zn-4R","Zn-5R","Zn-6R","Zn-7R","Zp-10R","Zp-11R","Zp-12R","Zp-13R","Zp-14R","Zp-15R","Zp-16R","Zp-17R","Zp-20R","Zp-21R","Zp-22R","Zp-23R","Zp-24R","Zp-8R","Zp-9R"],[0.07658640812388,0.19907412725581,0.169565317281625,0.41722189036438,0.0194251711620641,0.998265384659567,0.999471586275028,0.996151725405505,0,0,0.342822023871064,0.136150072945744,0.474342626736051,0.517351099879387,0.573547149437699,0.00368745117860589,0.721106018311105,0.299886276176582,0.332295603625544,0.83038678475836,0.483905801263642,0.303327934199492,0.521958392206867,0.000206772963286312,0.353497294751473]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>core.abundance.90<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":1},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

``` r
core_taxa_id90 = subset(taxonomy, rownames(taxonomy) %in% core90_name)
DT::datatable(core_taxa_id90)
```

<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-ccafefc9f9f7ba425490" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-ccafefc9f9f7ba425490">{"x":{"filter":"none","vertical":false,"data":[["ASV_1"],["Bacteria"],["Cyanobacteria"],["Cyanobacteria"],["Nostocales"],["Nostocaceae"],["ASV_1"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Kingdom<\/th>\n      <th>Phylum<\/th>\n      <th>Class<\/th>\n      <th>Order<\/th>\n      <th>Family<\/th>\n      <th>Genus<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

Visualitation of the core microbiome with Core line plots considering
the compositional (relative) abundances.

``` r
det = c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences = seq(.05, 1, .05)

plot_core(ps2prun.ra2, prevalences = prevalences, 
          detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)") + 
  theme_bw()
```

    ## Warning: Transformation introduced infinite values in continuous x-axis
    ## Transformation introduced infinite values in continuous x-axis

![](CoreMicrobiota_files/figure-gfm/core_plot-1.png)<!-- -->

``` r
prevalences = seq(.05, 1, .05)
detections = round(10^seq(log10(1e-3), log10(.2), length = 10), digits = 4)

p.core = plot_core(ps2prun.ra2, 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(5, "Spectral")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .5) + theme_classic() +
xlab("Detection Threshold (Relative Abundance (%))")

p.core 
```

![](CoreMicrobiota_files/figure-gfm/core_plot-2.png)<!-- -->
