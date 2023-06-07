Analyses of bacterial diversity
================

### Zamia Coralloid root Cyanobacteria community

### Contributors

Juan Carlos Villarreal A.

Adriel M. Sierra

To study the coralloid root cyanobacteria community of two Zamia
species, we used the processed *phyloseq* object with the whole
bacterial community identified using the 16S marker. From this dataset,
we subset all ASVs whose phylum corresponds to Cyanobacteria.

The dataset includes a total of forty-four coralloid roots of *Z.
pseudoparasitica* from two localities, and ten of *Zamia nana* from one
locality.

``` r
knitr::opts_knit$set(root.dir = "/Users/adrielsierra/Documents/GitHub/cycads_rbcLX_roots")
root16ps = readRDS("./Data/Zamia_root_bacteria16s_ps_mai23.rds")
ps_cyano = subset_taxa(root16ps, Phylum=="Cyanobacteria")
ps_cyano
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 101 taxa and 53 samples ]
    ## sample_data() Sample Data:       [ 53 samples by 5 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 101 taxa by 7 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 101 tips and 100 internal nodes ]

``` r
knitr::opts_knit$set(root.dir = "/Users/adrielsierra/Documents/GitHub/cycads_rbcLX_roots")
rbclx.ps = readRDS("./Data/RBCLX_ps_rootZamia.rds")
rbclx.ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 175 taxa and 25 samples ]
    ## sample_data() Sample Data:       [ 25 samples by 3 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 175 taxa by 6 taxonomic ranks ]

Analyses of Cyanobacterial diversity using the 16S and rbcL-X markers

### Compute and visualize alpha diversity

Alpha diversity was assessed using four indexes: i) the Chao index which
shows the observed richness, ii) the exponential of Shannon which
weights the ASVs by their frequency, iii) the Simpson’s diversity index
which considers the presence and relative abundances of ASVs present in
a sample, and iv) the Pielou index which measures diversity along with
species richness.

To the data frame we add an extra column with read count by sample, the
Zamia host species, and locality respective to each sample.

``` r
richness.df.cyano = estimate_richness(ps_cyano, split = TRUE, measures = c("Observed","Shannon", "Simpson"))
richness.df.cyano$Pielou = richness.df.cyano$Shannon/log(specnumber(otu_table(ps_cyano)))
richness.df.cyano$readcount = sample_sums(ps_cyano)
richness.df.cyano$Host = as.factor(ps_cyano@sam_data$Species)
richness.df.cyano$Locality = as.factor(ps_cyano@sam_data$Locality)
```

``` r
richness.rbclx.df = estimate_richness(rbclx.ps, split = TRUE, measures = c("Observed","Shannon", "Simpson"))
```

    ## Warning in estimate_richness(rbclx.ps, split = TRUE, measures = c("Observed", : The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

``` r
richness.rbclx.df$Pielou = richness.rbclx.df$Shannon/log(specnumber(t(otu_table(rbclx.ps))))
richness.rbclx.df$readcount = sample_sums(rbclx.ps)
richness.rbclx.df$Host = as.factor(rbclx.ps@sam_data$Species_name)
richness.rbclx.df$Locality = as.factor(rbclx.ps@sam_data$Site)
```

### Alpha diversity table : Richness indexes estimates of the cyanobacteria community identified with the 16S marker in coralloid roots of the two species *Z. nana* and *Z. pseudoparasitica*.

``` r
DT::datatable(richness.df.cyano) %>%
    formatRound(columns=c('Shannon', 'Simpson', 'Pielou'), digits=4)
```

<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-eb0b299c97b0fe631641" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-eb0b299c97b0fe631641">{"x":{"filter":"none","vertical":false,"data":[["X10S","X11S","X12S","X13S","X14S","X15S","X16S","X17S","X18S","X19S","X1S","X20S","X22S","X23S","X25S","X26S","X27S","X28S","X2S","X3S","X4S","X5S","X6S","X7S","X8S","X9S","ZAP02","ZAP10","ZAP12","ZAP14","ZAP18","ZAP25","ZAP26","ZAP27","ZAP28","ZAP29","ZAP30","ZAP31","ZAP32","ZAP33","ZAP34","ZAP36","ZAP37","ZAP39","ZAP40","ZAP41","ZAP43","ZAP44","ZAP45","ZAP46","ZAP47","ZAP51","ZAP55"],[6,8,6,5,9,11,5,4,7,4,4,7,8,4,1,5,4,8,4,5,5,0,0,0,4,10,3,2,5,1,2,5,2,6,1,1,3,3,3,7,5,5,6,6,2,7,2,7,2,1,2,2,1],[1.23151616046822,1.62870065785748,1.20551777217355,1.24525964504716,1.29082594328878,1.29742206113786,1.23824383328713,1.05875835050612,1.26625081514199,1.32711964607567,0.901097602859054,1.64754896367754,1.27582569961181,1.06204654051338,0,0.731843181930048,0.915560059790578,1.13134075398397,0.830485787497943,1.1517006792885,1.14982288606581,0,0,0,1.0771568754573,1.72566918982964,0.591426563052181,0.15159848960351,1.18246148887696,0,0.0141100040118955,0.216856722017325,0.185538809311542,1.31802335544431,0,0,0.360722799185723,0.872607596044421,0.0733341600444549,0.134881525005117,0.201442340275103,0.095268130025434,0.595736156238222,0.0781771709741379,0.0100904150601421,0.429292653407039,0.0162049676310515,0.0570675286855453,0.138558987317699,0,0.00330475011514022,0.0260027687599942,0],[0.645939221043004,0.764602531464018,0.641384845598984,0.656908584075497,0.651997561338661,0.640520518406902,0.649900357607682,0.599915487201957,0.661275472266613,0.721738271604938,0.554057875388821,0.779386660584234,0.652995420244125,0.60940815555836,0,0.333645764442138,0.476704518756363,0.499153002957055,0.500561043382718,0.65220491306183,0.650488468168212,1,1,1,0.645513284611707,0.79979419566131,0.341304278160551,0.0674849625898577,0.645318559556787,0,0.003890495155497,0.0790719156468261,0.0871541460770209,0.673366152650308,0,0,0.190255469627867,0.538106625532808,0.0242231493011235,0.0431110837340427,0.0769376957827524,0.0332810400898567,0.27468588480505,0.0229918999203644,0.00264239114748399,0.188286917646158,0.00456827208200838,0.0170983945227124,0.060273029781531,0,0.000742528206467918,0.00794101592738272,0],[0.687322255926904,0.783239454061191,0.672812279146443,0.7737233200651,0.58748020416452,0.541067024842895,0.769364151124303,0.763732710887482,0.650724194925814,0.957314466029841,0.650004521500845,0.846672681406891,0.613542469956219,0.766104638595965,null,0.454719735552405,0.660436978947935,0.544059898442752,0.599068863576024,0.71559186619798,0.714425127668847,-0,-0,-0,0.777004441240845,0.749448606733467,0.538339657359186,0.218710389157242,0.734704631810626,null,0.020356432814885,0.134740657183447,0.267675920086205,0.735602840716202,null,null,0.32834404175748,0.794281663372178,0.0667516291241933,0.0693153921164339,0.125163163312366,0.0591934173349696,0.332486679417346,0.0436315098743801,0.0145573917677784,0.220612782977389,0.0233788264390842,0.0293269083946399,0.199898363873842,null,0.00476774660249002,0.037514065539426,null],[2999,167,3269,2443,4579,5137,2299,3445,2249,225,206,2861,1719,2854,1441,2715,565,3063,4605,4873,16692,0,0,0,2420,2481,395,429,95,159,7183,1746,219,994,72,53,492,662,1717,14034,9959,7343,4341,20631,6803,2251,2184,9412,6848,14398,5385,8529,699],["Zamia pseudoparasitica","Zamia nana","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia nana","Zamia nana","Zamia nana","Zamia pseudoparasitica","Zamia nana","Zamia nana","Zamia nana","Zamia nana","Zamia nana","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica"],["Donoso","El Valle, Cerro Gaital","Donoso","Donoso","Donoso","Donoso","Donoso","El Valle, Cerro Gaital","El Valle, Cerro Gaital","El Valle, Cerro Gaital","Donoso","El Valle, Cerro Gaital","El Valle, Cerro Gaital","El Valle, Cerro Gaital","El Valle, Cerro Gaital","El Valle, Cerro Gaital","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Observed<\/th>\n      <th>Shannon<\/th>\n      <th>Simpson<\/th>\n      <th>Pielou<\/th>\n      <th>readcount<\/th>\n      <th>Host<\/th>\n      <th>Locality<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"targets":2,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 4, 3, \",\", \".\", null);\n  }"},{"targets":3,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 4, 3, \",\", \".\", null);\n  }"},{"targets":4,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 4, 3, \",\", \".\", null);\n  }"},{"className":"dt-right","targets":[1,2,3,4,5]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":["options.columnDefs.0.render","options.columnDefs.1.render","options.columnDefs.2.render"],"jsHooks":[]}</script>

### Alpha diversity table : Richness indexes estimates of the cyanobacteria community identified with the rbcL-X marker in coralloid roots of the two species *Z. nana* and *Z. pseudoparasitica*.

``` r
DT::datatable(richness.rbclx.df)%>%
    formatRound(columns=c('Shannon', 'Simpson', 'Pielou'), digits=4)
```

<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-463b650862ea34ea3afc" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-463b650862ea34ea3afc">{"x":{"filter":"none","vertical":false,"data":[["Zn.1R","Zn.25R","Zn.26R","Zn.27R","Zn.28R","Zn.3R","Zn.4R","Zn.5R","Zn.6R","Zn.7R","Zp.10R","Zp.11R","Zp.12R","Zp.13R","Zp.14R","Zp.15R","Zp.16R","Zp.17R","Zp.20R","Zp.21R","Zp.22R","Zp.23R","Zp.24R","Zp.8R","Zp.9R"],[4,33,26,33,14,4,4,4,5,5,31,71,27,24,24,5,18,30,25,9,19,25,27,19,24],[0.891441786532632,1.47159416005061,1.55955007031733,1.16509253071246,0.299346842431142,0.0145280699749611,0.00489491124997962,0.0294248000856802,0.399060246470943,1.1956296032135,1.21487742364917,1.55293979109531,1.09923218875263,1.07389739717025,0.955114660962273,0.0374454731270453,0.757578511952533,1.33978851912582,1.36915160985682,0.614837171790981,0.95658685177063,1.13054293310144,1.11621186522196,0.914251196022902,1.09222073171069],[0.549273833659806,0.714167365734258,0.706479952376432,0.616548825389151,0.109116314993521,0.0034899516841157,0.0010432912115812,0.00767308278726386,0.175467631947463,0.68312418510185,0.59115904731978,0.676079792355861,0.600088591326476,0.587472066446833,0.542287002228107,0.0103971675958732,0.424025452588679,0.658541671576696,0.678991196788473,0.297788475068478,0.563797595500054,0.593401885277932,0.600696836296099,0.317915371657297,0.575605508745944],[0.643039322335913,0.420875440473352,0.478669079358952,0.333216076393613,0.113429457962849,0.0104797872532821,0.00353093209296859,0.02122550658138,0.247950072126366,0.742886441270195,0.353780416555546,0.364310676825488,0.333521419124507,0.337910386197793,0.300534450310964,0.0232661805949464,0.262104177389263,0.393916720725884,0.425350862956286,0.279824455876227,0.324879156450182,0.351222909677708,0.33867327498378,0.310500982527314,0.343675969639291],[71514,79551,38287,105238,104016,85283,93896,93853,73607,89205,90140,93740,86399,85548,91586,77595,92333,76841,68986,63312,89516,109232,97686,87052,80694],["Zamia nana","Zamia nana","Zamia nana","Zamia nana","Zamia nana","Zamia nana","Zamia nana","Zamia nana","Zamia nana","Zamia nana","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica"],["Gaital","Gaital","Gaital","Gaital","Gaital","Gaital","Gaital","Gaital","Gaital","Gaital","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Observed<\/th>\n      <th>Shannon<\/th>\n      <th>Simpson<\/th>\n      <th>Pielou<\/th>\n      <th>readcount<\/th>\n      <th>Host<\/th>\n      <th>Locality<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"targets":2,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 4, 3, \",\", \".\", null);\n  }"},{"targets":3,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 4, 3, \",\", \".\", null);\n  }"},{"targets":4,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 4, 3, \",\", \".\", null);\n  }"},{"className":"dt-right","targets":[1,2,3,4,5]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":["options.columnDefs.0.render","options.columnDefs.1.render","options.columnDefs.2.render"],"jsHooks":[]}</script>

### Correlation of richness with read count

Here, we evaluated if sequencing depth bias toward over or
underestimation of cyanobacterial community diversity by correlating the
observed species richness and the read count by sample (coverage) with
the Kendall method.

**16S**

``` r
cor.16s.cyano = ggscatter(richness.df.cyano, x = "readcount", y = "Observed", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "kendall",
          xlab = "Read count", ylab = "Observed species")
cor.16s.cyano
```

![](Cyanobacterial_Community_files/figure-gfm/cor.16S-1.png)<!-- -->

**rbcL-X**

``` r
cor.rbclx = ggscatter(richness.rbclx.df, x = "readcount", y = "Observed", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "kendall",
          xlab = "Read count", ylab = "Observed species")
cor.rbclx
```

![](Cyanobacterial_Community_files/figure-gfm/cor.rbclx-1.png)<!-- -->

### Rarefaction curves

For the two Zamia species separately, species accumulation curves (100
permutations) were generated from sample-based rarefactions using the
Coleman method (Coleman et al. 1982) with the function specaccum in the
vegan package (Oksanen et al. 2016).

**16S**

``` r
cyano.zpseudo <- subset_samples(ps_cyano, Species %in% c("Zamia pseudoparasitica"))
cyano.znana <- subset_samples(ps_cyano, Species %in% c("Zamia nana"))

rarezpseudo.cyano<-specaccum(otu_table(cyano.zpseudo), method = "coleman", permutations = 100, conditioned =TRUE, gamma = "coleman")
plot(rarezpseudo.cyano, main="Zamia pseudoparasitica", font.main= 4, add = F, ci = 2, ci.type = "polygon", col = "black",  lwd=1.9, ci.lty=0, ci.col = rgb(blue=0.1, green=0, red=0, alpha =0.3), xlab = "Number of Samples", ylab="Number of ASVs")
```

![](Cyanobacterial_Community_files/figure-gfm/rare16S-1.png)<!-- -->

``` r
rareznana.cyano<-specaccum(otu_table(cyano.znana), method = "coleman", permutations = 100, conditioned =TRUE, gamma = "coleman")
plot(rareznana.cyano, main= "Zamia nana", font.main= 4, add = F, ci = 2, ci.type = "polygon", col = "black",  lwd=1.9, ci.lty=0, ci.col = rgb(blue=0.1, green=0, red=0, alpha =0.3), xlab = "Number of Samples", ylab="Number of ASVs")
```

![](Cyanobacterial_Community_files/figure-gfm/rare16S-2.png)<!-- -->

**rbcL-X**

``` r
ps_pseudo <- subset_samples(rbclx.ps, Species_name %in% c("Zamia pseudoparasitica"))
ps_nana <- subset_samples(rbclx.ps, Species_name %in% c("Zamia nana"))

rarezpseudo.rbclx <-specaccum(t(otu_table(ps_pseudo)), method = "coleman", permutations = 100, conditioned =TRUE, gamma = "coleman")
plot(rarezpseudo.rbclx, main="Zamia pseudoparasitica", font.main= 4, add = F, ci = 2, ci.type = "polygon", col = "black",  lwd=1.9, ci.lty=0, ci.col = rgb(blue=0.1, green=0, red=0, alpha =0.3), xlab = "Number of Samples", ylab="Number of ASVs")
```

![](Cyanobacterial_Community_files/figure-gfm/rare.rbclx-1.png)<!-- -->

``` r
rareznana.rbclx <-specaccum(t(otu_table(ps_nana)), method = "coleman", permutations = 100, conditioned =TRUE, gamma = "coleman")
plot(rareznana.rbclx, main= "Zamia nana", font.main= 4, add = F, ci = 2, ci.type = "polygon", col = "black",  lwd=1.9, ci.lty=0, ci.col = rgb(blue=0.1, green=0, red=0, alpha =0.3), xlab = "Number of Samples", ylab="Number of ASVs")
```

![](Cyanobacterial_Community_files/figure-gfm/rare.rbclx-2.png)<!-- -->

### Extrapolate species richness

Estimated species richness was also extrapolated for each species with
three non-parametric methods. The Chao Index with standard error (Chao
1987) was estimated using the small sampling corrections to reduce
estimation bias by multiplying (N-1)/N to the basic Chao equation
(Gotelli & Colwell 2001) as implemented by the function specpool in the
vegan R package. First-order Jackknife estimators and the Bootstrap
estimator with standard error were also used, as it is considered a more
accurate measure of species richness from datasets with a high number of
singletons and doubletons. Pearson’s Chi-squared Test was used to test
for significant differences between the observed species and
extrapolated richness with the three non-parametric methods.

**16S**

``` r
pool16s = specpool(ps_cyano@otu_table, ps_cyano@sam_data$Species, smallsample = TRUE)

DT::datatable(pool16s)%>%
    formatRound(colnames(pool16s),digits=2)
```

<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-08d3fae1b807027aa9fa" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-08d3fae1b807027aa9fa">{"x":{"filter":"none","vertical":false,"data":[["Zamia nana","Zamia pseudoparasitica"],[25,81],[105.222222222222,2274.48863636364],[68.2120875570371,2258.50242280753],[41.8888888888889,146.477272727273],[8.78831997034897,14.170650468865],[55.3055555555556,209.499471458774],[31.7957136541635,105.654014454995],[4.06788465475545,5.99717976754863],[9,44]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Species<\/th>\n      <th>chao<\/th>\n      <th>chao.se<\/th>\n      <th>jack1<\/th>\n      <th>jack1.se<\/th>\n      <th>jack2<\/th>\n      <th>boot<\/th>\n      <th>boot.se<\/th>\n      <th>n<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"targets":1,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":2,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":3,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":4,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":5,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":6,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":7,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":8,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":9,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"className":"dt-right","targets":[1,2,3,4,5,6,7,8,9]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":["options.columnDefs.0.render","options.columnDefs.1.render","options.columnDefs.2.render","options.columnDefs.3.render","options.columnDefs.4.render","options.columnDefs.5.render","options.columnDefs.6.render","options.columnDefs.7.render","options.columnDefs.8.render"],"jsHooks":[]}</script>

``` r
cyano_pseudo_acc = poolaccum(cyano.zpseudo@otu_table, permutations = 999)
plot(cyano_pseudo_acc)
```

![](Cyanobacterial_Community_files/figure-gfm/extra.plot-1.png)<!-- -->

``` r
chisq.test(cyano_pseudo_acc$chao, cyano_pseudo_acc$S, correct = TRUE)
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  cyano_pseudo_acc$chao
    ## X-squared = 3520182, df = 40918, p-value < 2.2e-16

``` r
chisq.test(cyano_pseudo_acc$jack1, cyano_pseudo_acc$S, correct = TRUE)
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  cyano_pseudo_acc$jack1
    ## X-squared = 30702, df = 40918, p-value = 1

``` r
chisq.test(cyano_pseudo_acc$boot, cyano_pseudo_acc$S, correct = TRUE)
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  cyano_pseudo_acc$boot
    ## X-squared = 20988, df = 40918, p-value = 1

``` r
cyano_nana_acc = poolaccum(cyano.znana@otu_table, permutations = 999)
plot(cyano_nana_acc)
```

![](Cyanobacterial_Community_files/figure-gfm/extra.plot-2.png)<!-- -->

``` r
chisq.test(cyano_nana_acc$chao, cyano_nana_acc$S,correct = TRUE)
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  cyano_nana_acc$chao
    ## X-squared = 65438, df = 5988, p-value < 2.2e-16

``` r
chisq.test(cyano_nana_acc$jack1, cyano_nana_acc$S, correct = TRUE)
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  cyano_nana_acc$jack1
    ## X-squared = 4823.3, df = 5988, p-value = 1

``` r
chisq.test(cyano_nana_acc$boot, cyano_nana_acc$S, correct = TRUE)
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  cyano_nana_acc$boot
    ## X-squared = 3283.2, df = 5988, p-value = 1

**rbcL-X**

``` r
poolrbclx = specpool(t(rbclx.ps@otu_table), rbclx.ps@sam_data$Species_name, smallsample = TRUE)

DT::datatable(poolrbclx)%>%
    formatRound(colnames(poolrbclx), digits=2)
```

<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-9b1ef33e7b226099bed5" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-9b1ef33e7b226099bed5">{"x":{"filter":"none","vertical":false,"data":[["Zamia nana","Zamia pseudoparasitica"],[71,141],[158.12,950.2],[39.1931095984996,368.268814862187],[110.6,236.2],[20.3823453017556,56.8534959347268],[138.688888888889,319.771428571429],[87.7130973674,178.225902611439],[10.837785770461,27.4206735336516],[10,15]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Species<\/th>\n      <th>chao<\/th>\n      <th>chao.se<\/th>\n      <th>jack1<\/th>\n      <th>jack1.se<\/th>\n      <th>jack2<\/th>\n      <th>boot<\/th>\n      <th>boot.se<\/th>\n      <th>n<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"targets":1,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":2,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":3,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":4,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":5,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":6,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":7,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":8,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":9,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"className":"dt-right","targets":[1,2,3,4,5,6,7,8,9]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":["options.columnDefs.0.render","options.columnDefs.1.render","options.columnDefs.2.render","options.columnDefs.3.render","options.columnDefs.4.render","options.columnDefs.5.render","options.columnDefs.6.render","options.columnDefs.7.render","options.columnDefs.8.render"],"jsHooks":[]}</script>

``` r
ps_pseudo_acc = poolaccum(t(ps_pseudo@otu_table), permutations = 999)
plot(ps_pseudo_acc)
```

![](Cyanobacterial_Community_files/figure-gfm/extraplot2-1.png)<!-- -->

``` r
chisq.test(ps_pseudo_acc$chao, ps_pseudo_acc$S,correct = TRUE)
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  ps_pseudo_acc$chao
    ## X-squared = 1089363, df = 11976, p-value < 2.2e-16

``` r
chisq.test(ps_pseudo_acc$jack1, ps_pseudo_acc$S, correct = TRUE)
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  ps_pseudo_acc$jack1
    ## X-squared = 99270, df = 11976, p-value < 2.2e-16

``` r
chisq.test(ps_pseudo_acc$boot, ps_pseudo_acc$S, correct = TRUE)
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  ps_pseudo_acc$boot
    ## X-squared = 65818, df = 11976, p-value < 2.2e-16

``` r
ps_nana_acc = poolaccum(t(ps_nana@otu_table), permutations = 999)
plot(ps_nana_acc)
```

![](Cyanobacterial_Community_files/figure-gfm/extraplot2-2.png)<!-- -->

``` r
chisq.test(ps_nana_acc$chao, ps_nana_acc$S, correct = TRUE)
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  ps_nana_acc$chao
    ## X-squared = 153847, df = 6986, p-value < 2.2e-16

``` r
chisq.test(ps_nana_acc$jack1, ps_nana_acc$S,  correct = TRUE)
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  ps_nana_acc$jack1
    ## X-squared = 15025, df = 6986, p-value < 2.2e-16

``` r
chisq.test(ps_nana_acc$boot, ps_nana_acc$S, correct = TRUE)
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  ps_nana_acc$boot
    ## X-squared = 12163, df = 6986, p-value < 2.2e-16

### Complete richness dataframe with the two markers

First we joined in a single the two dataframes with richness estimates
with the 16S marker and rbcL-X marker.

``` r
cyano.div.com = data.frame(rbind(richness.df.cyano, richness.rbclx.df))
cyano.div.com$Marker = as.factor(c(rep("16S",53), rep("rbcL-X",25)))
```

The assumptions of normality was tested with the Shapiro-Wilk Test for
the richness estimates from the dataset generated with the rbcL-X marker

``` r
shapiro.test(richness.rbclx.df$Observed)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  richness.rbclx.df$Observed
    ## W = 0.82983, p-value = 0.0007467

``` r
shapiro.test(richness.rbclx.df$Shannon)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  richness.rbclx.df$Shannon
    ## W = 0.89509, p-value = 0.01436

``` r
shapiro.test(richness.rbclx.df$Simpson)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  richness.rbclx.df$Simpson
    ## W = 0.81165, p-value = 0.0003579

``` r
shapiro.test(richness.rbclx.df$Pielou)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  richness.rbclx.df$Pielou
    ## W = 0.90601, p-value = 0.02485

The homogeneity of variance were tested with the Levene test

``` r
test_names <- c( "Observed", "Shannon", "Simpson", "Pielou")
for (test in test_names) {
  test_formula <- paste(test,"~","Marker")
  print(paste("Testing homogeneity of variance:", test))
  print(leveneTest(as.formula(test_formula), data=cyano.div.com))
}
```

    ## [1] "Testing homogeneity of variance: Observed"
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value    Pr(>F)    
    ## group  1  32.433 2.212e-07 ***
    ##       76                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "Testing homogeneity of variance: Shannon"
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  5.4782 0.02188 *
    ##       76                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "Testing homogeneity of variance: Simpson"
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value   Pr(>F)   
    ## group  1   6.988 0.009962 **
    ##       76                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "Testing homogeneity of variance: Pielou"
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value    Pr(>F)    
    ## group  1  17.354 8.723e-05 ***
    ##       70                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Statistical test for significance between host species and study localities

Here, the significant differences between species and markers were
assessed by fitting an Analysis of Variance Model (ANOVA). Two-way ANOVA
model use *Zamia* species host and molecular marker as predictors of
alpha diversity richness indexes in the coralloid roots.

``` r
for (test in test_names) {
  test_formula <- paste(test,"~","Host*Marker")
  print(paste("Alpha diversity:", test))
  print(summary(aov(as.formula(test_formula), data=cyano.div.com)))
}
```

    ## [1] "Alpha diversity: Observed"
    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Host         1      0       0   0.001  0.98119    
    ## Marker       1   4712    4712  73.559 1.05e-12 ***
    ## Host:Marker  1    588     588   9.173  0.00338 ** 
    ## Residuals   74   4740      64                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "Alpha diversity: Shannon"
    ##             Df Sum Sq Mean Sq F value  Pr(>F)   
    ## Host         1  0.869  0.8692   3.203 0.07760 . 
    ## Marker       1  0.786  0.7861   2.897 0.09296 . 
    ## Host:Marker  1  2.676  2.6764   9.862 0.00242 **
    ## Residuals   74 20.082  0.2714                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "Alpha diversity: Simpson"
    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## Host         1  0.073  0.0734   0.816 0.3693  
    ## Marker       1  0.054  0.0535   0.595 0.4430  
    ## Host:Marker  1  0.494  0.4936   5.486 0.0219 *
    ## Residuals   74  6.658  0.0900                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "Alpha diversity: Pielou"
    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## Host         1  0.251  0.2514   3.785 0.0559 .
    ## Marker       1  0.405  0.4051   6.097 0.0161 *
    ## Host:Marker  1  0.441  0.4406   6.633 0.0122 *
    ## Residuals   68  4.518  0.0664                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 6 observations deleted due to missingness

### Alpha diversity Figures : Alpha diversity (richness indexes) for the complete bacterial community in coralloid roots of the two species *Z. nana* and *Z. pseudoparasitica*.

``` r
cyano.div.com$Locality = factor(cyano.div.com$Locality, levels= c("El Valle, Cerro Gaital", "Donoso", "El Cope"))


Figsupp_A_cyano = 
  ggplot(cyano.div.com, aes(x=Host, y=Observed, fill=Marker)) + 
  ggtitle("Alpha Diversity", subtitle="Cyanobacteria community")+
  scale_fill_brewer(palette="Pastel2")+
  geom_boxplot(width = 0.3)+
  #scale_fill_grey()+
  theme_classic()+
  theme(text = element_text(size=20))+
  labs(x=" ", y = "Observed species")

Figsupp_A_cyano
```

![](Cyanobacterial_Community_files/figure-gfm/divplot.marker-1.png)<!-- -->

``` r
Figsupp_B_cyano = 
  ggplot(cyano.div.com, aes(x=Host, y=Shannon, fill=Marker)) + 
  #ggtitle("Alpha Diversity", subtitle="Complete bacterial community")+
  scale_fill_brewer(palette="Pastel2")+
  geom_boxplot(width = 0.3)+
  #scale_fill_grey()+
  theme_classic()+ theme(text = element_text(size=20))+
  labs(x=" ", y = "Shannon")

Figsupp_B_cyano
```

![](Cyanobacterial_Community_files/figure-gfm/divplot.marker-2.png)<!-- -->

``` r
Figsupp_C_cyano = 
  ggplot(cyano.div.com, aes(x=Host, y=Simpson, fill=Marker)) + 
  #ggtitle("Alpha Diversity", subtitle="Complete bacterial community")+
  scale_fill_brewer(palette="Pastel2")+
  geom_boxplot(width = 0.3)+
  #scale_fill_grey()+
  theme_classic()+ theme(text = element_text(size=20))+
  labs(x=" ", y = "Simpson")

Figsupp_C_cyano
```

![](Cyanobacterial_Community_files/figure-gfm/divplot.marker-3.png)<!-- -->

``` r
Figsupp_D_cyano = 
  ggplot(cyano.div.com, aes(x=Host, y=Pielou, fill=Marker)) + 
  #ggtitle("Alpha Diversity", subtitle="Complete bacterial community")+
  scale_fill_brewer(palette="Pastel2")+
  geom_boxplot(width = 0.3)+
  #scale_fill_grey()+
  theme_classic()+ theme(text = element_text(size=20))+
  labs(x=" ", y = "Pielou")

Figsupp_D_cyano
```

    ## Warning: Removed 6 rows containing non-finite values (`stat_boxplot()`).

![](Cyanobacterial_Community_files/figure-gfm/divplot.marker-4.png)<!-- -->
