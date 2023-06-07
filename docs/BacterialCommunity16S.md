Analyses of bacterial diversity using 16S
================

### Zamia Coralloid root bacterial community

### Contributors

Juan Carlos Villarreal A.

Adriel M. Sierra

To study the overall coralloid root bacterial community of Zamia
species, we used the processed phyloseq object with the whole bacterial
community.

The dataset include a total of fourty four coralloid roots of Z.
pseudoparasitica from two localities, and ten of Zamia nana from one
locality

``` r
knitr::opts_knit$set(root.dir = "/Users/adrielsierra/Documents/GitHub/cycads_rbcLX_roots")
root16ps = readRDS("./Data/Zamia_root_bacteria16s_ps_mai23.rds")
root16ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1561 taxa and 53 samples ]
    ## sample_data() Sample Data:       [ 53 samples by 5 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1561 taxa by 7 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 1561 tips and 1559 internal nodes ]

### Compute and visualize alpha diversity

Alpha diversity was assessed using four indexes: i) the Chao index which
shows the observed richness, ii) the exponential of Shannon which
weights the ASVs by their frequency, iii) the Simpson’s diversity index
which considers the presence and relative abundances of ASVs present in
a sample, and iv) the Pielou index which measures diversity along with
species richness.

``` r
richness.df = estimate_richness(root16ps, split = TRUE, measures = c("Observed","Shannon", "Simpson"))

richness.df$Pielou = (richness.df$Shannon)/log(specnumber(otu_table(root16ps)))
```

Now, we add extra column with read count by sample, the Zamia host
species and locality respective to each sample

``` r
richness.df$readcount = sample_sums(root16ps)

richness.df$Host = as.factor(root16ps@sam_data$Species)

richness.df$Locality = as.factor(root16ps@sam_data$Locality)
```

### Alpha diversity table : Alpha diversity (richness indexes) for the complete bacterial community in coralloid roots of the two species Z. nana and Z. pseudoparasitica.

``` r
 DT::datatable(richness.df) %>%
formatRound(columns = c("Shannon", "Simpson", "Pielou"), digits=4)
```

<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-ffa3c8687ef01085dc24" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-ffa3c8687ef01085dc24">{"x":{"filter":"none","vertical":false,"data":[["X10S","X11S","X12S","X13S","X14S","X15S","X16S","X17S","X18S","X19S","X1S","X20S","X22S","X23S","X25S","X26S","X27S","X28S","X2S","X3S","X4S","X5S","X6S","X7S","X8S","X9S","ZAP02","ZAP10","ZAP12","ZAP14","ZAP18","ZAP25","ZAP26","ZAP27","ZAP28","ZAP29","ZAP30","ZAP31","ZAP32","ZAP33","ZAP34","ZAP36","ZAP37","ZAP39","ZAP40","ZAP41","ZAP43","ZAP44","ZAP45","ZAP46","ZAP47","ZAP51","ZAP55"],[62,39,27,44,32,57,36,28,47,17,70,71,60,53,42,38,58,46,16,35,36,0,33,36,66,25,57,80,48,7,32,44,5,116,15,34,32,87,24,43,70,40,118,39,20,75,6,27,35,1,3,22,38],[3.37213181445108,2.50137248188099,2.65601028354699,3.24123880312947,2.54475609880836,3.16488774039833,2.84840615278171,2.39062203117206,3.13690980698085,2.38451885526106,3.81444253831758,3.53364808289095,3.65048899100866,3.34860882705927,3.05967362223999,2.77304190187043,3.74287810507723,3.04160730963147,1.5784451458262,2.55063252961205,1.75862894480902,0,3.01993074988582,2.24036372186678,2.97668971077958,2.74353305570356,3.52328947909838,3.65411918972826,3.18116835430628,1.53412546271815,1.01627950508218,2.63285524370825,1.25674877853834,4.09693574634844,2.20250908045689,3.16972458081712,2.98361412288347,3.69052481406956,2.08122792348536,1.06947575274652,2.09010176395178,1.65444043279578,3.61607794915107,0.94751453796996,1.03707736949188,3.34189012784467,0.804292923720309,1.03545522512246,2.05394677646532,0,0.0213274080309622,1.80192732397652,2.9259365927408],[0.941467183602364,0.84181695854838,0.900603704104191,0.942103472468879,0.866497066627394,0.914213905061429,0.916154361775748,0.850716715308751,0.923049905439316,0.875207649476213,0.967559651736702,0.95129717885771,0.964838575162363,0.945238991814304,0.923482255218224,0.870314974253723,0.972536616433153,0.911113037611698,0.693454449421305,0.873936805750878,0.73323034005336,1,0.941976086820012,0.829183578738408,0.894522813974589,0.922069607695016,0.960727840535533,0.950727903386563,0.931456086894586,0.70897924951979,0.326392900518485,0.817548430460068,0.648279083861519,0.972308304157078,0.857255924137388,0.946089495163019,0.920533323910211,0.956433890020833,0.741700028190788,0.335360226853234,0.615476327291631,0.530265893386622,0.918554227095699,0.323210054564159,0.409285552852205,0.918986974190606,0.439277321971743,0.547069956562814,0.736451953387529,0,0.00591137986938295,0.750101553425328,0.916688883207982],[0.817063729901839,0.68277068150977,0.805868248802912,0.856521241450946,0.734261400804554,0.782796836037511,0.794862871300714,0.7174295646837,0.814750532054298,0.841630532435381,0.897833779019981,0.828973365305495,0.891593009706894,0.843416277687086,0.818604621252006,0.762330232221636,0.921790576879331,0.794435151962955,0.56930374604968,0.717407265195728,0.490754751129262,-0,0.863699190348446,0.625185400256876,0.710485682965808,0.852326465813851,0.871443186270909,0.833887946005195,0.821752305122543,0.788384532278082,0.293236280427831,0.695751402125721,0.780861920070993,0.861861368279816,0.813319147318774,0.898865716846165,0.860889059802051,0.826377237510516,0.654874975247822,0.284344283699275,0.491962835044526,0.448494035485681,0.757978829832225,0.25863207159922,0.346184930691933,0.774035744850595,0.448884427588276,0.314170654440727,0.577706244492102,null,0.019413043392058,0.582951334757237,0.804362141338497],[9682,4942,8160,8095,8109,11383,9972,6818,8570,2706,9540,10742,10657,11947,9130,6752,8218,8924,6114,9775,19154,0,12683,29539,22676,5067,3432,5204,3991,333,8748,4190,398,10751,1133,2088,2074,8360,3517,16863,15527,10591,14224,24880,8956,8445,3065,17139,13758,14398,5399,22425,5073],["Zamia pseudoparasitica","Zamia nana","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia nana","Zamia nana","Zamia nana","Zamia pseudoparasitica","Zamia nana","Zamia nana","Zamia nana","Zamia nana","Zamia nana","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica","Zamia pseudoparasitica"],["Donoso","El Valle, Cerro Gaital","Donoso","Donoso","Donoso","Donoso","Donoso","El Valle, Cerro Gaital","El Valle, Cerro Gaital","El Valle, Cerro Gaital","Donoso","El Valle, Cerro Gaital","El Valle, Cerro Gaital","El Valle, Cerro Gaital","El Valle, Cerro Gaital","El Valle, Cerro Gaital","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","Donoso","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope","El Cope"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Observed<\/th>\n      <th>Shannon<\/th>\n      <th>Simpson<\/th>\n      <th>Pielou<\/th>\n      <th>readcount<\/th>\n      <th>Host<\/th>\n      <th>Locality<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"targets":2,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 4, 3, \",\", \".\", null);\n  }"},{"targets":3,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 4, 3, \",\", \".\", null);\n  }"},{"targets":4,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 4, 3, \",\", \".\", null);\n  }"},{"className":"dt-right","targets":[1,2,3,4,5]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":["options.columnDefs.0.render","options.columnDefs.1.render","options.columnDefs.2.render"],"jsHooks":[]}</script>

The assumptions of normality was tested with the Shapiro-Wilk Test

``` r
shapiro.test(richness.df$Observed)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  richness.df$Observed
    ## W = 0.94207, p-value = 0.01244

``` r
shapiro.test(richness.df$Shannon)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  richness.df$Shannon
    ## W = 0.93602, p-value = 0.007031

``` r
shapiro.test(richness.df$Simpson)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  richness.df$Simpson
    ## W = 0.74905, p-value = 3.715e-08

``` r
shapiro.test(richness.df$Pielou)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  richness.df$Pielou
    ## W = 0.81468, p-value = 1.338e-06

The homogeneity of variance were tested with the Levene test

``` r
test_names <- c( "Observed", "Shannon", "Simpson", "Pielou")
for (test in test_names) {
  test_formula <- paste(test,"~","Host")
  print(paste("Testing homogeneity of variance:", test))
  print(leveneTest(as.formula(test_formula), data=richness.df))
}
```

    ## [1] "Testing homogeneity of variance: Observed"
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  1.2698 0.2651
    ##       51               
    ## [1] "Testing homogeneity of variance: Shannon"
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  5.0997 0.02824 *
    ##       51                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "Testing homogeneity of variance: Simpson"
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  3.4914 0.06744 .
    ##       51                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "Testing homogeneity of variance: Pielou"
    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  3.6555 0.06162 .
    ##       50                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Statistical test for significance between host species and study localities

Here, the significant differences between species and markers were
assessed using the Kruskal-Wallis nonparametric test. Post hoc pairwise
comparisons were performed with the exact sum of Wilcoxon ranks.

``` r
for (test in test_names) {
  test_formula <- paste(test,"~","Host")
  print(paste("Alpha diversity:", test))
  print(kruskal.test(as.formula(test_formula), data=richness.df))
}
```

    ## [1] "Alpha diversity: Observed"
    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Observed by Host
    ## Kruskal-Wallis chi-squared = 0.81082, df = 1, p-value = 0.3679
    ## 
    ## [1] "Alpha diversity: Shannon"
    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Shannon by Host
    ## Kruskal-Wallis chi-squared = 2.1572, df = 1, p-value = 0.1419
    ## 
    ## [1] "Alpha diversity: Simpson"
    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Simpson by Host
    ## Kruskal-Wallis chi-squared = 2.0881, df = 1, p-value = 0.1485
    ## 
    ## [1] "Alpha diversity: Pielou"
    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Pielou by Host
    ## Kruskal-Wallis chi-squared = 2.434, df = 1, p-value = 0.1187

``` r
for (test in test_names) {
  test_formula <- paste(test,"~","Locality")
  print(paste("Alpha diversity:", test))
  print(kruskal.test(as.formula(test_formula), data=richness.df))
}
```

    ## [1] "Alpha diversity: Observed"
    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Observed by Locality
    ## Kruskal-Wallis chi-squared = 0.91613, df = 2, p-value = 0.6325
    ## 
    ## [1] "Alpha diversity: Shannon"
    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Shannon by Locality
    ## Kruskal-Wallis chi-squared = 4.3259, df = 2, p-value = 0.115
    ## 
    ## [1] "Alpha diversity: Simpson"
    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Simpson by Locality
    ## Kruskal-Wallis chi-squared = 7.5496, df = 2, p-value = 0.02294
    ## 
    ## [1] "Alpha diversity: Pielou"
    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Pielou by Locality
    ## Kruskal-Wallis chi-squared = 4.1152, df = 2, p-value = 0.1278

### Alpha diversity Figures : Alpha diversity (richness indexes) for the complete bacterial community in coralloid roots of the two species Z. nana and Z. pseudoparasitica.

``` r
Figsupp_A = 
  ggplot(richness.df, aes(x=Host, y=Observed, fill=Locality)) + 
  ggtitle("Alpha Diversity", subtitle="Complete bacterial community")+
  scale_fill_brewer(palette="Pastel2")+
  geom_boxplot(width = 0.3)+
  #scale_fill_grey()+
  theme_classic()+
  theme(text = element_text(size=12))+
  labs(x=" ", y = "Observed species")

Figsupp_A
```

![](BacterialCommunity16S_files/figure-gfm/divplot-1.png)<!-- -->

``` r
Figsupp_B = 
  ggplot(richness.df, aes(x=Host, y=Shannon, fill=Locality)) + 
  #ggtitle("Alpha Diversity", subtitle="Complete bacterial community")+
  scale_fill_brewer(palette="Pastel2")+
  geom_boxplot(width = 0.3)+
  #scale_fill_grey()+
  theme_classic()+ theme(text = element_text(size=12))+
  labs(x=" ", y = "Shannon")

Figsupp_B
```

![](BacterialCommunity16S_files/figure-gfm/divplot-2.png)<!-- -->

``` r
Figsupp_C = 
  ggplot(richness.df, aes(x=Host, y=Simpson, fill=Locality)) + 
  #ggtitle("Alpha Diversity", subtitle="Complete bacterial community")+
  scale_fill_brewer(palette="Pastel2")+
  geom_boxplot(width = 0.3)+
  #scale_fill_grey()+
  theme_classic()+ theme(text = element_text(size=12))+
  labs(x=" ", y = "Simpson")
  
Figsupp_C
```

![](BacterialCommunity16S_files/figure-gfm/divplot-3.png)<!-- -->

``` r
Figsupp_D = 
  ggplot(richness.df, aes(x=Host, y=Pielou, fill=Locality)) + 
  #ggtitle("Alpha Diversity", subtitle="Complete bacterial community")+
  scale_fill_brewer(palette="Pastel2")+
  geom_boxplot(width = 0.3)+
  #scale_fill_grey()+
  theme_classic()+ theme(text = element_text(size=12))+
  labs(x=" ", y = "Pielou")

Figsupp_D
```

    ## Warning: Removed 1 rows containing non-finite values (`stat_boxplot()`).

![](BacterialCommunity16S_files/figure-gfm/divplot-4.png)<!-- -->
