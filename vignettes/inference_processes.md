---
title: "Inference of germline mutational processes"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Inference of germline mutational processes}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



The vignette describes inference of spatially-varying mutational processes. The guideline starts with matrix of regional mutation rates of mutation types in genomic regions, infers mutational components, estimates their biological relevance and groups in strand-independend/strand-dependent mutational processes. Finally, the guideline provides quality estimates of reconstructed processes using bootstrap and reflection correlations. 

## Loading the data


```r
library(spacemut)
```

Matrix of regional mutation rates is loaded from the local web server:


```r
info <- read.table("http://pklab.med.harvard.edu/ruslan/spacemut/TOPMed.30kb.txt",header=TRUE)
dim(info)
## [1] 87496   195
```


```r
head(info[,1:6])
```


| chr|  start|    end|     AAA_C|     AAC_C|     AAG_C|
|---:|------:|------:|---------:|---------:|---------:|
|   1| 810000| 840000| 0.0131313| 0.0148515| 0.0196078|
|   1| 840000| 870000| 0.0171359| 0.0096852| 0.0055659|
|   1| 870000| 900000| 0.0183968| 0.0081081| 0.0107527|
|   1| 900000| 930000| 0.0305882| 0.0191571| 0.0225000|
|   1| 930000| 960000| 0.0284698| 0.0124481| 0.0469613|
|   1| 960000| 990000| 0.0198413| 0.0186047| 0.0272109|

Retain only columns of regional mutation rates:

```r
rate <- info[,-c(1:3)]
```

## Inference of mutational components and processes

Function `extract.comp` decomposes observed regional mutation rates in a small number of mutational components using independent component analysis. Below we estimate 13 mutational components:


```r
comp.info <- extract.comp(rate,13)
```

Output object `comp.info` consists of matrices of components spectra  and intensities. Matrix of components spectra:


```r
icM <- comp.info[[1]]
head(icM[,1:4])
```


|      |     comp.1|     comp.2|     comp.3|     comp.4|
|:-----|----------:|----------:|----------:|----------:|
|AAA_C | -0.1518822|  0.1471248|  0.0935838|  1.5667148|
|AAC_C | -0.2755600| -0.1149889| -0.1862629|  1.7804012|
|AAG_C | -0.3187923|  1.0127495| -0.1973558|  1.9464271|
|AAT_C | -0.5558670| -0.3287671| -0.0797052|  0.9466719|
|CAA_C | -0.0351220| -0.0821468| -0.1965866|  1.2996496|
|CAC_C |  0.3389197|  0.1762089| -0.2458765| -0.2286360|

Matrix of components intensities:


```r
icS <- comp.info[[2]]
head(icS[,1:4])
```


|     comp.1|     comp.2|     comp.3|     comp.4|
|----------:|----------:|----------:|----------:|
|  0.1796166|  0.0721438| -0.4479047|  0.1779447|
|  0.5102889| -0.1320399| -0.1451787|  0.1336934|
| -0.0918585| -0.0342882| -0.0810809| -0.0979801|
|  0.0956659|  0.1482672| -0.2599748|  0.1147958|
|  0.1351631|  0.3616662| -0.3003602|  0.3959587|
|  0.1522517|  0.0206364| -0.0299601|  0.1796010|


Of note, spectra and intensities have negative values due to Z-score transformation of mutation rates in the method. Positive (negative) values in spectra and intensities indicate higher (lower) values compared to genome-wide average.
Reflection matrix of correlations between spectrum and reverse complementary spectrum of each pair of components is used to separate extracted components in strand symmetric, asymmetric and noise:


```r
plot.reflection.matrix(icM)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)


For example, `comp. 3` is symmetic, since it has spectrum reverse complementary to itself:


```r
reflection.scatter(3,3,icM)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)


On the other hand `comp. 1` is not symmetrical:

```r
reflection.scatter(1,1,icM)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)

However, `comp. 1` is reflected to `comp. 2` spectrum representing a single strand-dependent process:


```r
reflection.scatter(1,2,icM)
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png)

Function `reflection.test` formally classifies components in symmetric/asymmetric, noise and combine them in mutaitonal processes.


```r
ref.prop <- reflection.test(icM)
```

Object `ref.prop` contains classification of components, including symmetric/asymmteric `type` and reflected component `ref.comp` for each component:

```r
head(ref.prop[[1]])
```


|       |type       |ref.comp |
|:------|:----------|:--------|
|comp.1 |asymmetric |2        |
|comp.2 |asymmetric |1        |
|comp.3 |symmetric  |3        |
|comp.4 |symmetric  |4        |
|comp.5 |symmetric  |5        |
|comp.6 |asymmetric |7        |

A symmetric mutational component corresponds to a strand-independent mutational process, while a pair of two reflected components corresponds to a single strand-dependent mutational process. Annotation of mutational processes includes components that correspond to it (`comp.X` `comp.Y`):


```r
head(ref.prop[[2]])
```


|          |type               |comp. X |comp. Y |
|:---------|:------------------|:-------|:-------|
|process.1 |strand-dependent   |1       |2       |
|process.2 |strand-independent |3       |3       |
|process.3 |strand-independent |4       |4       |
|process.4 |strand-independent |5       |5       |
|process.5 |strand-dependent   |6       |7       |
|process.6 |strand-independent |8       |8       |

Spectrum visualization of symmetric `comp. 3` with equal rates of complementary mutations: 


```r
draw.signature(icM[,3])
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20-1.png)

Spectrum visualization of asymmetric `comp. 1` with imbalanced rates of complementary mutations: 


```r
draw.signature(icM[,1])
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21-1.png)

Routine `plot.intensities` enables exploration of spatial variation in component intensity along the genome. Intensity of maternal process, estimated as sum of intensities of `comp. 4` and `comp. 5`, along left arm of chromosome 8 is shown below:

```r
plot.intensities(icS[,4]+icS[,5], info, chr=8, start=0,end=4e+7,  span.wind=30)
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22-1.png)

Finally, robustness of components can be evaluated using statistical bootstrap of genomic regions or reflection correlation. Rounine  `spectra.bootstrap` estimates Spearman correlations between a component in original and `n.boot` bootstrapped inferences:


```r
boot <- spectra.bootstrap(icM,rate,n.boot=500,n.cores=20)
```

Visualization of summarized statistics:

```r
visualize.bootstrap(boot,icM=icM)
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24-1.png)


## Estimate number of components to extract
Number 13 of initially extracted components may seems rather arbitrary. Reflection property suggests a natural criterion to evaluate number of components to extract that maximumizes number of components having reflection property. More specifically, routine `select.ncomponents` extracts a range of components from `n.min` and `n.max` and evaluates optimal choice:


```r
stat.extract <- select.ncomponents(rate,n.min=2,n.max=35,cutoff=0.8,n.cores=20)

head(stat.extract)
```


|  n| comp| proc|
|--:|----:|----:|
|  2|    2|    2|
|  3|    3|    3|
|  4|    4|    4|

Number of components and processes having reflection property depending on number of extracted components:

```r
show.optimal.comp(stat.extract)
```

![plot of chunk unnamed-chunk-27](figure/unnamed-chunk-27-1.png)




