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
|AAA_C |  0.1561095| -0.0811611|  0.0150918|  0.2402698|
|AAC_C | -0.0874669| -0.1869301| -0.2633562|  0.0571817|
|AAG_C |  1.0181631| -0.2543222| -0.2757838|  0.2118468|
|AAT_C | -0.3327339| -0.5852982| -0.1747418| -0.0527771|
|CAA_C | -0.0730356|  0.0184412| -0.2400683|  0.2539164|
|CAC_C |  0.1277553|  0.3380247| -0.2356787|  1.5230192|

Matrix of components intensities:


```r
icS <- comp.info[[2]]
head(icS[,1:4])
```


|     comp.1|     comp.2|     comp.3|     comp.4|
|----------:|----------:|----------:|----------:|
|  0.0761317|  0.1822939| -0.4542016|  0.0451050|
| -0.1306416|  0.5153931| -0.1456649|  0.0761229|
| -0.0335879| -0.0986314| -0.0813891| -0.0086307|
|  0.1540814|  0.1043745| -0.2628011|  0.1553104|
|  0.3756479|  0.1449107| -0.3151124|  0.0182006|
|  0.0235961|  0.1611414| -0.0358629|  0.1210491|


Of note, spectra and intensities have negative values due to Z-score transformation of mutation rates in the method. Positive (negative) values in spectra and intensities indicate higher (lower) values compared to genome-wide average.
Reflection matrix of correlations between spectrum and reverse complementary spectrum of each pair of components is used to separate extracted components in strand symmetric, asymmetric and noise:


```r
plot.reflection.matrix(icM)
```

![plot of chunk unnamed-chunk-71](figure/unnamed-chunk-71-1.png)


For example, `comp. 3` is symmetic, since it has spectrum reverse complementary to itself:


```r
reflection.scatter(3,3,icM)
```

![plot of chunk unnamed-chunk-72](figure/unnamed-chunk-72-1.png)


On the other hand `comp. 1` is not symmetrical:

```r
reflection.scatter(1,1,icM)
```

![plot of chunk unnamed-chunk-73](figure/unnamed-chunk-73-1.png)

However, `comp. 1` is reflected to `comp. 2` spectrum representing a single strand-dependent process:


```r
reflection.scatter(1,2,icM)
```

![plot of chunk unnamed-chunk-74](figure/unnamed-chunk-74-1.png)

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
|comp.4 |asymmetric |5        |
|comp.5 |asymmetric |4        |
|comp.6 |symmetric  |6        |

A symmetric mutational component corresponds to a strand-independent mutational process, while a pair of two reflected components corresponds to a single strand-dependent mutational process. Annotation of mutational processes includes components that correspond to it (`comp.X` `comp.Y`):


```r
head(ref.prop[[2]])
```


|          |type               |comp. X |comp. Y |
|:---------|:------------------|:-------|:-------|
|process.1 |strand-dependent   |1       |2       |
|process.2 |strand-independent |3       |3       |
|process.3 |strand-dependent   |4       |5       |
|process.4 |strand-independent |6       |6       |
|process.5 |strand-dependent   |10      |7       |
|process.6 |strand-independent |8       |8       |

Spectrum visualization of symmetric `comp. 3` with equal rates of complementary mutations: 


```r
draw.signature(icM[,3])
```

![plot of chunk unnamed-chunk-80](figure/unnamed-chunk-80-1.png)

Spectrum visualization of asymmetric `comp. 1` with imbalanced rates of complementary mutations: 


```r
draw.signature(icM[,1])
```

![plot of chunk unnamed-chunk-81](figure/unnamed-chunk-81-1.png)

Routine `plot.intensities` enables exploration of spatial variation in component intensity along the genome. Intensity of maternal process, estimated as sum of intensities of `comp. 4` and `comp. 5`, along left arm of chromosome 8 is shown below:

```r
plot.intensities(icS[,4]+icS[,5], info, chr=8, start=0,end=4e+7,  span.wind=30)
```

![plot of chunk unnamed-chunk-82](figure/unnamed-chunk-82-1.png)

Finally, robustness of components can be evaluated using statistical bootstrap of genomic regions or reflection correlation. Rounine  `spectra.bootstrap` estimates Spearman correlations between a component in original and `n.boot` bootstrapped inferences:


```r
boot <- spectra.bootstrap(icM,rate,n.boot=500,n.cores=20)
```

Visualization of summarized statistics:

```r
visualize.bootstrap(boot,icM=icM)
```

![plot of chunk unnamed-chunk-84](figure/unnamed-chunk-84-1.png)


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

![plot of chunk unnamed-chunk-87](figure/unnamed-chunk-87-1.png)




