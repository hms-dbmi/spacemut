# vrnmf_germline

## Loading the data
The inference requires only a matrix of window-specific mutation rates. Below we use a matrix calculated based on TOPMed dataset:

```r
library(vrnmf)
library(spacemut)

rate.info <- read.table("http://pklab.med.harvard.edu/ruslan/spacemut/tracks/TOPMed_10kb.txt",header=TRUE)
```


```r
head(rate.info[,1:6])
```


| chr|  start|    end|     AAA_C|     AAC_C|     AAG_C|
|---:|------:|------:|---------:|---------:|---------:|
|   1| 770000| 870000| 0.0056657| 0.0357143| 0.0074014|
|   1| 780000| 880000| 0.0082418| 0.0070821| 0.0168406|
|   1| 790000| 890000| 0.0078431| 0.0091439| 0.0314421|
|   1| 800000| 900000| 0.0165351| 0.0293756| 0.0296436|
|   1| 810000| 910000| 0.0082497| 0.0000000| 0.0207556|
|   1| 820000| 920000| 0.0093458| 0.0198020| 0.0185262|

Retain only columns of regional mutation rates:

```r
rate <- rate.info[,-c(1:3)]
```

## Vrnmf inference of spectra of mutational components

First, matrix is preprocessed and required statistics of covariance matrix are estimated using `vol_preprocess` routine:

```r
vol.info <- vol_preprocess(rate)
```


Volume-regularized NMF is applied to the co-occurence mutation matrix stored in the object `vol.info` to infer `n.comp = 14` components using volume weight `wvol = 7.9e-3` (see the manuscript for the motivation behind the choice of parameters). Alternating optimization was executed for `n.iter = 3e+3` iterations with `vol.iter` and `c.iter` updates of matrices inside each iteration:


```r
vr <- volnmf_main(vol.info, n.comp = 14, wvol = 7.9e-3,
                  n.iter = 3e+3, vol.iter = 1e+1, c.iter = 1e+1, 
                  verbose = FALSE)
```

Non-normalized spectra of mutational components are available in the matrix `vr$C`. To reflect the results in the paper, mutational spectra are then renormalized to the genome-average vector of mutation rates:


```r
S <- vr$C*vol.info$col.factors/colMeans(rate)
rownames(S) <- colnames(rate)
colnames(S) <- paste("comp.",1:ncol(S))
```

Note that order of components (columns of `S`) is arbitrary and may not be consistent with that used in the manuscript. 
Spectra can now be visualized:

```r
draw.signature(S[,2])
```

![plot of chunk unnamed-chunk-56](figure/unnamed-chunk-56-1.png)

Additionally, we can explore reflection properties of components:

```r
plot.reflection.matrix(S)
```

![plot of chunk unnamed-chunk-57](figure/unnamed-chunk-57-1.png)

And visualize reflections of individual components:

```r
reflection.scatter(3,3,S)
```

![plot of chunk unnamed-chunk-58](figure/unnamed-chunk-58-1.png)

## Estimation of intensities of mutational components

Given known spectra of mutational components, intensities can then be obtained using non-negative least squares. A fast per-window estimation is implemented in a routine `infer_intensities`. However, we will additionally model offset spectrum, that would reflect all unaccounted mutational forces, using a routine `factor_intensities` (it takes a few to ten minutes to optimize):


```r
intensities.info <- factor_intensities(as.matrix(vr$C), t(as.matrix(vol.info$X.process)), 
                                  fit.nmf = FALSE, fit.factor = TRUE, 
                                  n.iter = 1e+2, verbose = TRUE) 
```

![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-1.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-2.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-3.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-4.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-5.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-6.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-7.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-8.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-9.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-10.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-11.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-12.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-13.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-14.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-15.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-16.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-17.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-18.png)![plot of chunk unnamed-chunk-59](figure/unnamed-chunk-59-19.png)

```r
str(intensities.info)
```


```r
intensities <- cbind(rate.info[,1:3],t(intensities.info$intensities))
intensities <-t(intensities.info$intensities)
```


```r
plot.intensities(intensities[,2], rate.info, chr=16, start=0,end=4e+7,  span.wind=100)
```

![plot of chunk unnamed-chunk-61](figure/unnamed-chunk-61-1.png)

