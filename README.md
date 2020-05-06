# spacemut
Spatial analysis of genome mutation patterns

Instructions on inference of mutational processes from spatial patterns are described in https://github.com/hms-dbmi/spacemut/blob/master/vignettes/inference_processes.md


## Installation instructions

```{r setup}
install.packages("devtools")
devtools::install_github("hms-dbmi/spacemut")
library(spacemut)
```

## gnomAD data

[Download](http://pklab.med.harvard.edu/ruslan/spacemut/gnomAD_100kb_rates.txt) matrix of mutation rates per 100kb.

[Download](http://pklab.med.harvard.edu/ruslan/spacemut/gnomAD_100kb_spectra.txt) spectra of 12 mutational components.

[Download](http://pklab.med.harvard.edu/ruslan/spacemut/gnomAD_100kb_intensities.txt) genome intensities of 12 mutational components.

[Download](http://pklab.med.harvard.edu/ruslan/spacemut/gnomAD_100kb_processes_annotation.txt) annotation of 8 processes grouping 12 mutational components.
