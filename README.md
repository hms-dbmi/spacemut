# spacemut
Spatial analysis of germline mutation patterns.

This page containts data and code from our [manuscript](https://www.biorxiv.org/content/10.1101/2020.01.10.893024v1.abstract). We developed a computational approach that employs co-variation of mutation frequencies along the genome to extract mutational processes operating in the human germline. Toward that goal, we implementated volume-regularized NMF (vrnmf), a powerful technique to reconstruct non-negative sources from their observed mixtures. Analysis of large-scale whole genome sequencing TOPMed dataset reveals 4 strand-independent and 5 strand-dependent processes. For 7 of which we provided biological associations or interpretations.  

## Data

Below we provide description of data files for TOPMed and gnomAD datset. 

### TOPMed dataset

[TOPMed_10kb.txt](http://pklab.med.harvard.edu/ruslan/spacemut/tracks_update/TOPMed_10kb.txt) contains estimated mutation frequency of 192 mutation type across 263,870 non-intersecting 10 kb genomic window. Each row of a matrix contains information on a window: chromosome, start and end positions and 192 mutation frequencies.

#### Standard deviation normalization

Before inference of mutational components each mutation rate was normalized to its standard deviation (sd) across genomic windows. Below are presented tracks of mutational components for sd-normalized mutation rates (these tracks are used in the paper):

[TOPMed_10kb_spectra_sd.txt](http://pklab.med.harvard.edu/ruslan/spacemut/tracks_update/TOPMed_10kb_spectra_sdnorm.txt) contains spectra of mutational components. Each column contains a vector of relative probabilities to generate 192 sd-normalized mutation types for a mutational component. 

[TOPMed_10kb_intensities_sd.txt](http://pklab.med.harvard.edu/ruslan/spacemut/tracks_update/TOPMed_10kb_intensities_sdnorm.txt) contains intensities of mutational components.  Each row represents a genomic window and contains information on its chromosome, start and end positions, as well as window-average intensities of 14 sd-normalized mutational components. 

#### Original scales

Below are spectra and intensities of mutational components converted in the original scales of mutation rates (by multiplying to mutation type-specific standard deviation):

[TOPMed_10kb_spectra.txt](http://pklab.med.harvard.edu/ruslan/spacemut/tracks_update/TOPMed_10kb_spectra.txt) contains spectra of mutational components. Each column contains a vector of relative probabilities to generate 192 mutation types for a mutational component. 

[TOPMed_10kb_intensities.txt](http://pklab.med.harvard.edu/ruslan/spacemut/tracks_update/TOPMed_10kb_intensities.txt) contains intensities of mutational components.  Each row represents a genomic window and contains information on its chromosome, start and end positions, as well as window-average intensities of 14 mutational components. 


### gnomAD dataset

[gnomAD_100kb.txt](http://pklab.med.harvard.edu/ruslan/spacemut/tracks_update/gnomAD_100kb.txt) contains estimated mutation frequency of 192 mutation type across 26,625 non-intersecting 100 kb genomic window. 

#### Standard deviation normalization

[gnomAD_100kb_spectra_sd.txt](http://pklab.med.harvard.edu/ruslan/spacemut/tracks_update/gnomAD_100kb_spectra_sdnorm.txt) contains spectra of 12 sd-normalized mutational components. 

[gnomAD_100kb_intensities_sd.txt](http://pklab.med.harvard.edu/ruslan/spacemut/tracks_update/gnomAD_100kb_intensities_sdnorm.txt) contains intensities of 12 sd-normalized mutational components. 

#### Original scales

[gnomAD_100kb_spectra.txt](http://pklab.med.harvard.edu/ruslan/spacemut/tracks_update/gnomAD_100kb_spectra.txt) contains spectra of 12 mutational components. 

[gnomAD_100kb_intensities.txt](http://pklab.med.harvard.edu/ruslan/spacemut/tracks_update/gnomAD_100kb_intensities.txt) contains intensities of 12 mutational components. 

(gnomAD dataset reveals only 12 components compared to 14 components of TOPMed dataset due to more than order of magnitude smaller sample size) 

## Vrnmf inference of mutational processes from TOPMed dataset

Inference of mutational components from TOPMed dataset can be reproduced with a few simple steps using [_vrnmf_ R package](https://github.com/kharchenkolab/vrnmf). 

Following installation of _vrnmf_ and _spacemut_ R packages (below), one can proceed to a guide through [application of vrnmf to inference of mutational processes in germline](https://github.com/hms-dbmi/spacemut/blob/master/vignettes/vrnmf_germline.md).

### Installation instructions

```{r setup}
install.packages("devtools")
devtools::install_github("kharchenkolab/vrnmf")
library(vrnmf)
```
A detailed guide through _vrnmf_ can be found [here](https://github.com/kharchenkolab/vrnmf/blob/master/vignettes/volume_regularized_NMF.md).

We also recommend to install _spacemut_ R package that provides a set of routines to visualize and interpret inferred mutational components.

```{r setup}
install.packages("devtools")
devtools::install_github("hms-dbmi/spacemut")
library(spacemut)
```

## Available code

Code of simulations of mutational processes is available here:
http://pklab.med.harvard.edu/ruslan/spacemut/simulations_topmed.R

Other code is available upon the request.
