# Deconvolution of cellular heterogeneity in brain transcriptomes (2018)

This repository contains the R scripts necessary to replicate the figures in this manuscript.

## IHC image analysis

IHCimageAnalysis contains the R scripts neuroProp.R, astroProp.R, microgliaProp.R, endoProp.R and oligoProp.R which can be used to process the IHC images and output proportion estimates. This image analysis is primarily performed using the EBImage package.

The IHC images can be found at [neurons](https://www.dropbox.com/s/k1lq99b8hmpgy7o/NeuN.zip?dl=0), [astrocytes](https://www.dropbox.com/s/z67ft05x1xnzofv/GFAP.zip?dl=0), [microglia](https://www.dropbox.com/s/4ihwqhmh3gtap1r/iba1.zip?dl=0), [oligodendrocytes](https://www.dropbox.com/s/yvf0wb22zy3lwam/Olig2.zip?dl=0) and [endothelial cells](https://www.dropbox.com/s/lyykkg62u0j4984/PECAM.zip?dl=0). Careful, these zip files are over 10GB each.



## Cell type deconvolution analysis

CellTypeDeconvAnalysis contains the R script DeconvAnalysis.Rmd which can be used to generate all of the figures from the manuscript. This folder also contains all of the necassary data.


![alt text](https://github.com/ellispatrick/CortexCellDeconv/raw/master/CellTypeDeconvAnalysis/figures/Figure2b_proportionZhang-1.png "Figure 2b")



