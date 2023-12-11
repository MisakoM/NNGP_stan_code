# Stan code for NNGP

This stan code is for estimating the effect of environmental factors on local extinction of Red list vascular plants following phylogenetic Gaussian Processes (GPs) implemented by nearest neighbour GP (NNGP) approximation using cmdstanr in R.
Full methods are described in XXXX(2023).


# Data list
 ・splist : species list data for calculating the phylogenetic distance using V.Phylomaker (Jin & Quian 2019 Ecography)
 
 ・spdata : data format for spdata, which includes the headers and their discription. The distribution data for Red List vascular plants are not allowed to be published for conservation reasons. Therefore, this code does not include the actual data, but only provides the format with column names retained. 

# Reference
Jin, Y., & Qian, H. (2019). V. PhyloMaker: an R package that can generate very large phylogenies for vascular plants. Ecography, 42(8), 1353-1359. https://doi.org/10.1111/ecog.04434

XXXX(2023)
