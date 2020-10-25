# SahulHuman
R code for reproducing age-structured models of ancient humans entering Sahul.

The repository's code and associated data are to reproduce the models presented in the following companion papers:

- Bradshaw, CJA, S Ulm, AN Williams, MI Bird, RG Roberts, Z Jacobs, F Laviano, LS Weyrich, T Friedrich, K Norman, F Saltré. 2019. Minimum founding populations for the first peopling of Sahul. <i>Nature Ecology and Evolution</i> 3: 1057-1063 http://doi.org/10.1038/s41559-019-0902-6

- Bird, MI, SA Condie, S O’Connor, D O’Grady, C Reepmeyer, S Ulm, M Zega, F Saltré, CJA Bradshaw. 2019. Early human settlement of Sahul was not an accident. <i>Scientific Reports</i> 9: 8220 http://doi.org/10.1038/s41598-019-42946-9

The first R code file ('humanpopmodelgithub.R) is a stochastic projection of the human population entering Sahul required to estimate minimum viable population size. It requires the following R libraries:

- boot
- tcltk
- sp
- rgdal
- raster

The code also requires the source code: 'matrixOperators.r' (provided in this repository), as well as the following data files:

- 'world2013lifetable.csv' (modern human demographic data from Bradshaw & Brook. 2014. Human population reduction is not a quick fix for environmental problems. <i>Proceedings of the National Academy of Sciences of the USA</i> 111: 16610–16615. doi:10.1073/pnas.1410465111)
- 'ClimateSahul_Npp.csv' (hindcasted net primary production values for northern Sahul produced by the LOVECLIM global circulation model)


The second R code file ('human arrival population model_github.R') is a variant of the above, only this time applied to the island-hopping scenario presented in the Bird et al. analysis.

This code requires the following R libraries:

- boot
- tcltk
- plotly
- sp
- rgdal
- raster

and the same data files as above ('world2013lifetable.csv' & 'ClimateSahul_Npp.csv').


The repository also includes a global sensitivity analysis ('ancienthumanfound_gsa_sim_func_v2.R) derived from the paper: Prowse, TAA, CJA Bradshaw, et al. 2016. An efficient protocol for the global sensitivity analysis of stochastic ecological models. <i>Ecosphere</i> 7: e01238 http://dx.doi.org/10.1002/ecs2.1238. This requires the following R libraries:

- iterators
- snow
- doSNOW
- foreach
- lhs
- data.table

The sensitivity analysis also relies on the source code: 'matrixOperators.r' (provided in this repository), as well as the following data files:

- 'ClimateSahul_Npp.csv' (hindcasted net primary production values for northern Sahul produced by the LOVECLIM global circulation model)
- 'world2013lifetable.csv' (modern human demographic data from Bradshaw & Brook. 2014. Human population reduction is not a quick fix for environmental problems. Proceedings of the National Academy of Sciences of the USA 111: 16610–16615 http://dx.doi.org/10.1073/pnas.1410465111)

