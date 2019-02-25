# SahulHuman
R code for reproducing age-structured models of ancient humans entering Sahul.

The code includes a global sensitivity analysis ('ancienthumanfound_gsa_sim_func.R) that uses the following R libraries:

- iterators
- snow
- doSNOW
- foreach
- lhs
- data.table

The sensitivity analysis also requires the source code: 'matrixOperators.r' (provided in this repository), as well as the following data files:

- 'ClimateSahul_Npp.csv' (hindcasted net primary production values for northern Sahul produced by the LOVECLIM global circulation model)
- 'world2013lifetable.csv' (modern human demographic data from Bradshaw & Brook. 2014. Human population reduction is not a quick fix for environmental problems. Proceedings of the National Academy of Sciences of the USA 111: 16610–16615. doi:10.1073/pnas.1410465111)

