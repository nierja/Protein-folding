#!/usr/bin/env Rscript
# given output "angles.csv" and "pdb_angles", performs
# Kolmogorov–Smirnov test
args = commandArgs(trailingOnly=TRUE)

library(cramer)

# import and assign data
angles_exp <- read.csv(args[1], FALSE, sep = ";")
angles_pdb <- read.csv(args[2], FALSE, sep = ";")

bond_pdb <- angles_pdb[[1]]
bond_experimental <- angles_exp[[1]]

torsion_pdb <- angles_pdb[[2]]
torsion_experimental <- angles_exp[[2]]

# perform Kolmogorov-Smirnov test and store the output
#sink("output_bond.txt", append = TRUE)
print(ks.test(bond_pdb, bond_experimental))
#sink()
#sink("output_torsion.txt", append = TRUE)
print(ks.test(torsion_pdb, torsion_experimental))
#sink()

