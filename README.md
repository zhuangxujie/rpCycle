# rpCycle
To calculated the relative proportion of cells in each cell cycle phase for 414 samples based on gene expression profile. 

@param expr the gene expression data set. A data frame with row names as symbols and columns as samples.
@param marker a list of cell cycle phase marker symbol genes.
@return a matrix, row is samples, col is cell cycle phase.


An example:
source("R/rpCycle.R")

load("data/expr.rda")

res <- rpCycle(expr)

res
