load("data/genemarker.rda")

#' The cell cycle phase prportion in samples
#' @description calculated the relative proportion of cells in each cell cycle phase for samples based on gene expression profile.
#'
#' @param expr the gene expression data set. A data frame with row names as symbols and columns as samples.
#' @param marker a list of cell cycle phase marker symbol genes.
#'
#' @return a matrix, row is samples, col is cell cycle phase.
#' @export
#'
#' @examples
rpCycle <- function(expr, marker = genemarker) {

  if (length(unique(rownames(expr))) != nrow(expr)) {
    stop("The gene symbols in expression are not unique!")
  }

  # anti-log if max < 50 in mixture file
  if(max(expr) < 50) {
    expr <- 2^expr
  }

  # normalize expression matrix
  normalize <- function(x) {
    return((x - min(x)) / (max(x) - min(x)))
  }

  expr.pro <- apply(expr, 2, normalize) * 10

  # update marker gene list, delete genes that not in expression
  marker.new <- lapply(marker, function(x) {
    x <- intersect(rownames(expr.pro), x)
  })

  # get sample by cell cycle phase matrix
  avg <- as.data.frame(do.call(cbind,
                               lapply(marker.new, function(x){
                                 apply(expr.pro[x, ], 2, mean)
                               })))

  # calculate cor for marker and all cell cycle phase
  marker.end <- list()
  corScore <- list()
  for (i in 1:length(marker.new)){
    phase <- c()
    score <- c()
    for (j in 1:length(marker.new[[i]])){
      temp <- cor.test(as.numeric(expr.pro[marker.new[[i]][j], ]),
                       as.numeric(avg[, i]))
      if (!is.na(temp$estimate)) { # & temp$estimate > 0
        phase <- c(phase, marker.new[[i]][j])
        score <- c(score, temp$estimate)
      }
    }
    marker.end <- c(marker.end, list(phase))
    corScore <- c(corScore, list(score))
  }
  names(corScore) <- names(marker.new)
  names(marker.end) <- names(marker.new)

  # get cell cycle phase score in every sample
  res <- matrix(0L, length(marker.end), ncol(expr.pro))
  for (i in 1:length(marker.end)) {
    for (j in 1:ncol(expr.pro)) {
      res[i, j] <- sum(as.numeric(expr.pro[marker.end[[i]], j]) * (1 + exp(as.numeric(corScore[[i]]))))
    }
  }
  rownames(res) <- names(marker.end)
  colnames(res) <- colnames(expr.pro)

  # make sure the RES is greater than 0
  out <- res - apply(res, 1, function(x){floor(min(x))})

  # get every cell phase proportion in samples
  out <- apply(out, 2, function(x){x / sum(x)})

  # matrix-row as samples, col as cell cycle phase
  return(t(out))
}
