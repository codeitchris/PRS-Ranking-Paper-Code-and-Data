### Instructions
## AA Matrix of the comparion graph
## WW Matrix of winner indices: encodes observed winner information
## Output reference guide: matrix RR
## -----------------------------------------------------------------    
## RR[1,] theta.hat
## RR[2,] rank based on theta hat
## RR{3,] & RR[4,] two-sided CI for rank                               ## RR[1,] to RR[6,] are based on the Vanilla spectral method
## RR[5,] left-sided CI for rank
## RR[6,] uniform left-sided CI for rank
## -----------------------------------------------------------------    
## RR[7,] theta.hat
## RR[8,] rank based on theta hat
## RR{9,] & RR[10,] two-sided CI for rank                              ## RR[7,] to RR[12,] are based on the two-stage theta estimators
## RR[11,] left-sided CI for rank
## RR[12,] uniform left-sided CI for rank

### Read AUC results, compare and form graph ###
library(readr)
library(dplyr)

current_dir <- getwd()


### calculation ###
output_dir <- file.path(current_dir, "spectral_ranking")

aCSV <- read_csv(file.path("YOUR PATH TO AA0 OF APPLIED OR METHOD PAPERS"))
AA0 <- as.matrix(aCSV)

wCSV <- read_csv(file.path("YOUR PATH TO WW0 OF APPLIED OR METHOD PAPERS"))
WW0 <- as.matrix(wCSV)

AA2 <- AA0
WW2 <- WW0

## Vanilla spectral method ##

B <- 2000
n <- ncol(AA2)
L2 <- nrow(AA2)
fAvec2 <- numeric(L2) + 2 ## weight for vanilla spectral method ##



## -----------------------------------------------------------------    ## compute matrix P ##
dval2 <- 2 * max(colSums(AA2)) ## value of d ##
P2 <- matrix(0, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    if (j != i) {
      P2[i, j] <- sum(AA2[, i] * AA2[, j] * WW2[, j] / fAvec2) / dval2
    }
  }
  P2[i, i] <- 1 - sum(P2[i, ])
}


## -----------------------------------------------------------------    ## solve theta and pi ##
tmp.P2 <- t(t(P2) - diag(n)) %*% (t(P2) - diag(n))
tmp.svd2 <- svd(tmp.P2)
pihat2 <- abs(tmp.svd2$v[, n])
print(pihat2)
thetahat2 <- log(pihat2) - mean(log(pihat2))
print(thetahat2)

## -----------------------------------------------------------------    ## output ##
RR2 <- matrix(0, 12, n)
colnames(RR2) <- aCSV[0,]
RR2[1, ] <- thetahat2
RR2[2, ] <- n + 1 - rank(thetahat2)

Vmatrix2 <- matrix(0, L2, n)
tauhatvec2 <- numeric(n)
tmp.pimatrix2 <- t(AA2) * pihat2
tmp.pivec2 <- colSums(tmp.pimatrix2)
tmp.var2 <- numeric(n)
for (oo in 1:n) {
  tauhatvec2[oo] <- sum(AA2[, oo] * (1 - pihat2[oo] / tmp.pivec2) * pihat2[oo] / fAvec2, na.rm = TRUE) / dval2
  print(tauhatvec2[oo])
  tmp.var2[oo] <- sum(AA2[, oo] * (tmp.pivec2 - pihat2[oo]) / fAvec2 / fAvec2) * pihat2[oo] / dval2 / dval2 / tauhatvec2[oo] / tauhatvec2[oo]
  Vmatrix2[, oo] <- (AA2[, oo] * WW2[, oo] * tmp.pivec2 - AA2[, oo] * pihat2[oo]) / fAvec2
}
sigmahatmatrix2 <- matrix(tmp.var2, n, n) + t(matrix(tmp.var2, n, n))

## -----------------------------------------------------------------    ## Weighted bootstrap ##
Wmatrix2 <- matrix(rnorm(L2 * B), L2, B)
tmp.Vtau2 <- (t(Vmatrix2) / tauhatvec2) %*% Wmatrix2

## -----------------------------------------------------------------    ##
R.left.m2 <- numeric(n)
R.right.m2 <- numeric(n)
R.left.one.m2 <- numeric(n)

for (ooo in 1:n) {
  print(ooo)
  tmpGMmatrix02 <- matrix(rep(tmp.Vtau2[ooo, ], n) - c(t(tmp.Vtau2)), B, n)
  tmpGMmatrix2 <- abs(t(t(tmpGMmatrix02) / sqrt(sigmahatmatrix2[ooo, ])) / dval2)
  tmpGMmatrixone2 <- t(t(tmpGMmatrix02) / sqrt(sigmahatmatrix2[ooo, ])) / dval2
  tmp.GMvecmax2 <- apply(tmpGMmatrix2, 1, max)
  tmp.GMvecmaxone2 <- apply(tmpGMmatrixone2, 1, max)
  cutval2 <- quantile(tmp.GMvecmax2, 0.95)
  cutvalone2 <- quantile(tmp.GMvecmaxone2, 0.95)
  tmp.theta.sd2 <- sqrt(sigmahatmatrix2[ooo, ])
  tmp.theta.sd2 <- tmp.theta.sd2[-ooo]
  R.left.m2[ooo] <- 1 + sum(1 * (((thetahat2[-ooo] - thetahat2[ooo]) / tmp.theta.sd2) > cutval2))
  R.right.m2[ooo] <- n - sum(1 * (((thetahat2[-ooo] - thetahat2[ooo]) / tmp.theta.sd2) < (-cutval2)))
  R.left.one.m2[ooo] <- 1 + sum(1 * (((thetahat2[-ooo] - thetahat2[ooo]) / tmp.theta.sd2) > cutvalone2))
}

## two-sided CI for rank ##
RR2[3, ] <- R.left.m2
RR2[4, ] <- R.right.m2
## left-sided CI for rank ##
RR2[5, ] <- R.left.one.m2


## -----------------------------------------------------------------    ## Weighted bootstrap ##
Wmatrix2 <- matrix(rnorm(L2 * B), L2, B)
tmp.Vtau2 <- (t(Vmatrix2) / tauhatvec2) %*% Wmatrix2
Mval <- n
GMvecmax2 <- numeric(B) - 1
GMvecmaxone2 <- numeric(B) - Inf
tmpTMval2 <- -1
tmpTMvalone2 <- -Inf
## -----------------------------------------------------------------    ##
for (ooo in 1:n) {
  tmpGMmatrix02 <- matrix(rep(tmp.Vtau2[ooo, ], n) - c(t(tmp.Vtau2)), B, n)
  tmpGMmatrix2 <- abs(t(t(tmpGMmatrix02) / sqrt(sigmahatmatrix2[ooo, ])) / dval2)
  tmpGMmatrixone2 <- t(t(tmpGMmatrix02) / sqrt(sigmahatmatrix2[ooo, ])) / dval2
  tmp.GMvecmax2 <- apply(tmpGMmatrix2, 1, max)
  tmp.GMvecmaxone2 <- apply(tmpGMmatrixone2, 1, max)
  GMvecmax2 <- c(GMvecmax2, tmp.GMvecmax2)
  GMvecmaxone2 <- c(GMvecmaxone2, tmp.GMvecmaxone2)
}
## -----------------------------------------------------------------    ##
GMmaxmatrixone2 <- matrix(GMvecmaxone2, B)
GMmaxone2 <- apply(GMmaxmatrixone2, 1, max)
cutvalone2 <- quantile(GMmaxone2, 0.95)
## -----------------------------------------------------------------    ##
R.left.one2 <- numeric(n)
for (oooo in 1:n) {
  tmp.theta.sd2 <- sqrt(sigmahatmatrix2[oooo, ])
  tmp.theta.sd2 <- tmp.theta.sd2[-oooo]
  R.left.one2[oooo] <- 1 + sum(1 * (((thetahat2[-oooo] - thetahat2[oooo]) / tmp.theta.sd2) > cutvalone2))
}
RR2[6, ] <- R.left.one2 ## uniform left-sided CI for rank ##



## two-stage theta estimators ##

fMLE2 <- numeric(L2)
ftwo2 <- rowSums(AA2 * t(matrix(exp(thetahat2), n, nrow(AA2))))
fMLE2 <- ftwo2 ## weight for two-stage spectral method ##


## -----------------------------------------------------------------    ## compute matrix P ##
dval2 <- 2 * max(colSums(AA2))
PMLE2 <- matrix(0, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    if (j != i) {
      PMLE2[i, j] <- sum(AA2[, i] * AA2[, j] * WW2[, j] / fMLE2, na.rm = TRUE) / dval2
    }
  }
  PMLE2[i, i] <- 1 - sum(PMLE2[i, ])
}


## -----------------------------------------------------------------    ## solve theta and pi ##
tmp.PMLE2 <- t(t(PMLE2) - diag(n)) %*% (t(PMLE2) - diag(n))
print(PMLE2)

tmp.svd.MLE2 <- svd(tmp.PMLE2)
pihatMLE2 <- abs(tmp.svd.MLE2$v[, n])
thetahatMLE2 <- log(pihatMLE2) - mean(log(pihatMLE2))
RR2[7, ] <- thetahatMLE2 ## theta hat ##
RR2[8, ] <- n + 1 - rank(thetahatMLE2) ## rank based on theta hat ##


## -----------------------------------------------------------------    ## compute sample xihat ##
pihat2 <- pihatMLE2
thetahat2 <- thetahatMLE2
P2 <- PMLE2
fAvec2 <- fMLE2
Vmatrix2 <- matrix(0, L2, n)
tauhatvec2 <- numeric(n)
tmp.pimatrix2 <- t(AA2) * pihat2
tmp.pivec2 <- colSums(tmp.pimatrix2)
tmp.var2 <- numeric(n)
for (oo in 1:n) {
  tauhatvec2[oo] <- sum(AA2[, oo] * (1 - pihat2[oo] / tmp.pivec2) * pihat2[oo] / fAvec2, na.rm = TRUE) / dval2
  tmp.var2[oo] <- sum(AA2[, oo] * (tmp.pivec2 - pihat2[oo]) / fAvec2 / fAvec2, na.rm = TRUE) * pihat2[oo] / dval2 / dval2 / tauhatvec2[oo] / tauhatvec2[oo]
  Vmatrix2[, oo] <- (AA2[, oo] * WW2[, oo] * tmp.pivec2 - AA2[, oo] * pihat2[oo]) / fAvec2
}
sigmahatmatrix2 <- matrix(tmp.var2, n, n) + t(matrix(tmp.var2, n, n))
## -----------------------------------------------------------------    ## Weighted bootstrap ##
Wmatrix2 <- matrix(rnorm(L2 * B), L2, B)
Vmatrix2[is.na(Vmatrix2)] <- 0
Vmatrix2[is.nan(Vmatrix2)] <- 0
tmp.Vtau2 <- (t(Vmatrix2) / tauhatvec2) %*% Wmatrix2

## -----------------------------------------------------------------    ##
R.left.m2 <- numeric(n)
R.right.m2 <- numeric(n)
R.left.one.m2 <- numeric(n)
for (ooo in 1:n) {
  print(ooo)
  tmpGMmatrix02 <- matrix(rep(tmp.Vtau2[ooo, ], n) - c(t(tmp.Vtau2)), B, n)
  tmpGMmatrix2 <- abs(t(t(tmpGMmatrix02) / sqrt(sigmahatmatrix2[ooo, ])) / dval2)
  tmpGMmatrixone2 <- t(t(tmpGMmatrix02) / sqrt(sigmahatmatrix2[ooo, ])) / dval2
  tmp.GMvecmax2 <- apply(tmpGMmatrix2, 1, max)
  tmp.GMvecmaxone2 <- apply(tmpGMmatrixone2, 1, max)
  cutval2 <- quantile(tmp.GMvecmax2, 0.95)
  cutvalone2 <- quantile(tmp.GMvecmaxone2, 0.95)
  tmp.theta.sd2 <- sqrt(sigmahatmatrix2[ooo, ])
  tmp.theta.sd2 <- tmp.theta.sd2[-ooo]
  R.left.m2[ooo] <- 1 + sum(1 * (((thetahat2[-ooo] - thetahat2[ooo]) / tmp.theta.sd2) > cutval2), na.rm = TRUE)
  R.right.m2[ooo] <- n - sum(1 * (((thetahat2[-ooo] - thetahat2[ooo]) / tmp.theta.sd2) < (-cutval2)), na.rm = TRUE)
  R.left.one.m2[ooo] <- 1 + sum(1 * (((thetahat2[-ooo] - thetahat2[ooo]) / tmp.theta.sd2) > cutvalone2), na.rm = TRUE)
}

## two-sided CI for rank ##
RR2[9, ] <- R.left.m2
RR2[10, ] <- R.right.m2
## left-sided CI for rank ##
RR2[11, ] <- R.left.one.m2

## -----------------------------------------------------------------    ## Weighted bootstrap ##
Wmatrix2 <- matrix(rnorm(L2 * B), L2, B)
tmp.Vtau2 <- (t(Vmatrix2) / tauhatvec2) %*% Wmatrix2
Mval <- n
GMvecmax2 <- numeric(B) - 1
GMvecmaxone2 <- numeric(B) - Inf
tmpTMval2 <- -1
tmpTMvalone2 <- -Inf

## -----------------------------------------------------------------    ##
for (ooo in 1:n) {
  tmpGMmatrix02 <- matrix(rep(tmp.Vtau2[ooo, ], n) - c(t(tmp.Vtau2)), B, n)
  tmpGMmatrix2 <- abs(t(t(tmpGMmatrix02) / sqrt(sigmahatmatrix2[ooo, ])) / dval2)
  tmpGMmatrixone2 <- t(t(tmpGMmatrix02) / sqrt(sigmahatmatrix2[ooo, ])) / dval2
  tmp.GMvecmax2 <- apply(tmpGMmatrix2, 1, max)
  tmp.GMvecmaxone2 <- apply(tmpGMmatrixone2, 1, max)
  GMvecmax2 <- c(GMvecmax2, tmp.GMvecmax2)
  GMvecmaxone2 <- c(GMvecmaxone2, tmp.GMvecmaxone2)
}
GMmaxmatrixone2 <- matrix(GMvecmaxone2, B)
GMmaxone2 <- apply(GMmaxmatrixone2, 1, max)
cutvalone2 <- quantile(GMmaxone2, 0.95)
R.left.one2 <- numeric(n)
for (oooo in 1:n) {
  tmp.theta.sd2 <- sqrt(sigmahatmatrix2[oooo, ])
  tmp.theta.sd2 <- tmp.theta.sd2[-oooo]
  R.left.one2[oooo] <- 1 + sum(1 * (((thetahat2[-oooo] - thetahat2[oooo]) / tmp.theta.sd2) > cutvalone2))
}
RR2[12, ] <- R.left.one2 ## uniform left-sided CI for rank ##


colnames(RR2) <- colnames(AA0)
## RR[1,] theta.hat
## RR[2,] rank based on theta hat
## RR{3,] & RR[4,] two-sided CI for rank                               ## RR[1,] to RR[6,] are based on the Vanilla spectral method
## RR[5,] left-sided CI for rank
## RR[6,] uniform left-sided CI for rank
## -----------------------------------------------------------------    
## RR[7,] theta.hat
## RR[8,] rank based on theta hat
## RR{9,] & RR[10,] two-sided CI for rank                              ## RR[7,] to RR[12,] are based on the two-stage theta estimators
## RR[11,] left-sided CI for rank
## RR[12,] uniform left-sided CI for rank
rownames(RR2) <-c("theta.hat","rank","leftOf2CI",
                  "rightOf2CI","left-sidedCI","unifleft-sidedCI",
                  "2theta.hat","secondRank","2leftOf2CI",
                  "2rightOf2CI","2left-sidedCI","2unifleft-sidedCI")


dates<-c(2007,2019.12,2015.10,2020.12,2020.12,2020.12,2021.10,2017.06,
         2017.05,2022.10,2019.04,2019.04,2019.11,2020.05)


library(ggplot2)
dfT <- as.data.frame(t(RR2))
dfT |>
  ggplot(aes(x =rank, y = reorder(rownames(dfT),dates, decreasing=TRUE) )) +
  geom_errorbar(aes(xmin = `leftOf2CI`, xmax = `rightOf2CI`),color="black",size=.2) +
  geom_point(color = "red", size = 1) +
  labs(title  = "Method Paper Ranking",
       y = "Method",
       x = "Rank"  ) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) + 
  theme_linedraw() + geom_text(aes(label = rank), vjust = -.5,color="blue")


library(gt)
library(tidyverse)


# Traits table
table <- read_csv(file.path("YOUR PATH TO TRAIT COMPARISON FROM APPLIED PAPERS"))

myTable<-table |> 
  gt(id = "mygt") |>
  tab_spanner(
    label = 'Methods',
    columns=colnames(RR2)
  )|>
  cols_label(traits~"Traits")|>
  tab_header(title="Number of Method Comparisons per Trait for Applied Papers")|>
  opt_row_striping(row_striping=TRUE) |>
  tab_options(latex.use_longtable = TRUE)|>
  gt_split(row_every_n = 50)

myTable
gtsave(grp_pull(myTable,1), filename = "LongTable1.png",vwidth = 2000)

tbl<-rm_spanners(grp_pull(myTable,2), spanners = everything(), levels = NULL)|>
  rm_header()
gtsave(tbl, filename = "LongTable2.png",vwidth = 2000)

tbl<-rm_spanners(grp_pull(myTable,3), spanners = everything(), levels = NULL)|>
  rm_header()
gtsave(tbl, filename = "LongTable3.png",vwidth = 2000)

# Head to head table
table <- read_csv(file.path("YOUR PATH TO COMPARISON DF APPLIED OR METHOD PAPERS"))
rownames(table) = colnames(table)
secondTable<-table |>
  gt(rownames_to_stub = TRUE) |>
  opt_row_striping(row_striping=TRUE) |>
  tab_header(title="Head to Head Record of Methods for Applied Papers")
secondTable
gtsave(secondTable, filename = "h2hcomp.png",vwidth = 2500)

# Violin Plot
library(ggstatsplot)
library(gghalves) 
library(tidyverse)
library(paletteer)

table <- read_csv(file.path("YOUR PATH TO VIOLON DATA FROM APPLIED PAPERS"))
myPalette <- colorBlindness::paletteMartin
ggbetweenstats(
  data = table,
  x = Method,
  y = Ranking,
  violin.args = list(width = 0, linewidth = 0),
  point.args = list(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6)),
  ggtheme = ggplot2::theme_bw(),
  results.subtitle = FALSE,
  pairwise.display = "none",
  package = "colorBlindness",
  palette = "paletteMartin",
  title = "Normalized Ranking Over Individual Traits"
) +
  geom_half_violin(data = table, aes(x = Method, y = Ranking,fill=Method),
                   side = "right", # Add half-violins to the right
                   position = position_nudge(x = 0.15), # Adjust position
                   alpha = 0.6) +
   scale_fill_paletteer_d("colorBlindness::paletteMartin")



