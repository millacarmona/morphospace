##################################################################

test_that(desc = "testing phy_prcomp, general behavior", code = {
  data(tails)

  Y <- geomorph::two.d.array(expected_shapes(tails$shapes, tails$data$species))
  tree <- tails$tree
  ppca <- phy_prcomp(Y, tree)

  result1 <- nrow(ppca$x) == nrow(Y)
  result2 <- ncol(ppca$x) == min(ncol(Y), nrow(Y) - 1)

  #note: phyl.pca gives as many vectors as min(ncol(Y), nrow(Y) - 1); however, if
  #one make computations by hand, the amount is min(ncol(Y), nrow(Y)). Idk the reason.

  result3 <- round(sum(apply(Y,2,var)),10) == round(sum(apply(ppca$x,2,var)),10)
  result4 <- sum(apply(ppca$x,2,stats::sd)) > sum(ppca$sdev)

  result5 <- all(round(ppca$center,10) == round(apply(Y, 2, phytools::fastAnc, tree = tree)[1,],10))

  expect_true(all(result1,result2,result3,result4,result5))
})

##################################################################

test_that(desc = "testing phyalign_comp, general behavior", code = {
  data(tails)

  Y <- geomorph::two.d.array(expected_shapes(tails$shapes, tails$data$species))
  tree <- tails$tree
  paca <- phyalign_comp(Y, tree)

  result1 <- nrow(paca$x) == nrow(Y)
  result2 <- ncol(paca$x) == min(ncol(Y), nrow(Y))

  result3 <- round(sum(apply(Y,2,var)),10) == round(sum(apply(paca$x,2,var)),10)
  result4 <- sum(apply(paca$x,2,stats::sd)) > sum(paca$sdev)

  result5 <- all(round(paca$center,10) == round(apply(Y, 2, phytools::fastAnc, tree = tree)[1,],10))

  expect_true(all(result1,result2,result3,result4,result5))
})

##################################################################

test_that(desc = "testing bg_prcomp, general behavior", code = {
  data(tails)

  Y <- geomorph::two.d.array(tails$shapes)
  X <- tails$data$species
  bgpca <- bg_prcomp(x = Y, groups = X)
  model_ndims <- ncol(model.matrix(~X)) - 1

  result1 <- nrow(bgpca$x) == nrow(Y)
  result2 <- ncol(bgpca$x) == min(model_ndims, nrow(Y))

  result3 <- round(sum(apply(Y,2,var)),10) >= round(sum(apply(bgpca$x,2,var)),10)

  result4 <- all(round(bgpca$grcenter,10) == round( apply(X=Y,MARGIN=2,FUN=tapply,X,mean)  ,10))
  result5 <- all(round(bgpca$center,10) == round(colMeans(Y),10))
  result6 <- all(round(colMeans(bgpca$x),10) == 0)

  expect_true(all(result1,result2,result3,result4,result5,result6))
})


test_that(desc = "testing bg_prcomp, general behavior (with loocv)", code = {
  data(tails)

  Y <- geomorph::two.d.array(tails$shapes)
  X <- tails$data$species
  bgpca <- bg_prcomp(x = Y, groups = X, LOOCV = TRUE)
  model_ndims <- ncol(model.matrix(~X)) - 1

  result1 <- nrow(bgpca$x) == nrow(Y)
  result2 <- ncol(bgpca$x) == min(model_ndims, nrow(Y))

  result3 <- round(sum(apply(Y,2,var)),10) >= round(sum(apply(bgpca$x,2,var)),10)

  result4 <- all(round(bgpca$grcenter,10) == round( apply(X=Y,MARGIN=2,FUN=tapply,X,mean)  ,10))
  result5 <- all(round(bgpca$center,10) == round(colMeans(Y),10))
  result6 <- all(round(colMeans(bgpca$x),10) == 0)

  expect_true(all(result1,result2,result3,result4,result5,result6))
})

##################################################################

test_that(desc = "testing pls2b, general behavior", code = {
  data(tails)

  Y <- geomorph::two.d.array(tails$shapes)
  X <- tails$sizes
  pls <- pls2b(y = Y, x = X)
  model_ndimsX <- ncol(model.matrix(~X)) - 1
  model_ndimsY <- ncol(model.matrix(~Y)) - 1

  result1 <- nrow(pls$yscores) == nrow(Y)
  result2 <- ncol(pls$yscores) == min(model_ndimsX, model_ndimsY, nrow(Y))

  result3 <- round(sum(apply(Y,2,var)),10) >= round(sum(apply(pls$yscores,2,var)),10)

  result4 <- all(round(pls$ycenter,10) == round(colMeans(Y),10))
  result5 <- all(round(colMeans(pls$yscores),10) == 0)


  result6 <- nrow(pls$xscores) == nrow(cbind(X))
  result7 <- ncol(pls$xscores) == min(model_ndimsX, model_ndimsY, nrow(cbind(X)))

  result8 <- round(sum(apply(cbind(X),2,var)),10) >= round(sum(apply(pls$xscores,2,var)),10)

  result9 <- all(round(pls$xcenter,10) == round(colMeans(cbind(X)),10))
  result10 <- all(round(colMeans(pls$xscores),10) == 0)

  expect_true(all(result1,result2,result3,result4,result5,
                  result6,result7,result8,result9,result10))
})


test_that(desc = "testing (phylogenetic) pls2b, general behavior (with loocv)", code = {
  data(tails)

  Y <- geomorph::two.d.array(tails$shapes)
  X <- tails$sizes
  pls <- pls2b(y = Y, x = X, LOOCV = TRUE)
  model_ndimsX <- ncol(model.matrix(~X)) - 1
  model_ndimsY <- ncol(model.matrix(~Y)) - 1

  result1 <- nrow(pls$yscores) == nrow(Y)
  result2 <- ncol(pls$yscores) == min(model_ndimsX, model_ndimsY, nrow(Y))

  result3 <- round(sum(apply(Y,2,var)),10) >= round(sum(apply(pls$yscores,2,var)),10)

  result4 <- all(round(pls$ycenter,10) == round(colMeans(Y),10))
  result5 <- all(round(colMeans(pls$yscores),10) == 0)


  result6 <- nrow(pls$xscores) == nrow(cbind(X))
  result7 <- ncol(pls$xscores) == min(model_ndimsX, model_ndimsY, nrow(cbind(X)))

  result8 <- round(sum(apply(cbind(X),2,var)),10) >= round(sum(apply(pls$xscores,2,var)),10)

  result9 <- all(round(pls$xcenter,10) == round(colMeans(cbind(X)),10))
  result10 <- all(round(colMeans(pls$xscores),10) == 0)

  expect_true(all(result1,result2,result3,result4,result5,
                  result6,result7,result8,result9,result10))
})


##################################################################

test_that(desc = "testing burnaby, general behavior", code = {
  data(tails)

  Y <- geomorph::two.d.array(tails$shapes)
  X <- tails$sizes
  burn1 <- burnaby(x = Y, vars = X)

  result1 <- nrow(burn1$x) == nrow(Y)
  result2 <- round(sum(apply(Y,2,var)),10) >= round(sum(burn1$sdev^2),10)

  result3 <- all(round(burn1$center,10) == round(colMeans(Y),10))
  result4 <- all(round(colMeans(burn1$x),10) == 0)

  result5 <- ncol(burn1$x) == min(ncol(Y) - 1, nrow(Y))

  effdims <- 1:(((ncol(Y)-4)/2)-1)
  allo_shapes <- expected_shapes(shapes = geomorph::arrayspecs(Y, k =2, p = 9), x = cbind(X), returnarray = FALSE)
  nallo_scores <- proj_eigen(x = allo_shapes, vectors = burn1$rotation, center = burn1$center)[,effdims]

  result6 <- all(round(apply(nallo_scores, 2, var),10) == 0)
  result7 <- round(Morpho::angle.calc(burn1$rotation[,1], lm(allo_shapes ~ X)$coef[2,])*57.2958, 3) == 90


  axmat <- lm(Y ~ X)$coef[2,]
  cent <- colMeans(lm(Y ~ X)$fitted)
  burn2 <- burnaby(x = Y, axmat = axmat, center = cent)

  result8 <- nrow(burn2$x) == nrow(Y)
  result9 <- round(sum(apply(Y,2,var)),10) >= round(sum(burn2$sdev^2),10)

  result10 <- all(round(burn2$center,10) == round(colMeans(Y),10))
  result11 <- all(round(colMeans(burn2$x),10) == 0)

  result12 <- ncol(burn2$x) == min(ncol(Y) - 1, nrow(Y))

  effdims <- 1:(((ncol(Y)-4)/2)-1)
  allo_shapes <- expected_shapes(shapes = geomorph::arrayspecs(Y, k =2, p = 9), x = cbind(X), returnarray = FALSE)
  nallo_scores <- proj_eigen(x = allo_shapes, vectors = burn2$rotation, center = burn2$center)[,effdims]

  result13 <- all(round(apply(nallo_scores, 2, var),10) == 0)
  result14 <- round(Morpho::angle.calc(burn2$rotation[,1], lm(allo_shapes ~ X)$coef[2,])*57.2958, 3) == 90

  expect_true(all(result1,result2,result3,result4,result5,result6,result7,
                  result8,result9,result10,result11,result12,result13,result14))
})


test_that(desc = "testing (phylogenetic) burnaby, general behavior", code = {
  data(tails)

  Y <- geomorph::two.d.array(expected_shapes(tails$shapes, tails$data$species))
  X <- cbind(tapply(tails$sizes, tails$data$species, mean))
  tree <- tails$tree
  pburn1 <- burnaby(x = Y, vars = X, tree = tree)

  result1 <- nrow(pburn1$x) == nrow(Y)
  result2 <- round(sum(apply(Y,2,var)),10) >= round(sum(pburn1$sdev^2),10)

  result3 <- all(round(pburn1$center,10) == round(apply(Y, 2, phytools::fastAnc, tree = tree)[1,],10))
  result4 <- ncol(pburn1$x) == min(ncol(Y) - 1, nrow(Y))

  effdims <- 1:(((ncol(Y)-4)/2)-1)
  allo_shapes <- expected_shapes(shapes = geomorph::arrayspecs(Y, k =2, p = 9), x = cbind(X), tree = tree, returnarray = FALSE)
  nallo_scores <- proj_eigen(x = allo_shapes, vectors = pburn1$rotation, center = pburn1$center)[,effdims]

  result5 <- all(round(apply(nallo_scores, 2, var),10) == 0)
  result6 <- round(Morpho::angle.calc(pburn1$rotation[,1], lm(allo_shapes ~ X)$coef[2,])*57.2958, 3) == 90


  allo_shapes <- expected_shapes(shapes = geomorph::arrayspecs(Y, k =2, p = 9), x = cbind(X),
                                 tree = tree, returnarray = FALSE)
  axmat <- adapt_model(geomorph::procD.pgls(allo_shapes ~ X, phy = tree))$coef[2,]
  cent <- colMeans(adapt_model(geomorph::procD.pgls(allo_shapes ~ X, phy = tree))$fitted)
  pburn2 <- suppressWarnings(burnaby(x = Y, axmat = axmat, center = cent))

  result7 <- nrow(pburn2$x) == nrow(Y)
  result8 <- round(sum(apply(Y,2,var)),10) >= round(sum(pburn2$sdev^2),10)


  result9 <- all(round(pburn2$center,10) == round(colMeans(mvMORPH::mvgls(allo_shapes ~ X, tree = tree, model = "BM")$fitted),10))
  result10 <- ncol(pburn2$x) == min(ncol(Y) - 1, nrow(Y))

  effdims <- 1:(((ncol(Y)-4)/2)-1)
  nallo_scores <- proj_eigen(x = allo_shapes, vectors = pburn2$rotation, center = pburn2$center)[,effdims]

  result11 <- all(round(apply(nallo_scores, 2, var),10) == 0)
  result12 <- round(Morpho::angle.calc(pburn2$rotation[,1], lm(allo_shapes ~ X)$coef[2,])*57.2958, 3) == 90

  expect_true(all(result1,result2,result3,result4,result5,result6,result7,result8,result9,result10,
                  result11,result12))
})


