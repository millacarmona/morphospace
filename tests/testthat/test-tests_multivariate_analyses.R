##################################################################

test_that(desc = "testing phy_prcomp, general behavior", code = {
  data(tails)

  Y <- geomorph::two.d.array(consensus(tails$shapes, tails$data$species))
  tree <- tails$tree
  ppca <- phy_prcomp(Y, tree)

  result1 <- nrow(ppca$x) == nrow(Y)
  result2 <- ncol(ppca$x) == min(ncol(Y), nrow(Y) - 1)

  result3 <- round(sum(apply(Y,2,var)),10) == round(sum(apply(ppca$x,2,var)),10)
  result4 <- sum(apply(ppca$x,2,stats::sd)) > sum(ppca$sdev)

  result5 <- all(round(ppca$center,10) == round(apply(Y, 2, phytools::fastAnc, tree = tree)[1,],10))

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


test_that(desc = "testing pls2b, general behavior (with loocv)", code = {
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

test_that(desc = "testing pls2b, general behavior", code = {
  data(tails)

  Y <- geomorph::two.d.array(consensus(tails$shapes, tails$data$species))
  X <- tapply(tails$sizes, tails$data$species, mean)
  tree <- tails$tree
  ppls <- pls2b(y = Y, x = X, tree)
  model_ndimsX <- ncol(model.matrix(~X)) - 1
  model_ndimsY <- ncol(model.matrix(~Y)) - 1

  result1 <- nrow(ppls$yscores) == nrow(Y)
  result2 <- ncol(ppls$yscores) == min(model_ndimsX, model_ndimsY, nrow(Y))

  result3 <- round(sum(apply(Y,2,var)),10) >= round(sum(apply(ppls$yscores,2,var)),10)

  result4 <- all(round(ppls$ycenter,10) == round(apply(Y, 2, phytools::fastAnc, tree = tree)[1,],10))


  result5 <- nrow(ppls$xscores) == nrow(cbind(X))
  result6 <- ncol(ppls$xscores) == min(model_ndimsX, model_ndimsY, nrow(cbind(X)))

  result7 <- round(sum(apply(cbind(X),2,var)),10) >= round(sum(apply(ppls$xscores,2,var)),10)

  result8 <- all(round(ppls$xcenter,10) == round(apply(cbind(X), 2, phytools::fastAnc, tree = tree)[1,],10))


  expect_true(all(result1,result2,result3,result4,result5,
                  result6,result7,result8))
})


test_that(desc = "testing pls2b, general behavior (with loocv)", code = {
  data(tails)

  Y <- geomorph::two.d.array(consensus(tails$shapes, tails$data$species))
  X <- tapply(tails$sizes, tails$data$species, mean)
  tree <- tails$tree
  ppls <- pls2b(y = Y, x = X, tree, LOOCV = TRUE)
  model_ndimsX <- ncol(model.matrix(~X)) - 1
  model_ndimsY <- ncol(model.matrix(~Y)) - 1

  result1 <- nrow(ppls$yscores) == nrow(Y)
  result2 <- ncol(ppls$yscores) == min(model_ndimsX, model_ndimsY, nrow(Y))

  result3 <- round(sum(apply(Y,2,var)),10) >= round(sum(apply(ppls$yscores,2,var)),10)

  result4 <- all(round(ppls$ycenter,10) == round(apply(Y, 2, phytools::fastAnc, tree = tree)[1,],10))


  result5 <- nrow(ppls$xscores) == nrow(cbind(X))
  result6 <- ncol(ppls$xscores) == min(model_ndimsX, model_ndimsY, nrow(cbind(X)))

  result7 <- round(sum(apply(cbind(X),2,var)),10) >= round(sum(apply(ppls$xscores,2,var)),10)

  result8 <- all(round(ppls$xcenter,10) == round(apply(cbind(X), 2, phytools::fastAnc, tree = tree)[1,],10))


  expect_true(all(result1,result2,result3,result4,result5,
                  result6,result7,result8))
})

