example(transfer_labels, echo = TRUE)
test_that("outputted NMF factors are all nonnegative", {
  source_factors <- t(res$nmf$h)
  expect_equal(sum(source_factors >= 0), length(source_factors))
})

# test_that("predicted probabilities are in [0,1]",{
#
#
#   vis_anno_target <- res$targets
#   multinom_mod <- res$multinom
#   projections <- reducedDim(vis_anno_target, "nmf_projections")
#
#   preds <- predict(multinom_mod, newx=projections, type='class',
#                                     na.action=na.exclude)
#   expect_equal(sum(preds >= 0 & preds <= 1), length(preds))
#   expect_equal(as.vector(rowSums(preds)), rep(1, nrow(preds)))
# })
