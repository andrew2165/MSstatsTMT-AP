context("groupComparision")

test_that("groupComparision works", {
  
  output<-groupComparisonTMT(data = MSstatsTMTAP::quant.pd.msstats, moderated = TRUE)
  expect_equal(output$ComparisonResult, MSstatsTMTAP::test.pairwise$ComparisonResult, tolerance=1e-5)
  
})

test_that("groupComparision handle wrong input", {
  
  expect_error(groupComparisonTMT())
  expect_error(groupComparisonTMT(data = MSstatsTMTAP::quant.pd.msstats$ProteinLevelData))
  expect_error(groupComparisonTMT(data = MSstatsTMTAP::quant.pd.msstats$FeatureLevelData))
  expect_error(groupComparisonTMT(data = MSstatsTMTAP::quant.pd.msstats,
                                   contrast.matrix = matrix(c(1,2,3,4),nrow = 2)))
  expect_error(groupComparisonTMT(data = MSstatsTMTAP::quant.pd.msstats,remove_norm_channel = "abc"))
  expect_error(groupComparisonTMT(data = MSstatsTMTAP::quant.pd.msstats,
                                   moderated = "abc"))
  
})

