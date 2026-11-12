context("PDtoMSstatsTMTFormat")
context("MaxQtoMSstatsTMTFormat")
context("SpectroMinetoMSstatsTMTFormat")

test_that("PDtoMSstatsTMTFormat works", {

    expect_error(PDtoMSstatsTMTFormat(input = MSstatsTMTAP::raw.pd[, !colnames(MSstatsTMTAP::raw.pd) == "Protein.Accessions"], # missing columns in input
                                    annotation = MSstatsTMTAP::annotation.pd))

    expect_error(PDtoMSstatsTMTFormat(input = MSstatsTMTAP::raw.pd,
                                      annotation = MSstatsTMTAP::annotation.pd[, !colnames(MSstatsTMTAP::annotation.pd) == "Condition"])) # missing columns in annotation

})

test_that("MaxQtoMSstatsTMTFormat works", {

    expect_error(MaxQtoMSstatsTMTFormat(evidence = MSstatsTMTAP::evidence, proteinGroups = MSstatsTMTAP::proteinGroups,
                                        annotation = MSstatsTMTAP::annotation.pd)) # wrong annotation file

})

test_that("SpectroMinetoMSstatsTMTFormat works", {

    expect_error(SpectroMinetoMSstatsTMTFormat(input = MSstatsTMTAP::raw.mine,
                                        annotation = MSstatsTMTAP::annotation.mine,
                                        summaryforMultipleRows = average)) # wrong argument value

})
