bait_example = function() {
    # Each of the last four samples differs from the first in just one key.
    samples = data.frame(
        Mixture = factor(c("m1", "m2", "m1", "m1", "m1")),
        TechRepMixture = c(1L, 1L, 2L, 1L, 1L),
        Run = factor(c("r1", "r1", "r1", "r2", "r1")),
        Channel = factor(c("126", "126", "126", "126", "127")),
        Condition = factor(c("A", "B", "A", "A", "B")),
        BioReplicate = c("a1", "b1", "a1", "a1", "b2")
    )
    proteins = samples[rep(seq_len(nrow(samples)), each = 3L), ]
    proteins$Protein = factor(rep(c("bait|id", "same", "higher"), 5L))
    proteins$Abundance = log2(c(256, 1024, 1024, 512, 2048, 4096,
                               128, 1024, 64, 64, 1024, 32, 32, 1024, NA))
    proteins$QC = rep(c(TRUE, FALSE, NA), 5L)
    list(FeatureLevelData = data.frame(Intensity = c(12, 20)),
         ProteinLevelData = proteins)
}

test_that("bait normalization calculates log2 ratios using every key", {
    quant = bait_example()
    quant$ProteinLevelData = quant$ProteinLevelData[15:1, ]
    adjusted = normalizeToBaitTMT(quant, "bait|id")
    proteins = adjusted$ProteinLevelData
    audit = attr(adjusted, "BaitNormalization")

    expect_identical(proteins$Abundance, rev(c(2, 2, 2, 3, 3, -1, 4, -1, 5, NA)))
    expect_equal(proteins$Abundance[8] - proteins$Abundance[10], 0)
    expect_equal(proteins$Abundance[7] - proteins$Abundance[9], 1)
    keep = which(quant$ProteinLevelData$Protein != "bait|id")
    attr(proteins, "MSstatsTMTAbundanceScale") = NULL
    expect_identical(proteins[names(proteins) != "Abundance"],
                     quant$ProteinLevelData[keep, names(proteins) != "Abundance"])
    expect_identical(adjusted$FeatureLevelData, quant$FeatureLevelData)
    expect_identical(audit$bait, "bait|id")
    expect_identical(audit$matching_keys, c("Mixture", "TechRepMixture", "Run", "Channel"))
    expect_identical(audit$missing_bait, "error")
    expect_identical(audit$original_abundance, quant$ProteinLevelData$Abundance)
    expect_identical(audit$retained_rows, keep)
    expect_identical(audit$bait_abundance$Abundance, c(5, 6, 7, 9, 8))
    expect_equal(nrow(audit$dropped_samples), 0)
})

test_that("data.table inputs are unchanged and the audit survives serialization", {
    quant = lapply(bait_example(), data.table::as.data.table)
    before = data.table::copy(quant)
    adjusted = normalizeToBaitTMT(quant, "bait|id")

    expect_identical(quant, before)
    expect_s3_class(adjusted$ProteinLevelData, "data.table")
    expect_identical(adjusted$FeatureLevelData, quant$FeatureLevelData)
    data.table::set(adjusted$ProteinLevelData, i = 1L, j = "Abundance", value = -10)
    expect_identical(quant, before)

    restored = unserialize(serialize(adjusted, NULL))
    expect_identical(attr(restored, "BaitNormalization"),
                     attr(adjusted, "BaitNormalization"))
    expect_error(normalizeToBaitTMT(restored, "same"), "already been applied")
})

test_that("repacking the tables preserves the normalization and plotting guards", {
    for (quant in list(bait_example(), lapply(bait_example(), data.table::as.data.table))) {
        adjusted = normalizeToBaitTMT(quant, "bait|id")
        repacked = adjusted[c("FeatureLevelData", "ProteinLevelData")]
        repacked = unserialize(serialize(repacked, NULL))

        expect_null(attr(repacked, "BaitNormalization"))
        expect_error(normalizeToBaitTMT(repacked, "same"), "already been applied")
        expect_error(dataProcessPlotsTMT(repacked, "ProfilePlot"), "original protein summary")
        expect_error(dataProcessPlotsTMT(repacked, "QCPlot"), "original protein summary")
    }
})

test_that("invalid input and entirely absent exact bait fail explicitly", {
    quant = bait_example()
    expect_error(normalizeToBaitTMT(quant$ProteinLevelData, "bait|id"),
                 "FeatureLevelData and ProteinLevelData")
    expect_error(normalizeToBaitTMT(quant["ProteinLevelData"], "bait|id"),
                 "FeatureLevelData and ProteinLevelData")
    expect_error(normalizeToBaitTMT(quant, c("same", "higher")), "one nonempty")
    expect_error(normalizeToBaitTMT(quant, NA_character_), "one nonempty")
    expect_error(normalizeToBaitTMT(quant, " "), "one nonempty")
    expect_error(normalizeToBaitTMT(quant, "bait", "drop"),
                 "absent from ProteinLevelData")
    quant$ProteinLevelData$Abundance = as.character(quant$ProteinLevelData$Abundance)
    expect_error(normalizeToBaitTMT(quant, "bait|id"), "numeric vector")
    quant$ProteinLevelData$Abundance = NULL
    expect_error(normalizeToBaitTMT(quant, "bait|id"), "Missing ProteinLevelData columns: Abundance")
})

test_that("invalid matching keys and ambiguous observations fail", {
    for (key in c("Mixture", "TechRepMixture", "Run", "Channel")) {
        quant = bait_example()
        quant$ProteinLevelData[[key]][1] = NA
        expect_error(normalizeToBaitTMT(quant, "bait|id"),
                     paste0("Missing or blank values in ProteinLevelData\\$", key))
    }
    for (row in c(1L, 2L)) {
        quant = bait_example()
        quant$ProteinLevelData = quant$ProteinLevelData[c(1:15, row), ]
        expect_error(normalizeToBaitTMT(quant, "bait|id"), "Duplicate protein observation")
    }
    quant = bait_example()
    quant$ProteinLevelData$Condition[1] = "B"
    expect_error(normalizeToBaitTMT(quant, "bait|id"), "Conflicting Condition or BioReplicate")
    quant = bait_example()
    quant$ProteinLevelData$BioReplicate[1] = "other"
    expect_error(normalizeToBaitTMT(quant, "bait|id"), "Conflicting Condition or BioReplicate")
    quant$ProteinLevelData$BioReplicate[1] = NA
    expect_error(normalizeToBaitTMT(quant, "bait|id"), "Missing or blank BioReplicate")
})

test_that("missing bait errors by default and drop reports whole sample losses", {
    quant = bait_example()
    quant$ProteinLevelData$Abundance[c(4, 7)] = c(NA, Inf)
    quant$ProteinLevelData = quant$ProteinLevelData[-1, ]
    expect_error(normalizeToBaitTMT(quant, "bait|id"),
                 "m1 / 1 / r1 / 126 / absent.*m2 / 1 / r1 / 126 / non-finite")
    expect_warning(adjusted <- normalizeToBaitTMT(quant, "bait|id", "drop"),
                   "3 biological sample/channel.*Dropped 6 client row")
    audit = attr(adjusted, "BaitNormalization")
    expect_identical(audit$missing_bait, "drop")
    expect_identical(audit$dropped_samples$Reason, c("absent", "non-finite", "non-finite"))
    expect_identical(audit$dropped_samples$BaitAbundance, c(NA_real_, NA_real_, Inf))
    expect_identical(audit$dropped_samples$ClientRows, rep(2L, 3L))
    expect_identical(adjusted$ProteinLevelData$Abundance, c(4, -1, 5, NA_real_))
})

test_that("numeric channel matching and annotation checks ignore data.table rounding", {
    previous_rounding = data.table::getNumericRounding()
    data.table::setNumericRounding(2L)
    tryCatch({
        quant = bait_example()
        p = quant$ProteinLevelData
        p$Channel = as.double(as.character(p$Channel))
        p$Channel[3] = 126 + 1e-12
        quant$ProteinLevelData = p

        expect_error(normalizeToBaitTMT(quant, "bait|id"), "1 biological sample/channel.*absent")
        expect_warning(adjusted <- normalizeToBaitTMT(quant, "bait|id", "drop"),
                       "Dropped 1 client row")
        audit = attr(adjusted, "BaitNormalization")
        expect_identical(audit$dropped_samples$Channel, p$Channel[3])
        expect_identical(audit$dropped_samples$Reason, "absent")
        expect_identical(adjusted$ProteinLevelData$Abundance[1], 2)

        p$Channel[3] = p$Channel[1]
        p$BioReplicate = rep(1, nrow(p))
        p$BioReplicate[1] = 1 + 1e-12
        quant$ProteinLevelData = p
        expect_error(normalizeToBaitTMT(quant, "bait|id"), "Conflicting Condition or BioReplicate")
        expect_equal(data.table::getNumericRounding(), 2L)
    }, finally = data.table::setNumericRounding(previous_rounding))
})

test_that("drop requires usable client data and rejects all non-finite baits", {
    quant = bait_example()
    quant$ProteinLevelData$Abundance[seq(1, 15, 3)] = c(NA, NaN, Inf, -Inf, NA)
    expect_error(normalizeToBaitTMT(quant, "bait|id"),
                 "5 biological sample/channel")
    expect_error(normalizeToBaitTMT(quant, "bait|id", "drop"),
                 "No usable client abundances remain")
    quant = bait_example()
    quant$ProteinLevelData$Abundance[quant$ProteinLevelData$Protein != "bait|id"] = NA
    expect_error(normalizeToBaitTMT(quant, "bait|id"),
                 "No usable client abundances remain")
    quant$ProteinLevelData$Abundance[2] = Inf
    expect_error(normalizeToBaitTMT(quant, "bait|id"), "Infinite adjusted client abundance")
})

test_that("special channels stay on the original scale in the audit only", {
    quant = bait_example()
    p = quant$ProteinLevelData
    levels(p$Condition) = c("A", "B", "Norm", "Empty")
    p$Condition[1:6] = rep(c("Norm", "Empty"), each = 3L)
    p$BioReplicate[1:6] = NA
    p$Abundance[4] = NA
    p = p[-1, ]  # The reference has no bait row at all.
    quant$ProteinLevelData = p
    adjusted = normalizeToBaitTMT(quant, "bait|id")
    audit = attr(adjusted, "BaitNormalization")
    expect_identical(audit$special_channels, p[1:5, ])
    expect_equal(nrow(audit$dropped_samples), 0)
    expect_false(any(adjusted$ProteinLevelData$Condition %in% c("Norm", "Empty")))
    expect_identical(adjusted$FeatureLevelData, quant$FeatureLevelData)
    expect_error(dataProcessPlotsTMT(adjusted, "ProfilePlot"), "original protein summary")
    expect_error(dataProcessPlotsTMT(adjusted, "QCPlot"), "original protein summary")
})

test_that("group comparison on bait ratios agrees with manual subtraction", {
    quant = MSstatsTMT::quant.pd.msstats
    p = quant$ProteinLevelData
    # Treat each mixture/channel as one biological subject across technical runs.
    p$BioReplicate = paste(p$Mixture, p$Channel)
    quant$ProteinLevelData = p
    bait = "P04406"
    adjusted = normalizeToBaitTMT(quant, bait)

    key = with(p, paste(Mixture, TechRepMixture, Run, Channel))
    is_bait = p$Protein == bait
    manual = quant
    manual$ProteinLevelData = p[!is_bait, ]
    manual$ProteinLevelData$Abundance = p$Abundance[!is_bait] -
        p$Abundance[is_bait][match(key[!is_bait], key[is_bait])]
    adjusted_proteins = adjusted$ProteinLevelData
    attr(adjusted_proteins, "MSstatsTMTAbundanceScale") = NULL
    expect_identical(adjusted_proteins, manual$ProteinLevelData)

    contrast = matrix(c(-1, 0, 0, 1), nrow = 1,
                      dimnames = list("1-0.125", c("0.125", "0.5", "0.667", "1")))
    actual = groupComparisonTMT(adjusted, contrast, remove_norm_channel = FALSE,
                                remove_empty_channel = FALSE, use_log_file = FALSE,
                                verbose = FALSE)$ComparisonResult
    expected = groupComparisonTMT(manual, contrast, use_log_file = FALSE,
                                  verbose = FALSE)$ComparisonResult
    expect_equal(actual, expected)
    expect_equal(nrow(actual), 9)
    expect_false(bait %in% actual$Protein)
    expect_true(all(is.finite(actual$log2FC)))
    expect_true(all(is.finite(actual$pvalue)))
})
