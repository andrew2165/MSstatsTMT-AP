#' Normalize protein summaries to an affinity-purification bait
#'
#' Subtract the bait's log2 abundance from each client protein's log2
#' abundance within the same Mixture, TechRepMixture, Run, and Channel.
#' The result is a log2 client-to-bait ratio. Apply this after protein
#' summarization and any reference-channel alignment, before group comparison.
#'
#' @param data Output of \code{\link{proteinSummarization}}.
#' @param bait One exact identifier from \code{ProteinLevelData$Protein}.
#' @param missing_bait Either \code{"error"} (default) or \code{"drop"}.
#'   If a biological channel containing client rows has no bait row or has a
#'   non-finite bait abundance, stop or warn and drop all clients in that
#'   channel, respectively. A bait identifier absent from the entire input
#'   always causes an error.
#'   If dropping channels removes an entire condition, pairwise group
#'   comparisons omit contrasts involving that condition.
#'
#' @details
#'   Bait rows and channels annotated \code{"Norm"} or \code{"Empty"} are
#'   excluded from the returned protein table. The original special-channel
#'   rows, when present, are kept in the audit for quality control. Missing
#'   client abundances and negative log2 ratios are retained. Analyze bait-free
#'   biological controls separately.
#'
#'   Use the original summarization output with
#'   \code{\link{dataProcessPlotsTMT}}. Feature-level measurements retain
#'   their original scale. Contrasts from \code{\link{groupComparisonTMT}}
#'   describe differences in client recovery per recovered bait. Binding
#'   affinity, absolute stoichiometry, and interaction specificity require
#'   additional evidence. The model estimates variation across replicates without
#'   explicitly propagating uncertainty in the bait summaries. Assess bait
#'   peptide consistency, interference, missingness, imputation, and saturation
#'   before use. A finite bait summary alone is insufficient for this assessment.
#'
#'   Saving and loading the result with \code{saveRDS()} and \code{readRDS()}
#'   preserves its normalization markers. Rewrapping the same tables keeps the
#'   protein table marker but discards the outer audit. Rebuilding the protein
#'   table, for example with \code{data.table::rbindlist()}, can also discard
#'   its marker. If both markers are lost, the plotting and repeat-normalization
#'   checks can no longer identify the data as bait-normalized.
#'
#' @return A list with unchanged \code{FeatureLevelData} and bait-adjusted
#'   \code{ProteinLevelData}. The \code{BaitNormalization} attribute marks
#'   the object as normalized and contains \code{bait}, \code{matching_keys},
#'   \code{missing_bait}, \code{original_abundance} in input row order,
#'   \code{retained_rows} giving the input row indices of the output protein
#'   rows, \code{bait_abundance} with the original bait rows' sample annotations
#'   and abundances, \code{dropped_samples} with sample annotations, bait
#'   abundance, reason, and client-row count for each sample dropped due to
#'   unusable bait, and \code{special_channels} with the original protein rows
#'   for the reference and empty channels. To change the bait or policy, start
#'   from the original summary.
#'
#' @export
#' @seealso \code{\link{proteinSummarization}}, \code{\link{groupComparisonTMT}}
#' @examples
#' data(quant.pd.msstats)
#' adjusted <- normalizeToBaitTMT(quant.pd.msstats, bait = "P04406")
#' head(adjusted$ProteinLevelData)
#' head(attr(adjusted, "BaitNormalization")$bait_abundance)
normalizeToBaitTMT = function(data, bait, missing_bait = c("error", "drop")) {
    missing_bait = match.arg(missing_bait)
    if (!is.list(data) ||
        !is.data.frame(data[["FeatureLevelData"]]) ||
        !is.data.frame(data[["ProteinLevelData"]])) {
        stop("Supply a proteinSummarization() result with FeatureLevelData and ProteinLevelData.",
             call. = FALSE)
    }
    if (.hasBaitNormalization(data)) {
        stop("Bait normalization has already been applied. Use the original protein summary.",
             call. = FALSE)
    }
    if (!is.character(bait) || length(bait) != 1L ||
        is.na(bait) || !nzchar(trimws(bait))) {
        stop("Supply one nonempty bait protein identifier.", call. = FALSE)
    }

    keys = c("Mixture", "TechRepMixture", "Run", "Channel")
    sample_cols = c(keys, "Condition", "BioReplicate")
    required = c(sample_cols, "Protein", "Abundance")
    proteins = data[["ProteinLevelData"]]
    missing_cols = setdiff(required, names(proteins))
    if (length(missing_cols)) {
        stop("Missing ProteinLevelData columns: ", paste(missing_cols, collapse = ", "),
             call. = FALSE)
    }
    if (!is.numeric(proteins$Abundance) || !is.null(dim(proteins$Abundance))) {
        stop("ProteinLevelData$Abundance must be a numeric vector of log2 abundances.",
             call. = FALSE)
    }
    input = data.table::as.data.table(proteins)
    for (column in c(keys, "Protein", "Condition")) {
        value = as.character(input[[column]])
        if (anyNA(value) || any(!nzchar(trimws(value)))) {
            stop("Missing or blank values in ProteinLevelData$", column,
                 ". Check the sample annotation and protein identifiers.", call. = FALSE)
        }
    }
    is_bait = as.character(input$Protein) == bait
    if (!any(is_bait)) {
        stop("Bait protein ", sQuote(bait), " is absent from ProteinLevelData. ",
             "Check the exact Protein identifier and summarization filters.",
             call. = FALSE)
    }
    observation_cols = c(sample_cols, "Protein")
    observations = input[, observation_cols, with = FALSE]
    # Use exact integer codes so data.table's numeric rounding cannot merge identifiers.
    for (column in observation_cols) {
        value = observations[[column]]
        if (is.double(value) && !inherits(value, "integer64")) {
            value = unclass(value)
            data.table::set(observations, j = column, value = match(value, unique(value)))
        }
    }
    duplicate = anyDuplicated(observations, by = c(keys, "Protein"))
    if (duplicate) {
        stop("Duplicate protein observation for ", sQuote(input$Protein[duplicate]),
             " at input row ", duplicate, ". Resolve duplicate rows for each ",
             "Protein, Mixture, TechRepMixture, Run, and Channel.", call. = FALSE)
    }

    is_special = input$Condition %in% c("Norm", "Empty")
    subject = as.character(input$BioReplicate[!is_special])
    if (anyNA(subject) || any(!nzchar(trimws(subject)))) {
        stop("Missing or blank BioReplicate values in biological channels.",
             call. = FALSE)
    }
    sample_rows = which(!duplicated(observations, by = sample_cols))
    samples = observations[sample_rows, keys, with = FALSE]
    if (anyDuplicated(samples, by = keys)) {
        stop("Conflicting Condition or BioReplicate annotations for the same ",
             "Mixture, TechRepMixture, Run, and Channel. Check the sample annotation.",
             call. = FALSE)
    }

    bait_rows = input[is_bait, c(sample_cols, "Abundance"), with = FALSE]
    bait_match = observations[is_bait][samples, on = keys, which = TRUE]
    sample_bait = bait_rows$Abundance[bait_match]
    client_rows = which(!is_bait & !is_special)
    client_sample = samples[observations[client_rows], on = keys, which = TRUE]
    missing_sample = !is.finite(sample_bait)
    drop_client = missing_sample[client_sample]

    client_count = tabulate(client_sample, nrow(samples))
    dropped = which(missing_sample & client_count > 0L)
    dropped_samples = as.data.frame(input[sample_rows[dropped], sample_cols, with = FALSE])
    dropped_samples$BaitAbundance = sample_bait[dropped]
    dropped_samples$Reason = rep("non-finite", length(dropped))
    dropped_samples$Reason[is.na(bait_match[dropped])] = "absent"
    dropped_samples$ClientRows = client_count[dropped]
    if (length(dropped)) {
        shown = seq_len(min(length(dropped), 10L))
        channels = apply(dropped_samples[shown, c(keys, "Reason")], 1L,
                         paste, collapse = " / ")
        details = paste(channels, collapse = "; ")
        if (length(dropped) > length(shown)) {
            details = paste0(details, "; ...")
        }
        problem = paste0("Bait ", sQuote(bait), " is absent or non-finite in ",
                         length(dropped), " biological sample/channel combination(s) ",
                         "(Mixture / TechRepMixture / Run / Channel / reason): ", details, ".")
        if (missing_bait == "error") {
            stop(problem, call. = FALSE)
        }
    }

    retained = client_rows[!drop_client]
    adjusted = as.double(input$Abundance[retained]) -
        sample_bait[client_sample[!drop_client]]
    if (any(is.infinite(adjusted))) {
        stop("Infinite adjusted client abundance. Check the client and bait summaries.",
             call. = FALSE)
    }
    if (!any(is.finite(adjusted))) {
        stop(if (length(dropped)) paste(problem, "No usable client abundances remain.")
             else "No usable client abundances remain.", call. = FALSE)
    }
    if (length(dropped)) {
        warning(problem, " Dropped ", sum(drop_client), " client row(s). ",
                "See the BaitNormalization audit for the full list.", call. = FALSE)
    }

    data$ProteinLevelData = proteins[retained, , drop = FALSE]
    data$ProteinLevelData$Abundance = adjusted
    attr(data$ProteinLevelData, "MSstatsTMTAbundanceScale") = "log2_client_bait"
    attr(data, "BaitNormalization") = list(
        bait = bait, matching_keys = keys, missing_bait = missing_bait,
        original_abundance = input$Abundance, retained_rows = retained,
        bait_abundance = as.data.frame(bait_rows),
        dropped_samples = dropped_samples,
        special_channels = as.data.frame(proteins[which(is_special), , drop = FALSE])
    )
    data
}

.hasBaitNormalization = function(data) {
    !is.null(attr(data, "BaitNormalization", exact = TRUE)) ||
        (is.list(data) && identical(
            attr(data[["ProteinLevelData"]], "MSstatsTMTAbundanceScale", exact = TRUE),
            "log2_client_bait"))
}
