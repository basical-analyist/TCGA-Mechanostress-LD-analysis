# ============================================================
# TCGA PanCancer Atlas 2018 — compact analysis
#
# Outputs
#   1) BRCA scatter A: LD buffering vs Turnover
#      point color = Mechanostress
#   2) BRCA scatter B: Mechanostress vs Adjusted LD excess
#      Adjusted LD excess = Z(residuals(lm(LD_buffering ~ Turnover)))
#   3) Four-cohort Mechanostress correlation heatmaps
#      color only / rho + raw P / rho + raw-P stars
#
# Heatmap order: CRC, LUAD, BRCA, SKCM
# Scoring: hierarchical mean Z-score, using every complete sample
# in each PanCancer Atlas RNA-seq sample list (no tumor filter).
# ============================================================


# 0. Packages and settings -----------------------------------

cran_pkgs <- c("httr", "jsonlite", "ggplot2", "writexl", "scales")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(ggplot2)
  library(writexl)
  library(scales)
})

make_output_dir <- function(base_dir) {
  if (!dir.exists(base_dir)) return(NULL)
  candidate <- file.path(base_dir, "TCGA_compact_scatter_heatmap_results")
  candidate_cache <- file.path(candidate, "cache")
  suppressWarnings(dir.create(candidate_cache, recursive = TRUE, showWarnings = FALSE))
  if (!dir.exists(candidate_cache)) return(NULL)
  write_test <- tempfile(pattern = "write_test_", tmpdir = candidate_cache)
  writable <- suppressWarnings(file.create(write_test))
  if (isTRUE(writable)) unlink(write_test)
  if (isTRUE(writable)) candidate else NULL
}

# First use the current working directory. If it is a protected drive root
# (for example D:/), automatically fall back to Documents or the user home.
out_dir <- make_output_dir(getwd())
if (is.null(out_dir)) {
  fallback_roots <- unique(c(file.path(path.expand("~"), "Documents"), path.expand("~")))
  for (root in fallback_roots) {
    out_dir <- make_output_dir(root)
    if (!is.null(out_dir)) break
  }
}
if (is.null(out_dir)) {
  stop("A writable output directory could not be created. Set getwd() to a writable folder.")
}
cache_dir <- file.path(out_dir, "cache")
message("Output directory: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))

api_base <- "https://www.cbioportal.org/api"
rho_color_limit <- 0.5

# Requested order is preserved by this named vector.
studies <- c(
  CRC  = "coadread_tcga_pan_can_atlas_2018",
  LUAD = "luad_tcga_pan_can_atlas_2018",
  BRCA = "brca_tcga_pan_can_atlas_2018",
  SKCM = "skcm_tcga_pan_can_atlas_2018"
)


# 1. Gene sets ------------------------------------------------

gene_sets <- list(
  Mechanostress = list(
    Focal_adhesion = c("ITGB1", "PTK2", "PXN", "TLN1", "VCL"),
    Actomyosin_tension = c("RHOA", "ROCK1", "MYH9", "MYL9", "ACTN1"),
    Membrane_cytoskeleton = c("FLNA", "EZR", "MSN", "SPTAN1", "ANK3")
  ),
  LD_buffering = list(
    Neutral_lipid_synthesis = c("DGAT1", "DGAT2", "GPAM", "AGPAT2", "LPIN1", "SOAT1"),
    LD_organization = c("PLIN2", "PLIN3", "CIDEC"),
    Lipid_channeling = c("SCD", "FASN", "ACSL3")
  ),
  Turnover = list(
    Autophagosome = c("ATG5", "ATG7", "BECN1", "MAP1LC3B", "ULK1"),
    Lysosome = c("LAMP1", "LAMP2", "CTSD", "LIPA", "ATP6V1A", "ATP6V0A1"),
    Trafficking_lipolysis = c("RAB7A", "STX17", "VAMP8", "PNPLA2", "LIPE")
  )
)

required_genes <- unique(unlist(gene_sets, use.names = FALSE))
stopifnot(length(required_genes) == 43L)


# 2. Compact helpers -----------------------------------------

zscore <- function(x, label = "value") {
  x <- as.numeric(x)
  if (any(!is.finite(x)) || !is.finite(sd(x)) || sd(x) == 0) {
    stop(label, " contains invalid values or has zero variance.")
  }
  as.numeric(scale(x))
}

api_call <- function(path, query = list(), body = NULL, attempts = 3L) {
  url <- paste0(api_base, path)
  last_error <- "unknown error"

  for (attempt in seq_len(attempts)) {
    response <- tryCatch(
      if (is.null(body)) {
        GET(url, query = query, accept_json(), timeout(600))
      } else {
        # Explicit serialization is required for top-level JSON arrays
        # such as the HUGO symbols sent to /genes/fetch.
        json_body <- toJSON(
          body, auto_unbox = TRUE, na = "null", null = "null", digits = 22
        )
        raw_body <- charToRaw(enc2utf8(as.character(json_body)))
        POST(url, query = query, body = raw_body, encode = "raw",
             accept_json(), content_type_json(), timeout(600))
      },
      error = identity
    )

    if (inherits(response, "response") && status_code(response) < 300) {
      return(fromJSON(content(response, "text", encoding = "UTF-8"),
                      flatten = TRUE))
    }

    last_error <- if (inherits(response, "response")) {
      paste0(
        "HTTP ", status_code(response), ": ",
        substr(content(response, "text", encoding = "UTF-8"), 1, 300)
      )
    } else {
      conditionMessage(response)
    }
    message("API attempt ", attempt, " failed: ", last_error)
  }
  stop("cBioPortal API failed after ", attempts, " attempts: ", last_error)
}

# Resolve the 43 HUGO symbols once.
gene_info <- api_call(
  "/genes/fetch",
  query = list(geneIdType = "HUGO_GENE_SYMBOL", projection = "SUMMARY"),
  body = unname(as.character(required_genes))
)
gene_map <- setNames(as.integer(gene_info$entrezGeneId),
                     toupper(gene_info$hugoGeneSymbol))
missing_map <- setdiff(required_genes, names(gene_map))
if (length(missing_map)) stop("Entrez mapping failed: ", paste(missing_map, collapse = ", "))
entrez_ids <- unname(gene_map[required_genes])

fetch_cohort <- function(cohort, study_id) {
  cache_file <- file.path(cache_dir, paste0(cohort, "_43genes_all_samples.rds"))
  if (file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (is.list(cached) && is.matrix(cached$expression) &&
        identical(rownames(cached$expression), required_genes) &&
        ncol(cached$expression) >= 50L &&
        all(is.finite(cached$expression))) {
      message(cohort, ": using cached expression matrix (N=",
              ncol(cached$expression), ").")
      return(cached)
    }
  }

  profile_id <- paste0(study_id, "_rna_seq_v2_mrna")
  sample_ids <- unique(as.character(unlist(api_call(
    paste0("/sample-lists/", profile_id, "/sample-ids")
  ), use.names = FALSE)))

  chunks <- split(entrez_ids, ceiling(seq_along(entrez_ids) / 15))
  molecular_data <- do.call(rbind, lapply(seq_along(chunks), function(i) {
    message(cohort, ": downloading expression chunk ", i, "/", length(chunks))
    x <- api_call(
      paste0("/molecular-profiles/", profile_id, "/molecular-data/fetch"),
      query = list(projection = "SUMMARY"),
      body = list(sampleIds = sample_ids,
                  entrezGeneIds = as.integer(chunks[[i]]))
    )
    as.data.frame(x, stringsAsFactors = FALSE)
  }))

  expr <- matrix(
    NA_real_, nrow = length(required_genes), ncol = length(sample_ids),
    dimnames = list(required_genes, sample_ids)
  )
  lookup <- setNames(required_genes, as.character(entrez_ids))
  r <- match(unname(lookup[as.character(molecular_data$entrezGeneId)]),
             required_genes)
  c <- match(as.character(molecular_data$sampleId), sample_ids)
  ok <- !is.na(r) & !is.na(c)
  expr[cbind(r[ok], c[ok])] <- as.numeric(molecular_data$value[ok])

  complete <- colSums(!is.finite(expr)) == 0L
  expr <- expr[, complete, drop = FALSE]
  if (ncol(expr) < 50L) stop(cohort, ": fewer than 50 complete samples.")

  negative_n <- sum(expr < 0)
  if (negative_n) expr[expr < 0] <- 0
  q99 <- unname(quantile(expr, 0.99, na.rm = TRUE))
  transformation <- if (min(expr) >= 0 && q99 > 50) {
    expr <- log2(expr + 1)
    "log2(expression + 1)"
  } else {
    "No additional log2 transform"
  }

  result <- list(
    expression = expr,
    qc = data.frame(
      Cohort = cohort,
      Sample_list_N = length(sample_ids),
      Complete_sample_N = ncol(expr),
      Incomplete_sample_N_removed = sum(!complete),
      Negative_value_N_floored = negative_n,
      Transformation = transformation,
      Sample_rule = "All complete RNA-seq sample-list samples; no tumor filter",
      stringsAsFactors = FALSE
    )
  )
  if (!dir.exists(cache_dir)) {
    stop("Cache directory disappeared before saving: ", cache_dir)
  }
  saveRDS(result, cache_file)
  result
}

score_hierarchical <- function(expr, cohort) {
  missing <- setdiff(required_genes, rownames(expr))
  if (length(missing)) stop(cohort, " missing genes: ", paste(missing, collapse = ", "))

  # Gene-wise cohort Z-score.
  expr_z <- t(scale(t(expr[required_genes, , drop = FALSE])))
  if (any(!is.finite(expr_z))) stop(cohort, ": invalid gene-wise Z-scores.")

  # Gene mean within submodule -> equal mean of three submodules -> final Z-score.
  major_scores <- sapply(names(gene_sets), function(major) {
    submodule_scores <- sapply(gene_sets[[major]], function(genes) {
      colMeans(expr_z[genes, , drop = FALSE])
    })
    zscore(rowMeans(submodule_scores), paste(cohort, major))
  })

  scores <- data.frame(
    sample_id = colnames(expr), Cohort = cohort,
    major_scores, check.names = FALSE, stringsAsFactors = FALSE
  )
  fit <- lm(LD_buffering ~ Turnover, data = scores)
  scores$Adjusted_LD_predicted <- as.numeric(fitted(fit))
  scores$Adjusted_LD_residual <- as.numeric(residuals(fit))
  scores$Adjusted_LD_excess <- zscore(
    scores$Adjusted_LD_residual, paste(cohort, "Adjusted LD excess")
  )
  scores
}

cor_row <- function(scores, cohort, test, x, y) {
  keep <- complete.cases(scores[[x]], scores[[y]])
  ct <- suppressWarnings(cor.test(scores[[x]][keep], scores[[y]][keep],
                                  method = "spearman", exact = FALSE,
                                  alternative = "two.sided"))
  data.frame(
    Cohort = cohort, Test = test, X = x, Y = y, N = sum(keep),
    Spearman_rho = unname(ct$estimate), Raw_P_value = ct$p.value,
    stringsAsFactors = FALSE
  )
}

format_p <- function(p) {
  if (p < 2.2e-16) return("< 2.2e-16")
  if (p < 0.0001) return("< 0.0001")
  formatC(p, format = "fg", digits = 3)
}

stat_label <- function(x) {
  p_text <- format_p(x$Raw_P_value)
  paste0("Spearman \u03c1 = ", sprintf("%.3f", x$Spearman_rho),
         if (startsWith(p_text, "<")) paste0("    P ", p_text)
         else paste0("    P = ", p_text),
         "    N = ", x$N)
}


# 3. Download once and calculate all scores -----------------

cohort_data <- Map(fetch_cohort, names(studies), unname(studies))
names(cohort_data) <- names(studies)
scores <- Map(function(x, nm) score_hierarchical(x$expression, nm),
              cohort_data, names(cohort_data))
names(scores) <- names(studies)
qc_table <- do.call(rbind, lapply(cohort_data, `[[`, "qc"))
rownames(qc_table) <- NULL


# 4. BRCA scatter plots --------------------------------------

brca <- scores$BRCA
scatter_stats <- rbind(
  cor_row(brca, "BRCA", "Turnover_vs_LD_buffering", "Turnover", "LD_buffering"),
  cor_row(brca, "BRCA", "Mechanostress_vs_Adjusted_LD_excess",
          "Mechanostress", "Adjusted_LD_excess")
)

scatter_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold", colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.line = element_line(linewidth = 0.6, colour = "black"),
    axis.line.x.top = element_line(linewidth = 0.6),
    axis.line.y.right = element_line(linewidth = 0.6),
    axis.ticks.x.top = element_line(linewidth = 0.5),
    axis.ticks.y.right = element_line(linewidth = 0.5),
    axis.text.x.top = element_blank(), axis.text.y.right = element_blank(),
    axis.title.x.top = element_blank(), axis.title.y.right = element_blank(),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 1, size = 8.5,
                                 margin = margin(t = 2, b = 7)),
    plot.tag = element_text(face = "bold", size = 15)
  )

common_scatter <- list(
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.85, colour = "black"),
  scale_x_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)),
  scale_y_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)),
  scatter_theme
)

scatter_A <- ggplot(brca, aes(LD_buffering, Turnover, colour = Mechanostress)) +
  geom_point(size = 1.7, alpha = 0.72) +
  common_scatter +
  scale_colour_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, name = "Mechanostress\nscore"
  ) +
  labs(
    tag = "A", title = "PanCancer Atlas 2018 \u2014 hierarchical z-score",
    subtitle = stat_label(scatter_stats[1, ]),
    x = "LD-buffering score", y = "Turnover-competence score"
  )

scatter_B <- ggplot(brca, aes(Mechanostress, Adjusted_LD_excess)) +
  geom_point(size = 1.7, alpha = 0.62, colour = "#6F6F6F") +
  common_scatter +
  labs(
    tag = "B", title = "PanCancer Atlas 2018 \u2014 hierarchical z-score",
    subtitle = stat_label(scatter_stats[2, ]),
    x = "Mechanostress score",
    y = "Adjusted LD excess\n[Z(residuals: LD buffering ~ Turnover)]"
  )

scatter_individuals <- list(
  "01A_LD_buffering_vs_Turnover" = scatter_A,
  "01B_Mechanostress_vs_Adjusted_LD_excess" = scatter_B
)
for (nm in names(scatter_individuals)) {
  print(scatter_individuals[[nm]])
  stem <- file.path(out_dir, nm)
  ggsave(paste0(stem, ".pdf"), scatter_individuals[[nm]],
         device = grDevices::cairo_pdf,
         width = 6.8, height = 5.4, units = "in")
  ggsave(paste0(stem, ".png"), scatter_individuals[[nm]],
         width = 6.8, height = 5.4, units = "in", dpi = 600, bg = "white")
}


# 5. Four-cohort correlation heatmaps ------------------------

heatmap_rows <- c(
  LD_buffering = "LD buffering",
  Turnover = "Turnover competence",
  Adjusted_LD_excess = "Adjusted LD excess\n[Z(regression residual)]"
)

heatmap_stats <- do.call(rbind, lapply(names(studies), function(cohort) {
  do.call(rbind, lapply(names(heatmap_rows), function(endpoint) {
    x <- cor_row(scores[[cohort]], cohort,
                 paste0("Mechanostress_vs_", endpoint),
                 "Mechanostress", endpoint)
    x$Heatmap_row <- unname(heatmap_rows[endpoint])
    x
  }))
}))
rownames(heatmap_stats) <- NULL

heatmap_stats$BH_FDR_global_12_tests <- p.adjust(
  heatmap_stats$Raw_P_value, method = "BH"
)
heatmap_stats$Raw_P_significance <- cut(
  heatmap_stats$Raw_P_value,
  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
  labels = c("****", "***", "**", "*", "ns"), right = FALSE
)
heatmap_stats$Cohort <- factor(heatmap_stats$Cohort, levels = names(studies))
heatmap_stats$Heatmap_row <- factor(
  heatmap_stats$Heatmap_row, levels = rev(unname(heatmap_rows))
)
heatmap_stats$Label_P <- paste0(
  "\u03c1 = ", sprintf("%.3f", heatmap_stats$Spearman_rho), "\nP = ",
  ifelse(heatmap_stats$Raw_P_value < 0.0001,
         formatC(heatmap_stats$Raw_P_value, format = "e", digits = 1),
         formatC(heatmap_stats$Raw_P_value, format = "f", digits = 4))
)
heatmap_stats$Label_stars <- paste0(
  "\u03c1 = ", sprintf("%.3f", heatmap_stats$Spearman_rho), "\n",
  heatmap_stats$Raw_P_significance
)
heatmap_stats$Label_colour <- ifelse(
  abs(heatmap_stats$Spearman_rho) >= 0.65 * rho_color_limit, "white", "black"
)

heatmap_base <- ggplot(
  heatmap_stats, aes(Cohort, Heatmap_row, fill = Spearman_rho)
) +
  geom_tile(colour = "black", linewidth = 0.55, width = 0.98, height = 0.98) +
  scale_fill_gradient2(
    name = "Spearman\n\u03c1", low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-rho_color_limit, rho_color_limit),
    oob = squish
  ) +
  scale_x_discrete(position = "top", drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  coord_fixed(ratio = 0.82) +
  labs(
    title = "Correlation with Mechanostress score",
    subtitle = "TCGA PanCancer Atlas 2018 | hierarchical mean Z-score",
    x = NULL, y = NULL,
    caption = paste0("Adjusted LD excess = Z(residuals from within-cohort ",
                     "LD buffering ~ Turnover regression)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", size = 12, colour = "black"),
    axis.text.y = element_text(face = "bold", size = 11, colour = "black"),
    plot.title = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold"),
    plot.caption = element_text(size = 8.5, hjust = 0),
    plot.margin = margin(10, 15, 10, 15)
  )

heatmaps <- list(
  NO_STATS = heatmap_base,
  RHO_RAW_P = heatmap_base +
    geom_text(aes(label = Label_P, colour = Label_colour),
              size = 3.25, lineheight = 0.95, fontface = "bold",
              show.legend = FALSE) + scale_colour_identity(),
  RHO_P_STARS = heatmap_base +
    geom_text(aes(label = Label_stars, colour = Label_colour),
              size = 3.4, lineheight = 0.95, fontface = "bold",
              show.legend = FALSE) + scale_colour_identity() +
    labs(caption = paste0(
      "Adjusted LD excess = Z(regression residual); ",
      "**** P<0.0001, *** P<0.001, ** P<0.01, * P<0.05, ns P>=0.05"
    ))
)

for (nm in names(heatmaps)) {
  print(heatmaps[[nm]])
  stem <- file.path(out_dir, paste0("02_4cohort_heatmap_", nm))
  ggsave(paste0(stem, ".pdf"), heatmaps[[nm]],
         device = grDevices::cairo_pdf,
         width = 8.2, height = 4.9, units = "in")
  ggsave(paste0(stem, ".png"), heatmaps[[nm]],
         width = 8.2, height = 4.9, units = "in", dpi = 600, bg = "white")
}


# 6. Tables ---------------------------------------------------

gene_set_table <- do.call(rbind, lapply(names(gene_sets), function(major) {
  do.call(rbind, lapply(names(gene_sets[[major]]), function(submodule) {
    data.frame(Major_score = major, Submodule = submodule,
               Gene = gene_sets[[major]][[submodule]], stringsAsFactors = FALSE)
  }))
}))
rownames(gene_set_table) <- NULL

# Convert factors back to text for clean spreadsheet/CSV export.
heatmap_export <- heatmap_stats
heatmap_export$Cohort <- as.character(heatmap_export$Cohort)
heatmap_export$Heatmap_row <- as.character(heatmap_export$Heatmap_row)
heatmap_export$Raw_P_significance <- as.character(heatmap_export$Raw_P_significance)

write.csv(scatter_stats, file.path(out_dir, "Scatter_correlations.csv"), row.names = FALSE)
write.csv(heatmap_export, file.path(out_dir, "Heatmap_correlations.csv"), row.names = FALSE)
write.csv(qc_table, file.path(out_dir, "Cohort_QC.csv"), row.names = FALSE)
for (nm in names(scores)) {
  write.csv(scores[[nm]], file.path(out_dir, paste0("Scores_", nm, ".csv")),
            row.names = FALSE)
}

write_xlsx(
  c(
    list(Scatter_correlations = scatter_stats,
         Heatmap_correlations = heatmap_export,
         Cohort_QC = qc_table,
         Gene_sets = gene_set_table),
    setNames(scores, paste0("Scores_", names(scores)))
  ),
  file.path(out_dir, "TCGA_compact_scatter_heatmap_results.xlsx")
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done. Results: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
