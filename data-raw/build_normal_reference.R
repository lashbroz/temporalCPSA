## Build normal_reference.rds for ageTMP
##
## The normal/reference data used for AD-TMP trajectory analysis were
## previously stored in the paper-analysis object `for_tadj.RData`.
##
## Source/provenance:
## - DEveLopmental Trajectory Atlas (DELTA), dorsolateral prefrontal cortex
##   (DLPFC)
## - PMID: 30518843
## - URL: http://amp.pharm.mssm.edu/DELTA
##
## This script converts the legacy study object into documented package data.

legacy_file <- Sys.getenv(
  "AGETMP_FOR_TADJ",
  "/Users/lashbn01/Dropbox/HOPE_AYA/for_tadj.RData"
)
normal_cluster_file <- Sys.getenv(
  "AGETMP_NORMAL_PROTEIN_CLUSTERS",
  "/Users/lashbn01/Dropbox/HOPE_AYA/normal_protein_clusters_scaled.RData"
)

if (!file.exists(legacy_file)) {
  stop("Legacy for_tadj.RData file not found: ", legacy_file, call. = FALSE)
}

loaded_objects <- load(legacy_file)

protein_clusters <- NULL
if (file.exists(normal_cluster_file)) {
  normal_cluster_env <- new.env(parent = emptyenv())
  cluster_objects <- load(normal_cluster_file, envir = normal_cluster_env)
  if ("n.con.mfclusters" %in% cluster_objects) {
    raw_cluster <- factor(
      normal_cluster_env$n.con.mfclusters$col.clusters[[4]]$consensusClass
    )
    remapped_cluster <- raw_cluster
    levels(remapped_cluster) <- list("3" = "1", "4" = "2", "1" = "3", "2" = "4")
    protein_clusters <- list(
      consensus_k4_raw = raw_cluster,
      consensus_k4 = remapped_cluster,
      remap = c("1" = "3", "2" = "4", "3" = "1", "4" = "2"),
      source_file = normal_cluster_file,
      source_object = "n.con.mfclusters$col.clusters[[4]]$consensusClass"
    )
  }
}

normal_reference <- list(
  protein = list(
    matrix = breen.prot,
    sample_metadata = breen.prot.meta,
    clusters = protein_clusters
  ),
  transcript = list(
    matrix = breen.trans,
    sample_metadata = breen.trans.meta
  ),
  brainspan = list(
    annotation = bs.anno,
    matrix = bs.data,
    sample_metadata = bs.meta
  ),
  provenance = list(
    source_name = "DEveLopmental Trajectory Atlas (DELTA), DLPFC",
    pmid = "30518843",
    url = "http://amp.pharm.mssm.edu/DELTA",
    legacy_file = legacy_file,
    normal_cluster_file = if (file.exists(normal_cluster_file)) normal_cluster_file else NA_character_,
    legacy_objects = loaded_objects,
    build_date = as.character(Sys.Date())
  )
)

out_dir <- file.path("temporalCPSA", "inst", "extdata")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

saveRDS(normal_reference, file.path(out_dir, "normal_reference.rds"))
