suppressWarnings(library(data.table))
rm(list = ls())
box::use(
    util/LogUtil[log_msg],
    util/GeneUtil[string_to_adj_binary],
    src/DoIntegPIN_JMK[do_integrate_pin, compute_stochastic_matrix],
    src/CompSR_AS[comp_sr_optimized],
    future[...],
    future.apply[...]
)

#' Batch Calculate Signaling Entropy with Parallel Support
#'
#' Iterates through each sample in the expression matrix using parallel processing,
#' integrates with the STRING network, and calculates the Signaling Entropy Rate.
#'
#' @param string_path Character. Path to the STRING interaction file.
#' @param exp_path Character. Path to the expression matrix.
#' @param out_adj_path Character. Optional path to save/load precomputed adjacency matrix.
#' @param result_output_path Character. Path for the consolidated output file.
#' @param string_score_cutoff Numeric. Confidence threshold for STRING edges.
#' @param species Integer or Character. NCBI taxonomy ID or species name.
#' @param rewrite_adj Logical. If TRUE, recomputes the adjacency matrix.
#' @param n_cores Integer. Number of CPU cores for parallel processing. Default is 1 (sequential).
#'
#' @return A data frame of Sample IDs and Signaling Entropy Rates.
#' @export
signalling_entropy_parallel <- function(string_path, 
                               exp_path, 
                               out_adj_path = NULL,
                               result_output_path = "all_samples_entropy.txt",
                               string_score_cutoff = 0,
                               species = 9606,
                               rewrite_adj = FALSE,
                               n_cores = 1) {
    
    log_msg("INFO", "Initializing parallel data loading...")

    # 1. Adjacency Matrix Preparation (Single-threaded)
    if (!is.null(out_adj_path) && file.exists(out_adj_path) && !rewrite_adj) {
        log_msg("INFO", "Loading cached adjacency matrix.")
        adj_dt <- suppressWarnings(data.table::fread(out_adj_path))
        adj_mat <- as.matrix(adj_dt[, -1, with = FALSE])
        rownames(adj_mat) <- adj_dt[[1]]
    } else {
        log_msg("INFO", "Building adjacency matrix from STRING.")
        adj_mat <- string_to_adj_binary(string_path, cutoff = string_score_cutoff, 
                                        convert_to_entrez = TRUE, species = species)
        if (!is.null(out_adj_path)) {
            data.table::fwrite(data.table::as.data.table(adj_mat, keep.rownames = "GeneID"), 
                               file = out_adj_path, sep = "\t")
        }
    }

    # 2. Expression Matrix Loading
    exp_dt <- suppressWarnings(data.table::fread(exp_path))
    exp_mat <- as.matrix(exp_dt[, -1, with = FALSE])
    rownames(exp_mat) <- exp_dt[[1]]
    sample_names <- colnames(exp_mat)
    
    # 3. Parallel Backend Setup
    if (n_cores > 1) {
        log_msg("INFO", paste("Setting up parallel backend with", n_cores, "cores."))
        future::plan(future::multisession, workers = n_cores)
        # Ensure cleanup of workers on exit
        on.exit(future::plan(future::sequential), add = TRUE)
        apply_func <- future.apply::future_lapply
    } else {
        apply_func <- lapply
    }

    log_msg("INFO", paste("Processing", length(sample_names), "samples..."))

    # 4. Parallelized Iteration
    # future_lapply automatically exports necessary objects (adj_mat, etc.) to workers
    entropy_results_list <- apply_func(sample_names, function(sample_id) {

        current_exp <- exp_mat[, sample_id, drop = FALSE]
        integrated <- do_integrate_pin(adj_mat, current_exp)
        
        prob_mat <- compute_stochastic_matrix(as.vector(integrated$e), integrated$a)
        entropy_res <- comp_sr_optimized(prob_mat)
        
        return(data.frame(SampleID = sample_id, 
                          SignallingEntropyRate = entropy_res$sr, 
                          stringsAsFactors = FALSE))
    })

    # 5. Consolidation and Output
    final_entropy_df <- do.call(rbind, entropy_results_list)
    write.table(final_entropy_df, file = result_output_path, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    log_msg("INFO", paste("Task complete. Results saved to:", result_output_path))
    return(final_entropy_df)
}

signalling_entropy_parallel(string_path = "/data/pub/zhousha/Totipotent20251031/data/STRING/9606.protein.links.v12.0.txt",
          exp_path = "/data/pub/zhousha/Totipotent20251031/RNAseqML/matrix/human_all_tpm.tsv",
          out_adj_path = "/data/pub/zhousha/Totipotent20251031/data/STRING/9606_adj_matrix.tsv",
          rewrite_adj = TRUE,
          n_cores = 1,
          result_output_path = "/data/pub/zhousha/Totipotent20251031/RNAseqML/signalling/signalling_entropy_local.txt")