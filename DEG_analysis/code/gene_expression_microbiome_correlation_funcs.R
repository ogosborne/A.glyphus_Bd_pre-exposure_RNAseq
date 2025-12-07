# Function to run bootstrapped CCA
bootstrap_cca <- function(ex, mb, perm_obj, K = 1, nboot = 100, seed = 123){
  set.seed(seed)
  n_genes <- ncol(ex)
  n_mb <- ncol(mb)
  # matrices to store selections and weights
  ex_sel <- matrix(0, nrow = n_genes, ncol = nboot)
  mb_sel <- matrix(0, nrow = n_mb, ncol = nboot)
  ex_wts <- matrix(0, nrow = n_genes, ncol = nboot)
  mb_wts <- matrix(0, nrow = n_mb, ncol = nboot)
  # bootstraps
  for(b in 1:nboot){
    idx <- sample(1:nrow(ex), replace = TRUE)  # resample samples
    ex_b <- ex[idx, ]
    mb_b <- mb[idx, ]
    # run CCA on bootstrap sample
    cca_b <- PMA::CCA(ex_b, mb_b,
                      typex = "standard", typez = "standard",
                      K = K, niter = 500,
                      penaltyx = perm_obj$bestpenaltyx,
                      penaltyz = perm_obj$bestpenaltyz,
                      v = perm_obj$v.init,
                      trace = FALSE)
    # record selection (non-zero) and weights
    ex_sel[, b] <- as.numeric(cca_b$u != 0)
    mb_sel[, b] <- as.numeric(cca_b$v != 0)
    ex_wts[, b] <- cca_b$u
    mb_wts[, b] <- cca_b$v
  }
  # get summaries
  ex_df <- data.frame(
    feature = colnames(ex),
    selection_freq = rowMeans(ex_sel),
    mean_abs_weight = rowMeans(abs(ex_wts)),
    mean_signed_weight = rowMeans(ex_wts)
  )
  mb_df <- data.frame(
    feature = colnames(mb),
    selection_freq = rowMeans(mb_sel),
    mean_abs_weight = rowMeans(abs(mb_wts)),
    mean_signed_weight = rowMeans(mb_wts)
  )
  return(list(ex = ex_df, mb = mb_df,
              ex_wts = ex_wts, mb_wts = mb_wts))
}

# function to get all-by-all correlations for selected expression + microbiome features
allcors_ex_mb <- function(dat, only_mb_ex = TRUE, TPs = c("TP1","TP2")){
  out <- list()
  for(tp in TPs){
    my.ex <- dat[[tp]]$ex
    my.mb <- dat[[tp]]$mb
    my.cors <- Hmisc::rcorr(my.ex, my.mb, type = "pearson")
    if(only_mb_ex == FALSE){
      my.cors$Q <- my.cors$P
      my.cors$Q[upper.tri(my.cors$Q)] <- p.adjust(my.cors$Q[upper.tri(my.cors$Q)], method = "fdr")
      my.cors$Q[lower.tri(my.cors$Q)] <- p.adjust(my.cors$Q[lower.tri(my.cors$Q)], method = "fdr")
    } else {
      # keep only ex vs mb
      for(d in c("r", "n", "P")){
        my.cors[[d]] <- my.cors[[d]][colnames(my.ex), colnames(my.mb)]
      }
      # adjust p values
      my.cors$Q <- my.cors$P
      my.cors$Q[,] <- p.adjust(my.cors$P, method = "fdr")
    }
    out[[tp]] <- my.cors
  }
  out
}

# function to get taxon names from phyloseq tax table
get_taxon_names <- function(tax_table_obj, inc.asv = FALSE) {
  # Convert to data.frame if needed
  tt <- as.data.frame(tax_table_obj)
  rank_names <- colnames(tt)
  # Detect genus/species columns (case-insensitive)
  genus_col   <- grep("^Genus$", rank_names, ignore.case = TRUE, value = TRUE)
  species_col <- grep("^Species$", rank_names, ignore.case = TRUE, value = TRUE)
  genus   <- if (length(genus_col)) tt[[genus_col]] else rep(NA, nrow(tt))
  species <- if (length(species_col)) tt[[species_col]] else rep(NA, nrow(tt))
  # Build names: prefer Genus species → Genus sp. → next available rank + " sp."
  base_name <- ifelse(!is.na(genus) & genus != "" & !is.na(species) & species != "",
                      paste(genus, species),
                      ifelse(!is.na(genus) & genus != "",
                             paste0(genus, " sp."),
                             apply(tt, 1, function(x) {
                               # find most specific available name
                               last_non_na <- tail(na.omit(x[x != ""]), 1)
                               ifelse(length(last_non_na) == 0, NA, paste0(last_non_na, " sp."))
                             })
                      ))
  
  # Optionally append ASV / row name
  if (inc.asv) {
    asv_names <- sub("^[BF]_","", rownames(tt))
    base_name <- ifelse(!is.na(base_name) & base_name != "",
                        paste0(base_name, " (", asv_names, ")"),
                        asv_names)
  }
  return(base_name)
}

# function to get gene names
get_gene_names <- function(gene_list, gene_info) {
  # Check required columns exist
  required_cols <- c("gene", "Preferred_name", "Description")
  if (!all(required_cols %in% colnames(gene_info))) {
    stop("gene_info must contain columns: 'gene', 'Preferred_name', and 'Description'")
  }
  sapply(gene_list, function(g) {
    row <- gene_info[gene_info$gene == g, ]
    if (nrow(row) == 0) {
      return(g)  # gene not found
    }
    pref <- row[["Preferred_name"]][1]
    desc <- row[["Description"]][1]
    if (!is.na(pref) && pref != "-") {
      return(pref)
    } else if (!is.na(desc) && desc != "-") {
      return(desc)
    } else {
      return(g)
    }
  }, USE.NAMES = TRUE)
}

# function to filter df by significance and abs correlation 
filter_by_significance <- function(df, alpha = 0.05, 
                                   sig_filter = c("none", "TP1", "TP2", "either", "both"),
                                   min_abs_cor = 0) {
  sig_filter <- match.arg(sig_filter)
  if (sig_filter == "none") return(df)
  # per-timepoint logicals for sig + abs cor
  sig_TP1 <- (df$TP1_Q < alpha) & (abs(df$TP1_r) >= min_abs_cor)
  sig_TP2 <- (df$TP2_Q < alpha) & (abs(df$TP2_r) >= min_abs_cor)
  keep <- switch(sig_filter,
                 "TP1" = sig_TP1,
                 "TP2" = sig_TP2,
                 "either" = sig_TP1 | sig_TP2,
                 "both" = sig_TP1 & sig_TP2)
  df[keep, , drop = FALSE]
}

# function to convert cors output by allcors_ex_mb to df
cors2df <- function(cors, gene_info, taxonomy, 
                    only_mb_ex = TRUE, 
                    TPs = c("TP1", "TP2"), 
                    symmetric = TRUE, 
                    inc.asv = TRUE,
                    sig_filter = c("none", "TP1", "TP2", "either", "both"),
                    alpha = 0.05, 
                    min_abs_cor = 0,
                    show_progress = TRUE) {
  
  sig_filter <- match.arg(sig_filter)
  # helper function: check if a feature is a gene or microbe
  get_type <- function(f, genes, microbes) {
    if (f %in% genes) return("gene")
    # detect fungi/bacteria by prefix
    if (startsWith(f, "F_")) return("fungi")
    if (startsWith(f, "B_")) return("bacteria")
    # fallback if not found
    if (f %in% microbes) return("microbe")
    return(NA)
  }
  # extract gene and microbe names
  genes <- gene_info$gene
  microbes <- rownames(taxonomy)
  # extract correlation matrices for both timepoints
  r1 <- cors[[TPs[1]]]$r
  q1 <- cors[[TPs[1]]]$Q
  r2 <- cors[[TPs[2]]]$r
  q2 <- cors[[TPs[2]]]$Q
  # dimensions
  n <- nrow(r1)
  p <- ncol(r1)
  # set iteration ranges
  if (symmetric) {
    i_range <- 1:(n - 1)
    j_fun <- function(i) (i + 1):n
    ncors <- n * (n - 1) / 2
  } else {
    i_range <- 1:n
    j_fun <- function(i) 1:p
    ncors <- n * p
  }
  my.df <- data.frame(matrix(NA, 
                             ncol = 11, 
                             nrow = ncors,
                             dimnames = list(NULL, 
                                             c("feature1", "feature2", "cor_type", 
                                               "feature1_type", "feature2_type", 
                                               "feature1_name", "feature2_name", 
                                               "TP1_r", "TP1_Q", "TP2_r", "TP2_Q")
                             )))
  # initialize progress bar
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = length(i_range), style = 3)
  }
  cnt <- 1
  for (ii in seq_along(i_range)) {
    i <- i_range[ii]
    for (j in j_fun(i)) {
      f1 <- rownames(r1)[i]
      f2 <- colnames(r1)[j]
      if (symmetric && f1 == f2) next
      # get feature names
      my.df[cnt, "feature1"] <- f1
      my.df[cnt, "feature2"] <- f2
      # get cors
      my.df[cnt, "TP1_r"] <- r1[i, j]
      my.df[cnt, "TP1_Q"] <- q1[i, j]
      my.df[cnt, "TP2_r"] <- r2[i, j]
      my.df[cnt, "TP2_Q"] <- q2[i, j]
      # get feature type
      f1_type <- get_type(f1, genes, microbes)
      f2_type <- get_type(f2, genes, microbes)
      my.df[cnt, "feature1_type"] <- f1_type
      my.df[cnt, "feature2_type"] <- f2_type
      # get cor type
      if (f1_type == "gene" && f2_type == "gene") {
        cor_type <- "gene-gene"
      } else if (f1_type %in% c("fungi","bacteria") && f2_type %in% c("fungi","bacteria")) {
        cor_type <- "microbe-microbe"
      } else {
        cor_type <- "gene-microbe"
      }
      my.df[cnt, "cor_type"] <- cor_type
      # get feature 1 name
      if (f1_type == "gene") {
        f1_name <- get_gene_names(f1, gene_info)
      } else {
        f1_name <- get_taxon_names(taxonomy[f1, , drop = FALSE], inc.asv = inc.asv)
      }
      # get feature 2 name
      if (f2_type == "gene") {
        f2_name <- get_gene_names(f2, gene_info)
      } else {
        f2_name <- get_taxon_names(taxonomy[f2, , drop = FALSE], inc.asv = inc.asv)
      }
      my.df[cnt, "feature1_name"] <- f1_name
      my.df[cnt, "feature2_name"] <- f2_name
      # increment counter
      cnt <- cnt + 1
    }
    if (show_progress) setTxtProgressBar(pb, ii)
  }
  if (show_progress) close(pb)
  # remove empty rows
  my.df <- my.df[complete.cases(my.df$feature1), ]
  # optionally keep only gene–microbe correlations
  if (only_mb_ex) {
    my.df <- my.df[my.df$cor_type == "gene-microbe", ]
  }
  # filtering by significance
  my.df <- filter_by_significance(df = my.df, alpha = alpha, min_abs_cor = min_abs_cor, sig_filter = sig_filter)
  # remove NA correlations (optional cleanup)
  my.df <- my.df[!is.na(my.df$TP1_r) & !is.na(my.df$TP2_r), ]
  return(my.df)
}

# function to compute correlations between the requested feature1 vs feature2 per group per TP
compute_group_cors <- function(CCAdat, feature_pairs, TPs = c("TP1","TP2")) {
  # feature_pairs: data.frame with columns feature1 and feature2
  out <- list()
  for(tp in TPs) {
    ex_all <- CCAdat[[tp]]$ex_res    # samples x genes
    mb_all <- CCAdat[[tp]]$mb_res    # samples x microbes
    trt <- CCAdat[[tp]]$trt
    groups <- sort(unique(trt))
    out[[tp]] <- list()
    # get features
    f1_list <- as.character(feature_pairs$feature1)
    f2_list <- as.character(feature_pairs$feature2)
    for(g in groups) {
      # subset rows for this group
      idx <- which(trt == g)
      ex_grp <- ex_all[idx, , drop = FALSE]
      mb_grp <- mb_all[idx, , drop = FALSE]
      # allocate result vectors
      R <- numeric(length(f1_list))
      N <- integer(length(f1_list))
      P <- numeric(length(f1_list))
      # loop pairs
      for(i in seq_along(f1_list)) {
        f1 <- f1_list[i]; f2 <- f2_list[i]
        # check features present
        x <- ex_grp[, f1]
        y <- mb_grp[, f2]
        # count non-missing paired observations
        keep <- !is.na(x) & !is.na(y)
        n_obs <- sum(keep)
        N[i] <- n_obs
        if (n_obs > 2) {
          # Pearson correlation and p-value
          test <- suppressWarnings(cor.test(x[keep], y[keep], method = "pearson"))
          R[i] <- test$estimate
          P[i] <- test$p.value
        } else {
          R[i] <- NA
          P[i] <- NA
        }
      }
      # FDR correction per group
      P_adj <- p.adjust(P, method = "fdr")
      dfg <- data.frame(feature1 = f1_list,
                        feature2 = f2_list,
                        r = R,
                        n = N,
                        p = P,
                        p_adj = P_adj,
                        stringsAsFactors = FALSE)
      out[[tp]][[g]] <- dfg
    }
  }
  return(out)
}

# function to run Fisher's Z tests between each group to identify correlations that significantly differ between groups
run_pairwise_cocor_tests <- function(group_corrs, alpha = 0.05) {
  TPs <- names(group_corrs)
  out <- list()
  for (tp in TPs) {
    cors_by_group <- group_corrs[[tp]]
    groups <- names(cors_by_group)
    # assume identical feature1/2 in all groups
    feature_pairs <- cors_by_group[[1]][, c("feature1", "feature2")]
    df <- data.frame()
    # all pairwise group comparisons
    combs <- combn(groups, 2, simplify = FALSE)
    for (k in seq_len(nrow(feature_pairs))) {
      f1 <- feature_pairs$feature1[k]
      f2 <- feature_pairs$feature2[k]
      pairwise_res <- data.frame()
      for (c in combs) {
        g1 <- c[1]
        g2 <- c[2]
        r1 <- cors_by_group[[g1]]$r[k]
        r2 <- cors_by_group[[g2]]$r[k]
        n1 <- cors_by_group[[g1]]$n[k]
        n2 <- cors_by_group[[g2]]$n[k]
        # Skip invalid r values
        if (is.na(r1) || is.na(r2) || abs(r1) >= 1 || abs(r2) >= 1) next
        # Run independent group test using cocor
        cocor_res <- cocor::cocor.indep.groups(
          r1.jk = r1,
          r2.hm = r2,
          n1 = n1,
          n2 = n2,
          alternative = "two.sided",
          test = "fisher1925"
        )
        p_val <- cocor_res@fisher1925$p.value
        z_val <- cocor_res@fisher1925$statistic
        pairwise_res <- rbind(pairwise_res, data.frame(
          feature1 = f1,
          feature2 = f2,
          group1 = g1,
          group2 = g2,
          z_diff = z_val,
          p_value = p_val,
          stringsAsFactors = FALSE
        ))
      }
      df <- rbind(df, pairwise_res)
    }
    # FDR correction per TP
    df$p_adj <- p.adjust(df$p_value, method = "fdr")
    out[[tp]] <- df
  }
  return(out)
}

# function to test whether the sign of gene-microbe correlations is significantly associated between groups with Fisher's exact tests
sign_concordance_pairwise <- function(group_corrs) {
  # Flatten groups into "Group_TP" format
  all_groups <- list()
  for (tp in names(group_corrs)) {
    for (grp in names(group_corrs[[tp]])) {
      nm <- paste0(grp, "_", tp)   # <-- Group_TP
      all_groups[[nm]] <- group_corrs[[tp]][[grp]]
    }
  }
  group_names <- names(all_groups)
  K <- length(group_names)
  rows <- list()
  get_sign <- function(rvec) {
    ifelse(rvec > 0, "pos", ifelse(rvec < 0, "neg", NA))
  }
  for(a in 1:(K-1)) {
    for(b in (a+1):K) {
      g1 <- group_names[a]
      g2 <- group_names[b]
      d1 <- all_groups[[g1]][, c("feature1","feature2","r")]
      d2 <- all_groups[[g2]][, c("feature1","feature2","r")]
      colnames(d1)[3] <- "r1"
      colnames(d2)[3] <- "r2"
      merged <- merge(d1, d2, by = c("feature1","feature2"))
      s1 <- get_sign(merged$r1)
      s2 <- get_sign(merged$r2)
      tab <- table(s1, s2)
      pp <- tab["pos","pos"]
      nn <- tab["neg","neg"]
      pn <- tab["pos","neg"]
      np <- tab["neg","pos"]
      concord <- (pp + nn) / sum(pp, nn, pn, np)
      ft <- fisher.test(tab)
      rows[[length(rows)+1]] <- data.frame(
        group1 = g1, group2 = g2, n = nrow(merged),
        pos_pos = pp, neg_neg = nn,
        pos_neg = pn, neg_pos = np,
        concordance = concord,
        fisher_stat = as.numeric(ft$estimate),
        fisher_p = ft$p.value,
        fisher_p_adj = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  res <- do.call(rbind, rows)
  res$fisher_p_adj <- p.adjust(res$fisher_p, method = "fdr")
  return(res)
}

# function to test whether the value of gene-microbe correlations is significantly correlated between groups with permutation-based correlation tests
permtest_group_similarity <- function(group_corrs, B = 9999,
                                      method = "Pearson") {
  # combine TP1 and TP2 groups into a single named list
  all_groups <- list()
  for (tp in names(group_corrs)) {
    for (g in names(group_corrs[[tp]])) {
      nm <- paste(g, tp, sep = "_")
      all_groups[[nm]] <- group_corrs[[tp]][[g]]
    }
  }
  group_names <- names(all_groups)
  K <- length(all_groups)
  out <- data.frame(
    group1 = character(),
    group2 = character(),
    r_obs = numeric(),
    p_perm = numeric(),
    stringsAsFactors = FALSE
  )
  for (i in 1:(K-1)) {
    for (j in (i+1):K) {
      g1 <- group_names[i]
      g2 <- group_names[j]
      r1 <- all_groups[[g1]]$r
      r2 <- all_groups[[g2]]$r
      # remove any NA pairs before permutation test
      keep <- complete.cases(r1, r2)
      r1 <- r1[keep]
      r2 <- r2[keep]
      # permutation test
      pc <- PermCor::perm_test(
        x = r1, 
        y = r2, 
        B = B,
        method = method,
        alternative = "two.sided"
      )
      out <- rbind(out, data.frame(
        group1 = g1,
        group2 = g2,
        r_obs = pc$estimate,
        p_perm = pc$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
  out$p_adj <- p.adjust(out$p_perm, method = "fdr")
  return(out)
}

# function to plot treatment-treatment correlation of gene-microbe correlations
plot_group_cors_matrix <- function(group_corrs, perm_results = NULL, alpha = 0.05, aes) {
  # Flatten groups with TP prefix
  all_groups <- list()
  for(tp in names(group_corrs)) {
    for(g in names(group_corrs[[tp]])) {
      nm <- paste(g, tp, sep = "_")
      all_groups[[nm]] <- group_corrs[[tp]][[g]]
    }
  }
  group_names <- names(all_groups)
  K <- length(group_names)
  clean_names <- gsub("_", " ", group_names)  # for diagonal labels
  # Helper to get permutation results for a pair
  get_perm_cell <- function(g1, g2) {
    if(is.null(perm_results)) return(list(r_val=NA, p_val=NA, label="NA"))
    idx <- which((perm_results$group1 == g1 & perm_results$group2 == g2) |
                   (perm_results$group1 == g2 & perm_results$group2 == g1))[1]
    r_val <- perm_results$r_obs[idx]
    p_val <- perm_results$p_adj[idx]
    rho_sign <- perm_results$concordance[idx]
    # Format P-value + r
    if(is.na(p_val)) lbl <- "NA"
    else if(p_val < 0.001) lbl <- paste0("r = ", round(r_val,2), "\nP <0.001")
    else if(p_val >= alpha) lbl <- paste0("r = ", round(r_val,2), "\nP = n.s")
    else lbl <- paste0("r = ", round(r_val,2), "\nP = ", round(p_val,3))
    # Add SC (sign concordance) to label if available
    if(!is.na(rho_sign)) lbl <- paste0(lbl, "\nSC = ", round(rho_sign,2))
    return(list(r_val=r_val, p_val=p_val, label=lbl))
  }
  plot_list <- list()
  for(i in 1:K) {
    for(j in 1:K) {
      g1 <- group_names[i]
      g2 <- group_names[j]
      if(i == j) {
        # Diagonal: density + group label
        r_vec <- all_groups[[g1]]$r
        r_vec <- r_vec[!is.na(r_vec)]
        p <- ggplot2::ggplot(data.frame(r=r_vec), ggplot2::aes(x=r)) +
          ggplot2::geom_density(fill=aes[g1,"col"], alpha=0.5) +
          ggplot2::xlim(-1,1) +
          ggplot2::theme_minimal() +
          ggplot2::xlab(clean_names[i]) +
          ggplot2::ylab("") +
          ggplot2::theme(axis.text.y=ggplot2::element_blank(),
                axis.ticks.y=ggplot2::element_blank(),
                axis.text.x=ggplot2::element_text(size=6),
                axis.title.x=ggplot2::element_text(size=8),
                plot.margin=ggplot2::margin(1,1,1,1))
      } else if(i > j) {
        # Lower triangle: scatter plot with regression line, no axis titles
        df1 <- all_groups[[g1]]
        df2 <- all_groups[[g2]]
        merged <- merge(df1, df2, by=c("feature1","feature2"), suffixes=c("_1","_2"))
        merged <- merged[complete.cases(merged[,c("r_1","r_2")]), ]
        show_x <- i == K  # bottom row
        show_y <- j == 1  # leftmost column
        p <- ggplot2::ggplot(merged, aes(x=r_1, y=r_2)) +
          ggplot2::geom_point(alpha=0.5, size=0.8) +
          ggplot2::geom_smooth(method="lm", se=FALSE, color="darkred") +
          ggplot2::xlim(-1,1) + 
          ggplot2::ylim(-1,1) +
          ggplot2::theme_minimal() +
          ggplot2::labs(x="r1", y="r2") +   
          ggplot2::theme(plot.margin=margin(1,1,1,1),
                axis.title.x=if(show_x) ggplot2::element_text(size=7) else ggplot2::element_blank(),
                axis.text.x=if(show_x) ggplot2::element_text(size=6) else ggplot2::element_blank(),
                axis.title.y=if(show_y) ggplot2::element_text(size=7) else ggplot2::element_blank(),
                axis.text.y=if(show_y) ggplot2::element_text(size=6) else ggplot2::element_blank())
      } else {
        # Upper triangle: heatmap with r_obs and formatted P
        cell <- get_perm_cell(g1, g2)
        p <- ggplot2::ggplot(data.frame(x=1, y=1, fill=cell$r_val, label=cell$label),
                             ggplot2::aes(x=x, y=y, fill=fill)) +
          ggplot2::geom_tile(color="white") +
          ggplot2::geom_text(ggplot2::aes(label=label), size=3) +
          ggplot2::scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, limits=c(-1,1)) +
          ggplot2::theme_void() +
          ggplot2::theme(legend.position="none")
      }
      plot_list[[paste(i,j)]] <- p
    }
  }
  final_plot <- patchwork::wrap_plots(plot_list, ncol=K)
  return(final_plot)
}
