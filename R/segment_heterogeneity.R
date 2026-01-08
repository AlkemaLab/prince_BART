# This is required for R CMD check to suppress visible binding notes for NSE variables
if (getRversion() >= "2.15.1") utils::globalVariables(c("segment_label", "is_overall", "ci_lower", "ci_upper", "d"))
#' Effect Heterogeneity by Segments
#'
#' Partition the sample into covariate-defined segments and summarize how
#' estimated complier effects vary across segments, following the heterogeneity
#' analysis in the Prince BART paper.
#'
#' @param princebart A fitted object from \code{prince_BART}.
#' @param data Optional data.frame of covariates used for segmentation. If NULL,
#'   uses \code{princebart$data$X}.
#' @param vars Optional character vector of covariate names to include in the
#'   segmentation model. If NULL, uses all covariates in \code{data}.
#' @param rpart_control An \code{rpart.control} object for tree fitting.
#' @param min_compliers_bucket Minimum effective complier count for pruning
#'   (sum of posterior mean complier probabilities within a segment). Default is 100.
#' @param ci_level Credible interval level (e.g., 0.9 for 90% credible intervals).
#'   Default is 0.9.
#' @param plot Logical; if TRUE, return ggplot2 objects for segment effects
#'   (and the max--min contrast plot if \code{contrast = TRUE}).
#' @param contrast Logical; if TRUE, compute posterior draws for the difference
#'   between the segments with the highest and lowest estimated mean effects.
#'
#' @details
#' This function fits a shallow \code{rpart} regression tree to posterior mean
#' conditional complier effects \eqn{\mathrm{CATE}_C(x)} and assigns each unit to
#' a terminal-node segment. Tree fitting is weighted by each unit's posterior mean
#' complier probability, so that splits emphasize regions of the covariate space
#' with more compliers.
#'
#' For each segment, the function aggregates posterior draws of
#' \eqn{\mathrm{CATE}_C(x)} using complier-probability weights to obtain a
#' segment-specific average effect among compliers (an \emph{MCATE_C} estimand in
#' the paper). These segment-level summaries are mixed (sample-based) estimands:
#' they average conditional effects over the empirical covariate distribution within
#' each segment rather than defining new population-level causal parameters.
#'
#' The main segment summaries are returned in \code{res$effects} and visualized in
#' \code{res$plot$effect} (if \code{plot = TRUE}).
#'
#' Optionally, segments with small effective numbers of compliers can be pruned/merged
#' for stability. If \code{contrast = TRUE}, heterogeneity is summarized by the
#' posterior distribution of the difference between the segments with the largest and
#' smallest estimated mean effects. Numerical results are in \code{res$contrast$summary}
#' and the histogram is in \code{res$plot$diff} (if \code{plot = TRUE}).
#'
#' @return A list with components:
#' \item{tree}{Fitted (and optionally pruned) rpart tree.}
#' \item{segment_data}{Input data with \code{cate} (posterior mean \eqn{\mathrm{CATE}_C(x)}),
#'   \code{w} (posterior mean complier probability), and \code{segment}.}
#' \item{effects}{Segment-level effect summaries (MCATE_C): posterior mean, sd, credible
#'   interval bounds, and \code{p_gt0}, plus \code{n} (segment size) and \code{n_complier}
#'   (estimated complier count).}
#' \item{draws}{List of posterior draw matrices by segment.}
#' \item{contrast}{If \code{contrast = TRUE}, posterior comparison of the highest- vs
#'   lowest-effect segments, containing \code{$summary} and posterior draws in \code{$draws}.}
#' \item{plot}{If \code{plot = TRUE}, a list of ggplot objects including \code{$effect}
#'   (segment effects) and, when \code{contrast = TRUE}, \code{$diff} (difference distribution).}
#'
#' @export

segment_heterogeneity <- function(
  princebart,
  data = NULL,
  vars = NULL,
  rpart_control = rpart::rpart.control(cp = 0, maxdepth = 3L),
  min_compliers_bucket = 100,
  ci_level = 0.9,
  plot = FALSE,
  contrast = FALSE
) {
  if (!inherits(princebart, "princebart")) {
    stop("princebart must be a fitted princebart object")
  }
  if (is.null(princebart$probs) || length(dim(princebart$probs)) != 4) {
    stop("princebart$probs must be a 4D array")
  }

  if (is.null(data)) {
    if (is.null(princebart$data$X)) {
      stop("data must be provided when princebart$data$X is missing")
    }
    data <- princebart$data$X
  }
  data <- as.data.frame(data)

  n_units <- dim(princebart$probs)[4]
  if (nrow(data) != n_units) {
    stop("data and princebart$probs must have the same number of units")
  }

  if (is.null(vars)) {
    vars <- setdiff(names(data), c("cate", "w", "segment", "e"))
  }
  if (length(vars) == 0) {
    stop("vars must contain at least one covariate name")
  }
  if (!all(vars %in% names(data))) {
    stop("vars contains names not found in data")
  }

  alpha <- (1 - ci_level) / 2

  prob <- princebart$probs
  p_n <- prob[, , "p_n", ]
  p_a <- prob[, , "p_a", ]
  p_c <- 1 - p_n - p_a
  my1c <- prob[, , "m_y1c", ]
  my0c <- prob[, , "m_y0c", ]
  cate_draws <- my1c - my0c

  data$cate <- apply(cate_draws, 3, mean, na.rm = TRUE)
  data$w <- apply(p_c, 3, mean, na.rm = TRUE)

  rpart_data <- data[, unique(c("cate", vars, "w")), drop = FALSE]
  seg_formula <- paste("cate ~", paste(vars, collapse = " + ")) |>
    stats::as.formula()

  tree <- rpart::rpart(
    seg_formula,
    data = rpart_data,
    weights = rpart_data$w,
    control = rpart_control
  )

  if (min_compliers_bucket > 0) {
    tree <- prune_weighted_tree(rpart_data, tree, min_compliers_bucket)
  }

  data$segment <- assign_segments(tree)

  segment_draws <- function(include) {
    if (!any(include)) {
      return(matrix(NA_real_, nrow = dim(prob)[1], ncol = dim(prob)[2]))
    }
    numer <- apply(
      cate_draws[, , include, drop = FALSE] * p_c[, , include, drop = FALSE],
      1:2,
      sum
    )
    denom <- apply(p_c[, , include, drop = FALSE], 1:2, sum)
    est <- numer / denom
    est[denom == 0] <- NA_real_
    est
  }

  summarize_draws <- function(draws) {
    v <- as.vector(draws)
    v <- v[is.finite(v)]
    if (length(v) == 0) {
      return(list(mean = NA_real_, sd = NA_real_, ci90 = c(NA_real_, NA_real_),
                  p_gt0 = NA_real_))
    }
    list(
      mean = mean(v),
      sd = stats::sd(v),
      ci = stats::quantile(v, c(alpha, 1 - alpha), names = FALSE),
      p_gt0 = mean(v > 0)
    )
  }

  effects <- list()
  draws <- list()

  add_segment_result <- function(seg_name, include) {
    est <- segment_draws(include)
    sumry <- summarize_draws(est)
    effects[[seg_name]] <<- data.frame(
      segment = seg_name,
      mean = sumry$mean,
      sd = sumry$sd,
      ci_lower = sumry$ci[1],
      ci_upper = sumry$ci[2],
      p_gt0 = sumry$p_gt0,
      n = sum(include),
      n_complier = sum(data$w[include]),
      stringsAsFactors = FALSE
    )
    draws[[seg_name]] <<- est
  }

  add_segment_result("overall", rep(TRUE, nrow(data)))
  for (seg in levels(data$segment)) {
    add_segment_result(seg, data$segment == seg)
  }

  effects <- do.call(rbind, effects)
  class(effects) <- c("princebart_segment_effects", class(effects))

  contrast_res <- NULL
  if (contrast) {
    seg_means <- tapply(data$cate, data$segment, mean, na.rm = TRUE)
    seg_min <- names(which.min(seg_means))
    seg_max <- names(which.max(seg_means))
    est_min <- segment_draws(data$segment == seg_min)
    est_max <- segment_draws(data$segment == seg_max)
    diff_draws <- as.vector(est_max - est_min)
    diff_draws <- diff_draws[is.finite(diff_draws)]
    if (length(diff_draws) == 0) {
      diff_draws <- NA_real_
    }
    contrast_res <- list(
      segment_low = seg_min,
      segment_high = seg_max,
      draws = diff_draws,
      summary = data.frame(
        mean = mean(diff_draws, na.rm = TRUE),
        sd = stats::sd(diff_draws, na.rm = TRUE),
        ci_lower = stats::quantile(diff_draws, alpha, na.rm = TRUE),
        ci_upper = stats::quantile(diff_draws, 1 - alpha, na.rm = TRUE),
        p_gt0 = mean(diff_draws > 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    )
  }
  wrap_label <- function(x, width = 40) {
    vapply(
      x,
      function(s) paste(strwrap(s, width = width), collapse = "\n"),
      character(1)
    )
  }

  plot_res <- NULL
  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("ggplot2 is required for plot = TRUE")
    }
    plot_df <- effects
    plot_df$segment_label <- as.character(plot_df$segment)
    plot_df$is_overall <- plot_df$segment == "overall"

    effect_plot <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        y = stats::reorder(wrap_label(segment_label, 40), mean),
        x = mean,
        color = !is_overall
      )
    ) +
      ggplot2::geom_vline(xintercept = 0, color = "grey") +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
        width = 0.5
      ) +
      ggplot2::geom_point() +
      ggplot2::labs(
        x = paste("Subgroup-specific effect size estimate ("
          , ci_level * 100, "% CI)", sep = ""
        ),
        y = "",
        title = ""
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")

    plot_res <- list(effect = effect_plot)

    if (contrast && !is.null(contrast_res)) {
      diff_df <- data.frame(d = contrast_res$draws)
      diff_plot <- ggplot2::ggplot(diff_df, ggplot2::aes(x = d)) +
        ggplot2::geom_histogram(bins = 30, fill = "grey70", color = "white") +
        ggplot2::geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
        ggplot2::labs(
          x = "Difference in effect across subgroups",
          y = "Count",
          title = ""
        ) +
        ggplot2::theme_minimal()
      plot_res$diff <- diff_plot
    }
  }

  out <- list(
    tree = tree,
    segment_data = data,
    effects = effects,
    draws = draws,
    contrast = contrast_res,
    plot = plot_res
  )
  class(out) <- c("princebart_segment", "list")
  out
}

assign_segments <- function(tree) {
  party_tree <- partykit::as.party(tree)
  list_rules <- utils::getFromNamespace(".list.rules.party", "partykit")
  rules <- list_rules(party_tree)
  rules <- simplify_and_pretty_rules(rules)
  nodes <- stats::predict(party_tree, type = "node")
  factor(nodes, labels = rules)
}

#' @export
print.princebart_segment_effects <- function(x, digits = 3, ...) {
  print.data.frame(x, digits = digits, row.names = FALSE, ...)
  invisible(x)
}

prune_weighted_tree <- function(data, tree, min_weighted_bucket = 100) {
  if (is.null(tree$cptable) || nrow(tree$cptable) == 0) {
    return(tree)
  }
  cp_seq <- sort(unique(tree$cptable[, "CP"]))
  for (cp in cp_seq) {
    tree <- rpart::prune(tree, cp = cp)
    seg_weights <- weighted_segment_sizes(data, tree)
    if (all(seg_weights > min_weighted_bucket, na.rm = TRUE)) {
      break
    }
  }
  tree
}

weighted_segment_sizes <- function(data, tree) {
  seg <- assign_segments(tree)
  tapply(data$w, seg, sum)
}

simplify_and_pretty_rules <- function(rules, digits = 1) {

  pretty_num <- function(x) {
    if (is.na(x)) return(NA_character_)
    if (abs(x) >= 1000) return(paste0(round(x/1000, digits), "K"))
    formatC(x, format = "f", digits = digits)
  }

  simplify_one <- function(rule) {
    parts <- strsplit(rule, " & ", fixed = TRUE)[[1]]

    # Parse: var op value  (value numeric)
    rx <- "^\\s*([A-Za-z.][A-Za-z0-9._]*)\\s*(<=|>=|<|>)\\s*([-+]?[0-9]*\\.?[0-9]+)\\s*$"
    m <- regexec(rx, parts, perl = TRUE)
    g <- regmatches(parts, m)

    ok  <- lengths(g) == 4
    raw <- parts[!ok]
    g   <- g[ok]
    if (length(g) == 0) return(rule)

    var <- vapply(g, "[[", "", 2)
    op  <- vapply(g, "[[", "", 3)
    val <- as.numeric(vapply(g, "[[", "", 4))

    out <- character(0)
    for (v in unique(var)) {
      idx <- var == v

      # strongest lower bound: keep the maximum among >= or >
      lb <- idx & op %in% c(">=", ">")
      if (any(lb)) {
        best <- max(val[lb])
        out <- c(out, paste(v, ">=", pretty_num(best)))
      }

      # strongest upper bound: keep the minimum among < or <=
      ub <- idx & op %in% c("<", "<=")
      if (any(ub)) {
        best <- min(val[ub])
        out <- c(out, paste(v, "<", pretty_num(best)))
      }
    }

    paste(c(out, raw), collapse = " & ")
  }

  vapply(rules, simplify_one, character(1))
}
