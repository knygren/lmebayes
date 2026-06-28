#' Plot Block~2 sweep-history diagnostics (cross-chain mean or SD)
#'
#' Plots cross-chain Block~2 hyperparameter summaries stored on an object of
#' class \code{"two_block_sweep_history"} (typically
#' \code{fit$sweep_history$pilot} or \code{fit$sweep_history$main} from
#' \code{\link{lmerb}} / \code{\link{glmerb}}; see
#' \code{\link[glmbayesCore]{print.two_block_sweep_history}}).
#'
#' @param hist Object of class \code{"two_block_sweep_history"} (see
#'   \code{\link[glmbayesCore]{print.two_block_sweep_history}}).
#' @param coef_focus List of length-2 character vectors
#'   \code{c(re_component, covariate)}, matching rows of \code{hist$table}.
#'   Example: \code{list(c("(Intercept)", "(Intercept)"), c("violent_i", "(Intercept)"))}.
#' @param what One or both of \code{"sd"} and \code{"mean"} (cross-chain
#'   summary after each inner sweep). Default both (\code{sd} first).
#' @param engine \code{"base"} for one panel per coefficient (default), or
#'   \code{"ggplot"} for a single faceted figure (requires \pkg{ggplot2}).
#' @param stage_label Character label for titles; defaults to \code{hist$stage}.
#' @return \code{hist} invisibly.
#' @seealso \code{\link{lmerb}}, \code{\link{glmerb}},
#'   \code{\link[glmbayesCore]{print.two_block_sweep_history}}
#' @export
plot_sweep_history_diag <- function(
    hist,
    coef_focus,
    what = c("sd", "mean"),
    engine = c("base", "ggplot"),
    stage_label = hist$stage
) {
  if (!inherits(hist, "two_block_sweep_history")) {
    stop(
      "'hist' must be a two_block_sweep_history object ",
      "(e.g. fit$sweep_history$main).",
      call. = FALSE
    )
  }
  if (!is.list(coef_focus) || !length(coef_focus)) {
    stop("'coef_focus' must be a non-empty list of c(re_component, covariate) pairs.",
         call. = FALSE)
  }

  what <- match.arg(what, c("sd", "mean"), several.ok = TRUE)
  what <- what[order(match(what, c("sd", "mean")))]
  engine <- match.arg(engine)
  stage_label <- as.character(stage_label)[1L]
  if (!nzchar(stage_label)) {
    stage_label <- if (!is.null(hist$stage)) hist$stage else "stage"
  }

  if (identical(engine, "ggplot") && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("'engine = \"ggplot\"' requires the ggplot2 package.", call. = FALSE)
  }

  sh_tab <- hist$table
  sh_tab$re_component <- as.character(sh_tab$re_component)
  sh_tab$covariate <- as.character(sh_tab$covariate)
  sh_sweeps <- sh_tab[sh_tab$sweep > 0L, , drop = FALSE]
  if (!nrow(sh_sweeps)) {
    warning("No sweep rows in sweep history for stage ", stage_label, call. = FALSE)
    return(invisible(hist))
  }

  sh_plot <- do.call(rbind, lapply(coef_focus, function(cc) {
    if (length(cc) < 2L) {
      stop("Each element of 'coef_focus' must be c(re_component, covariate).",
           call. = FALSE)
    }
    sh_sweeps[
      sh_sweeps$re_component == as.character(cc[1L]) &
        sh_sweeps$covariate == as.character(cc[2L]),
      ,
      drop = FALSE
    ]
  }))
  rownames(sh_plot) <- NULL

  plot_one <- function(re_comp, cov, metric, ylab) {
    sub <- sh_sweeps[
      sh_sweeps$re_component == re_comp & sh_sweeps$covariate == cov,
      ,
      drop = FALSE
    ]
    sub <- sub[order(sub$sweep), , drop = FALSE]
    if (!nrow(sub)) {
      warning("No sweep rows for ", re_comp, " | ", cov, call. = FALSE)
      return(invisible(FALSE))
    }
    y <- if (metric == "sd") sub$sd else sub$mean
    if (!any(is.finite(y))) {
      warning("No finite ", metric, " values for ", re_comp, " | ", cov,
              call. = FALSE)
      return(invisible(FALSE))
    }
    graphics::plot(
      sub$sweep, y,
      type = "b", pch = 16,
      xlab = "Inner sweep", ylab = ylab,
      main = paste(re_comp, cov, sep = " | ")
    )
    if (metric == "mean") {
      mode_val <- sh_tab$mean[
        sh_tab$re_component == re_comp &
          sh_tab$covariate == cov &
          sh_tab$sweep == 0L
      ]
      if (length(mode_val) == 1L && is.finite(mode_val)) {
        graphics::abline(h = mode_val, lty = 2, col = "gray40")
      }
    }
    invisible(TRUE)
  }

  for (metric in what) {
    ylab <- if (metric == "sd") "Cross-chain SD" else "Cross-chain mean"
    cat(sprintf(
      "\n=== %s sweep history (%s; %s) ===\n\n",
      stage_label, ylab, engine
    ))

    if (identical(engine, "base")) {
      op <- graphics::par(
        mfrow = c(length(coef_focus), 1L),
        mar = c(4, 4, 2.5, 1),
        oma = c(0, 0, 2, 0)
      )
      n_plotted <- 0L
      for (cc in coef_focus) {
        if (isTRUE(plot_one(as.character(cc[1L]), as.character(cc[2L]), metric, ylab))) {
          n_plotted <- n_plotted + 1L
        }
      }
      graphics::mtext(
        sprintf("%s Block 2 fixef: cross-chain %s by inner sweep", stage_label, metric),
        outer = TRUE, line = 0.5, cex = 0.95
      )
      if (metric == "mean") {
        graphics::mtext(
          "Dashed line = ICM mode (sweep 0)",
          outer = TRUE, line = -1.5, cex = 0.85
        )
      }
      graphics::par(op)
      if (n_plotted == 0L) {
        warning(
          "No panels drawn for ", metric, ". Check coef_focus against ",
          "unique(hist$table[, c('re_component', 'covariate')]).",
          call. = FALSE
        )
      }
    } else if (nrow(sh_plot)) {
      sh_plot$coef <- interaction(
        sh_plot$re_component, sh_plot$covariate, sep = " | "
      )
      y_var <- if (metric == "sd") "sd" else "mean"
      aes_y <- switch(
        y_var,
        sd = sh_plot$sd,
        mean = sh_plot$mean
      )
      p <- ggplot2::ggplot(
        sh_plot,
        ggplot2::aes(sweep, aes_y, group = coef, colour = coef)
      ) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::facet_wrap(~ coef, scales = "free_y") +
        ggplot2::labs(
          x = "Inner sweep",
          y = ylab,
          title = sprintf(
            "%s Block 2 fixef - cross-chain %s by sweep",
            stage_label, metric
          )
        ) +
        ggplot2::theme(legend.position = "none")
      if (metric == "mean") {
        mode_df <- do.call(rbind, lapply(coef_focus, function(cc) {
          sh_tab[
            sh_tab$re_component == as.character(cc[1L]) &
              sh_tab$covariate == as.character(cc[2L]) &
              sh_tab$sweep == 0L,
            ,
            drop = FALSE
          ]
        }))
        if (nrow(mode_df)) {
          mode_df$coef <- interaction(
            mode_df$re_component, mode_df$covariate, sep = " | "
          )
          mode_df$yint <- mode_df$mean
          p <- p + ggplot2::geom_hline(
            ggplot2::aes(yintercept = yint, linetype = "ICM mode"),
            data = mode_df,
            colour = "gray40"
          ) +
            ggplot2::scale_linetype_manual(
              name = NULL, values = c("ICM mode" = "dashed")
            )
        }
      }
      print(p)
    }
  }

  invisible(hist)
}
