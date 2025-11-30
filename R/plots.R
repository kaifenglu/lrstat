#' @title Plot method for assess_phregr objects
#' @description Generate line plots of observed and simulated paths of
#' standardized score processes of a Cox model.
#'
#' @param x An object of class \code{assess_phregr}.
#' @param nsim The number of simulation samples used in the plot.
#' @param ... Ensures that all arguments starting from "..." are named.
#'
#' @return A list of ggplot2 objects for the line plots, one for each
#' covariate.
#'
#' @keywords internal
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @method plot assess_phregr
#' @export
plot.assess_phregr <- function(x, nsim = 20, ...) {
  if (!inherits(x, "assess_phregr")) {
    stop("x must be of class 'assess_phregr'")
  }

  format_pvalue <- function(p) {
    # Handle the case of p < 0.0001
    ifelse(p < 0.0001,
           "<.0001",
           # Handle the case of p > 0.9999
           ifelse(p > 0.9999,
                  ">.9999",
                  # For all other cases, format to 4 decimal places
                  sprintf("%.4f", p)))
  }

  g <- list()
  for (j in 1:length(x$covariates)) {
    g[[j]] <- ggplot2::ggplot(
      data = rbind(data.frame(t = x$time,
                              u = x$score_t[,j],
                              sim = 0,
                              type = "observed"),
                   data.frame(t = rep(x$time, nsim),
                              u = do.call(c, lapply(1:nsim, function(i)
                                x$score_t_list[[i]][,j])),
                              sim = rep(1:nsim, each=length(x$time)),
                              type = "simulated")),
      ggplot2::aes(x = .data$t, y = .data$u, group = .data$sim,
                   linetype =.data$type, linewidth = .data$type,
                   color = .data$type)) +
      ggplot2::scale_linewidth_manual(values = c("observed" = 0.75,
                                                 "simulated" = 0.5)) +
      ggplot2::geom_step() +
      ggplot2::annotate("text",
                        x = max(x$time) * 0.8, y = -Inf,
                        label = paste0("Pr > MaxAbsVal: ",
                                       format_pvalue(x$p_value[j]), "\n",
                                       "(", x$resample, " simulations)"),
                        hjust = 0, vjust = -0.5,
                        size = 3) +
      ggplot2::labs(
        x = "Time", y = "Standardized Score Process",
        title = paste("Checking Proportional Hazards Assumption for",
                      x$covariates[j]),
        subtitle = paste("Observed Path and First",
                         nsim, "Simulated Paths")) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
                     plot.title = ggplot2::element_text(hjust = 0.5,
                                                        face = "bold"),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5))
  }

  g
}
