# ggplot2 autoplot for krls objects. Discoverable via library(ggplot2);
# ggplot2 is in Suggests:.

utils::globalVariables(c(".data"))

autoplot.krls <-
function(object, ...)
  {
    if (!requireNamespace("ggplot2", quietly = TRUE))
      stop("ggplot2 is required for autoplot(). install.packages(\"ggplot2\")")
    if (!inherits(object, "krls"))
      stop("object is not of class 'krls'")
    if (is.null(object$derivatives))
      stop("\n autoplot() requires marginal effects; refit with krls(..., derivative = TRUE)\n")

    derivs <- as.matrix(object$derivatives)
    colnames(derivs) <- colnames(object$X)
    long <- do.call(rbind, lapply(colnames(derivs), function(j)
      data.frame(term = j, dy_dx = derivs[, j], stringsAsFactors = FALSE)))
    long$term <- factor(long$term, levels = colnames(derivs))

    ame <- data.frame(
      term     = factor(colnames(object$X), levels = colnames(object$X)),
      estimate = as.numeric(object$avgderivatives),
      stringsAsFactors = FALSE
    )

    ggplot2::ggplot(long, ggplot2::aes(x = .data$dy_dx)) +
      ggplot2::geom_histogram(bins = 30, fill = "grey80", color = "grey40") +
      ggplot2::geom_vline(data = ame,
                          ggplot2::aes(xintercept = .data$estimate),
                          color = "darkblue", linewidth = 0.7) +
      ggplot2::geom_vline(xintercept = 0, color = "grey60",
                          linetype = "dashed") +
      ggplot2::facet_wrap(~ term, scales = "free") +
      ggplot2::labs(x = expression(partialdiff * hat(y) / partialdiff * x),
                    y = "Number of observations",
                    title = "Pointwise marginal effects",
                    subtitle = "Blue line: average marginal effect (AME). Dashed line at 0 for reference.") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9,
                                                            color = "grey30"))
  }
