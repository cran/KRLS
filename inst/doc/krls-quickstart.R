## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5,
  dpi = 96
)
set.seed(1)

## ----setup--------------------------------------------------------------------
library(KRLS)

## ----data---------------------------------------------------------------------
set.seed(1)
n  <- 200
x1 <- runif(n, -2, 2)
x2 <- runif(n, -2, 2)
# True surface: nonlinear in x1, modulated by x2.
y  <- sin(x1) + 0.5 * x2 * (x1 > 0) + rnorm(n, sd = 0.3)
df <- data.frame(y = y, x1 = x1, x2 = x2)

## ----fit----------------------------------------------------------------------
fit <- krls(y ~ x1 + x2, data = df, print.level = 0)

## ----tidy---------------------------------------------------------------------
library(generics)
tidy(fit)
glance(fit)

## ----autoplot, eval = requireNamespace("ggplot2", quietly = TRUE)-------------
library(ggplot2)
autoplot(fit)

## ----augment------------------------------------------------------------------
aug <- augment(fit, data = df)
head(aug, 3)

