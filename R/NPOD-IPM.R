#' @title IPM
#' @description
#' The IPM function solves a convex optimization problem using an interior-point Newton method.
#' It iteratively updates lambda, `lam`, while ensuring feasibility constraints.
#' A damped step-size approach is used to maintain numerical stability.
#' @param Psi A matrix of non-negative values.
#' @return A list with the objective function value and the optimal lambda.
#' @keywords internal
#' @importFrom Matrix Diagonal
IPM <- function(Psi) {
  Psi <- abs(Psi)

  # Check Input Specifications
  row <- nrow(Psi)
  col <- ncol(Psi)

  # Check that Psi has only non-negative entries.
  if (min(Psi) < 0) {
    stop("The matrix Psi has a negative element.")
  }

  # Check that Psi * ecol > 0.
  ecol <- matrix(1, nrow = col, ncol = 1) # Ensure ecol is a column vector
  Plam <- Psi %*% ecol
  if (min(Plam) <= 1.e-15) {
    stop("The vector Psi * ecol has a non-positive entry.")
  }

  # Stopping tolerance
  eps <- 1.e-8

  # Initialization
  sig <- 0
  erow <- rep(1, row)
  lam <- ecol
  w <- 1 / Plam
  Ptw <- t(Psi) %*% w
  shrink <- 2 * max(Ptw)
  lam <- lam * shrink
  Plam <- Plam * shrink
  w <- w / shrink
  Ptw <- Ptw / shrink
  y <- ecol - Ptw
  R <- erow - w * Plam
  normR <- max(abs(R))
  gap <- abs((sum(log(w)) + sum(log(Plam)))) / (1 + max(abs(sum(log(Plam)))))
  mu <- (t(lam) %*% y) / col
  mu <- mu[1]

  # Reporting structures
  normNiter <- c(normR)
  muiter <- c(mu)
  iter <- 0

  # Newton Iteration
  while (mu > eps || normR > eps || gap > eps) {
    iter <- iter + 1
    smu <- sig * mu

    # Form the matrix H
    inner <- lam / y
    H <- Psi %*% Matrix::Diagonal(x = inner, n = col) %*% t(Psi) + Matrix::Diagonal(x = Plam / w, n = row)
    H <- as.matrix(H)
    # Cholesky factorization of H
    cholH <- try(chol(H), silent = TRUE)
    if ("try-error" %in% class(cholH)) {
      lam <- lam / sum(lam)
      return(list(obj = sum(log(Psi %*% lam)), lam = lam))
    }
    # Right-hand side for dw
    smuyinv <- smu * (ecol / y)
    rhsdw <- erow / w - Psi %*% smuyinv

    # Compute dw
    dw <- solve(cholH, solve(t(cholH), rhsdw))

    # Compute dy
    dy <- -t(Psi) %*% dw

    # Compute dlam
    dlam <- smuyinv - lam - inner * dy

    # Damped Newton Step

    alfpri <- -1 / min(min(dlam / lam), -0.5)
    alfdual <- -1 / min(min(dy / y), -0.5)
    alfdual <- min(alfdual, -1 / min(min(dw / w), -0.5))
    alfpri <- min(1, 0.99995 * alfpri)
    alfdual <- min(1, 0.99995 * alfdual)

    # Newton Update
    if (alfpri == 1) {
      lam <- lam + dlam
    } else {
      lam <- lam + alfpri * dlam
    }


    if (alfdual == 1) {
      w <- w + dw
      y <- y + dy
    } else {
      w <- w + alfdual * dw
      y <- y + alfdual * dy
    }

    # Update mu and intermediary vectors and stopping criteria values
    mu <- (t(lam) %*% y) / col
    mu <- mu[1]
    Plam <- Psi %*% lam
    R <- erow - w * Plam
    Ptw <- Ptw - alfdual * dy
    normR <- max(abs(R))
    gap <- abs((sum(log(w)) + sum(log(Plam)))) / (1 + max(abs(sum(log(Plam)))))

    # Update sig
    if (mu < eps && normR > eps) {
      sig <- 1
    } else {
      c1 <- 2
      c2 <- 1.e+2
      sig <- min(c(0.3, max(c((1 - alfpri)^c1, (1 - alfdual)^c1, (normR - mu) / (normR + c2 * mu)))))
    }

    # Update reporting information
    muiter <- c(muiter, mu)
    normNiter <- c(normNiter, normR)
  }

  lam <- lam / row
  obj <- sum(log(Psi %*% lam))
  lam <- lam / sum(lam)

  list(fobj = obj, lambda = lam)
}
