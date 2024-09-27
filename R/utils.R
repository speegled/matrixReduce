#' @keywords internal
trace <- function(A) {
  sum(diag(A))
}

S1 <- matrix(c(1,0,0,-1), nrow = 2)
S2 <- matrix(c(0,-1,1,0), nrow = 2)
Sn <- function(n) {
  matrix(c(1,0,n,1), nrow = 2)
}





get_evec <- function(A, eig, debug = FALSE) {
  if(debug) {
    browser()
  }
  if(any(abs(eigen(A)$value - eigen(A)$value) > .00000001)) {
    warning("non integer eigenvalues")
  }
  if(A[1,2] != 0) {
    x <- rep(A[1,2], 2)
    y <- eig  - A[1,1]
  } else if (A[2,1] != 0) {
    y <- rep(A[2,1], 2)
    x <- eig - A[2,2]
  } else {
    if(A[1,1] != 0 && A[2,2] != 0) {
      x <- c(1, 0)
      y <- c(0, 1)
    }
  }
  evec <- round(matrix(c(x, y), ncol = 2, byrow = T))
  evec[,1] <- evec[,1]/pracma::gcd(evec[1,1], evec[2,1])
  evec[,2] <- evec[,2]/pracma::gcd(evec[2,2], evec[1,2])
  evec
}



rfunc <- function(b, a, del) {
  #browser()
  if(abs(a) > sqrt(del)) {
    r <- b %% (2 * a)
    if(r > a) {
      r <- r - 2 * a
    }
  } else {
    r <- b %% (2 * a)
    while(r < sqrt(del) - 2 * abs(a)) {
      r <- r + 2 * abs(a)
    }
  }
  return(r)
}
