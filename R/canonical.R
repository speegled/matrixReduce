

#' Determine Whether Matrix is Reduced
#'
#' @param A 2x2 integer matrix
#' @export
is_reduced <- function(A) {
  del <- round(trace(A)^2 - 4 * det(A))
  a <- A[1,1]
  b <- A[1,2]
  c <- A[2,1]
  d <- A[2,2]
  if(del < 0) {
    if(abs(d-a) <= c && c <= -b) {
      if(abs(d - a) == c || c == -b) {
        return(d >= a)
      }
      return(TRUE)
    }
    return(FALSE)
  }
  if(del > 0 && !is_square(del)) {
    return(c > 0 && abs(sqrt(del) - 2 * c) < (d - a) && (d - a) < sqrt(del))
  }
  if(is_square(del)) {
    if(a != d) {
      return(0 <= b && b <= a - d)
    } else{
      return(b >= 0)
    }
  }
  stop("fall through in is_reduced")
}

is_square <- function(x) {
  if(x < 0) {
    return(FALSE)
  }
  abs(sqrt(x) - round(sqrt(x))) < .00000001
}

#' Find Canonical Matrix
#'
#' @param A a 2x2 integer matrix
#' @param return_p Boolean whether to return P such that return value is P^{-1} AP
#' @param debug Boolean
#'
#' @export
ff <- function(A, return_p = FALSE, debug = FALSE) {
  PP <- matrix(c(1, 0, 0, 1), nrow = 2)
  if(debug) {
    browser()
  }
  i <- 0
  del <- round(trace(A)^2 - 4 * det(A))
  if(del < 0) {
    if(A[2,1] < 0) {
      A <- S1 %*% A %*% solve(S1)
      PP <- PP %*% solve(S1)
    }
    while(!is_reduced(A) && i < 100) {
      if(A[2,1] > -A[1,2] || (A[2,1] == -A[1,2] && -A[2,1] <= A[2,2] - A[1,1] && A[2,2] - A[1,1] < 0 )) {
        A <- S2 %*% A %*% solve(S2)
        PP <- PP %*% solve(S2)
      }
      if(!is_reduced(A)) {
        n <- ceiling((A[2,2] - A[1,1])/(2 * A[2,1]) - 1/2)
        A <- Sn(n) %*% A %*% solve(Sn(n))
        PP <- PP %*% solve(Sn(n))
      }
      i <- i + 1
    }
  }
  if(del > 0 && !is_square(del)) {
    while(!is_reduced(A) && i < 100) {
      n <- (rfunc(A[1,1] - A[2,2], A[1,2], del) + A[2,2] - A[1,1])/(2 * A[1,2])
      if(A[1,2] < 0) {
        P1 <- matrix(c(-n, -1, 1, 0), nrow = 2)
        A <- round(P1 %*% A %*% solve(P1))
        PP <- PP %*% solve(P1)
      } else {
        P1 <- matrix(c(-n, 1, 1, 0), nrow = 2)
        A <- round(P1 %*% A %*% solve(P1))
        PP <- PP %*% solve(P1)
      }
      i <- i + 1
    }
    O <- paste(A, collapse =" ")
    while(!anyDuplicated(O) && i < 100) {
      i <- i + 1
      n <- (rfunc(A[1,1] - A[2,2], A[1,2], del) + A[2,2] - A[1,1])/(2 * A[1,2])
      if(A[1,2] < 0) {
        P1 <- matrix(c(-n, -1, 1, 0), nrow = 2)
        A <- round(P1 %*% A %*% solve(P1))
        PP <- PP %*% solve(P1)
      } else {
        P1 <- matrix(c(-n, 1, 1, 0), nrow = 2)
        A <- round(P1 %*% A %*% solve(P1))
        PP <- PP %*% solve(P1)
      }
      O <- c(O, paste(round(A), collapse = " "))
    }
    len <- length(O) - 1
    A <- list(len)
    for(i in 1:len) {
      A[[i]] <- matrix(as.integer(unlist(strsplit(O[i], split = " "))), nrow = 2)
    }
  }
  if(is_square(del) && !is_reduced(A)) {
    if(round(A[1,1] - A[2,2]) < 0) {
      S <- matrix(c(0, 1, 1, 0), nrow = 2)
      A <- S %*% A %*% solve(S)
      PP <- PP %*% solve(S)
    }
    eig <- eigen(A)$value
    eig <- ifelse(abs(eig - Re(eig)) < .00001, Re(eig), eig)
    if(abs(eig[1] - eig[2]) > .000001) {
      ee <- get_evec(A, eig, debug)
      zw <- numbers::extGCD(ee[1,1], ee[2,1])
      zw[3] <- -zw[3]
      P <- matrix(c(ee[,1], zw[3:2]), nrow = 2)
      A <- round(solve(P) %*% A %*% P)
      PP <- PP %*% P
      if(A[1,2] > 0) {
        n <- floor(A[1,2]/(A[1,1] - A[2,2]))
        A <- round(Sn(n) %*% A %*% Sn(-n))
        PP <- PP %*% Sn(-n)
      }
      if(A[1,2] < 0) {
        n <- floor(A[1,2]/(A[1,1] - A[2,2]))
        A <- round(Sn(n) %*% A %*% Sn(-n))
        PP <- PP %*% Sn(-n)
      }
    } else {
      ee <- get_evec(A, eig, debug)
      zw <- numbers::extGCD(ee[1,1], ee[2,1])
      zw[3] <- -zw[3]
      P <- matrix(c(ee[,1], zw[3:2]), nrow = 2)
      A <- round(solve(P) %*% A %*% P)
      PP <- PP %*% P
      if(A[1,2] < 0) {
        A <- round(S1 %*% A %*% S1)
        PP <- PP %*% S1
      }
    }
  }
  if(i == 100) {
    warning("fall through")
  }
  if(return_p) {
    return(list(A = A, P = PP))
  }
  A
}
