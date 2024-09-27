#' Plot Image
#'
#' @param A 2x2 integer matrix
#' @param dd data frame with names x and y containing two points of parallelogram
#' @param det Boolean changes the form of plot by printing determinant
#'
#' @export
plot_image <- function(A, dd = data.frame(x = c(1, 1), y = c(1, -1)), det = FALSE) {
  dd <- rbind(dd, -dd) %>%
    distinct()
  dd$piece = 1
  #dd <- rename(dd, piece = p)
  dd <- mutate(dd, point = as.character(1:n()))
  dd
  hull <- dd %>%
    group_by(piece) %>%
    slice(chull(x, y))

  # Find the convex hull of the points being plotted
  hull_11 <- dd %>%
    filter(piece == 1) %>%
    slice(chull(x, y)) %>%
    select(x, y)
  hull_12 <- dd %>%
    filter(piece == 2) %>%
    slice(chull(x, y)) %>%
    select(x, y)
  # Define the scatterplot
  #A <- matrix(c(0,-10,1,0), nrow = 2)
  dd2 <- purrr::map_df(1:nrow(dd), function(x) {
    val <- A %*% c(dd$x[x], dd$y[x])
    data.frame(x = val[1], y = val[2], piece = dd$piece[x], color = dd$point[x])
  })
  p <- ggplot(dd, aes(x, y)) + geom_point(shape = 21, mapping = aes(fill = point), size = 2) + geom_point(shape = 21, data = dd2, aes(fill = color), size = 2)
  # Overlay the convex hull
  #p + geom_polygon(mapping = aes(x, y, group = piece), data = hull, alpha = 0.5)

  hull2 <- dd2 %>%
    group_by(piece) %>%
    slice(chull(x, y))
  if(det) {
    p + geom_polygon(mapping = aes(x, y, group = piece), data = hull, alpha = 0.5)  +
      geom_polygon(mapping = aes(x, y, group = piece), data = hull2, alpha = 0.2) +
      ggtitle(label = det(A))
    # scale_y_continuous(breaks = integer_breaks() ) +
    # scale_x_continuous(breaks = integer_breaks() )

  } else {
    p + geom_polygon(mapping = aes(x, y, group = piece), data = hull, alpha = 0.5)  +
      geom_polygon(mapping = aes(x, y, group = piece), data = hull2, alpha = 0.2) +
      annotate("text", x=1, y=-.6, label= "K") +
      annotate("text", x = 4, y = -1.8, label = "A(K)") +
      theme_minimal() +
      ggtitle(label = "Strictly Expansive Example")
  }
}

yget_u <- function(u, AA) {
  P <- AA %*% matrix(c(1, u + 1, 0, 1), nrow = 2)
  AAP_inv <- solve(P)
  l1 <- max(abs(AAP_inv %*% matrix(c(1, u), ncol = 1)))
  l2 <- max(abs(AAP_inv %*% matrix(c(1, u + 2), ncol = 1)))
  max(l1, l2)
}

xget_u <- function(u, AA) {
  P <- AA %*% matrix(c(u + 1, 1, 1, 0), nrow = 2)
  AAP_inv <- solve(P)
  l1 <- max(abs(AAP_inv %*% matrix(c(u, 1), ncol = 1)))
  l2 <- max(abs(AAP_inv %*% matrix(c(u + 2, 1), ncol = 1)))
  max(l1, l2)
}

#' Eigenvector Transform
#'
#' This is only useful to me for testing hypotheses about matrices with complex eigenvalues.
#'
#' @param A 2x2 integer matrix with comples eigenvalues
transform_eigenvectors <- function(A) {
  ee <- eigen(A)$vec
  mag <- abs(ee)
  if(mag[1] > mag[2]) {
    ee[,1] * (Conj(ee[2,1])/abs(ee[2,1]))
  } else {
    ee[,1] * (Conj(ee[1,1])/abs(ee[1,1]))
  }
}

mult_mat_df <- function(A, dd) {
  ret <- data.frame(t(A %*% t(dd)))
  names(ret) <- c("x", "y")
  ret
}


