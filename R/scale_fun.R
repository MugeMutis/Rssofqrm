#' @importFrom pcaPP l1median
scale_fun <- function(data){
  n <- dim(data)[1]
  p <- dim(data)[2]

  data <- as.matrix(data)
    if(p == 1){
      cen.data <- apply(data, 2, median)
    }else{
      cen.data <- l1median(data)
    }


  scl.data <- (data - matrix(cen.data, nrow = n, ncol = p, byrow = TRUE))

  return(scl.data)
}
