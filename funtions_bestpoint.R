#Funciones

unif.scale <- function(x){
  return((x -min(x))/(max(x)-min(x)))
}

# Funcion de escalamiento que utiliza la desviacion estandar poblacional
scale.p <- function (x, center = TRUE, scale = TRUE) 
{
  x <- as.matrix(x)
  nc <- ncol(x)
  if (is.logical(center)) {
    if (center) {
      center <- colMeans(x, na.rm = TRUE)
      x <- sweep(x, 2L, center, check.margin = FALSE)
    }
  }
  else {
    if (!is.numeric(center)) 
      center <- as.numeric(center)
    if (length(center) == nc) 
      x <- sweep(x, 2L, center, check.margin = FALSE)
    else stop("length of 'center' must equal the number of columns of 'x'")
  }
  if (is.logical(scale)) {
    if (scale) {
      f <- function(v) {
        v <- v[!is.na(v)]
        sqrt(sum(v^2)/length(v))
      }
      scale <- apply(x, 2L, f)
      x <- sweep(x, 2L, scale, `/`, check.margin = FALSE)
    }
  }
  else {
    if (!is.numeric(scale)) 
      scale <- as.numeric(scale)
    if (length(scale) == nc) 
      x <- sweep(x, 2L, scale, `/`, check.margin = FALSE)
    else stop("length of 'scale' must equal the number of columns of 'x'")
  }
  if (is.numeric(center)) 
    attr(x, "scaled:center") <- center
  if (is.numeric(scale)) 
    attr(x, "scaled:scale") <- scale
  x
}

# Grafico de correlaciones simbolicas

plot.sym.core <- function (prin.corre,gtitle) 
{
  v <- c("green", "red", "blue", "cyan", "brown", "yellow", 
         "pink", "purple", "orange", "gray")
  #msg <- paste("Correlation Circle")
  graphics::plot(-1:1, -1:1, type = "n", xlab = "C1", 
                 ylab = "C2",asp = 1,main = gtitle)
  graphics::abline(h = 0, lty = 3)
  graphics::abline(v = 0, lty = 3)
  graphics::symbols(0, 0, circles = 1, inches = FALSE, add = TRUE)
  c1 <- 1
  c2 <- 2
  n <- dim(prin.corre)[1]
  f <- dim(prin.corre)[2]
  CRTI <- prin.corre[, colnames(prin.corre) != "varname"]
  vars <- prin.corre$varname
  for (k in 1:n) {
    x1 <- min(CRTI[k, c1], CRTI[k, c2])
    x2 <- max(CRTI[k, c1], CRTI[k, c2])
    y1 <- min(CRTI[k, c2 + 1], CRTI[k, c2 + 2])
    y2 <- max(CRTI[k, c2 + 1], CRTI[k, c2 + 2])
    
    RSDA:::plotX.slice(x1, y1, x2, y2, v, vars, k)
  }
}

# Función para calcular angulos entre vectores
# thetasfun <- function(x){
#   theta1 <- atan2(x[3],x[1])
#   theta2 <- atan2(x[4],x[2])
#   theta <- max(abs(theta2),abs(theta1))-min(abs(theta2),abs(theta1))
#   return(theta)
# }

thetasfun <- function(x){
  theta1 <- atan2(x[3],x[1])
  theta2 <- atan2(x[4],x[2])
  
  #theta1 <- ifelse(theta1 >= 0, theta1, 2*pi + theta1)
  #theta2 <- ifelse(theta2 >= 0, theta2, 2*pi + theta2)
  #theta <- max(theta1,theta2)-min(theta1,theta2)
  theta <- abs(theta2-theta1)
  theta <- ifelse(theta <= pi, theta, 2*pi - theta)
  return(theta)
}


# Funcion para el metodo de centros
centers.pca.j <- function (sym.data) 
{
  N <- sym.data$N
  M <- sym.data$M
  seq.min <- seq(from = 1, by = 2, length.out = M)
  seq.max <- seq(from = 2, by = 2, length.out = M)
  sym.var.names <- sym.data$sym.var.names
  sym.data.vertex <- RSDA:::vertex.interval.new.j(sym.data)
  sym.data.vertex.matrix <- sym.data.vertex$vertex
  dim.vertex <- dim(sym.data.vertex.matrix)[1]
  tot.individuals <- N + dim.vertex
  min.interval <- as.matrix(sym.data$data[, seq.min])
  max.interval <- as.matrix(sym.data$data[, seq.max])
  centers.vector <- as.numeric(sapply(seq.min,FUN = function(j){
    (sym.data$data[,j] + sym.data$data[,j+1])/2
  }))
  M.x <- matrix(centers.vector, nrow = N)
  colnames(M.x) <- sym.var.names
  M.x <- scale.p(M.x)
  mean.var <- attr(M.x, "scaled:center")
  desv.var <- attr(M.x, "scaled:scale")
  sym.data.vertex.matrix.cent <- sym.data.vertex.matrix
  for (j in 1:M) {
    sym.data.vertex.matrix.cent[, j] <- (sym.data.vertex.matrix.cent[,j] - 
                                           mean.var[j])/desv.var[j]
  }
  M.x <- rbind(M.x, sym.data.vertex.matrix.cent)
  pca.opt <- FactoMineR::PCA(X = M.x, 
                             scale.unit = FALSE, 
                             ind.sup = (N + 1):tot.individuals, 
                             ncp = M, 
                             graph = FALSE)
  pca.opt.sym <- RSDA:::sym.interval.pca.limits.new.j(sym.data, 
                                                      pca.opt$ind.sup$coord, 
                                                      sym.data.vertex$num.vertex)
  
  V <- pca.opt$svd$U
  minm.c <- matrix(NA,nrow = N,ncol = M)
  maxm.c <- matrix(NA,nrow = N,ncol = M)
  
  for (j in 1:M) {
    minm.c[,j] <- (min.interval[,j] - mean.var[j])/desv.var[j]
    maxm.c[,j] <- (max.interval[,j] - mean.var[j])/desv.var[j]
  }
  
  vmpos <- pmax(V,0)
  vmneg <- pmin(V,0)
  
  rmin <- sapply(1:M,FUN = function(j){
    1/N*(colSums(maxm.c*vmneg[,j]) + colSums(minm.c*vmpos[,j]))
  })
  
  rmax <- sapply(1:M,FUN = function(j){
    1/N*(colSums(minm.c*vmneg[,j]) + colSums(maxm.c*vmpos[,j]))
  })
  
  Rm <- cbind(rmin, rmax)                 
  Rm <- as.data.frame(Rm[, c(matrix(1:(2*M), nrow = 2, byrow = TRUE))])
  rownames(Rm) <- sym.data$sym.var.names
  colnames(Rm) <- apply(expand.grid(c("min","max"), paste0("Dim",1:M)),
                        1, paste, collapse=".")
  Rm <- tibble::rownames_to_column(Rm,
                                   var = "varname")
  Rm <- tibble::as_tibble(Rm)
  
  return(list(Sym.Components = pca.opt.sym, 
              pca.opt = pca.opt,
              IPrinCore = Rm,
              scaled.var = list("mean" = mean.var, "desv" = desv.var)))
}

#Funcion de cos2 individuos
pca.ind.cos2.fun <- function (x, M, N, sym.var.names, 
                              sym.data.vertex.matrix, 
                              tot.individuals, num.dimen.aux) 
{
  M.x <- matrix(x, nrow = N)
  colnames(M.x) <- sym.var.names
  M.x <- rbind(M.x, sym.data.vertex.matrix)
  pca.min <- PCA(X = M.x, scale.unit = TRUE, ind.sup = (N + 
                                                          1):tot.individuals, ncp = M, graph = FALSE)
  
  # out <- list(pca.max = pca.max,
  #             out = -sum(pca.min$ind.sup$cos2[,1:num.dimen.aux]^2))
  return(-sum(pca.min$ind.sup$cos2[,1:num.dimen.aux]))
}

# pca.ind.cos2.fun <- function (x, M, N, sym.var.names, 
#                               sym.data.vertex.matrix, 
#                               tot.individuals, num.dimen.aux) 
# {
#   M.x <- matrix(x, nrow = N)
#   colnames(M.x) <- sym.var.names
#   M.x <- scale.p(M.x)
#   mean.var <- attr(M.x, "scaled:center")
#   desv.var <- attr(M.x, "scaled:scale")
#   sym.data.vertex.matrix.cent <- sym.data.vertex.matrix
#   for (j in 1:M) {
#     sym.data.vertex.matrix.cent[, j] <- (sym.data.vertex.matrix.cent[, 
#                                                                      j] - mean.var[j])/desv.var[j]
#   }
#   M.x <- rbind(M.x, sym.data.vertex.matrix.cent)
#   pca.opt <- FactoMineR::PCA(X = M.x, 
#                              scale.unit = FALSE, 
#                              ind.sup = (N + 1):tot.individuals, 
#                              ncp = M, 
#                              graph = FALSE)
#   # out <- list(pca.max = pca.max,
#   #             out = -sum(pca.min$ind.sup$cos2[,1:num.dimen.aux]^2))
#   return(-sum(pca.opt$ind.sup$cos2[,1:num.dimen.aux]^2))
# }

#Funcion optimizadora individuos
optim.pca.ind.cos2.fun <- function (sym.data, num.dimension) 
{
  N <- sym.data$N
  M <- sym.data$M
  num.dimen.aux <- num.dimension
  seq.min <- seq(from = 1, by = 2, length.out = M)
  seq.max <- seq(from = 2, by = 2, length.out = M)
  sym.var.names <- sym.data$sym.var.names
  sym.data.vertex <- RSDA:::vertex.interval.new.j(sym.data)
  sym.data.vertex.matrix <- sym.data.vertex$vertex
  dim.vertex <- dim(sym.data.vertex.matrix)[1]
  tot.individuals <- N + dim.vertex
  min.interval <- as.matrix(sym.data$data[, seq.min])
  max.interval <- as.matrix(sym.data$data[, seq.max])
  init.point <- as.vector(as.matrix(RSDA:::centers.interval.j(sym.data)$centers))
  res.opt <- nloptr::lbfgs(init.point, 
                           pca.ind.cos2.fun, 
                           lower = as.vector(min.interval), 
                           upper = as.vector(max.interval), 
                           nl.info = FALSE, 
                           control = list(xtol_rel = 1e-08, 
                                          maxeval = 20000), 
                           N = N, 
                           M = M, 
                           sym.var.names = sym.var.names, 
                           sym.data.vertex.matrix = sym.data.vertex.matrix, 
                           tot.individuals = tot.individuals, 
                           num.dimen.aux = num.dimen.aux)
  M.x <- matrix(res.opt$par, nrow = N)
  colnames(M.x) <- sym.var.names
  M.x <- scale.p(M.x)
  mean.var <- attr(M.x, "scaled:center")
  desv.var <- attr(M.x, "scaled:scale")
  sym.data.vertex.matrix.cent <- sym.data.vertex.matrix
  for (j in 1:M) {
    sym.data.vertex.matrix.cent[, j] <- (sym.data.vertex.matrix.cent[, 
                                                                     j] - mean.var[j])/desv.var[j]
  }
  M.x <- rbind(M.x, sym.data.vertex.matrix.cent)
  pca.opt <- FactoMineR::PCA(X = M.x, 
                             scale.unit = FALSE, 
                             ind.sup = (N + 1):tot.individuals, 
                             ncp = M, 
                             graph = FALSE)
  pca.opt.sym <- RSDA:::sym.interval.pca.limits.new.j(sym.data, 
                                                      pca.opt$ind.sup$coord, 
                                                      sym.data.vertex$num.vertex)
  V <- pca.opt$svd$U
  minm.c <- matrix(NA,nrow = N,ncol = M)
  maxm.c <- matrix(NA,nrow = N,ncol = M)
  
  for (j in 1:M) {
    minm.c[,j] <- (min.interval[,j] - mean.var[j])/desv.var[j]
    maxm.c[,j] <- (max.interval[,j] - mean.var[j])/desv.var[j]
  }
  
  vmpos <- pmax(V,0)
  vmneg <- pmin(V,0)
  
  rmin <- sapply(1:M,FUN = function(j){
    1/N*(colSums(maxm.c*vmneg[,j]) + colSums(minm.c*vmpos[,j]))
  })
  
  rmax <- sapply(1:M,FUN = function(j){
    1/N*(colSums(minm.c*vmneg[,j]) + colSums(maxm.c*vmpos[,j]))
  })
  
  Rm <- cbind(rmin, rmax)                 
  Rm <- as.data.frame(Rm[, c(matrix(1:(2*M), nrow = 2, byrow = TRUE))])
  rownames(Rm) <- sym.data$sym.var.names
  colnames(Rm) <- apply(expand.grid(c("min","max"), paste0("Dim",1:M)),
                        1, paste, collapse=".")
  Rm <- tibble::rownames_to_column(Rm,
                                   var = "varname")
  Rm <- tibble::as_tibble(Rm)
  
  return(list(Sym.Components = pca.opt.sym, 
              pca.opt = pca.opt, 
              res.opt = res.opt,
              IPrinCore = Rm,
              scaled.var = list("mean" = mean.var, "desv" = desv.var)))
}

# Funcion de cos2 variables
pca.var.cos2.fun <- function (x, M, N, sym.var.names, 
                              num.dimen.aux,
                              min.interval,max.interval) 
{
  M.x <- matrix(x, nrow = N)
  colnames(M.x) <- sym.var.names
  M.x <- scale.p(M.x)
  mean.var <- attr(M.x, "scaled:center")
  desv.var <- attr(M.x, "scaled:scale")
  pca.opt <- FactoMineR::PCA(X = M.x, 
                             scale.unit = FALSE, 
                             ncp = M, 
                             graph = FALSE)
  V <- pca.opt$svd$U
  minm.c <- matrix(NA,nrow = N,ncol = M)
  maxm.c <- matrix(NA,nrow = N,ncol = M)
  
  for (j in 1:M) {
    minm.c[,j] <- (min.interval[,j] - mean.var[j])/desv.var[j]
    maxm.c[,j] <- (max.interval[,j] - mean.var[j])/desv.var[j]
  }
  
  vmpos <- pmax(V,0)
  vmneg <- pmin(V,0)
  
  rmin <- sapply(1:M,FUN = function(j){
    1/N*(colSums(maxm.c*vmneg[,j]) + colSums(minm.c*vmpos[,j]))
  })
  
  rmax <- sapply(1:M,FUN = function(j){
    1/N*(colSums(minm.c*vmneg[,j]) + colSums(maxm.c*vmpos[,j]))
  })
  
  Rm <- cbind(rmin, rmax)                 
  Rm <- Rm[, c(matrix(1:(2*M), nrow = 2, byrow = TRUE))]
  
  return(-sum(Rm[,1:(num.dimen.aux*2)]^2))
}



# Función optimizadora de variables
optim.pca.var.cos2.fun <- function (sym.data, num.dimension) 
{
  N <- sym.data$N
  M <- sym.data$M
  num.dimen.aux <- num.dimension
  seq.min <- seq(from = 1, by = 2, length.out = M)
  seq.max <- seq(from = 2, by = 2, length.out = M)
  sym.var.names <- sym.data$sym.var.names
  sym.data.vertex <- RSDA:::vertex.interval.new.j(sym.data)
  sym.data.vertex.matrix <- sym.data.vertex$vertex
  dim.vertex <- dim(sym.data.vertex.matrix)[1]
  tot.individuals <- N + dim.vertex
  min.interval <- as.matrix(sym.data$data[, seq.min])
  max.interval <- as.matrix(sym.data$data[, seq.max])
  init.point <- as.vector(as.matrix(RSDA:::centers.interval.j(sym.data)$centers))
  res.opt <- nloptr::lbfgs(init.point, 
                           pca.var.cos2.fun, 
                           lower = as.vector(min.interval), 
                           upper = as.vector(max.interval), 
                           nl.info = FALSE, 
                           control = list(xtol_rel = 1e-08, 
                                          maxeval = 20000), 
                           N = N, 
                           M = M, 
                           sym.var.names = sym.var.names, 
                           num.dimen.aux = num.dimen.aux,
                           min.interval = min.interval,
                           max.interval = max.interval)
  M.x <- matrix(res.opt$par, nrow = N)
  colnames(M.x) <- sym.var.names
  M.x <- scale.p(M.x)
  mean.var <- attr(M.x, "scaled:center")
  desv.var <- attr(M.x, "scaled:scale")
  sym.data.vertex.matrix.cent <- sym.data.vertex.matrix
  for (j in 1:M) {
    sym.data.vertex.matrix.cent[, j] <- (sym.data.vertex.matrix.cent[, 
                                                                     j] - mean.var[j])/desv.var[j]
  }
  M.x <- rbind(M.x, sym.data.vertex.matrix.cent)
  pca.opt <- FactoMineR::PCA(X = M.x, 
                             scale.unit = FALSE, 
                             ind.sup = (N + 1):tot.individuals, 
                             ncp = M, 
                             graph = FALSE)
  
  V <- pca.opt$svd$U
  minm.c <- matrix(NA,nrow = N,ncol = M)
  maxm.c <- matrix(NA,nrow = N,ncol = M)
  for (j in 1:M) {
    minm.c[,j] <- (min.interval[,j] - mean.var[j])/desv.var[j]
    maxm.c[,j] <- (max.interval[,j] - mean.var[j])/desv.var[j]
  }
  vmpos <- pmax(V,0)
  vmneg <- pmin(V,0)
  
  rmin <- sapply(1:M,FUN = function(j){
    1/N*(colSums(maxm.c*vmneg[,j]) + colSums(minm.c*vmpos[,j]))
  })
  
  rmax <- sapply(1:M,FUN = function(j){
    1/N*(colSums(minm.c*vmneg[,j]) + colSums(maxm.c*vmpos[,j]))
  })
  
  Rm <- cbind(rmin, rmax)                 
  Rm <- as.data.frame(Rm[, c(matrix(1:(2*M), nrow = 2, byrow = TRUE))])
  rownames(Rm) <- sym.var.names
  colnames(Rm) <- apply(expand.grid(c("min","max"), paste0("Dim",1:M)),
                        1, paste, collapse=".")
  Rm <- tibble::rownames_to_column(Rm,
                                   var = "varname")
  Rm <- tibble::as_tibble(Rm)
  
  pca.opt.sym <- RSDA:::sym.interval.pca.limits.new.j(sym.data, 
                                                      pca.opt$ind.sup$coord, 
                                                      sym.data.vertex$num.vertex)
  return(list(Sym.Components = pca.opt.sym, 
              pca.opt = pca.opt, 
              res.opt = res.opt,
              IPrinCore = Rm,
              scaled.var = list("mean" = mean.var, "desv" = desv.var)))
}

# Funcion de distancia

pca.distance.fun <- function (x, N, M, sym.var.names, sym.data.vertex.matrix, tot.individuals) 
{
  M.x <- matrix(x, nrow = N)
  colnames(M.x) <- sym.var.names
  M.x <- scale.p(M.x)
  mean.var <- attr(M.x, "scaled:center")
  desv.var <- attr(M.x, "scaled:scale")
  sym.data.vertex.matrix.cent <- sym.data.vertex.matrix
  for (j in 1:M) {
    sym.data.vertex.matrix.cent[, j] <- (sym.data.vertex.matrix.cent[, 
                                                                     j] - mean.var[j])/desv.var[j]
  }
  M.x <- rbind(M.x, sym.data.vertex.matrix.cent)
  pca.opt <- FactoMineR::PCA(X = M.x, scale.unit = FALSE, ind.sup = (N + 
                                                                       1):tot.individuals, ncp = M, graph = FALSE)
  opt.dist.pca <- pca.opt$ind.sup$dist * pca.opt$ind.sup$dist
  return(sum(opt.dist.pca))
}

# Funcion de optimizadora de distancia

optim.pca.distance.fun <- function(sym.data) 
{
  N <- sym.data$N
  M <- sym.data$M
  seq.min <- seq(from = 1, by = 2, length.out = M)
  seq.max <- seq(from = 2, by = 2, length.out = M)
  sym.var.names <- sym.data$sym.var.names
  sym.data.vertex <- RSDA:::vertex.interval.new.j(sym.data)
  sym.data.vertex.matrix <- sym.data.vertex$vertex
  dim.vertex <- dim(sym.data.vertex.matrix)[1]
  tot.individuals <- N + dim.vertex
  min.interval <- as.matrix(sym.data$data[, seq.min])
  max.interval <- as.matrix(sym.data$data[, seq.max])
  init.point <- as.vector(as.matrix(RSDA:::centers.interval.j(sym.data)$centers))
  res.opt <- nloptr::lbfgs(init.point, pca.distance.fun, 
                           lower = as.vector(min.interval), 
                           upper = as.vector(max.interval), nl.info = FALSE, 
                           control = list(xtol_rel = 1e-08, maxeval = 20000), N = N, 
                           M = M, sym.var.names = sym.var.names, sym.data.vertex.matrix = sym.data.vertex.matrix, 
                           tot.individuals = tot.individuals)
  M.x <- matrix(res.opt$par, nrow = N)
  colnames(M.x) <- sym.var.names
  M.x <- scale.p(M.x)
  mean.var <- attr(M.x, "scaled:center")
  desv.var <- attr(M.x, "scaled:scale")
  sym.data.vertex.matrix.cent <- sym.data.vertex.matrix
  for (j in 1:M) {
    sym.data.vertex.matrix.cent[, j] <- (sym.data.vertex.matrix.cent[, 
                                                                     j] - mean.var[j])/desv.var[j]
  }
  M.x <- rbind(M.x, sym.data.vertex.matrix.cent)
  pca.opt <- FactoMineR::PCA(X = M.x, scale.unit = FALSE, ind.sup = (N + 
                                                                       1):tot.individuals, ncp = M, graph = FALSE)
  pca.opt.sym <- RSDA:::sym.interval.pca.limits.new.j(sym.data, pca.opt$ind.sup$coord, 
                                                      sym.data.vertex$num.vertex)
  V <- pca.opt$svd$U
  minm.c <- matrix(NA,nrow = N,ncol = M)
  maxm.c <- matrix(NA,nrow = N,ncol = M)
  
  for (j in 1:M) {
    minm.c[,j] <- (min.interval[,j] - mean.var[j])/desv.var[j]
    maxm.c[,j] <- (max.interval[,j] - mean.var[j])/desv.var[j]
  }
  
  vmpos <- pmax(V,0)
  vmneg <- pmin(V,0)
  
  rmin <- sapply(1:M,FUN = function(j){
    1/N*(colSums(maxm.c*vmneg[,j]) + colSums(minm.c*vmpos[,j]))
  })
  
  rmax <- sapply(1:M,FUN = function(j){
    1/N*(colSums(minm.c*vmneg[,j]) + colSums(maxm.c*vmpos[,j]))
  })
  
  Rm <- cbind(rmin, rmax)                 
  Rm <- as.data.frame(Rm[, c(matrix(1:(2*M), nrow = 2, byrow = TRUE))])
  rownames(Rm) <- sym.data$sym.var.names
  colnames(Rm) <- apply(expand.grid(c("min","max"), paste0("Dim",1:M)),
                        1, paste, collapse=".")
  Rm <- tibble::rownames_to_column(Rm,
                                   var = "varname")
  Rm <- tibble::as_tibble(Rm)
  
  return(list(Sym.Components = pca.opt.sym, pca.opt = pca.opt, 
              res.opt = res.opt,
              IPrinCore = Rm,
              scaled.var = list("mean" = mean.var, "desv" = desv.var)))
}

# Funcion de varianza

pca.variance.fun <- function (x, M, N, sym.var.names, sym.data.vertex.matrix, tot.individuals, 
                              num.dimen.aux) 
{
  M.x <- matrix(x, nrow = N)
  colnames(M.x) <- sym.var.names
  M.x <- scale.p(M.x)
  mean.var <- attr(M.x, "scaled:center")
  desv.var <- attr(M.x, "scaled:scale")
  sym.data.vertex.matrix.cent <- sym.data.vertex.matrix
  for (j in 1:M) {
    sym.data.vertex.matrix.cent[, j] <- (sym.data.vertex.matrix.cent[, 
                                                                     j] - mean.var[j])/desv.var[j]
  }
  M.x <- rbind(M.x, sym.data.vertex.matrix.cent)
  pca.opt <- FactoMineR::PCA(X = M.x, scale.unit = FALSE, ind.sup = (N + 
                                                                       1):tot.individuals, ncp = M, graph = FALSE)
  #out <- list(pca.max = pca.max, out = -sum(pca.max$eig[(1:num.dimen.aux)]))
  return(-sum(pca.opt$eig[(1:num.dimen.aux)]))
}

# Funcion de optimizadora de varianza

optim.pca.variance.fun <- function (sym.data, num.dimension) 
{
  N <- sym.data$N
  M <- sym.data$M
  num.dimen.aux <- num.dimension
  seq.min <- seq(from = 1, by = 2, length.out = M)
  seq.max <- seq(from = 2, by = 2, length.out = M)
  sym.var.names <- sym.data$sym.var.names
  sym.data.vertex <- RSDA:::vertex.interval.new.j(sym.data)
  sym.data.vertex.matrix <- sym.data.vertex$vertex
  dim.vertex <- dim(sym.data.vertex.matrix)[1]
  tot.individuals <- N + dim.vertex
  min.interval <- as.matrix(sym.data$data[, seq.min])
  max.interval <- as.matrix(sym.data$data[, seq.max])
  init.point <- as.vector(as.matrix(RSDA:::centers.interval.j(sym.data)$centers))
  res.opt <- nloptr::lbfgs(init.point, pca.variance.fun, 
                           lower = as.vector(min.interval), 
                           upper = as.vector(max.interval), nl.info = FALSE, 
                           control = list(xtol_rel = 1e-08, maxeval = 20000), N = N, 
                           M = M, sym.var.names = sym.var.names, sym.data.vertex.matrix = sym.data.vertex.matrix, 
                           tot.individuals = tot.individuals, num.dimen.aux = num.dimen.aux)
  M.x <- matrix(res.opt$par, nrow = N)
  colnames(M.x) <- sym.var.names
  M.x <- scale.p(M.x)
  mean.var <- attr(M.x, "scaled:center")
  desv.var <- attr(M.x, "scaled:scale")
  sym.data.vertex.matrix.cent <- sym.data.vertex.matrix
  for (j in 1:M) {
    sym.data.vertex.matrix.cent[, j] <- (sym.data.vertex.matrix.cent[, 
                                                                     j] - mean.var[j])/desv.var[j]
  }
  M.x <- rbind(M.x, sym.data.vertex.matrix.cent)
  pca.opt <- FactoMineR::PCA(X = M.x, scale.unit = FALSE, ind.sup = (N + 
                                                                       1):tot.individuals, ncp = M, graph = FALSE)
  pca.opt.sym <- RSDA:::sym.interval.pca.limits.new.j(sym.data, pca.opt$ind.sup$coord, 
                                                      sym.data.vertex$num.vertex)
  V <- pca.opt$svd$U
  minm.c <- matrix(NA,nrow = N,ncol = M)
  maxm.c <- matrix(NA,nrow = N,ncol = M)
  
  for (j in 1:M) {
    minm.c[,j] <- (min.interval[,j] - mean.var[j])/desv.var[j]
    maxm.c[,j] <- (max.interval[,j] - mean.var[j])/desv.var[j]
  }
  
  vmpos <- pmax(V,0)
  vmneg <- pmin(V,0)
  
  rmin <- sapply(1:M,FUN = function(j){
    1/N*(colSums(maxm.c*vmneg[,j]) + colSums(minm.c*vmpos[,j]))
  })
  
  rmax <- sapply(1:M,FUN = function(j){
    1/N*(colSums(minm.c*vmneg[,j]) + colSums(maxm.c*vmpos[,j]))
  })
  
  Rm <- cbind(rmin, rmax)                 
  Rm <- as.data.frame(Rm[, c(matrix(1:(2*M), nrow = 2, byrow = TRUE))])
  thetas <- apply(Rm,1,thetasfun)
  rownames(Rm) <- sym.data$sym.var.names
  colnames(Rm) <- apply(expand.grid(c("min","max"), paste0("Dim",1:M)),
                        1, paste, collapse=".")
  Rm <- tibble::rownames_to_column(Rm,
                                   var = "varname")
  Rm <- tibble::as_tibble(Rm)
  
  return(list(Sym.Components = pca.opt.sym, pca.opt = pca.opt, 
              res.opt = res.opt,
              IPrinCore = Rm,
              scaled.var = list("mean" = mean.var, "desv" = desv.var)))
}

# Funcion de angulos
pca.theta.fun <- function (x, M, N, sym.var.names,
                           min.interval,max.interval) 
{
  M.x <- matrix(x, nrow = N)
  colnames(M.x) <- sym.var.names
  M.x <- scale.p(M.x)
  mean.var <- attr(M.x, "scaled:center")
  desv.var <- attr(M.x, "scaled:scale")
  pca.opt <- FactoMineR::PCA(X = M.x, 
                             scale.unit = FALSE, 
                             ncp = M, 
                             graph = FALSE)
  V <- pca.opt$svd$U
  minm.c <- matrix(NA,nrow = N,ncol = M)
  maxm.c <- matrix(NA,nrow = N,ncol = M)
  
  for (j in 1:M) {
    minm.c[,j] <- (min.interval[,j] - mean.var[j])/desv.var[j]
    maxm.c[,j] <- (max.interval[,j] - mean.var[j])/desv.var[j]
  }
  
  vmpos <- pmax(V,0)
  vmneg <- pmin(V,0)
  
  rmin <- sapply(1:M,FUN = function(j){
    1/N*(colSums(maxm.c*vmneg[,j]) + colSums(minm.c*vmpos[,j]))
  })
  
  rmax <- sapply(1:M,FUN = function(j){
    1/N*(colSums(minm.c*vmneg[,j]) + colSums(maxm.c*vmpos[,j]))
  })
  
  Rm <- cbind(rmin, rmax)                 
  Rm <- Rm[, c(matrix(1:(2*M), nrow = 2, byrow = TRUE))]
  thetas <- apply(Rm,1,thetasfun)
  
  return(sum(thetas^2))
}

# Función optimizadora de angulos
optim.pca.theta.fun <- function (sym.data) 
{
  N <- sym.data$N
  M <- sym.data$M
  #num.dimen.aux <- num.dimension
  seq.min <- seq(from = 1, by = 2, length.out = M)
  seq.max <- seq(from = 2, by = 2, length.out = M)
  sym.var.names <- sym.data$sym.var.names
  sym.data.vertex <- RSDA:::vertex.interval.new.j(sym.data)
  sym.data.vertex.matrix <- sym.data.vertex$vertex
  dim.vertex <- dim(sym.data.vertex.matrix)[1]
  tot.individuals <- N + dim.vertex
  min.interval <- as.matrix(sym.data$data[, seq.min])
  max.interval <- as.matrix(sym.data$data[, seq.max])
  init.point <- as.vector(as.matrix(RSDA:::centers.interval.j(sym.data)$centers))
  res.opt <- nloptr::lbfgs(init.point, 
                           pca.theta.fun, 
                           lower = as.vector(min.interval), 
                           upper = as.vector(max.interval), 
                           nl.info = FALSE, 
                           control = list(xtol_rel = 1e-08, 
                                          maxeval = 20000), 
                           N = N, 
                           M = M, 
                           sym.var.names = sym.var.names, 
                           #num.dimen.aux = num.dimen.aux,
                           min.interval = min.interval,
                           max.interval = max.interval)
  M.x <- matrix(res.opt$par, nrow = N)
  colnames(M.x) <- sym.var.names
  M.x <- scale.p(M.x)
  mean.var <- attr(M.x, "scaled:center")
  desv.var <- attr(M.x, "scaled:scale")
  sym.data.vertex.matrix.cent <- sym.data.vertex.matrix
  for (j in 1:M) {
    sym.data.vertex.matrix.cent[, j] <- (sym.data.vertex.matrix.cent[, 
                                                                     j] - mean.var[j])/desv.var[j]
  }
  M.x <- rbind(M.x, sym.data.vertex.matrix.cent)
  pca.opt <- FactoMineR::PCA(X = M.x, 
                             scale.unit = FALSE, 
                             ind.sup = (N + 1):tot.individuals, 
                             ncp = M, 
                             graph = FALSE)
  
  V <- pca.opt$svd$U
  minm.c <- matrix(NA,nrow = N,ncol = M)
  maxm.c <- matrix(NA,nrow = N,ncol = M)
  for (j in 1:M) {
    minm.c[,j] <- (min.interval[,j] - mean.var[j])/desv.var[j]
    maxm.c[,j] <- (max.interval[,j] - mean.var[j])/desv.var[j]
  }
  vmpos <- pmax(V,0)
  vmneg <- pmin(V,0)
  
  rmin <- sapply(1:M,FUN = function(j){
    1/N*(colSums(maxm.c*vmneg[,j]) + colSums(minm.c*vmpos[,j]))
  })
  
  rmax <- sapply(1:M,FUN = function(j){
    1/N*(colSums(minm.c*vmneg[,j]) + colSums(maxm.c*vmpos[,j]))
  })
  
  Rm <- cbind(rmin, rmax)                 
  Rm <- as.data.frame(Rm[, c(matrix(1:(2*M), nrow = 2, byrow = TRUE))])
  rownames(Rm) <- sym.var.names
  colnames(Rm) <- apply(expand.grid(c("min","max"), paste0("Dim",1:M)),
                        1, paste, collapse=".")
  Rm <- tibble::rownames_to_column(Rm,
                                   var = "varname")
  Rm <- tibble::as_tibble(Rm)
  
  pca.opt.sym <- RSDA:::sym.interval.pca.limits.new.j(sym.data, 
                                                      pca.opt$ind.sup$coord, 
                                                      sym.data.vertex$num.vertex)
  return(list(Sym.Components = pca.opt.sym, 
              pca.opt = pca.opt, 
              res.opt = res.opt,
              IPrinCore = Rm,
              scaled.var = list("mean" = mean.var, "desv" = desv.var)))
}




sym.pca.complete <- function (sym.data, method = c("classic", "tops", "centers","optimized.var.cos2",
                                                   "optimized.angles",
                                                   "principal.curves", "optimized.distance", "optimized.variance", 
                                                   "optimized.ind.cos2","fixed"), d = 3,fixed.matrix = NULL, ...) 
{
  all_interval <- all(sapply(sym.data, function(x) any(class(x) %in% 
                                                         "symbolic_interval")))
  if (!all_interval) {
    stop("All variables have to be of the same type")
  }
  sym.data <- RSDA:::to.v2(sym.data)
  method <- match.arg(method)
  if (method == "classic") {
    if ((sym.data$sym.var.types[1] != "$C") && (sym.data$sym.var.types[1] != 
                                                "$I")) {
      stop("Variables have to be continuos or Interval")
    }
    if (sym.data$sym.var.types[1] == "$C") {
      res <- FactoMineR::PCA(sym.data$data, scale.unit = TRUE, 
                             ncp = sym.data$M, graph = FALSE)
    }
    else if (sym.data$sym.var.types[1] == "$I") {
      nn <- sym.data$N
      mm <- sym.data$M
      centers <- matrix(0, nn, mm)
      centers <- as.data.frame(centers)
      rownames(centers) <- sym.data$sym.obj.names
      colnames(centers) <- sym.data$sym.var.names
      for (i in 1:nn) {
        for (j in 1:mm) {
          centers[i, j] <- (sym.var(sym.data, j)$var.data.vector[i, 
                                                                 1] + sym.var(sym.data, j)$var.data.vector[i, 
                                                                                                           2])/2
        }
      }
      res <- FactoMineR::PCA(centers, scale.unit = TRUE, 
                             ncp = sym.data$M, graph = FALSE)
    }
    return(res)
  }
  if (method == "centers") {
    res <- centers.pca.j(sym.data)
    class(res) <- c("symbolic_pca", class(res))
    return(res)
  }
  if (method == "tops") {
    res <- vertex.pca.j(sym.data)
    class(res) <- c("symbolic_pca", class(res))
    return(res)
  }
  if (method == "principal.curves") {
    res <- sym.interval.pc(sym.data, "vertex", 150, FALSE, 
                           FALSE, TRUE)
    class(res) <- c("symbolic_pca_curves", class(res))
    return(res)
  }
  if (method == "optimized.distance") {
    res <- optim.pca.distance.fun(sym.data)
    class(res) <- c("symbolic_pca_optimized", class(res))
    return(res)
  }
  if (method == "optimized.variance") {
    res <- optim.pca.variance.fun(sym.data, num.dimension = d)
    class(res) <- c("symbolic_pca_optimized", class(res))
    return(res)
  }
  if (method == "optimized.ind.cos2") {
    res <- optim.pca.ind.cos2.fun(sym.data, num.dimension = d)
    class(res) <- c("symbolic_pca_optimized", class(res))
    return(res)
  }
  if (method == "optimized.var.cos2") {
    res <- optim.pca.var.cos2.fun(sym.data, num.dimension = d)
    class(res) <- c("symbolic_pca_optimized", class(res))
    return(res)
  }
  if (method == "optimized.angles") {
    res <- optim.pca.theta.fun(sym.data)
    class(res) <- c("symbolic_pca_optimized", class(res))
    return(res)
  }
  if (method == "fixed") {
    res <- fixed.pca.j.new(sym.data, fixed.matrix)
    return(res)
  }
  return(TRUE)
}

fitness <- function(i, m = M, n = N,
                    var.names = sym.var.names,
                    vertex.matrix = sym.data.vertex.matrix, 
                    nind = tot.individuals, 
                    d = 2,
                    min.v = min.interval,
                    max.v = max.interval){
  #M.x <- matrix(x, nrow = N)
  M.x <- matrix(x[i,], nrow = N)
  colnames(M.x) <- var.names
  M.x <- scale.p(M.x)
  mean.var <- attr(M.x, "scaled:center")
  desv.var <- attr(M.x, "scaled:scale")
  
  sym.data.vertex.matrix.cent <- vertex.matrix
  for (j in 1:m) {
    sym.data.vertex.matrix.cent[, j] <- (sym.data.vertex.matrix.cent[, 
                                                                     j] - mean.var[j])/desv.var[j]
  }
  M.x <- rbind(M.x, sym.data.vertex.matrix.cent)
  pca.opt <- FactoMineR::PCA(X = M.x, 
                             scale.unit = FALSE, 
                             ind.sup = (n + 1):nind, 
                             ncp = m, 
                             graph = FALSE)
  V <- pca.opt$svd$U
  minm.c <- matrix(NA,nrow = N,ncol = m)
  maxm.c <- matrix(NA,nrow = N,ncol = m)
  
  for (j in 1:m) {
    minm.c[,j] <- (min.v[,j] - mean.var[j])/desv.var[j]
    maxm.c[,j] <- (max.v[,j] - mean.var[j])/desv.var[j]
  }
  
  vmpos <- pmax(V,0)
  vmneg <- pmin(V,0)
  
  rmin <- sapply(1:m,FUN = function(j){
    1/n*(colSums(maxm.c*vmneg[,j]) + colSums(minm.c*vmpos[,j]))
  })
  
  rmax <- sapply(1:m,FUN = function(j){
    1/n*(colSums(minm.c*vmneg[,j]) + colSums(maxm.c*vmpos[,j]))
  })
  
  Rm <- cbind(rmin, rmax)                 
  Rm <- Rm[, c(matrix(1:(2*M), nrow = 2, byrow = TRUE))]
  thetas <- apply(Rm,1,thetasfun)
  
  cr1 <- sum(pca.opt$ind.sup$dist^2)
  cr2 <- -sum(pca.opt$eig[(1:d)])
  cr3 <- -sum(pca.opt$ind.sup$cos2[,1:d]^2)
  cr4 <- -sum(Rm[,1:(d*2)]^2)
  cr5 <- sum(thetas^2) 
  return(c(cr1,cr2,cr3,cr4,cr5))
}

caRamel2 <- function (nobj, nvar, minmax, bounds, func, popsize, archsize, 
                      maxrun, prec, repart_gene = c(5, 5, 5, 5), gpp = NULL, blocks = NULL, 
                      pop = NULL, funcinit = NULL, objnames = NULL, listsave = NULL, 
                      write_gen = FALSE, carallel = 1, numcores = NULL, graph = TRUE, 
                      sensitivity = FALSE, verbose = TRUE, worklist = NULL) 
{
  start_time <- Sys.time()
  if (nobj <= 1) {
    if (verbose) 
      message("the number of objectives must be greater than one!")
    return(list(success = FALSE, message = "the number of objectives must be greater than one!"))
  }
  if (nvar < 1) {
    if (verbose) 
      message("the number of variables must be greater than zero!")
    return(list(success = FALSE, message = "the number of variables must be greater than zero!"))
  }
  if (length(minmax) != nobj) {
    if (verbose) 
      message("the dimension of 'minmax' is incorrect!")
    return(list(success = FALSE, message = "the dimension of 'minmax' is incorrect!"))
  }
  if ((dim(bounds)[1] != nvar) | (dim(bounds)[2] != 2)) {
    if (verbose) 
      message("the dimension of 'bounds' is incorrect!")
    return(list(success = FALSE, message = "the dimension of 'bounds' is incorrect!"))
  }
  if (!"function" %in% class(func)) {
    if (verbose) 
      message("'func' is not an R function!")
    return(list(success = FALSE, message = "'func' is not an R function!"))
  }
  if (popsize < 1) {
    if (verbose) 
      message("'popsize' must be strictly positive!")
    return(list(success = FALSE, message = "'popsize' must be strictly positive!"))
  }
  if (archsize < 1) {
    if (verbose) 
      message("'archsize' must be strictly positive!")
    return(list(success = FALSE, message = "'archsize' must be strictly positive!"))
  }
  if (maxrun < 1) {
    if (verbose) 
      message("'maxrun' must be strictly positive!")
    return(list(success = FALSE, message = "'maxrun' must be strictly positive!"))
  }
  if (length(repart_gene) != 4) {
    if (verbose) 
      message("the dimension of'repart_gene' must be 4!")
    return(list(success = FALSE, message = "the dimension of 'repart_gene' must be 4!"))
  }
  if (sum(repart_gene <= 0) > 0) {
    if (verbose) 
      message("parameter values for each rule of 'repart_gene' must be strictly positive!")
    return(list(success = FALSE, message = "parameter values for each rule of 'repart_gene' must be strictly positive!"))
  }
  if (!is.null(gpp)) {
    if (gpp < 1) {
      if (verbose) 
        message("gpp must be greater than zero!")
      return(list(success = FALSE, message = "gpp be greater than zero!"))
    }
  }
  initialise_calc <- 0
  if (!is.null(funcinit)) {
    if (!"function" %in% class(funcinit)) {
      if (verbose) 
        message("'funcinit' is not an R function!")
      return(list(success = FALSE, message = "'funcinit' is not an R function!"))
    }
    initialise_calc <- 1
  }
  if (is.null(objnames)) 
    objnames <- paste("Obj", as.character(c(1:nobj)), sep = "")
  writefile <- FALSE
  if (!is.null(listsave)) {
    if (!"list" %in% class(listsave)) {
      if (verbose) 
        message("'listsave' is not an R list!")
      return(list(success = FALSE, message = "'listsave' is not an R list!"))
    }
    writefile <- TRUE
    if (is.null(listsave$pmt)) {
      if (verbose) 
        message(" 'listsave$pmt' must be defined!")
      return(list(success = FALSE, message = " 'listsave$pmt' must be defined!"))
    }
    if (is.null(listsave$obj)) {
      if (verbose) 
        message(" 'listsave$obj' must be defined!")
      return(list(success = FALSE, message = " 'listsave$obj' must be defined!"))
    }
    if (is.null(listsave$evol)) {
      if (verbose) 
        message(" 'listsave$evol' must be defined!")
      return(list(success = FALSE, message = " 'listsave$evol' must be defined!"))
    }
    ecrit_total_pop <- 0
    if (!is.null(listsave$totalpop)) {
      ecrit_total_pop <- 1
    }
  }
  if (write_gen == TRUE) {
    if (writefile == FALSE) {
      if (verbose) 
        message(" 'listsave' must be defined to use write_gen!")
      return(list(success = FALSE, message = " 'listsave' must be defined to use write_gen!"))
    }
    listsave$RadPmt <- gsub(pattern = ".txt", replacement = "", 
                            listsave$pmt)
    listsave$RadObj <- gsub(pattern = ".txt", replacement = "", 
                            listsave$obj)
    if (ecrit_total_pop == 1) 
      listsave$RadPop <- gsub(pattern = ".txt", replacement = "", 
                              listsave$totalpop)
  }
  carallel <- as.integer(carallel)
  if (!(carallel %in% c(0, 1, 2))) {
    if (verbose) 
      message("'carallel' value is not 0, 1 or 2")
    return(list(success = FALSE, message = "wrong value for paramemter 'carallel'"))
  }
  save_crit <<- c()
  sp <- (bounds[, 2] - bounds[, 1])/(2 * sqrt(3))
  gsearch <- ceiling((nvar/10) * nobj/log(nobj + 1))
  precis <- matrix(data = prec, nrow = nobj)
  nrun <- 0
  ngen <- 0
  if (is.null(gpp)) 
    gpp <- ceiling(nvar * (nobj + 1) * 4/sum(repart_gene))
  jacobian <- NA
  if (carallel == 1) {
    if (is.null(numcores)) 
      numcores <- detectCores()
    cl <- makeCluster(numcores)
    if (initialise_calc == 1) 
      funcinit(cl, numcores)
  }
  if (!is.null(pop)) {
    pop <- as.matrix(pop)
    if (length(pop[1, ]) < (nvar + nobj)) {
      x <<- pop[, 1:nvar]
      if (carallel == 1) {
        newfeval <- NULL
        clusterExport(cl = cl, varlist = c("x", "worklist"), 
                      envir = environment())
        res <- parLapply(cl, 1:dim(x)[1], func)
        for (j in 1:dim(x)[1]) {
          newfeval <- rbind(newfeval, as.numeric(res[[j]][1:nobj]))
        }
      }
      else if (carallel == 0) {
        newfeval <- matrix(data = 0, nrow = dim(x)[1], 
                           ncol = nobj)
        for (i in 1:dim(x)[1]) {
          res <- func(i)
          newfeval[i, ] <- res[1:nobj]
        }
      }
      else {
        newfeval <- func(pop[, 1:nvar])
      }
      pop <- cbind(pop[, 1:nvar], newfeval)
      nrun <- nrun + length(pop[, 1])
    }
  }
  if (verbose) {
    message(paste("Beginning of caRamel optimization <--", 
                  date()))
    message(paste("Number of variables :", as.character(nvar)))
    message(paste("Number of functions :", as.character(nobj)))
    pb <- txtProgressBar(min = 0, max = 1, initial = 0, title = "caRamel progress :", 
                         label = "caRamel progress :", style = 3)
  }
  while (nrun < maxrun) {
    ngen <- ngen + 1
    if (is.null(pop)) {
      A <- t(matrix(data = bounds[, 2] - bounds[, 1], ncol = popsize, 
                    nrow = nvar))
      B <- t(matrix(data = bounds[, 1], ncol = popsize, 
                    nrow = nvar))
      x <- A * matrix(runif(popsize * nvar), ncol = nvar) + 
        B
      probj <- matrix(data = NaN, nrow = popsize, ncol = nobj)
    }
    else {
      vamax <- (ngen%%gpp) == 0
      param <- as.matrix(pop[, 1:nvar])
      dim(param) <- c(dim(pop)[1], nvar)
      if (dim(param)[1] < 4) {
        if (carallel == 1) 
          stopCluster(cl)
        if (verbose) {
          close(pb)
          message("Optimization failed")
        }
        return(list(success = FALSE, message = "The number of feasible points is not sufficient! Try to increase the size of the population..."))
      }
      crit <- pop[, (nvar + 1):(nvar + nobj)]
      dim(crit) <- c(dim(pop)[1], nobj)
      Xp <- newXval2(param, crit, minmax, sp, bounds, repart_gene, 
                     blocks, vamax)
      x <- Xp$x
      probj <- Xp$pcrit
    }
    additional_eval <- NULL
    if (carallel == 1) {
      newfeval <- NULL
      clusterExport(cl = cl, varlist = c("x", "worklist"), 
                    envir = environment())
      res <- parLapply(cl, 1:dim(x)[1], func)
      for (j in 1:dim(x)[1]) {
        newfeval <- rbind(newfeval, as.numeric(res[[j]][1:nobj]))
        additional_eval <- rbind(additional_eval, res[[j]][(nobj + 
                                                              1):length(res[[1]])])
      }
      nadditional <- ncol(additional_eval)
      if (length(res[[1]]) < (nobj + 1)) {
        additional_eval <- NULL
        nadditional <- 0
      }
    }
    else if (carallel == 0) {
      newfeval <- matrix(data = 0, nrow = dim(x)[1], ncol = nobj)
      x <<- x
      for (i in 1:dim(x)[1]) {
        res <- func(i)
        newfeval[i, ] <- res[1:nobj]
        additional_eval <- rbind(additional_eval, res[(nobj + 
                                                         1):length(res)])
      }
      nadditional <- ncol(additional_eval)
      if (length(res) < (nobj + 1)) {
        additional_eval <- NULL
        nadditional <- 0
      }
    }
    else {
      newfeval <- func(x)
      additional_eval <- NULL
      nadditional <- 0
    }
    nrun <- nrun + dim(x)[1]
    if (verbose) 
      setTxtProgressBar(pb, min(nrun, maxrun)/maxrun)
    detect_nan <- is.na(newfeval)
    set_ok <- !rowSums(detect_nan)
    if (length(set_ok[set_ok == TRUE]) == 0) {
      if (carallel == 1) 
        stopCluster(cl)
      if (verbose) {
        message("Optimization failed")
        close(pb)
      }
      return(list(success = FALSE, message = "No feasible points! Try to increase the size of the population..."))
    }
    newfeval <- newfeval[set_ok, ]
    dim(newfeval) <- c(sum(set_ok, na.rm = TRUE), nobj)
    x <- x[set_ok, ]
    dim(x) <- c(sum(set_ok, na.rm = TRUE), nvar)
    additional_eval <- additional_eval[set_ok, ]
    probj <- probj[set_ok, ]
    pop1 <- rbind(pop, cbind(x, newfeval, additional_eval))
    matobj <- pop1[, (nvar + 1):(nvar + nobj)]
    dim(matobj) <- c(dim(pop1)[1], nobj)
    ind <- decrease_pop(matobj, minmax, prec, archsize, popsize)
    arch <- matrix(pop1[ind$arch, ], nrow = length(ind$arch), 
                   ncol = nobj + nvar + nadditional)
    pop <- pop1[c(ind$arch, ind$pop), ]
    dim(pop) <- c(length(c(ind$arch, ind$pop)), nvar + nobj)
    param_arch <- arch[, 1:nvar]
    crit_arch <- matrix(arch[, (nvar + 1):(nvar + nobj)], 
                        nrow = length(ind$arch), ncol = nobj)
    if (nadditional > 0) {
      additional_eval <- matrix(arch[, (nvar + nobj + 1):(nvar + 
                                                            nobj + nadditional)], nrow = length(ind$arch), 
                                ncol = nadditional)
    }
    a <- c(lapply(c(1:nobj), function(i) max(crit_arch[, 
                                                       i])))
    maxcrit <- as.data.frame(a, col.names = objnames)
    a <- c(lapply(c(1:nobj), function(i) min(crit_arch[, 
                                                       i])))
    mincrit <- as.data.frame(a, col.names = objnames)
    crit <- mincrit
    crit[minmax] <- maxcrit[minmax]
    save_crit <- cbind(save_crit, c(nrun, t(crit)))
    if (graph == TRUE) 
      plot_population(MatObj = crit_arch, nobj, ngen, nrun, 
                      objnames, MatEvol = save_crit, popsize)
    if (writefile == TRUE) {
      if (write_gen == TRUE) {
        listsave$pmt <- paste(listsave$RadPmt, "_gen", 
                              ngen, ".txt", sep = "")
        listsave$obj <- paste(listsave$RadObj, "_gen", 
                              ngen, ".txt", sep = "")
        if (ecrit_total_pop == 1) {
          listsave$totalpop <- paste(listsave$RadPop, 
                                     "_gen", ngen, ".txt", sep = "")
        }
      }
      write.table(param_arch, listsave$pmt, row.names = FALSE, 
                  col.names = FALSE)
      write.table(cbind(crit_arch, additional_eval), listsave$obj, 
                  row.names = FALSE, col.names = FALSE)
      write.table(t(save_crit), listsave$evol, row.names = FALSE, 
                  col.names = FALSE)
      if (ecrit_total_pop == 1) {
        write.table(pop, listsave$totalpop, row.names = FALSE, 
                    col.names = FALSE)
      }
    }
  }
  if (verbose) 
    close(pb)
  if (sensitivity == TRUE) {
    if (verbose) 
      message("Computing the sensitivity of the Pareto front...")
    dx <- 1e-04
    dxinv <- 1/dx
    jacobian <- list()
    xopt <- param_arch
    nfront <- dim(crit_arch)[1]
    dim(xopt) <- c(nfront, nvar)
    x <- xopt
    for (k in 1:nobj) {
      nameId <- paste("Jacobian_", toString(k), sep = "")
      jacobian[[nameId]] <- matrix(data = 0, nrow = nfront, 
                                   ncol = nvar)
      dim(jacobian[[nameId]]) <- c(nfront, nvar)
    }
    for (j in 1:nvar) {
      x[, j] <- x[, j] + dx
      if (carallel == 1) {
        newfeval <- NULL
        clusterExport(cl = cl, varlist = c("x", "worklist"), 
                      envir = environment())
        res <- parLapply(cl, 1:nfront, func)
        for (e in 1:nfront) {
          newfeval <- rbind(newfeval, as.numeric(res[[e]][1:nobj]))
        }
      }
      else if (carallel == 0 | carallel == 2) {
        newfeval <- matrix(data = 0, nrow = nfront, ncol = nobj)
        x <<- x
        for (e in 1:nfront) {
          res <- func(e)
          newfeval[e, ] <- res[1:nobj]
        }
      }
      for (k in 1:nobj) {
        jacobian[[k]][, j] <- (newfeval[, k] - crit_arch[, 
                                                         k]) * dxinv
      }
      x[, j] <- xopt[, j]
    }
    nrun <- nrun + dim(x)[1]
  }
  if (carallel == 1) 
    stopCluster(cl)
  end_time <- Sys.time()
  if (verbose) {
    message(paste("Done in", as.character(end_time - start_time), 
                  units(end_time - start_time), "-->", date()))
    message(paste("Size of the Pareto front :", as.character(dim(crit_arch)[1])))
    message(paste("Number of calls :", as.character(nrun)))
  }
  return(list(success = TRUE, parameters = param_arch, objectives = cbind(crit_arch, 
                                                                          additional_eval), derivatives = jacobian, save_crit = t(save_crit), 
              total_pop = pop, gpp = gpp))
}

newXval2 <- function (param, crit, isperf, sp, bounds, repart_gene, blocks, 
                      fireworks) 
{
  nobj <- dim(crit)[2]
  npar <- dim(param)[2]
  nvec <- dim(param)[1]
  a <- 3/8
  obj <- matrix(0, nrow = dim(crit)[1], ncol = nobj)
  for (i in 1:nobj) {
    obj[, i] <- (val2rank(crit[, i], 3) - a)/(nvec + 1 - 
                                                2 * a)
  }
  obj[, !isperf] <- 1 - obj[, !isperf]
  Fo <- dominate(obj)
  param_arch <- as.matrix(param[Fo == 1, ])
  if (dim(param_arch)[2] < npar) {
    param_arch <- t(param_arch)
  }
  obj_arch <- obj[Fo == 1, ]
  crit_arch <- crit[Fo == 1, ]
  n_inter <- repart_gene[1]
  n_extra <- repart_gene[2]
  n_cov <- repart_gene[3]
  n_recomb <- repart_gene[4]
  xnew <- NULL
  project_crit <- NULL
  if (n_inter > 0 | n_extra > 0) {
    simplices <- suppressMessages(delaunayn2(obj,options = "Qt Qc Qx Qbb"))
    nf <- apply((matrix(Fo[simplices], nrow = dim(simplices)[1], 
                        ncol = dim(simplices)[2]) == 1), 1, sum)
    ix <- which(nf > 0)
    simplices <- simplices[ix, ]
    simplices <- matrix(data = simplices, ncol = (nobj + 
                                                    1))
    nbsimp <- length(ix)
    na <- nobj * (nobj + 1)/2
    volume <- matrix(0, nrow = nbsimp, ncol = 1)
    oriedge <- NULL
    ledge <- NULL
    for (s in 1:nbsimp) {
      P <- simplices[s, ]
      S <- obj[P, ]
      volume[s] <- vol_splx(S)
      direc <- Dimprove(S, Fo[P])
      if (!is.null(direc$oriedge)) {
        oriedge <- rbind(oriedge, matrix(P[direc$oriedge], 
                                         nrow = dim(direc$oriedge)[1]))
        ledge <- c(ledge, direc$ledge)
      }
    }
    unik <- !duplicated(oriedge)
    oriedge <- oriedge[unik, ]
    ledge <- ledge[unik]
  }
  if (n_inter > 0) {
    carain <- Cinterp(param, crit, simplices, volume, n_inter)
    xnew <- rbind(xnew, carain$xnew)
    project_crit <- rbind(project_crit, carain$pcrit)
  }
  if (n_extra > 0) {
    caraex <- Cextrap(param, crit, oriedge, ledge, n_extra)
    xnew <- rbind(xnew, caraex$xnew)
    project_crit <- rbind(project_crit, caraex$pcrit)
  }
  if (n_cov > 0) {
    iref <- sort(unique(c(simplices)))
    xref <- as.matrix(param[iref, ])
    xcov <- Cusecovar(xref, sqrt(2), n_cov)
    critcov <- matrix(NaN, nrow = n_cov, ncol = nobj)
    xnew <- rbind(xnew, xcov)
    project_crit <- rbind(project_crit, critcov)
  }
  if (dim(param_arch)[1] > 1 & n_recomb > 0) {
    xrecomb <- Crecombination(param_arch, blocks, n_recomb)
    critrec <- matrix(NaN, nrow = n_recomb, ncol = nobj)
    xnew <- rbind(xnew, xrecomb)
    project_crit <- rbind(project_crit, critrec)
  }
  if (fireworks) {
    if (sum(Fo == 1) == 1) {
      obj_arch <- matrix(obj_arch, 1, nobj)
      crit_arch <- matrix(crit_arch, 1, nobj)
    }
    sp <- matrix(sp, nrow = 1, ncol = length(sp))
    m <- apply(obj_arch, 2, max)
    ipp <- apply(obj_arch, 2, which.max)
    m <- max(apply(obj_arch, 1, min))
    maximin <- which.max(apply(obj_arch, 1, min))
    ipp <- c(ipp, maximin)
    for (i in seq_len(length(ipp))) {
      xcloud <- matrix(param_arch[ipp[i], ], nrow = npar, 
                       ncol = npar, byrow = TRUE)
      ccloud <- matrix(crit_arch[ipp[i], ], nrow = npar, 
                       ncol = nobj, byrow = TRUE)
      devi <- c(rnorm(npar) * sp)
      xcloud <- xcloud + diag(x = devi, nrow = length(devi), 
                              ncol = length(devi))
      xcloud <- xcloud[sp > 0, ]
      ccloud <- ccloud[sp > 0, ]
      xnew <- rbind(xnew, xcloud)
      project_crit <- rbind(project_crit, ccloud)
    }
  }
  minp <- matrix(bounds[, 1], nrow = dim(xnew)[1], ncol = npar, 
                 byrow = TRUE)
  maxp <- matrix(bounds[, 2], nrow = dim(xnew)[1], ncol = npar, 
                 byrow = TRUE)
  ll <- which(xnew > maxp, arr.ind = TRUE)
  nelem <- dim(ll)[1]
  if (nelem > 0) {
    for (i in 1:nelem) {
      rnum <- ll[i, 1]
      cnum <- ll[i, 2]
      xnew[rnum, cnum] <- maxp[rnum, cnum]
    }
  }
  ll <- which(xnew < minp, arr.ind = TRUE)
  nelem <- dim(ll)[1]
  if (nelem > 0) {
    for (i in 1:nelem) {
      rnum <- ll[i, 1]
      cnum <- ll[i, 2]
      xnew[rnum, cnum] <- minp[rnum, cnum]
    }
  }
  return(list(xnew = xnew, pcrit = project_crit))
}

delaunayn2 <- function (p, options = NULL, output.options = NULL, full = FALSE) 
{
  tmp_stdout <- tempfile("Rf")
  tmp_stderr <- tempfile("Rf")
  on.exit(unlink(c(tmp_stdout, tmp_stderr)))
  if (is.data.frame(p)) {
    p <- as.matrix(p)
  }
  storage.mode(p) <- "double"
  if (any(is.na(p))) {
    stop("The first argument should not contain any NAs")
  }
  default.options <- "Qt Qc Qx"
  if (ncol(p) < 4) {
    default.options <- "Qt Qc Qz"
  }
  if (is.null(options)) {
    options <- default.options
  }
  options <- tryCatch(geometry:::qhull.options(options, output.options, 
                                               supported_output.options <- c("Fa", "Fn"), full = full), 
                      error = function(e) {
                        stop(e)
                      })
  if (!grepl("Qt", options) & !grepl("QJ", options)) {
    options <- paste(options, "Qt")
  }
  out <- .Call("C_delaunayn", p, as.character(options), tmp_stdout, 
               tmp_stderr, PACKAGE = "geometry")
  if (nrow(out$tri) > 0) {
    missing.points <- length(setdiff(seq(1, nrow(p)), unique(as.vector(out$tri))))
    if (missing.points > 0) {
      warning(paste0(missing.points, " points missing from triangulation.\nIt is possible that setting the 'options' argument of delaunayn may help.\nFor example:\noptions = \"", 
                     default.options, " Qbb\"\noptions = \"", default.options, 
                     " QbB\"\nIf these options do not work, try shifting the centre of the points\nto the origin by subtracting the mean coordinates from every point."))
    }
  }
  out[which(sapply(out, is.null))] <- NULL
  if (is.null(out$areas) & is.null(out$neighbours)) {
    attr(out$tri, "delaunayn") <- attr(out$tri, "delaunayn")
    return(out$tri)
  }
  class(out) <- "delaunayn"
  out$p <- p
  return(out)
}

evalfun <- function(x, m = M, n = N,
                    var.names = sym.var.names,
                    vertex.matrix = sym.data.vertex.matrix, 
                    nind = tot.individuals, 
                    d = 2,
                    min.v = min.interval,
                    max.v = max.interval,
                    symobj = sym.data){
  M.x <- matrix(x, nrow = N)
  colnames(M.x) <- var.names
  M.x <- scale.p(M.x)
  mean.var <- attr(M.x, "scaled:center")
  desv.var <- attr(M.x, "scaled:scale")
  
  sym.data.vertex.matrix.cent <- vertex.matrix
  for (j in 1:m) {
    sym.data.vertex.matrix.cent[, j] <- (sym.data.vertex.matrix.cent[, 
                                                                     j] - mean.var[j])/desv.var[j]
  }
  M.x <- rbind(M.x, sym.data.vertex.matrix.cent)
  pca.opt <- FactoMineR::PCA(X = M.x, 
                             scale.unit = FALSE, 
                             ind.sup = (n + 1):nind, 
                             ncp = m, 
                             graph = FALSE)
  
  pca.opt.sym <- RSDA:::sym.interval.pca.limits.new.j(sym.data, 
                                                      pca.opt$ind.sup$coord, 
                                                      sym.data.vertex$num.vertex)
  
  V <- pca.opt$svd$U
  minm.c <- matrix(NA,nrow = N,ncol = m)
  maxm.c <- matrix(NA,nrow = N,ncol = m)
  
  for (j in 1:m) {
    minm.c[,j] <- (min.v[,j] - mean.var[j])/desv.var[j]
    maxm.c[,j] <- (max.v[,j] - mean.var[j])/desv.var[j]
  }
  
  vmpos <- pmax(V,0)
  vmneg <- pmin(V,0)
  
  rmin <- sapply(1:m,FUN = function(j){
    1/n*(colSums(maxm.c*vmneg[,j]) + colSums(minm.c*vmpos[,j]))
  })
  
  rmax <- sapply(1:m,FUN = function(j){
    1/n*(colSums(minm.c*vmneg[,j]) + colSums(maxm.c*vmpos[,j]))
  })
  
  Rm <- cbind(rmin, rmax)                 
  Rm <- as.data.frame(Rm[, c(matrix(1:(2*M), nrow = 2, byrow = TRUE))])
  
  rownames(Rm) <- sym.var.names
  colnames(Rm) <- apply(expand.grid(c("min","max"), paste0("Dim",1:M)),
                        1, paste, collapse=".")
  Rm <- tibble::rownames_to_column(Rm,
                                   var = "varname")
  Rm <- tibble::as_tibble(Rm)
  
  return(list(Sym.Components = pca.opt.sym,
              pca.opt = pca.opt,
              IPrinCore = Rm,
              scaled.var = list("mean" = mean.var, "desv" = desv.var)))
}

