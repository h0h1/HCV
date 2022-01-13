dissimilarity <- function(feature_domain, dist_method){

  # calculate the dissimilarity matrix
  if(dist_method == 'correlation'){
    dist_matrix <- 1 - cor(t(feature_domain))
    return(dist_matrix)
  }

  if(dist_method == 'abscor'){
    dist_matrix <- 1 - abs(cor(t(feature_domain)))
    return(dist_matrix)
  }
  else{
    dist_matrix <- as.matrix(dist(feature_domain, method = dist_method))
  }
  return(dist_matrix)
}

LanceWilliams_algorithm <- function(n_i, n_j, n_h){

  # Dynamic update of the dissimilarity matrix

  n_k <- n_i + n_j # merge i-cluster and j-cluster, denoted by k-cluster

  dict <- vector(mode = 'list', length = 8) # Build the dictionary to store the LanceWilliams coefficients
  names(dict) <- c('single', 'complete', 'average', 'centroid', 'ward', 'median', 'weight')
  dict[[1]] <- c(0.5,0.5,0,-0.5)
  dict[[2]] <- c(0.5,0.5,0,0.5)
  dict[[3]] <- c(n_i/n_k, n_j/n_k, 0, 0)
  dict[[4]] <- c(n_i/n_k, n_j/n_k, -n_i*n_j/n_k/n_k, 0)
  dict[[5]] <- c((n_i+n_h)/(n_i+n_j+n_h), (n_j+n_h)/(n_i+n_j+n_h), -n_h/(n_i+n_j+n_h), 0) 
  dict[[6]] <- c(0.5, 0.5, -0.25, 0)
  dict[[7]] <- c(0.5, 0.5, 0, 0)
  return(dict)
}

DistanceBetweenCluster <- function(linkage, alpha_i, alpha_j, beta, gamma, d_hi, d_hj, d_ij){

  # return the distance between cluster by given linkage

  #if(linkage == 'ward.D2'){
  #  return((alpha_i * (d_hi ** 2) + alpha_j * (d_hj ** 2) + beta * (d_ij ** 2) + gamma * (d_hi - d_hj) ** 2) ** 0.5)
  #}
  #else{
    return(alpha_i * d_hi + alpha_j * d_hj + beta * d_ij + gamma * abs(d_hi - d_hj))
  #}
}

#' Adjacency Matrix from Tessellation 
#' 
#' This function deals with spatial data having a point-level geometry domain.  It converts the spatial proximity into an adjacency matrix based on Voronoi tessellation or Delaunay triangulation.
#'
#' @param geometry_domain \emph{n} by \emph{d} matrix (NA not allowed) of geographical coordinates for \emph{n} points in \emph{d}-dimensional space.
#' @return A matrix with 0-1 values indicating the adjacency between the \emph{n} input points.
#' @export 
#' @author ShengLi Tzeng and Hao-Yun Hsu.  
#' @references  Gallier, J. (2011). Dirichlet–Voronoi Diagrams and Delaunay Triangulations. In Geometric Methods and Applications (pp. 301-319). Springer, New York, NY.
#' @examples
#' require(fields)
#' require(alphahull)
#' pts <- ChicagoO3$x
#' rownames(pts) <- LETTERS[1:20]
#' Vcells <- delvor(pts)
#' plot(Vcells,wlines='vor',pch='.')
#' text(pts,rownames(pts))
#' Amat <- tessellation_adjacency_matrix(pts)  
#
tessellation_adjacency_matrix <- function(geometry_domain){

  n <- nrow(geometry_domain)
  adj_mat <- matrix(0, nrow=n, ncol=n)
  cell <- geometry::delaunayn(geometry_domain, output.options = T)$tri
  ndim <- dim(cell)
  for(k in 1:ndim[1]){
    pair <- cell[k, ]
    for(i in 1:(ndim[2]-1)){
      for(j in (i+1):ndim[2]){
        adj_mat[pair[i], pair[j]] <- 1        
      }
    }
  }
  adj_mat <- pmin( adj_mat + t(adj_mat),1) 
  rownames(adj_mat)=colnames(adj_mat)=rownames(geometry_domain)
  return(adj_mat)
}

FindMinimum <- function(dist_matrix, adj_mat){

  # Find the minimum value among the lower triangle part in a dissimilarity matrix

  min_dist = c(-1,-1, Inf) # The first two are the (i, j) terns and third is the minimum value
  weighted_matrix <- dist_matrix / adj_mat
  weighted_matrix[adj_mat==0]=Inf
  n <- nrow(weighted_matrix)
  for(i in 2:n){
    for(j in 1:(i-1)){
      dist <- weighted_matrix[i, j]
      if(dist <= min_dist[3]){
        min_dist[3] = dist
        min_dist[1] = i
        min_dist[2] = j
      }
    }
  }
  return(min_dist)
}

AGNES <- function(dist_matrix, adj_mat, linkage, iterate){
  count <- 0
  n <- nrow(dist_matrix)
  dslist <- list() # data structure list
  height <- rep(0, n - 1)
  merge <- matrix(0, ncol = 2, nrow = n - 1)
  n_stat <- nrow(dist_matrix)
  cuttree <- matrix(-1, nrow = n - iterate + 2, ncol = n_stat)
  cluster <- vector(mode = 'list', length = n)
  key <- which(c('single', 'complete', 'average', 'centroid',  'ward',
                 'median', 'weight') == linkage)				 
  for(i in 1:n){
    cluster[[i]] <- i
  }
  while(iterate <= n + 1){
    min_dist <- FindMinimum(dist_matrix, adj_mat)
    n_i <- length(cluster[[min_dist[1]]])
    n_j <- length(cluster[[min_dist[2]]])
    for(h in 1:n){
      n_h <- length(cluster[[h]])
      coef <- LanceWilliams_algorithm(n_i, n_j, n_h)[[key]]

      if(h != min_dist[1]){
        dist_matrix[h, min_dist[1]] <- DistanceBetweenCluster(
          linkage, coef[1], coef[2],coef[3],coef[4],
          dist_matrix[h, min_dist[1]], dist_matrix[h, min_dist[2]],
          dist_matrix[min_dist[1], min_dist[2]]
        )

        dist_matrix[min_dist[1], h] <- dist_matrix[h, min_dist[1]]

        adj_mat[h, min_dist[1]] <- (adj_mat[h, min_dist[1]] | adj_mat[h, min_dist[2]])
        adj_mat[min_dist[1], h] <- adj_mat[h, min_dist[1]]
      }
    }


    dist_matrix <- dist_matrix[,-min_dist[2], drop=F]
    dist_matrix <- dist_matrix[-min_dist[2],, drop=F]
    adj_mat <- adj_mat[,-min_dist[2], drop=F]
    adj_mat <- adj_mat[-min_dist[2],, drop=F]

    n <- n - 1
    count <- count + 1

    dslist[[count]] <- c(cluster[[min_dist[1]]], cluster[[min_dist[2]]])
    height[count] <- min_dist[[3]]

    ##### Generate the merge attribute #####
    if(length(cluster[[min_dist[1]]]) == 1){
      merge[count, 1] <- -cluster[[min_dist[1]]]
    }
    if(length(cluster[[min_dist[2]]]) == 1){
      merge[count, 2] <- -cluster[[min_dist[2]]]
    }

    if(length(cluster[[min_dist[1]]]) > 1){
      for(s in count:1){
        if(all(cluster[[min_dist[1]]] %in% dslist[[s]])){
          merge[count, 1] <- s
        }
      }
    }
    if(length(cluster[[min_dist[2]]]) > 1){
      for(s in count:1){
        if(all(cluster[[min_dist[2]]] %in% dslist[[s]])){
          merge[count, 2] <- s
        }
      }
    }


    cluster[[min_dist[1]]] <- c(cluster[[min_dist[1]]], cluster[[min_dist[2]]])
    cluster <- cluster[-min_dist[2]]
    cuttree[n - iterate + 3,] <- Clusterlabels(cluster, n_stat)

  }
  return(list(cuttree, merge, height))
}

Clusterlabels <- function(cluster, n){

  # return the label for the iteration

  labels <- rep(-1, n)
  for(i in 1:length(cluster)){
    for(item in cluster[[i]]){
      labels[item] <- i
    }
  }
  return(labels)
}

merge2order <- function (m)
{
    LC = list()
    RC = list()
    n = NROW(m)
    for (i in 1:n) {
        if (m[i, 1] < 0)
            LC[[i]] = abs(m[i, 1])
        else {
            LC[[i]] = c(LC[[m[i, 1]]], RC[[m[i, 1]]])
        }
        if (m[i, 2] < 0)
            RC[[i]] = abs(m[i, 2])
        else {
            RC[[i]] = c(LC[[m[i, 2]]], RC[[m[i, 2]]])
        }
    }
    return(c(LC[[n]], RC[[n]]))
}

#' Hierarchical Clustering from Vertex-links
#' 
#' This function implements the hierarchical clustering for spatial data. It modified typically used hierarchical agglomerative clustering algorithms for introducing the spatial homogeneity, by considering geographical locations as vertices and converting spatial adjacency into whether a shared edge exists between a pair of vertices.  
#' 
#' @param geometry_domain  one of the three formats: (i) \emph{n} by \emph{d} matrix (NA not allowed), (ii)  a \code{SpatialPolygonsDataFrame} object defining polygons, (iii) a matrix with 0-1 value adjacency (with \code{adjacency=TRUE})
#' @param feature_domain  either (i) \emph{n} by \emph{p} matrix (NA allowed) for \emph{n} samples with \emph{p} attributes, or (ii) \emph{n} by \emph{n} matrix (NA not allowed) with dissimilarity between \emph{n} samples (with \code{diss = 'precomputed'})
#' @param linkage the agglomeration method to be used, one of "ward", "single", "complete", "average" (= UPGMA), "weight" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). Default is \code{'ward'}.
#' @param  diss  character indicating if \code{feature_domain} is a dissimilarity matrix: 'none' for not dissimilarity, and 'precomputed' for dissimilarity. Default is 'none'.
#' @param adjacency logical indicating if \code{geometry_domain } is a adjacency matrix. Default is FALSE.
#' @param dist_method the distance measure to be used when \code{feature_domain} is not a dissimilarity matrix (\code{diss ='none'}), one of "euclidean", "correlation", "abscor", "maximum", "manhattan", "canberra", "binary" or "minkowski".  Default is \code{'euclidean'}.
#' 
#' @return An object of class \code{hclust} which describes the tree produced by the clustering process. See the documentation in \code{hclust}.
#' @references  Carvalho, A. X. Y., Albuquerque, P. H. M., de Almeida Junior, G. R., and Guimaraes, R. D. (2009). Spatial hierarchical clustering. Revista Brasileira de Biometria, 27(3), 411-442.
#' @export 
#' @author ShengLi Tzeng and Hao-Yun Hsu.  
#' @details \code{geometry_domain} can be a user-specifid adjacency matrix, an \emph{n} by \emph{d} matrix with geographical coordinates for point-level data, or a \code{SpatialPolygonsDataFrame} object defining polygons for areal data. If an adjacency matrix is given, the user should use \code{adjacency=TRUE}.
#' @seealso \code{\link{hclust}}
#' @examples
#' set.seed(0)
#' pcase <- synthetic_data(3,30,0.02,100,2,2)
#' HCVobj <- HCV(pcase$geo,  pcase$feat)
#smi <- getCluster(pcase$geo,pcase$feat,HCVobj,method="SMI")
#' smi <- getCluster(HCVobj,method="SMI")
#' par(mfrow=c(2,2))
#' labcolor  <-  (pcase$labels+1)%%3+1
#' plot(pcase$feat,  col  =  labcolor,  pch=19,  xlab  =  'First  attribute', 
#'   ylab  =  'Second  attribute',  main  =  'Feature  domain')
#' plot(pcase$geo,  col  =  labcolor,  pch=19,  xlab  =  'First  attribute', 
#'   ylab  =  'Second  attribute',  main  =  'Geometry  domain')
#' plot(pcase$feat,  col=factor(smi),pch=19,  xlab  =  'First  attribute', 
#'   ylab  =  'Second  attribute',main  =  'Feature  domain')
#' plot(pcase$geo,  col=factor(smi),pch=19,  xlab  =  'First  attribute', 
#'   ylab  =  'Second  attribute',main  =  'Geometry  domain')
#'
HCV <- function(geometry_domain, feature_domain,
                                linkage='ward',  diss = 'none',
                                adjacency=FALSE,dist_method = 'euclidean'){
  n <- nrow(geometry_domain)
  dimension <- ncol(geometry_domain)
  if(dimension < 2){
    cat('Error : The input dimension of geometry domain should be greater than 1')
    return()
  }

  ##### distance matrix preprocessing #####
  if(diss == 'precomputed'){
    dist_matrix <- feature_domain
  }
  else{
    dist_matrix <- dissimilarity(feature_domain, dist_method)
  }
  dist_matrix <- as.matrix(dist_matrix)

  ##### adjacency matrix preprocessing #####
  if(adjacency == T){
    adj_mat <- geometry_domain*1
  }
  else{
    if("SpatialPolygonsDataFrame"%in%class(geometry_domain)) #areal data
    adj_mat <- rgeos::gTouches(geometry_domain,byid=T)*1 else{       #point-level data	
     if(is.null(rownames(geometry_domain)))
	 {
	  rownames(geometry_domain) <- rownames(feature_domain)
      warning(" No rownames for geometry_domain and/or feature_domain. \n Be sure that the two datasets have the same order for samples.")	 
	 } 
	 adj_mat <- tessellation_adjacency_matrix(geometry_domain)	 
	}
  }
  adj_mat <- as.matrix(adj_mat)
 
  if(is.null(rownames(dist_matrix))) 
     rnames <- rownames(feature_domain) else
	 rnames <- rownames(dist_matrix)
  if( (!is.null(rownames(adj_mat))) & (!is.null(rnames)) ) {
     if(is.null(colnames(adj_mat))) 
	    colnames(adj_mat) <- rownames(adj_mat)
	 adj_mat <- adj_mat[rnames,rnames] 
  } else 
  warning(" No rownames for geometry_domain and/or feature_domain. \n Be sure that the two datasets have the same order for samples.")

  result <- AGNES(dist_matrix, adj_mat, linkage, 3)
  clust <- list()

  rname <- rownames(dist_matrix)
  if(is.null(rname)) rname <- rep('x',1:nrow(dist_matrix))
  clust$call <- match.call()
  clust$hmatrix <- result[[1]]     
  clust$labels <- rname
  clust$merge <- result[[2]]
  clust$height <- result[[3]]
  clust$method <- linkage
  clust$dist.method <- dist_method
  clust$adjacency.matrix <- Matrix(adj_mat, sparse=TRUE)
  clust$dist.matrix <- dist_matrix
  clust$order <- merge2order(result[[2]])
  clust$connected <- 'Connected'
  clust$cut_for_connect <- 1
  class(clust) <- 'hclust'

  for(i in 1:length(result[[3]])){

    # height
    if(result[[3]][i]==Inf){
      cat('The graph is not connected\n')
      cat('Use cutree(output,', n-i+1, ') to find the components')
      clust$connected <- 'Not connected'
      clust$cut_for_connect <- n-i+1
      break
    }
  }

  return(clust)
}


#' Generating Point-level Data Having Several Groups 
#'
#' Generation of synthetic point-level data based on a method proposed by Lin et al. (2005).
#' 
#' @param k integer specifying the number of groups.
#' @param f positive number controlling the concentration of generated samples toward large groups.
#' @param r positive number controlling the variance of individual attributes on the feature domain.
#' @param n integer specifying the total number of sampled points.
#' @param feature integer specifying the number of attributes for the feature domain.
#' @param geometry integer specifying the number of attributes for the geometry domain.
#' @param homogeneity logical indicating whether to force the centers of the feature domain to be the same as those of the geometry domain. Default is TRUE.
#' @return A list with two matrices and a vector of labels. One matrix is for the feature domain and the other is for the geometry domain, both of which have  \emph{n} sampled points. The vector of labels indicates which cluster each sample belongs to.
#' @export 
#' @author ShengLi Tzeng and Hao-Yun Hsu.  
#' @references Lin, C. R., Liu, K. H., and Chen, M. S. (2005). Dual clustering: integrating data clustering over optimization and constraint domains. IEEE Transactions on Knowledge and Data Engineering, 17(5), 628-637.
#' @examples
#' set.seed(0)
#' pcase <- synthetic_data(3,30,0.02,100,2,2)
#' par(mfrow=c(1,2))
#' labcolor <- (pcase$labels+1)%%3+1
#' plot(pcase$feat, col = labcolor, pch=19, xlab = 'First attribute', 
#'   ylab = 'Second attribute', main = 'Feature domain')
#' plot(pcase$geo, col = labcolor, pch=19, xlab = 'First attribute', 
#'   ylab = 'Second attribute', main = 'Geometry domain')
#' 
synthetic_data <- function(k, f, r, n, feature, geometry, homogeneity = TRUE){
  geometry_domain <- matrix(0, ncol = geometry, nrow = n)
  feature_domain <- matrix(0, ncol = feature, nrow = n)
  labels <- rep(-1, n)
  geometry_domain_center = runif(geometry * k, 0, 1)
  geometry_domain_center = matrix(geometry_domain_center, ncol = geometry)
  feature_domain_center = runif(feature * k, 0, 1)
  feature_domain_center = matrix(feature_domain_center, ncol = feature)
  if(homogeneity){
    feature_domain_center <- geometry_domain_center
  }
  for(z in 1:n){
    prob <- rep(0, k)
    p <- runif(geometry, 0, 1)
    geometry_domain[z,] <- p
    for(i in 1:k){
      q <- (1 / (sum((p - geometry_domain_center[i,]) ** 2))**0.5) ** f
      sums <- 0
      for(j in 1:k){
        sums <- sums + (1 / (sum((p - geometry_domain_center[j,]) ** 2))**0.5) ** f
      }
      prob[i] <- q / sums
    }
    labels[z] <- sample(1:k, replace = T, size = 1, prob = prob)
    mean <- feature_domain_center[labels[z],]
    cov <- diag(feature) * r
    feature_domain[z,] <- MASS::mvrnorm(1, mean, cov)
  }
  rownames(geometry_domain) <- 1:n
  rownames(feature_domain) <- 1:n
  list <- list(geo = geometry_domain, feat = feature_domain, labels = labels)
  return(list)
}

#' Drawing a Thematic Map with a Quantitative Feature
#'
#' Plot the polygons in a \code{SpatialPolygonsDataFrame} object, and turn the values of a quantitative feature into colors over individual polygons.
#' @param map  \code{SpatialPolygonsDataFrame} object consisting of data and polygons. 
#' @param feat numberic vector having the same elements as the number of polygons in the input \code{map}.
#' @param color vector of distinct colors for converting values of \code{feat}.
#' @param main character specifying the main title. 
#' @param bar_title character specifying the text over the color bar. 
#' @param zlim length-2 numeric vector specifying the range of values to be converted.
#' @return  A colored map.
#' @seealso \code{\link{SpatialPolygonsDataFrame}}
#' @examples
#' require(sp)
#' grd  <-  GridTopology(c(1,1),  c(1,1),  c(5,5))
#' polys  <-  as(grd,  "SpatialPolygons") 
#' centroids  <-  coordinates(polys)
#' gdomain  <-  SpatialPolygonsDataFrame(polys,  data=data.frame(x=centroids[,1],  
#'   y=centroids[,2], row.names=row.names(polys)))
#' feat <- gdomain$x*5+gdomain$y^2
#' plotMap(gdomain,feat)
#'
plotMap <- function(map, feat, color = topo.colors(10), main = "",
                    bar_title = "rank", zlim = NULL ) {
  layout(t(1:2),widths = c(6,1))
  par(mar = c(4,4,1,0.5))
  nc <- length(color) 
  if(is.null(zlim)) zlim=c(min(feat)-0.1^5,max(feat)+0.1^5)
  feat[feat>zlim[2]] <- zlim[2]
  feat[feat>zlim[2]] <- zlim[2]
  grp <- cut(feat, seq(zlim[1],zlim[2], l=nc) )  
  plot(map, col = color[grp], main = main)  
  par(mar = c(5,1,5,2.5))
  scaleby <- (max(zlim) - min(zlim)) / nc
  z <- seq(from = min(zlim), to = max(zlim), by = scaleby)
  image(y = z,
        z = t(z),
        col = color,
        axes = FALSE,
        main = bar_title,
        cex.main = 1)
  axis(4, cex.axis = 0.8, mgp = c(0,.5,0),las=1)  
}

TreeStructure <- function(hclust){

  ##### Return the fatherNode and the childNode #####
  merge <- hclust$merge
  n <- nrow(merge) + 1
  leftNode <- list()
  rightNode <- list()
  fatherNode <- list()
  height <- hclust$height

  for(i in 1:(n-1)){
    node <- merge[i,1]

    if(node < 0)
      leftNode[[i]] <- -node
    else
      leftNode[[i]] <- fatherNode[[node]]

    node <- merge[i,2]

    if(node < 0)
      rightNode[[i]] <- -node
    else
      rightNode[[i]] <- fatherNode[[node]]

    fatherNode[[i]] <- sort(c(leftNode[[i]], rightNode[[i]]))
  }

  return(list(leftNode=leftNode, rightNode=rightNode, fatherNode=fatherNode,
              height=height))
}


newCrit <- function( HCVobj, geodesic, k){
  n <- HCVobj$labels
  labels <- as.numeric(cutree(HCVobj, k))
  ek <- rep(0, k)
  WSS <- rep(0, k)
  for(i in 1:k){
    Ak <- geodesic[labels == i,labels == i]
	Dk <- HCVobj$dist.matrix[labels == i,labels == i]
    ek[i] <- mean(Ak[lower.tri(Ak)])^2 
	WSS[i] <- mean(Dk[lower.tri(Dk)]^2)
  }  
  edgesum <- sum(ek)
  alpha <- (k*ek - edgesum) / edgesum
  fa <- 1 / (1 + exp(-alpha))
  total <- sum(WSS * fa + ek * (1 - fa)) / k  
}

SMI <- function( HCVobj, max_k){
  index <- vector(length = max_k)
  geodesic <- dst(as.matrix(HCVobj$adjacency.matrix))  
  #e <- mean(geodesic[lower.tri(geodesic)])   
  Dk <- HCVobj$dist.matrix
  index[1] <- mean(Dk[lower.tri(Dk)]^2)
  gd <- geodesic[lower.tri(geodesic)]   
  fd <- Dk[lower.tri(Dk)]
  geodesic <- geodesic / mean(gd) * mean(fd)     
  
  for(i in 2:max_k){
    index[i] <- newCrit(HCVobj, geodesic, i)
  }
  #diff <- -diff(index)
  #ratio <- diff / index[-(max_k)]
  ratio <-  index[-(max_k)]/index[-1]
  bestCluster <- which.max(ratio)+1
  res <- list()
  res$bestCluster <- bestCluster
  res$index <- index
  res$ratio <- ratio
  return(res)
}



FindComponents <- function(adj, label){

  k <- max(label)
  n <- length(label)
  neighbor <- vector(mode = 'list', length = k)

  for(i in 1:n){
    neighbor[[label[i]]] <- c(neighbor[[label[i]]], i)
  }

  components <- list()
  count <- 0
  for(j in 1:k){
    invisible(utils::capture.output(
      hclust <- HCV(adj[neighbor[[j]], neighbor[[j]]],
                                    matrix(1,nrow=length(neighbor[[j]]))
                                    , adjacency = T))
    )

    comp <- as.numeric(cutree(hclust, hclust$cut_for_connect))

    for(l in 1:max(comp)){

      count <- count + 1
      components[[count]] <- neighbor[[j]][comp == l]

    }
  }

  assignments <- rep(0, n)
  for(i in 1:count){
    for(j in 1:length(components[[i]])){
      assignments[components[[i]][j]] <- i
    }
  }
  return(list(labels = assignments, neighbors = neighbor))
}

affinityMat <- function(hclust, kernel = 'none', KNN = 7){
  n <- nrow(hclust$merge) + 1
  Tree <- TreeStructure(hclust)
  leftNode <- Tree$leftNode
  rightNode <- Tree$rightNode
  height <- Tree$height
  aveheight <- mean(height)

  affinityMatrix <- matrix(Inf, ncol=n, nrow=n)
  diag(affinityMatrix) <- 0

  for(i in 1:(n-1)){
    for(l in leftNode[[n-i]]){
      for(r in rightNode[[n-i]]){
        affinityMatrix[l, r] <- min(affinityMatrix[l, r], height[n-i])
        affinityMatrix[r, l] <- min(affinityMatrix[r, l], height[n-i])
      }
    }
  }
  if(kernel == 'Gaussian'){
    affinityMatrix <- exp(-affinityMatrix**2 / (aveheight**2))
  }  
  diag(affinityMatrix) <- 0
  return(affinityMatrix)
}


#'  Determining Appropriate Clusters for HCV Objects
#'
#' The funciton provides two methods to determine an  appropriate number of clusters for an \code{HCV} object, and reports individual cluster members. One of the method is a novel internal index named Spatial Mixture Index (SMI), considering both the within-cluster sum of squared difference of geographical attributes and non-geographical attributes. The other is an M3C-based method taking account of the stability of clusters. 
#@param  geometry_domain \emph{n} by \emph{d} matrix (NA not allowed) of geographical coordinates for \emph{n} points or the centers of \emph{n} polygons in \emph{d}-dimensional space.
#@param  feature_domain \emph{n} by \emph{p} matrix (NA allowed) for \emph{n} samples with \emph{p} attributes on the feature domain.
#' @param  HCVobj an object resulting from calling the \code{HCV} function.
#' @param method character indicating the method to determine an appropriate number of clusters. Default 'SMI' is faster, while 'M3C' is more precise but slower.
#' @param Kmax integer for the upper bound of the potential number of clusters to be considered.
#' @param niter integer for the number of resampling, only used in \code{method='M3C'}.
#' @param criterion character indicating whether to use 'PAC' or 'entropy' as the objective function. Default is 'PAC'. Only used in \code{method='M3C'}. See the reference for details.
#' @export
#' @seealso \code{\link{M3C}}
#' @author ShengLi Tzeng and Hao-Yun Hsu.  
#' @references John, Christopher R., et al. (2020). M3C: Monte Carlo reference-based consensus clustering. Scientific reports, 10(1), 1-14.
#' @examples
#' set.seed(0)
#' pcase  <-  synthetic_data(3,30,0.02,100,2,2)
#' HCVobj <- HCV(pcase$geo,  pcase$feat)
#smi=getCluster(pcase$geo,pcase$feat,HCVobj,method="SMI")
#' smi <- getCluster(HCVobj,method="SMI")
#' par(mfrow=c(2,2))
#' labcolor  <-  (pcase$labels+1)%%3+1
#' plot(pcase$feat,  col  =  labcolor,  pch=19,  xlab  =  'First  attribute', 
#'   ylab  =  'Second  attribute',  main  =  'Feature  domain')
#' plot(pcase$geo,  col  =  labcolor,  pch=19,  xlab  =  'First  attribute', 
#'   ylab  =  'Second  attribute',  main  =  'Geometry  domain')
#' plot(pcase$feat,  col=factor(smi),pch=19,  xlab  =  'First  attribute', 
#'   ylab  =  'Second  attribute',main  =  'Feature  domain')
#' plot(pcase$geo,  col=factor(smi),pch=19,  xlab  =  'First  attribute', 
#'   ylab  =  'Second  attribute',main  =  'Geometry  domain')
#'
getCluster <- function( HCVobj, method=c('SMI','M3C'), Kmax=10, niter=25, criterion='PAC'){
# only M3C needs niter, criterion
 method <- match.arg(method, c("SMI","M3C"))

 if(method=="M3C") {
  mes <- affinityMat(HCVobj, kernel = 'none')
  mes <- cmdscale(mes)
  mes <- t(mes)
  colnames(mes) <- paste0('x',1:ncol(mes))
  mccces <-  M3C::M3C(mes, removeplots = TRUE, iters=niter, maxK=Kmax, objective=criterion, silent=TRUE)
  glabel <- FindComponents(as.matrix(HCVobj$adjacency.matrix), mccces$assignments)
  grps <- glabel$labels
  return(grps)
 }

 if(method=="SMI") {
  k_find <- SMI(HCVobj, Kmax)$bestCluster
  mes <- affinityMat(HCVobj, kernel = 'none')
  glabel <- cluster::pam(as.dist(mes),k=k_find[[1]])
  return(glabel$clustering)
 }
}

