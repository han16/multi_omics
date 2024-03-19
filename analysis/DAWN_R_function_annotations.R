#translate into synonyms
#>covarianceSelection::symbol_synonyms
function (vec, verbose = T) 
{
  dbCon <- org.Hs.eg.db::org.Hs.eg_dbconn()
  sqlQuery <- "SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;"
  aliasSymbol <- DBI::dbGetQuery(dbCon, sqlQuery) # Send query, retrieve results and then clear result set
  sapply(1:length(vec), function(i) {
    if (verbose & i%%max(floor(length(vec)/10), 1) == 0) 
      cat("*")
    res <- aliasSymbol[which(aliasSymbol[, 2] %in% vec[i]), 
                       5]
    if (length(res) == 0) # no matched gene 
      return(NA)     
    if (length(res) > 1) {
      len_vec <- sapply(res, function(x) {
        length(which(aliasSymbol[, 2] %in% x))
      })
      res_final <- res[which(len_vec == max(len_vec))]
      if (length(res_final) > 1) {
        res_final <- sort(res_final, decreasing = F)[1]
      }
      res <- res_final
    }
    res
  })
}

#######################################################################################
#average non-unique genes
#> covarianceSelection:::average_same_columns
function (dat) 
{
  col_vec <- colnames(dat)
  tab_vec <- table(col_vec)
  idx <- which(tab_vec > 1)  # find non-unique genes, i.e. genes with multiple appearances 
  remove_idx <- numeric(0)
  if (length(idx) > 0) {
    for (i in idx) {
      col_idx <- which(col_vec == names(tab_vec)[i])
      dat[, col_idx[1]] <- rowMeans(dat[, col_idx])
      remove_idx <- c(remove_idx, col_idx[-1])
    }
    dat <- dat[, -remove_idx]
  }
  dat
}

#################################
#' Matching
#'
#' Outputs a result so that \code{vec1 == vec2[matching(vec1, vec2)]}.
#'
#' @param vec1 numeric vector
#' @param vec2 numeric vector
#'
#' @return vector of indices
#' @export
matching <- function(vec1, vec2){
  stopifnot(length(vec1) == length(vec2))
  stopifnot(length(vec1) == length(unique(vec1)))
  stopifnot(length(vec2) == length(unique(vec2)))
  stopifnot(is.vector(vec1), is.vector(vec2))
  stopifnot(all(vec2 %in% vec1), all(vec1 %in% vec2))
  
  order(vec2)[rank(vec1)]
}
############################


#' Gene synonyms
#' 
#' This function is an in-house way to map genes to their most commonly-known form
#' based on the \code{org.Hs.eg.db} package. If there are conflicts (i.e., multiple
#' genes that claim to be the official alias of a gene), we take the most "common" gene (i.e.,
#' how often it appears as the official gene). If there are still ties, we take the 
#' first gene alphabetically.
#' 
#' If no genes are found for a particular gene character, return \code{character(0)}.
#' 
#' At the time of this function's writing, we are using version 3.8.2 of the \code{org.Hs.eg.db} package.
#' 
#' @param vec a character vector
#' @param verbose a boolean
#' 
#' @return a character vec, in order with respect to the input \code{vec}
#' @export
#' 
#' @source The code for this comes from the user Duff at \url{https://www.biostars.org/p/14971/}.
symbol_synonyms <- function(vec, verbose = T){
  dbCon <- org.Hs.eg.db::org.Hs.eg_dbconn()
  sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
  aliasSymbol <- DBI::dbGetQuery(dbCon, sqlQuery)
  
  sapply(1:length(vec), function(i){
    if(verbose & i %% max(floor(length(vec)/10),1) == 0) cat('*')
    
    res <- aliasSymbol[which(aliasSymbol[,2] %in% vec[i]), 5]
    if(length(res) == 0) return(NA)
    
    #if there are more than one, take the most common one
    if(length(res) > 1){
      len_vec <- sapply(res, function(x){
        length(which(aliasSymbol[,2] %in% x))
      })
      
      #if there is still a tie, take the first alphabetic one
      res_final <- res[which(len_vec == max(len_vec))]
      if(length(res_final) > 1){
        res_final <- sort(res_final, decreasing = F)[1]
      }
      
      res <- res_final
    }
    
    res
  })
}

average_same_columns <- function(dat){
  col_vec <- colnames(dat)
  tab_vec <- table(col_vec)
  idx <- which(tab_vec > 1)
  remove_idx <- numeric(0)
  
  if(length(idx) > 0){
    for(i in idx){
      col_idx <- which(col_vec == names(tab_vec)[i])
      dat[,col_idx[1]] <- rowMeans(dat[,col_idx])
      remove_idx <- c(remove_idx, col_idx[-1])
    }
    
    dat <- dat[,-remove_idx]
  }
  
  dat
}

####################

#' Binning of subjects
#'
#' This must be the in the format of ID.REGION.(whatever), for example,
#' "HSB100.PFC.6", where the ID is prefixed with "HSB". The splits are
#' a vector of integers between 1 and 15 in increasing order.
#'
#' @param vec a vector of valid subjects
#' @param splits a vector of integers
#'
#' @return a table
#' @export
binning <- function(vec, splits = c(2,5,8,15)){
  stopifnot(all(diff(splits) > 0), min(splits) >= 1, max(splits) <= 15)
  
  name_mat <- sapply(vec, .split_name)
  time <- .time_to_splits(as.numeric(name_mat[3,]), splits = splits)
  
  region <- factor(as.character(name_mat[2,]), levels = c("PFC", "VIIAS", "SHA", "MDCBC"))
  
  table(region, time)
}

.time_to_splits <- function(vec, splits = c(2,5,8,15)){
  vec[vec <= splits[1]] <- 1
  for(j in 2:length(splits)){
    vec[intersect(which(vec >= splits[j-1]+1), which(vec <= splits[j]))] <- j
  }
  
  factor(vec, levels = c(1,2,3,4))
}

####################

.split_name <- function(str){
  id <- unlist(strsplit(str,"\\."))[1]
  subregion <- unlist(strsplit(str,"\\."))[2]
  
  dat <- covarianceSelection::region_subregion
  if(subregion %in% dat$region){
    region <- subregion
  } else {
    stopifnot(subregion %in% dat$subregion)
    region <- dat$region[which(dat$subregion == subregion)]
  }
  
  dat <- covarianceSelection::brainspan_id
  stopifnot(id %in% dat$Braincode)
  time <- dat$Stage[which(dat$Braincode == id)]
  
  c(id, region, time)
}

#' Split the dataset into smaller datasets
#'
#' According to the naming convention set in \code{.split_name}
#' and \code{region_subregion}, split the dataset according to the
#' rows The output is a list where each data frame originated
#' from the same ID and same brain region.
#'
#' @param dat data frame
#'
#' @return list of data frames
#' @export
extractor <- function(dat){
  stopifnot(is.data.frame(dat))
  id <- sapply(rownames(dat), .split_name)
  id <- as.factor(apply(id, 2, function(x){paste0(x[1], ".", x[2], ".", x[3])}))
  
  base::split(dat, f = id, drop = F)
}

############################

#' Regroup a list of matrices into another list of matrices
#'
#' \code{dat_list} is a list of matrices, and \code{idx_list} is
#' a list of integer vectors, where each integer from 1 to \code{length(dat_list)}
#' appears exactly once. Outputs a list of matrices of length \code{length(idx_list)}
#' matrices.
#'
#' @param dat_list a list of matrices
#' @param idx_list a list of integer vectors
#'
#' @return a list of matrices
#' @export
regroup <- function(dat_list, idx_list){
  .is.listOfMatrix(dat_list, "dat_list")
  .is.listofNumeric(idx_list, "idx_list")
  .regroup_check(length(dat_list), idx_list)
  
  dat <- do.call(cbind, dat_list)
  idx.vec <- .populate(sapply(dat_list, ncol), idx_list)
  
  lis <- vector("list", length(idx_list))
  
  for(i in 1:length(idx_list)){
    idx <- which(idx.vec == i)
    lis[[i]] <- dat[,idx]
  }
  
  names(lis) <- names(idx_list)
  lis
}

.populate <- function(ncol.vec, idx_list){
  idx.vec <- .partition_renamed(idx_list)
  rep(idx.vec, times = ncol.vec)
}

.partition_renamed <- function(idx_list){
  idx.vec <- unlist(idx_list)
  n <- max(idx.vec)
  
  vec <- rep(0, n)
  
  for(i in 1:length(idx_list)){
    vec[idx_list[[i]]] <- i
  }
  
  vec
}

.regroup_check <- function(len, idx_list){
  vec <- unlist(idx_list)
  
  .is.nonNegInteger(vec)
  if(length(vec) != len) stop(paste("idx_list is missing elements",
                                    "compared to dat_list"))
  if(!all(sort(vec) == 1:len)) stop(paste("idx_list does not have",
                                          "all consecutive integers"))
  
  TRUE
}












#https://github.com/linnykos/covarianceSelection/blob/master/covarianceSelection/R/graphicalModel.R

#' Graphical model estimate
#'
#' Estimated using neighbhorhood selection, cross validation to select lambda
#'
#' @param dat the matrix with \code{n} rows and \code{d} columns
#' @param primary_idx index vector that is a subset of \code{1:ncol(dat)}
#' @param lambda either a character vector (\code{"lambda.1se"} or \code{"lambda.min"}) or a numeric positive scalar
#' @param verbose boolean
#' @param tol numeric
#'
#' @return a list that contains an \code{d} by \code{d} \code{sparseMatrix} encoding the
#' estimated adjacency matrix and a numeric vector \code{lambda_vec}
#' @export
graphicalModel <- function(dat, primary_idx, lambda = "lambda.1se", verbose = F, tol = 1e-6){
  n <- nrow(dat); d <- ncol(dat)
  
  if(verbose) print("Starting to estimate coefficients")
  coef_list <- .compute_reg_coefficients_cv(dat, primary_idx, lambda = lambda, verbose = verbose)
  coef_mat <- sapply(coef_list, function(x){ x$vec })
  lambda_vec <- sapply(coef_list, function(x){ x$lambda })
  if(abs(diff(range(lambda_vec))) <= tol) lambda_vec <- min(lambda_vec)
  stopifnot(nrow(coef_mat) == ncol(dat), ncol(coef_mat) == length(primary_idx))
  
  adj_mat <- cbind(coef_mat, matrix(0, nrow = nrow(coef_mat), ncol = nrow(coef_mat) - ncol(coef_mat)))
  adj_mat <- .symmetrize(adj_mat)
  adj_mat[which(abs(adj_mat) >= tol)] <- 1
  adj_mat[which(abs(adj_mat) <= tol)] <- 0
  adj_mat <- Matrix::Matrix(adj_mat, sparse = T)
  
  list(adj_mat = adj_mat, lambda_vec = lambda_vec)
}

#' Graphical model estimate for a range of lambda values
#'
#' @param dat the matrix with \code{n} rows and \code{d} columns
#' @param primary_idx index vector that is a subset of \code{1:ncol(dat)}
#' @param lambda_min minimum value of \code{lambda} when using \code{graphicalModel}
#' @param lambda_max maximum value of \code{lambda} when using \code{graphicalModel}
#' @param lambda_length number of \code{lambda} values to try (exponential growth)
#' @param verbose boolean
#' @param tol numeric
#'
#' @return a list, each being an output for \code{graphicalModel} for a different value of \code{lambda}
#' @export
#' 
graphicalModel_range <- function(dat,  primary_idx, lambda_min, lambda_max, lambda_length = 15, verbose = F, tol = 1e-6){
  lambda_seq <- seq(lambda_min, lambda_max, length.out = lambda_length)
  
  lapply(lambda_seq, function(x){
    if(verbose) print(x)
    graphicalModel(dat,  primary_idx, lambda = x, verbose = verbose, tol = tol)
  })
}

##################

.compute_reg_coefficients_cv <- function(dat, primary_idx = 1:ncol(dat), lambda = "lambda.1se", verbose = F){
  d <- ncol(dat)
  stopifnot(all(primary_idx %in% 1:d))
  
  func <- function(i){
    if(verbose & i %% floor(length(primary_idx)/10) == 0) cat('*')
    
    vec <- rep(0, d)
    y_vec <- dat[,i]
    
    if(i %in% primary_idx){
      x_mat <- dat[,-i]; idx_vec <- c(1:d)[-i]
    } else{
      x_mat <- dat[,primary_idx]; idx_vec <- primary_idx
    }
    
    if(is.numeric(lambda)){
      res <- glmnet::glmnet(x = x_mat, y = y_vec, intercept = F, lambda = lambda)  
      # glmnet: https://glmnet.stanford.edu/articles/glmnet.html; use default Gaussian assumption 
      vec[idx_vec] <- as.numeric(res$beta) # for beta https://www.rdocumentation.org/packages/glmnet/versions/4.1-8/topics/glmnet
    } else {
      # because of argument lambda = "lambda.1se", it always do CV step 
      res <- glmnet::cv.glmnet(x = x_mat, y = y_vec, intercept = F) # cv.glmnet is for cross-validation 
      vec[idx_vec] <- as.numeric(stats::coef(res, s = lambda))[-1] # get model coefficients 
      #lambda.1se is the value of 𝜆 that gives the most regularized model such that the cross-validated error is within one standard error of the minimum
     # https://glmnet.stanford.edu/articles/glmnet.html
      }
    
    list(vec = vec, lambda = ifelse(is.numeric(lambda), lambda, res[[which(names(res) == lambda)]]))
  }
  
  i <- 0 #debugging purposes only
  foreach::"%dopar%"(foreach::foreach(i = primary_idx), func(i))
}
################

.l2norm <- function(vec){
  sqrt(sum((vec)^2))
}
###################

.symmetrize <- function(mat){
  stopifnot(ncol(mat) == nrow(mat))
  (mat + t(mat))/2
}


#####################
#> covarianceSelection::screen
function (dat, pv, p_thres = 0.1, num_genes = 3500) 
{
  stopifnot(ncol(dat) == length(pv))
  stopifnot(is.numeric(dat), is.matrix(dat))
  n <- nrow(dat)
  d <- ncol(dat)
  cor_mat <- abs(stats::cor(dat))
  primary <- which(pv < p_thres)
  candidates <- c(1:d)[-primary]
  if (length(candidates) == 0 | length(primary) >= num_genes) {
    return(list(primary = primary, secondary = NA, cor_thres = NA))
  }
  cor_vec <- apply(cor_mat[candidates, primary], 1, max)
  idx <- order(cor_vec, decreasing = T)[1:(num_genes - length(primary))]
  cor_thres <- min(cor_vec[idx])
  secondary <- candidates[idx]
  stopifnot(length(intersect(primary, secondary)) == 0)
  list(primary = sort(primary), secondary = sort(secondary), 
       cor_thres = cor_thres)
}