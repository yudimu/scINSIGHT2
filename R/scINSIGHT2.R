#' Main unction in scINSIGHT2 package
#'
#' @description
#' To estimate latent factors, corresponding loadings and coefficients of covariates.
#'
#' @param object A Seurat object containing count matrix Y and covariate information X in meta.data.
#' @param B Number of different initialization seeds to be selected from. The candidates are 1, ..., B. (default 1)
#' @param libsize library size of all the cells if provided.
#' @param p_candidate The vector of candidate p values.
#' @param maxIter The number of maximum iteration. (default 5000)
#' @param alpha Step size.(default 0.01)
#' @param tol1 The stopping rule for RMSE. (default 1e-5)
#' @param tol2 The stopping rule for mean RMSE. (default 1.5e-5)
#' @param num.cores Number of cores used for estimation. (default 1)
#' @param out.dir Output directory.
#'
#' @useDynLib scINSIGHT2
#' @return A list of estimates.
#' @import Rcpp
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat SplitObject
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat SelectIntegrationFeatures
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @importFrom rlist list.rbind
#' @importFrom parallel mclapply
#' @importFrom fastDummies dummy_cols
#' @importFrom stats median
#' @importFrom stringr str_sub
#'
#' @export

scINSIGHT2_estimate = function(object,
                               libsize=NULL,
                               p_candidate=NULL,
                               B=1,
                               maxIter=5000,
                               alpha=0.01,
                               tol1=1e-5,
                               tol2=1.5e-5,
                               num.cores = 1,
                               out.dir=NULL)
{

  #Create out.dir
  if(!is.null(out.dir)){

    if (!dir.exists(out.dir)) {
      # Create the directory if it does not exist
      if(str_sub(out.dir, -1)!="/"){
        out.dir = paste0(out.dir,"/")
      }
      dir.create(out.dir, recursive = TRUE)

    } else {
      message("Directory already exists: ", out.dir)
    }
  }

  #track time
  start_time = Sys.time()

  #input setting
  if (is.null(libsize)){
    #Calculate the log library size factor
    libsize_Y = rowSums(Y)
    s_Y = libsize_Y/median(libsize_Y)
    logs = log(s_Y)
  }else{
    logs = log(libsize/median(libsize))
  }

  if(is.null(libsize)){
    #Find highly variable genes
    data.list <- SplitObject(object, split.by = "orig.ident")
    # normalize and identify variable features for each dataset independently
    data.list <- lapply(X = data.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })

    features_select <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 2000)

    #Subset the original subject to only 2000 genes
    Seurat_obj = subset(object, features = features_select)
  }else{
    Seurat_obj = object
  }

  #Input count matrix, covariates, individual index and log library size factor
  Y = t(as.matrix(Seurat_obj@assays$RNA$counts))
  X = Seurat_obj@meta.data
  X = X[,!(colnames(X) %in% c("orig.ident", "nCount_RNA", "nFeature_RNA"))]
  individual = Seurat_obj$orig.ident

  #set dummy variables for individual
  ind = data.frame(ind=as.factor(individual))
  ind = dummy_cols(ind, select_columns = 'ind', remove_first_dummy = T)[,-1]
  X = cbind(X, ind)
  X = as.matrix(X)


  seeds = 1:B
  pairs = expand.grid(seeds, p_candidate)
  colnames(pairs) = c('seed', 'p')
  pair_list = split(pairs, seq(nrow(pairs)))


  output_seed = mclapply(1:length(pair_list), function(i){
    p = pair_list[[i]]$p
    seed = pair_list[[i]]$seed
    print(paste0('Run model p=', p, ' seed=', seed))

    res = estimation(Y,
                     X,
                     logs=logs,
                     p=pair_list[[i]]$p,
                     seed=seed,
                     maxIter=maxIter,
                     alpha=alpha,
                     tol1=tol1,
                     tol2=tol2,
                     out.dir = out.dir)
    gc()
    return(res)
  }, mc.cores = num.cores, mc.preschedule = T)

  names(output_seed) = pair_list
  names(output_seed) = gsub("[^0-9]+", "", names(output_seed))

  #save U
  U_all = lapply(output_seed, function(df){ df[['U']]})
  saveRDS(U_all, file = paste0(out.dir, '_U_all.rds'))

  #select hyper-parameters
  n_cell = nrow(Y)
  selected = selection(U_all = U_all, individual = individual, n_cell = n_cell, p_candidate = p_candidate, seeds = seeds)
  p_final = selected$p_final
  seed_final = selected$seed_final

  print(paste0('The selected p is ', p_final, ' and the corresponding optimal seed is ', seed_final, '.'))

  #save the final output
  final = readRDS(paste0(out.dir, 'Results_p=', p_final, '_seed=', seed_final, '.rds'))

  #normalize U
  gllvm_u = final$U

  new = data.frame(u = gllvm_u, individual = individual)
  new1 = new %>% group_split(individual, .keep = F)

  U_new = list()

  for (i in 1:length(unique(individual))){
    U_new[[i]] = unname(as.matrix(as.data.frame(new1[i][[1]])))
  }

  #normalization on U and clustering
  norm_clust = norm_clust_strict_weighted(U_new)
  U_normalized = list.rbind(norm_clust$U)
  final_clusters = norm_clust$clusters

  #end time
  end_time = Sys.time()


  #save final list
  final$U_norm = U_normalized
  final$time = end_time - start_time
  final$clusters = final_clusters
  final$seed = seed_final
  final = final[!names(final) %in% c("beta", "likelihood")]


  return(final)
}























