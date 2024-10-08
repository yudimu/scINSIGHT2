#' Selection of latent dimension p and seed d
#'
#' @description
#' Select latent dimension p and the corresponding initialization seed d to finalize the estimation
#'
#' @param U_all The estimation of latent factors U with each candidate p and seed pair.
#' @param individual A vector indicating different samples.
#' @param n_cell Total number of cells in the count matrix.
#' @param p_candidate The candidate p values.
#' @param seeds Candidate initialization seeds.
#'
#' @return A list of selected latent dimension p and seed.


selection = function(U_all, individual, n_cell, p_candidate, seeds) {
  print(table(individual))
  print(seeds)
  # Initialize number of seeds
  n_res = length(seeds)

  # Helper function to convert clustering assignments to connectivity matrix
  clust2connect = function(assignment) {
    n = length(assignment)
    mat = matrix(0, nrow = n, ncol = n)
    for (l in unique(assignment)) {
      id = which(assignment == l)
      mat[id, id] = 1
    }
    return(mat)
  }

  consmat = list()

  # Loop through each candidate p value
  for (j in 1:length(p_candidate)) {
    # Generate connectivity matrices for each seed
    consmat_all = lapply(seeds, function(k) {
      name = as.character(paste0(k, p_candidate[j]))
      pick = names(U_all)[grep(name, names(U_all))]
      u_pick = U_all[pick]

      new = data.frame(u_pick[[1]], individual = individual)
      new1 = new %>% group_split(individual, .keep = F)
      new = list()

      for (i in 1:length(unique(individual))) {
        new[[i]] = unname(as.matrix(as.data.frame(new1[i][[1]])))
      }

      clust = norm_clust_strict_weighted(new)
      assign = unlist(clust$clusters)

      if (n_cell > 20000) {
        n_sample = min(2^16 - 1, floor(n_cell * 0.2))
        index = sample(1:n_cell, n_sample)
        assign = assign[index]
      }

      return(clust2connect(assign) / n_res)
    })

    names(consmat_all) = seeds
    assign(paste0('consmat_all_p', p_candidate[j]), consmat_all)
    consmat[[j]] = Reduce('+', consmat_all)
  }

  names(consmat) = p_candidate
  saveRDS(consmat, paste0(out.dir, 'cosmat.rds'))


  entropy = c()

  # Calculate entropy for each p candidate
  for (i in p_candidate) {
    mat = consmat[[paste(i)]]
    entropy = c(entropy, -sum((mat * log(mat, 2)), na.rm = T))
  }

  # Select the best p value based on entropy
  if (length(p_candidate) <= 2) {
    p_final = p_candidate[which.min(entropy)]
  } else {

    selected_numbers = numeric(0)
    index = numeric(0)

    # Loop to find local minima in entropy values
    for (i in 2:(length(entropy) - 1)) {
      if (entropy[i] < entropy[i - 1] && entropy[i] < entropy[i + 1]) {
        selected_numbers = c(selected_numbers, entropy[i])
        index = c(index, i)
      }
    }

    if (length(selected_numbers) == 0) {
      p_final = p_candidate[which.min(entropy)]
    } else {
      p_final = p_candidate[min(index)]
    }

  }

  print(entropy)

  # Select the best seed based on Frobenius norm difference
  if (length(seeds) == 1) {
    seed_final = seeds
  } else {
    consmat_pfinal = lapply(seeds, function(k){
      name = as.character(paste0(k, p_final))
      pick = names(U_all)[grep(name, names(U_all))]
      u_pick = U_all[pick]

      new = data.frame(u_pick[[1]], individual = individual)
      new1 = new %>% group_split(individual, .keep = F)
      new = list()

      for (i in 1:length(unique(individual))){
        new[[i]] = unname(as.matrix(as.data.frame(new1[i][[1]])))
      }

      clust = norm_clust_strict_weighted(new)
      assign = unlist(clust$clusters)

      if(n_cell > 20000){
        n_sample = min(2^16-1, floor(n_cell*0.2))
        index = sample(1:n_cell, n_sample)
        assign = assign[index]
      }
      return(clust2connect(assign)/n_res)
    })


    consmat_test = Reduce('+', consmat_pfinal)

    frobenius_norm_difference = sapply(seeds, function(x)  norm(consmat_pfinal[[x]]-consmat_test,  type = "F"))
    print('frobenius_norm_difference:')
    print(frobenius_norm_difference)
    seed_final = seeds[which.min(frobenius_norm_difference)]
  }

  return(list(p_final = p_final, seed_final = seed_final))
}
