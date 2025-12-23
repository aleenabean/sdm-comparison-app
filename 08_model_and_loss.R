# SINR-style MLP (sigmoid outputs) and single-positive multi-label loss
# with random background locations

library(torch)

source("paths.R")
check_project_dirs()

# MLP with residual connections, sigmoid outputs

sinr_mlp <- nn_module(
  "sinr_mlp",
  
  initialize = function(input_dim,
                        num_species,
                        hidden_dim = 256L,
                        dropout = 0.2) {
    
    self$fc1 <- nn_linear(input_dim, hidden_dim)
    self$fc2 <- nn_linear(hidden_dim, hidden_dim)
    self$fc3 <- nn_linear(hidden_dim, hidden_dim)
    self$out <- nn_linear(hidden_dim, num_species)
    
    self$dropout_p <- dropout
  },
  
  forward = function(x) {
    # x: [batch_size, input_dim]
    
    h1 <- self$fc1(x)
    h1 <- nnf_relu(h1)
    
    h2 <- self$fc2(h1)
    h2 <- nnf_relu(h2)
    h2 <- h2 + h1  # residual
    
    h3 <- self$fc3(h2)
    h3 <- nnf_relu(h3)
    
    h3 <- nnf_dropout(h3, p = self$dropout_p, training = self$training)
    
    logits <- self$out(h3)         # [B, S]
    probs  <- logits$sigmoid()     # [B, S] in (0,1)
    
    probs
  }
)

# -------------------------------------------------------------------
# 2. SINR-style loss: single-positive multi-label + background
# -------------------------------------------------------------------
#
# sinr_spl_loss:
#   p_data   : [B, S] probabilities at observed locations
#   y        : [B] 1-based species labels (1..S)
#   p_rand   : [B_r, S] probabilities at random background locations (or NULL)
#   lambda_pos : weight for positive term (e.g. 50)
#
# Loss (to minimize) ~ -[ λ Σ log(p_pos) + Σ log(1-p_neg_data) + Σ log(1-p_rand) ]

sinr_spl_loss <- function(p_data,
                          y,
                          p_rand = NULL,
                          lambda_pos = 50) {
  # Ensure float / long
  if (p_data$dtype != torch_float()) {
    p_data <- p_data$to(dtype = torch_float())
  }
  if (y$dtype != torch_long()) {
    y <- y$to(dtype = torch_long())
  }
  
  B <- p_data$size(1)
  S <- p_data$size(2)
  
  # ---- Build targets for data: one-hot (single-positive) ------------
  # y is 1-based; we create a one-hot [B, S] on CPU, then move to device.
  
  y_cpu <- y$to(device = "cpu")
  y_vec <- as.integer(as_array(y_cpu))
  
  if (any(is.na(y_vec))) {
    stop("sinr_spl_loss: label tensor contains NA values.")
  }
  if (any(y_vec < 1L | y_vec > S)) {
    stop(sprintf(
      "sinr_spl_loss: labels out of range [1, %d]. Observed [%d, %d]",
      S, min(y_vec), max(y_vec)
    ))
  }
  
  target_mat <- matrix(0, nrow = B, ncol = S)
  target_mat[cbind(seq_len(B), y_vec)] <- 1
  
  target_data <- torch_tensor(
    target_mat,
    dtype  = torch_float(),
    device = p_data$device
  )  # [B, S]
  
  # background targets (all zeros)
  if (!is.null(p_rand)) {
    if (p_rand$dtype != torch_float()) {
      p_rand <- p_rand$to(dtype = torch_float())
    }
    target_rand <- torch_zeros_like(p_rand)  # [B_r, S]
    
    # Concatenate data + random along batch dimension
    probs_all   <- torch_cat(list(p_data, p_rand), dim = 1L)
    targets_all <- torch_cat(list(target_data, target_rand), dim = 1L)
    
    B_total <- probs_all$size(1)
  } else {
    probs_all   <- p_data
    targets_all <- target_data
    B_total <- B
  }
  
  #compute LAN-like objective
  eps <- 1e-7
  p <- probs_all$clamp(min = eps, max = 1 - eps)
  
  # positive term: only where target == 1
  pos_term <- (targets_all * p$log())$sum()
  
  # negative term: where target == 0
  neg_term <- ((1 - targets_all) * (1 - p)$log())$sum()
  
  total_obj <- lambda_pos * pos_term + neg_term
  
  # loss to minimize
  loss <- -total_obj / B  # normalize by number of positive samples
  
  loss
}
