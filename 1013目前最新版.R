## ============
##  EM 與 naive 
## ============

start_time <- Sys.time()
options(error = function() { traceback(3) })

## ---- 套件 ----
library(brglm2)
library(dplyr)
library(ggplot2)
library(data.table)
library(purrr)
library(furrr)
library(future)
library(future.apply)
library(writexl)

## ---- 載入函式 ----
source("C:/Users/kimi1/OneDrive/文件/論文/diagnostic_funcs.R")  
source("C:/Users/kimi1/OneDrive/文件/論文/hetero_EM_logistic.R")
source("C:/Users/kimi1/OneDrive/文件/論文/diag_all_naive.R")

## ---- AUC 與 cut_vec ----
AUC_fast <- function(y, s) {  
  y <- as.integer(y); stopifnot(all(y %in% 0:1))
  pos <- s[y == 1]; neg <- s[y == 0]
  n1 <- length(pos); n0 <- length(neg)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[1:n1]) - n1*(n1+1)/2) / (n1*n0)
}

make_cut_vec <- function(p, n_main, contam_n) {  
  k <- p + 1
  c(
    CookD    = 4 / (n_main - contam_n),
    DFFITS   = 2 * sqrt(k / (n_main - contam_n)),
    GDFFITS  = 3 * sqrt(k / (n_main - contam_n)),
    GSDFBETA = 9 * k / (n_main - contam_n - 3*p),
    mCDstar  = 3 / sqrt((n_main - p) / (n_main - contam_n)),
    GCD_GSPR = 1,
    StdR2_CS = 5, StdR2_MF = 5, StdR2_NK = 5,
    GD = 4 / (n_main - p), MD = 4 / (n_main - p)
  )
}

## ---- RC-like 起始值 ----
get_rc_init <- function(Y, W, sigma2_mat, use_brglm2 = TRUE,
                        max_abs_beta = 50, verbose = TRUE) {
  stopifnot(is.matrix(W), is.matrix(sigma2_mat),
            all(dim(W) == dim(sigma2_mat)), length(Y) == nrow(W))
  nmW <- colnames(W)
  nzv <- which(apply(W, 2, function(z) sd(z, na.rm=TRUE) < 1e-8))
  if (length(nzv)) {
    if (verbose) message("[RC-init] drop near-zero var: ", paste(nmW[nzv], collapse=", "))
    W <- W[, -nzv, drop=FALSE]
    sigma2_mat <- sigma2_mat[, -nzv, drop=FALSE]
    nmW <- colnames(W)
    if (ncol(W) == 0L) {
      mu0 <- pmin(pmax(mean(Y), 1e-6), 1 - 1e-6)
      return(list(beta = c("(Intercept)" = qlogis(mu0)), type = "SAFE"))
    }
  }
  df  <- data.frame(Y, W)
  frm <- as.formula(paste("Y ~", paste(nmW, collapse=" + ")))
  base_fit <- NULL
  if (use_brglm2 && requireNamespace("brglm2", quietly = TRUE)) {
    base_fit <- try(brglm2::brglm(frm, family = binomial("logit"), data = df, type = "AS_mixed"), silent = TRUE)
  }
  if (inherits(base_fit, "try-error") || is.null(base_fit) || any(!is.finite(coef(base_fit)))) {
    base_fit <- try(glm(frm, family = binomial(), data = df, control = glm.control(maxit = 200, epsilon = 1e-8)), silent = TRUE)
  }
  if (inherits(base_fit, "try-error") || any(!is.finite(coef(base_fit)))) {
    if (verbose) message("[RC-init] SAFE (base fit failed)")
    mu0 <- pmin(pmax(mean(Y), 1e-6), 1 - 1e-6)
    return(list(beta = c("(Intercept)" = qlogis(mu0), setNames(rep(0, ncol(W)), nmW)), type = "SAFE"))
  }
  b0 <- coef(base_fit)[c("(Intercept)", nmW)]
  varW <- pmax(apply(W, 2, var), 1e-12)
  ed   <- pmax(colMeans(sigma2_mat), 0)
  kappa_hat <- pmax((varW - ed) / varW, 1e-3)
  b0[nmW] <- b0[nmW] / kappa_hat[nmW]
  if (any(abs(b0) > max_abs_beta)) {
    if (verbose) message("[RC-init] clip |beta| to ", max_abs_beta)
    b0 <- pmin(pmax(b0, -max_abs_beta), max_abs_beta)
  }
  list(beta = b0, type = "RC_LIKE")
}

check_nonclassical_error <- function(shared, beta_true, tag = "RAW", idx = NULL) {
  # idx=NULL 表示用全體 main；否則用 keep 子集檢查
  Y <- shared$main_dat$Y
  X <- shared$X_main
  W <- shared$W_main %||% as.matrix(shared$main_dat[, grep("^W\\d+$", names(shared$main_dat)), drop=FALSE])
  if (!is.null(idx)) { Y <- Y[idx]; X <- X[idx, , drop=FALSE]; W <- W[idx, , drop=FALSE] }
  
  Delta <- as.data.frame(W - X)
  
  # (1) Δ 與 Y 關聯（古典：兩組均值差 ~ 0）
  by_Y <- sapply(Delta, function(d) {
    est <- try(t.test(d ~ Y)$estimate, silent = TRUE)
    if (inherits(est, "try-error")) return(NA_real_)
    unname(diff(est))  # mean in group 2 - group 1
  })
  
  # (2) Δ 與 X 同欄相關（古典：corr ~ 0）
  cor_DX <- sapply(seq_along(Delta), function(j) suppressWarnings(cor(Delta[[j]], X[, j])))
  
  # (3) Δ 與 η 相關（古典：corr ~ 0）
  eta_true <- as.numeric(beta_true[1] + X %*% beta_true[-1])
  cor_Deta <- sapply(Delta, function(d) suppressWarnings(cor(d, eta_true)))
  
  cat(sprintf("\n[check Δ=W−X] tag=%s | n=%d | p=%d\n", tag, nrow(X), ncol(X)))
  print(round(by_Y,   3))
  print(round(cor_DX, 3))
  print(round(cor_Deta, 3))
  invisible(list(by_Y = by_Y, cor_DX = cor_DX, cor_Deta = cor_Deta))
}
`%||%` <- function(a,b) if (!is.null(a)) a else b

## ========== 共享資料 + 同步汙染 ==========
make_shared_data <- function(kappa_tar, n_total, n_val,
                             p, beta_true, contam_n, sigma_fac,
                             r_target = 0.5,
                             q_prop = 1, tau_dir = 0, mix_prob = 0,
                             bound_fac = 10,
                             flip_prop = 0,
                             seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  ## 1) 生成 X, W, Δ²
  rho <- sqrt(r_target)
  Z   <- matrix(rnorm(n_total * (p + 1)), n_total, p + 1)
  X   <- sqrt(1 - rho^2) * Z[, 1:p] + rho * Z[, p + 1]
  colnames(X) <- paste0("X", 1:p)
  
  sigma_X2   <- apply(X, 2, var)
  bar_delta2 <- sigma_X2 * (1 / kappa_tar - 1)
  
  delta2_mat <- matrix(0, n_total, p)
  for (j in 1:p) delta2_mat[, j] <- runif(n_total, 0.5 * bar_delta2[j], 1.5 * bar_delta2[j])
  
  W <- X + matrix(rnorm(n_total * p, sd = sqrt(as.vector(delta2_mat))), n_total, p)
  colnames(W) <- paste0("W", 1:p)
  
  ## 2) 產生 Y（用 X 與 beta_true）
  eta <- beta_true[1] + X %*% beta_true[-1]
  Y   <- rbinom(n_total, 1, plogis(eta))
  
  ## 3) 主樣本
  full_dat <- data.frame(Y, W)
  main_dat <- full_dat[(n_val + 1):n_total, ]
  n_main   <- n_total - n_val
  
  ## 4) 共同汙染
  # 首先，確定總汙染點的索引
  contam_idx_total <- sample(seq_len(n_main), contam_n)
  
  # 根據 flip_prop，將總汙染點分為兩組
  n_flip <- round(contam_n * flip_prop)
  n_shift <- contam_n - n_flip
  
  if (n_flip > 0) {
    contam_idx_flip <- sample(contam_idx_total, n_flip)
  } else {
    contam_idx_flip <- integer(0)
  }
  
  if (n_shift > 0) {
    contam_idx_shift <- setdiff(contam_idx_total, contam_idx_flip)
  } else {
    contam_idx_shift <- integer(0)
  }
  
  # 取得 W 矩陣和 Y 向量
  Wmat <- as.matrix(main_dat[, paste0("W", 1:p)])
  
  # =======================================================
  # (A) 執行 Y-flipping 汙染 (只對 contam_idx_flip 操作)
  # =======================================================
  if (n_flip > 0) {
    main_dat$Y[contam_idx_flip] <- 1 - main_dat$Y[contam_idx_flip]
  }
  
  # =======================================================
  # (B) 執行 W-shifting 汙染 (只對 contam_idx_shift 操作)
  # =======================================================
  if (n_shift > 0) {
    beta_s  <- beta_true[-1]
    b0_true <- beta_true[1]
    eta0  <- as.numeric(b0_true + Wmat %*% beta_s)
    s_eta <- sd(eta0)
    k <- sigma_fac
    B <- bound_fac * s_eta
    
    # 注意：迴圈只跑在 shift 的索引上
    for (i in contam_idx_shift) {
      sgn_base <- if (main_dat$Y[i] == 0) +1 else -1
      sgn <- if (runif(1) < mix_prob) sample(c(-1,+1), 1) else sgn_base
      d_move_nom <- sgn * (k * s_eta)
      eta_i <- eta0[i]
      d_max <- sgn * (B - sgn * eta_i)
      d_move <- sgn * min(abs(d_move_nom), abs(d_max))
      
      q <- max(1, floor(q_prop * p))
      j_sub <- sample.int(p, q, replace = FALSE)
      b_sub <- beta_s[j_sub]
      b_norm_sub <- sqrt(sum(b_sub^2)); if (b_norm_sub < 1e-12) next
      u_sub <- b_sub / b_norm_sub
      if (tau_dir > 0) {
        z <- rnorm(length(j_sub)); u_sub <- u_sub + tau_dir*z; u_sub <- u_sub / sqrt(sum(u_sub^2))
      }
      t_step <- d_move / b_norm_sub
      Wmat[i, j_sub] <- Wmat[i, j_sub] + t_step * u_sub
    }
  } # end of n_shift > 0
  
  # 最後，將可能被修改過的 W 矩陣寫回 main_dat
  for (j in seq_len(p)) main_dat[[paste0("W", j)]] <- Wmat[, j]
  
  # 返回的 contam_idx 應為所有汙染點的總和
  good_idx <- setdiff(seq_len(n_main), contam_idx_total)
  list(
    main_dat    = main_dat,
    delta2_main = delta2_mat[(n_val + 1):n_total, , drop=FALSE],
    contam_idx  = contam_idx_total,
    G.ind       = good_idx,
    n_main      = n_main,
    X_main      = X[(n_val + 1):n_total, , drop = FALSE],
    W_main      = W[(n_val + 1):n_total, , drop = FALSE],
    contam_idx_flip  = contam_idx_flip,
    contam_idx_shift = contam_idx_shift
  )
}

## ========== 用同一份資料跑 EM ==========
run_EM_from_shared <- function(shared, p, beta_true, kappa_tar, cut_vec,
                               em_diag_func = diag_all_1stepEM,
                               verbose = TRUE){
  
  main_dat    <- shared$main_dat
  delta2_main <- shared$delta2_main
  contam_idx  <- shared$contam_idx
  G.ind       <- shared$G.ind
  n_main      <- shared$n_main
  
  if (length(unique(main_dat$Y)) < 2L) return(NULL)
  p1 <- mean(main_dat$Y[G.ind] == 1)
  if (p1 <= 0.1 || p1 >= 0.9) return(NULL)
  
  ## --- EM（RC 起始）
  Wmat <- as.matrix(main_dat[, paste0("W", 1:p)])
  ini  <- get_rc_init(main_dat$Y, Wmat, sigma2_mat = delta2_main, verbose = verbose)
  beta_init <- ini$beta
  init_type <- ini$type
  
  em_out <- hetero_em_logit(
    Y = main_dat$Y, W = Wmat, sigma2_mat = delta2_main,
    init_beta = beta_init,
    B0 = 120, Bmax = 800, growth = 1.6, B_vcov = 1200,
    maxit = 300, tol = 1e-4, verbose = verbose, vcov = "sandwich",
    seed = 123, use_sobol = TRUE, antithetic = TRUE
  )
  
  ## --- 多輪小步修剪 ---
  .z <- function(x) { m <- stats::median(x, na.rm=TRUE); s <- stats::mad(x, center=m, constant=1.4826, na.rm=TRUE); if (!is.finite(s) || s<=0) return(rep(0,length(x))); (x-m)/s }
  keep <- seq_len(n_main); dropped <- integer(0)
  trim_budget_frac <- 0.05; trim_budget <- floor(trim_budget_frac * n_main)
  max_round <- 2
  em_curr <- em_out
  
  for (iter in seq_len(max_round)) {
    
    ## A) 針對當前 keep 子樣本，重算這一輪的 cut（門檻）
    contam_idx_iter <- match(intersect(contam_idx, keep), keep)
    cut_vec_iter <- make_cut_vec(
      p = p,
      n_main  = length(keep),
      contam_n = length(contam_idx_iter)
    )
    
    ## B) 這一輪的診斷，用當前樣本 + 當前 cut_vec_iter
    infl_sub <- em_diag_func(
      Y = main_dat$Y[keep],
      W = as.matrix(main_dat[keep, paste0("W", 1:p), drop = FALSE]),
      sigma2_mat = delta2_main[keep, , drop = FALSE],
      cut_vec = cut_vec_iter,                     # ← 原本是 cut_vec
      G.ind = match(intersect(G.ind, keep), keep),
      em_out = em_curr,
      dbg = FALSE
    )
    
    ## C) 立旗一律用「這一輪」的 cut_vec_iter
    idx_keep <- keep[infl_sub$id]
    flag_rcs <- infl_sub$R_CS      >  cut_vec_iter["StdR2_CS"]
    flag_mdf <- abs(infl_sub$MDF)  >  cut_vec_iter["GDFFITS"]
    flag_mcd <- infl_sub$mCDstar   >  cut_vec_iter["mCDstar"]
    hit_cnt  <- flag_rcs + flag_mdf + flag_mcd
    cand_local <- which(hit_cnt >= 3L)
    if (!length(cand_local)) break
    
    score_local <- .z(infl_sub$R_CS) + .z(abs(infl_sub$MDF)) + .z(infl_sub$mCDstar)
    ord_local   <- order(score_local[cand_local], decreasing = TRUE)
    
    remain <- trim_budget - length(dropped)
    if (remain <= 0L) break
    drop_n <- min(remain, max(0L, floor(0.5 * length(cand_local))))
    if (drop_n <= 0L) break
    
    cand_global <- idx_keep[cand_local]
    to_drop     <- cand_global[ord_local][seq_len(drop_n)]
    
    # 3) 模擬修剪，再做守門檢查
    keep_after <- setdiff(keep, to_drop)
    yk_after   <- main_dat$Y[keep_after]
    p1k_after  <- mean(yk_after == 1)
    
    if (length(unique(yk_after)) < 2L || p1k_after <= 0.05 || p1k_after >= 0.95) {
      # 不通過：不要修剪，直接結束修剪流程
      break
    }
    
    # 4) 通過 → 正式修剪並重估 em_curr
    keep    <- keep_after
    dropped <- c(dropped, to_drop)
    
    SigmaX_keep <- .est_SigmaX_full(
      W = as.matrix(main_dat[keep, paste0("W", 1:p), drop = FALSE]),
      sigma2_mat = delta2_main[keep, , drop = FALSE]
    )
    ini2 <- get_rc_init(
      Y = main_dat$Y[keep],
      W = as.matrix(main_dat[keep, paste0("W", 1:p), drop = FALSE]),
      sigma2_mat = delta2_main[keep, , drop = FALSE],
      verbose = FALSE
    )
    em_curr <- hetero_em_logit(
      Y = main_dat$Y[keep],
      W = as.matrix(main_dat[keep, paste0("W", 1:p)]),
      sigma2_mat = delta2_main[keep, , drop = FALSE],
      init_beta = ini2$beta,
      B0 = 150, Bmax = 600, growth = 1.5, B_vcov = 1000,
      maxit = 250, tol = 1e-4, vcov = "sandwich", verbose = FALSE,
      SigmaX = SigmaX_keep, use_sobol = TRUE, antithetic = TRUE
    )
  }
  em_out <- em_curr
  # 修剪完（取得 keep）後，加這行看看 keep 上 Δ 是否更接近古典
  if (exists("check_nonclassical_error", inherits = TRUE)) {
    check_nonclassical_error(shared, beta_true, tag = "EM_KEEP", idx = keep)
  }
  
  # 1) 映射到 keep 的座標 + 新 n
  n_main_post <- length(keep)
  contam_idx_post <- match(intersect(contam_idx, keep), keep)
  
  # 2) 用新 n / 新 contam_n 重算 cut（關鍵）
  cut_vec_post <- make_cut_vec(p = p, n_main = n_main_post, contam_n = length(contam_idx_post))
  
  # === 修剪結束後：用最終 keep + 最終 em_out 重算一次診斷 ===
  infl <- em_diag_func(
    Y = main_dat$Y[keep],
    W = as.matrix(main_dat[keep, paste0("W", 1:p), drop = FALSE]),
    sigma2_mat = delta2_main[keep, , drop = FALSE],
    cut_vec = cut_vec_post,
    G.ind = match(intersect(G.ind, keep), keep),
    em_out = em_out,
    dbg = FALSE
  )
  n_main    <- n_main_post
  contam_idx <- contam_idx_post
  
  ## 只取斜率並做效能統計
  nm <- colnames(Wmat)
  b <- tryCatch(em_out$beta, error = function(e) NULL)
  if (is.null(b) || any(!is.finite(b)) || max(abs(b)) > 15) {
    beta_EM <- rep(NA_real_, length(nm))
  } else {
    slope_names <- setdiff(names(b), "(Intercept)")
    beta_EM <- as.numeric(b[slope_names][match(nm, slope_names)])
  }
  get_bias_mse <- function(est, beta_true) {
    d <- est - beta_true
    safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    c(bias = safe_mean(d), mse = safe_mean(d^2))
  }
  perf_df <- data.frame(method="EM_trim", t(get_bias_mse(beta_EM, beta_true[-1])))
  perf_df$kappa <- kappa_tar
  perf_df$init_type <- "RC_LIKE"
  perf_df$estimator <- "TrimEM"
  
  ## 指標與 AUC（你的寫法）
  if (!identical(infl$id, seq_len(n_main))) infl <- infl[order(infl$id), , drop = FALSE]
  flag_mat <- as.matrix(infl[, grepl("^flag_", names(infl)), drop = FALSE]); flag_mat[is.na(flag_mat)] <- FALSE
  is_bad <- seq_len(n_main) %in% contam_idx
  TP <- colSums(flag_mat[ is_bad, , drop = FALSE]); FP <- colSums(flag_mat[!is_bad, , drop = FALSE])
  CIR <- TP / sum(is_bad); SR  <- FP / (n_main - sum(is_bad))
  
  labels <- as.integer(is_bad)
  scores <- list(
    CookD = infl$CookD, DFFITS = abs(infl$DFFITS),
    GDF = abs(infl$GDF), MDF = abs(infl$MDF),
    GDB = abs(infl$GDB), MDB = abs(infl$MDB),
    GD = infl$GD, MD = infl$MD,
    GCD_GSPR = abs(infl$GCD_GSPR), mCDstar = infl$mCDstar,
    R_CS = infl$R_CS_raw, R_MF = infl$R_MF_raw, R_NK = infl$R_NK_raw
  )
  auc_by_method <- vapply(scores, function(s) {
    ok <- is.finite(s) & !is.na(labels)
    if (sum(ok) < 2 || length(unique(labels[ok])) < 2 || length(unique(s[ok])) < 2) return(NA_real_)
    AUC_fast(labels[ok], s[ok])
  }, numeric(1))
  
  diag_df <- data.frame(method = sub("^flag_", "", colnames(flag_mat)), CIR = CIR, SR = SR, kappa = kappa_tar)
  diag_df$init_type <- "RC_LIKE"
  diag_df$estimator <- "TrimEM"
  diag_df$AUC <- as.numeric(auc_by_method[match(diag_df$method, names(auc_by_method))])
  
  list(perf = perf_df, diag = diag_df, infl = infl, Q_trace = em_out$Q_trace, keep = keep)
}

run_EM_rawdiag_from_shared <- function(shared, p, beta_true, kappa_tar, cut_vec,
                                       em_diag_func = diag_all_1stepEM, verbose = FALSE) {
  main_dat    <- shared$main_dat
  delta2_main <- shared$delta2_main
  n_main      <- shared$n_main
  G.ind       <- shared$G.ind
  
  # 基本檢查
  if (length(unique(main_dat$Y)) < 2L) return(NULL)
  p1 <- mean(main_dat$Y[G.ind] == 1)
  if (p1 <= 0.1 || p1 >= 0.9) return(NULL)
  
  # --- 一次性 EM（不修剪）
  Wmat <- as.matrix(main_dat[, paste0("W", 1:p)])
  ini  <- get_rc_init(main_dat$Y, Wmat, sigma2_mat = delta2_main, verbose = verbose)
  em_out <- hetero_em_logit(
    Y = main_dat$Y, W = Wmat, sigma2_mat = delta2_main,
    init_beta = ini$beta,
    B0 = 120, Bmax = 800, growth = 1.6, B_vcov = 1200,
    maxit = 300, tol = 1e-4, verbose = verbose, vcov = "sandwich",
    seed = 123, use_sobol = TRUE, antithetic = TRUE
  )
  
  # --- 直接診斷（不修剪）
  infl <- em_diag_func(
    Y = main_dat$Y,
    W = Wmat,
    sigma2_mat = delta2_main,
    cut_vec = cut_vec,                    # 用原始 n_main 的 cut
    G.ind = G.ind,
    em_out = em_out,
    dbg = FALSE
  )
  
  # --- AUC/CIR 用真實污染索引（與 run_EM_from_shared 相同邏輯）
  if (!identical(infl$id, seq_len(n_main))) infl <- infl[order(infl$id), , drop = FALSE]
  flag_mat <- as.matrix(infl[, grepl("^flag_", names(infl)), drop = FALSE]); flag_mat[is.na(flag_mat)] <- FALSE
  
  is_bad <- seq_len(n_main) %in% shared$contam_idx   # ← 直接用 contam_idx
  TP <- colSums(flag_mat[ is_bad, , drop = FALSE])
  FP <- colSums(flag_mat[!is_bad, , drop = FALSE])
  CIR <- TP / sum(is_bad)
  SR  <- FP / (n_main - sum(is_bad))
  
  labels <- as.integer(is_bad)
  scores <- list(
    CookD = infl$CookD, DFFITS = abs(infl$DFFITS),
    GDF = abs(infl$GDF), MDF = abs(infl$MDF),
    GDB = abs(infl$GDB), MDB = abs(infl$MDB),
    GD = infl$GD, MD = infl$MD,
    GCD_GSPR = abs(infl$GCD_GSPR), mCDstar = infl$mCDstar,
    R_CS = infl$R_CS_raw, R_MF = infl$R_MF_raw, R_NK = infl$R_NK_raw
  )
  auc_by_method <- vapply(scores, function(s) {
    ok <- is.finite(s) & !is.na(labels)
    if (sum(ok) < 2 || length(unique(labels[ok])) < 2 || length(unique(s[ok])) < 2) return(NA_real_)
    AUC_fast(labels[ok], s[ok])
  }, numeric(1))
  
  diag_df <- data.frame(method = sub("^flag_", "", colnames(flag_mat)),
                        CIR = CIR, SR = SR, kappa = kappa_tar)
  diag_df$init_type <- "RC_LIKE"
  diag_df$estimator <- "EM_raw"              # 標記：EM-不修剪
  diag_df$AUC <- as.numeric(auc_by_method[match(diag_df$method, names(auc_by_method))])
  
  # 偏誤（可選）：用一次性 EM 的斜率
  nm <- colnames(Wmat)
  b  <- tryCatch(em_out$beta, error = function(e) NULL)
  if (is.null(b) || any(!is.finite(b)) || max(abs(b)) > 15) {
    beta_EM <- rep(NA_real_, length(nm))
  } else {
    slope_names <- setdiff(names(b), "(Intercept)")
    beta_EM <- as.numeric(b[slope_names][match(nm, slope_names)])
  }
  get_bias_mse <- function(est, beta_true) {
    d <- est - beta_true
    safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    c(bias = safe_mean(d), mse = safe_mean(d^2))
  }
  perf_df <- data.frame(method = "EM_raw", t(get_bias_mse(beta_EM, beta_true[-1])))
  perf_df$kappa <- kappa_tar
  perf_df$init_type <- "RC_LIKE"
  perf_df$estimator <- "EM_raw"
  
  list(perf = perf_df, diag = diag_df, infl = infl)
}

## ========== 用同一份資料跑 naive ==========
run_naive_from_shared <- function(shared, p, beta_true, kappa_tar, cut_vec, idx_use = NULL) {
  
  # 先從 shared 解包（一定要先做，下面才有 main_dat / contam_idx 可用）
  main_dat    <- shared$main_dat
  contam_idx  <- shared$contam_idx
  n_main      <- shared$n_main
  
  # 若指定子集（例如 EM 修剪後的 keep），再子集並用新 n 重算 cut
  if (!is.null(idx_use)) {
    main_dat <- main_dat[idx_use, , drop = FALSE]
    n_main   <- nrow(main_dat)
    
    # 把舊 contam 索引映射到子集座標
    contam_map <- match(contam_idx, idx_use)
    contam_idx <- contam_map[!is.na(contam_map)]
    
    # 用新 n / 新 contam_n 重算 cut
    cut_vec <- make_cut_vec(p = p, n_main = n_main, contam_n = length(contam_idx))
  }
  
  # 基本檢查（在子集之後）
  if (length(unique(main_dat$Y)) < 2L) return(NULL)
  p1 <- mean(main_dat$Y[setdiff(seq_len(n_main), contam_idx)] == 1)
  if (p1 <= 0.1 || p1 >= 0.9) return(NULL)
  
  # 診斷（naive，不校正測誤）
  infl <- diag_all_naive(
    Y = main_dat$Y,
    W = as.matrix(main_dat[, paste0("W", 1:p)]),
    G.ind = setdiff(seq_len(n_main), contam_idx),
    cut_vec = cut_vec,
    dbg = FALSE
  )
  
  # 配適一個 naive GLM 以做偏誤/MSE（保持你原本流程）
  keep <- seq_len(n_main)
  glm2 <- try(glm(Y ~ ., family = binomial(),
                  data = main_dat[keep, c("Y", paste0("W", 1:p))],
                  control = glm.control(maxit = 200, epsilon = 1e-8)), silent = TRUE)
  if (inherits(glm2, "try-error") || any(!is.finite(coef(glm2)))) {
    if (requireNamespace("brglm2", quietly = TRUE)) {
      glm2 <- try(brglm2::brglm(Y ~ ., data = main_dat[keep, c("Y", paste0("W", 1:p))],
                                family = binomial("logit"), type = "AS_mixed"), silent = TRUE)
    }
  }
  
  nm <- paste0("W", 1:p)
  b <- tryCatch(coef(glm2), error = function(e) NULL)
  if (is.null(b) || any(!is.finite(b)) || max(abs(b)) > 25) {
    beta_naive <- rep(NA_real_, length(nm))
  } else {
    slope_names <- setdiff(names(b), "(Intercept)")
    beta_naive <- as.numeric(b[slope_names][match(nm, slope_names)])
  }
  
  get_bias_mse <- function(est, beta_true) {
    d <- est - beta_true
    safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    c(bias = safe_mean(d), mse = safe_mean(d^2))
  }
  perf_df <- data.frame(method = "GLM_naive", t(get_bias_mse(beta_naive, beta_true[-1])))
  perf_df$kappa <- kappa_tar
  perf_df$init_type <- "NAIVE"
  perf_df$estimator <- "NaiveGLM"
  
  # CIR/SR/AUC（保持你原本的寫法，R 系列用 raw 分數做 AUC）
  if (!identical(infl$id, seq_len(n_main))) infl <- infl[order(infl$id), , drop = FALSE]
  flag_mat <- as.matrix(infl[, grepl("^flag_", names(infl)), drop = FALSE]); flag_mat[is.na(flag_mat)] <- FALSE
  is_bad <- seq_len(n_main) %in% contam_idx
  TP <- colSums(flag_mat[ is_bad, , drop = FALSE]); FP <- colSums(flag_mat[!is_bad, , drop = FALSE])
  CIR <- TP / sum(is_bad); SR  <- FP / (n_main - sum(is_bad))
  
  labels <- as.integer(is_bad)
  scores <- list(
    CookD = infl$CookD, DFFITS = abs(infl$DFFITS),
    GDF = abs(infl$GDF), MDF = abs(infl$MDF),
    GDB = abs(infl$GDB), MDB = abs(infl$MDB),
    GD = infl$GD, MD = infl$MD,
    GCD_GSPR = abs(infl$GCD_GSPR), mCDstar = infl$mCDstar,
    R_CS = infl$R_CS_raw, R_MF = infl$R_MF_raw, R_NK = infl$R_NK_raw
  )
  auc_by_method <- vapply(scores, function(s) {
    ok <- is.finite(s) & !is.na(labels)
    if (sum(ok) < 2 || length(unique(labels[ok])) < 2 || length(unique(s[ok])) < 2) return(NA_real_)
    AUC_fast(labels[ok], s[ok])
  }, numeric(1))
  
  diag_df <- data.frame(method = sub("^flag_", "", colnames(flag_mat)), CIR = CIR, SR = SR, kappa = kappa_tar)
  diag_df$init_type <- "NAIVE"
  diag_df$estimator <- "NaiveGLM"
  diag_df$AUC <- as.numeric(auc_by_method[match(diag_df$method, names(auc_by_method))])
  
  list(perf = perf_df, diag = diag_df, infl = infl)
}

## ========== 全域設定 ==========
p_list <- c(3)
n_list <- c(100, 200)
r_target <- 0.5
kappa_set <- c(0.75, 0.85, 0.95)
contam_rate_list <- c(0.1)
sigma_fac_list <- c(1, 2, 3)
flip_prop_list <- c(0, 0.5, 1.0)
R <- 25   

## ---- 平行策略 ----
plan(multisession, workers = max(1, parallel::detectCores() - 1L))

## ========== 主迴圈：每一回合「先生成一次共享資料」，再各跑 EM 與 naive ==========
total_jobs <- length(p_list) * length(n_list) * length(contam_rate_list) *
  length(sigma_fac_list) * R * length(kappa_set)

set.seed(2025)
result_perf <- list(); result_diag <- list(); result_meta <- list(); id <- 1L

for (p in p_list) for (n_total in n_list) for (contam_rate in contam_rate_list) for (sigma_fac in sigma_fac_list) for (flip_prop in flip_prop_list) {
  contam_n <- floor(n_total * contam_rate)
  n_val <- 0L
  n_main <- n_total - n_val
  beta_true <- c(1, rep(2, p))
  cut_vec   <- make_cut_vec(p, n_main, contam_n)
  
  setting_res <- future_lapply(seq_len(R), function(r) {
    lapply(kappa_set, function(kappa_tar) {
      ## --- 生成同一份資料 + 汙染
      shared <- make_shared_data(
        kappa_tar = kappa_tar, n_total = n_total, n_val = n_val,
        p = p, beta_true = beta_true, contam_n = contam_n,
        sigma_fac = sigma_fac, r_target = r_target, flip_prop = flip_prop,
        ## 統一汙染參數
        q_prop = 1, tau_dir = 0, mix_prob = 0, bound_fac = 10,
        seed = NULL 
      )
      if (r == 1) check_nonclassical_error(shared, beta_true, tag = "RAW")
      
      # 把每個回合實際翻轉/位移的人數記下來
      meta <- data.table(
        p = p, n_total = n_total, contam_rate = contam_rate, sigma_fac = sigma_fac,
        flip_prop = flip_prop, kappa = kappa_tar, R = r,
        flipped = length(shared$contam_idx_flip),
        shifted = length(shared$contam_idx_shift)
      )
      
      ## --- EM 與 naive 都吃同一份 shared
      em_res    <- run_EM_from_shared(shared, p, beta_true, kappa_tar, cut_vec, em_diag_func = diag_all_1stepEM, verbose = FALSE)
      em_raw    <- run_EM_rawdiag_from_shared(shared, p, beta_true, kappa_tar, cut_vec, em_diag_func = diag_all_1stepEM, verbose = FALSE)
      naive_res <- run_naive_from_shared(shared, p, beta_true, kappa_tar, cut_vec, idx_use = em_res$keep)
      message(sprintf("p=%d, n=%d, R=%d, kappa=%.2f", p, n_total, r, kappa_tar))
      list(em = em_res, em_raw = em_raw, naive = naive_res, meta = meta)
    })
  }, future.seed = TRUE, future.scheduling = 1)
  
  flat <- unlist(setting_res, recursive = FALSE)
  flat <- Filter(function(z) !is.null(z$em) && !is.null(z$naive), flat)
  meta_part <- rbindlist(lapply(flat, function(z) z$meta), fill = TRUE)
  
  perf_part <- rbindlist(lapply(flat, function(z) rbind(z$em$perf, z$em_raw$perf, z$naive$perf)), fill = TRUE)
  diag_part <- rbindlist(lapply(flat, function(z) rbind(z$em$diag, z$em_raw$diag, z$naive$diag)), fill = TRUE)
  
  perf_part[, `:=`(p = p, n_total = n_total, contam_rate = contam_rate, sigma_fac = sigma_fac, flip_prop = flip_prop)]
  diag_part[, `:=`(p = p, n_total = n_total, contam_rate = contam_rate, sigma_fac = sigma_fac, flip_prop = flip_prop)]
  
  result_perf[[id]] <- perf_part
  result_diag[[id]] <- diag_part
  result_meta[[id]] <- meta_part
  id <- id + 1L
}

final_perf <- rbindlist(result_perf, fill = TRUE)
final_diag <- rbindlist(result_diag, fill = TRUE)
final_meta <- rbindlist(result_meta, fill = TRUE)

## 小結
## ---- 係數偏誤/MSE 的摘要（顯示 flip_prop 與翻轉/位移人數）----
summary_tbl <- final_perf %>%
  dplyr::group_by(method, estimator, kappa, p, n_total, contam_rate, sigma_fac, flip_prop) %>%
  dplyr::summarise(
    mean_bias = if (any(is.finite(bias))) mean(bias[is.finite(bias)]) else NA_real_,
    sd_bias   = if (sum(is.finite(bias)) > 1) sd(bias[is.finite(bias)]) else NA_real_,
    mean_mse  = if (any(is.finite(mse)))  mean(mse[is.finite(mse)])   else NA_real_,
    .groups   = "drop"
  ) %>%
  dplyr::mutate(
    contam_n = floor(n_total * contam_rate),
    flip_n   = as.integer(round(contam_n * flip_prop)),
    shift_n  = contam_n - flip_n
  ) %>%
  dplyr::relocate(flip_prop, flip_n, shift_n, contam_n, .after = sigma_fac)

tibble::as_tibble(summary_tbl) %>% print(n = Inf, width = Inf)

## ---- 診斷指標的摘要（同樣顯示 flip_prop 與翻轉/位移人數）----
summary_all <- final_diag %>%
  dplyr::group_by(estimator, method, kappa, p, n_total, contam_rate, sigma_fac, flip_prop) %>%
  dplyr::summarise(
    CIR_mean = mean(CIR, na.rm = TRUE),
    SR_mean  = mean(SR,  na.rm = TRUE),
    BA_mean  = mean((CIR + (1 - SR))/2, na.rm = TRUE),
    AUC_mean = mean(AUC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    contam_n = floor(n_total * contam_rate),
    flip_n   = as.integer(round(contam_n * flip_prop)),
    shift_n  = contam_n - flip_n
  ) %>%
  dplyr::relocate(flip_prop, flip_n, shift_n, contam_n, .after = sigma_fac)

# 先把每個設定下（p, n_total, contam_rate, sigma_fac, flip_prop）的平均翻轉/位移數算出來
meta_summary <- final_meta[, .(
  flipped_mean = mean(flipped),
  shifted_mean = mean(shifted)
), by = .(p, n_total, contam_rate, sigma_fac, flip_prop)]

# 合併到 summary_all（用相同的 key）
summary_all <- merge(
  as.data.frame(summary_all),    # dplyr tibble 轉 data.frame，避免型別小摩擦
  meta_summary,
  by = c("p","n_total","contam_rate","sigma_fac","flip_prop"),
  all.x = TRUE
)

tibble::as_tibble(summary_all) %>% print(n = Inf, width = Inf)

end_time <- Sys.time()
cat(sprintf("\n總運行時間：%.1f 秒\n", as.numeric(difftime(end_time, start_time, units = "secs"))))

out_path <- "C:/Users/kimi1/OneDrive/文件/論文/A run25.xlsx"
sheets <- list(
  Compare_Multi = as.data.frame(summary_all),
  Summary_Tbl   = as.data.frame(summary_tbl)
)
write_xlsx(x = sheets, path = out_path)