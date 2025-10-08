##################
#  diag_all_naive
##################

if (!requireNamespace("MASS", quietly = TRUE))
  stop("請先安裝套件 MASS 以使用 ginv()。")

# ───── 小工具 ───── #
.logistic <- function(z, clamp = 30) {
  z <- pmin(pmax(z, -clamp), clamp)
  1/(1 + exp(-z))
}
.safe_inv <- function(A) tryCatch(solve(A), error = function(e) MASS::ginv(A))
.pseudoR2_vec <- function(y, pi, type = c("CS","MF","NK"), eps = 1e-12) {
  type <- match.arg(type)
  y  <- as.numeric(y)
  if (!all(y %in% c(0,1))) stop("y 必須是 0/1")
  n  <- length(y)
  if (length(pi) != n) stop("y 與 pi 長度需相同")
  pi <- pmin(pmax(pi, eps), 1 - eps)
  mu <- pmin(pmax(mean(y), eps), 1 - eps)
  ll_full <- sum(y * log(pi) + (1 - y) * log(1 - pi))
  ll_null <- sum(y * log(mu) + (1 - y) * log(1 - mu))
  if (type == "CS") {
    return(1 - exp((2 / n) * (ll_null - ll_full)))
  } else if (type == "MF") {
    return(1 - (ll_full / ll_null))
  } else {
    cs <- 1 - exp((2 / n) * (ll_null - ll_full))
    denom <- 1 - exp((2 / n) * ll_null)
    if (!is.finite(denom) || denom <= 0) return(NA_real_)
    return(cs / denom)
  }
}

# ───── 未校正（naive）的一階近似診斷 ───── #
diag_all_naive <- function(Y, W, G.ind = NULL, cut_vec = NULL, dbg = FALSE) {
  stopifnot(is.matrix(W), length(Y) == nrow(W))
  n <- nrow(W); p <- ncol(W); k <- p + 1
  if (is.null(G.ind)) G.ind <- seq_len(n)
  B.ind <- setdiff(seq_len(n), G.ind)
  
  .fallback_cut <- function(cut_vec) {
    if (!is.null(cut_vec)) return(cut_vec)
    if (exists("cut_vec", envir = .GlobalEnv, inherits = FALSE))
      return(get("cut_vec", envir = .GlobalEnv))
    stop("cut_vec 未提供；請先用 make_cut_vec(...) 產生並傳入。")
  }
  cut_vec <- .fallback_cut(cut_vec)
  if (is.na(cut_vec["GCD_GSPR"])) cut_vec["GCD_GSPR"] <- 1
  
  Xg <- cbind(1, W[G.ind, , drop = FALSE])
  df_g <- data.frame(Y = Y[G.ind], Xg[, -1, drop = FALSE])
  colnames(df_g) <- c("Y", colnames(W))
  
  fit_g <- try(glm(Y ~ ., family = binomial(), data = df_g,
                   control = glm.control(maxit = 200, epsilon = 1e-8)), silent = TRUE)
  if (inherits(fit_g, "try-error") || any(!is.finite(coef(fit_g)))) {
    if (requireNamespace("brglm2", quietly = TRUE)) {
      fit_g <- try(brglm2::brglm(Y ~ ., data = df_g,
                                 family = binomial("logit"), type = "AS_mixed"),
                   silent = TRUE)
    }
  }
  if (inherits(fit_g, "try-error") || any(!is.finite(coef(fit_g))))
    stop("good-only GLM 失敗：請檢查資料或門檻設定。")
  
  beta_ref <- coef(fit_g)
  pig  <- .logistic(as.vector(Xg %*% beta_ref), clamp = 30)
  vg   <- pig * (1 - pig)
  J_clean <- crossprod(Xg, Xg * vg)
  J_clean_inv <- .safe_inv(J_clean)
  
  X     <- cbind(1, W)
  pi_ref<- .logistic(as.vector(X %*% beta_ref), clamp = 30)
  v_all <- pmax(pi_ref * (1 - pi_ref), 1e-8)
  
  quad_xJx <- rowSums((X %*% J_clean_inv) * X)
  hii_G_all <- pmin(pmax(v_all * quad_xJx, 1e-8), 1 - 1e-6)
  
  Zhat <- X * sqrt(v_all)
  M_hat_inv <- .safe_inv(crossprod(Zhat))
  quad_xhat <- rowSums((X %*% M_hat_inv) * X)
  hii_M_all <- pmin(pmax(v_all * quad_xhat, 1e-8), 1 - 1e-6)
  
  R_CS_all <- .pseudoR2_vec(Y[G.ind], pi_ref[G.ind], type = "CS")
  R_MF_all <- .pseudoR2_vec(Y[G.ind], pi_ref[G.ind], type = "MF")
  R_NK_all <- .pseudoR2_vec(Y[G.ind], pi_ref[G.ind], type = "NK")
  
  term_vec <- quad_xJx
  
  out <- data.frame(id = seq_len(n),
                    GD = NA, MD = NA, GDF = NA, MDF = NA,
                    GDB = NA, MDB = NA, GCD_GSPR = NA,
                    mCDstar = NA, R_CS = NA, R_MF = NA, R_NK = NA)
  
  Rdiff_CS <- Rdiff_MF <- Rdiff_NK <- numeric(n)
  
  df_all <- data.frame(Y, W)
  glm_fit <- try(glm(Y ~ ., data = df_all, family = binomial(),
                     control = glm.control(maxit = 200, epsilon = 1e-8)), silent = TRUE)
  if (inherits(glm_fit, "try-error") || any(!is.finite(coef(glm_fit)))) {
    if (requireNamespace("brglm2", quietly = TRUE)) {
      glm_fit <- try(brglm2::brglm(Y ~ ., data = df_all,
                                   family = binomial("logit"), type = "AS_mixed"),
                     silent = TRUE)
    }
  }
  CookD_trad  <- DFFITS_trad <- rep(NA_real_, n)
  if (!inherits(glm_fit, "try-error") && is.list(glm_fit)) {
    CookD_trad  <- suppressWarnings(cooks.distance(glm_fit))
    DFFITS_trad <- suppressWarnings(dffits(glm_fit))
  }
  
  out <- cbind(out, CookD = CookD_trad, DFFITS = DFFITS_trad)
  
  resid_vec <- (Y - pi_ref)
  scores_all <- X * resid_vec
  Zhat        <- X * sqrt(v_all)
  XtVtX_hat   <- crossprod(Zhat)   # modified 幾何：X' V X
  
  for (i in seq_len(n)) {
    vi    <- v_all[i]
    hii_G <- hii_G_all[i]
    hii_M <- hii_M_all[i]
    
    score_i <- scores_all[i, ]
    den   <- if (i %in% B.ind) (1 + hii_G) else (1 - hii_G)
    d_beta <- as.vector(J_clean_inv %*% score_i) / den
    if (!(i %in% B.ind)) d_beta <- -d_beta
    
    beta_new_i <- as.vector(beta_ref + d_beta)
    idx_R <- if (i %in% B.ind) c(G.ind, i) else setdiff(G.ind, i)
    
    pi_approx_i <- .logistic(as.vector(cbind(1, W[idx_R, , drop = FALSE]) %*% beta_new_i),
                             clamp = 30)
    
    # 先建 V 權重的 Modified 資訊
    eta_ref  <- as.vector(X %*% beta_ref)
    eta_ref  <- pmin(pmax(eta_ref, -30), 30)         # 同 EM 版，避免數值爆炸
    pi_ref   <- plogis(eta_ref)
    Vref_s   <- sqrt(pi_ref * (1 - pi_ref))
    
    # 之後
    GD_i <- as.numeric(t(d_beta) %*% J_clean   %*% d_beta / k)  # k = p+1
    MD_i <- as.numeric(t(d_beta) %*% XtVtX_hat %*% d_beta / k)
    
    R_CS_i <- .pseudoR2_vec(Y[idx_R], pi_approx_i, type = "CS")
    R_MF_i <- .pseudoR2_vec(Y[idx_R], pi_approx_i, type = "MF")
    R_NK_i <- .pseudoR2_vec(Y[idx_R], pi_approx_i, type = "NK")
    Rdiff_CS[i] <- abs(R_CS_all - R_CS_i)
    Rdiff_MF[i] <- abs(R_MF_all - R_MF_i)
    Rdiff_NK[i] <- abs(R_NK_all - R_NK_i)
    
    if (i %in% B.ind) {
      ri_G <- (Y[i] - pi_ref[i]) / sqrt(vi * (1 + hii_G));  hiiS_G <- hii_G / (1 + hii_G)
      ri_M <- (Y[i] - pi_ref[i]) / sqrt(vi * (1 + hii_M));  hiiS_M <- hii_M / (1 + hii_M)
    } else {
      ri_G <- (Y[i] - pi_ref[i]) / sqrt(vi * (1 - hii_G));  hiiS_G <- hii_G / (1 - hii_G)
      ri_M <- (Y[i] - pi_ref[i]) / sqrt(vi * (1 - hii_M));  hiiS_M <- hii_M / (1 - hii_M)
    }
    
    GDF_i <- ri_G * sqrt(hiiS_G)
    MDF_i <- ri_M * sqrt(hiiS_M)
    mCD_i <- abs(MDF_i) * sqrt((n - length(B.ind) - p) / p)
    GDB_i <- hiiS_G * ri_G^2
    MDB_i <- hiiS_M * ri_M^2
    GCD_i <- (1/(p+1)) * (ri_G^2) * term_vec[i]
    
    out[i, c("GD","MD","GDF","MDF","GDB","MDB","GCD_GSPR","mCDstar","R_CS","R_MF","R_NK")] <-
      c(GD_i, MD_i, GDF_i, MDF_i, GDB_i, MDB_i, GCD_i, mCD_i, Rdiff_CS[i], Rdiff_MF[i], Rdiff_NK[i])
  }
  
  # --- R-squared based diagnostics: raw & standardized (match EM) ---
  
  sd_safe  <- function(x) { x <- x[is.finite(x)]; if (length(x) >= 2) sd(x) else NA_real_ }
  mad_safe <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2) return(NA_real_)
    stats::mad(x, center = stats::median(x), constant = 1.4826, na.rm = TRUE)
  }
  
  # 共同尺度（與 EM 對齊）：優先用「好點」的 MAD；失敗→ overall SD；再失敗→ 1e-6
  .get_common_scale <- function(name, Rdiff, G.ind, B.ind, verbose = FALSE) {
    s <- mad_safe(Rdiff[G.ind])
    if (!is.finite(s) || s <= 0) s <- sd_safe(Rdiff)
    if (!is.finite(s) || s <= 0) s <- 1e-6
    if (verbose) {
      sdb <- sd_safe(Rdiff[B.ind]); sdg <- sd_safe(Rdiff[G.ind])
      message(sprintf("[naive][scale %s] sd_bad=%.4g, sd_good=%.4g, use_scale=%.6g",
                      name, sdb, sdg, s))
    }
    s
  }
  
  # Raw differences
  out$R_CS_raw <- abs(Rdiff_CS)
  out$R_MF_raw <- abs(Rdiff_MF)
  out$R_NK_raw <- abs(Rdiff_NK)
  
  # Standardized values（與 EM 同一步驟）
  s_cs <- .get_common_scale("CS", Rdiff_CS, G.ind, B.ind, verbose = FALSE)
  s_mf <- .get_common_scale("MF", Rdiff_MF, G.ind, B.ind, verbose = FALSE)
  s_nk <- .get_common_scale("NK", Rdiff_NK, G.ind, B.ind, verbose = FALSE)
  
  out$R_CS <- out$R_CS_raw / s_cs
  out$R_MF <- out$R_MF_raw / s_mf
  out$R_NK <- out$R_NK_raw / s_nk
  
  # 旗標
  GDFFITS_cut  <- cut_vec["GDFFITS"]; GSDFBETA_cut <- cut_vec["GSDFBETA"]
  out <- transform(out,
                   flag_CookD    = CookD   > cut_vec["CookD"],
                   flag_DFFITS   = abs(DFFITS) > cut_vec["DFFITS"],
                   flag_GD       = GD  > cut_vec["GD"],
                   flag_MD       = MD  > cut_vec["MD"],
                   flag_GDF      = abs(GDF) > GDFFITS_cut,
                   flag_MDF      = abs(MDF) > GDFFITS_cut,
                   flag_GDB      = abs(GDB) > GSDFBETA_cut,
                   flag_MDB      = abs(MDB) > GSDFBETA_cut,
                   flag_mCDstar  = mCDstar > cut_vec["mCDstar"],
                   flag_GCD_GSPR = abs(GCD_GSPR) > cut_vec["GCD_GSPR"],
                   flag_R_CS     = R_CS > cut_vec["StdR2_CS"],
                   flag_R_MF     = R_MF > cut_vec["StdR2_MF"],
                   flag_R_NK     = R_NK > cut_vec["StdR2_NK"])
  return(out)
}