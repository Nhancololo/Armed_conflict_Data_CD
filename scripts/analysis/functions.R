# BAMLSS spatial theme
theme_bamlss_spatial <- function(base_size = 12, base_family = "serif") {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      text = ggplot2::element_text(family = base_family),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(color = "gray40", size = base_size),
      strip.background = ggplot2::element_rect(fill = "gray95"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

# Global color palette
spatial_palette <- function() {
  c(
    primary = "#fc4e2a",
    secondary = "#fd8d3c",
    highlight = "#b10026",
    neutral = "#feb24c",
    light = "#f5f5f5"
  )
}

# Bounded logarithm
bounded_log <- function(x, minval = 1e-300) {
  log(pmax(x, minval))
}

# Build spatial penalty matrix
build_spatial_penalty <- function(shape, id_col) {
  if (!inherits(shape, "sf")) {
    stop("Shape must be an sf object.", call. = FALSE)
  }
  
  shape <- shape[order(shape[[id_col]]), ]
  nb <- spdep::poly2nb(shape, row.names = shape[[id_col]])
  adj_mat <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
  penalty_mat <- diag(rowSums(adj_mat)) - adj_mat
  rownames(penalty_mat) <- colnames(penalty_mat) <- rownames(adj_mat)
  
  list(
    nb = nb,
    penalty = penalty_mat,
    listw = spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  )
}

# make.link("log")$linkinv floors at double.eps, so the likelihood is flat
# in eta below log(double.eps); ZANBI and ZOINB reach that region
positive_log_linkinv <- function(eta) {
  pmax(exp(eta), .Machine$double.xmin)
}

# Numerical score and hessian factory on predictor scale
num_score_hess_factory <- function(d_fun, link, positive_log = FALSE) {
  lk <- tryCatch(stats::make.link(link), error = function(e) link_object(link))
  
  # Only values sitting exactly on the boundary are moved inside the open
  # parameter space; interior values are differentiated where they are.
  perturb <- function(par, id, delta) {
    z <- par
    x <- z[[id]]
    if (link == "logit") {
      x <- pmin(pmax(x, .Machine$double.xmin), 1 - .Machine$double.neg.eps)
    }
    if (link == "log") x <- pmax(x, .Machine$double.xmin)
    eta <- lk$linkfun(x)
    z[[id]] <- if (identical(link, "log") && isTRUE(positive_log)) {
      positive_log_linkinv(eta + delta)
    } else {
      lk$linkinv(eta + delta)
    }
    z
  }
  
  logdens <- function(y, par) {
    out <- tryCatch(
      d_fun(y, par, log = TRUE),
      error = function(e) rep(-Inf, length(y))
    )
    out <- as.numeric(out)
    out[!is.finite(out)] <- log(1e-300)
    out
  }
  
  score_fun <- function(y, par, id = NULL, ..., delta = 1e-5) {
    lp <- logdens(y, perturb(par, id,  delta))
    lm <- logdens(y, perturb(par, id, -delta))
    (lp - lm) / (2 * delta)
  }
  
  hess_fun <- function(y, par, id = NULL, ..., delta = 1e-4) {
    lp <- logdens(y, perturb(par, id,  delta))
    l0 <- logdens(y, par)
    lm <- logdens(y, perturb(par, id, -delta))
    h <- -(lp - 2 * l0 + lm) / (delta^2)
    # Where the observed information is not positive, the squared score
    # (outer product of gradients) is the working weight
    bad <- !is.finite(h) | h <= 1e-10
    if (any(bad)) {
      g <- (lp - lm) / (2 * delta)
      h[bad] <- g[bad]^2
    }
    h[!is.finite(h) | h <= 1e-10] <- 1e-8
    h
  }
  
  attr(score_fun, "numerical_derivative") <- TRUE
  attr(hess_fun, "numerical_derivative") <- TRUE
  list(score = score_fun, hess = hess_fun)
}

# Zero-one inflated negative binomial (ZOINB) family
family_zoinb <- function() {
  d_fun <- function(y, par, log = FALSE) {
    n <- max(length(y), length(par$mu), length(par$sigma),
             length(par$nu), length(par$tau))
    y <- rep_len(as.numeric(y), n)
    mu <- pmax(rep_len(as.numeric(par$mu), n), .Machine$double.xmin)
    sigma <- pmax(rep_len(as.numeric(par$sigma), n), .Machine$double.xmin)
    nu <- pmin(pmax(rep_len(as.numeric(par$nu), n),
                    .Machine$double.eps), 1 - .Machine$double.eps)
    tau <- pmin(pmax(rep_len(as.numeric(par$tau), n),
                     .Machine$double.eps), 1 - .Machine$double.eps)
    size <- 1 / sigma
    
    log_denom <- stats::pnbinom(
      1, mu = mu, size = size, lower.tail = FALSE, log.p = TRUE
    )
    
    ld <- rep(-Inf, n)
    valid_y <- is.finite(y) & y >= 0 & abs(y - round(y)) < 1e-8
    idx0 <- valid_y & y == 0
    idx1 <- valid_y & y == 1
    idx_gt1 <- valid_y & y > 1
    
    if (any(idx0)) ld[idx0] <- log(nu[idx0])
    if (any(idx1)) ld[idx1] <- log1p(-nu[idx1]) + log(tau[idx1])
    if (any(idx_gt1)) {
      ld[idx_gt1] <- log1p(-nu[idx_gt1]) + log1p(-tau[idx_gt1]) +
        stats::dnbinom(
          y[idx_gt1], mu = mu[idx_gt1], size = size[idx_gt1], log = TRUE
        ) -
        log_denom[idx_gt1]
    }
    if (log) ld else exp(ld)
  }
  
  # CDF from log-survival ratios, with the clipping of d_fun
  p_fun <- function(y, par, ...) {
    y <- as.numeric(y)
    n <- length(y)
    mu <- pmax(rep_len(as.numeric(par$mu), n), .Machine$double.xmin)
    sigma <- pmax(rep_len(as.numeric(par$sigma), n), .Machine$double.xmin)
    nu <- pmin(pmax(rep_len(as.numeric(par$nu), n),
                    .Machine$double.eps), 1 - .Machine$double.eps)
    tau <- pmin(pmax(rep_len(as.numeric(par$tau), n),
                     .Machine$double.eps), 1 - .Machine$double.eps)
    size <- 1 / sigma
    out <- rep(NA_real_, n)
    valid <- is.finite(y) & is.finite(mu) & is.finite(sigma) &
      is.finite(nu) & is.finite(tau) & mu > 0 & sigma > 0
    q <- floor(y)
    out[valid & q < 0] <- 0
    out[valid & q == 0] <- nu[valid & q == 0]
    idx1 <- valid & q == 1
    out[idx1] <- nu[idx1] + (1 - nu[idx1]) * tau[idx1]
    idx <- which(valid & q >= 2)
    if (length(idx)) {
      log_sf1 <- stats::pnbinom(1, mu = mu[idx], size = size[idx],
                                lower.tail = FALSE, log.p = TRUE)
      log_sfy <- stats::pnbinom(q[idx], mu = mu[idx], size = size[idx],
                                lower.tail = FALSE, log.p = TRUE)
      good <- is.finite(log_sf1) & is.finite(log_sfy)
      if (any(good)) {
        ii <- idx[good]
        log_ratio <- pmin(log_sfy[good] - log_sf1[good], 0)
        cdf_trunc <- -expm1(log_ratio)
        out[ii] <- nu[ii] + (1 - nu[ii]) * tau[ii] +
          (1 - nu[ii]) * (1 - tau[ii]) * cdf_trunc
      }
    }
    pmin(pmax(out, 0), 1)
  }
  
  sh_mu <- num_score_hess_factory(d_fun, "log", positive_log = TRUE)
  sh_sigma <- num_score_hess_factory(d_fun, "log", positive_log = TRUE)
  sh_nu <- num_score_hess_factory(d_fun, "logit")
  sh_tau <- num_score_hess_factory(d_fun, "logit")
  
  f <- list(
    family = "ZOINB",
    names = c("mu", "sigma", "nu", "tau"),
    links = c("log", "log", "logit", "logit"),
    map2par = function(eta) {
      out <- vector("list", length(eta))
      names(out) <- names(eta)
      link_vec <- stats::setNames(
        c("log", "log", "logit", "logit"),
        c("mu", "sigma", "nu", "tau")
      )
      for (nm in names(out)) {
        out[[nm]] <- if (identical(link_vec[[nm]], "log")) {
          positive_log_linkinv(eta[[nm]])
        } else {
          stats::make.link(link_vec[[nm]])$linkinv(eta[[nm]])
        }
      }
      out
    },
    d = d_fun,
    p = p_fun,
    loglik = function(y, par, ...) sum(d_fun(y, par, log = TRUE)),
    score = list(
      mu = function(y, p, ...) sh_mu$score(y, p, id = "mu"),
      sigma = function(y, p, ...) sh_sigma$score(y, p, id = "sigma"),
      nu = function(y, p, ...) sh_nu$score(y, p, id = "nu"),
      tau = function(y, p, ...) sh_tau$score(y, p, id = "tau")
    ),
    hess = list(
      mu = function(y, p, ...) sh_mu$hess(y, p, id = "mu"),
      sigma = function(y, p, ...) sh_sigma$hess(y, p, id = "sigma"),
      nu = function(y, p, ...) sh_nu$hess(y, p, id = "nu"),
      tau = function(y, p, ...) sh_tau$hess(y, p, id = "tau")
    )
  )
  class(f) <- "family.bamlss"
  f
}

# Zero-adjusted negative binomial, gamlss.dist::ZANBI parameterisation:
# NBI(mu, sigma) for the positive part and nu = P(Y = 0). log P(Y > 0) is
# computed directly; bamlss::ZANBI_bamlss gives +Inf for y > 0 when
# mu * sigma is small.
family_zanbi <- function() {
  links <- c(mu = "log", sigma = "log", nu = "logit")
  
  # log P(Y = 0) = -log(1 + mu * sigma) / sigma, which tends to -mu as
  # sigma -> 0; the log1p term is split when mu * sigma overflows.
  log_p0 <- function(mu, sigma) {
    ms <- mu * sigma
    l1p <- ifelse(is.finite(ms), log1p(ms), log(mu) + log(sigma))
    ifelse(ms > 0, -l1p / sigma, -mu)
  }
  
  log_pos <- function(mu, sigma) log(-expm1(log_p0(mu, sigma)))
  
  log_nb <- function(y, mu, sigma) {
    out <- stats::dpois(y, mu, log = TRUE)
    nb <- sigma > 0
    out[nb] <- stats::dnbinom(
      y[nb], size = 1 / sigma[nb], mu = mu[nb], log = TRUE
    )
    out
  }
  
  # log P(Y > q) for the untruncated NBI, q >= 0
  log_upper <- function(q, mu, sigma) {
    out <- stats::ppois(q, mu, lower.tail = FALSE, log.p = TRUE)
    nb <- sigma > 0
    out[nb] <- stats::pnbinom(
      q[nb], size = 1 / sigma[nb], mu = mu[nb],
      lower.tail = FALSE, log.p = TRUE
    )
    out
  }
  
  expand <- function(y, par) {
    n <- max(length(y), length(par$mu), length(par$sigma), length(par$nu))
    z <- list(
      y = rep_len(as.numeric(y), n),
      mu = rep_len(as.numeric(par$mu), n),
      sigma = rep_len(as.numeric(par$sigma), n),
      nu = rep_len(as.numeric(par$nu), n)
    )
    z$ok <- with(z, is.finite(mu) & mu > 0 & is.finite(sigma) &
                   sigma >= 0 & is.finite(nu) & nu >= 0 & nu <= 1)
    z
  }
  
  d_fun <- function(y, par, log = FALSE) {
    z <- expand(y, par)
    count <- is.finite(z$y) & z$y >= 0 & abs(z$y - round(z$y)) < 1e-8
    ld <- rep(-Inf, length(z$y))
    
    zero <- z$ok & count & z$y == 0
    ld[zero] <- log(z$nu[zero])
    
    pos <- z$ok & count & z$y > 0
    if (any(pos)) {
      mu <- z$mu[pos]
      sigma <- z$sigma[pos]
      ld[pos] <- log1p(-z$nu[pos]) + log_nb(z$y[pos], mu, sigma) -
        log_pos(mu, sigma)
    }
    if (log) ld else exp(ld)
  }
  
  p_fun <- function(y, par, ...) {
    z <- expand(y, par)
    q <- floor(z$y)
    ok <- z$ok & !is.na(q)
    out <- rep(NA_real_, length(q))
    out[ok & q < 0] <- 0
    out[ok & q == 0] <- z$nu[ok & q == 0]
    
    pos <- ok & q >= 1
    if (any(pos)) {
      mu <- z$mu[pos]
      sigma <- z$sigma[pos]
      tail_ratio <- exp(log_upper(q[pos], mu, sigma) - log_pos(mu, sigma))
      out[pos] <- z$nu[pos] + (1 - z$nu[pos]) * (1 - pmin(tail_ratio, 1))
    }
    out
  }
  
  sh_mu <- num_score_hess_factory(d_fun, "log", positive_log = TRUE)
  sh_sigma <- num_score_hess_factory(d_fun, "log", positive_log = TRUE)
  sh_nu <- num_score_hess_factory(d_fun, "logit")
  
  f <- list(
    family = "ZANBI",
    names = names(links),
    links = links,
    type = "discrete",
    map2par = function(eta) {
      out <- lapply(names(eta), function(nm) {
        if (identical(links[[nm]], "log")) {
          positive_log_linkinv(eta[[nm]])
        } else {
          stats::make.link(links[[nm]])$linkinv(eta[[nm]])
        }
      })
      stats::setNames(out, names(eta))
    },
    d = d_fun,
    p = p_fun,
    loglik = function(y, par, ...) sum(d_fun(y, par, log = TRUE)),
    initialize = list(
      mu = function(y, ...) {
        rep(if (any(y > 0)) mean(y[y > 0]) else 1, length(y))
      },
      sigma = function(y, ...) rep(1, length(y)),
      nu = function(y, ...) rep(mean(y == 0), length(y))
    ),
    score = list(
      mu = function(y, p, ...) sh_mu$score(y, p, id = "mu"),
      sigma = function(y, p, ...) sh_sigma$score(y, p, id = "sigma"),
      nu = function(y, p, ...) sh_nu$score(y, p, id = "nu")
    ),
    hess = list(
      mu = function(y, p, ...) sh_mu$hess(y, p, id = "mu"),
      sigma = function(y, p, ...) sh_sigma$hess(y, p, id = "sigma"),
      nu = function(y, p, ...) sh_nu$hess(y, p, id = "nu")
    )
  )
  class(f) <- "family.bamlss"
  f
}

# log K_v(x), modified Bessel function of the second kind; the uniform
# (Debye) expansion takes over where besselK() overflows (large order v).
log_bessel_k <- function(x, v) {
  n <- max(length(x), length(v))
  x <- rep_len(as.numeric(x), n)
  v <- abs(rep_len(as.numeric(v), n))
  out <- suppressWarnings(log(besselK(x, v, expon.scaled = TRUE)) - x)
  
  bad <- !is.finite(out) & is.finite(x) & is.finite(v) & x > 0 & v > 0
  if (any(bad)) {
    xb <- x[bad]
    vb <- v[bad]
    s <- sqrt(1 + (xb / vb)^2)
    t <- 1 / s
    u1 <- (3 * t - 5 * t^3) / 24
    u2 <- (81 * t^2 - 462 * t^4 + 385 * t^6) / 1152
    u3 <- (30375 * t^3 - 369603 * t^5 + 765765 * t^7 - 425425 * t^9) / 414720
    series <- 1 - u1 / vb + u2 / vb^2 - u3 / vb^3
    out[bad] <- 0.5 * log(pi / (2 * vb)) - vb * (s + log(xb / vb / (1 + s))) -
      0.5 * log(s) + log(series)
  }
  out
}

# Sichel family, gamlss.dist::SICHEL parameterisation (mean mu, dispersion
# sigma, shape nu), with c = K_{nu+1}(1/sigma) / K_nu(1/sigma) and
# alpha = sqrt(1 + 2 * sigma * mu / c) / sigma:
#   P(Y = y) = (mu / c)^y K_{nu+y}(alpha) / (y! (sigma * alpha)^(y + nu) *
#              K_nu(1 / sigma)).
# Evaluated in R on the log scale, since the compiled dSICHEL can crash the
# session inside the sampler. For sigma > 1e4 and nu > 0 the negative
# binomial limit is kept.
family_sichel <- function() {
  links <- c(mu = "log", sigma = "log", nu = "identity")
  
  expand <- function(y, par) {
    n <- max(length(y), length(par$mu), length(par$sigma), length(par$nu))
    z <- list(
      y = rep_len(as.numeric(y), n),
      mu = rep_len(as.numeric(par$mu), n),
      sigma = rep_len(as.numeric(par$sigma), n),
      nu = rep_len(as.numeric(par$nu), n)
    )
    z$ok <- with(z, is.finite(mu) & mu > 0 & is.finite(sigma) & sigma > 0 &
                   is.finite(nu))
    z
  }
  
  log_pmf <- function(y, mu, sigma, nu) {
    log_k0 <- log_bessel_k(1 / sigma, nu)
    cc <- exp(log_bessel_k(1 / sigma, nu + 1) - log_k0)
    alpha <- sqrt(1 + 2 * sigma * mu / cc) / sigma
    -lgamma(y + 1) - nu * log(sigma * alpha) +
      y * log(mu / (cc * sigma * alpha)) +
      log_bessel_k(alpha, nu + y) - log_k0
  }
  
  d_fun <- function(y, par, log = FALSE) {
    z <- expand(y, par)
    count <- is.finite(z$y) & z$y >= 0 & abs(z$y - round(z$y)) < 1e-8
    ld <- rep(-Inf, length(z$y))
    
    nb <- z$ok & count & z$sigma > 1e4 & z$nu > 0
    if (any(nb)) {
      ld[nb] <- stats::dnbinom(z$y[nb], size = z$nu[nb], mu = z$mu[nb],
                               log = TRUE)
    }
    si <- z$ok & count & !nb
    if (any(si)) {
      ld[si] <- suppressWarnings(
        log_pmf(round(z$y[si]), z$mu[si], z$sigma[si], z$nu[si])
      )
    }
    ld[is.na(ld)] <- -Inf
    if (log) ld else exp(ld)
  }
  
  p_fun <- function(y, par, ...) {
    z <- expand(y, par)
    q <- floor(z$y)
    out <- rep(NA_real_, length(q))
    out[z$ok & !is.na(q) & q < 0] <- 0
    for (i in which(z$ok & !is.na(q) & q >= 0)) {
      k <- seq.int(0, q[i])
      pk <- d_fun(k, list(mu = z$mu[i], sigma = z$sigma[i], nu = z$nu[i]))
      out[i] <- min(sum(pk), 1)
    }
    out
  }
  
  sh_mu <- num_score_hess_factory(d_fun, "log")
  sh_sigma <- num_score_hess_factory(d_fun, "log")
  sh_nu <- num_score_hess_factory(d_fun, "identity")
  
  f <- list(
    family = "SICHEL",
    names = names(links),
    links = links,
    type = "discrete",
    map2par = function(eta) {
      out <- lapply(names(eta), function(nm) {
        stats::make.link(links[[nm]])$linkinv(eta[[nm]])
      })
      stats::setNames(out, names(eta))
    },
    d = d_fun,
    p = p_fun,
    loglik = function(y, par, ...) sum(d_fun(y, par, log = TRUE)),
    initialize = list(
      mu = function(y, ...) (y + mean(y)) / 2,
      sigma = function(y, ...) {
        m <- mean(y)
        rep(max((stats::var(y) - m) / m^2, 0.1), length(y))
      },
      nu = function(y, ...) rep(-0.5, length(y))
    ),
    score = list(
      mu = function(y, p, ...) sh_mu$score(y, p, id = "mu"),
      sigma = function(y, p, ...) sh_sigma$score(y, p, id = "sigma"),
      nu = function(y, p, ...) sh_nu$score(y, p, id = "nu")
    ),
    hess = list(
      mu = function(y, p, ...) sh_mu$hess(y, p, id = "mu"),
      sigma = function(y, p, ...) sh_sigma$hess(y, p, id = "sigma"),
      nu = function(y, p, ...) sh_nu$hess(y, p, id = "nu")
    )
  )
  class(f) <- "family.bamlss"
  f
}


# Log-likelihood as the sum of log-densities. The default bamlss version
# replaces non-finite values by -100; here such states get zero density.
strict_loglik <- function(d_fun) {
  force(d_fun)
  function(y, par, ...) {
    ld <- suppressWarnings(as.numeric(d_fun(y, par, log = TRUE)))
    if (length(ld) != length(y) || anyNA(ld) || any(ld == Inf)) return(-Inf)
    sum(ld)
  }
}

# Observation-safe log-density for gamlss.dist families. Several of their
# densities recycle parameters wrongly when a subset (x == 0, sigma < 1e-4,
# sigma > 1e4) is smaller than n. Each subset pattern is evaluated on its
# own; what cannot be evaluated gets log-density -Inf.
safe_density <- function(d_fun) {
  force(d_fun)
  
  eval_group <- function(y, par) {
    out <- tryCatch(
      suppressWarnings(as.numeric(d_fun(y, par, log = TRUE))),
      error = function(e) NULL
    )
    if (length(out) == length(y)) return(out)
    vapply(seq_along(y), function(i) {
      z <- tryCatch(
        suppressWarnings(
          as.numeric(d_fun(y[i], lapply(par, `[`, i), log = TRUE))
        ),
        error = function(e) NA_real_
      )
      if (length(z) == 1L) z else NA_real_
    }, numeric(1))
  }
  
  function(y, par, log = FALSE, ...) {
    y <- as.numeric(unlist(y, use.names = FALSE))
    n <- length(y)
    par <- lapply(par, function(z) rep_len(as.numeric(z), n))
    unusable <- !is.finite(y)
    key <- as.integer(y == 0)
    for (z in par) {
      unusable <- unusable | !is.finite(z)
      key <- key * 4L + (z <= 1e-4) * 2L + (z >= 1e4)
    }
    
    out <- rep(-Inf, n)
    ok <- which(!unusable)
    for (idx in split(ok, key[ok])) {
      ld <- eval_group(y[idx], lapply(par, `[`, idx))
      ld[is.na(ld) | ld == Inf] <- -Inf
      out[idx] <- ld
    }
    if (log) out else exp(out)
  }
}

# Score and hessian that cannot stop a chain: missing, infinite or failed
# derivatives fall back to a central difference of the safe density. They
# only shape the proposal, not the target.
guard_derivative <- function(fun, backup, p) {
  force(fun)
  force(backup)
  function(y, par, ...) {
    n <- NROW(y)
    out <- tryCatch(
      suppressWarnings(as.numeric(fun(y, par, ...))),
      error = function(e) NULL
    )
    bad <- if (length(out) == n) !is.finite(out) else rep(TRUE, n)
    if (!any(bad)) return(out)
    alt <- backup(y, par, id = p)
    if (length(out) != n) return(alt)
    out[bad] <- alt[bad]
    out
  }
}

guard_derivatives <- function(fam) {
  for (p in fam$names) {
    use_positive_log <-
      fam$family %in% c("ZANBI", "ZOINB") &&
      identical(unname(fam$links[[p]]), "log")
    num <- num_score_hess_factory(
      fam$d, fam$links[[p]], positive_log = use_positive_log
    )
    for (what in c("score", "hess")) {
      fam[[what]][[p]] <- guard_derivative(fam[[what]][[p]], num[[what]], p)
    }
  }
  fam
}

# Small design with parameters that change from one observation to the next
# and two extreme second-parameter values that reach the 1e-4 / 1e4 branches.
family_probe_design <- function(fam) {
  y <- c(0, 0, 0, 1, 2, 0, 3, 5, 0, 8, 1, 0, 2, 13, 0, 4, 1, 0, 6, 2)
  i <- seq_along(y)
  eta <- lapply(seq_along(fam$names), function(j) {
    1.2 * sin(i * (1.3 + 0.7 * j))
  })
  names(eta) <- fam$names
  if (length(eta) > 1L) eta[[2L]][c(5L, 13L)] <- c(-11, 11)
  list(y = y, par = fam$map2par(eta))
}

# Compare the density on the whole design with the density observation by
# observation; a difference means it depends on batching.
probe_family_density <- function(fam, raw_d, safe_d) {
  des <- family_probe_design(fam)
  n <- length(des$y)
  eval_one <- function(f, i) {
    z <- tryCatch(
      suppressWarnings(
        as.numeric(f(des$y[i], lapply(des$par, `[`, i), log = TRUE))
      ),
      error = function(e) NA_real_
    )
    if (length(z) == 1L) z else NA_real_
  }
  eval_all <- function(f) {
    z <- tryCatch(
      suppressWarnings(as.numeric(f(des$y, des$par, log = TRUE))),
      error = function(e) NULL
    )
    if (length(z) == n) z else rep(NA_real_, n)
  }
  ref <- vapply(seq_len(n), function(i) eval_one(raw_d, i), numeric(1))
  use <- is.finite(ref)
  agrees <- function(z) {
    all(is.finite(z[use]) &
          abs(z[use] - ref[use]) <= 1e-8 * (1 + abs(ref[use])))
  }
  list(
    batch_sensitive = !agrees(eval_all(raw_d)),
    safe_consistent = agrees(eval_all(safe_d)),
    usable = mean(use)
  )
}

family_probe_cache <- new.env(parent = emptyenv())

# Gauss-Legendre quadrature
gauss_legendre <- function(n = 80L) {
  n <- as.integer(n)
  if (n < 8L) stop("`n_quad` must be at least 8.", call. = FALSE)
  
  i <- seq_len(n - 1L)
  off <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- off
  J[cbind(i + 1L, i)] <- off
  
  ee <- eigen(J, symmetric = TRUE)
  ord <- order(ee$values)
  list(nodes = ee$values[ord], weights = 2 * (ee$vectors[1L, ord]^2))
}

# RSB (Rescaled-Beta mixture) family
family_rsb <- function(a = 0.5, b = 0.5, n_quad = 30L) {
  if (!isTRUE(all.equal(as.numeric(a), 0.5)) ||
      !isTRUE(all.equal(as.numeric(b), 0.5))) {
    stop("The verified RSB implementation requires a=b=0.5.", call. = FALSE)
  }
  
  gl <- gauss_legendre(n_quad)
  tt <- (pi / 4) * (gl$nodes + 1)
  L <- tan(tt)^2
  
  log_eta <- numeric(length(L))
  small_L <- L < 40
  log_eta[small_L] <- log(expm1(L[small_L]))
  log_eta[!small_L] <- L[!small_L] + log1p(-exp(-L[!small_L]))
  
  quad_w <- gl$weights / 2
  quad_w <- quad_w / sum(quad_w)
  log_quad_w <- log(quad_w)
  
  log_pois_from_lograte <- function(y, log_rate) {
    y <- as.numeric(y)
    log_rate <- as.numeric(log_rate)
    out <- rep(-Inf, max(length(y), length(log_rate)))
    y <- rep(y, length.out = length(out))
    log_rate <- rep(log_rate, length.out = length(out))
    
    good <-
      is.finite(y) & y >= 0 & abs(y - round(y)) < 1e-8 & is.finite(log_rate)
    moderate <- good & log_rate < 700
    
    if (any(moderate)) {
      out[moderate] <- y[moderate] * log_rate[moderate] -
        exp(log_rate[moderate]) - lgamma(y[moderate] + 1)
    }
    out
  }
  
  log_integrated_pois <- function(y, lambda) {
    lambda <- pmax(as.numeric(lambda), 1e-300)
    n <- max(length(y), length(lambda))
    y <- rep(as.numeric(y), length.out = n)
    lambda <- rep(lambda, length.out = n)
    loglam <- log(lambda)
    
    lm <- matrix(-Inf, nrow = n, ncol = length(log_eta))
    for (j in seq_along(log_eta)) {
      lm[, j] <- log_pois_from_lograte(y, loglam + log_eta[j]) + log_quad_w[j]
    }
    
    mm <- apply(lm, 1L, max)
    ans <- rep(-Inf, n)
    finite_rows <- is.finite(mm)
    if (any(finite_rows)) {
      ans[finite_rows] <- mm[finite_rows] +
        log(rowSums(exp(lm[finite_rows, , drop = FALSE] - mm[finite_rows])))
    }
    ans
  }
  
  integrated_cdf <- function(q, lambda) {
    lambda <- pmax(as.numeric(lambda), 1e-300)
    n <- max(length(q), length(lambda))
    q <- rep(as.numeric(q), length.out = n)
    lambda <- rep(lambda, length.out = n)
    out <- numeric(n)
    
    for (j in seq_along(log_eta)) {
      lr <- log(lambda) + log_eta[j]
      pj <- numeric(n)
      valid_idx <- is.finite(lr) & lr < 700
      if (any(valid_idx)) {
        pj[valid_idx] <- stats::ppois(q[valid_idx], lambda = exp(lr[valid_idx]))
      }
      out <- out + quad_w[j] * pj
    }
    pmin(pmax(out, 0), 1)
  }
  
  d_fun <- function(y, par, log = FALSE) {
    lambda <- pmax(as.numeric(par$lambda), 1e-300)
    s <- pmin(pmax(as.numeric(par$s), 1e-12), 1 - 1e-12)
    n <- max(length(y), length(lambda), length(s))
    y <- rep(as.numeric(y), length.out = n)
    lambda <- rep(lambda, length.out = n)
    s <- rep(s, length.out = n)
    
    l0 <- stats::dpois(y, lambda = lambda, log = TRUE)
    li <- log_integrated_pois(y, lambda)
    a0 <- log1p(-s) + l0
    a1 <- log(s) + li
    
    mm <- pmax(a0, a1)
    ld <- mm + log(exp(a0 - mm) + exp(a1 - mm))
    
    if (log) ld else exp(ld)
  }
  
  p_fun <- function(y, par, ...) {
    lambda <- pmax(as.numeric(par$lambda), 1e-300)
    s <- pmin(pmax(as.numeric(par$s), 1e-12), 1 - 1e-12)
    n <- max(length(y), length(lambda), length(s))
    y <- rep(as.numeric(y), length.out = n)
    lambda <- rep(lambda, length.out = n)
    s <- rep(s, length.out = n)
    
    p0 <- stats::ppois(y, lambda = lambda)
    pi <- integrated_cdf(y, lambda)
    out <- (1 - s) * p0 + s * pi
    out[y < 0] <- 0
    pmin(pmax(out, 0), 1)
  }
  
  sh_lam <- num_score_hess_factory(d_fun, "log")
  sh_s <- num_score_hess_factory(d_fun, "logit")
  
  f <- list(
    family = "RSB",
    names = c("lambda", "s"),
    links = c("log", "logit"),
    map2par = function(eta) {
      out <- vector("list", length(eta))
      names(out) <- names(eta)
      link_vec <- stats::setNames(c("log", "logit"), c("lambda", "s"))
      for (nm in names(out)) {
        out[[nm]] <- stats::make.link(link_vec[[nm]])$linkinv(eta[[nm]])
      }
      out
    },
    d = d_fun,
    p = p_fun,
    loglik = function(y, par, ...) sum(d_fun(y, par, log = TRUE)),
    score = list(
      lambda = function(y, p, ...) sh_lam$score(y, p, id = "lambda"),
      s = function(y, p, ...) sh_s$score(y, p, id = "s")
    ),
    hess = list(
      lambda = function(y, p, ...) sh_lam$hess(y, p, id = "lambda"),
      s = function(y, p, ...) sh_s$hess(y, p, id = "s")
    )
  )
  class(f) <- "family.bamlss"
  f
}

# Validate a family for the model-search engine
validate_model_family <- function(fam, requested_name = NULL) {
  lab <- if (!is.null(requested_name) && nzchar(requested_name)) {
    requested_name
  } else {
    if (!is.null(fam$family)) as.character(fam$family)[1L] else "<family>"
  }
  
  if (!inherits(fam, "family.bamlss")) {
    stop(
      "Family '", lab, "' is not a valid family.bamlss object.", call. = FALSE
    )
  }
  
  params <- fam$names
  if (is.null(params) || !length(params) || any(!nzchar(params))) {
    stop(
      "Family '", lab, "' does not expose a valid `names` vector.",
      call. = FALSE
    )
  }
  
  if (is.null(fam$links) || length(fam$links) != length(params)) {
    stop(
      "Family '", lab,
      "' has an invalid `links` specification for parameters: ",
      paste(params, collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  required <- c("d", "map2par", "loglik")
  missing_core <- required[!vapply(
    required,
    function(z) !is.null(fam[[z]]) && is.function(fam[[z]]),
    logical(1)
  )]
  
  if (length(missing_core)) {
    stop(
      "Family '", lab, "' is incomplete for the model-search target. Missing: ",
      paste(missing_core, collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  score_functions_available <- !is.null(fam$score) && all(vapply(
    params, function(p) !is.null(fam$score[[p]]) && is.function(fam$score[[p]]),
    logical(1)
  ))
  hessian_functions_available <- !is.null(fam$hess) && all(vapply(
    params, function(p) !is.null(fam$hess[[p]]) && is.function(fam$hess[[p]]),
    logical(1)
  ))
  
  if (!score_functions_available || !hessian_functions_available) {
    bad_score <- params[!vapply(
      params,
      function(p) !is.null(fam$score[[p]]) && is.function(fam$score[[p]]),
      logical(1)
    )]
    bad_hess <- params[!vapply(
      params, function(p) !is.null(fam$hess[[p]]) && is.function(fam$hess[[p]]),
      logical(1)
    )]
    stop(
      "Family '", lab,
      "' lacks the score/hess functions the model search needs. ",
      if (length(bad_score)) {
        paste0("score: ", paste(bad_score, collapse = ", "), ". ")
      } else {
        ""
      },
      if (length(bad_hess)) {
        paste0("hess: ", paste(bad_hess, collapse = ", "), ". ")
      } else {
        ""
      },
      call. = FALSE
    )
  }
  
  fam
}

# Retrieve standardized distribution family
get_model_family <- function(family_name) {
  requested_name <- if (is.character(family_name) &&
                        length(family_name) == 1L) {
    family_name
  } else if (is.list(family_name) && !is.null(family_name$family)) {
    as.character(family_name$family)[1L]
  } else {
    "<custom>"
  }
  
  user_supplied <- FALSE
  if (is.list(family_name) && !is.null(family_name$family)) {
    fam <- family_name
    user_supplied <- TRUE
  } else if (identical(family_name, "ZANBI")) {
    fam <- family_zanbi()
  } else if (identical(family_name, "SICHEL")) {
    fam <- family_sichel()
  } else if (identical(family_name, "ZOINB")) {
    fam <- family_zoinb()
  } else if (identical(family_name, "RSB")) {
    fam <- family_rsb()
  } else {
    lookup_name <- family_name
    if (identical(lookup_name, "ZINB")) lookup_name <- "ZINBI"
    
    bamlss_func <- paste0(lookup_name, "_bamlss")
    ns_bamlss <- asNamespace("bamlss")
    ns_gamlss <- asNamespace("gamlss.dist")
    
    if (exists(bamlss_func, envir = ns_bamlss, inherits = FALSE)) {
      fam <- get(bamlss_func, envir = ns_bamlss, inherits = FALSE)()
    } else if (exists(lookup_name, envir = ns_gamlss, inherits = FALSE)) {
      g_fun <- get(lookup_name, envir = ns_gamlss, inherits = FALSE)
      fam <- tryCatch(
        bamlss::gF(g_fun),
        error = function(e) bamlss::bamlss.family(g_fun)
      )
    } else if (exists(lookup_name, envir = .GlobalEnv, inherits = FALSE)) {
      custom_fun <- get(lookup_name, envir = .GlobalEnv, inherits = FALSE)
      fam <- if (is.function(custom_fun)) custom_fun() else custom_fun
      user_supplied <- TRUE
    } else {
      stop("Family '", requested_name, "' not found.", call. = FALSE)
    }
  }
  
  native_loglik <- is.function(fam$loglik)
  ns_bamlss <- asNamespace("bamlss")
  if (exists("complete.bamlss.family", envir = ns_bamlss, inherits = FALSE)) {
    fam <- get(
      "complete.bamlss.family", envir = ns_bamlss, inherits = FALSE
    )(fam)
  }
  
  # Every family gets the observation-safe density and the strict sum as
  # log-likelihood; only a log-likelihood supplied by the user is kept.
  raw_d <- fam$d
  fam$d <- safe_density(raw_d)
  if (!(native_loglik && user_supplied)) fam$loglik <- strict_loglik(fam$d)
  fam <- validate_model_family(fam, requested_name)
  fam <- guard_derivatives(fam)
  
  probe_key <- paste(c(requested_name, fam$names, fam$links), collapse = "|")
  probe <- family_probe_cache[[probe_key]]
  if (is.null(probe)) {
    probe <- probe_family_density(fam, raw_d, fam$d)
    family_probe_cache[[probe_key]] <- probe
    if (probe$batch_sensitive) {
      message(
        "  Family ", requested_name, ": the density from gamlss.dist ",
        "changes with batching; observation-safe evaluation is used."
      )
    }
    if (!probe$safe_consistent || probe$usable < 0.5) {
      warning(
        "Family ", requested_name, ": the density is not finite or not ",
        "batch-invariant on the probe design (usable ",
        sprintf("%.0f%%", 100 * probe$usable), ").", call. = FALSE
      )
    }
  }
  # SICHEL is evaluated in plain R now, but results written with the
  # gamlss.dist density were affected by batching and are judged as before.
  attr(fam, "batch_sensitive") <- probe$batch_sensitive ||
    identical(requested_name, "SICHEL")
  fam
}

# Resolve an offset for one distribution family
resolve_family_offset <- function(offsets, family_name) {
  if (is.null(offsets)) return(NULL)
  if (!is.list(offsets)) return(offsets)
  
  nm <- names(offsets)
  if (!is.null(nm) && family_name %in% nm) return(offsets[[family_name]])
  
  family_keyed <- length(offsets) > 0L && all(vapply(
    offsets, function(z) is.null(z) || is.list(z), logical(1)
  ))
  
  if (family_keyed) return(NULL)
  offsets
}

# Obtain a link object compatible with BAMLSS
link_object <- function(link) {
  ns <- asNamespace("bamlss")
  if (exists("make.link2", envir = ns, inherits = FALSE)) {
    return(get("make.link2", envir = ns, inherits = FALSE)(link))
  }
  stats::make.link(link)
}

# Formula update
update_formula_terms <- function(base_f, add_terms) {
  if (length(add_terms) == 0 || add_terms == "" || add_terms == "1") {
    return(base_f)
  }
  if (is.null(base_f)) base_f <- as.formula("~ 1")
  
  base_str <- paste(deparse(base_f, width.cutoff = 500), collapse = " ")
  
  if (grepl("~\\s*1\\s*$", base_str)) {
    new_str <- sub("1\\s*$", add_terms, base_str)
  } else {
    new_str <- paste(base_str, "+", add_terms)
  }
  
  new_f <- as.formula(new_str)
  environment(new_f) <- environment(base_f)
  return(new_f)
}

# Parse a distributional formula with explicit BAMLSS specials
parse_distributional_formula <- function(formula, fam_obj) {
  ns <- asNamespace("bamlss")
  if (!exists("lin", envir = ns, inherits = FALSE)) {
    stop(
      "The installed bamlss version does not provide the `lin()` constructor. ",
      "Update bamlss before running the model search.",
      call. = FALSE
    )
  }
  bamlss::bamlss.formula(formula = formula, family = fam_obj, specials = "lin")
}

# Aggregate hierarchical component PIPs
aggregate_component_pips <- function(component_pip) {
  component_pip <- as.matrix(component_pip)
  rn <- rownames(component_pip)
  cn <- colnames(component_pip)
  mp <- candidate_map(rn)
  vars <- unique(mp$variable)
  
  assoc <- matrix(NA_real_, length(vars), length(cn), dimnames = list(vars, cn))
  nl <- matrix(NA_real_, length(vars), length(cn), dimnames = list(vars, cn))
  nl_given <- matrix(
    NA_real_, length(vars), length(cn), dimnames = list(vars, cn)
  )
  
  for (v in vars) {
    rl <- match(paste0("L::", v), rn)
    rnl <- match(paste0("NL::", v), rn)
    if (!is.na(rl)) assoc[v, ] <- component_pip[rl, ]
    if (!is.na(rnl)) {
      nl[v, ] <- component_pip[rnl, ]
      good <- is.finite(assoc[v, ]) & assoc[v, ] > 0
      nl_given[v, good] <- pmin(1, pmax(0, nl[v, good] / assoc[v, good]))
    }
  }
  list(
    association = assoc,
    nonlinear = nl,
    nonlinear_given_association = nl_given
  )
}

# Aggregate hierarchical indicator traces
aggregate_trace_indicators <- function(trace_inclusion) {
  if (!length(trace_inclusion)) return(list())
  
  first <- trace_inclusion[[1L]]
  rn <- rownames(first)
  cn <- colnames(first)
  mp <- candidate_map(rn)
  vars <- unique(mp$variable)
  
  make_one <- function(z, component = c("association", "nonlinear")) {
    component <- match.arg(component)
    out <- matrix(FALSE, length(vars), length(cn), dimnames = list(vars, cn))
    for (v in vars) {
      key <- paste0(if (component == "association") "L::" else "NL::", v)
      r <- match(key, rn)
      if (!is.na(r)) out[v, ] <- z[r, ]
    }
    out
  }
  
  list(
    association = lapply(trace_inclusion, make_one, component = "association"),
    nonlinear = lapply(trace_inclusion, make_one, component = "nonlinear")
  )
}

# ESS and MCSE for an indicator trace
indicator_trace_diagnostics <- function(trace_list) {
  if (!length(trace_list)) return(list(ess = numeric(0), mcse = numeric(0)))
  
  first <- trace_list[[1L]]
  X <- do.call(rbind, lapply(trace_list, function(z) as.numeric(z)))
  nm <- as.vector(outer(rownames(first), colnames(first), paste, sep = "::"))
  
  if (length(nm) == ncol(X)) colnames(X) <- nm
  ess <- vapply(
    seq_len(ncol(X)),
    function(j) effective_sample_size(X[, j]),
    numeric(1)
  )
  names(ess) <- colnames(X)
  
  pp <- colMeans(X)
  mcse <- sqrt(pp * (1 - pp) / ess)
  bm <- batch_means_stats(X)
  mcse <- pmax(mcse, bm$mcse, na.rm = TRUE)
  ess <- pmin(ess, bm$ess, na.rm = TRUE)
  mcse[!is.finite(mcse)] <- NA_real_
  ess[!is.finite(ess)] <- NA_real_
  list(ess = ess, mcse = mcse)
}

# Candidate component key
candidate_key <- function(variable, component = c("linear", "nonlinear")) {
  component <- match.arg(component)
  paste0(if (component == "linear") "L::" else "NL::", variable)
}

# Parse a candidate component key
candidate_specification <- function(key, data = NULL) {
  key <- as.character(key)[1L]
  if (startsWith(key, "L::")) {
    variable <- substring(key, 4L)
    component <- "linear"
  } else if (startsWith(key, "NL::")) {
    variable <- substring(key, 5L)
    component <- "nonlinear"
  } else {
    variable <- key
    if (!is.null(data) && variable %in% names(data)) {
      x <- data[[variable]]
      cond <- is.numeric(x) && sum(is.finite(x)) > 0L &&
        length(unique(x[is.finite(x)])) > 7L
      component <- if (cond) "legacy_smooth" else "linear"
    } else {
      component <- "linear"
    }
  }
  list(key = key, variable = variable, component = component)
}

# Expand candidate pool into hierarchical components
expand_candidate_covariates <- function(candidate_covariates, data) {
  unique(unlist(lapply(candidate_covariates, function(v) {
    sp <- candidate_specification(v, data)
    if (startsWith(v, "L::") || startsWith(v, "NL::")) return(v)
    if (!sp$variable %in% names(data)) {
      stop("Candidate variable not found: ", sp$variable, call. = FALSE)
    }
    
    x <- data[[sp$variable]]
    continuous <- is.numeric(x) && sum(is.finite(x)) > 0L &&
      length(unique(x[is.finite(x)])) > 7L
    
    if (continuous) {
      c(
        candidate_key(sp$variable, "linear"),
        candidate_key(sp$variable, "nonlinear")
      )
    } else {
      candidate_key(sp$variable, "linear")
    }
  }), use.names = FALSE))
}

# Candidate component map
candidate_map <- function(candidate_keys, data = NULL) {
  z <- lapply(candidate_keys, candidate_specification, data = data)
  data.frame(
    key = vapply(z, `[[`, character(1), "key"),
    variable = vapply(z, `[[`, character(1), "variable"),
    component = vapply(z, `[[`, character(1), "component"),
    stringsAsFactors = FALSE
  )
}

# Basis dimension for a nonlinear candidate
nonlinear_basis_dimension <- function(x, param_name = "mu", min_k = 3L) {
  if (!is.numeric(x)) return(NA_integer_)
  z <- x[is.finite(x)]
  n_unique <- length(unique(z))
  if (n_unique < (as.integer(min_k) + 1L)) return(NA_integer_)
  
  target_k <- if (param_name %in% c("mu", "lambda")) 10L else 4L
  kk <- min(as.integer(target_k), as.integer(n_unique - 1L))
  
  if (!is.finite(kk) || kk < as.integer(min_k)) return(NA_integer_)
  as.integer(kk)
}

# Check hierarchical allocation invariant
hierarchy_valid <- function(inclusion) {
  rn <- rownames(inclusion)
  if (is.null(rn)) return(TRUE)
  
  nl <- which(startsWith(rn, "NL::"))
  if (!length(nl)) return(TRUE)
  
  for (r in nl) {
    v <- substring(rn[r], 5L)
    rp <- match(paste0("L::", v), rn)
    if (is.na(rp)) return(FALSE)
    if (any(inclusion[r, ] & !inclusion[rp, ])) return(FALSE)
  }
  TRUE
}

# Hierarchy-admissible birth mask
birth_mask <- function(inclusion, admissibility) {
  mask <- !inclusion & admissibility
  rn <- rownames(inclusion)
  nl <- which(startsWith(rn, "NL::"))
  
  if (length(nl)) {
    for (r in nl) {
      rp <- match(paste0("L::", substring(rn[r], 5L)), rn)
      if (is.na(rp)) {
        mask[r, ] <- FALSE
      } else {
        mask[r, ] <- mask[r, ] & inclusion[rp, ]
      }
    }
  }
  mask
}

# Hierarchy-admissible death mask
death_mask <- function(inclusion, admissibility) {
  mask <- inclusion & admissibility
  rn <- rownames(inclusion)
  lr <- which(startsWith(rn, "L::"))
  
  if (length(lr)) {
    for (r in lr) {
      rn_nl <- paste0("NL::", substring(rn[r], 4L))
      rc <- match(rn_nl, rn)
      if (!is.na(rc)) mask[r, ] <- mask[r, ] & !inclusion[rc, ]
    }
  }
  mask
}

# Hierarchical model dimensions
model_dimensions <- function(inclusion) {
  rn <- rownames(inclusion)
  if (is.null(rn)) {
    return(list(
      distinct = sum(rowSums(inclusion) > 0L),
      associations = sum(inclusion),
      nonlinear = 0L,
      components = sum(inclusion),
      p = nrow(inclusion),
      K = ncol(inclusion)
    ))
  }
  
  mp <- candidate_map(rn)
  vars <- unique(mp$variable)
  lin_rows <- which(mp$component != "nonlinear")
  nl_rows <- which(mp$component == "nonlinear")
  
  var_active <- vapply(vars, function(v) {
    rr <- which(mp$variable == v)
    any(inclusion[rr, , drop = FALSE])
  }, logical(1))
  
  list(
    distinct = sum(var_active),
    associations = if (length(lin_rows)) {
      sum(inclusion[lin_rows, , drop = FALSE])
    } else {
      0L
    },
    nonlinear = if (length(nl_rows)) {
      sum(inclusion[nl_rows, , drop = FALSE])
    } else {
      0L
    },
    components = sum(inclusion),
    p = length(vars),
    K = ncol(inclusion)
  )
}

# Prepare candidate covariates
auto_format_covariates <- function(
    vars, data, param_name = "mu",
    c_linear = 1.5,
    a_lambda = 1,
    b_lambda = 1) {
  
  if (!length(vars)) return("1")
  if (!is.finite(c_linear) || c_linear <= 0) {
    stop("`c_linear` must be positive.", call. = FALSE)
  }
  if (!is.finite(a_lambda) ||
      a_lambda <= 0 ||
      !is.finite(b_lambda) ||
      b_lambda <= 0) {
    stop("`a_lambda` and `b_lambda` must be positive.", call. = FALSE)
  }
  
  formatted <- vapply(vars, function(key) {
    sp <- candidate_specification(key, data)
    v <- sp$variable
    if (!v %in% names(data)) {
      stop("Candidate variable not found: ", v, call. = FALSE)
    }
    
    x <- data[[v]]
    continuous <- is.numeric(x) && sum(is.finite(x)) > 0L &&
      length(unique(x[is.finite(x)])) > 7L
    
    if (sp$component == "linear") {
      sprintf('lin(~ %s, ridge = FALSE)', v)
    } else if (sp$component == "nonlinear") {
      if (!continuous) {
        stop(
          "Nonlinear candidate requested for a non-continuous variable: ", v,
          call. = FALSE
        )
      }
      kk <- nonlinear_basis_dimension(x, param_name = param_name)
      if (is.na(kk)) {
        stop(
          "Too few distinct values for a smooth of ",
          v, " -> ", param_name, ".", call. = FALSE
        )
      }
      sprintf(
        paste0(
          's(%s, bs = "tp", k = %d, m = c(2, 0), ',
          'xt = list(prior = "ig", a = %.17g, b = %.17g))'
        ),
        v, kk, a_lambda, b_lambda
      )
    } else if (continuous) {
      kk <- nonlinear_basis_dimension(x, param_name = param_name)
      if (is.na(kk)) {
        stop(
          "Too few distinct values for a smooth of ",
          v, " -> ", param_name, ".", call. = FALSE
        )
      }
      sprintf(
        paste0(
          's(%s, bs = "ps", k = %d, ',
          'xt = list(prior = "ig", a = %.17g, b = %.17g))'
        ),
        v, kk, a_lambda, b_lambda
      )
    } else {
      sprintf('lin(~ %s, ridge = FALSE)', v)
    }
  }, character(1))
  
  paste(formatted, collapse = " + ")
}

# Prepare candidate data
prepare_candidate_data <- function(data, vars) {
  mp <- candidate_map(vars, data)
  base_vars <- unique(mp$variable)
  missing_vars <- setdiff(base_vars, names(data))
  
  if (length(missing_vars)) {
    stop(
      "Candidate variables not found in data: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
  
  out <- data
  for (v in base_vars) {
    if (is.character(out[[v]])) out[[v]] <- factor(out[[v]], exclude = NULL)
  }
  out
}

# Construct a complete district-month count panel
build_district_month_panel <- function(
    data, space_col = "admin2", month_col = "year_month", y_col = "fatalities",
    cell_covariates, static_by_space = character(0),
    require_complete_panel = TRUE) {
  
  data <- as.data.frame(data)
  required <- unique(c(space_col, month_col, y_col, cell_covariates))
  miss <- setdiff(required, names(data))
  
  if (length(miss)) {
    stop(
      "District-month panel construction is missing required columns: ",
      paste(miss, collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(cell_covariates)) {
    stop("`cell_covariates` contains duplicated names.", call. = FALSE)
  }
  
  if (length(setdiff(static_by_space, cell_covariates))) {
    stop(
      "`static_by_space` must be a subset of `cell_covariates`.", call. = FALSE
    )
  }
  
  parse_month <- function(x) {
    if (inherits(x, "Date")) return(as.Date(format(x, "%Y-%m-01")))
    z <- as.character(x)
    out <- suppressWarnings(
      as.Date(ifelse(nchar(z) >= 7L, paste0(substr(z, 1L, 7L), "-01"), z))
    )
    out
  }
  
  month <- parse_month(data[[month_col]])
  if (anyNA(month)) {
    bad <- unique(as.character(data[[month_col]][is.na(month)]))
    stop(
      "Cannot parse district-month field `", month_col, "`; examples: ",
      paste(utils::head(bad, 5L), collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  space <- as.character(data[[space_col]])
  if (anyNA(space) || any(!nzchar(space))) {
    stop("District identifier contains missing/empty values.", call. = FALSE)
  }
  
  y <- suppressWarnings(as.numeric(data[[y_col]]))
  if (any(!is.finite(y)) || any(y < 0) || any(abs(y - round(y)) > 1e-8)) {
    stop(
      "District-month response must be finite, non-negative integer counts.",
      call. = FALSE
    )
  }
  
  work <- data
  work$.space_id <- space
  work$.month <- month
  work$.count <- as.integer(round(y))
  grp <- split(
    seq_len(nrow(work)), interaction(work$.space_id, work$.month, drop = TRUE)
  )
  
  one_value <- function(x, idx, v, gname) {
    z <- x[idx]
    keep <- !is.na(z)
    z <- z[keep]
    if (!length(z)) return(NA)
    if (is.factor(z)) z <- as.character(z)
    if (inherits(z, "Date")) z <- as.character(z)
    
    if (is.numeric(z)) {
      zz <- z[is.finite(z)]
      if (!length(zz)) return(NA_real_)
      ref <- zz[1L]
      tol <- sqrt(.Machine$double.eps) * max(1, abs(ref))
      if (any(abs(zz - ref) > tol)) {
        stop(
          "Covariate `", v, "` is not constant within district-month cell ",
          gname,
          "; event-level variables cannot be collapsed to one value.",
          call. = FALSE
        )
      }
      return(ref)
    }
    
    uz <- unique(as.character(z))
    if (length(uz) > 1L) {
      stop(
        "Covariate `", v, "` is not constant within district-month cell ",
        gname,
        "; event-level variables cannot be collapsed to one value.",
        call. = FALSE
      )
    }
    uz[1L]
  }
  
  observed <- lapply(names(grp), function(gname) {
    idx <- grp[[gname]]
    ans <- list(
      .space_id = work$.space_id[idx[1L]],
      .month = work$.month[idx[1L]],
      .count = sum(work$.count[idx], na.rm = TRUE)
    )
    for (v in cell_covariates) ans[[v]] <- one_value(work[[v]], idx, v, gname)
    as.data.frame(ans, stringsAsFactors = FALSE, check.names = FALSE)
  })
  
  observed <- do.call(rbind, observed)
  observed$.month <- as.Date(observed$.month)
  observed$.count <- as.integer(observed$.count)
  
  spaces <- sort(unique(work$.space_id))
  months <- seq(min(work$.month), max(work$.month), by = "month")
  
  grid <- expand.grid(
    .space_id = spaces,
    .month = months,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  grid$.month <- as.Date(grid$.month, origin = "1970-01-01")
  out <- merge(
    grid, observed, by = c(".space_id", ".month"), all.x = TRUE, sort = TRUE
  )
  out$.count[is.na(out$.count)] <- 0L
  
  for (v in static_by_space) {
    lookup <- lapply(spaces, function(s) {
      z <- work[[v]][work$.space_id == s]
      z <- z[!is.na(z)]
      if (is.factor(z)) z <- as.character(z)
      
      if (is.numeric(z)) {
        z <- z[is.finite(z)]
        if (!length(z)) return(NA_real_)
        ref <- z[1L]
        tol <- sqrt(.Machine$double.eps) * max(1, abs(ref))
        if (any(abs(z - ref) > tol)) {
          stop(
            "Variable `", v,
            "` was declared district-static but varies within district ", s,
            ".", call. = FALSE
          )
        }
        return(ref)
      }
      
      uz <- unique(as.character(z))
      if (length(uz) > 1L) {
        stop(
          "Variable `", v,
          "` was declared district-static but varies within district ", s, ".",
          call. = FALSE
        )
      }
      if (length(uz)) uz[1L] else NA_character_
    })
    
    names(lookup) <- spaces
    miss_v <- is.na(out[[v]])
    if (any(miss_v)) {
      vals <- unlist(lookup[out$.space_id[miss_v]], use.names = FALSE)
      if (is.numeric(work[[v]])) vals <- suppressWarnings(as.numeric(vals))
      out[[v]][miss_v] <- vals
    }
  }
  
  missing_by_var <- vapply(
    cell_covariates, function(v) sum(is.na(out[[v]])), integer(1)
  )
  
  if (isTRUE(require_complete_panel) && any(missing_by_var > 0L)) {
    bad <- missing_by_var[missing_by_var > 0L]
    stop(
      "After adding the zero-event months, covariates are missing in: ",
      paste(paste0(names(bad), "=", bad, " cell(s)"), collapse = ", "),
      ". Supply a complete district-month covariate panel (nothing is ",
      "imputed).",
      call. = FALSE
    )
  }
  
  names(out)[names(out) == ".space_id"] <- space_col
  names(out)[names(out) == ".count"] <- y_col
  out[[month_col]] <- format(out$.month, "%Y-%m")
  out$date_obj <- out$.month
  out$.month <- NULL
  out <- out[order(out[[space_col]], out$date_obj), , drop = FALSE]
  rownames(out) <- NULL
  
  if (anyDuplicated(out[c(space_col, month_col)])) {
    stop(
      "Internal error: district-month constructor produced duplicate cells.",
      call. = FALSE
    )
  }
  
  attr(out, "analysis_unit") <- "district_month"
  attr(out, "zero_count") <- sum(out[[y_col]] == 0L)
  attr(out, "complete_grid") <- !any(missing_by_var > 0L)
  out
}

# Validate the common analysis sample
check_analysis_sample <- function(
    data, formula_base, candidate_covariates, offset = NULL) {
  if (!is.data.frame(data)) data <- as.data.frame(data)
  n <- nrow(data)
  if (n < 1L) stop("Analysis data contain no rows.", call. = FALSE)
  
  cand_vars <- unique(candidate_map(candidate_covariates, data)$variable)
  missing_cand <- setdiff(cand_vars, names(data))
  if (length(missing_cand)) {
    stop(
      "Candidate variables not found in the analysis data: ",
      paste(missing_cand, collapse = ", "),
      call. = FALSE
    )
  }
  
  structural_vars <- unique(unlist(lapply(formula_base, function(f) {
    if (is.null(f)) character(0) else all.vars(f)
  }), use.names = FALSE))
  
  vars <- unique(c(intersect(structural_vars, names(data)), cand_vars))
  
  bad <- vapply(vars, function(v) {
    x <- data[[v]]
    if (is.numeric(x)) {
      sum(!is.finite(x))
    } else {
      sum(is.na(x))
    }
  }, integer(1))
  bad <- sort(bad[bad > 0L], decreasing = TRUE)
  
  if (length(bad)) {
    stop(
      "All candidate models must use the same rows, but these variables ",
      "have missing or non-finite values: ",
      paste(sprintf("%s=%d", names(bad), bad), collapse = ", "),
      ". Define a complete-case sample before fitting.",
      call. = FALSE
    )
  }
  
  if (!is.null(offset)) {
    offs <- if (is.list(offset)) offset else list(offset = offset)
    if (is.null(names(offs))) names(offs) <- paste0("offset", seq_along(offs))
    bad_offset <- character(0)
    
    for (nm in names(offs)) {
      z <- offs[[nm]]
      if (is.null(z)) next
      z <- suppressWarnings(as.numeric(z))
      
      if (!(length(z) %in% c(1L, n))) {
        bad_offset <- c(bad_offset, paste0(nm, "[length=", length(z), "]"))
      } else if (any(!is.finite(z))) {
        bad_offset <- c(
          bad_offset, paste0(nm, "[non-finite=", sum(!is.finite(z)), "]")
        )
      }
    }
    
    if (length(bad_offset)) {
      stop(
        "Offset validation failed: ", paste(bad_offset, collapse = ", "),
        ". Each offset must be finite and have length 1 or nrow(data).",
        call. = FALSE
      )
    }
  }
  
  invisible(TRUE)
}

# Model-state key
state_key <- function(inclusion) {
  paste0(as.integer(inclusion), collapse = "")
}

# Stable numerical utilities
symmetrize_matrix <- function(A) {
  (A + t(A)) / 2
}

make_positive_definite <- function(A, floor = 1e-8) {
  A <- symmetrize_matrix(as.matrix(A))
  if (!nrow(A)) return(A)
  ee <- eigen(A, symmetric = TRUE)
  vals <- pmax(ee$values, floor)
  ee$vectors %*% (diag(vals, nrow = length(vals)) %*% t(ee$vectors))
}

invert_positive_definite <- function(A, floor = 1e-8) {
  A <- symmetrize_matrix(as.matrix(A))
  if (!nrow(A)) return(A)
  ee <- eigen(A, symmetric = TRUE)
  vals <- pmax(ee$values, floor)
  ee$vectors %*% (diag(1 / vals, nrow = length(vals)) %*% t(ee$vectors))
}

logdet_positive_definite <- function(A, floor = 1e-12) {
  if (!nrow(A)) return(0)
  vals <- eigen(
    symmetrize_matrix(as.matrix(A)),
    symmetric = TRUE,
    only.values = TRUE
  )$values
  sum(log(pmax(vals, floor)))
}

draw_multivariate_normal <- function(mean, cov, floor = 1e-10) {
  mean <- as.numeric(mean)
  if (!length(mean)) return(numeric(0))
  ee <- eigen(symmetrize_matrix(as.matrix(cov)), symmetric = TRUE)
  vals <- pmax(ee$values, floor)
  as.numeric(mean + ee$vectors %*% (sqrt(vals) * stats::rnorm(length(mean))))
}

log_multivariate_normal_density <- function(x, mean, cov, floor = 1e-10) {
  x <- as.numeric(x)
  mean <- as.numeric(mean)
  d <- length(x)
  if (!d || length(mean) != d) return(-Inf)
  
  ee <- eigen(symmetrize_matrix(as.matrix(cov)), symmetric = TRUE)
  vals <- pmax(ee$values, floor)
  z <- crossprod(ee$vectors, x - mean)
  -0.5 * (d * log(2 * pi) + sum(log(vals)) + sum((z^2) / vals))
}

# Predictor-scale standardization of a Gaussian RJ proposal
standardize_gaussian_on_predictor <- function(
    mean, cov, X, target_rms = 1, floor = 1e-10) {
  
  mean <- as.numeric(mean)
  cov <- make_positive_definite(as.matrix(cov), floor)
  X <- as.matrix(X)
  d <- length(mean)
  
  if (!d || ncol(X) != d || any(dim(cov) != c(d, d))) {
    stop(
      "Incompatible dimensions in the proposal standardization.", call. = FALSE
    )
  }
  if (!is.finite(target_rms) || target_rms <= 0) {
    stop(
      "`proposal_link_rms` must be finite and strictly positive.", call. = FALSE
    )
  }
  
  eta_mean <- as.numeric(X %*% mean)
  XV <- X %*% cov
  eta_var <- rowSums(XV * X)
  eta_var[!is.finite(eta_var)] <- Inf
  eta_var <- pmax(eta_var, 0)
  
  second_moment <- eta_mean^2 + eta_var
  raw_rms <- sqrt(base::mean(second_moment))
  if (!is.finite(raw_rms)) raw_rms <- Inf
  
  scale <- if (raw_rms > target_rms) target_rms / raw_rms else 1
  scale <- as.numeric(scale)
  if (!is.finite(scale) || scale <= 0) {
    stop(
      "Failed to construct a positive predictor-scale proposal factor.",
      call. = FALSE
    )
  }
  
  mean2 <- scale * mean
  cov2 <- make_positive_definite((scale^2) * cov, floor)
  
  eta_mean2 <- as.numeric(X %*% mean2)
  XV2 <- X %*% cov2
  eta_var2 <- pmax(rowSums(XV2 * X), 0)
  scaled_rms <- sqrt(base::mean(eta_mean2^2 + eta_var2))
  
  list(
    mean = mean2,
    cov = cov2,
    scale = scale,
    raw_rms = raw_rms,
    scaled_rms = scaled_rms
  )
}

# BAMLSS internal function lookup
get_bamlss_internal <- function(name) {
  ns <- asNamespace("bamlss")
  if (!exists(name, envir = ns, inherits = FALSE)) {
    stop(
      "Required bamlss internal function not available: ", name, call. = FALSE
    )
  }
  get(name, envir = ns, inherits = FALSE)
}

# Dynamic-allocation admissibility matrix
build_admissibility_matrix <- function(
    formula_base, candidate_covariates, params, data,
    external_admissibility = NULL) {
  
  E <- matrix(
    TRUE, nrow = length(candidate_covariates), ncol = length(params),
    dimnames = list(candidate_covariates, params)
  )
  mp <- candidate_map(candidate_covariates, data)
  
  overlap <- character(0)
  for (k in params) {
    fv <- if (!is.null(formula_base[[k]])) {
      all.vars(formula_base[[k]])
    } else {
      character(0)
    }
    hit <- unique(mp$variable[mp$variable %in% fv])
    if (length(hit)) overlap <- c(overlap, paste0(hit, "::", k))
  }
  
  if (length(overlap)) {
    stop(
      "Dynamic candidate(s) also appear in structural predictors: ",
      paste(unique(overlap), collapse = ", "),
      ". Remove them from the structural formula or from the candidates.",
      call. = FALSE
    )
  }
  
  for (j in seq_along(candidate_covariates)) {
    v <- mp$variable[j]
    x <- data[[v]]
    n_unique <- if (is.numeric(x)) {
      length(unique(x[is.finite(x)]))
    } else {
      length(unique(x[!is.na(x)]))
    }
    
    if (n_unique < 2L) E[j, ] <- FALSE
    
    if (mp$component[j] == "nonlinear") {
      continuous <- is.numeric(x) && sum(is.finite(x)) > 0L && n_unique > 7L
      if (!continuous) {
        E[j, ] <- FALSE
      } else {
        for (kk in seq_along(params)) {
          E[j, kk] <- isTRUE(E[j, kk]) &&
            !is.na(nonlinear_basis_dimension(x, param_name = params[kk]))
        }
      }
    }
  }
  
  if (!is.null(external_admissibility)) {
    X <- as.matrix(external_admissibility)
    if (is.null(rownames(X)) || is.null(colnames(X))) {
      stop(
        "`external_admissibility` must have row and column names.",
        call. = FALSE
      )
    }
    
    param_alias <- c(lambda = "mu", s = "nu")
    missing_vars <- setdiff(unique(mp$variable), rownames(X))
    
    if (length(missing_vars)) {
      stop(
        "External admissibility is missing candidate variable(s): ",
        paste(missing_vars, collapse = ", "),
        call. = FALSE
      )
    }
    
    for (j in seq_len(nrow(mp))) {
      v <- mp$variable[j]
      for (kk in seq_along(params)) {
        p <- params[kk]
        p_external <- if (p %in% colnames(X)) p else unname(param_alias[p])
        
        if (is.na(p_external) || !p_external %in% colnames(X)) {
          stop(
            "External admissibility has no column for family parameter '", p,
            "'.",
            call. = FALSE
          )
        }
        E[j, kk] <- isTRUE(E[j, kk]) && isTRUE(as.logical(X[v, p_external]))
      }
    }
  }
  
  E
}

# Model prior
log_model_prior <- function(inclusion, a_V = 0.5, a_A = 0.2) {
  if (!nrow(inclusion) || !ncol(inclusion)) return(0)
  if (!hierarchy_valid(inclusion)) return(-Inf)
  
  dm <- model_dimensions(inclusion)
  if (dm$p <= 0L || dm$K <= 0L) return(0)
  
  -a_V * dm$distinct * log(dm$p) -
    a_A * dm$components * log(dm$p * dm$K)
}

# Exact prior ratio for a birth candidate
log_model_prior_ratio <- function(inclusion, row, col, a_V = 0.5, a_A = 0.2) {
  if (isTRUE(inclusion[row, col])) return(-Inf)
  plus <- inclusion
  plus[row, col] <- TRUE
  log_model_prior(plus, a_V, a_A) - log_model_prior(inclusion, a_V, a_A)
}

# Prior complexity summary for the allocation matrix
model_prior_summary <- function(
    p, K, a_V = 0.5, a_A = 0.2, admissibility = NULL) {
  K <- as.integer(K)
  
  if (is.null(admissibility)) {
    p <- as.integer(p)
    if (p <= 0L || K <= 0L) {
      return(list(
        p_empty = 1, expected_distinct = 0, expected_allocations = 0,
        expected_associations = 0, expected_nonlinear = 0
      ))
    }
    Kj <- rep(K, p)
    r <- (p * K)^(-a_A)
    row_extra <- p^(-a_V) * ((1 + r)^Kj - 1)
    Zrow <- 1 + row_extra
    ea <- sum(p^(-a_V) * Kj * r * (1 + r)^(pmax(Kj - 1, 0)) / Zrow)
    
    return(list(
      p_empty = prod(1 / Zrow),
      expected_distinct = sum(row_extra / Zrow),
      expected_allocations = ea,
      expected_associations = ea,
      expected_nonlinear = 0
    ))
  }
  
  E <- as.matrix(admissibility)
  rn <- rownames(E)
  if (is.null(rn)) {
    stop(
      "Hierarchical prior summary requires admissibility row names.",
      call. = FALSE
    )
  }
  
  mp <- candidate_map(rn)
  vars <- unique(mp$variable)
  p0 <- length(vars)
  if (p0 <= 0L || K <= 0L) {
    return(list(
      p_empty = 1, expected_distinct = 0, expected_allocations = 0,
      expected_associations = 0, expected_nonlinear = 0
    ))
  }
  
  r <- (p0 * K)^(-a_A)
  p_empty_rows <- e_dist <- e_assoc <- e_nl <- e_comp <- numeric(length(vars))
  
  for (j in seq_along(vars)) {
    v <- vars[j]
    rl <- match(paste0("L::", v), rn)
    rnl <- match(paste0("NL::", v), rn)
    states_by_k <- vector("list", K)
    
    for (k in seq_len(K)) {
      canL <- !is.na(rl) && isTRUE(E[rl, k])
      canNL <- !is.na(rnl) && isTRUE(E[rnl, k])
      states_by_k[[k]] <- if (!canL) 0L else if (canNL) 0:2 else 0:1
    }
    
    grid <- expand.grid(
      states_by_k, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    mat <- as.matrix(grid)
    assoc <- rowSums(mat >= 1L)
    nl <- rowSums(mat >= 2L)
    comp <- assoc + nl
    nonempty <- comp > 0L
    
    lw <- comp * log(r)
    lw[nonempty] <- lw[nonempty] - a_V * log(p0)
    mm <- max(lw)
    w <- exp(lw - mm)
    pr <- w / sum(w)
    
    p_empty_rows[j] <- sum(pr[!nonempty])
    e_dist[j] <- sum(pr * nonempty)
    e_assoc[j] <- sum(pr * assoc)
    e_nl[j] <- sum(pr * nl)
    e_comp[j] <- sum(pr * comp)
  }
  
  list(
    p_empty = prod(p_empty_rows),
    expected_distinct = sum(e_dist),
    expected_allocations = sum(e_comp),
    expected_associations = sum(e_assoc),
    expected_nonlinear = sum(e_nl)
  )
}

# Build formula list for one allocation state
build_distributional_formula <- function(
    inclusion, formula_base, data, fam_obj,
    c_linear = 2.0, a_lambda = 1, b_lambda = 1) {
  
  out <- formula_base
  for (p_name in fam_obj$names) {
    base_f <- out[[p_name]]
    if (is.null(base_f)) base_f <- stats::as.formula("~ 1")
    
    active <- rownames(inclusion)[which(inclusion[, p_name])]
    if (length(active)) {
      add <- auto_format_covariates(
        active, data, p_name,
        c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda
      )
      base_f <- update_formula_terms(base_f, add)
    }
    out[[p_name]] <- base_f
  }
  out[fam_obj$names]
}

# Selected components for one distributional parameter
selected_components_by_parameter <- function(
    inclusion,
    parameter,
    data) {
  
  inclusion <- as.matrix(inclusion)
  
  if (is.null(rownames(inclusion)) || is.null(colnames(inclusion))) {
    stop("`inclusion` must have row and column names.", call. = FALSE)
  }
  
  if (!parameter %in% colnames(inclusion)) {
    stop(
      "Distributional parameter `", parameter,
      "` is absent from the inclusion matrix.",
      call. = FALSE
    )
  }
  
  active_rows <- which(inclusion[, parameter])
  
  if (!length(active_rows)) {
    return(
      list(
        linear_rows = integer(0),
        linear_keys = character(0),
        linear_variables = character(0),
        nonlinear_rows = integer(0),
        nonlinear_keys = character(0),
        nonlinear_variables = character(0)
      )
    )
  }
  
  keys <- rownames(inclusion)[active_rows]
  mp <- candidate_map(keys, data)
  allowed_components <- c("linear", "nonlinear")
  bad_components <- setdiff(unique(mp$component), allowed_components)
  
  if (length(bad_components)) {
    stop(
      "Selected posterior allocation contains unsupported component type(s): ",
      paste(bad_components, collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  is_linear <- mp$component == "linear"
  is_nonlinear <- mp$component == "nonlinear"
  
  list(
    linear_rows = active_rows[is_linear],
    linear_keys = keys[is_linear],
    linear_variables = mp$variable[is_linear],
    nonlinear_rows = active_rows[is_nonlinear],
    nonlinear_keys = keys[is_nonlinear],
    nonlinear_variables = mp$variable[is_nonlinear]
  )
}

# Build a blocked post-selection formula
build_selected_posterior_formula <- function(
    inclusion,
    formula_base,
    data,
    fam_obj,
    a_lambda = 1,
    b_lambda = 1) {
  
  inclusion <- as.matrix(inclusion)
  
  if (is.null(rownames(inclusion)) || is.null(colnames(inclusion))) {
    stop("`inclusion` must have row and column names.", call. = FALSE)
  }
  
  if (!identical(colnames(inclusion), fam_obj$names)) {
    stop(
      "Selected allocation columns do not match family parameters.",
      call. = FALSE
    )
  }
  
  out <- formula_base
  
  for (p_name in fam_obj$names) {
    base_f <- out[[p_name]]
    if (is.null(base_f)) {
      base_f <- stats::as.formula("~ 1")
    }
    
    part <- selected_components_by_parameter(
      inclusion = inclusion,
      parameter = p_name,
      data = data
    )
    
    additions <- character(0)
    
    if (length(part$linear_variables)) {
      linear_variables <- unique(part$linear_variables)
      inner_formula <- paste(
        deparse(stats::reformulate(linear_variables), width.cutoff = 500L),
        collapse = " "
      )
      additions <- c(
        additions,
        paste0("lin(", inner_formula, ", ridge = FALSE)")
      )
    }
    
    if (length(part$nonlinear_keys)) {
      nonlinear_terms <- auto_format_covariates(
        vars = part$nonlinear_keys,
        data = data,
        param_name = p_name,
        a_lambda = a_lambda,
        b_lambda = b_lambda
      )
      additions <- c(additions, nonlinear_terms)
    }
    
    if (length(additions)) {
      base_f <- update_formula_terms(base_f, paste(additions, collapse = " + "))
    }
    
    out[[p_name]] <- base_f
  }
  
  out[fam_obj$names]
}

# Elliptical slice sampler for a selected linear block
GMCMC_ess_linear_block <- function(
    family,
    theta,
    id,
    eta,
    y,
    data,
    weights = NULL,
    offset = NULL,
    ...) {
  
  if (length(id) != 2L) {
    stop(
      "Elliptical-slice update needs a two-element GMCMC term id.",
      call. = FALSE
    )
  }
  
  p_name <- id[[1L]]
  term_name <- id[[2L]]
  block <- theta[[p_name]][[term_name]]
  
  if (!is.numeric(block) || !length(block) || any(!is.finite(block))) {
    stop(
      "Elliptical-slice block contains invalid coefficient values.",
      call. = FALSE
    )
  }
  
  beta0 <- as.numeric(block)
  names(beta0) <- names(block)
  
  if (is.null(names(beta0)) || any(!nzchar(names(beta0)))) {
    stop("Elliptical-slice block requires named coefficients.", call. = FALSE)
  }
  
  if (any(!startsWith(names(beta0), "b"))) {
    stop(
      "Elliptical-slice linear block contains non-beta parameters.",
      call. = FALSE
    )
  }
  
  prior_sd <- data$xt$ess_prior_sd
  
  if (is.null(prior_sd) ||
      length(prior_sd) != 1L ||
      !is.finite(prior_sd) ||
      prior_sd <= 0) {
    stop(
      "Elliptical-slice block is missing a valid `xt$ess_prior_sd`.",
      call. = FALSE
    )
  }
  
  prior_sd <- as.numeric(prior_sd)
  current_fit <- attr(block, "fitted.values")
  
  if (is.null(current_fit)) {
    current_fit <- data$fit.fun(data$X, beta0)
  }
  current_fit <- as.numeric(current_fit)
  
  if (length(current_fit) != length(y) || any(!is.finite(current_fit))) {
    stop(
      "Elliptical-slice block produced invalid fitted values.", call. = FALSE
    )
  }
  
  eta_work <- eta
  if (!is.null(offset)) {
    for (nm in names(offset)) {
      eta_work[[nm]] <- eta_work[[nm]] + as.numeric(offset[[nm]])
    }
  }
  
  eta_base <- eta_work
  eta_base[[p_name]] <- as.numeric(eta_work[[p_name]]) - current_fit
  
  evaluate_beta <- function(beta) {
    fit <- as.numeric(data$fit.fun(data$X, beta))
    if (length(fit) != length(y) || any(!is.finite(fit))) {
      return(list(loglik = -Inf, fit = fit))
    }
    
    eta_candidate <- eta_base
    eta_candidate[[p_name]] <- eta_base[[p_name]] + fit
    par_candidate <- family$map2par(eta_candidate)
    
    ll <- tryCatch(
      as.numeric(family$loglik(y, par_candidate))[1L],
      error = function(e) -Inf
    )
    
    list(loglik = if (is.finite(ll)) ll else -Inf, fit = fit)
  }
  
  target_log_prior <- as.numeric(data$prior(block))[1L]
  expected_log_prior <- sum(
    stats::dnorm(beta0, mean = 0, sd = prior_sd, log = TRUE)
  )
  prior_tolerance <- 1e-8 * (1 + abs(expected_log_prior))
  
  if (!is.finite(target_log_prior) || !is.finite(expected_log_prior) ||
      abs(target_log_prior - expected_log_prior) > prior_tolerance) {
    stop(
      "Elliptical-slice sampling needs the spherical Gaussian block prior.",
      call. = FALSE
    )
  }
  
  current <- evaluate_beta(beta0)
  if (!is.finite(current$loglik)) {
    stop(
      "Current elliptical-slice state has non-finite log-likelihood.",
      call. = FALSE
    )
  }
  
  auxiliary <- stats::rnorm(length(beta0), mean = 0, sd = prior_sd)
  names(auxiliary) <- names(beta0)
  
  log_slice <- current$loglik + log(stats::runif(1L))
  angle <- stats::runif(1L, min = 0, max = 2 * pi)
  angle_min <- angle - 2 * pi
  angle_max <- angle
  contractions <- 0L
  
  repeat {
    beta_candidate <- beta0 * cos(angle) + auxiliary * sin(angle)
    names(beta_candidate) <- names(beta0)
    
    candidate <- evaluate_beta(beta_candidate)
    
    if (is.finite(candidate$loglik) && candidate$loglik >= log_slice) {
      break
    }
    
    contractions <- contractions + 1L
    
    if (angle < 0) {
      angle_min <- angle
    } else {
      angle_max <- angle
    }
    
    bracket_width <- angle_max - angle_min
    collapse_tolerance <-
      .Machine$double.eps * max(1, abs(angle_min), abs(angle_max))
    
    if (!is.finite(bracket_width) || bracket_width <= collapse_tolerance) {
      stop(
        "Elliptical-slice angular bracket collapsed numerically.", call. = FALSE
      )
    }
    
    angle <- stats::runif(1L, min = angle_min, max = angle_max)
  }
  
  out <- block
  out[] <- beta_candidate
  attr(out, "fitted.values") <- candidate$fit
  
  old_hessian <- attr(block, "hess")
  if (!is.null(old_hessian)) {
    attr(out, "hess") <- old_hessian
  }
  
  list(parameters = out, alpha = 0)
}

# Clean auxiliary ESS samples
remove_ess_auxiliary_samples <- function(model) {
  chains <- extract_sample_chains(model)
  
  cleaned <- lapply(chains, function(ch) {
    M <- as.matrix(ch)
    drop <- grepl("\\.ess_shrink$", colnames(M))
    
    if (!any(drop)) {
      return(coda::mcmc(M))
    }
    
    M <- M[, !drop, drop = FALSE]
    coda::mcmc(M)
  })
  
  model$samples <- coda::mcmc.list(cleaned)
  model
}

# Prior for a block of selected linear effects
make_selected_linear_block_prior <- function(
    c_linear,
    n_beta,
    term_key) {
  
  c_linear <- as.numeric(c_linear)[1L]
  n_beta <- as.integer(n_beta)[1L]
  term_key <- as.character(term_key)[1L]
  
  if (!is.finite(c_linear) || c_linear <= 0) {
    stop("`c_linear` must be finite and strictly positive.", call. = FALSE)
  }
  if (!is.finite(n_beta) || n_beta < 1L) {
    stop("`n_beta` must be a positive integer.", call. = FALSE)
  }
  
  f <- function(parameters) {
    nm <- names(parameters)
    if (is.null(nm)) return(-Inf)
    
    idx <- grep("^b", nm)
    beta <- suppressWarnings(as.numeric(parameters[idx]))
    
    if (length(beta) != n_beta || any(!is.finite(beta))) {
      return(-Inf)
    }
    
    sum(stats::dnorm(beta, mean = 0, sd = c_linear, log = TRUE))
  }
  
  attr(f, "var_prior") <- "proper_candidate_prior"
  attr(f, "candidate_term_key") <- term_key
  attr(f, "candidate_component") <- "linear_block"
  
  f
}

# Find a blocked linear term
find_selected_linear_block <- function(
    frame,
    parameter,
    variables) {
  
  sc <- frame$x[[parameter]]$smooth.construct
  if (is.null(sc) || !length(sc)) {
    stop(
      "Parameter `", parameter,
      "` contains no smooth.construct terms.",
      call. = FALSE
    )
  }
  
  variables <- sort(unique(as.character(variables)))
  
  hits <- names(sc)[
    vapply(
      names(sc),
      function(label) {
        obj <- sc[[label]]
        object_terms <- if (is.null(obj$term)) {
          character(0)
        } else {
          sort(unique(as.character(obj$term)))
        }
        
        is_linear <- startsWith(label, "lin(") || inherits(obj, "linear.smooth")
        isTRUE(is_linear) && identical(object_terms, variables)
      },
      logical(1)
    )
  ]
  
  if (length(hits) != 1L) {
    stop(
      "Could not identify exactly one blocked linear term for parameter `",
      parameter, "`. Expected variables: ", paste(variables, collapse = ", "),
      "; matches: ",
      if (length(hits)) paste(hits, collapse = ", ") else "<none>", ".",
      call. = FALSE
    )
  }
  
  hits[[1L]]
}

# Design columns generated by one linear candidate
selected_linear_design_columns <- function(
    data,
    variable) {
  
  f <- stats::reformulate(variable)
  X <- stats::model.matrix(f, data = data, na.action = stats::na.pass)
  
  keep <- colnames(X) != "(Intercept)"
  X <- X[, keep, drop = FALSE]
  
  if (!ncol(X)) {
    stop(
      "Linear candidate `", variable, "` produced no design columns.",
      call. = FALSE
    )
  }
  
  paste0("b.", colnames(X))
}

# Install priors for the blocked selected posterior
install_selected_posterior_priors <- function(
    frame,
    inclusion,
    formula_base,
    data,
    fam_obj,
    offset,
    runtime,
    c_linear = 1.0,
    c0 = 2.5,
    cplus = 1.0,
    a_lambda = 1,
    b_lambda = 1) {
  
  inclusion <- as.matrix(inclusion)
  meta <- list()
  candidate_term_keys <- character(0)
  
  for (p_name in fam_obj$names) {
    k <- match(p_name, colnames(inclusion))
    if (is.na(k)) {
      stop(
        "Parameter `", p_name, "` is absent from the inclusion matrix.",
        call. = FALSE
      )
    }
    
    part <- selected_components_by_parameter(
      inclusion = inclusion, parameter = p_name, data = data
    )
    
    if (length(part$linear_variables)) {
      variables <- unique(part$linear_variables)
      label <- find_selected_linear_block(
        frame = frame, parameter = p_name, variables = variables
      )
      
      obj <- frame$x[[p_name]]$smooth.construct[[label]]
      X <- as.matrix(obj$X)
      design_names <- colnames(X)
      
      pm <- term_parameter_map(frame = frame, p_name = p_name, label = label)
      term_key <- paste(p_name, label, sep = "::")
      q <- ncol(X)
      
      if (!is.finite(q) || q < 1L) {
        stop("Blocked linear term has no coefficient columns.", call. = FALSE)
      }
      if (!is.finite(c_linear) || c_linear <= 0) {
        stop("`c_linear` must be finite and positive.", call. = FALSE)
      }
      
      obj$prior <- make_selected_linear_block_prior(
        c_linear = c_linear,
        n_beta = q,
        term_key = term_key
      )
      
      if (is.null(obj$xt)) {
        obj$xt <- list()
      }
      
      obj$xt$ess_prior_sd <- as.numeric(c_linear)
      obj$xt$pS <- NULL
      obj$xt$pm <- NULL
      obj$propose <- GMCMC_ess_linear_block
      
      frame$x[[p_name]]$smooth.construct[[label]] <- obj
      candidate_term_keys <- c(candidate_term_keys, term_key)
      
      local_beta_names <- pm$local_names[!startsWith(pm$local_names, "tau2")]
      
      for (jj in seq_along(part$linear_variables)) {
        variable <- part$linear_variables[[jj]]
        expected_columns <- selected_linear_design_columns(
          data = data, variable = variable
        )
        column_index <- match(expected_columns, design_names)
        
        block <- list(
          row = part$linear_rows[[jj]],
          col = k,
          variable = variable,
          candidate_key = part$linear_keys[[jj]],
          component = "linear",
          parameter = p_name,
          raw_label = label,
          label = label,
          term_key = term_key,
          smooth = FALSE,
          c_linear = c_linear,
          a_lambda = a_lambda,
          b_lambda = b_lambda,
          X = X[, column_index, drop = FALSE],
          design_colnames = design_names[column_index],
          S = list(),
          parameter_names = pm$beta_names[column_index],
          beta_names = pm$beta_names[column_index],
          tau_names = character(0),
          local_parameter_names = local_beta_names[column_index],
          fixed = pm$fixed,
          fxsp = pm$fxsp,
          n_penalties = 0L
        )
        meta[[paste(part$linear_rows[[jj]], k, sep = ":")]] <- block
      }
    }
    
    if (length(part$nonlinear_rows)) {
      for (jj in seq_along(part$nonlinear_rows)) {
        block <- get_candidate_block(
          row = part$nonlinear_rows[[jj]], col = k,
          inclusion_template = inclusion, formula_base = formula_base,
          data = data, fam_obj = fam_obj, offset = offset, runtime = runtime,
          c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda
        )
        
        sc <- frame$x[[p_name]]$smooth.construct
        candidate_labels <- unique(c(block$label, block$raw_label))
        label <- candidate_labels[candidate_labels %in% names(sc)][1L]
        
        block <- refresh_block_from_frame(block, frame, label = label)
        obj <- frame$x[[p_name]]$smooth.construct[[block$label]]
        
        obj$prior <- make_candidate_prior(
          block, c_linear = c_linear, c0 = c0, cplus = cplus,
          a_lambda = a_lambda, b_lambda = b_lambda
        )
        
        frame$x[[p_name]]$smooth.construct[[block$label]] <- obj
        meta[[paste(part$nonlinear_rows[[jj]], k, sep = ":")]] <- block
        candidate_term_keys <- c(candidate_term_keys, block$term_key)
      }
    }
  }
  
  frame$candidate_metadata <- meta
  frame$candidate_term_keys <- unique(candidate_term_keys)
  frame
}

# Build the blocked selected-model BAMLSS transformer
make_selected_block_model_transform <- function(
    candidate_inclusion,
    formula_base,
    data,
    offset,
    c_linear = 1.0,
    c0 = 2.5,
    cplus = 1.0,
    a_lambda = 1,
    b_lambda = 1) {
  
  runtime <- new_runtime_state(cache_frames = FALSE, cache_blocks = TRUE)
  
  function(bf, ...) {
    if (!inherits(bf, "bamlss.frame")) {
      stop(
        "Blocked selected-model transformer requires a bamlss.frame object.",
        call. = FALSE
      )
    }
    
    raw_x <- bf$x
    bf$x <- bamlss::bamlss.engine.setup(bf$x, propose = "iwlsC_gp")
    
    bf$raw_design <- raw_x
    bf$model_offset <- offset
    
    bf <- install_structural_priors(bf)
    bf <- install_selected_posterior_priors(
      frame = bf,
      inclusion = candidate_inclusion,
      formula_base = formula_base,
      data = data,
      fam_obj = bf$family,
      offset = offset,
      runtime = runtime,
      c_linear = c_linear,
      c0 = c0,
      cplus = cplus,
      a_lambda = a_lambda,
      b_lambda = b_lambda
    )
    
    bf
  }
}

# Project a selected-model start into the blocked posterior frame
project_start_to_blocked_posterior <- function(
    start,
    source_frame,
    target_frame,
    fam_obj,
    a_lambda = 1,
    b_lambda = 1,
    tolerance = 1e-8) {
  
  source <- flatten_parameters(start)
  
  if (!length(source) ||
      is.null(names(source)) ||
      any(!nzchar(names(source))) ||
      any(!is.finite(source))) {
    stop(
      "Source selected-model start is not a finite named parameter vector.",
      call. = FALSE
    )
  }
  
  target <- generic_initial_parameters(target_frame, fam_obj)
  target <- seed_candidate_reference_state(
    theta = target,
    frame = target_frame,
    a_lambda = a_lambda,
    b_lambda = b_lambda,
    beta_value = 0
  )
  
  # structural and nonlinear parameters keep their names
  target <- merge_parameters(target, source)
  
  source_meta <- source_frame$candidate_metadata
  target_meta <- target_frame$candidate_metadata
  
  if (is.null(source_meta) || is.null(target_meta)) {
    stop("Selected-model frames are missing candidate metadata.", call. = FALSE)
  }
  
  target_linear <- Filter(
    function(z) identical(z$component, "linear"), target_meta
  )
  
  for (target_block in target_linear) {
    hit <- which(
      vapply(
        source_meta,
        function(source_block) {
          identical(source_block$component, "linear") &&
            identical(source_block$parameter, target_block$parameter) &&
            identical(source_block$variable, target_block$variable)
        },
        logical(1)
      )
    )
    
    if (length(hit) != 1L) {
      stop(
        "Could not identify exactly one source linear block for ",
        target_block$parameter, "::", target_block$variable, ".", call. = FALSE
      )
    }
    
    source_block <- source_meta[[hit]]
    source_columns <- colnames(source_block$X)
    target_columns <- colnames(target_block$X)
    column_map <- match(target_columns, source_columns)
    
    if (anyNA(column_map)) {
      stop(
        "Blocked start projection could not align design columns for ",
        target_block$parameter, "::", target_block$variable, ".", call. = FALSE
      )
    }
    
    source_X <- source_block$X[, column_map, drop = FALSE]
    same_design <- isTRUE(
      all.equal(
        unname(as.matrix(source_X)),
        unname(as.matrix(target_block$X)),
        tolerance = tolerance,
        check.attributes = FALSE
      )
    )
    
    if (!same_design) {
      stop(
        "Blocking changed the processed design matrix for ",
        target_block$parameter, "::", target_block$variable,
        ". The posterior target will not be changed silently.", call. = FALSE
      )
    }
    
    source_beta_names <- source_block$beta_names[column_map]
    target_beta_names <- target_block$beta_names
    
    if (length(source_beta_names) != length(target_beta_names) ||
        any(!source_beta_names %in% names(source)) ||
        any(!target_beta_names %in% names(target))) {
      stop(
        "Incompatible coefficient ownership in the blocked start for ",
        target_block$parameter, "::", target_block$variable, ".", call. = FALSE
      )
    }
    
    target[target_beta_names] <- source[source_beta_names]
  }
  
  if (any(!is.finite(target))) {
    stop("Blocked posterior start contains non-finite values.", call. = FALSE)
  }
  
  source_state <- merge_parameters(parameter_template(source_frame), source)
  source_eta <- state_predictors(
    theta = source_state, frame = source_frame,
    offset = source_frame$model_offset
  )
  target_eta <- state_predictors(
    theta = target, frame = target_frame, offset = target_frame$model_offset
  )
  
  if (!identical(names(source_eta), names(target_eta))) {
    stop(
      "Source and blocked posterior predictors have different parameter names.",
      call. = FALSE
    )
  }
  
  eta_difference <- vapply(
    names(source_eta),
    function(p_name) max(
      abs(as.numeric(source_eta[[p_name]]) - as.numeric(target_eta[[p_name]]))
    ),
    numeric(1)
  )
  
  if (any(!is.finite(eta_difference)) || any(eta_difference > tolerance)) {
    stop(
      "Blocked posterior reparameterization changed the linear predictor. ",
      "Maximum absolute differences: ",
      paste(
        sprintf("%s=%.6g", names(eta_difference), eta_difference),
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }
  
  source_prior <- term_prior_components(
    theta = source_state, frame = source_frame
  )
  target_prior <- term_prior_components(theta = target, frame = target_frame)
  
  if (!isTRUE(source_prior$finite) || !isTRUE(target_prior$finite)) {
    stop(
      "Prior-equivalence validation produced a non-finite prior.", call. = FALSE
    )
  }
  
  prior_difference <- c(
    structural = target_prior$structural_prior - source_prior$structural_prior,
    candidate = target_prior$candidate_prior - source_prior$candidate_prior
  )
  
  if (any(abs(prior_difference) > tolerance)) {
    stop(
      "Blocked posterior reparameterization changed the prior density. ",
      "Differences: ",
      paste(
        sprintf("%s=%.6g", names(prior_difference), prior_difference),
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }
  
  y <- extract_response(source_frame)
  source_loglik <- as.numeric(
    fam_obj$loglik(y, fam_obj$map2par(source_eta))
  )[1L]
  target_loglik <- as.numeric(
    fam_obj$loglik(y, fam_obj$map2par(target_eta))
  )[1L]
  
  if (!is.finite(source_loglik) ||
      !is.finite(target_loglik) ||
      abs(source_loglik - target_loglik) > tolerance) {
    stop(
      "Blocked posterior reparameterization changed the likelihood: source=",
      source_loglik, ", blocked=", target_loglik, ".", call. = FALSE
    )
  }
  
  attr(target, "blocking_validation") <- list(
    maximum_eta_difference = max(eta_difference),
    structural_prior_difference = prior_difference[["structural"]],
    candidate_prior_difference = prior_difference[["candidate"]],
    loglik_difference = target_loglik - source_loglik
  )
  
  target
}

# Runtime cache
new_runtime_state <- function(
    cache_frames = TRUE,
    cache_blocks = TRUE,
    max_frame_cache = Inf) {
  e <- new.env(parent = emptyenv())
  e$frame_cache <- if (isTRUE(cache_frames)) {
    new.env(hash = TRUE, parent = emptyenv())
  } else {
    NULL
  }
  e$block_cache <- if (isTRUE(cache_blocks)) {
    new.env(hash = TRUE, parent = emptyenv())
  } else {
    NULL
  }
  
  max_frame_cache <- as.numeric(max_frame_cache)[1L]
  if (!is.finite(max_frame_cache) && !is.infinite(max_frame_cache)) {
    stop("`max_frame_cache` must be a positive number or Inf.", call. = FALSE)
  }
  if (is.finite(max_frame_cache) && max_frame_cache < 1) {
    stop(
      "`max_frame_cache` must be at least 1 when frame caching is enabled.",
      call. = FALSE
    )
  }
  e$max_frame_cache <- if (is.null(e$frame_cache)) 0 else max_frame_cache
  e$frame_cache_order <- character(0)
  e$frame_cache_evictions <- 0L
  e$frame_cache_peak <- 0L
  
  e$frame_hits <- 0L
  e$frame_misses <- 0L
  e$block_hits <- 0L
  e$block_misses <- 0L
  e$initial_optimizer_fits <- 0L
  e$initial_optimizer_warning_count <- 0L
  e$initial_optimizer_warning_messages <- integer(0)
  e$within_attempts <- 0L
  e$within_successes <- 0L
  e$within_failures <- 0L
  e$within_error_messages <- integer(0)
  
  e$within_linear_ess_attempts <- numeric(0)
  e$within_linear_ess_successes <- numeric(0)
  e$within_linear_ess_failures <- numeric(0)
  e$within_linear_ess_moves <- numeric(0)
  e$within_linear_ess_coefficients_moved <- numeric(0)
  e$within_linear_ess_max_abs_change <- numeric(0)
  e$within_linear_ess_contractions <- numeric(0)
  e$within_linear_ess_max_contractions <- numeric(0)
  e$within_linear_ess_error_messages <- integer(0)
  
  e$proposal_failures <- 0L
  e$proposal_error_messages <- integer(0)
  e$lambda_laplace_draws <- 0L
  e$lambda_prior_draws <- 0L
  e$lambda_laplace_failures <- 0L
  e$beta_local_draws <- 0L
  e$beta_prior_draws <- 0L
  e$predictor_scaled_local_draws <- 0L
  e$predictor_scaled_prior_draws <- 0L
  e$predictor_scale_sum <- 0
  e$predictor_scale_min <- 1
  e
}

# Increment a named runtime statistic
increment_runtime_stat <- function(runtime, field, key, value = 1) {
  key <- as.character(key)[1L]
  value <- as.numeric(value)[1L]
  
  if (!nzchar(key)) {
    stop("Runtime statistic key must be non-empty.", call. = FALSE)
  }
  if (!is.finite(value)) {
    stop("Runtime statistic increment must be finite.", call. = FALSE)
  }
  
  x <- runtime[[field]]
  if (is.null(x)) x <- numeric(0)
  
  if (key %in% names(x)) {
    x[[key]] <- x[[key]] + value
  } else {
    x <- c(x, stats::setNames(value, key))
  }
  
  runtime[[field]] <- x
  invisible(NULL)
}

# Update a named runtime maximum
update_runtime_maximum <- function(runtime, field, key, value) {
  key <- as.character(key)[1L]
  value <- as.numeric(value)[1L]
  
  if (!nzchar(key)) {
    stop("Runtime maximum key must be non-empty.", call. = FALSE)
  }
  if (!is.finite(value)) return(invisible(NULL))
  
  x <- runtime[[field]]
  if (is.null(x)) x <- numeric(0)
  
  if (key %in% names(x)) {
    x[[key]] <- max(x[[key]], value)
  } else {
    x <- c(x, stats::setNames(value, key))
  }
  
  runtime[[field]] <- x
  invisible(NULL)
}

# Ratio of named runtime statistics
runtime_stat_ratio <- function(numerator, denominator) {
  keys <- union(names(numerator), names(denominator))
  if (!length(keys)) return(numeric(0))
  
  num <- stats::setNames(rep(0, length(keys)), keys)
  den <- stats::setNames(rep(0, length(keys)), keys)
  
  num[names(numerator)] <- as.numeric(numerator)
  den[names(denominator)] <- as.numeric(denominator)
  
  out <- num / den
  out[den <= 0] <- NA_real_
  out
}

# Active linear candidate blocks
active_linear_candidate_blocks <- function(frame, parameter) {
  meta <- frame$candidate_metadata
  if (is.null(meta) || !length(meta)) return(list())
  
  keep <- vapply(
    meta,
    function(block) {
      identical(block$component, "linear") &&
        identical(block$parameter, parameter) &&
        !isTRUE(block$smooth)
    },
    logical(1)
  )
  
  blocks <- meta[keep]
  if (length(blocks) > 1L) {
    rows <- vapply(
      blocks, function(block) as.integer(block$row)[1L], integer(1)
    )
    blocks <- blocks[order(rows)]
  }
  
  blocks
}

# Elliptical slice update for active linear candidates
update_linear_candidates_ess <- function(
    state,
    frame,
    fam_obj,
    runtime,
    parameter,
    a_V = 0.5,
    a_A = 0.2) {
  
  blocks <- active_linear_candidate_blocks(frame = frame, parameter = parameter)
  if (!length(blocks)) return(state)
  
  increment_runtime_stat(runtime, "within_linear_ess_attempts", parameter)
  
  beta_names <- unlist(
    lapply(blocks, function(block) block$beta_names), use.names = FALSE
  )
  if (!length(beta_names) || any(!nzchar(beta_names))) {
    stop(
      "Active linear candidate blocks contain invalid coefficient names for ",
      parameter, ".", call. = FALSE
    )
  }
  if (anyDuplicated(beta_names)) {
    stop(
      "Duplicated coefficient names in the active linear blocks for ",
      parameter, ".", call. = FALSE
    )
  }
  
  design <- do.call(cbind, lapply(blocks, function(block) {
    X <- as.matrix(block$X)
    if (ncol(X) != length(block$beta_names)) {
      stop(
        "Candidate design and coefficient dimensions disagree for ",
        block$candidate_key, " -> ", parameter, ".", call. = FALSE
      )
    }
    X
  }))
  
  if (ncol(design) != length(beta_names)) {
    stop(
      "Joint linear candidate design has an invalid coefficient dimension for ",
      parameter, ".", call. = FALSE
    )
  }
  
  prior_sd <- vapply(
    blocks, function(block) as.numeric(block$c_linear)[1L], numeric(1)
  )
  if (any(!is.finite(prior_sd)) || any(prior_sd <= 0)) {
    stop(
      "Active linear candidate blocks contain an invalid prior scale for ",
      parameter, ".", call. = FALSE
    )
  }
  
  prior_tolerance <- sqrt(.Machine$double.eps) * max(1, abs(prior_sd[[1L]]))
  if (any(abs(prior_sd - prior_sd[[1L]]) > prior_tolerance)) {
    stop(
      "Joint elliptical-slice updates need one prior scale within ", parameter,
      ".", call. = FALSE
    )
  }
  prior_sd <- prior_sd[[1L]]
  
  theta <- merge_parameters(parameter_template(frame), state$theta)
  missing_beta <- setdiff(beta_names, names(theta))
  if (length(missing_beta)) {
    stop(
      "Current RJMCMC state is missing active linear coefficient(s): ",
      paste(missing_beta, collapse = ", "), ".", call. = FALSE
    )
  }
  
  beta_current <- as.numeric(theta[beta_names])
  if (any(!is.finite(beta_current))) {
    stop(
      "Current linear candidate coefficients are non-finite for ", parameter,
      ".", call. = FALSE
    )
  }
  
  y <- extract_response(frame)
  if (nrow(design) != length(y)) {
    stop(
      "Joint linear candidate design has incompatible row dimension for ",
      parameter, ".", call. = FALSE
    )
  }
  
  eta_current <- state_predictors(
    theta = theta, frame = frame, offset = frame$model_offset
  )
  if (!parameter %in% names(eta_current)) {
    stop(
      "Distributional predictor is unavailable for ", parameter, ".",
      call. = FALSE
    )
  }
  
  current_contribution <- as.numeric(design %*% beta_current)
  eta_base <- eta_current
  eta_base[[parameter]] <- eta_current[[parameter]] - current_contribution
  
  evaluate_loglik <- function(beta) {
    contribution <- as.numeric(design %*% beta)
    eta_candidate <- eta_base
    eta_candidate[[parameter]] <- eta_base[[parameter]] + contribution
    par_candidate <- tryCatch(
      fam_obj$map2par(eta_candidate), error = function(e) NULL
    )
    if (is.null(par_candidate)) return(-Inf)
    loglik <- tryCatch(
      as.numeric(fam_obj$loglik(y, par_candidate))[1L], error = function(e) -Inf
    )
    if (!is.finite(loglik)) return(-Inf)
    loglik
  }
  
  current_loglik <- evaluate_loglik(beta_current)
  if (!is.finite(current_loglik)) {
    stop(
      "Current RJMCMC state has non-finite conditional log-likelihood for ",
      parameter, ".", call. = FALSE
    )
  }
  
  auxiliary <- stats::rnorm(length(beta_current), mean = 0, sd = prior_sd)
  log_slice <- current_loglik + log(stats::runif(1L))
  angle <- stats::runif(1L, min = 0, max = 2 * pi)
  lower <- angle - 2 * pi
  upper <- angle
  contractions <- 0L
  
  repeat {
    beta_proposed <- beta_current * cos(angle) + auxiliary * sin(angle)
    proposed_loglik <- evaluate_loglik(beta_proposed)
    
    if (is.finite(proposed_loglik) && proposed_loglik >= log_slice) break
    
    contractions <- contractions + 1L
    if (angle < 0) lower <- angle else upper <- angle
    
    width <- upper - lower
    collapse_tolerance <- .Machine$double.eps * max(1, abs(lower), abs(upper))
    
    if (!is.finite(width) || width <= collapse_tolerance) {
      increment_runtime_stat(runtime, "within_linear_ess_failures", parameter)
      record_runtime_message(
        runtime, "within_linear_ess_error_messages",
        paste0(parameter, ": elliptical-slice bracket collapsed numerically")
      )
      return(state)
    }
    angle <- stats::runif(1L, min = lower, max = upper)
  }
  
  theta_proposed <- theta
  theta_proposed[beta_names] <- beta_proposed
  proposed_state <- state
  proposed_state$theta <- theta_proposed
  
  proposed_state$logPost <- target_log_posterior(
    proposed_state, frame, fam_obj,
    offset = frame$model_offset,
    a_V = a_V, a_A = a_A
  )
  
  if (!is.finite(proposed_state$logPost)) {
    increment_runtime_stat(runtime, "within_linear_ess_failures", parameter)
    record_runtime_message(
      runtime, "within_linear_ess_error_messages",
      paste0(parameter, ": elliptical-slice proposal outside the support")
    )
    return(state)
  }
  
  delta <- beta_proposed - beta_current
  movement_tolerance <-
    sqrt(.Machine$double.eps) * pmax(1, abs(beta_current), abs(beta_proposed))
  moved <- abs(delta) > movement_tolerance
  
  increment_runtime_stat(runtime, "within_linear_ess_successes", parameter)
  increment_runtime_stat(
    runtime, "within_linear_ess_contractions", parameter, contractions
  )
  update_runtime_maximum(
    runtime, "within_linear_ess_max_contractions", parameter, contractions
  )
  
  if (any(moved)) {
    increment_runtime_stat(runtime, "within_linear_ess_moves", parameter)
    increment_runtime_stat(
      runtime, "within_linear_ess_coefficients_moved", parameter, sum(moved)
    )
    update_runtime_maximum(
      runtime, "within_linear_ess_max_abs_change", parameter, max(abs(delta))
    )
  }
  
  proposed_state
}

# Record a counted runtime message
record_runtime_message <- function(runtime, field, message_text) {
  msg <- trimws(as.character(message_text)[1L])
  if (!nzchar(msg)) msg <- "<empty message>"
  
  x <- runtime[[field]]
  if (is.null(x)) x <- integer(0)
  if (msg %in% names(x)) {
    x[[msg]] <- x[[msg]] + 1L
  } else {
    x <- c(x, stats::setNames(1L, msg))
  }
  runtime[[field]] <- x
  invisible(NULL)
}

# Record a trans-dimensional proposal failure
record_proposal_failure <- function(runtime, message_text) {
  runtime$proposal_failures <- runtime$proposal_failures + 1L
  record_runtime_message(runtime, "proposal_error_messages", message_text)
  invisible(NULL)
}

# Compact runtime message table
runtime_message_table <- function(x, max_rows = 10L) {
  if (is.null(x) || !length(x)) {
    return(data.frame(
      message = character(0),
      count = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  ord <- order(as.integer(x), decreasing = TRUE)
  ord <- head(ord, as.integer(max_rows))
  data.frame(
    message = names(x)[ord],
    count = as.integer(x[ord]),
    stringsAsFactors = FALSE
  )
}

# Build an undecorated BAMLSS frame
build_raw_frame <- function(
    inclusion, formula_base, data, fam_obj, offset,
    c_linear = 1.5, a_lambda = 1, b_lambda = 1) {
  
  form <- build_distributional_formula(
    inclusion, formula_base, data, fam_obj,
    c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda
  )
  
  parsed_form <- parse_distributional_formula(form, fam_obj)
  bf <- bamlss::bamlss.frame(
    parsed_form, data = data, family = fam_obj, offset = offset
  )
  
  raw_x <- bf$x
  bf$x <- bamlss::bamlss.engine.setup(bf$x, propose = "iwlsC_gp")
  bf$raw_design <- raw_x
  bf$model_formula <- form
  bf$model_state_key <- state_key(inclusion)
  bf$model_offset <- if (!is.null(bf$offset)) bf$offset else offset
  bf
}

# Exact proper prior for a dynamic candidate block
log_candidate_prior_values <- function(
    beta, tau2, block,
    c_linear = 1.5, c0 = 316.227766, cplus = 316.227766,
    a_lambda = 1, b_lambda = 1,
    floor = 1e-12) {
  
  beta <- as.numeric(beta)
  d <- ncol(block$X)
  
  if (length(beta) != d || any(!is.finite(beta))) return(-Inf)
  
  if (!is.finite(c_linear) || c_linear <= 0 ||
      !is.finite(c0) || c0 <= 0 || !is.finite(cplus) || cplus <= 0 ||
      !is.finite(a_lambda) || a_lambda <= 0 ||
      !is.finite(b_lambda) || b_lambda <= 0) {
    return(-Inf)
  }
  
  if (!isTRUE(block$smooth) || !length(block$S)) {
    return(sum(stats::dnorm(beta, mean = 0, sd = c_linear, log = TRUE)))
  }
  
  S <- lapply(block$S, function(z) symmetrize_matrix(as.matrix(z)))
  if (any(vapply(S, function(z) nrow(z) != d || ncol(z) != d, logical(1)))) {
    return(-Inf)
  }
  
  tau2 <- as.numeric(tau2)
  if (!length(tau2) || any(!is.finite(tau2)) || any(tau2 <= 0)) return(-Inf)
  if (length(tau2) == 1L && length(S) > 1L) tau2 <- rep(tau2, length(S))
  if (length(tau2) != length(S)) return(-Inf)
  lambda <- 1 / tau2
  
  Ksum <- Reduce(`+`, S)
  ee <- eigen(symmetrize_matrix(Ksum), symmetric = TRUE)
  tol <- max(floor, max(abs(ee$values), 1) * 1e-10)
  pos <- which(ee$values > tol)
  nul <- which(ee$values <= tol)
  
  Q <- matrix(0, d, d)
  if (length(nul)) {
    P0 <- tcrossprod(ee$vectors[, nul, drop = FALSE])
    Q <- Q + P0 / (c0^2)
  }
  if (length(pos)) {
    Pp <- tcrossprod(ee$vectors[, pos, drop = FALSE])
    Q <- Q + Pp / (cplus^2)
  }
  for (j in seq_along(S)) Q <- Q + lambda[j] * S[[j]]
  Q <- symmetrize_matrix(Q)
  
  ev <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  qtol <- max(floor, max(abs(ev), 1) * 1e-12)
  if (any(!is.finite(ev)) || any(ev < -qtol)) return(-Inf)
  
  ev <- pmax(ev, floor)
  logdetQ <- sum(log(ev))
  quad <- as.numeric(crossprod(beta, Q %*% beta))
  if (!is.finite(logdetQ) || !is.finite(quad)) return(-Inf)
  
  log_beta <- 0.5 * logdetQ - 0.5 * d * log(2 * pi) - 0.5 * quad
  log_tau2 <- sum(
    stats::dgamma(lambda, shape = a_lambda, rate = b_lambda, log = TRUE) -
      2 * log(tau2)
  )
  as.numeric(log_beta + log_tau2)
}

# Candidate prior function for a BAMLSS model term
make_candidate_prior <- function(
    block,
    c_linear = 1.5, c0 = 316.227766, cplus = 316.227766,
    a_lambda = 1, b_lambda = 1) {
  force(block)
  force(c_linear)
  force(c0)
  force(cplus)
  force(a_lambda)
  force(b_lambda)
  
  block_local <- block
  c_linear_local <- as.numeric(c_linear)[1L]
  c0_local <- as.numeric(c0)[1L]
  cplus_local <- as.numeric(cplus)[1L]
  a_lambda_local <- as.numeric(a_lambda)[1L]
  b_lambda_local <- as.numeric(b_lambda)[1L]
  
  local_get_par <- function(x, name) {
    idx <- grep(paste0("^", name), names(x))
    if (length(idx)) x[idx] else numeric(0)
  }
  
  f <- function(parameters) {
    beta <- tryCatch(
      as.numeric(local_get_par(parameters, "b")),
      error = function(e) numeric(0)
    )
    tau2 <- tryCatch(
      as.numeric(local_get_par(parameters, "tau2")),
      error = function(e) numeric(0)
    )
    log_candidate_prior_values(
      beta = beta, tau2 = tau2, block = block_local,
      c_linear = c_linear_local, c0 = c0_local, cplus = cplus_local,
      a_lambda = a_lambda_local, b_lambda = b_lambda_local
    )
  }
  
  attr(f, "var_prior") <- "proper_candidate_prior"
  attr(f, "candidate_term_key") <- block_local$term_key
  attr(f, "candidate_key") <- block_local$candidate_key
  attr(f, "candidate_component") <- block_local$component
  f
}

# Resolve term-local BAMLSS parameter ownership
term_parameter_map <- function(frame, p_name, label) {
  sc <- frame$x[[p_name]]$smooth.construct
  if (is.null(sc) || is.null(sc[[label]])) {
    stop(
      "Cannot resolve term-local parameters: processed term '", label,
      "' is absent from parameter ", p_name, ".",
      call. = FALSE
    )
  }
  obj <- sc[[label]]
  
  plist <- tryCatch(
    bamlss::parameters(frame, list = TRUE),
    error = function(e) {
      stop(
        "BAMLSS term-local parameter extraction failed: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  term_par <- NULL
  if (!is.null(plist[[p_name]]) && !is.null(plist[[p_name]]$s)) {
    term_par <- plist[[p_name]]$s[[label]]
  }
  
  if (is.null(term_par)) {
    available <- if (!is.null(plist[[p_name]]) && !is.null(plist[[p_name]]$s)) {
      paste(names(plist[[p_name]]$s), collapse = ", ")
    } else {
      "<none>"
    }
    stop(
      "BAMLSS parameters(list=TRUE) has no smooth-term parameter block for ",
      p_name, "::", label, ". Available smooth terms: ", available, ".",
      call. = FALSE
    )
  }
  
  term_par <- unlist(term_par, recursive = TRUE, use.names = TRUE)
  local_names <- names(term_par)
  if (is.null(local_names) ||
      any(!nzchar(local_names)) ||
      anyDuplicated(local_names)) {
    stop(
      "Term-local BAMLSS parameter names are missing or duplicated for ",
      p_name, "::", label, ".",
      call. = FALSE
    )
  }
  
  flat <- parameter_template(frame)
  flat_names <- names(flat)
  prefix <- paste0(p_name, ".s.", label, ".")
  direct <- paste0(prefix, local_names)
  
  if (all(direct %in% flat_names) && !anyDuplicated(direct)) {
    global_names <- direct
  } else {
    owned <- flat_names[startsWith(flat_names, prefix)]
    suffix <- substring(owned, nchar(prefix) + 1L)
    idx <- match(local_names, suffix)
    
    if (length(owned) != length(local_names) ||
        anyNA(idx) ||
        anyDuplicated(idx)) {
      stop(
        "Could not map term-local parameters of ", p_name, "::", label,
        " to the flat BAMLSS state without ambiguity. Local names: ",
        paste(local_names, collapse = ", "), "; flat prefix matches: ",
        if (length(owned)) paste(owned, collapse = ", ") else "<none>", ".",
        call. = FALSE
      )
    }
    global_names <- owned[idx]
  }
  
  tau_idx <- startsWith(local_names, "tau2")
  beta_idx <- !tau_idx
  
  if (sum(beta_idx) != ncol(obj$X)) {
    stop(
      "Term-local coefficient count does not match the processed design for ",
      p_name, "::", label, ": beta=", sum(beta_idx),
      ", ncol(X)=", ncol(obj$X),
      "; local parameters=", paste(local_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  
  list(
    local_names = local_names,
    global_names = global_names,
    beta_names = global_names[beta_idx],
    tau_names = global_names[tau_idx],
    fixed = isTRUE(obj$fixed),
    fxsp = isTRUE(obj$fxsp),
    n_penalties = if (is.null(obj$S)) 0L else length(obj$S)
  )
}

# Refresh a candidate block from the processed frame
refresh_block_from_frame <- function(block, frame, label = NULL) {
  p_name <- block$parameter
  sc <- frame$x[[p_name]]$smooth.construct
  if (is.null(sc)) {
    stop(
      "Processed BAMLSS frame has no smooth.construct list for parameter ",
      p_name, ".", call. = FALSE
    )
  }
  
  candidates <- unique(c(
    if (!is.null(label)) as.character(label) else character(0),
    if (!is.null(block$label)) as.character(block$label) else character(0),
    if (!is.null(block$raw_label)) {
      as.character(block$raw_label)
    } else {
      character(0)
    }
  ))
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[candidates %in% names(sc)]
  
  if (!length(hit)) {
    stop(
      "Candidate term ", block$candidate_key, " -> ", p_name,
      " is absent from the processed BAMLSS frame. Expected one of: ",
      if (length(candidates)) paste(candidates, collapse = ", ") else "<none>",
      "; available: ",
      if (length(names(sc))) paste(names(sc), collapse = ", ") else "<none>",
      ".",
      call. = FALSE
    )
  }
  
  label <- hit[1L]
  obj <- sc[[label]]
  
  X <- as.matrix(obj$X)
  y <- frame$y
  n_expected <- if (is.matrix(y) || is.data.frame(y)) nrow(y) else length(y)
  
  if (nrow(X) != n_expected && !is.null(obj$binning$match.index)) {
    X <- X[obj$binning$match.index, , drop = FALSE]
  }
  if (nrow(X) != n_expected) {
    stop(
      "Candidate design block has incompatible row dimension for ",
      p_name, "::", label, ": nrow(X)=", nrow(X),
      ", response rows=", n_expected, ".",
      call. = FALSE
    )
  }
  if (!ncol(X) || any(!is.finite(X))) {
    stop(
      "Candidate processed design is empty or non-finite for ",
      p_name, "::", label, ".", call. = FALSE
    )
  }
  
  S <- obj$S
  if (is.null(S)) S <- list()
  S <- lapply(S, function(z) {
    if (is.function(z)) {
      np <- attr(z, "npar")
      if (is.null(np)) np <- ncol(X)
      z(c("b" = rep(0, np)))
    } else {
      z
    }
  })
  S <- lapply(S, as.matrix)
  
  pm <- term_parameter_map(frame, p_name, label)
  
  block$label <- label
  block$term_key <- paste(p_name, label, sep = "::")
  block$X <- X
  block$design_colnames <- colnames(X)
  
  if (is.null(block$design_colnames) ||
      length(block$design_colnames) != ncol(X)) {
    block$design_colnames <- paste0("basis_", seq_len(ncol(X)))
  }
  
  block$S <- S
  block$parameter_names <- pm$global_names
  block$beta_names <- pm$beta_names
  block$tau_names <- pm$tau_names
  block$local_parameter_names <- pm$local_names
  block$fixed <- pm$fixed
  block$fxsp <- pm$fxsp
  block$n_penalties <- pm$n_penalties
  
  if (length(block$beta_names) != ncol(block$X)) {
    stop(
      "Candidate term ownership is inconsistent for ", block$term_key,
      ": beta names=", length(block$beta_names),
      ", ncol(X)=", ncol(block$X), ".",
      call. = FALSE
    )
  }
  
  if (!isTRUE(block$smooth)) {
    if (length(block$tau_names) || length(block$S)) {
      stop(
        "Unpenalized linear candidate unexpectedly owns smoothing parameters: ",
        block$term_key,
        " | local=", paste(block$local_parameter_names, collapse = ", "),
        " | tau_n=", length(block$tau_names),
        " | penalties=", length(block$S),
        " | fixed=", block$fixed, " | fxsp=", block$fxsp, ".",
        call. = FALSE
      )
    }
  } else {
    if (!length(block$S)) {
      stop(
        "Smooth candidate has no penalty matrix: ", block$term_key,
        call. = FALSE
      )
    }
    if (!length(block$tau_names) &&
        !isTRUE(block$fixed) &&
        !isTRUE(block$fxsp)) {
      stop(
        "Smooth candidate has no term-local tau2 parameter: ", block$term_key,
        " | local=", paste(block$local_parameter_names, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }
  
  block
}

# Extract a proper candidate design block
get_candidate_block <- function(
    row, col, inclusion_template, formula_base, data, fam_obj, offset, runtime,
    c_linear = 2.0, a_lambda = 1, b_lambda = 1) {
  
  key <- paste(row, col, sep = ":")
  if (!is.null(runtime$block_cache) &&
      exists(key, runtime$block_cache, inherits = FALSE)) {
    runtime$block_hits <- runtime$block_hits + 1L
    return(get(key, runtime$block_cache, inherits = FALSE))
  }
  runtime$block_misses <- runtime$block_misses + 1L
  
  zero <- inclusion_template
  zero[,] <- FALSE
  one <- zero
  one[row, col] <- TRUE
  
  bf0 <- build_raw_frame(
    zero, formula_base, data, fam_obj, offset,
    c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda
  )
  bf1 <- build_raw_frame(
    one, formula_base, data, fam_obj, offset,
    c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda
  )
  
  p_name <- colnames(inclusion_template)[col]
  
  sc0_raw <- bf0$raw_design[[p_name]]$smooth.construct
  sc1_raw <- bf1$raw_design[[p_name]]$smooth.construct
  nr0 <- if (is.null(sc0_raw)) character(0) else names(sc0_raw)
  nr1 <- if (is.null(sc1_raw)) character(0) else names(sc1_raw)
  new_raw <- setdiff(nr1, nr0)
  
  if (length(new_raw) != 1L) {
    stop(
      "Could not identify a unique raw candidate block for ",
      rownames(inclusion_template)[row], " -> ", p_name,
      ". New raw terms: ", paste(new_raw, collapse = ", "),
      call. = FALSE
    )
  }
  
  sc0_proc <- bf0$x[[p_name]]$smooth.construct
  sc1_proc <- bf1$x[[p_name]]$smooth.construct
  np0 <- if (is.null(sc0_proc)) character(0) else names(sc0_proc)
  np1 <- if (is.null(sc1_proc)) character(0) else names(sc1_proc)
  new_proc <- setdiff(np1, np0)
  
  if (length(new_proc) != 1L) {
    stop(
      "Could not identify a unique processed candidate block for ",
      rownames(inclusion_template)[row], " -> ", p_name,
      ". New processed terms: ", paste(new_proc, collapse = ", "),
      call. = FALSE
    )
  }
  
  candidate_key <- rownames(inclusion_template)[row]
  spec <- candidate_specification(candidate_key, data)
  block <- list(
    row = row, col = col, variable = spec$variable,
    candidate_key = candidate_key,
    component = spec$component, parameter = p_name,
    raw_label = new_raw, label = new_proc,
    term_key = paste(p_name, new_proc, sep = "::"),
    smooth = spec$component %in% c("nonlinear", "legacy_smooth"),
    c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda
  )
  
  block <- refresh_block_from_frame(block, bf1, label = new_proc)
  
  if (!is.null(runtime$block_cache)) {
    assign(key, block, envir = runtime$block_cache)
  }
  block
}

# Install candidate priors in active terms
# bamlss::GMCMC_slice ignores the offset; it is added here
slice_with_offset <- function(
    family, theta, id, eta, y, data, weights = NULL, offset = NULL, ...) {
  if (!is.null(offset)) {
    for (j in names(offset)) eta[[j]] <- eta[[j]] + offset[[j]]
  }
  slice <- get("GMCMC_slice", envir = asNamespace("bamlss"))
  slice(family, theta, id, eta, y, data, ...)
}

# N(0, sd^2) prior for every coefficient of a structural parametric term
structural_gaussian_prior <- function(sd) {
  force(sd)
  f <- function(parameters) {
    nm <- names(parameters)
    # same selection as bamlss get.par(., "b"); names are the design columns
    b <- if (is.null(nm)) parameters else parameters[
      !grepl("tau2|edf|lasso|bw", nm)
    ]
    b <- suppressWarnings(as.numeric(b))
    if (!length(b) || any(!is.finite(b))) return(-Inf)
    sum(stats::dnorm(b, mean = 0, sd = sd, log = TRUE))
  }
  attr(f, "var_prior") <- "structural_gaussian"
  f
}

# Intercepts and year effects get the prior set by
# options(rbst.structural_sd = ...) and slice updates; unset keeps the bamlss
# defaults
install_structural_priors <- function(frame) {
  sd <- getOption("rbst.structural_sd")
  if (is.null(sd)) return(frame)
  sd <- as.numeric(sd)[1L]
  if (!is.finite(sd) || sd <= 0) {
    stop("Option `rbst.structural_sd` must be a positive number.",
         call. = FALSE)
  }
  for (p_name in names(frame$x)) {
    obj <- frame$x[[p_name]]$smooth.construct$model.matrix
    if (is.null(obj)) next
    obj$prior <- structural_gaussian_prior(sd)
    obj$propose <- slice_with_offset
    frame$x[[p_name]]$smooth.construct$model.matrix <- obj
  }
  frame
}

install_candidate_priors <- function(
    frame, inclusion, formula_base, data, fam_obj, offset, runtime,
    c_linear = 2.0, c0 = 316.227766, cplus = 316.227766,
    a_lambda = 1, b_lambda = 1) {
  
  frame <- install_structural_priors(frame)
  active <- which(inclusion, arr.ind = TRUE)
  meta <- list()
  if (!nrow(active)) {
    frame$candidate_metadata <- meta
    frame$candidate_term_keys <- character(0)
    return(frame)
  }
  
  for (ii in seq_len(nrow(active))) {
    r <- active[ii, 1L]
    k <- active[ii, 2L]
    block <- get_candidate_block(
      r, k, inclusion, formula_base, data, fam_obj, offset, runtime,
      c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda
    )
    
    p_name <- block$parameter
    sc <- frame$x[[p_name]]$smooth.construct
    candidate_labels <- unique(c(block$label, block$raw_label))
    label <- candidate_labels[candidate_labels %in% names(sc)][1L]
    
    if (length(label) != 1L || is.na(label) || !nzchar(label)) {
      stop(
        "Active candidate term not found in processed BAMLSS frame for ",
        block$candidate_key, " -> ", p_name,
        ". Expected one of: ", paste(candidate_labels, collapse = ", "),
        ". Available terms: ",
        if (is.null(sc)) "<none>" else paste(names(sc), collapse = ", "),
        call. = FALSE
      )
    }
    
    block <- refresh_block_from_frame(block, frame, label = label)
    obj <- frame$x[[p_name]]$smooth.construct[[block$label]]
    
    obj$prior <- make_candidate_prior(
      block,
      c_linear = c_linear, c0 = c0, cplus = cplus,
      a_lambda = a_lambda, b_lambda = b_lambda
    )
    
    prior_term_key <- attr(obj$prior, "candidate_term_key")
    if (!identical(prior_term_key, block$term_key)) {
      stop(
        "Candidate-prior factory binding mismatch: expected ", block$term_key,
        ", got ", if (is.null(prior_term_key)) "<missing>" else prior_term_key,
        ".",
        call. = FALSE
      )
    }
    
    frame$x[[p_name]]$smooth.construct[[block$label]] <- obj
    meta[[paste(r, k, sep = ":")]] <- block
  }
  
  frame$candidate_metadata <- meta
  frame$candidate_term_keys <- unique(
    vapply(meta, `[[`, character(1), "term_key")
  )
  frame
}

# Obtain a cached BAMLSS frame
get_model_frame <- function(
    inclusion, formula_base, data, fam_obj, offset, runtime,
    c_linear = 2.0, a_lambda = 1, b_lambda = 1,
    c0 = 316.227766, cplus = 316.227766) {
  
  key <- state_key(inclusion)
  if (!is.null(runtime$frame_cache) &&
      exists(key, runtime$frame_cache, inherits = FALSE)) {
    runtime$frame_hits <- runtime$frame_hits + 1L
    runtime$frame_cache_order <- c(
      runtime$frame_cache_order[runtime$frame_cache_order != key],
      key
    )
    return(get(key, runtime$frame_cache, inherits = FALSE))
  }
  runtime$frame_misses <- runtime$frame_misses + 1L
  
  bf <- build_raw_frame(
    inclusion, formula_base, data, fam_obj, offset,
    c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda
  )
  
  bf <- install_candidate_priors(
    bf, inclusion, formula_base, data, fam_obj, offset, runtime,
    c_linear = c_linear, c0 = c0, cplus = cplus,
    a_lambda = a_lambda, b_lambda = b_lambda
  )
  
  if (!is.null(runtime$frame_cache)) {
    if (is.finite(runtime$max_frame_cache)) {
      while (length(runtime$frame_cache_order) >= runtime$max_frame_cache) {
        evict_key <- runtime$frame_cache_order[1L]
        if (exists(evict_key, runtime$frame_cache, inherits = FALSE)) {
          rm(list = evict_key, envir = runtime$frame_cache)
        }
        runtime$frame_cache_order <- runtime$frame_cache_order[-1L]
        runtime$frame_cache_evictions <- runtime$frame_cache_evictions + 1L
      }
    }
    assign(key, bf, envir = runtime$frame_cache)
    runtime$frame_cache_order <- c(runtime$frame_cache_order, key)
    runtime$frame_cache_peak <- max(
      runtime$frame_cache_peak,
      length(runtime$frame_cache_order)
    )
  }
  bf
}

# Extract a named BAMLSS parameter template
parameter_template <- function(frame) {
  p <- bamlss::parameters(frame, list = FALSE)
  if (is.list(p)) p <- unlist(p, recursive = TRUE, use.names = TRUE)
  nm <- names(p)
  p <- stats::setNames(as.numeric(p), nm)
  
  if (is.null(names(p)) || any(!nzchar(names(p)))) {
    stop(
      "BAMLSS parameter template does not contain stable parameter names.",
      call. = FALSE
    )
  }
  if (anyDuplicated(names(p))) {
    dup <- unique(names(p)[duplicated(names(p))])
    stop(
      "BAMLSS parameter template contains duplicated names: ",
      paste(dup, collapse = ", "),
      call. = FALSE
    )
  }
  p
}

# Coerce a BAMLSS parameter container to a named numeric vector
flatten_parameters <- function(x) {
  if (is.null(x) || !length(x)) return(numeric(0))
  
  if (is.atomic(x) && !is.list(x)) {
    nm <- names(x)
    z <- suppressWarnings(as.numeric(x))
    if (!is.null(nm)) names(z) <- nm
    return(z)
  }
  
  z <- tryCatch(
    unlist(x, recursive = TRUE, use.names = TRUE),
    error = function(e) NULL
  )
  if (is.null(z) || !length(z)) return(numeric(0))
  
  nm <- names(z)
  z <- suppressWarnings(as.numeric(z))
  if (!is.null(nm)) names(z) <- nm
  z
}

# Copy common parameters between nested models
merge_parameters <- function(template, source = NULL) {
  out <- flatten_parameters(template)
  if (is.null(names(out)) || any(!nzchar(names(out)))) {
    stop("Parameter template must be a named numeric vector.", call. = FALSE)
  }
  if (is.null(source) || !length(source)) return(out)
  
  src <- flatten_parameters(source)
  if (!length(src)) return(out)
  sn <- names(src)
  
  if (is.null(sn) || any(!nzchar(sn))) return(out)
  
  for (nm in names(out)) {
    j <- match(nm, sn)
    if (!is.na(j)) {
      out[[nm]] <- src[[j]]
      next
    }
    
    js <- which(endsWith(sn, nm))
    if (length(js) == 1L) out[[nm]] <- src[[js]]
  }
  
  storage.mode(out) <- "double"
  out
}

# Seed a finite reference state for active candidate terms
seed_candidate_reference_state <- function(
    theta, frame, a_lambda = 1, b_lambda = 1, beta_value = 0) {
  
  theta <- flatten_parameters(theta)
  if (!length(theta) || is.null(names(theta)) ||
      any(!nzchar(names(theta))) || any(!is.finite(theta))) {
    stop(
      "The reference state needs a finite named parameter vector.",
      call. = FALSE
    )
  }
  
  if (!is.finite(a_lambda) || a_lambda <= 0 ||
      !is.finite(b_lambda) || b_lambda <= 0) {
    stop(
      "The reference state needs positive `a_lambda` and `b_lambda`.",
      call. = FALSE
    )
  }
  
  if (!is.finite(beta_value)) {
    stop("Candidate reference-state beta value must be finite.", call. = FALSE)
  }
  
  meta <- frame$candidate_metadata
  if (is.null(meta) || !length(meta)) return(theta)
  
  tau2_ref <- b_lambda / a_lambda
  
  for (key in names(meta)) {
    block <- meta[[key]]
    beta_names <- block$beta_names
    tau_names <- block$tau_names
    
    if (length(beta_names) != ncol(block$X) ||
        any(!beta_names %in% names(theta))) {
      stop(
        "Candidate reference-state seed found inconsistent beta ownership for ",
        block$term_key, ".",
        call. = FALSE
      )
    }
    theta[beta_names] <- rep(beta_value, length(beta_names))
    
    if (isTRUE(block$smooth)) {
      if (length(block$S) && !length(tau_names)) {
        stop(
          "Penalized smooth without tau2 in the reference state for ",
          block$term_key, ".",
          call. = FALSE
        )
      }
      if (length(tau_names)) {
        if (any(!tau_names %in% names(theta))) {
          stop(
            "Candidate reference-state seed cannot locate tau2 parameters for ",
            block$term_key, ".",
            call. = FALSE
          )
        }
        theta[tau_names] <- rep(tau2_ref, length(tau_names))
      }
    } else if (length(tau_names)) {
      stop(
        "Candidate reference-state seed found tau2 owned by a linear term: ",
        block$term_key, ".",
        call. = FALSE
      )
    }
  }
  
  if (any(!is.finite(theta))) {
    stop(
      "Candidate reference-state seed produced non-finite parameters.",
      call. = FALSE
    )
  }
  
  storage.mode(theta) <- "double"
  theta
}

# Extract response vector from a BAMLSS frame
extract_response <- function(frame) {
  y <- frame$y
  if (is.data.frame(y) && ncol(y) == 1L) y <- y[[1L]]
  as.numeric(y)
}

# Add external offsets on the predictor scale
add_offset <- function(eta, offset, n) {
  if (is.null(offset)) return(eta)
  
  if (!is.list(offset)) {
    if (length(eta) != 1L) {
      stop(
        "Unnamed numeric offset is only valid for a one-parameter family.",
        call. = FALSE
      )
    }
    offset <- stats::setNames(list(offset), names(eta)[1L])
  }
  
  for (p_name in intersect(names(eta), names(offset))) {
    z <- as.numeric(offset[[p_name]])
    if (length(z) == 1L) z <- rep(z, n)
    if (length(z) != n) {
      stop("Offset length mismatch for parameter ", p_name, call. = FALSE)
    }
    eta[[p_name]] <- eta[[p_name]] + z
  }
  eta
}

# Synchronize a parameter vector with BAMLSS term states
state_coefficients <- function(theta, frame) {
  template <- parameter_template(frame)
  theta <- merge_parameters(template, theta)
  
  if (!is.numeric(theta) || is.list(theta) ||
      is.null(names(theta)) || any(!nzchar(names(theta))) ||
      any(!is.finite(theta))) {
    stop(
      "State parameters must be a finite named numeric vector.", call. = FALSE
    )
  }
  
  x_theta <- frame$x
  for (p_name in names(x_theta)) {
    sc <- x_theta[[p_name]]$smooth.construct
    if (is.null(sc)) next
    for (label in names(sc)) {
      pm <- term_parameter_map(frame, p_name, label)
      local_pars <- theta[pm$global_names]
      names(local_pars) <- pm$local_names
      x_theta[[p_name]]$smooth.construct[[label]]$state$parameters <- local_pars
    }
  }
  x_theta
}

# Compute distributional predictors for a model state
state_predictors <- function(theta, frame, offset = NULL, x_theta = NULL) {
  if (is.null(x_theta)) x_theta <- state_coefficients(theta, frame)
  y_len <- length(extract_response(frame))
  
  eta <- list()
  for (p_name in names(x_theta)) {
    sc <- x_theta[[p_name]]$smooth.construct
    eta_p <- rep(0, y_len)
    if (!is.null(sc)) {
      for (label in names(sc)) {
        obj <- sc[[label]]
        pars <- obj$state$parameters
        
        # the parametric block is named "(Intercept)", "year_factor..."
        b_idx <- which(!startsWith(names(pars), "tau2"))
        if (length(b_idx) != ncol(obj$X)) {
          stop("Term '", label, "' of ", p_name, " has ", length(b_idx),
               " coefficients for ", ncol(obj$X), " design columns.",
               call. = FALSE)
        }
        eta_p <- eta_p + as.numeric(obj$X %*% as.numeric(pars[b_idx]))
      }
    }
    eta[[p_name]] <- eta_p
  }
  
  add_offset(eta, if (is.null(offset)) frame$model_offset else offset, y_len)
}

# Pointwise log density
evaluate_log_density <- function(y, par, fam_obj) {
  out <- tryCatch(
    fam_obj$d(y, par, log = TRUE),
    error = function(e) rep(NA_real_, length(y))
  )
  out <- as.numeric(out)
  if (length(out) == 1L && length(y) > 1L) out <- rep(out, length(y))
  out
}

# Direct BAMLSS term-prior decomposition
term_prior_components <- function(
    theta, frame, details = FALSE, x_theta = NULL) {
  theta <- flatten_parameters(theta)
  if (!length(theta) || is.null(names(theta)) ||
      any(!nzchar(names(theta))) || any(!is.finite(theta))) {
    return(list(
      structural_prior = -Inf,
      candidate_prior = -Inf,
      finite = FALSE,
      table = if (isTRUE(details)) data.frame() else NULL,
      error = "theta is not a finite named numeric vector"
    ))
  }
  
  sync_error <- NULL
  if (is.null(x_theta)) {
    x_theta <- tryCatch(
      state_coefficients(theta, frame),
      error = function(e) {
        sync_error <<- conditionMessage(e)
        NULL
      }
    )
  }
  if (is.null(x_theta)) {
    return(list(
      structural_prior = -Inf,
      candidate_prior = -Inf,
      finite = FALSE,
      table = if (isTRUE(details)) data.frame() else NULL,
      error = paste0("BAMLSS state synchronization failed: ", sync_error)
    ))
  }
  
  candidate_keys <- if (!is.null(frame$candidate_term_keys)) {
    frame$candidate_term_keys
  } else {
    character(0)
  }
  
  structural <- 0
  candidate <- 0
  rows <- if (isTRUE(details)) list() else NULL
  pos <- 0L
  global_error <- NULL
  
  for (p_name in names(x_theta)) {
    sc <- x_theta[[p_name]]$smooth.construct
    if (is.null(sc) || !length(sc)) next
    
    for (sj in names(sc)) {
      obj <- sc[[sj]]
      term_key <- paste(p_name, sj, sep = "::")
      kind <- if (term_key %in% candidate_keys) "candidate" else "structural"
      err <- NULL
      val <- NA_real_
      
      tryCatch({
        if (is.null(obj$prior) || !is.function(obj$prior)) {
          stop("term prior is missing or non-functional")
        }
        if (is.null(obj$state) || is.null(obj$state$parameters)) {
          stop("term state parameters are missing after synchronization")
        }
        
        z <- as.numeric(obj$prior(obj$state$parameters))
        if (!length(z)) stop("term prior returned length zero")
        val <- sum(z)
        if (!is.finite(val)) stop("term prior returned a non-finite value")
      }, error = function(e) {
        err <<- conditionMessage(e)
      })
      
      if (is.finite(val)) {
        if (identical(kind, "candidate")) {
          candidate <- candidate + val
        } else {
          structural <- structural + val
        }
      } else {
        if (identical(kind, "candidate")) {
          candidate <- -Inf
        } else {
          structural <- -Inf
        }
        if (is.null(global_error)) {
          global_error <- paste0(
            term_key, ": ", if (is.null(err)) "non-finite prior" else err
          )
        }
      }
      
      if (isTRUE(details)) {
        pos <- pos + 1L
        rows[[pos]] <- data.frame(
          parameter = p_name,
          term = sj,
          term_key = term_key,
          kind = kind,
          log_prior = val,
          finite = is.finite(val),
          error = if (is.null(err)) "" else err,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  tab <- if (isTRUE(details)) {
    if (length(rows)) do.call(rbind, rows) else data.frame()
  } else {
    NULL
  }
  
  list(
    structural_prior = as.numeric(structural),
    candidate_prior = as.numeric(candidate),
    finite = is.finite(structural) && is.finite(candidate),
    table = tab,
    error = global_error
  )
}

# Decompose the joint log posterior
posterior_components <- function(
    state, frame, fam_obj, offset = NULL,
    a_V = 0.5, a_A = 0.2,
    details = FALSE) {
  
  theta <- state$theta
  if (is.list(theta)) {
    theta <- tryCatch(
      merge_parameters(parameter_template(frame), theta),
      error = function(e) NULL
    )
  }
  
  bad_theta <- is.null(theta) || !is.numeric(theta) ||
    is.null(names(theta)) || any(!nzchar(names(theta))) ||
    any(!is.finite(theta))
  
  if (bad_theta) {
    return(list(
      logLik = -Inf,
      structural_prior = -Inf,
      candidate_prior = -Inf,
      model_prior = log_model_prior(state$inclusion, a_V, a_A),
      total = -Inf,
      prior_table = if (isTRUE(details)) data.frame() else NULL,
      likelihood = if (isTRUE(details)) list(error = "invalid theta") else NULL,
      theta = theta,
      error = "invalid theta"
    ))
  }
  
  sync_error <- NULL
  x_theta <- tryCatch(
    state_coefficients(theta, frame),
    error = function(e) {
      sync_error <<- conditionMessage(e)
      NULL
    }
  )
  
  if (is.null(x_theta)) {
    mp <- log_model_prior(state$inclusion, a_V, a_A)
    return(list(
      logLik = -Inf,
      structural_prior = -Inf,
      candidate_prior = -Inf,
      model_prior = mp,
      total = -Inf,
      prior_table = if (isTRUE(details)) data.frame() else NULL,
      likelihood = if (isTRUE(details)) {
        list(
          error = sync_error,
          pointwise_nonfinite = length(extract_response(frame)),
          pointwise_min = NA_real_,
          pointwise_max = NA_real_,
          y_range = range(extract_response(frame), na.rm = TRUE),
          eta_ranges = NULL,
          parameter_ranges = NULL
        )
      } else NULL,
      theta = theta,
      error = sync_error
    ))
  }
  
  y <- extract_response(frame)
  like_error <- NULL
  eta <- par <- NULL
  ll <- NA_real_
  
  tryCatch({
    eta <- state_predictors(theta, frame, offset, x_theta = x_theta)
    par <- fam_obj$map2par(eta)
    ll <- as.numeric(fam_obj$loglik(y, par))[1L]
    if (!is.finite(ll)) stop("family log-likelihood is non-finite")
  }, error = function(e) {
    like_error <<- conditionMessage(e)
  })
  
  pr <- term_prior_components(
    theta, frame, details = details, x_theta = x_theta
  )
  mp <- log_model_prior(state$inclusion, a_V, a_A)
  
  total <- if (is.finite(ll) && isTRUE(pr$finite) && is.finite(mp)) {
    ll + pr$structural_prior + pr$candidate_prior + mp
  } else {
    -Inf
  }
  
  like_details <- NULL
  if (isTRUE(details)) {
    ld <- if (!is.null(par)) {
      evaluate_log_density(y, par, fam_obj)
    } else {
      rep(NA_real_, length(y))
    }
    
    par_ranges <- if (!is.null(par)) {
      lapply(par, function(z) {
        z <- as.numeric(z)
        c(min = suppressWarnings(min(z, na.rm = TRUE)),
          max = suppressWarnings(max(z, na.rm = TRUE)),
          nonfinite = sum(!is.finite(z)))
      })
    } else NULL
    
    eta_ranges <- if (!is.null(eta)) {
      lapply(eta, function(z) {
        z <- as.numeric(z)
        c(min = suppressWarnings(min(z, na.rm = TRUE)),
          max = suppressWarnings(max(z, na.rm = TRUE)),
          nonfinite = sum(!is.finite(z)))
      })
    } else NULL
    
    like_details <- list(
      error = like_error,
      pointwise_nonfinite = sum(!is.finite(ld)),
      pointwise_min = if (any(is.finite(ld))) {
        min(ld[is.finite(ld)])
      } else {
        NA_real_
      },
      pointwise_max = if (any(is.finite(ld))) {
        max(ld[is.finite(ld)])
      } else {
        NA_real_
      },
      y_range = range(y, na.rm = TRUE),
      eta_ranges = eta_ranges,
      parameter_ranges = par_ranges
    )
  }
  
  list(
    logLik = if (is.finite(ll)) ll else -Inf,
    structural_prior = pr$structural_prior,
    candidate_prior = pr$candidate_prior,
    model_prior = mp,
    total = as.numeric(total),
    prior_table = pr$table,
    likelihood = like_details,
    theta = theta,
    error = if (!is.null(like_error)) like_error else pr$error
  )
}

# Print posterior diagnostics
print_posterior_diagnostics <- function(comp, label = "state") {
  fmt <- function(z) {
    if (length(z) && is.finite(z[1L])) {
      format(z[1L], digits = 8)
    } else {
      as.character(z[1L])
    }
  }
  
  message("  [Posterior diagnostics: ", label, "]")
  message("    logLik           = ", fmt(comp$logLik))
  message("    structural_prior = ", fmt(comp$structural_prior))
  message("    candidate_prior  = ", fmt(comp$candidate_prior))
  message("    model_prior      = ", fmt(comp$model_prior))
  message("    TOTAL            = ", fmt(comp$total))
  
  if (!is.null(comp$likelihood)) {
    if (!is.null(comp$likelihood$error) && nzchar(comp$likelihood$error)) {
      message("    likelihood error = ", comp$likelihood$error)
    }
    message(
      "    pointwise non-finite log densities = ",
      comp$likelihood$pointwise_nonfinite
    )
    message(
      "    y range = [", paste(comp$likelihood$y_range, collapse = ", "), "]"
    )
    
    if (!is.null(comp$likelihood$eta_ranges)) {
      for (nm in names(comp$likelihood$eta_ranges)) {
        z <- comp$likelihood$eta_ranges[[nm]]
        message(sprintf(
          "    eta[%s] range = [%g, %g], nonfinite=%d",
          nm, z["min"], z["max"], as.integer(z["nonfinite"])
        ))
      }
    }
    if (!is.null(comp$likelihood$parameter_ranges)) {
      for (nm in names(comp$likelihood$parameter_ranges)) {
        z <- comp$likelihood$parameter_ranges[[nm]]
        message(sprintf(
          "    par[%s] range = [%g, %g], nonfinite=%d",
          nm, z["min"], z["max"], as.integer(z["nonfinite"])
        ))
      }
    }
  }
  
  tab <- comp$prior_table
  if (!is.null(tab) && nrow(tab)) {
    bad <- tab[!tab$finite, , drop = FALSE]
    if (nrow(bad)) {
      message("    non-finite prior terms:")
      for (i in seq_len(nrow(bad))) {
        message(
          "      - ", bad$term_key[i], " [", bad$kind[i], "]: ",
          if (nzchar(bad$error[i])) {
            bad$error[i]
          } else {
            as.character(bad$log_prior[i])
          }
        )
      }
    }
  }
  invisible(comp)
}

# Exact offset-correct joint log posterior
target_log_posterior <- function(
    state, frame, fam_obj, offset = NULL, a_V = 0.5, a_A = 0.2) {
  
  posterior_components(
    state, frame, fam_obj, offset = offset,
    a_V = a_V, a_A = a_A, details = FALSE
  )$total
}

# Build a generic family-aware initial parameter vector
generic_initial_parameters <- function(frame, fam_obj) {
  start0 <- parameter_template(frame)
  
  itau <- grep("tau2", names(start0), fixed = TRUE)
  if (length(itau)) {
    bad_tau <- !is.finite(start0[itau]) | start0[itau] <= 0
    if (any(bad_tau)) start0[itau[bad_tau]] <- 1
  }
  
  y_raw <- frame$y
  if (is.data.frame(y_raw) && ncol(y_raw) == 1L) y_raw <- y_raw[[1L]]
  if (is.matrix(y_raw) && ncol(y_raw) == 1L) y_raw <- y_raw[, 1L]
  y_vec <- as.numeric(y_raw)
  
  links <- fam_obj$links
  if (is.null(names(links))) names(links) <- fam_obj$names
  
  for (p_name in fam_obj$names) {
    int_name <- paste0(p_name, ".p.(Intercept)")
    if (!int_name %in% names(start0)) next
    
    eta0 <- numeric(0)
    init_fun <- if (!is.null(fam_obj$initialize)) {
      fam_obj$initialize[[p_name]]
    } else {
      NULL
    }
    
    if (!is.null(init_fun) && is.function(init_fun)) {
      p0 <- tryCatch(init_fun(y_raw), error = function(e) NULL)
      if (!is.null(p0)) {
        p0 <- suppressWarnings(as.numeric(p0))
        p0 <- p0[is.finite(p0)]
        if (length(p0)) {
          lk <- tryCatch(link_object(links[[p_name]]), error = function(e) NULL)
          if (!is.null(lk) && is.function(lk$linkfun)) {
            eta0 <- tryCatch(
              suppressWarnings(as.numeric(lk$linkfun(p0))),
              error = function(e) numeric(0)
            )
            eta0 <- eta0[is.finite(eta0)]
          }
        }
      }
    }
    
    eta_center <- if (length(eta0)) stats::median(eta0) else 0
    off_raw <- NULL
    
    if (!is.null(frame$model_offset) && is.list(frame$model_offset) &&
        !is.null(frame$model_offset[[p_name]])) {
      off_raw <- suppressWarnings(as.numeric(frame$model_offset[[p_name]]))
    } else if (!is.null(frame$model_offset) && !is.list(frame$model_offset) &&
               length(fam_obj$names) == 1L) {
      off_raw <- suppressWarnings(as.numeric(frame$model_offset))
    }
    
    if (!is.null(off_raw) && length(off_raw)) {
      if (length(off_raw) == 1L && length(eta0) > 1L) {
        off_raw <- rep(off_raw, length(eta0))
      }
      if (length(eta0) && length(off_raw) == length(eta0) &&
          all(is.finite(off_raw))) {
        start0[[int_name]] <- stats::median(eta0 - off_raw)
      } else {
        z <- off_raw[is.finite(off_raw)]
        start0[[int_name]] <- eta_center - if (length(z)) mean(z) else 0
      }
    } else {
      start0[[int_name]] <- eta_center
    }
  }
  
  if (!is.numeric(start0) ||
      is.null(names(start0)) ||
      any(!is.finite(start0))) {
    stop(
      "Initialization did not give a finite named parameter vector.",
      call. = FALSE
    )
  }
  start0
}

# Verify prior decomposition on the same synchronized BAMLSS state
same_state_prior_check <- function(theta, frame, details = FALSE) {
  x_theta <- state_coefficients(theta, frame)
  pr <- term_prior_components(
    theta, frame, details = details, x_theta = x_theta
  )
  
  engine_prior <- tryCatch({
    ns <- asNamespace("bamlss")
    if (exists("get.log.prior", envir = ns, inherits = FALSE)) {
      as.numeric(
        get("get.log.prior", envir = ns, inherits = FALSE)(x_theta)
      )[1L]
    } else {
      pr$structural_prior + pr$candidate_prior
    }
  }, error = function(e) NA_real_)
  
  target_prior <- pr$structural_prior + pr$candidate_prior
  delta <- target_prior - engine_prior
  tol <- sqrt(.Machine$double.eps) * (1 + abs(engine_prior))
  
  list(
    consistent = is.finite(engine_prior) && isTRUE(pr$finite) &&
      is.finite(delta) && abs(delta) <= tol,
    engine_prior = engine_prior,
    target_prior = target_prior,
    delta = delta,
    tolerance = tol,
    prior = pr,
    x_theta = x_theta
  )
}

# Initialize the empty model
initialize_model_state <- function(
    inclusion, frame, fam_obj, runtime,
    optimizer_maxit = 400L,
    optimizer_eps = .Machine$double.eps^0.25,
    optimizer_nu = TRUE,
    a_V = 0.5, a_A = 0.2) {
  
  runtime$initial_optimizer_fits <- runtime$initial_optimizer_fits + 1L
  start0 <- generic_initial_parameters(frame, fam_obj)
  
  make_state <- function(theta) {
    st <- list(inclusion = inclusion, theta = theta, logPost = -Inf)
    st$logPost <- target_log_posterior(
      st, frame, fam_obj, offset = frame$model_offset,
      a_V = a_V, a_A = a_A
    )
    st
  }
  
  st_start <- make_state(start0)
  
  if (is.finite(st_start$logPost)) {
    chk <- same_state_prior_check(st_start$theta, frame)
    if (!isTRUE(chk$consistent)) {
      stop(
        "Prior decomposition disagrees with BAMLSS get.log.prior().",
        call. = FALSE
      )
    }
  }
  
  opt_error <- NULL
  opt <- tryCatch(
    suppressWarnings(
      bamlss::opt_bfit(
        frame$x, frame$y, fam_obj,
        weights = frame$weights,
        offset = frame$model_offset,
        start = start0,
        maxit = as.integer(optimizer_maxit),
        eps = optimizer_eps,
        nu = optimizer_nu,
        verbose = FALSE
      )
    ),
    error = function(e) {
      opt_error <<- conditionMessage(e)
      NULL
    }
  )
  
  st_opt <- NULL
  if (!is.null(opt) && !is.null(opt$parameters)) {
    theta_opt <- tryCatch(
      merge_parameters(start0, opt$parameters),
      error = function(e) NULL
    )
    if (!is.null(theta_opt) && is.numeric(theta_opt) && !is.list(theta_opt) &&
        all(is.finite(theta_opt))) {
      st_opt <- make_state(theta_opt)
    }
  }
  
  candidates <- list(family_start = st_start)
  if (!is.null(st_opt)) candidates$optimizer_export <- st_opt
  
  finite <- vapply(candidates, function(z) is.finite(z$logPost), logical(1))
  if (!any(finite)) {
    dg <- tryCatch(
      term_prior_components(st_start$theta, frame, details = TRUE),
      error = function(e) NULL
    )
    tab <- if (!is.null(dg) && !is.null(dg$table) && nrow(dg$table)) {
      paste(utils::capture.output(print(
        dg$table[, c("term_key", "log_prior", "error")], row.names = FALSE
      )), collapse = "\n")
    } else {
      "<prior table unavailable>"
    }
    stop(
      "No finite initial state: the log-posterior is not finite at the start.",
      "\nTerm priors at the start (finite priors point to the likelihood):\n",
      tab, call. = FALSE
    )
  }
  
  cand_f <- candidates[finite]
  lp <- vapply(cand_f, function(z) z$logPost, numeric(1))
  chosen_name <- names(cand_f)[which.max(lp)]
  st <- cand_f[[chosen_name]]
  
  st
}

# sam_GMCMC names the parametric block "<par>.p.<coef>" while the frame
# template names it "<par>.s.model.matrix.<coef>"; rename the draw to match
align_draw_names <- function(draw, template_names) {
  nm <- names(draw)
  for (i in which(!nm %in% template_names)) {
    cand <- sub("^([^.]+)[.]p[.]", "\\1.s.model.matrix.", nm[i])
    if (cand != nm[i] && cand %in% template_names) nm[i] <- cand
  }
  names(draw) <- nm
  draw
}

# Inverse of align_draw_names: set.starting.values only reads the parametric
# block as "<par>.p.<coef>" and the block has no tau2 in its state
to_sampler_start <- function(start) {
  if (!is.numeric(start) || is.null(names(start))) return(start)
  nm <- names(start)
  mm <- grepl("^[^.]+[.]s[.]model[.]matrix[.]", nm)
  nm[mm] <- sub("^([^.]+)[.]s[.]model[.]matrix[.]", "\\1.p.", nm[mm])
  names(start) <- nm
  start[!(mm & grepl("tau2", nm, fixed = TRUE))]
}

# Stochastic within-model update
within_model_update <- function(
    state,
    frame,
    fam_obj,
    runtime,
    within_sweeps = 5L,
    a_V = 0.5,
    a_A = 0.2) {
  
  within_sweeps <- max(1L, as.integer(within_sweeps))
  runtime$within_attempts <- runtime$within_attempts + 1L
  sampler_error <- NULL
  
  sampled <- tryCatch(
    bamlss::sam_GMCMC(
      frame$x,
      frame$y,
      fam_obj,
      start = to_sampler_start(flatten_parameters(state$theta)),
      weights = frame$weights,
      offset = frame$model_offset,
      n.iter = within_sweeps,
      burnin = 0,
      thin = 1,
      verbose = FALSE,
      step = max(2L, within_sweeps + 1L),
      propose = "iwlsC_gp"
    ),
    error = function(e) {
      sampler_error <<- conditionMessage(e)
      NULL
    }
  )
  
  fail <- function(message) {
    runtime$within_failures <- runtime$within_failures + 1L
    record_runtime_message(runtime, "within_error_messages", message)
    state$logPost <- target_log_posterior(
      state, frame, fam_obj, offset = frame$model_offset,
      a_V = a_V, a_A = a_A
    )
    state
  }
  
  if (is.null(sampled)) return(
    fail(
      if (!is.null(sampler_error)) sampler_error else "sam_GMCMC returned NULL"
    )
  )
  
  matrix_error <- NULL
  sampled_matrix <- tryCatch({
    if (inherits(sampled, "mcmc.list")) {
      as.matrix(sampled[[length(sampled)]])
    } else {
      as.matrix(sampled)
    }
  }, error = function(e) {
    matrix_error <<- conditionMessage(e)
    NULL
  })
  
  if (is.null(sampled_matrix) ||
      !nrow(sampled_matrix) ||
      !ncol(sampled_matrix)) {
    if (is.null(matrix_error)) matrix_error <- "sam_GMCMC returned no draws"
    return(fail(matrix_error))
  }
  
  draw <- suppressWarnings(
    as.numeric(sampled_matrix[nrow(sampled_matrix), , drop = TRUE])
  )
  names(draw) <- colnames(sampled_matrix)
  
  if (!length(draw) || is.null(names(draw)) || any(!is.finite(draw))) {
    return(fail("sam_GMCMC returned a non-finite or unnamed terminal draw"))
  }
  
  template <- parameter_template(frame)
  draw <- align_draw_names(draw, names(template))
  theta <- tryCatch({
    current <- merge_parameters(template, state$theta)
    merge_parameters(current, draw)
  }, error = function(e) NULL
  )
  
  if (is.null(theta) || any(!is.finite(theta))) {
    return(
      fail("sam_GMCMC draw could not be mapped to a finite parameter state")
    )
  }
  
  if (nzchar(Sys.getenv("RBST_TRACE"))) {
    ll0 <- posterior_components(
      state, frame, fam_obj, offset = frame$model_offset,
      a_V = a_V, a_A = a_A
    )$logLik
    llv <- if ("logLik" %in% colnames(sampled_matrix)) {
      sampled_matrix[, "logLik"]
    } else NA_real_
    dl <- diff(c(ll0, llv))
    w <- if (all(is.na(dl))) 1L else which.min(dl)
    acc <- grep("accepted$", colnames(sampled_matrix), value = TRUE)
    hit <- acc[sampled_matrix[w, acc] > 0.5]
    message(sprintf(
      "[trace]   sweeps ll=%s | worst sweep %d dll=%.4g accepted: %s",
      paste(sprintf("%.6g", c(ll0, llv)), collapse = ","), w, dl[w],
      paste(sub("[.]accepted$", "", hit), collapse = ",")
    ))
  }
  
  updated_state <- state
  updated_state$theta <- theta
  
  updated_state$logPost <- target_log_posterior(
    updated_state, frame, fam_obj, offset = frame$model_offset,
    a_V = a_V, a_A = a_A
  )
  
  if (!is.finite(updated_state$logPost)) {
    return(
      fail("sam_GMCMC terminal draw is outside the specified posterior support")
    )
  }
  
  for (parameter in fam_obj$names) {
    updated_state <- update_linear_candidates_ess(
      state = updated_state,
      frame = frame,
      fam_obj = fam_obj,
      runtime = runtime,
      parameter = parameter,
      a_V = a_V,
      a_A = a_A
    )
  }
  
  runtime$within_successes <- runtime$within_successes + 1L
  updated_state
}

# Numerical predictor-scale score and information
numerical_derivatives <- function(y, eta, fam_obj, eps = 1e-4) {
  K <- length(eta)
  nm <- names(eta)
  n <- length(y)
  score <- matrix(0, n, K, dimnames = list(NULL, nm))
  W <- vector("list", K)
  names(W) <- nm
  for (k in seq_len(K)) W[[k]] <- vector("list", K)
  for (k in seq_len(K)) names(W[[k]]) <- nm
  
  ld0 <- evaluate_log_density(y, fam_obj$map2par(eta), fam_obj)
  for (k in seq_len(K)) {
    ep <- eta
    em <- eta
    ep[[k]] <- ep[[k]] + eps
    em[[k]] <- em[[k]] - eps
    ldp <- evaluate_log_density(y, fam_obj$map2par(ep), fam_obj)
    ldm <- evaluate_log_density(y, fam_obj$map2par(em), fam_obj)
    score[, k] <- (ldp - ldm) / (2 * eps)
    W[[k]][[k]] <- -(ldp - 2 * ld0 + ldm) / (eps^2)
  }
  list(score = score, W = W)
}

# Score and full cross-parameter observed information
score_information <- function(
    state, frame, fam_obj, full_cross_information = TRUE,
    fd_eps = 1e-4, info_floor = 1e-8) {
  
  y <- extract_response(frame)
  eta <- state_predictors(state$theta, frame, frame$model_offset)
  par <- fam_obj$map2par(eta)
  K <- length(fam_obj$names)
  nm <- fam_obj$names
  n <- length(y)
  
  score <- matrix(NA_real_, n, K, dimnames = list(NULL, nm))
  W <- vector("list", K)
  names(W) <- nm
  for (k in seq_len(K)) {
    W[[k]] <- vector("list", K)
    names(W[[k]]) <- nm
  }
  
  need_num_score <- rep(TRUE, K)
  need_num_hess <- rep(TRUE, K)
  
  for (k in seq_len(K)) {
    p_name <- nm[k]
    sf <- fam_obj$score[[p_name]]
    hf <- fam_obj$hess[[p_name]]
    
    if (is.function(sf) && !isTRUE(attr(sf, "numerical_derivative"))) {
      z <- tryCatch(as.numeric(sf(y, par)), error = function(e) NULL)
      if (!is.null(z) && length(z) == n && all(is.finite(z))) {
        score[, k] <- z
        need_num_score[k] <- FALSE
      }
    }
    if (is.function(hf) && !isTRUE(attr(hf, "numerical_derivative"))) {
      z <- tryCatch(as.numeric(hf(y, par)), error = function(e) NULL)
      if (!is.null(z) && length(z) == n && all(is.finite(z))) {
        W[[k]][[k]] <- z
        need_num_hess[k] <- FALSE
      }
    }
  }
  
  if (any(need_num_score | need_num_hess)) {
    ld0 <- evaluate_log_density(y, fam_obj$map2par(eta), fam_obj)
    for (k in which(need_num_score | need_num_hess)) {
      ep <- eta
      em <- eta
      ep[[k]] <- ep[[k]] + fd_eps
      em[[k]] <- em[[k]] - fd_eps
      ldp <- evaluate_log_density(y, fam_obj$map2par(ep), fam_obj)
      ldm <- evaluate_log_density(y, fam_obj$map2par(em), fam_obj)
      if (need_num_score[k]) score[, k] <- (ldp - ldm) / (2 * fd_eps)
      if (need_num_hess[k]) W[[k]][[k]] <- -(ldp - 2 * ld0 + ldm) / (fd_eps^2)
    }
  }
  
  for (k in seq_len(K)) {
    zsc <- as.numeric(score[, k])
    zsc[!is.finite(zsc)] <- 0
    score[, k] <- zsc
    
    z <- as.numeric(W[[k]][[k]])
    if (length(z) != n) z <- rep(info_floor, n)
    z[!is.finite(z)] <- info_floor
    W[[k]][[k]] <- pmax(z, info_floor)
  }
  
  if (isTRUE(full_cross_information) && K > 1L) {
    for (k in seq_len(K - 1L)) for (l in (k + 1L):K) {
      epp <- eta
      epm <- eta
      emp <- eta
      emm <- eta
      epp[[k]] <- epp[[k]] + fd_eps
      epp[[l]] <- epp[[l]] + fd_eps
      epm[[k]] <- epm[[k]] + fd_eps
      epm[[l]] <- epm[[l]] - fd_eps
      emp[[k]] <- emp[[k]] - fd_eps
      emp[[l]] <- emp[[l]] + fd_eps
      emm[[k]] <- emm[[k]] - fd_eps
      emm[[l]] <- emm[[l]] - fd_eps
      lpp <- evaluate_log_density(y, fam_obj$map2par(epp), fam_obj)
      lpm <- evaluate_log_density(y, fam_obj$map2par(epm), fam_obj)
      lmp <- evaluate_log_density(y, fam_obj$map2par(emp), fam_obj)
      lmm <- evaluate_log_density(y, fam_obj$map2par(emm), fam_obj)
      z <- -(lpp - lpm - lmp + lmm) / (4 * fd_eps^2)
      z[!is.finite(z)] <- 0
      lim <- 0.99 * sqrt(W[[k]][[k]] * W[[l]][[l]])
      z <- pmax(pmin(z, lim), -lim)
      W[[k]][[l]] <- z
      W[[l]][[k]] <- z
    }
  } else if (K > 1L) {
    for (k in seq_len(K)) for (l in seq_len(K)) if (k != l) {
      W[[k]][[l]] <- rep(0, n)
    }
  }
  
  list(score = score, W = W, eta = eta, par = par)
}

# Candidate prior precision conditional on tau2
candidate_prior_precision <- function(
    block, tau2 = NULL,
    c0 = 316.227766, cplus = 316.227766,
    c_linear = 2.0, floor = 1e-10) {
  
  d <- ncol(block$X)
  if (!isTRUE(block$smooth) || !length(block$S)) {
    return(diag(1 / c_linear^2, d))
  }
  
  S <- lapply(block$S, function(z) symmetrize_matrix(as.matrix(z)))
  Ksum <- Reduce(`+`, S)
  ee <- eigen(symmetrize_matrix(Ksum), symmetric = TRUE)
  tol <- max(floor, max(abs(ee$values), 1) * 1e-10)
  pos <- which(ee$values > tol)
  nul <- which(ee$values <= tol)
  
  Q <- matrix(0, d, d)
  if (length(nul)) {
    Q <- Q + tcrossprod(ee$vectors[, nul, drop = FALSE]) / (c0^2)
  }
  if (length(pos)) {
    Q <- Q + tcrossprod(ee$vectors[, pos, drop = FALSE]) / (cplus^2)
  }
  
  if (is.null(tau2) || !length(tau2)) {
    tau2 <- rep(1, length(S))
  } else if (any(!is.finite(tau2)) || any(tau2 <= 0)) {
    stop(
      "Candidate tau2 must be finite and strictly positive (got ",
      paste(format(tau2, digits = 3), collapse = ", "), "; called from ",
      paste(deparse(sys.call(-1L))[1L]), ").", call. = FALSE
    )
  }
  tau2 <- as.numeric(tau2)
  if (length(tau2) == 1L && length(S) > 1L) tau2 <- rep(tau2, length(S))
  if (length(tau2) != length(S)) {
    stop(
      "tau2 dimension does not match candidate penalty blocks.", call. = FALSE
    )
  }
  
  lambda <- 1 / tau2
  for (j in seq_along(S)) Q <- Q + lambda[j] * S[[j]]
  make_positive_definite(Q, floor)
}

# Candidate prior covariance
candidate_prior_covariance <- function(
    block, c0 = 316.227766, cplus = 316.227766,
    a_lambda = 1, b_lambda = 1, c_linear = 2.0,
    floor = 1e-10) {
  
  d <- ncol(block$X)
  if (!isTRUE(block$smooth) || !length(block$S)) return(diag(c_linear^2, d))
  lam_bar <- a_lambda / b_lambda
  if (!is.finite(lam_bar) || lam_bar <= 0) {
    stop("Invalid Gamma prior mean for lambda.", call. = FALSE)
  }
  tau_ref <- rep(1 / lam_bar, length(block$S))
  Q <- candidate_prior_precision(
    block, tau2 = tau_ref,
    c0 = c0, cplus = cplus, c_linear = c_linear,
    floor = floor
  )
  invert_positive_definite(Q, floor)
}

# Local log Bayes factor
local_log_bayes_factor <- function(Ueff, Ieff, Sigma0, floor = 1e-10) {
  d <- length(Ueff)
  if (!d) return(-Inf)
  Ieff <- make_positive_definite(Ieff, floor)
  Sigma0 <- make_positive_definite(Sigma0, floor)
  A <- diag(d) + Sigma0 %*% Ieff
  B <- Ieff + invert_positive_definite(Sigma0, floor)
  val <- -0.5 * logdet_positive_definite(A, floor) +
    0.5 * as.numeric(
      crossprod(Ueff, invert_positive_definite(B, floor) %*% Ueff)
    )
  if (!is.finite(val)) -Inf else val
}

# Validate the trans-dimensional birth geometry
validate_birth_geometry <- function(
    state, frame, data, fam_obj, formula_base,
    admissibility, offset, runtime,
    c_linear = 1.5,
    c0 = 316.227766,
    cplus = 316.227766,
    a_lambda = 1,
    b_lambda = 1,
    proposal_ridge = 1e-6,
    full_cross_information = TRUE,
    fd_eps = 1e-4,
    epsilon = 0.25,
    omega = 0.25,
    a_V = 0.5,
    a_A = 0.2) {
  
  admissible <- which(admissibility, arr.ind = TRUE)
  admissible_birth <- which(
    birth_mask(state$inclusion, admissibility), arr.ind = TRUE
  )
  if (!nrow(admissible_birth)) {
    stop(
      "No hierarchy-admissible birth exists from the initial state.",
      call. = FALSE
    )
  }
  
  deriv_error <- NULL
  deriv <- tryCatch(
    score_information(
      state, frame, fam_obj,
      full_cross_information = full_cross_information,
      fd_eps = fd_eps
    ),
    error = function(e) {
      deriv_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(deriv)) {
    stop("Score/information validation failed: ", deriv_error, call. = FALSE)
  }
  
  if (any(!is.finite(deriv$score))) {
    stop(
      "Score/information validation produced non-finite values.", call. = FALSE
    )
  }
  
  rows <- vector("list", nrow(admissible))
  cand_map <- candidate_map(rownames(state$inclusion), data)
  base_template <- parameter_template(frame)
  
  for (ii in seq_len(nrow(admissible))) {
    r <- admissible[ii, 1L]
    k <- admissible[ii, 2L]
    component <- rownames(state$inclusion)[r]
    parameter <- colnames(state$inclusion)[k]
    initial_birth_admissible <- any(admissible_birth[, 1L] == r &
                                      admissible_birth[, 2L] == k)
    
    block_err <- frame_err <- target_err <- NULL
    block <- tryCatch(
      get_candidate_block(
        r, k, state$inclusion, formula_base, data, fam_obj, offset, runtime,
        c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda
      ),
      error = function(e) {
        block_err <<- conditionMessage(e)
        NULL
      }
    )
    
    block_valid <- !is.null(block)
    if (block_valid) {
      checks <- c(
        nrow(block$X) == nrow(data),
        ncol(block$X) > 0L,
        all(is.finite(block$X)),
        length(block$beta_names) == ncol(block$X),
        length(block$label) == 1L,
        nzchar(block$label)
      )
      block_valid <- all(checks)
    }
    
    frame_available <- target_finite <- FALSE
    test_frame <- NULL
    test_inc <- state$inclusion
    
    if (block_valid) {
      test_inc[r, k] <- TRUE
      if (identical(cand_map$component[r], "nonlinear")) {
        v <- cand_map$variable[r]
        rL <- which(cand_map$variable == v & cand_map$component == "linear")
        if (length(rL) == 1L && isTRUE(admissibility[rL, k])) {
          test_inc[rL, k] <- TRUE
        }
      }
      
      if (hierarchy_valid(test_inc)) {
        test_frame <- tryCatch(
          get_model_frame(
            test_inc, formula_base, data, fam_obj, offset, runtime,
            c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda,
            c0 = c0, cplus = cplus
          ),
          error = function(e) NULL
        )
      }
      frame_available <- !is.null(test_frame)
    }
    
    if (frame_available) {
      theta_test <- tryCatch({
        z <- merge_parameters(parameter_template(test_frame), state$theta)
        seed_candidate_reference_state(
          z, test_frame, a_lambda = a_lambda, b_lambda = b_lambda,
          beta_value = 0
        )
      }, error = function(e) NULL)
      
      if (!is.null(theta_test) && all(is.finite(theta_test))) {
        test_state <- list(
          inclusion = test_inc, theta = theta_test, logPost = NA_real_
        )
        pc <- tryCatch(
          posterior_components(
            test_state, test_frame, fam_obj,
            offset = test_frame$model_offset,
            a_V = a_V, a_A = a_A,
            details = FALSE
          ),
          error = function(e) NULL
        )
        if (!is.null(pc) && is.finite(pc$total)) {
          target_finite <- TRUE
        }
      }
    }
    
    err <- paste(c(block_err, frame_err, target_err), collapse = " | ")
    geometry_valid <-
      isTRUE(block_valid) && isTRUE(frame_available) && isTRUE(target_finite)
    
    rows[[ii]] <- data.frame(
      Component = component,
      Parameter = parameter,
      Initial_birth_admissible = initial_birth_admissible,
      Term_beta_n = if (!is.null(block)) {
        length(block$beta_names)
      } else {
        NA_integer_
      },
      Term_tau_n = if (!is.null(block)) {
        length(block$tau_names)
      } else {
        NA_integer_
      },
      Term_penalties = if (!is.null(block)) length(block$S) else NA_integer_,
      Term_fixed = if (!is.null(block)) isTRUE(block$fixed) else NA,
      Term_fxsp = if (!is.null(block)) isTRUE(block$fxsp) else NA,
      Block_valid = block_valid,
      Frame_valid = frame_available,
      Target_valid = target_finite,
      Valid = geometry_valid,
      Error = err,
      stringsAsFactors = FALSE
    )
  }
  
  tab <- do.call(rbind, rows)
  if (any(!tab$Valid)) {
    bad <- tab[!tab$Valid, , drop = FALSE]
    print(bad, row.names = FALSE)
    stop("Candidate validation failed. Examine the table above.", call. = FALSE)
  }
  
  geometry_error <- NULL
  geometry <- tryCatch(
    evaluate_proposal_geometry(
      state, frame, data, fam_obj, formula_base, admissibility, offset, runtime,
      epsilon = epsilon, omega = omega,
      a_V = a_V, a_A = a_A,
      c_linear = c_linear, c0 = c0, cplus = cplus,
      a_lambda = a_lambda, b_lambda = b_lambda,
      proposal_ridge = proposal_ridge,
      full_cross_information = full_cross_information,
      fd_eps = fd_eps
    ),
    error = function(e) {
      geometry_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(geometry) || !length(geometry$q_b)) {
    stop(
      "Initial birth-geometry validation failed: ", geometry_error,
      call. = FALSE
    )
  }
  
  attr(tab, "initial_q_b") <- geometry$q_b
  invisible(tab)
}

# Evaluate efficient geometry and informed birth proposal
evaluate_proposal_geometry <- function(
    state, frame, data, fam_obj, formula_base, admissibility, offset, runtime,
    epsilon = 0.25, omega = 0.25,
    a_V = 0.5, a_A = 0.2,
    c_linear = 2.0, c0 = 316.227766, cplus = 316.227766,
    a_lambda = 1, b_lambda = 1,
    proposal_ridge = 1e-6,
    full_cross_information = TRUE,
    fd_eps = 1e-4) {
  
  inactive <- which(birth_mask(state$inclusion, admissibility), arr.ind = TRUE)
  if (!nrow(inactive)) return(NULL)
  
  deriv <- score_information(
    state, frame, fam_obj,
    full_cross_information = full_cross_information,
    fd_eps = fd_eps
  )
  K <- ncol(state$inclusion)
  
  XA <- vector("list", K)
  active_meta <- vector("list", K)
  dims <- integer(K)
  for (k in seq_len(K)) {
    rr <- which(state$inclusion[, k] & admissibility[, k])
    if (length(rr)) {
      blocks_k <- lapply(rr, function(r) {
        meta_key <- paste(r, k, sep = ":")
        bl <- frame$candidate_metadata[[meta_key]]
        if (is.null(bl)) {
          stop(
            "Active candidate metadata are missing from the current frame for ",
            rownames(state$inclusion)[r], " -> ", colnames(state$inclusion)[k],
            ".",
            call. = FALSE
          )
        }
        bl
      })
      XA[[k]] <- do.call(cbind, lapply(blocks_k, `[[`, "X"))
      active_meta[[k]] <- blocks_k
    } else {
      XA[[k]] <- matrix(0, nrow(data), 0L)
      active_meta[[k]] <- list()
    }
    dims[k] <- ncol(XA[[k]])
  }
  
  starts <- cumsum(c(1L, head(dims, -1L)))
  dA <- sum(dims)
  UA <- numeric(dA)
  IAA <- matrix(0, dA, dA)
  QA <- matrix(0, dA, dA)
  prior_grad_A <- numeric(dA)
  if (dA) {
    for (k in seq_len(K)) if (dims[k]) {
      ik <- starts[k]:(starts[k] + dims[k] - 1L)
      UA[ik] <- as.numeric(crossprod(XA[[k]], deriv$score[, k]))
      for (l in seq_len(K)) if (dims[l]) {
        il <- starts[l]:(starts[l] + dims[l] - 1L)
        IAA[ik, il] <- crossprod(XA[[k]], deriv$W[[k]][[l]] * XA[[l]])
      }
      
      local_start <- 1L
      for (bl in active_meta[[k]]) {
        dd <- ncol(bl$X)
        ig <- (starts[k] + local_start - 1L):(starts[k] + local_start + dd - 2L)
        beta <- as.numeric(state$theta[bl$beta_names])
        tau2 <- if (bl$smooth && length(bl$tau_names)) {
          as.numeric(state$theta[bl$tau_names])
        } else {
          numeric(0)
        }
        tau_ok <- !length(tau2) || (all(is.finite(tau2)) && all(tau2 > 0))
        if (length(beta) == dd && all(is.finite(beta)) && tau_ok) {
          Qbl <- candidate_prior_precision(
            bl, tau2 = tau2, c0 = c0, cplus = cplus,
            c_linear = c_linear, floor = proposal_ridge
          )
          QA[ig, ig] <- QA[ig, ig, drop = FALSE] + Qbl
          prior_grad_A[ig] <- -as.numeric(Qbl %*% beta)
        }
        local_start <- local_start + dd
      }
    }
    IAA_lik <- make_positive_definite(
      IAA + diag(proposal_ridge, dA), proposal_ridge
    )
    invIAA <- invert_positive_definite(IAA_lik, proposal_ridge)
    invIAA_UA <- invIAA %*% UA
    
    HAA <- make_positive_definite(
      IAA + QA + diag(proposal_ridge, dA), proposal_ridge
    )
    invHAA <- invert_positive_definite(HAA, proposal_ridge)
    GA <- UA + prior_grad_A
    invHAA_GA <- invHAA %*% GA
  } else {
    invIAA <- matrix(0, 0, 0)
    invIAA_UA <- numeric(0)
    invHAA <- matrix(0, 0, 0)
    invHAA_GA <- numeric(0)
  }
  
  candidates <- vector("list", nrow(inactive))
  log_w <- rep(-Inf, nrow(inactive))
  for (ii in seq_len(nrow(inactive))) {
    r <- inactive[ii, 1L]
    k <- inactive[ii, 2L]
    block <- get_candidate_block(
      r, k, state$inclusion, formula_base, data, fam_obj, offset, runtime,
      c_linear, a_lambda, b_lambda
    )
    Xc <- block$X
    Uc <- as.numeric(crossprod(Xc, deriv$score[, k]))
    Icc <- crossprod(Xc, deriv$W[[k]][[k]] * Xc)
    
    if (dA) {
      IcA <- matrix(0, nrow = ncol(Xc), ncol = dA)
      for (l in seq_len(K)) if (dims[l]) {
        il <- starts[l]:(starts[l] + dims[l] - 1L)
        IcA[, il] <- crossprod(Xc, deriv$W[[k]][[l]] * XA[[l]])
      }
      Ueff <- Uc - as.numeric(IcA %*% invIAA_UA)
      Ieff <- Icc - IcA %*% invIAA %*% t(IcA)
      Gperp <- Uc - as.numeric(IcA %*% invHAA_GA)
      Hbase <- Icc - IcA %*% invHAA %*% t(IcA)
    } else {
      Ueff <- Uc
      Ieff <- Icc
      Gperp <- Uc
      Hbase <- Icc
    }
    Ieff <- make_positive_definite(
      Ieff + diag(proposal_ridge, ncol(Xc)), proposal_ridge
    )
    Hbase <- make_positive_definite(
      Hbase + diag(proposal_ridge, ncol(Xc)), proposal_ridge
    )
    
    Gcond <- Uc
    Hcond <- make_positive_definite(
      Icc + diag(proposal_ridge, ncol(Xc)), proposal_ridge
    )
    
    Sigma0 <- candidate_prior_covariance(
      block, c0 = c0, cplus = cplus,
      a_lambda = a_lambda, b_lambda = b_lambda,
      c_linear = c_linear
    )
    logBF <- local_log_bayes_factor(Ueff, Ieff, Sigma0, floor = proposal_ridge)
    lpr <- log_model_prior_ratio(state$inclusion, r, k, a_V, a_A)
    log_w[ii] <- omega * (logBF + lpr)
    candidates[[ii]] <- list(
      row = r, col = k, block = block,
      Ueff = Ueff, Ieff = Ieff,
      Gperp = Gperp, Hbase = Hbase,
      Gcond = Gcond, Hcond = Hcond,
      logBF = logBF, log_prior_ratio = lpr
    )
  }
  
  finite <- is.finite(log_w)
  if (!any(finite)) {
    w <- rep(1 / length(log_w), length(log_w))
  } else {
    mm <- max(log_w[finite])
    w <- exp(pmin(log_w - mm, 700))
    w[!is.finite(w)] <- 0
    if (sum(w) <= 0) w[] <- 1
    w <- w / sum(w)
  }
  q <- (1 - epsilon) * w + epsilon / length(w)
  q <- q / sum(q)
  
  key <- vapply(
    candidates, function(z) paste(z$row, z$col, sep = ":"), character(1)
  )
  names(candidates) <- key
  names(q) <- key
  list(inactive_idx = inactive, q_b = q, candidates = candidates)
}

# Smoothing-variance proposal density
log_tau2_proposal <- function(tau2, a_lambda = 1, b_lambda = 1) {
  tau2 <- as.numeric(tau2)
  if (!length(tau2)) return(0)
  if (any(!is.finite(tau2)) || any(tau2 <= 0)) return(-Inf)
  sum(
    stats::dgamma(
      1 / tau2, shape = a_lambda, rate = b_lambda, log = TRUE
    ) - 2 * log(tau2)
  )
}

# Stable two-term log-sum-exp
logsumexp_pair <- function(a, b) {
  if (!is.finite(a) && !is.finite(b)) return(-Inf)
  m <- max(a, b)
  m + log(exp(a - m) + exp(b - m))
}

# Local marginal log posterior for spline precision on log(lambda)
local_log_lambda_profile <- function(
    z, geom_entry, block,
    c_linear = 2.0, c0 = 316.227766, cplus = 316.227766,
    a_lambda = 1, b_lambda = 1,
    proposal_ridge = 1e-6) {
  
  z <- as.numeric(z)
  if (!length(z) || any(!is.finite(z))) return(-Inf)
  lambda <- exp(z)
  if (any(!is.finite(lambda)) || any(lambda <= 0)) return(-Inf)
  tau2 <- exp(-z)
  if (any(!is.finite(tau2)) || any(tau2 <= 0)) return(-Inf)
  
  Q <- tryCatch(
    candidate_prior_precision(
      block, tau2 = tau2,
      c0 = c0, cplus = cplus, c_linear = c_linear,
      floor = proposal_ridge
    ),
    error = function(e) NULL
  )
  if (is.null(Q)) return(-Inf)
  
  H <- make_positive_definite(
    geom_entry$Hbase + Q + diag(proposal_ridge, ncol(block$X)),
    proposal_ridge
  )
  V <- invert_positive_definite(H, proposal_ridge)
  G <- as.numeric(geom_entry$Gperp)
  
  log_marg_beta <- 0.5 * logdet_positive_definite(Q, proposal_ridge) -
    0.5 * logdet_positive_definite(H, proposal_ridge) +
    0.5 * as.numeric(crossprod(G, V %*% G))
  
  log_prior_z <- sum(
    stats::dgamma(lambda, shape = a_lambda, rate = b_lambda, log = TRUE) + z
  )
  
  out <- log_marg_beta + log_prior_z
  if (is.finite(out)) as.numeric(out) else -Inf
}

# Conditional-maximization Laplace proposal for log(lambda)
laplace_log_lambda <- function(
    geom_entry, block,
    c_linear = 2.0, c0 = 316.227766, cplus = 316.227766,
    a_lambda = 1, b_lambda = 1,
    proposal_ridge = 1e-6,
    maxit = 25L, reltol = 1e-6) {
  
  nt <- length(block$S)
  if (!isTRUE(block$smooth) || nt < 1L) {
    return(
      list(
        available = FALSE, mean = numeric(0), cov = matrix(0, 0, 0),
        convergence = NA_integer_
      )
    )
  }
  
  lam0 <- a_lambda / b_lambda
  if (!is.finite(lam0) || lam0 <= 0) lam0 <- 1
  z0 <- rep(log(lam0), nt)
  
  fn <- function(z) {
    lp <- local_log_lambda_profile(
      z, geom_entry, block,
      c_linear = c_linear, c0 = c0, cplus = cplus,
      a_lambda = a_lambda, b_lambda = b_lambda,
      proposal_ridge = proposal_ridge
    )
    if (is.finite(lp)) -lp else .Machine$double.xmax / 100
  }
  
  opt <- tryCatch(
    stats::optim(
      z0, fn, method = "BFGS",
      control = list(maxit = as.integer(maxit), reltol = reltol)
    ),
    error = function(e) NULL
  )
  if (is.null(opt) ||
      opt$convergence != 0L ||
      any(!is.finite(opt$par)) ||
      !is.finite(opt$value)) {
    return(
      list(
        available = FALSE, mean = z0, cov = diag(1, nt),
        convergence = if (is.null(opt)) NA_integer_ else opt$convergence
      )
    )
  }
  
  HH <- tryCatch(stats::optimHess(opt$par, fn), error = function(e) NULL)
  if (is.null(HH) || any(!is.finite(HH))) {
    return(
      list(
        available = FALSE, mean = opt$par, cov = diag(1, nt),
        convergence = opt$convergence
      )
    )
  }
  HH <- symmetrize_matrix(HH)
  ev <- tryCatch(
    eigen(HH, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, nt)
  )
  if (any(!is.finite(ev)) || any(ev <= proposal_ridge)) {
    return(
      list(
        available = FALSE, mean = opt$par, cov = diag(1, nt),
        convergence = opt$convergence
      )
    )
  }
  
  list(
    available = TRUE,
    mean = as.numeric(opt$par),
    cov = invert_positive_definite(HH, proposal_ridge),
    convergence = opt$convergence,
    objective = opt$value
  )
}

# Draw or evaluate the RJMCMC parameter-block proposal
proposal_log_density <- function(
    geom_entry, block, beta = NULL, tau2 = NULL, draw = FALSE,
    c_linear = 2.0, c0 = 316.227766, cplus = 316.227766,
    a_lambda = 1, b_lambda = 1,
    proposal_ridge = 1e-6,
    lambda_prior_mix = 0.10,
    beta_prior_mix = 0.15,
    proposal_link_rms = 1.0,
    lambda_opt_maxit = 25L,
    runtime = NULL) {
  
  if (!is.finite(lambda_prior_mix) ||
      lambda_prior_mix <= 0 ||
      lambda_prior_mix > 1) {
    stop("`lambda_prior_mix` must lie in (0,1].", call. = FALSE)
  }
  if (!is.finite(beta_prior_mix) || beta_prior_mix <= 0 || beta_prior_mix > 1) {
    stop("`beta_prior_mix` must lie in (0,1].", call. = FALSE)
  }
  if (!is.finite(proposal_link_rms) || proposal_link_rms <= 0) {
    stop(
      "`proposal_link_rms` must be finite and strictly positive.", call. = FALSE
    )
  }
  
  prop_geom <- geom_entry
  d_block <- ncol(block$X)
  if (!is.null(geom_entry$Gcond) && length(geom_entry$Gcond) == d_block &&
      all(is.finite(geom_entry$Gcond))) {
    prop_geom$Gperp <- as.numeric(geom_entry$Gcond)
  }
  Hc <- if (!is.null(geom_entry$Hcond)) as.matrix(geom_entry$Hcond) else NULL
  if (!is.null(Hc) &&
      all(dim(Hc) == c(d_block, d_block)) &&
      all(is.finite(Hc))) {
    prop_geom$Hbase <- Hc
  }
  
  lap <- list(available = FALSE, mean = numeric(0), cov = matrix(0, 0, 0))
  lambda_component <- "none"
  
  if (isTRUE(block$smooth)) {
    nt <- max(1L, length(block$S))
    lap <- laplace_log_lambda(
      prop_geom, block,
      c_linear = c_linear, c0 = c0, cplus = cplus,
      a_lambda = a_lambda, b_lambda = b_lambda,
      proposal_ridge = proposal_ridge,
      maxit = lambda_opt_maxit
    )
    
    if (!isTRUE(lap$available) && !is.null(runtime) && isTRUE(draw)) {
      runtime$lambda_laplace_failures <- runtime$lambda_laplace_failures + 1L
    }
    
    if (isTRUE(draw)) {
      use_lap <- isTRUE(lap$available) && stats::runif(1) > lambda_prior_mix
      if (use_lap) {
        z <- draw_multivariate_normal(lap$mean, lap$cov, proposal_ridge)
        lambda_component <- "laplace"
        if (!is.null(runtime)) {
          runtime$lambda_laplace_draws <- runtime$lambda_laplace_draws + 1L
        }
      } else {
        lam <- stats::rgamma(nt, shape = a_lambda, rate = b_lambda)
        if (any(!is.finite(lam)) || any(lam <= 0)) {
          return(list(valid = FALSE, beta = beta, tau2 = tau2, log_q = -Inf,
                      mean = numeric(0), cov = matrix(0, 0, 0),
                      lambda_component = "prior", beta_component = "none",
                      laplace_available = isTRUE(lap$available)))
        }
        z <- log(lam)
        lambda_component <- "prior"
        if (!is.null(runtime)) {
          runtime$lambda_prior_draws <- runtime$lambda_prior_draws + 1L
        }
      }
      tau2 <- exp(-z)
      if (any(!is.finite(tau2)) || any(tau2 <= 0)) {
        return(list(valid = FALSE, beta = beta, tau2 = tau2, log_q = -Inf,
                    mean = numeric(0), cov = matrix(0, 0, 0),
                    lambda_component = lambda_component,
                    beta_component = "none",
                    laplace_available = isTRUE(lap$available)))
      }
    } else {
      tau2 <- as.numeric(tau2)
      if (length(tau2) == 1L && nt > 1L) tau2 <- rep(tau2, nt)
      if (length(tau2) != nt || any(!is.finite(tau2)) || any(tau2 <= 0)) {
        return(list(valid = FALSE, beta = beta, tau2 = tau2, log_q = -Inf,
                    mean = numeric(0), cov = matrix(0, 0, 0),
                    lambda_component = "evaluation", beta_component = "none",
                    laplace_available = isTRUE(lap$available)))
      }
      z <- -log(tau2)
      lambda_component <- "evaluation"
    }
  } else {
    tau2 <- numeric(0)
    z <- numeric(0)
  }
  
  Q <- candidate_prior_precision(
    block, tau2 = tau2,
    c0 = c0, cplus = cplus, c_linear = c_linear,
    floor = proposal_ridge
  )
  H <- make_positive_definite(
    prop_geom$Hbase + Q + diag(proposal_ridge, ncol(block$X)),
    proposal_ridge
  )
  V_local_raw <- invert_positive_definite(H, proposal_ridge)
  m_local_raw <- as.numeric(V_local_raw %*% prop_geom$Gperp)
  
  local_std <- standardize_gaussian_on_predictor(
    m_local_raw, V_local_raw, block$X,
    target_rms = proposal_link_rms,
    floor = proposal_ridge
  )
  m_local <- local_std$mean
  V_local <- local_std$cov
  
  V_prior_raw <- invert_positive_definite(Q, proposal_ridge)
  m_prior <- rep(0, length(m_local))
  prior_std <- standardize_gaussian_on_predictor(
    m_prior, V_prior_raw, block$X,
    target_rms = proposal_link_rms,
    floor = proposal_ridge
  )
  V_prior <- prior_std$cov
  beta_component <- "evaluation"
  
  if (isTRUE(draw)) {
    use_local_beta <- stats::runif(1) > beta_prior_mix
    if (use_local_beta) {
      beta <- draw_multivariate_normal(m_local, V_local, proposal_ridge)
      beta_component <- "local"
      if (!is.null(runtime)) {
        runtime$beta_local_draws <- runtime$beta_local_draws + 1L
        if (local_std$scale < 1) {
          runtime$predictor_scaled_local_draws <-
            runtime$predictor_scaled_local_draws + 1L
        }
        runtime$predictor_scale_sum <-
          runtime$predictor_scale_sum + local_std$scale
        runtime$predictor_scale_min <- min(
          runtime$predictor_scale_min, local_std$scale
        )
      }
    } else {
      beta <- draw_multivariate_normal(m_prior, V_prior, proposal_ridge)
      beta_component <- "prior"
      if (!is.null(runtime)) {
        runtime$beta_prior_draws <- runtime$beta_prior_draws + 1L
        if (prior_std$scale < 1) {
          runtime$predictor_scaled_prior_draws <-
            runtime$predictor_scaled_prior_draws + 1L
        }
        runtime$predictor_scale_sum <-
          runtime$predictor_scale_sum + prior_std$scale
        runtime$predictor_scale_min <- min(
          runtime$predictor_scale_min, prior_std$scale
        )
      }
    }
  }
  
  if (is.null(beta) ||
      length(beta) != length(m_local) ||
      any(!is.finite(beta))) {
    return(list(valid = FALSE, beta = beta, tau2 = tau2, log_q = -Inf,
                mean = m_local, cov = V_local,
                lambda_component = lambda_component,
                beta_component = beta_component,
                laplace_available = isTRUE(lap$available)))
  }
  
  logq_beta_local <- log_multivariate_normal_density(
    beta, m_local, V_local, proposal_ridge
  )
  logq_beta_prior <- log_multivariate_normal_density(
    beta, m_prior, V_prior, proposal_ridge
  )
  logq_beta <- logsumexp_pair(
    log1p(-beta_prior_mix) + logq_beta_local,
    log(beta_prior_mix) + logq_beta_prior
  )
  logq <- logq_beta
  
  if (isTRUE(block$smooth)) {
    lambda <- exp(z)
    if (any(!is.finite(lambda)) || any(lambda <= 0)) {
      return(list(valid = FALSE, beta = beta, tau2 = tau2, log_q = -Inf,
                  mean = m_local, cov = V_local,
                  lambda_component = lambda_component,
                  beta_component = beta_component,
                  laplace_available = isTRUE(lap$available)))
    }
    
    logq_prior_z <- sum(
      stats::dgamma(lambda, shape = a_lambda, rate = b_lambda, log = TRUE) + z
    )
    
    if (isTRUE(lap$available)) {
      logq_lap_z <- log_multivariate_normal_density(
        z, lap$mean, lap$cov, proposal_ridge
      )
      logq_z <- logsumexp_pair(
        log1p(-lambda_prior_mix) + logq_lap_z,
        log(lambda_prior_mix) + logq_prior_z
      )
    } else {
      logq_z <- logq_prior_z
    }
    
    logq <- logq_beta + logq_z + sum(z)
  }
  
  list(
    valid = is.finite(logq),
    beta = as.numeric(beta),
    tau2 = as.numeric(tau2),
    log_q = as.numeric(logq),
    mean = m_local,
    cov = V_local,
    prior_cov = V_prior,
    local_predictor_scale = local_std$scale,
    local_predictor_rms_before = local_std$raw_rms,
    local_predictor_rms_after = local_std$scaled_rms,
    prior_predictor_scale = prior_std$scale,
    prior_predictor_rms_before = prior_std$raw_rms,
    prior_predictor_rms_after = prior_std$scaled_rms,
    lambda_component = lambda_component,
    beta_component = beta_component,
    geometry_mode = if (!is.null(geom_entry$Gcond) &&
                        !is.null(geom_entry$Hcond)) {
      "conditional_identity_map"
    } else {
      "profiled_geometry"
    },
    laplace_available = isTRUE(lap$available),
    lambda_mean = if (isTRUE(lap$available)) lap$mean else numeric(0),
    lambda_cov = if (isTRUE(lap$available)) lap$cov else matrix(0, 0, 0)
  )
}

# Assign a proposed candidate block to a parameter vector
assign_candidate_parameters <- function(theta, qdraw, block) {
  theta <- flatten_parameters(theta)
  if (is.null(names(theta)) || any(!nzchar(names(theta)))) return(NULL)
  
  beta_names <- block$beta_names
  tau_names <- block$tau_names
  if (length(beta_names) != ncol(block$X) ||
      length(qdraw$beta) != length(beta_names) ||
      any(!beta_names %in% names(theta)) ||
      any(!is.finite(qdraw$beta))) {
    return(NULL)
  }
  
  theta[beta_names] <- as.numeric(qdraw$beta)
  
  if (isTRUE(block$smooth)) {
    if (length(tau_names) != length(qdraw$tau2) ||
        any(!tau_names %in% names(theta)) ||
        (length(qdraw$tau2) &&
         (any(!is.finite(qdraw$tau2)) || any(qdraw$tau2 <= 0)))) {
      return(NULL)
    }
    if (length(tau_names)) theta[tau_names] <- as.numeric(qdraw$tau2)
  } else if (length(tau_names)) {
    return(NULL)
  }
  
  storage.mode(theta) <- "double"
  theta
}

# Extract candidate values for the reverse q_psi density
extract_candidate_parameters <- function(theta, block) {
  theta <- flatten_parameters(theta)
  if (is.null(names(theta))) return(NULL)
  
  beta_names <- block$beta_names
  tau_names <- block$tau_names
  if (length(beta_names) != ncol(block$X) ||
      any(!beta_names %in% names(theta)) ||
      any(!tau_names %in% names(theta))) {
    return(NULL)
  }
  
  beta <- as.numeric(theta[beta_names])
  tau2 <- if (isTRUE(block$smooth) && length(tau_names)) {
    as.numeric(theta[tau_names])
  } else numeric(0)
  
  if (any(!is.finite(beta)) ||
      (length(tau2) && (any(!is.finite(tau2)) || any(tau2 <= 0)))) {
    return(NULL)
  }
  
  list(beta = beta, tau2 = tau2)
}

# State-dependent birth/death/swap attempt probabilities
move_probabilities <- function(inclusion, admissibility, swap_weight = 1) {
  nb <- sum(birth_mask(inclusion, admissibility))
  nd <- sum(death_mask(inclusion, admissibility))
  w <- c(
    birth = as.numeric(nb > 0L),
    death = as.numeric(nd > 0L),
    swap = as.numeric(nb > 0L && nd > 0L) * max(0, as.numeric(swap_weight))
  )
  if (sum(w) <= 0) return(c(birth = 0, death = 0, swap = 0))
  w / sum(w)
}

# Exact MH-Green birth move
birth_move <- function(
    state, geometry,
    formula_base, data, fam_obj, admissibility, offset, runtime,
    a_V = 0.5, a_A = 0.2,
    c_linear = 2.0, c0 = 316.227766, cplus = 316.227766,
    a_lambda = 1, b_lambda = 1,
    proposal_ridge = 1e-6,
    lambda_prior_mix = 0.10,
    beta_prior_mix = 0.15,
    proposal_link_rms = 1.0,
    lambda_opt_maxit = 25L,
    swap_weight = 1) {
  
  if (is.null(geometry) || !length(geometry$q_b)) {
    return(list(state = state, accepted = FALSE, failed = FALSE))
  }
  pick <- sample(seq_along(geometry$q_b), 1L, prob = geometry$q_b)
  ge <- geometry$candidates[[pick]]
  r <- ge$row
  k <- ge$col
  q_b <- geometry$q_b[pick]
  
  plus_inc <- state$inclusion
  plus_inc[r, k] <- TRUE
  frame_error <- NULL
  frame_plus <- tryCatch(
    get_model_frame(
      plus_inc, formula_base, data, fam_obj, offset, runtime,
      c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda,
      c0 = c0, cplus = cplus
    ),
    error = function(e) {
      frame_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(frame_plus)) {
    cand_lab <- paste0(
      rownames(state$inclusion)[r], " -> ", colnames(state$inclusion)[k]
    )
    msg <- paste0(
      "birth frame failure [", cand_lab, "] from model ",
      state_key(state$inclusion), ": ",
      if (is.null(frame_error)) "<unknown error>" else frame_error
    )
    record_proposal_failure(runtime, msg)
    stop(
      "Frame construction failed for an admissible birth, which would remove ",
      "part of the proposal support: ", msg, ".",
      call. = FALSE
    )
  }
  
  template_plus <- parameter_template(frame_plus)
  theta_plus <- merge_parameters(template_plus, state$theta)
  candidate_key <- paste(r, k, sep = ":")
  block_plus <- frame_plus$candidate_metadata[[candidate_key]]
  if (is.null(block_plus)) {
    stop(
      "Birth frame is missing term-local metadata for ",
      rownames(state$inclusion)[r], " -> ", colnames(state$inclusion)[k], ".",
      call. = FALSE
    )
  }
  qdraw <- proposal_log_density(
    ge, block_plus, draw = TRUE,
    c_linear = c_linear, c0 = c0, cplus = cplus,
    a_lambda = a_lambda, b_lambda = b_lambda,
    proposal_ridge = proposal_ridge,
    lambda_prior_mix = lambda_prior_mix,
    beta_prior_mix = beta_prior_mix,
    proposal_link_rms = proposal_link_rms,
    lambda_opt_maxit = lambda_opt_maxit,
    runtime = runtime
  )
  theta_plus <- assign_candidate_parameters(theta_plus, qdraw, block_plus)
  if (is.null(theta_plus) || !qdraw$valid) {
    record_proposal_failure(
      runtime, "birth: coefficient/smoothing proposal was invalid"
    )
    return(list(state = state, accepted = FALSE, failed = TRUE))
  }
  
  prop <- list(inclusion = plus_inc, theta = theta_plus, logPost = -Inf)
  prop_comp <- posterior_components(
    prop, frame_plus, fam_obj, offset = frame_plus$model_offset,
    a_V = a_V, a_A = a_A, details = FALSE
  )
  prop$logPost <- prop_comp$total
  if (!is.finite(prop$logPost)) {
    fail_part <- if (!is.finite(prop_comp$logLik)) {
      "likelihood"
    } else if (!is.finite(prop_comp$structural_prior)) {
      "structural-prior"
    } else if (!is.finite(prop_comp$candidate_prior)) {
      "candidate-prior"
    } else if (!is.finite(prop_comp$model_prior)) {
      "model-prior"
    } else "unknown"
    cand_lab <- paste0(
      rownames(state$inclusion)[r], " -> ", colnames(state$inclusion)[k]
    )
    record_proposal_failure(
      runtime,
      paste0("birth non-finite [", cand_lab, "]: ", fail_part)
    )
    return(
      list(state = state, accepted = FALSE, failed = TRUE, row = r, col = k)
    )
  }
  
  rho_M <- move_probabilities(
    state$inclusion, admissibility, swap_weight = swap_weight
  )
  rho_plus <- move_probabilities(
    plus_inc, admissibility, swap_weight = swap_weight
  )
  death_plus <- death_mask(plus_inc, admissibility)
  q_d_plus <- 1 / sum(death_plus)
  log_alpha <- (prop$logPost - state$logPost) +
    log(rho_plus["death"]) + log(q_d_plus) -
    log(rho_M["birth"]) - log(q_b) - qdraw$log_q
  accept <- is.finite(log_alpha) && log(stats::runif(1)) < min(0, log_alpha)
  list(
    state = if (accept) prop else state,
    accepted = accept,
    failed = FALSE,
    row = r, col = k,
    log_alpha = log_alpha,
    delta_logPost = prop$logPost - state$logPost,
    lambda_component = qdraw$lambda_component,
    laplace_available = qdraw$laplace_available
  )
}
# Exact MH-Green death move
death_move <- function(
    state, frame_current,
    formula_base, data, fam_obj, admissibility, offset, runtime,
    epsilon = 0.25, omega = 0.25,
    a_V = 0.5, a_A = 0.2,
    c_linear = 2.0, c0 = 316.227766, cplus = 316.227766,
    a_lambda = 1, b_lambda = 1,
    proposal_ridge = 1e-6,
    full_cross_information = TRUE,
    fd_eps = 1e-4,
    lambda_prior_mix = 0.10,
    beta_prior_mix = 0.15,
    proposal_link_rms = 1.0,
    lambda_opt_maxit = 25L,
    swap_weight = 1) {
  
  active <- which(death_mask(state$inclusion, admissibility), arr.ind = TRUE)
  if (!nrow(active)) {
    return(list(state = state, accepted = FALSE, failed = FALSE))
  }
  pick <- sample(seq_len(nrow(active)), 1L)
  r <- active[pick, 1L]
  k <- active[pick, 2L]
  q_d <- 1 / nrow(active)
  
  minus_inc <- state$inclusion
  minus_inc[r, k] <- FALSE
  frame_minus <- tryCatch(get_model_frame(
    minus_inc, formula_base, data, fam_obj, offset, runtime,
    c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda,
    c0 = c0, cplus = cplus
  ), error = function(e) NULL)
  if (is.null(frame_minus)) {
    record_proposal_failure(
      runtime, "death: reduced-model frame construction failed"
    )
    return(list(state = state, accepted = FALSE, failed = TRUE))
  }
  
  template_minus <- parameter_template(frame_minus)
  theta_minus <- merge_parameters(template_minus, state$theta)
  prop <- list(inclusion = minus_inc, theta = theta_minus, logPost = -Inf)
  prop$logPost <- target_log_posterior(
    prop, frame_minus, fam_obj, offset = frame_minus$model_offset, a_V = a_V,
    a_A = a_A
  )
  if (!is.finite(prop$logPost)) {
    record_proposal_failure(
      runtime, "death: reduced model had non-finite exact posterior"
    )
    return(list(state = state, accepted = FALSE, failed = TRUE))
  }
  
  reverse_geom_error <- NULL
  reverse_geom <- tryCatch(
    evaluate_proposal_geometry(
      prop, frame_minus, data, fam_obj, formula_base, admissibility, offset,
      runtime,
      epsilon = epsilon, omega = omega,
      a_V = a_V, a_A = a_A,
      c_linear = c_linear, c0 = c0, cplus = cplus,
      a_lambda = a_lambda, b_lambda = b_lambda,
      proposal_ridge = proposal_ridge,
      full_cross_information = full_cross_information,
      fd_eps = fd_eps
    ),
    error = function(e) {
      reverse_geom_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(reverse_geom) && !is.null(reverse_geom_error)) {
    record_proposal_failure(
      runtime, paste0("death: reverse geometry failed: ", reverse_geom_error)
    )
    return(
      list(state = state, accepted = FALSE, failed = TRUE, row = r, col = k)
    )
  }
  target_key <- paste(r, k, sep = ":")
  q_b_reverse <- if (!is.null(reverse_geom) &&
                     target_key %in% names(reverse_geom$q_b)) {
    reverse_geom$q_b[[target_key]]
  } else NA_real_
  ge <- if (!is.null(reverse_geom)) {
    reverse_geom$candidates[[target_key]]
  } else {
    NULL
  }
  if (!is.finite(q_b_reverse) || q_b_reverse <= 0 || is.null(ge)) {
    record_proposal_failure(
      runtime, "death: reverse birth candidate/probability unavailable"
    )
    return(list(state = state, accepted = FALSE, failed = TRUE))
  }
  
  active_block <- frame_current$candidate_metadata[[target_key]]
  if (is.null(active_block)) {
    stop(
      "Death move is missing term-local metadata for active candidate ",
      rownames(state$inclusion)[r], " -> ", colnames(state$inclusion)[k], ".",
      call. = FALSE
    )
  }
  vals <- extract_candidate_parameters(state$theta, active_block)
  if (is.null(vals)) {
    record_proposal_failure(
      runtime, "death: active candidate parameters could not be extracted"
    )
    return(list(state = state, accepted = FALSE, failed = TRUE))
  }
  qrev <- proposal_log_density(
    ge, ge$block, beta = vals$beta, tau2 = vals$tau2, draw = FALSE,
    c_linear = c_linear, c0 = c0, cplus = cplus,
    a_lambda = a_lambda, b_lambda = b_lambda,
    proposal_ridge = proposal_ridge,
    lambda_prior_mix = lambda_prior_mix,
    beta_prior_mix = beta_prior_mix,
    proposal_link_rms = proposal_link_rms,
    lambda_opt_maxit = lambda_opt_maxit,
    runtime = NULL
  )
  if (!qrev$valid) {
    record_proposal_failure(
      runtime,
      "death: reverse coefficient/smoothing proposal density was non-finite"
    )
    return(list(state = state, accepted = FALSE, failed = TRUE))
  }
  
  rho_M <- move_probabilities(
    state$inclusion, admissibility, swap_weight = swap_weight
  )
  rho_minus <- move_probabilities(
    minus_inc, admissibility, swap_weight = swap_weight
  )
  log_alpha <- (prop$logPost - state$logPost) +
    log(rho_minus["birth"]) + log(q_b_reverse) + qrev$log_q -
    log(rho_M["death"]) - log(q_d)
  accept <- is.finite(log_alpha) && log(stats::runif(1)) < min(0, log_alpha)
  list(
    state = if (accept) prop else state,
    accepted = accept,
    failed = FALSE,
    row = r, col = k,
    log_alpha = log_alpha,
    delta_logPost = prop$logPost - state$logPost,
    lambda_component = qrev$lambda_component,
    laplace_available = qrev$laplace_available
  )
}

# Exact MH-Green swap move
swap_move <- function(
    state, frame_current,
    formula_base, data, fam_obj, admissibility, offset, runtime,
    epsilon = 0.25, omega = 0.25,
    a_V = 0.5, a_A = 0.2,
    c_linear = 2.0, c0 = 316.227766, cplus = 316.227766,
    a_lambda = 1, b_lambda = 1,
    proposal_ridge = 1e-6,
    full_cross_information = TRUE,
    fd_eps = 1e-4,
    lambda_prior_mix = 0.10,
    beta_prior_mix = 0.15,
    proposal_link_rms = 1.0,
    lambda_opt_maxit = 25L,
    swap_weight = 1) {
  
  active <- which(death_mask(state$inclusion, admissibility), arr.ind = TRUE)
  if (!nrow(active) || !sum(birth_mask(state$inclusion, admissibility))) {
    return(list(state = state, accepted = FALSE, failed = FALSE))
  }
  
  pick_old <- sample(seq_len(nrow(active)), 1L)
  r_old <- active[pick_old, 1L]
  k_old <- active[pick_old, 2L]
  old_key <- paste(r_old, k_old, sep = ":")
  q_remove_f <- 1 / nrow(active)
  
  common_inc <- state$inclusion
  common_inc[r_old, k_old] <- FALSE
  
  frame_common <- tryCatch(
    get_model_frame(
      common_inc, formula_base, data, fam_obj, offset, runtime,
      c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda,
      c0 = c0, cplus = cplus
    ),
    error = function(e) NULL
  )
  if (is.null(frame_common)) {
    record_proposal_failure(
      runtime, "swap: common-model frame construction failed"
    )
    return(list(state = state, accepted = FALSE, failed = TRUE,
                row_out = r_old, col_out = k_old))
  }
  
  template_common <- parameter_template(frame_common)
  theta_common <- merge_parameters(template_common, state$theta)
  common_state <- list(
    inclusion = common_inc, theta = theta_common, logPost = NA_real_
  )
  
  geom_common <- tryCatch(
    evaluate_proposal_geometry(
      common_state, frame_common, data, fam_obj, formula_base, admissibility,
      offset, runtime,
      epsilon = epsilon, omega = omega,
      a_V = a_V, a_A = a_A,
      c_linear = c_linear, c0 = c0, cplus = cplus,
      a_lambda = a_lambda, b_lambda = b_lambda,
      proposal_ridge = proposal_ridge,
      full_cross_information = full_cross_information,
      fd_eps = fd_eps
    ),
    error = function(e) NULL
  )
  if (is.null(geom_common) || !length(geom_common$q_b) ||
      !old_key %in% names(geom_common$q_b)) {
    record_proposal_failure(
      runtime, "swap: common-model informed geometry was unavailable"
    )
    return(list(state = state, accepted = FALSE, failed = TRUE,
                row_out = r_old, col_out = k_old))
  }
  
  add_keys_f <- setdiff(names(geom_common$q_b), old_key)
  if (!length(add_keys_f)) {
    return(list(state = state, accepted = FALSE, failed = FALSE,
                row_out = r_old, col_out = k_old))
  }
  q_add_f_vec <- geom_common$q_b[add_keys_f]
  q_add_f_vec <- q_add_f_vec / sum(q_add_f_vec)
  new_key <- sample(add_keys_f, 1L, prob = q_add_f_vec)
  q_add_f <- unname(q_add_f_vec[[new_key]])
  ge_new <- geom_common$candidates[[new_key]]
  if (is.null(ge_new) || !is.finite(q_add_f) || q_add_f <= 0) {
    record_proposal_failure(
      runtime, "swap: forward add candidate/probability was invalid"
    )
    return(list(state = state, accepted = FALSE, failed = TRUE,
                row_out = r_old, col_out = k_old))
  }
  
  r_new <- ge_new$row
  k_new <- ge_new$col
  swap_inc <- common_inc
  swap_inc[r_new, k_new] <- TRUE
  
  frame_swap <- tryCatch(
    get_model_frame(
      swap_inc, formula_base, data, fam_obj, offset, runtime,
      c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda,
      c0 = c0, cplus = cplus
    ),
    error = function(e) NULL
  )
  if (is.null(frame_swap)) {
    record_proposal_failure(
      runtime, "swap: proposed-model frame construction failed"
    )
    return(list(state = state, accepted = FALSE, failed = TRUE,
                row_out = r_old, col_out = k_old,
                row_in = r_new, col_in = k_new))
  }
  
  template_swap <- parameter_template(frame_swap)
  theta_swap <- merge_parameters(template_swap, theta_common)
  block_new <- frame_swap$candidate_metadata[[new_key]]
  if (is.null(block_new)) {
    stop(
      "Swap frame is missing term-local metadata for incoming candidate ",
      rownames(state$inclusion)[r_new], " -> ",
      colnames(state$inclusion)[k_new], ".",
      call. = FALSE
    )
  }
  
  q_new <- proposal_log_density(
    ge_new, block_new, draw = TRUE,
    c_linear = c_linear, c0 = c0, cplus = cplus,
    a_lambda = a_lambda, b_lambda = b_lambda,
    proposal_ridge = proposal_ridge,
    lambda_prior_mix = lambda_prior_mix,
    beta_prior_mix = beta_prior_mix,
    proposal_link_rms = proposal_link_rms,
    lambda_opt_maxit = lambda_opt_maxit,
    runtime = runtime
  )
  theta_swap <- assign_candidate_parameters(theta_swap, q_new, block_new)
  if (is.null(theta_swap) || !q_new$valid) {
    record_proposal_failure(
      runtime, "swap: forward coefficient/smoothing proposal was invalid"
    )
    return(list(state = state, accepted = FALSE, failed = TRUE,
                row_out = r_old, col_out = k_old,
                row_in = r_new, col_in = k_new))
  }
  
  prop <- list(inclusion = swap_inc, theta = theta_swap, logPost = -Inf)
  prop_comp <- posterior_components(
    prop, frame_swap, fam_obj, offset = frame_swap$model_offset,
    a_V = a_V, a_A = a_A, details = FALSE
  )
  prop$logPost <- prop_comp$total
  if (!is.finite(prop$logPost)) {
    fail_part <- if (!is.finite(prop_comp$logLik)) {
      "likelihood"
    } else if (!is.finite(prop_comp$structural_prior)) {
      "structural-prior"
    } else if (!is.finite(prop_comp$candidate_prior)) {
      "candidate-prior"
    } else if (!is.finite(prop_comp$model_prior)) {
      "model-prior"
    } else "unknown"
    in_lab <- paste0(
      rownames(state$inclusion)[r_new], " -> ", colnames(state$inclusion)[k_new]
    )
    out_lab <- paste0(
      rownames(state$inclusion)[r_old], " -> ", colnames(state$inclusion)[k_old]
    )
    record_proposal_failure(
      runtime,
      paste0(
        "swap non-finite [out=", out_lab, "; in=", in_lab, "]: ", fail_part
      )
    )
    return(list(state = state, accepted = FALSE, failed = TRUE,
                row_out = r_old, col_out = k_old,
                row_in = r_new, col_in = k_new))
  }
  
  ge_old <- geom_common$candidates[[old_key]]
  active_block_old <- frame_current$candidate_metadata[[old_key]]
  if (is.null(active_block_old)) {
    stop(
      "Swap move is missing term-local metadata for outgoing candidate ",
      rownames(state$inclusion)[r_old], " -> ",
      colnames(state$inclusion)[k_old], ".",
      call. = FALSE
    )
  }
  old_vals <- extract_candidate_parameters(state$theta, active_block_old)
  if (is.null(old_vals)) {
    record_proposal_failure(
      runtime, "swap: removed candidate parameters could not be extracted"
    )
    return(list(state = state, accepted = FALSE, failed = TRUE,
                row_out = r_old, col_out = k_old,
                row_in = r_new, col_in = k_new))
  }
  
  q_old <- proposal_log_density(
    ge_old, ge_old$block,
    beta = old_vals$beta, tau2 = old_vals$tau2, draw = FALSE,
    c_linear = c_linear, c0 = c0, cplus = cplus,
    a_lambda = a_lambda, b_lambda = b_lambda,
    proposal_ridge = proposal_ridge,
    lambda_prior_mix = lambda_prior_mix,
    beta_prior_mix = beta_prior_mix,
    proposal_link_rms = proposal_link_rms,
    lambda_opt_maxit = lambda_opt_maxit,
    runtime = NULL
  )
  if (!q_old$valid) {
    record_proposal_failure(
      runtime,
      "swap: reverse coefficient/smoothing proposal density was non-finite"
    )
    return(list(state = state, accepted = FALSE, failed = TRUE,
                row_out = r_old, col_out = k_old,
                row_in = r_new, col_in = k_new))
  }
  
  add_keys_r <- setdiff(names(geom_common$q_b), new_key)
  q_add_r_vec <- geom_common$q_b[add_keys_r]
  q_add_r_vec <- q_add_r_vec / sum(q_add_r_vec)
  q_add_r <- if (old_key %in% names(q_add_r_vec)) {
    unname(q_add_r_vec[[old_key]])
  } else {
    NA_real_
  }
  if (!is.finite(q_add_r) || q_add_r <= 0) {
    record_proposal_failure(
      runtime, "swap: reverse add probability was invalid"
    )
    return(list(state = state, accepted = FALSE, failed = TRUE,
                row_out = r_old, col_out = k_old,
                row_in = r_new, col_in = k_new))
  }
  
  q_remove_r <- 1 / sum(death_mask(swap_inc, admissibility))
  rho_M <- move_probabilities(
    state$inclusion, admissibility, swap_weight = swap_weight
  )
  rho_swap <- move_probabilities(
    swap_inc, admissibility, swap_weight = swap_weight
  )
  
  log_alpha <- (prop$logPost - state$logPost) +
    log(rho_swap["swap"]) + log(q_remove_r) + log(q_add_r) + q_old$log_q -
    log(rho_M["swap"]) - log(q_remove_f) - log(q_add_f) - q_new$log_q
  
  accept <- is.finite(log_alpha) && log(stats::runif(1)) < min(0, log_alpha)
  
  list(
    state = if (accept) prop else state,
    accepted = accept,
    failed = FALSE,
    row_out = r_old, col_out = k_old,
    row_in = r_new, col_in = k_new,
    log_alpha = log_alpha,
    delta_logPost = prop$logPost - state$logPost,
    q_forward = q_remove_f * q_add_f,
    q_reverse = q_remove_r * q_add_r,
    lambda_component = q_new$lambda_component,
    laplace_available = q_new$laplace_available
  )
}

# Effective sample size for a scalar RJMCMC statistic
effective_sample_size <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3L || stats::var(x) <= .Machine$double.eps) return(NA_real_)
  
  if (requireNamespace("coda", quietly = TRUE)) {
    z <- tryCatch(
      as.numeric(coda::effectiveSize(coda::mcmc(x)))[1L],
      error = function(e) NA_real_
    )
    if (is.finite(z) && z > 0) return(min(z, n))
  }
  
  ac <- tryCatch(
    as.numeric(
      stats::acf(x, plot = FALSE, lag.max = min(n - 1L, floor(n / 2L)))$acf
    )[-1L],
    error = function(e) numeric(0)
  )
  if (!length(ac)) return(NA_real_)
  
  s <- 0
  j <- 1L
  while (j <= length(ac)) {
    pair <- ac[j]
    if (j + 1L <= length(ac)) pair <- pair + ac[j + 1L]
    if (!is.finite(pair) || pair <= 0) break
    s <- s + pair
    j <- j + 2L
  }
  tau <- max(1, 1 + 2 * s)
  min(n, n / tau)
}

# RJMCMC model-space diagnostics
model_trace_diagnostics <- function(
    trace_inclusion, trace_model_size, elapsed_sec, trace_logpost = NULL) {
  if (!length(trace_inclusion)) return(list())
  
  keys <- vapply(trace_inclusion, state_key, character(1))
  occ <- sort(table(keys) / length(keys), decreasing = TRUE)
  top_occ <- data.frame(
    model_key = names(head(occ, 10L)),
    posterior_occupancy = as.numeric(head(occ, 10L)),
    stringsAsFactors = FALSE
  )
  
  inc_mat <- do.call(rbind, lapply(trace_inclusion, function(z) as.numeric(z)))
  first <- trace_inclusion[[1L]]
  inc_names <- as.vector(
    outer(rownames(first), colnames(first), paste, sep = "::")
  )
  if (length(inc_names) == ncol(inc_mat)) colnames(inc_mat) <- inc_names
  
  ess_inc <- vapply(
    seq_len(ncol(inc_mat)), function(j) effective_sample_size(inc_mat[, j]),
    numeric(1)
  )
  names(ess_inc) <- colnames(inc_mat)
  pip_vec <- colMeans(inc_mat)
  pip_mcse <- sqrt(pip_vec * (1 - pip_vec) / ess_inc)
  pip_mcse[!is.finite(pip_mcse)] <- NA_real_
  
  sz <- as.integer(trace_model_size)
  ess_size <- effective_sample_size(sz)
  
  if (length(sz) >= 2L) {
    lev <- sort(unique(sz))
    tr <- table(
      factor(head(sz, -1L), levels = lev),
      factor(tail(sz, -1L), levels = lev)
    )
    P <- matrix(0, nrow(tr), ncol(tr), dimnames = dimnames(tr))
    rs <- rowSums(tr)
    positive_rows <- rs > 0
    if (any(positive_rows)) {
      P[positive_rows, ] <-
        tr[positive_rows, , drop = FALSE] / rs[positive_rows]
    }
    if (any(!positive_rows)) {
      P[!positive_rows, !positive_rows] <- diag(sum(!positive_rows))
    }
    ev <- tryCatch(
      Mod(eigen(t(P), only.values = TRUE)$values),
      error = function(e) numeric(0)
    )
    ev <- sort(ev[is.finite(ev)], decreasing = TRUE)
    second_eigen <- if (length(ev) >= 2L) ev[2L] else NA_real_
  } else {
    P <- matrix(1, 1, 1, dimnames = list(as.character(sz), as.character(sz)))
    second_eigen <- NA_real_
  }
  
  ess_per_second <- NA_real_
  if (is.finite(ess_size) && elapsed_sec > 0) {
    ess_per_second <- ess_size / elapsed_sec
  }
  
  postburn_switches <- if (length(keys) >= 2L) {
    sum(tail(keys, -1L) != head(keys, -1L))
  } else 0L
  
  split_drift <- function(x, scale = c("sd", "iqr")) {
    scale <- match.arg(scale)
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 20L) return(NA_real_)
    h <- floor(n / 2L)
    a <- x[seq_len(h)]
    b <- x[(h + 1L):n]
    denom <- if (scale == "sd") stats::sd(x) else stats::IQR(x)
    if (!is.finite(denom) || denom <= sqrt(.Machine$double.eps)) {
      return(if (isTRUE(all.equal(mean(a), mean(b)))) 0 else Inf)
    }
    loc_a <- if (scale == "sd") mean(a) else stats::median(a)
    loc_b <- if (scale == "sd") mean(b) else stats::median(b)
    abs(loc_a - loc_b) / denom
  }
  
  list(
    retained_model_draws = length(trace_inclusion),
    unique_models_visited_postburn = length(unique(keys)),
    postburn_model_switches = as.integer(postburn_switches),
    cumulative_unique_models = cumsum(!duplicated(keys)),
    top_model_occupancy = top_occ,
    ess_model_size = ess_size,
    ess_model_size_per_second = ess_per_second,
    ess_inclusion = ess_inc,
    pip_mcse = pip_mcse,
    model_size_transition_matrix = P,
    model_size_second_eigenvalue = second_eigen,
    model_size_split_drift = split_drift(sz, "sd"),
    logpost_split_drift = split_drift(trace_logpost, "iqr")
  )
}

# Assess RJMCMC kernel operation and posterior model-space mixing
assess_model_space_mixing <- function(
    model_diag, association_pip, nonlinear_pip,
    association_diag, nonlinear_diag,
    attempted, accepted, failed, runtime, visited_models_all,
    min_retained = 1000L, min_switches = 50L,
    min_effective_size = 50, max_pip_mcse = 0.05,
    max_split_drift = 0.50, max_second_eigenvalue = 0.995,
    min_within_success = 0.90, max_proposal_failure_rate = 0.02) {
  
  attempted_total <- sum(attempted)
  failure_count <- max(sum(failed), runtime$proposal_failures)
  failure_rate <- failure_count / max(1L, attempted_total)
  within_rate <- runtime$within_successes / max(1L, runtime$within_attempts)
  reverse_accepts <- unname(accepted["death"]) + unname(accepted["swap"])
  
  linear_ess_attempts <- sum(runtime$within_linear_ess_attempts)
  linear_ess_successes <- sum(runtime$within_linear_ess_successes)
  linear_ess_failures <- sum(runtime$within_linear_ess_failures)
  linear_ess_moves <- sum(runtime$within_linear_ess_moves)
  
  linear_ess_operational <- linear_ess_attempts == 0L || (
    linear_ess_successes == linear_ess_attempts &&
      linear_ess_failures == 0L &&
      linear_ess_moves == linear_ess_successes
  )
  
  kernel_operational <-
    is.finite(within_rate) &&
    within_rate >= min_within_success &&
    is.finite(failure_rate) &&
    failure_rate <= max(0.05, max_proposal_failure_rate) &&
    isTRUE(linear_ess_operational) &&
    unname(accepted["birth"]) >= 1L &&
    reverse_accepts >= 1L &&
    visited_models_all >= 2L
  
  summarize_indicator <- function(prob, diag_obj) {
    if (is.null(prob) || !length(prob)) {
      return(list(n_uncertain = 0L, min_ess = NA_real_, max_mcse = NA_real_))
    }
    p <- as.numeric(prob)
    nm <- as.vector(outer(rownames(prob), colnames(prob), paste, sep = "::"))
    names(p) <- nm
    uncertain <- is.finite(p) & p > 0.025 & p < 0.975
    if (!any(uncertain)) {
      return(list(n_uncertain = 0L, min_ess = NA_real_, max_mcse = 0))
    }
    keys <- names(p)[uncertain]
    ee <- unname(diag_obj$ess[keys])
    mm <- unname(diag_obj$mcse[keys])
    ee <- ee[is.finite(ee)]
    mm <- mm[is.finite(mm)]
    list(
      n_uncertain = length(keys),
      min_ess = if (length(ee)) min(ee) else NA_real_,
      max_mcse = if (length(mm)) max(mm) else Inf
    )
  }
  ainfo <- summarize_indicator(association_pip, association_diag)
  ninfo <- summarize_indicator(nonlinear_pip, nonlinear_diag)
  
  retained <- model_diag$retained_model_draws
  switch_req <- max(
    as.integer(min_switches), as.integer(ceiling(0.01 * retained))
  )
  ess_candidates <- c(model_diag$ess_model_size, ainfo$min_ess, ninfo$min_ess)
  ess_candidates <- ess_candidates[is.finite(ess_candidates)]
  exploration_ess <- if (length(ess_candidates)) {
    min(ess_candidates)
  } else {
    NA_real_
  }
  
  diagnostic_indicators <- c(
    kernel_operational = kernel_operational,
    retained_draws = is.finite(retained) && retained >= min_retained,
    model_switches = is.finite(model_diag$postburn_model_switches) &&
      model_diag$postburn_model_switches >= switch_req,
    effective_size = is.finite(exploration_ess) &&
      exploration_ess >= min_effective_size,
    association_mcse = is.finite(ainfo$max_mcse) &&
      ainfo$max_mcse <= max_pip_mcse,
    nonlinear_mcse = is.finite(ninfo$max_mcse) &&
      ninfo$max_mcse <= max_pip_mcse,
    model_size_stability = is.finite(model_diag$model_size_split_drift) &&
      model_diag$model_size_split_drift <= max_split_drift,
    logpost_stability = is.finite(model_diag$logpost_split_drift) &&
      model_diag$logpost_split_drift <= max_split_drift,
    transition_spectrum = !is.finite(model_diag$model_size_second_eigenvalue) ||
      model_diag$model_size_second_eigenvalue <= max_second_eigenvalue,
    proposal_failures = is.finite(failure_rate) &&
      failure_rate <= max_proposal_failure_rate,
    within_model_updates = is.finite(within_rate) &&
      within_rate >= min_within_success,
    linear_ess_updates = isTRUE(linear_ess_operational)
  )
  
  reasons <- names(diagnostic_indicators)[!diagnostic_indicators]
  kernel_reasons <- character(0)
  if (!is.finite(within_rate) || within_rate < min_within_success) {
    kernel_reasons <- c(kernel_reasons, "within_model_updates")
  }
  if (!is.finite(failure_rate) ||
      failure_rate > max(0.05, max_proposal_failure_rate)) {
    kernel_reasons <- c(kernel_reasons, "proposal_failures")
  }
  if (!isTRUE(linear_ess_operational)) {
    kernel_reasons <- c(kernel_reasons, "linear_ess_updates")
  }
  if (unname(accepted["birth"]) < 1L) {
    kernel_reasons <- c(kernel_reasons, "no_accepted_birth")
  }
  if (reverse_accepts < 1L) {
    kernel_reasons <- c(kernel_reasons, "no_accepted_death_or_swap")
  }
  if (visited_models_all < 2L) {
    kernel_reasons <- c(kernel_reasons, "fewer_than_two_models_visited")
  }
  
  list(
    kernel_operational = isTRUE(kernel_operational),
    mixing_within_reference = isTRUE(all(diagnostic_indicators)),
    diagnostic_indicators = diagnostic_indicators,
    diagnostics_outside_reference = reasons,
    kernel_diagnostics_outside_reference = unique(kernel_reasons),
    retained_draws_reference = as.integer(min_retained),
    model_switches_reference = switch_req,
    effective_size_reference = min_effective_size,
    pip_mcse_reference = max_pip_mcse,
    split_drift_reference = max_split_drift,
    second_eigenvalue_reference = max_second_eigenvalue,
    proposal_failure_rate = failure_rate,
    within_model_success_rate = within_rate,
    uncertain_association_count = ainfo$n_uncertain,
    uncertain_nonlinear_count = ninfo$n_uncertain,
    min_effective_size_observed = exploration_ess,
    max_association_pip_mcse = ainfo$max_mcse,
    max_nonlinear_pip_mcse = ninfo$max_mcse
  )
}

# Bayesian spatial model with dynamic covariate allocation
# Search snapshots: written through a temporary file, so an interrupted run
# leaves no partial file.
write_search_snapshot <- function(path, key, n_iter, iteration, state,
                                  runtime, elapsed) {
  snap <- list(
    key = key, n_iter = n_iter, iteration = iteration, state = state,
    runtime = runtime, elapsed = elapsed,
    seed = get0(".Random.seed", globalenv(), inherits = FALSE)
  )
  tmp <- paste0(path, ".tmp")
  saveRDS(snap, tmp, compress = FALSE)
  file.rename(tmp, path)
  invisible(path)
}

read_search_snapshot <- function(path, key, n_iter, label) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  snap <- tryCatch(readRDS(path), error = function(e) NULL)
  usable <- !is.null(snap) && identical(snap$key, key) &&
    identical(snap$n_iter, n_iter) && snap$iteration < n_iter
  if (!usable) {
    message("  [Resume] ", basename(path), " does not match; search restarts.")
    return(NULL)
  }
  message(
    "  [Resume] ", label, ": continuing the model search from iteration ",
    snap$iteration + 1L, " of ", n_iter, "."
  )
  snap
}

fit_spatial_model_dynamic <- function(
    label, formula_base, candidate_covariates, family_name, data, offset = NULL,
    candidate_admissibility = NULL,
    n_iter = 10000, n_burn = 4000, thin = 5, model_thin = 1L,
    optimizer_maxit = 400L,
    optimizer_eps = .Machine$double.eps^0.25,
    optimizer_nu = TRUE,
    within_sweeps = 5L,
    a_V = 0.5, a_A = 0.2,
    c_linear = 1.5,
    c0 = 316.227766, cplus = 316.227766,
    a_lambda = 1, b_lambda = 1,
    epsilon_birth = 0.25,
    omega_birth = 0.25,
    proposal_ridge = 1e-6,
    full_cross_information = TRUE,
    fd_eps = 1e-4,
    swap_weight = 1,
    lambda_prior_mix = 0.10,
    beta_prior_mix = 0.15,
    proposal_link_rms = 1.0,
    lambda_opt_maxit = 25L,
    cache_frames = TRUE,
    cache_blocks = TRUE,
    max_frame_cache = Inf,
    max_representative_states = Inf,
    store_theta = TRUE,
    resume_file = NULL, resume_key = NULL, resume_every = 500L) {
  
  n_iter <- as.integer(n_iter)
  n_burn <- as.integer(n_burn)
  thin <- as.integer(thin)
  model_thin <- as.integer(model_thin)
  within_sweeps <- as.integer(within_sweeps)
  
  data <- prepare_candidate_data(data, candidate_covariates)
  fam_obj <- get_model_family(family_name)
  check_analysis_sample(
    data = data,
    formula_base = formula_base,
    candidate_covariates = candidate_covariates,
    offset = offset
  )
  original_candidate_covariates <- unique(
    candidate_map(candidate_covariates, data)$variable
  )
  candidate_covariates <- expand_candidate_covariates(
    candidate_covariates, data
  )
  params <- fam_obj$names
  inclusion <- matrix(
    FALSE,
    nrow = length(candidate_covariates), ncol = length(params),
    dimnames = list(candidate_covariates, params)
  )
  admissibility <- build_admissibility_matrix(
    formula_base, candidate_covariates, params, data,
    external_admissibility = candidate_admissibility
  )
  runtime <- new_runtime_state(
    cache_frames = cache_frames,
    cache_blocks = cache_blocks,
    max_frame_cache = max_frame_cache
  )
  
  max_representative_states <- as.numeric(max_representative_states)[1L]
  if (!is.finite(max_representative_states) &&
      !is.infinite(max_representative_states)) {
    stop(
      "`max_representative_states` must be a positive number or Inf.",
      call. = FALSE
    )
  }
  if (is.finite(max_representative_states) && max_representative_states < 1) {
    stop("`max_representative_states` must be at least 1.", call. = FALSE)
  }
  
  prior_summary <- model_prior_summary(
    length(original_candidate_covariates), ncol(inclusion),
    a_V = a_V, a_A = a_A,
    admissibility = admissibility
  )
  
  frame0 <- get_model_frame(
    inclusion, formula_base, data, fam_obj, offset, runtime,
    c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda,
    c0 = c0, cplus = cplus
  )
  
  if (!is.null(frame0$family)) fam_obj <- frame0$family
  fam_label <- if (is.character(family_name) && length(family_name) == 1L) {
    family_name
  } else if (!is.null(fam_obj$family)) {
    as.character(fam_obj$family)[1L]
  } else "<custom>"
  fam_obj <- validate_model_family(fam_obj, fam_label)
  
  current <- initialize_model_state(
    inclusion, frame0, fam_obj, runtime,
    optimizer_maxit = optimizer_maxit,
    optimizer_eps = optimizer_eps,
    optimizer_nu = optimizer_nu,
    a_V = a_V, a_A = a_A
  )
  
  birth_geometry_check <- validate_birth_geometry(
    state = current,
    frame = frame0,
    data = data,
    fam_obj = fam_obj,
    formula_base = formula_base,
    admissibility = admissibility,
    offset = offset,
    runtime = runtime,
    c_linear = c_linear,
    c0 = c0,
    cplus = cplus,
    a_lambda = a_lambda,
    b_lambda = b_lambda,
    proposal_ridge = proposal_ridge,
    full_cross_information = full_cross_information,
    fd_eps = fd_eps,
    epsilon = epsilon_birth,
    omega = omega_birth,
    a_V = a_V,
    a_A = a_A
  )
  
  model_keep_iter <- seq_len(n_iter)
  model_keep_iter <- model_keep_iter[
    model_keep_iter > n_burn & ((model_keep_iter - n_burn) %% model_thin == 0L)
  ]
  theta_keep_iter <- seq_len(n_iter)
  theta_keep_iter <- theta_keep_iter[
    theta_keep_iter > n_burn & ((theta_keep_iter - n_burn) %% thin == 0L)
  ]
  
  trace_logpost <- rep(NA_real_, length(model_keep_iter))
  trace_model_size <- rep(NA_integer_, length(model_keep_iter))
  trace_inclusion <- vector("list", length(model_keep_iter))
  retained_states <- if (isTRUE(store_theta)) {
    vector("list", length(theta_keep_iter))
  } else {
    NULL
  }
  representative_states <- list()
  model_pos <- 0L
  theta_pos <- 0L
  
  move_names <- c("birth", "death", "swap")
  attempted <- stats::setNames(integer(3), move_names)
  accepted <- stats::setNames(integer(3), move_names)
  failed <- stats::setNames(integer(3), move_names)
  log_alpha_trace <- stats::setNames(vector("list", 3), move_names)
  delta_logpost_trace <- stats::setNames(vector("list", 3), move_names)
  
  cand_attempt_birth <- matrix(
    0L, nrow(inclusion), ncol(inclusion), dimnames = dimnames(inclusion)
  )
  cand_accept_birth <- cand_attempt_birth
  cand_attempt_death <- cand_attempt_birth
  cand_accept_death <- cand_attempt_birth
  cand_attempt_swap_in <- cand_attempt_birth
  cand_accept_swap_in <- cand_attempt_birth
  cand_attempt_swap_out <- cand_attempt_birth
  cand_accept_swap_out <- cand_attempt_birth
  
  visited_keys <- character(n_iter)
  model_switches <- 0L
  print_step <- max(1L, floor(n_iter / 10L))
  elapsed_before <- 0
  first_iter <- 1L
  
  # State that has to survive an interruption. The frame and block caches
  # are left out: they are rebuilt on demand.
  snapshot_names <- c(
    "current", "trace_logpost", "trace_model_size", "trace_inclusion",
    "retained_states", "representative_states", "model_pos", "theta_pos",
    "attempted", "accepted", "failed", "log_alpha_trace",
    "delta_logpost_trace", "cand_attempt_birth", "cand_accept_birth",
    "cand_attempt_death", "cand_accept_death", "cand_attempt_swap_in",
    "cand_accept_swap_in", "cand_attempt_swap_out", "cand_accept_swap_out",
    "visited_keys", "model_switches"
  )
  cache_names <- c("frame_cache", "block_cache", "frame_cache_order")
  
  snap <- read_search_snapshot(resume_file, resume_key, n_iter, label)
  if (!is.null(snap)) {
    for (nm in snapshot_names) {
      assign(nm, snap$state[[nm]], envir = environment())
    }
    for (nm in names(snap$runtime)) runtime[[nm]] <- snap$runtime[[nm]]
    runtime$frame_cache_order <- character(0)
    elapsed_before <- snap$elapsed
    first_iter <- snap$iteration + 1L
    if (!is.null(snap$seed)) assign(".Random.seed", snap$seed, globalenv())
  }
  t0 <- proc.time()[["elapsed"]] - elapsed_before
  # RBST_TRACE=1 prints the posterior components at every step
  trace_on <- nzchar(Sys.getenv("RBST_TRACE"))
  trace_step <- function(tag, st, fr) {
    pc <- posterior_components(
      st, fr, fam_obj, offset = fr$model_offset, a_V = a_V, a_A = a_A
    )
    th <- flatten_parameters(st$theta)
    big <- names(th)[which.max(abs(th))]
    message(sprintf(
      "[trace] %d %s ll=%.4g str=%.4g cand=%.4g |A|=%d max|theta|=%.3g %s",
      m, tag, pc$logLik, pc$structural_prior, pc$candidate_prior,
      sum(st$inclusion), max(abs(th)), big
    ))
  }
  
  for (m in seq(first_iter, n_iter)) {
    frame_current <- get_model_frame(
      current$inclusion, formula_base, data, fam_obj, offset, runtime,
      c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda,
      c0 = c0, cplus = cplus
    )
    
    if (trace_on) trace_step("start", current, frame_current)
    current <- within_model_update(
      current, frame_current, fam_obj, runtime,
      within_sweeps = within_sweeps,
      a_V = a_V, a_A = a_A
    )
    if (trace_on) trace_step("within", current, frame_current)
    
    key_before <- state_key(current$inclusion)
    rho <- move_probabilities(
      current$inclusion, admissibility, swap_weight = swap_weight
    )
    
    if (sum(rho) > 0) {
      move <- sample(names(rho), 1L, prob = rho)
      attempted[[move]] <- attempted[[move]] + 1L
      
      if (identical(move, "birth")) {
        geometry <- tryCatch(
          evaluate_proposal_geometry(
            current, frame_current, data, fam_obj, formula_base, admissibility,
            offset, runtime,
            epsilon = epsilon_birth, omega = omega_birth,
            a_V = a_V, a_A = a_A,
            c_linear = c_linear, c0 = c0, cplus = cplus,
            a_lambda = a_lambda, b_lambda = b_lambda,
            proposal_ridge = proposal_ridge,
            full_cross_information = full_cross_information,
            fd_eps = fd_eps
          ),
          error = function(e) stop(
            paste0("Birth-proposal geometry failed at iteration ", m),
            call. = FALSE
          )
        )
        
        out <- birth_move(
          current, geometry,
          formula_base, data, fam_obj, admissibility, offset, runtime,
          a_V = a_V, a_A = a_A,
          c_linear = c_linear, c0 = c0, cplus = cplus,
          a_lambda = a_lambda, b_lambda = b_lambda,
          proposal_ridge = proposal_ridge,
          lambda_prior_mix = lambda_prior_mix,
          beta_prior_mix = beta_prior_mix,
          proposal_link_rms = proposal_link_rms,
          lambda_opt_maxit = lambda_opt_maxit,
          swap_weight = swap_weight
        )
        
        if (!is.null(out$row) && !is.null(out$col)) {
          rc <- cbind(out$row, out$col)
          cand_attempt_birth[rc] <- cand_attempt_birth[rc] + 1L
          if (isTRUE(out$accepted)) {
            cand_accept_birth[rc] <- cand_accept_birth[rc] + 1L
          }
        }
        
      } else if (identical(move, "death")) {
        out <- death_move(
          current, frame_current,
          formula_base, data, fam_obj, admissibility, offset, runtime,
          epsilon = epsilon_birth, omega = omega_birth,
          a_V = a_V, a_A = a_A,
          c_linear = c_linear, c0 = c0, cplus = cplus,
          a_lambda = a_lambda, b_lambda = b_lambda,
          proposal_ridge = proposal_ridge,
          full_cross_information = full_cross_information,
          fd_eps = fd_eps,
          lambda_prior_mix = lambda_prior_mix,
          beta_prior_mix = beta_prior_mix,
          proposal_link_rms = proposal_link_rms,
          lambda_opt_maxit = lambda_opt_maxit,
          swap_weight = swap_weight
        )
        
        if (!is.null(out$row) && !is.null(out$col)) {
          rc <- cbind(out$row, out$col)
          cand_attempt_death[rc] <- cand_attempt_death[rc] + 1L
          if (isTRUE(out$accepted)) {
            cand_accept_death[rc] <- cand_accept_death[rc] + 1L
          }
        }
        
      } else {
        out <- swap_move(
          current, frame_current,
          formula_base, data, fam_obj, admissibility, offset, runtime,
          epsilon = epsilon_birth, omega = omega_birth,
          a_V = a_V, a_A = a_A,
          c_linear = c_linear, c0 = c0, cplus = cplus,
          a_lambda = a_lambda, b_lambda = b_lambda,
          proposal_ridge = proposal_ridge,
          full_cross_information = full_cross_information,
          fd_eps = fd_eps,
          lambda_prior_mix = lambda_prior_mix,
          beta_prior_mix = beta_prior_mix,
          proposal_link_rms = proposal_link_rms,
          lambda_opt_maxit = lambda_opt_maxit,
          swap_weight = swap_weight
        )
        
        if (!is.null(out$row_out) && !is.null(out$col_out)) {
          rc <- cbind(out$row_out, out$col_out)
          cand_attempt_swap_out[rc] <- cand_attempt_swap_out[rc] + 1L
          if (isTRUE(out$accepted)) {
            cand_accept_swap_out[rc] <- cand_accept_swap_out[rc] + 1L
          }
        }
        if (!is.null(out$row_in) && !is.null(out$col_in)) {
          rc <- cbind(out$row_in, out$col_in)
          cand_attempt_swap_in[rc] <- cand_attempt_swap_in[rc] + 1L
          if (isTRUE(out$accepted)) {
            cand_accept_swap_in[rc] <- cand_accept_swap_in[rc] + 1L
          }
        }
      }
      
      if (isTRUE(out$accepted)) accepted[[move]] <- accepted[[move]] + 1L
      if (isTRUE(out$failed)) failed[[move]] <- failed[[move]] + 1L
      if (!is.null(out$log_alpha) && is.finite(out$log_alpha)) {
        log_alpha_trace[[move]] <- c(log_alpha_trace[[move]], out$log_alpha)
      }
      if (!is.null(out$delta_logPost) && is.finite(out$delta_logPost)) {
        delta_logpost_trace[[move]] <- c(
          delta_logpost_trace[[move]], out$delta_logPost
        )
      }
      if (trace_on) {
        message(sprintf(
          "[trace] %d %s accepted=%s failed=%s delta=%.4g log_alpha=%.4g",
          m, move, isTRUE(out$accepted), isTRUE(out$failed),
          if (is.null(out$delta_logPost)) NA_real_ else out$delta_logPost,
          if (is.null(out$log_alpha)) NA_real_ else out$log_alpha
        ))
      }
      current <- out$state
    }
    
    key_after <- state_key(current$inclusion)
    if (!identical(key_before, key_after)) model_switches <- model_switches + 1L
    visited_keys[m] <- key_after
    
    if (m > n_burn && isTRUE(store_theta)) {
      candidate_rep <- list(
        inclusion = current$inclusion,
        theta = current$theta,
        logPost = current$logPost
      )
      old_rep <- representative_states[[key_after]]
      
      if (!is.null(old_rep)) {
        if (!is.finite(old_rep$logPost) || current$logPost > old_rep$logPost) {
          representative_states[[key_after]] <- candidate_rep
        }
      } else if (!is.finite(max_representative_states) ||
                 length(representative_states) < max_representative_states) {
        representative_states[[key_after]] <- candidate_rep
      } else if (length(representative_states)) {
        rep_logpost <- vapply(
          representative_states,
          function(z) {
            val <- suppressWarnings(as.numeric(z$logPost)[1L])
            if (length(val) && is.finite(val)) val else -Inf
          },
          numeric(1)
        )
        worst <- which.min(rep_logpost)
        if (current$logPost > rep_logpost[worst]) {
          representative_states[[worst]] <- NULL
          representative_states[[key_after]] <- candidate_rep
        }
      }
    }
    
    if (m > n_burn && ((m - n_burn) %% model_thin == 0L)) {
      model_pos <- model_pos + 1L
      trace_logpost[model_pos] <- current$logPost
      trace_model_size[model_pos] <-
        model_dimensions(current$inclusion)$associations
      trace_inclusion[[model_pos]] <- current$inclusion
    }
    
    if (isTRUE(store_theta) && m > n_burn && ((m - n_burn) %% thin == 0L)) {
      theta_pos <- theta_pos + 1L
      retained_states[[theta_pos]] <- list(
        inclusion = current$inclusion,
        theta = current$theta,
        logPost = current$logPost
      )
    }
    
    if (m == 1L || m %% print_step == 0L) {
      message(sprintf(
        paste0(
          "  -> Model search %s: %d/%d | |A|=%d | |NL|=%d | logPost=%.2f | ",
          "birth %.1f%% (fail %d) | death %.1f%% | swap %.1f%% | ",
          "frame builds %d / cache hits %d"
        ),
        label, m, n_iter,
        model_dimensions(current$inclusion)$associations,
        model_dimensions(current$inclusion)$nonlinear, current$logPost,
        100 * accepted["birth"] / max(1L, attempted["birth"]),
        failed["birth"],
        100 * accepted["death"] / max(1L, attempted["death"]),
        100 * accepted["swap"] / max(1L, attempted["swap"]),
        runtime$frame_misses, runtime$frame_hits
      ))
    }
    
    if (!is.null(resume_file) && m < n_iter && m %% resume_every == 0L) {
      write_search_snapshot(
        resume_file, resume_key, n_iter, m,
        state = mget(snapshot_names, envir = environment()),
        runtime = as.list(runtime)[
          setdiff(names(as.list(runtime)), cache_names)
        ],
        elapsed = proc.time()[["elapsed"]] - t0
      )
    }
  }
  
  elapsed_sec <- proc.time()[["elapsed"]] - t0
  # The caches hold one design per visited model; free them before the
  # post-selection chains are forked.
  for (cache in c(runtime$frame_cache, runtime$block_cache)) {
    rm(list = ls(cache, all.names = TRUE), envir = cache)
  }
  invisible(gc(verbose = FALSE))
  
  component_pip <- Reduce(
    "+", lapply(trace_inclusion, function(z) z * 1)
  ) / length(trace_inclusion)
  component_pip[!admissibility] <- NA_real_
  pip_summary <- aggregate_component_pips(component_pip)
  pip <- pip_summary$association
  nonlinear_pip <- pip_summary$nonlinear
  nonlinear_given_association <- pip_summary$nonlinear_given_association
  indicator_traces <- aggregate_trace_indicators(trace_inclusion)
  association_diag <- indicator_trace_diagnostics(indicator_traces$association)
  nonlinear_diag <- indicator_trace_diagnostics(indicator_traces$nonlinear)
  
  model_diag <- model_trace_diagnostics(
    trace_inclusion = trace_inclusion,
    trace_model_size = trace_model_size,
    elapsed_sec = elapsed_sec,
    trace_logpost = trace_logpost
  )
  mixing_assessment <- assess_model_space_mixing(
    model_diag = model_diag,
    association_pip = pip,
    nonlinear_pip = nonlinear_pip,
    association_diag = association_diag,
    nonlinear_diag = nonlinear_diag,
    attempted = attempted,
    accepted = accepted,
    failed = failed,
    runtime = runtime,
    visited_models_all = length(unique(visited_keys))
  )
  
  qfun <- function(z) {
    z <- z[is.finite(z)]
    if (!length(z)) {
      return(
        stats::setNames(rep(NA_real_, 5), c("q05", "q25", "q50", "q75", "q95"))
      )
    }
    stats::setNames(
      as.numeric(
        stats::quantile(z, probs = c(.05, .25, .5, .75, .95), names = FALSE)
      ),
      c("q05", "q25", "q50", "q75", "q95")
    )
  }
  
  diagnostics <- c(
    list(
      initial_optimizer_fits = runtime$initial_optimizer_fits,
      initial_optimizer_warning_count = runtime$initial_optimizer_warning_count,
      initial_optimizer_warning_messages = runtime_message_table(
        runtime$initial_optimizer_warning_messages
      ),
      unique_model_frames = runtime$frame_misses,
      frame_cache_hits = runtime$frame_hits,
      frame_cache_hit_rate =
        runtime$frame_hits / max(1, runtime$frame_hits + runtime$frame_misses),
      frame_cache_limit = runtime$max_frame_cache,
      frame_cache_peak = runtime$frame_cache_peak,
      frame_cache_evictions = runtime$frame_cache_evictions,
      retained_theta_states = if (is.null(retained_states)) {
        0L
      } else {
        length(retained_states)
      },
      representative_states_retained = length(representative_states),
      representative_state_limit = max_representative_states,
      unique_candidate_blocks = runtime$block_misses,
      block_cache_hits = runtime$block_hits,
      
      within_model_attempts = runtime$within_attempts,
      within_model_successes = runtime$within_successes,
      within_model_failures = runtime$within_failures,
      within_model_success_rate =
        runtime$within_successes / max(1L, runtime$within_attempts),
      within_model_error_messages =
        runtime_message_table(runtime$within_error_messages),
      within_linear_ess_attempts = runtime$within_linear_ess_attempts,
      within_linear_ess_successes = runtime$within_linear_ess_successes,
      within_linear_ess_failures = runtime$within_linear_ess_failures,
      within_linear_ess_success_rate = runtime_stat_ratio(
        runtime$within_linear_ess_successes, runtime$within_linear_ess_attempts
      ),
      within_linear_ess_moves = runtime$within_linear_ess_moves,
      within_linear_ess_move_rate = runtime_stat_ratio(
        runtime$within_linear_ess_moves, runtime$within_linear_ess_successes
      ),
      within_linear_ess_coefficients_moved =
        runtime$within_linear_ess_coefficients_moved,
      within_linear_ess_max_abs_change =
        runtime$within_linear_ess_max_abs_change,
      within_linear_ess_contractions = runtime$within_linear_ess_contractions,
      within_linear_ess_mean_contractions = runtime_stat_ratio(
        runtime$within_linear_ess_contractions,
        runtime$within_linear_ess_successes
      ),
      within_linear_ess_max_contractions =
        runtime$within_linear_ess_max_contractions,
      within_linear_ess_error_messages = runtime_message_table(
        runtime$within_linear_ess_error_messages
      ),
      
      proposal_failures = runtime$proposal_failures,
      proposal_error_messages =
        runtime_message_table(runtime$proposal_error_messages),
      elapsed_seconds = elapsed_sec,
      visited_models_all_iterations = length(unique(visited_keys)),
      accepted_model_switches = model_switches,
      move_attempted = attempted,
      move_accepted = accepted,
      move_failed = failed,
      move_acceptance_rate = accepted / pmax(attempted, 1L),
      log_alpha_quantiles = lapply(log_alpha_trace, qfun),
      delta_logpost_quantiles = lapply(delta_logpost_trace, qfun),
      birth_attempted = attempted["birth"],
      birth_accepted = accepted["birth"],
      birth_acceptance_rate = accepted["birth"] / max(1L, attempted["birth"]),
      death_attempted = attempted["death"],
      death_accepted = accepted["death"],
      death_acceptance_rate = accepted["death"] / max(1L, attempted["death"]),
      swap_attempted = attempted["swap"],
      swap_accepted = accepted["swap"],
      swap_acceptance_rate = accepted["swap"] / max(1L, attempted["swap"]),
      candidate_birth_attempted = cand_attempt_birth,
      candidate_birth_accepted = cand_accept_birth,
      candidate_death_attempted = cand_attempt_death,
      candidate_death_accepted = cand_accept_death,
      candidate_swap_in_attempted = cand_attempt_swap_in,
      candidate_swap_in_accepted = cand_accept_swap_in,
      candidate_swap_out_attempted = cand_attempt_swap_out,
      candidate_swap_out_accepted = cand_accept_swap_out,
      lambda_laplace_draws = runtime$lambda_laplace_draws,
      lambda_prior_draws = runtime$lambda_prior_draws,
      lambda_laplace_failures = runtime$lambda_laplace_failures,
      beta_local_draws = runtime$beta_local_draws,
      beta_prior_draws = runtime$beta_prior_draws,
      predictor_scaled_local_draws = runtime$predictor_scaled_local_draws,
      predictor_scaled_prior_draws = runtime$predictor_scaled_prior_draws,
      predictor_scale_mean = runtime$predictor_scale_sum /
        max(1L, runtime$beta_local_draws + runtime$beta_prior_draws),
      predictor_scale_min = runtime$predictor_scale_min,
      association_ess = association_diag$ess,
      association_pip_mcse = association_diag$mcse,
      nonlinear_ess = nonlinear_diag$ess,
      nonlinear_pip_mcse = nonlinear_diag$mcse,
      prior_complexity = prior_summary,
      candidate_block_validation = birth_geometry_check,
      initial_birth_probabilities = attr(birth_geometry_check, "initial_q_b"),
      kernel_operational = isTRUE(mixing_assessment$kernel_operational),
      mixing_within_reference = isTRUE(
        mixing_assessment$mixing_within_reference
      ),
      allocation_diagnostics_within_reference = isTRUE(
        mixing_assessment$mixing_within_reference
      ),
      mixing_diagnostic_indicators = mixing_assessment$diagnostic_indicators,
      mixing_diagnostics_outside_reference =
        mixing_assessment$diagnostics_outside_reference,
      kernel_diagnostics_outside_reference =
        mixing_assessment$kernel_diagnostics_outside_reference,
      mixing_references = mixing_assessment[c(
        "retained_draws_reference", "model_switches_reference",
        "effective_size_reference",
        "pip_mcse_reference", "split_drift_reference",
        "second_eigenvalue_reference"
      )],
      max_association_pip_mcse = mixing_assessment$max_association_pip_mcse,
      max_nonlinear_pip_mcse = mixing_assessment$max_nonlinear_pip_mcse,
      min_effective_size_observed =
        mixing_assessment$min_effective_size_observed,
      proposal_geometry = "predictor_scaled_conditional_identity_map"
    ),
    model_diag
  )
  
  list(
    model_name = label,
    family = family_name,
    trace_logpost = trace_logpost,
    trace_model_size = trace_model_size,
    trace_inclusion = trace_inclusion,
    component_inclusion_probabilities = component_pip,
    posterior_inclusion_probabilities = pip,
    posterior_nonlinear_probabilities = nonlinear_pip,
    posterior_nonlinear_given_association = nonlinear_given_association,
    retained_states = retained_states,
    representative_states = representative_states,
    terminal_state = current,
    admissibility = admissibility,
    diagnostics = diagnostics,
    formula_base = formula_base,
    candidate_covariates = candidate_covariates,
    original_candidate_covariates = original_candidate_covariates,
    candidate_map = candidate_map(candidate_covariates, data),
    hyperparameters = list(
      a_V = a_V, a_A = a_A,
      c_linear = c_linear, c0 = c0, cplus = cplus,
      a_lambda = a_lambda, b_lambda = b_lambda,
      epsilon_birth = epsilon_birth,
      omega_birth = omega_birth,
      proposal_ridge = proposal_ridge,
      within_sweeps = within_sweeps,
      model_thin = model_thin,
      theta_thin = thin,
      swap_weight = swap_weight,
      lambda_prior_mix = lambda_prior_mix,
      beta_prior_mix = beta_prior_mix,
      proposal_link_rms = proposal_link_rms,
      lambda_opt_maxit = lambda_opt_maxit,
      max_frame_cache = max_frame_cache,
      max_representative_states = max_representative_states
    )
  )
}

# Build family-specific structural formula list
family_formula_base <- function(
    fam_obj, base_formula, auxiliary_formulas = list()) {
  params <- fam_obj$names
  form_base <- list()
  form_base[[params[1L]]] <- base_formula
  if (length(params) > 1L) for (p in params[-1L]) {
    form_base[[p]] <-
      if (!is.null(auxiliary_formulas[[p]])) auxiliary_formulas[[p]] else ~ 1
  }
  form_base
}

# Initialize the selected model from the model-search output
selected_model_start <- function(
    search_fit, selected_inclusion, formula_base, data, fam_obj, offset = NULL,
    c_linear = 1.5, c0 = 316.227766, cplus = 316.227766,
    a_lambda = 1, b_lambda = 1) {
  
  states <- search_fit$retained_states
  reps <- search_fit$representative_states
  if (!is.null(reps) && length(reps)) states <- c(unname(reps), states)
  
  if (is.null(states) || !length(states)) {
    stop(
      "Post-selection sampling needs stored model-search states ",
      "(store_theta = TRUE).",
      call. = FALSE
    )
  }
  
  selected_inclusion <- as.matrix(selected_inclusion)
  if (!is.logical(selected_inclusion)) {
    storage.mode(selected_inclusion) <- "logical"
  }
  
  valid_state <- vapply(
    states,
    function(st) {
      !is.null(st) &&
        !is.null(st$inclusion) &&
        !is.null(st$theta) &&
        is.matrix(st$inclusion) &&
        identical(dim(st$inclusion), dim(selected_inclusion)) &&
        identical(dimnames(st$inclusion), dimnames(selected_inclusion))
    },
    logical(1)
  )
  states <- states[valid_state]
  
  if (!length(states)) {
    stop(
      "No stored model-search state matches the layout of the selected ",
      "allocation.",
      call. = FALSE
    )
  }
  
  runtime <- new_runtime_state(cache_frames = FALSE, cache_blocks = TRUE)
  frame_selected <- get_model_frame(
    selected_inclusion, formula_base, data, fam_obj, offset, runtime,
    c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda,
    c0 = c0, cplus = cplus
  )
  fam_selected <- if (!is.null(frame_selected$family)) {
    frame_selected$family
  } else {
    fam_obj
  }
  
  hp <- search_fit$hyperparameters
  if (is.null(hp$a_V) || is.null(hp$a_A)) {
    stop(
      "Model-search result is missing prior hyperparameters a_V/a_A.",
      call. = FALSE
    )
  }
  
  validate_start <- function(start, source, source_logpost = NA_real_) {
    start <- flatten_parameters(start)
    if (!is.numeric(start) || is.null(names(start)) ||
        any(!nzchar(names(start))) || any(!is.finite(start))) {
      return(NULL)
    }
    
    check_state <- list(
      inclusion = selected_inclusion,
      theta = start,
      logPost = -Inf
    )
    comp <- posterior_components(
      check_state, frame_selected, fam_selected,
      offset = frame_selected$model_offset,
      a_V = hp$a_V, a_A = hp$a_A,
      details = TRUE
    )
    if (!is.finite(comp$total)) return(NULL)
    
    message(sprintf(
      paste0(
        "  Post-selection sampler initialized from %s ",
        "(source logPost=%s; selected-model logPost=%.3f)."
      ),
      source,
      if (is.finite(source_logpost)) sprintf("%.3f", source_logpost) else "NA",
      comp$total
    ))
    
    list(
      start = start,
      source = source,
      source_logPost = source_logpost,
      selected_logPost = comp$total
    )
  }
  
  state_logpost <- vapply(states, function(st) {
    z <- suppressWarnings(as.numeric(st$logPost)[1L])
    if (length(z) && is.finite(z)) z else -Inf
  }, numeric(1))
  
  is_exact <- vapply(
    states,
    function(st) all(st$inclusion == selected_inclusion),
    logical(1)
  )
  covers_selected <- vapply(
    states,
    function(st) all(!selected_inclusion | st$inclusion),
    logical(1)
  )
  
  preferred <- if (any(is_exact)) which(is_exact) else which(covers_selected)
  
  if (length(preferred)) {
    j <- preferred[which.max(state_logpost[preferred])]
    chosen <- states[[j]]
    
    template <- generic_initial_parameters(frame_selected, fam_selected)
    template <- seed_candidate_reference_state(
      template, frame_selected,
      a_lambda = a_lambda, b_lambda = b_lambda,
      beta_value = 0
    )
    projected <- merge_parameters(template, chosen$theta)
    
    ans <- validate_start(
      projected,
      if (any(is_exact)) "exact stored allocation" else "stored superset",
      state_logpost[j]
    )
    if (!is.null(ans)) return(ans)
  }
  
  # If the median-probability allocation was never visited, the closest stored
  # state gives the shared parameters; other blocks start at zero.
  missing_selected <- vapply(
    states,
    function(st) sum(selected_inclusion & !st$inclusion),
    integer(1)
  )
  hamming <- vapply(
    states,
    function(st) sum(xor(selected_inclusion, st$inclusion)),
    integer(1)
  )
  ord <- order(missing_selected, hamming, -state_logpost)
  
  if (length(ord)) {
    j <- ord[1L]
    chosen <- states[[j]]
    
    template <- generic_initial_parameters(frame_selected, fam_selected)
    template <- seed_candidate_reference_state(
      template, frame_selected,
      a_lambda = a_lambda, b_lambda = b_lambda,
      beta_value = 0
    )
    projected <- merge_parameters(template, chosen$theta)
    
    ans <- validate_start(
      projected,
      "closest stored allocation (missing blocks start at zero)",
      state_logpost[j]
    )
    if (!is.null(ans)) return(ans)
  }
  
  init_runtime <- new_runtime_state(cache_frames = FALSE, cache_blocks = TRUE)
  init_state <- initialize_model_state(
    inclusion = selected_inclusion,
    frame = frame_selected,
    fam_obj = fam_selected,
    runtime = init_runtime,
    a_V = hp$a_V,
    a_A = hp$a_A
  )
  
  ans <- validate_start(
    init_state$theta,
    "selected-model family/optimizer initialization",
    init_state$logPost
  )
  if (is.null(ans)) {
    stop(
      "No finite initialization could be constructed for the selected ",
      "post-selection model.",
      call. = FALSE
    )
  }
  ans
}

# Posterior evidence for linear candidate components
candidate_effect_evidence <- function(
    search_fit, data, family_name, offset = NULL,
    c_linear = 2.0,
    a_lambda = 1, b_lambda = 1) {
  states <- search_fit$retained_states
  if (is.null(states) || !length(states)) return(data.frame())
  states <- Filter(
    function(st) !is.null(st) && !is.null(st$inclusion) && !is.null(st$theta),
    states
  )
  if (!length(states)) return(data.frame())
  fam_obj <- get_model_family(family_name)
  form_base <- search_fit$formula_base
  inc_template <- search_fit$terminal_state$inclusion
  keys <- rownames(inc_template)
  params <- colnames(inc_template)
  mp <- candidate_map(keys, data)
  lin_rows <- which(mp$component != "nonlinear")
  runtime <- new_runtime_state(cache_frames = FALSE, cache_blocks = TRUE)
  assoc_pip <- search_fit$posterior_inclusion_probabilities
  nl_pip <- search_fit$posterior_nonlinear_probabilities
  out <- list()
  pos <- 0L
  
  theta_value <- function(theta, nm) {
    th <- flatten_parameters(theta)
    if (nm %in% names(th)) return(as.numeric(th[[nm]]))
    hit <- which(endsWith(names(th), nm))
    if (length(hit) == 1L) return(as.numeric(th[[hit]]))
    NA_real_
  }
  finite_quantile <- function(x, p) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    as.numeric(stats::quantile(x, p, names = FALSE, type = 8))
  }
  
  for (r in lin_rows) for (k in seq_along(params)) {
    if (!isTRUE(search_fit$admissibility[r, k])) next
    v <- mp$variable[r]
    par <- params[k]
    block <- tryCatch(
      get_candidate_block(
        r, k, inc_template, form_base, data, fam_obj, offset, runtime,
        c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda
      ), error = function(e) NULL
    )
    if (is.null(block) || !length(block$beta_names)) next
    pip_assoc <- NA_real_
    if (v %in% rownames(assoc_pip)) pip_assoc <- assoc_pip[v, par]
    pip_nonlinear <- NA_real_
    if (!is.null(nl_pip) && v %in% rownames(nl_pip)) {
      pip_nonlinear <- nl_pip[v, par]
    }
    active <- vapply(
      states, function(st) isTRUE(st$inclusion[r, k]), logical(1)
    )
    n_active <- sum(active)
    for (bn in block$beta_names) {
      vals <- rep(0, length(states))
      mapped <- rep(TRUE, length(states))
      if (n_active) {
        av <- vapply(
          states[active], function(st) theta_value(st$theta, bn), numeric(1)
        )
        mapped[active] <- is.finite(av)
        vals[active] <- av
        vals[active & !mapped] <- NA_real_
      }
      cond <- vals[active & mapped]
      ma <- vals[mapped]
      pos <- pos + 1L
      out[[pos]] <- data.frame(
        Parameter = par,
        Covariate = v,
        Coefficient = bn,
        PIP_Association = pip_assoc,
        PIP_Nonlinear = pip_nonlinear,
        Active_Draws = n_active,
        Mapped_Active_Draws = length(cond),
        Mapping_Rate_Active = if (n_active) {
          length(cond) / n_active
        } else {
          NA_real_
        },
        Conditional_Mean = if (length(cond)) mean(cond) else NA_real_,
        Conditional_Median = if (length(cond)) {
          stats::median(cond)
        } else {
          NA_real_
        },
        Conditional_Q025 = finite_quantile(cond, .025),
        Conditional_Q975 = finite_quantile(cond, .975),
        Pr_Positive_Given_Included = if (length(cond)) {
          mean(cond > 0)
        } else {
          NA_real_
        },
        Pr_Negative_Given_Included = if (length(cond)) {
          mean(cond < 0)
        } else {
          NA_real_
        },
        Model_Averaged_Mean = if (length(ma)) mean(ma) else NA_real_,
        Model_Averaged_Median = if (length(ma)) stats::median(ma) else NA_real_,
        Model_Averaged_Q025 = finite_quantile(ma, .025),
        Model_Averaged_Q975 = finite_quantile(ma, .975),
        Retained_Theta_Draws = length(states),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(out)) return(data.frame())
  do.call(rbind, out)
}

# Key of what a checkpoint was computed from: family implementation and fit
# settings, deparsed without source references.
checkpoint_key_legacy <- function(fam_obj, settings) {
  ctrl <- c("keepNA", "keepInteger", "niceNames", "showAttributes")
  paste(c(
    fam_obj$family, fam_obj$names,
    deparse(fam_obj$d, control = ctrl),
    deparse(fam_obj$loglik, control = ctrl),
    deparse(settings, control = ctrl)
  ), collapse = "\n")
}

# The key ignores code text: only a change of the density used (revision
# below) or of the settings invalidates it. Legacy keys are still accepted.
family_density_revision <- function(fam_obj) {
  # SICHEL keeps the tag of the density it had before it moved to plain R:
  # the same probabilities, so an earlier search remains valid.
  fixed <- c(
    ZANBI = "log-space", ZOINB = "1", RSB = "1", SICHEL = "batch-safe"
  )
  nm <- as.character(fam_obj$family)[1L]
  if (nm %in% names(fixed)) return(fixed[[nm]])
  if (isTRUE(attr(fam_obj, "batch_sensitive"))) "batch-safe" else "1"
}

checkpoint_key <- function(fam_obj, settings) {
  ctrl <- c("keepNA", "keepInteger", "niceNames", "showAttributes")
  paste(c(
    "key2", fam_obj$family, fam_obj$names, family_density_revision(fam_obj),
    deparse(settings, control = ctrl)
  ), collapse = "\n")
}

checkpoint_keys <- function(fam_obj, settings) {
  ctrl <- c("keepNA", "keepInteger", "niceNames", "showAttributes")
  structure(
    c(
      checkpoint_key(fam_obj, settings),
      checkpoint_key_legacy(fam_obj, settings)
    ),
    family = as.character(fam_obj$family)[1L],
    settings = paste(deparse(settings, control = ctrl), collapse = "\n"),
    batch_sensitive = isTRUE(attr(fam_obj, "batch_sensitive"))
  )
}

# A legacy key holds if family, settings and density are unchanged; batch-
# sensitive families count only if written with the observation-safe density.
earlier_key_holds <- function(saved, key) {
  startsWith(saved, paste0(attr(key, "family"), "\n")) &&
    endsWith(saved, attr(key, "settings")) &&
    (!isTRUE(attr(key, "batch_sensitive")) ||
       grepl("eval_group", saved, fixed = TRUE))
}

read_checkpoint <- function(path, key, what, legacy_ok = TRUE) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  obj <- readRDS(path)
  saved <- attr(obj, "checkpoint_key")
  if (is.null(saved) && !legacy_ok) {
    message(
      "  [Checkpoint] ", basename(path), " has no key and was fitted with a ",
      "density that depends on batching; ", what, " is refitted."
    )
    return(NULL)
  }
  if (is.null(saved)) {
    message(
      "  [Checkpoint] ", basename(path), " was written without a key by an ",
      "older version of these functions and is reused as is; delete it to ",
      "refit ", what, " with the current code."
    )
    return(obj)
  }
  if (!saved %in% key && earlier_key_holds(saved, key)) {
    message(
      "  [Checkpoint] ", basename(path), " was written by an earlier ",
      "version of these functions with the same settings and is reused."
    )
    return(obj)
  }
  if (!saved %in% key) {
    message(
      "  [Checkpoint] ", basename(path), " belongs to a different family ",
      "implementation or settings; ", what, " is refitted."
    )
    return(NULL)
  }
  obj
}

write_checkpoint <- function(obj, path, key) {
  if (is.null(path)) return(invisible(NULL))
  attr(obj, "checkpoint_key") <- key[1L]
  saveRDS(obj, path)
  invisible(path)
}

# Fit models across distribution families
fit_all_spatial_families <- function(
    families, base_formula, candidate_pool, data, penalty_mat,
    auxiliary_formulas = list(), offsets = NULL,
    candidate_admissibility = NULL,
    n_iter = 10000, n_burn = 4000, thin = 5, model_thin = 1L,
    c_linear = 1.5,
    c0 = 316.227766, cplus = 316.227766,
    a_lambda = 1, b_lambda = 1,
    pip_threshold = 0.5,
    support_threshold = 0.30,
    continue_on_family_error = TRUE,
    run_posterior_sampling = TRUE,
    posterior_n_iter = 6000, posterior_n_burn = 2000, posterior_thin = 10,
    posterior_chains = 4L, posterior_cores = NULL,
    chain_seed = NULL,
    reference_min_chains = 2L,
    reference_min_draws = 400L,
    ess_reference = 100,
    split_drift_reference = 0.50,
    rhat_reference = 1.05,
    checkpoint_dir = NULL, resume_every = 500L,
    posterior_segment = 2000L,
    pooled_search_dirs = NULL, linearize = character(0),
    ...) {
  
  resume_every <- as.integer(resume_every)[1L]
  if (!is.finite(resume_every) || resume_every < 1L) {
    stop("`resume_every` must be a positive integer.", call. = FALSE)
  }
  if (!is.finite(pip_threshold) || pip_threshold <= 0 || pip_threshold >= 1) {
    stop("`pip_threshold` must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.finite(support_threshold) || support_threshold <= 0 ||
      support_threshold >= pip_threshold) {
    stop("`support_threshold` must lie in (0, pip_threshold).", call. = FALSE)
  }
  if (!missing(penalty_mat) && !is.null(penalty_mat)) {
    penalty_mat <- as.matrix(penalty_mat)
    if (nrow(penalty_mat) != ncol(penalty_mat) ||
        any(!is.finite(penalty_mat)) ||
        !isTRUE(all.equal(penalty_mat, t(penalty_mat), tolerance = 1e-8))) {
      stop(
        "`penalty_mat` must be a finite symmetric square matrix.", call. = FALSE
      )
    }
  }
  
  if (!is.null(checkpoint_dir) && !dir.exists(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  data <- prepare_candidate_data(data, candidate_pool)
  fitted_models <- list()
  model_formulas <- list()
  selected_vars_list <- list()
  supported_vars_list <- list()
  selection_types <- list()
  effect_evidence <- list()
  search_diagnostics <- list()
  search_results <- list()
  status_rows <- list()
  
  empty_by_param <- function(parameters) stats::setNames(
    replicate(length(parameters), character(0), simplify = FALSE), parameters
  )
  
  y_name <- all.vars(base_formula)[1L]
  search_settings <- list(
    base_formula = base_formula, auxiliary_formulas = auxiliary_formulas,
    candidate_pool = candidate_pool,
    candidate_admissibility = candidate_admissibility,
    offsets = offsets, n = nrow(data), y_total = sum(data[[y_name]]),
    n_iter = n_iter, n_burn = n_burn, thin = thin, model_thin = model_thin,
    c_linear = c_linear, c0 = c0, cplus = cplus,
    a_lambda = a_lambda, b_lambda = b_lambda, search_args = list(...),
    structural_sd = getOption("rbst.structural_sd"),
    mrf_theta = getOption("rbst.mrf_theta")
  )
  posterior_settings <- list(
    pip_threshold = pip_threshold, support_threshold = support_threshold,
    run_posterior_sampling = run_posterior_sampling,
    posterior_n_iter = posterior_n_iter, posterior_n_burn = posterior_n_burn,
    posterior_thin = posterior_thin, posterior_chains = posterior_chains,
    chain_seed = chain_seed
  )
  if (length(pooled_search_dirs)) {
    posterior_settings$pooled_search <- basename(pooled_search_dirs)
  }
  if (length(linearize)) posterior_settings$linearize <- sort(linearize)
  
  for (family_spec in families) {
    family_label <- "<custom>"
    if (is.character(family_spec) && length(family_spec) == 1L) {
      family_label <- family_spec
    } else if (is.list(family_spec) && !is.null(family_spec$family)) {
      family_label <- as.character(family_spec$family)[1L]
    }
    
    chk_file <- search_file <- resume_file <- NULL
    if (!is.null(checkpoint_dir)) {
      chk_file <- file.path(
        checkpoint_dir, paste0("checkpoint_", family_label, ".rds")
      )
      search_file <- file.path(
        checkpoint_dir, paste0("search_", family_label, ".rds")
      )
      resume_file <- file.path(
        checkpoint_dir, paste0("resume_", family_label, ".rds")
      )
    }
    
    # Keys are only available when the family itself can be built; if it
    # cannot, the error is reported by the fit below.
    keys <- tryCatch({
      fam_key_obj <- get_model_family(family_spec)
      list(
        search = checkpoint_keys(fam_key_obj, search_settings),
        full = checkpoint_keys(
          fam_key_obj, c(search_settings, posterior_settings)
        )
      )
    }, error = function(e) NULL)
    
    family_result <- if (is.null(keys)) NULL else {
      read_checkpoint(
        chk_file, keys$full, "the family",
        legacy_ok = !isTRUE(attr(fam_key_obj, "batch_sensitive"))
      )
    }
    
    if (!is.null(family_result)) {
      message(
        sprintf("\n=== Distribution family: %s (checkpoint) ===", family_label)
      )
      search_result_saved <- family_result$model_search
      
    } else {
      message(sprintf("\n=== Distribution family: %s ===", family_label))
      stage <- "family validation"
      search_result_saved <- NULL
      
      family_result <- tryCatch({
        fam_obj <- get_model_family(family_spec)
        parameters <- fam_obj$names
        form_base <- family_formula_base(
          fam_obj, base_formula, auxiliary_formulas
        )
        family_offset <- resolve_family_offset(offsets, family_label)
        
        stage <- "model search"
        search_result <- read_checkpoint(
          search_file, keys$search, "the model search",
          legacy_ok = !isTRUE(attr(fam_key_obj, "batch_sensitive"))
        )
        if (!is.null(search_result)) {
          message(
            "  [Checkpoint] Reusing the stored model search for ", family_label,
            "."
          )
        } else {
          search_result <- fit_spatial_model_dynamic(
            label = family_label,
            formula_base = form_base,
            candidate_covariates = candidate_pool,
            family_name = family_spec,
            data = data,
            offset = family_offset,
            candidate_admissibility = candidate_admissibility,
            n_iter = n_iter, n_burn = n_burn, thin = thin,
            model_thin = model_thin,
            c_linear = c_linear, c0 = c0, cplus = cplus,
            a_lambda = a_lambda, b_lambda = b_lambda,
            resume_file = resume_file, resume_key = keys$search[1L],
            resume_every = resume_every,
            ...
          )
          # Saved before any post-selection step, so a later failure does
          # not throw away the search.
          write_checkpoint(search_result, search_file, keys$search)
          if (!is.null(resume_file)) unlink(resume_file)
        }
        if (length(pooled_search_dirs)) {
          extra <- lapply(pooled_search_dirs, function(d) {
            read_checkpoint(
              file.path(d, paste0("search_", family_label, ".rds")),
              keys$search, "a replicate search"
            )
          })
          n_found <- sum(!vapply(extra, is.null, logical(1)))
          message(
            "  [Pooling] ", n_found + 1L, " independent searches for ",
            family_label, "."
          )
          search_result <- pool_search_fits(c(list(search_result), extra))
        }
        search_result_saved <- search_result
        
        diagnostics <- search_result$diagnostics
        kernel_operational <- isTRUE(diagnostics$kernel_operational)
        mixing_within_reference <- isTRUE(diagnostics$mixing_within_reference)
        allocation_audit <- selection_audit(
          search_result, threshold = pip_threshold
        )
        undecided <- allocation_audit$Decision == "Inconclusive"
        assoc_decision <- list(
          stable = !any(undecided & allocation_audit$Level == "association"),
          unstable = character(0)
        )
        nl_decision <- list(
          stable = !any(undecided & allocation_audit$Level == "nonlinear"),
          unstable = character(0)
        )
        
        threshold_stable <-
          isTRUE(assoc_decision$stable) && isTRUE(nl_decision$stable)
        
        diagnostics$selection_threshold_stable <- threshold_stable
        diagnostics$allocation_diagnostics_within_reference <-
          isTRUE(mixing_within_reference) && isTRUE(threshold_stable)
        diagnostics$unstable_threshold_decisions <- paste(
          allocation_audit$Key[undecided],
          allocation_audit$Parameter[undecided],
          sep = "::"
        )
        diagnostics$selection_audit <- allocation_audit
        outside <- character(0)
        if (!kernel_operational) {
          outside <- diagnostics$kernel_diagnostics_outside_reference
        }
        if (!mixing_within_reference) {
          outside <- c(
            outside, diagnostics$mixing_diagnostics_outside_reference
          )
        }
        if (!threshold_stable) {
          outside <- c(outside, "selection_threshold_instability")
        }
        diagnostics$diagnostics_outside_reference <- unique(outside)
        search_result$diagnostics <- diagnostics
        search_result_saved <- search_result
        
        if (length(diagnostics$diagnostics_outside_reference)) {
          warning(
            "Model-space diagnostics outside their reference values: ",
            paste(diagnostics$diagnostics_outside_reference, collapse = ", "),
            ". Results are kept for review.",
            call. = FALSE
          )
        }
        
        if (!isTRUE(run_posterior_sampling)) {
          list(
            completed = TRUE,
            status = "diagnostic_only",
            status_message = "Post-selection sampling disabled.",
            family = family_label,
            model = NULL,
            formulas = NULL,
            selections = empty_by_param(parameters),
            supported_selections = empty_by_param(parameters),
            selection_types = empty_by_param(parameters),
            effect_evidence = data.frame(),
            diagnostics = diagnostics,
            model_search = search_result
          )
        } else {
          pip <- search_result$posterior_inclusion_probabilities
          pip_nl <- search_result$posterior_nonlinear_probabilities
          component_keys <- rownames(
            search_result$component_inclusion_probabilities
          )
          selected_inclusion <- matrix(
            FALSE,
            nrow = length(component_keys), ncol = length(parameters),
            dimnames = list(component_keys, parameters)
          )
          selected_by_parameter <- list()
          supported_by_parameter <- list()
          effect_type_by_parameter <- list()
          
          for (param in parameters) {
            selected <- character(0)
            types <- character(0)
            for (v in rownames(pip)) {
              association_prob <- pip[v, param]
              nonlinear_prob <- if (!is.null(pip_nl) &&
                                    v %in% rownames(pip_nl)) {
                pip_nl[v, param]
              } else {
                NA_real_
              }
              linear_key <- paste0("L::", v)
              nonlinear_key <- paste0("NL::", v)
              if (is.finite(association_prob) &&
                  association_prob > pip_threshold &&
                  linear_key %in% component_keys) {
                selected_inclusion[linear_key, param] <- TRUE
                selected <- c(selected, v)
                nonlinear_selected <- is.finite(nonlinear_prob) &&
                  nonlinear_prob > pip_threshold &&
                  nonlinear_key %in% component_keys
                if (nonlinear_selected) {
                  selected_inclusion[nonlinear_key, param] <- TRUE
                }
                x <- data[[v]]
                is_continuous <- is.numeric(x) &&
                  sum(is.finite(x)) > 0L &&
                  length(unique(x[is.finite(x)])) > 7L
                types[v] <- if (nonlinear_selected) {
                  "nonlinear"
                } else if (is_continuous) {
                  "linear"
                } else {
                  "categorical"
                }
              }
            }
            selected_by_parameter[[param]] <- selected
            keep <- is.finite(pip[, param]) & pip[, param] >= support_threshold
            supported <- rownames(pip)[keep]
            supported_by_parameter[[param]] <- supported
            effect_type_by_parameter[[param]] <- types
            selection_label <- if (length(selected)) {
              paste0(
                selected, " [", unname(types[selected]), "]", collapse = ", "
              )
            } else {
              "<none>"
            }
            message(
              sprintf(
                "  -> Parameter [%s] selected covariates: %s", param,
                selection_label
              )
            )
          }
          
          # Sensitivity variant: the listed smooth effects enter linearly
          for (v in linearize) {
            nl_key <- paste0("NL::", v)
            if (!nl_key %in% rownames(selected_inclusion)) next
            was <- selected_inclusion[nl_key, ]
            selected_inclusion[nl_key, ] <- FALSE
            for (param in names(was)[was]) {
              effect_type_by_parameter[[param]][v] <- "linear"
            }
            if (any(was)) {
              message(
                "  -> Smooth effect of ", v, " set linear in: ",
                paste(names(was)[was], collapse = ", ")
              )
            }
          }
          
          if (!hierarchy_valid(selected_inclusion)) {
            stop(
              "Thresholded allocation violates the NL <= L hierarchy.",
              call. = FALSE
            )
          }
          selected_formulas <- build_distributional_formula(
            selected_inclusion, form_base, data, fam_obj,
            c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda
          )
          formula_text <- lapply(
            selected_formulas,
            function(x) paste(deparse(x, width.cutoff = 500), collapse = " ")
          )
          
          stage <- "posterior evidence summary"
          evidence <- candidate_effect_evidence(
            search_result, data = data, family_name = family_spec,
            offset = family_offset,
            c_linear = c_linear,
            a_lambda = a_lambda, b_lambda = b_lambda
          )
          
          stage <- "post-selection initialization"
          initial_state_info <- selected_model_start(
            search_fit = search_result,
            selected_inclusion = selected_inclusion,
            formula_base = form_base,
            data = data,
            fam_obj = fam_obj,
            offset = family_offset,
            c_linear = c_linear, c0 = c0, cplus = cplus,
            a_lambda = a_lambda, b_lambda = b_lambda
          )
          
          stage <- "post-selection sampler"
          fitted_model <- fit_spatial_model(
            label = family_label,
            formula_list = selected_formulas,
            family_name = family_spec,
            data = data,
            offset = family_offset,
            formula_base = form_base,
            candidate_inclusion = selected_inclusion,
            candidate_covariates = component_keys,
            c_linear = c_linear, c0 = c0, cplus = cplus,
            a_lambda = a_lambda, b_lambda = b_lambda,
            start = initial_state_info$start,
            optimizer = FALSE,
            sampler = bamlss::sam_GMCMC,
            n_iter = posterior_n_iter,
            n_burn = posterior_n_burn,
            thin = posterior_thin,
            chains = posterior_chains,
            cores = posterior_cores,
            chain_seed = chain_seed,
            chain_checkpoint = if (!is.null(checkpoint_dir)) {
              file.path(checkpoint_dir, paste0("post_", family_label))
            },
            segment_iter = posterior_segment,
            reference_min_chains = reference_min_chains,
            reference_min_draws = reference_min_draws,
            ess_reference = ess_reference,
            split_drift_reference = split_drift_reference,
            rhat_reference = rhat_reference
          )
          sampler_diagnostics <- attr(fitted_model, "sampler_diagnostics")
          sampler_diagnostics_within_reference <-
            !is.null(sampler_diagnostics) &&
            isTRUE(sampler_diagnostics$diagnostics_within_reference)
          
          attr(fitted_model, "posterior_start") <- initial_state_info
          attr(fitted_model, "selection_types") <- effect_type_by_parameter
          sampler_outside <- character(0)
          if (!is.null(sampler_diagnostics)) {
            sampler_outside <- sampler_diagnostics$diagnostics_outside_reference
          }
          attr(fitted_model, "diagnostic_summary") <- list(
            model_space_within_reference = isTRUE(
              diagnostics$allocation_diagnostics_within_reference
            ),
            posterior_sampler_within_reference =
              sampler_diagnostics_within_reference,
            model_space_diagnostics_outside_reference =
              diagnostics$diagnostics_outside_reference,
            sampler_diagnostics_outside_reference = sampler_outside
          )
          
          list(
            completed = TRUE, status = "completed", status_message = "",
            family = family_label,
            model = fitted_model,
            formulas = formula_text,
            selections = selected_by_parameter,
            supported_selections = supported_by_parameter,
            selection_types = effect_type_by_parameter,
            effect_evidence = evidence,
            diagnostics = diagnostics,
            model_search = search_result
          )
        }
      }, error = function(e) {
        list(
          completed = FALSE, family = family_label, stage = stage,
          error = conditionMessage(e)
        )
      })
      
      if (isTRUE(family_result$completed) && !is.null(keys)) {
        write_checkpoint(family_result, chk_file, keys$full)
      }
    }
    
    if (!isTRUE(family_result$completed)) {
      status_rows[[length(status_rows) + 1L]] <- data.frame(
        Family = family_label, Status = "failed", Stage = family_result$stage,
        Message = family_result$error, stringsAsFactors = FALSE
      )
      if (!is.null(search_result_saved)) {
        search_results[[family_label]] <- search_result_saved
        search_diagnostics[[family_label]] <- search_result_saved$diagnostics
      }
      if (!isTRUE(continue_on_family_error)) {
        stop(
          "Family '", family_label, "' failed during ", family_result$stage,
          ": ", family_result$error,
          call. = FALSE
        )
      }
      next
    }
    
    if (!is.null(family_result$model)) {
      fitted_models[[family_label]] <- family_result$model
    }
    if (!is.null(family_result$formulas)) {
      model_formulas[[family_label]] <- family_result$formulas
    }
    selected_vars_list[[family_label]] <- family_result$selections
    supported_vars_list[[family_label]] <- family_result$supported_selections
    selection_types[[family_label]] <- family_result$selection_types
    if (!is.null(family_result$effect_evidence) &&
        nrow(family_result$effect_evidence)) {
      effect_evidence[[family_label]] <- family_result$effect_evidence
    }
    search_diagnostics[[family_label]] <- family_result$diagnostics
    search_results[[family_label]] <- family_result$model_search
    stage_label <- "model search"
    if (identical(family_result$status, "completed")) stage_label <- "completed"
    status_rows[[length(status_rows) + 1L]] <- data.frame(
      Family = family_label,
      Status = family_result$status,
      Stage = stage_label,
      Message = family_result$status_message,
      stringsAsFactors = FALSE
    )
  }
  
  family_status <- if (length(status_rows)) do.call(rbind, status_rows) else {
    data.frame(
      Family = character(0), Status = character(0), Stage = character(0),
      Message = character(0)
    )
  }
  
  list(
    models = fitted_models,
    formulas = model_formulas,
    selections = selected_vars_list,
    supported_selections = supported_vars_list,
    selection_types = selection_types,
    effect_evidence = effect_evidence,
    diagnostics = search_diagnostics,
    model_search = search_results,
    family_status = family_status,
    pip_threshold = pip_threshold,
    support_threshold = support_threshold,
    run_posterior_sampling = isTRUE(run_posterior_sampling),
    diagnostic_references = list(
      reference_min_chains = as.integer(reference_min_chains),
      posterior_chains = as.integer(posterior_chains),
      concurrent_chains = as.integer(min(c(posterior_cores, posterior_chains))),
      chain_seed = chain_seed,
      reference_min_draws = as.integer(reference_min_draws),
      ess_reference = ess_reference,
      split_drift_reference = split_drift_reference,
      rhat_reference = rhat_reference
    )
  )
}

# Subset an external offset list
subset_offset <- function(offset, idx) {
  if (is.null(offset)) return(NULL)
  if (!is.list(offset)) return(offset[idx])
  lapply(offset, function(z) {
    if (length(z) == 1L) z else z[idx]
  })
}

# Convert a retained model state into a lightweight BAMLSS object
build_state_object <- function(state, formula_base, data, fam_obj, offset,
                               c_linear = 2.0, a_lambda = 1, b_lambda = 1) {
  form <- build_distributional_formula(
    state$inclusion, formula_base, data, fam_obj,
    c_linear, a_lambda, b_lambda)
  parsed_form <- parse_distributional_formula(form, fam_obj)
  bf <- bamlss::bamlss.frame(
    parsed_form, data = data, family = fam_obj, offset = offset
  )
  bf$parameters <- state$theta
  class(bf) <- unique(c("bamlss", class(bf)))
  bf
}

# Predict parameters from one retained model state
predict_state_parameters <- function(
    state, formula_base, train_data, newdata, fam_obj,
    offset_train = NULL, offset_new = NULL,
    c_linear = 2.0, a_lambda = 1, b_lambda = 1) {
  
  obj <- build_state_object(
    state, formula_base, train_data, fam_obj, offset_train,
    c_linear, a_lambda, b_lambda)
  eta <- lapply(fam_obj$names, function(p_name) {
    as.numeric(stats::predict(
      obj, newdata = newdata, model = p_name,
      type = "link", what = "parameters"
    ))
  })
  names(eta) <- fam_obj$names
  eta <- add_offset(eta, offset_new, nrow(newdata))
  fam_obj$map2par(eta)
}

# BMA predictive scores from the joint model-search chain
compute_bma_scores <- function(
    search_fit, formula_base, train_data, newdata, y_col,
    family_name = search_fit$family, offset_train = NULL, offset_new = NULL,
    max_states = 60L, tail_tol = 1e-10, zmax_cap = 100000L,
    c_linear = NULL, a_lambda = NULL, b_lambda = NULL) {
  
  states <- search_fit$retained_states
  if (!length(states)) {
    stop("Model-search fit has no retained joint states.", call. = FALSE)
  }
  if (length(states) > max_states) {
    keep <- unique(round(seq(1, length(states), length.out = max_states)))
    states <- states[keep]
  }
  
  hp <- search_fit$hyperparameters
  if (is.null(c_linear)) c_linear <- hp$c_linear
  if (is.null(a_lambda)) a_lambda <- hp$a_lambda
  if (is.null(b_lambda)) b_lambda <- hp$b_lambda
  
  fam_obj <- get_model_family(family_name)
  y <- as.integer(newdata[[y_col]])
  n <- length(y)
  
  pred_states <- lapply(states, function(st) {
    predict_state_parameters(
      st, formula_base, train_data, newdata, fam_obj, offset_train,
      offset_new, c_linear, a_lambda, b_lambda
    )
  })
  n_states <- length(pred_states)
  
  # model-averaged predictive mass at the values in z
  bma_mass <- function(z) {
    dens <- lapply(pred_states, function(pa) {
      as.numeric(fam_obj$d(z, pa, log = FALSE))
    })
    Reduce("+", dens) / n_states
  }
  p_obs <- pmax(bma_mass(y), .Machine$double.xmin)
  
  norm2 <- numeric(n)
  cdf <- numeric(n)
  compute_crps <- !identical(fam_obj$family, "RSB")
  crps_loss <- if (compute_crps) numeric(n) else NULL
  rps_loss <- numeric(n)
  k_top <- max(y)
  
  z_max <- 0L
  zmax_cap <- min(as.integer(zmax_cap), as.integer(max(1000L, max(y) * 3L)))
  tail_tolerance_reached <- FALSE
  repeat {
    pz <- bma_mass(rep(z_max, n))
    pz[!is.finite(pz) | pz < 0] <- 0
    norm2 <- norm2 + pz^2
    cdf <- pmin(1, cdf + pz)
    sq <- (cdf - as.numeric(y <= z_max))^2
    if (compute_crps) crps_loss <- crps_loss + sq
    if (z_max < k_top) rps_loss <- rps_loss + sq
    if (z_max >= max(y) && all(1 - cdf <= tail_tol)) {
      tail_tolerance_reached <- TRUE
      break
    }
    z_max <- z_max + 1L
    if (z_max > zmax_cap) break
  }
  
  data.frame(
    QS = mean(2 * p_obs - norm2),
    LS = mean(log(p_obs)),
    SPS = mean(p_obs / sqrt(pmax(norm2, .Machine$double.xmin))),
    CRPS = if (compute_crps && tail_tolerance_reached) {
      mean(crps_loss)
    } else {
      NA_real_
    },
    RPS_K = mean(rps_loss),
    TailMass = max(1 - cdf),
    BMA_States = n_states,
    Z_Max = z_max,
    Tail_Tolerance_Reached = tail_tolerance_reached
  )
}

# Prospective time blocks
make_prospective_time_blocks <- function(
    data, time_col, n_blocks = 3L, min_train_frac = 0.60) {
  tt <- sort(unique(data[[time_col]]))
  if (length(tt) < 4L) {
    stop(
      "Too few unique time points for prospective validation.", call. = FALSE
    )
  }
  first_test <- max(2L, floor(length(tt) * min_train_frac) + 1L)
  test_times <- tt[first_test:length(tt)]
  groups <- split(
    test_times, cut(seq_along(test_times), breaks = n_blocks, labels = FALSE)
  )
  lapply(groups, function(g) {
    test_idx <- which(data[[time_col]] %in% g)
    train_idx <- which(data[[time_col]] < min(g))
    list(train_idx = train_idx, test_idx = test_idx,
         train_end = max(data[[time_col]][train_idx]),
         test_start = min(g), test_end = max(g))
  })
}

# Prospective out-of-sample BMA validation
validate_prospective_bma <- function(
    families, base_formula, candidate_pool, data,
    auxiliary_formulas = list(), offsets = NULL,
    candidate_admissibility = NULL,
    time_col = "date_obj", n_blocks = 3L, min_train_frac = 0.60,
    n_iter = 1200, n_burn = 400, thin = 4,
    max_bma_states = 60L,
    c_linear = 1.5,
    c0 = 316.227766, cplus = 316.227766,
    a_lambda = 1, b_lambda = 1,
    ...) {
  
  data <- prepare_candidate_data(data, candidate_pool)
  splits <- make_prospective_time_blocks(
    data, time_col, n_blocks, min_train_frac
  )
  out <- list()
  pos <- 0L
  
  for (b in seq_along(splits)) {
    tr <- splits[[b]]$train_idx
    te <- splits[[b]]$test_idx
    if (!length(tr) || !length(te)) next
    dtr <- data[tr, , drop = FALSE]
    dte <- data[te, , drop = FALSE]
    
    for (fam in families) {
      fam_obj <- get_model_family(fam)
      form_base <- family_formula_base(
        fam_obj, base_formula, auxiliary_formulas
      )
      off_full <- resolve_family_offset(offsets, fam)
      off_tr <- subset_offset(off_full, tr)
      off_te <- subset_offset(off_full, te)
      
      fit <- fit_spatial_model_dynamic(
        label = sprintf("OOS_%s_B%d", fam, b),
        formula_base = form_base,
        candidate_covariates = candidate_pool,
        family_name = fam,
        data = dtr,
        offset = off_tr,
        candidate_admissibility = candidate_admissibility,
        n_iter = n_iter, n_burn = n_burn, thin = thin,
        c_linear = c_linear, c0 = c0, cplus = cplus,
        a_lambda = a_lambda, b_lambda = b_lambda,
        store_theta = TRUE,
        ...
      )
      mixing_within_reference <- isTRUE(fit$diagnostics$mixing_within_reference)
      sc <- compute_bma_scores(
        fit, form_base, dtr, dte,
        y_col = all.vars(base_formula)[1L], family_name = fam,
        offset_train = off_tr, offset_new = off_te,
        max_states = max_bma_states,
        c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda
      )
      oos_status <- "completed"
      if (!mixing_within_reference) {
        oos_status <- "completed_outside_reference_values"
      }
      pos <- pos + 1L
      out[[pos]] <- cbind(
        data.frame(Family = fam, Block = b, OOS_Status = oos_status,
                   Train_N = length(tr), Test_N = length(te),
                   Train_End = as.character(splits[[b]]$train_end),
                   Test_Start = as.character(splits[[b]]$test_start),
                   Test_End = as.character(splits[[b]]$test_end)),
        sc
      )
    }
  }
  dplyr::bind_rows(out)
}

# Build a BAMLSS transformer for the selected model
make_selected_model_transform <- function(
    candidate_inclusion, formula_base, data, offset,
    c_linear = 2.0, c0 = 316.227766, cplus = 316.227766,
    a_lambda = 1, b_lambda = 1) {
  
  runtime <- new_runtime_state(cache_frames = FALSE, cache_blocks = TRUE)
  
  function(bf, ...) {
    if (!inherits(bf, "bamlss.frame")) {
      stop(
        "Selected-model transformer requires a bamlss.frame object.",
        call. = FALSE
      )
    }
    
    raw_x <- bf$x
    bf$x <- bamlss::bamlss.engine.setup(bf$x, propose = "iwlsC_gp")
    bf$raw_design <- raw_x
    bf$model_offset <- offset
    
    bf <- install_candidate_priors(
      bf,
      inclusion = candidate_inclusion,
      formula_base = formula_base,
      data = data,
      fam_obj = bf$family,
      offset = offset,
      runtime = runtime,
      c_linear = c_linear, c0 = c0, cplus = cplus,
      a_lambda = a_lambda, b_lambda = b_lambda
    )
    
    slice_proposal <- slice_with_offset
    candidate_meta <- bf$candidate_metadata
    
    if (!is.null(candidate_meta) && length(candidate_meta)) {
      for (block in candidate_meta) {
        if (!identical(block$component, "linear")) {
          next
        }
        
        p_name <- block$parameter
        label <- block$label
        term <- bf$x[[p_name]]$smooth.construct[[label]]
        
        if (is.null(term)) {
          stop(
            "Selected posterior term `", p_name, "::", label, "` is missing.",
            call. = FALSE
          )
        }
        
        term$propose <- slice_proposal
        bf$x[[p_name]]$smooth.construct[[label]] <- term
      }
    }
    
    bf
  }
}

# Extract posterior samples by chain
extract_sample_chains <- function(model) {
  smp <- model$samples
  if (is.null(smp)) {
    stop("Fitted model contains no posterior sample object.", call. = FALSE)
  }
  chains <- if (inherits(smp, "mcmc.list")) {
    lapply(smp, as.matrix)
  } else {
    list(as.matrix(smp))
  }
  report_failed_chains(chains)
  valid <- vapply(
    chains,
    function(z) !is.null(z) &&
      nrow(z) > 0L &&
      ncol(z) > 0L &&
      !is.null(colnames(z)),
    logical(1)
  )
  if (!all(valid)) {
    stop(
      "Posterior samples contain an empty chain or missing parameter names.",
      call. = FALSE
    )
  }
  ref_names <- colnames(chains[[1L]])
  same_names <- vapply(chains, function(z) identical(colnames(z), ref_names),
                       logical(1))
  if (!all(same_names)) {
    stop(
      "Posterior sample chains have inconsistent parameter columns.",
      call. = FALSE
    )
  }
  bad <- nonfinite_sample_table(chains)
  if (nrow(bad)) {
    print(bad, row.names = FALSE)
    stop(
      "Posterior sample chains contain non-finite values in: ",
      paste(unique(bad$column), collapse = ", "), ".",
      call. = FALSE
    )
  }
  chains
}

# mclapply() turns an error inside a chain into a one-cell non-numeric chain
# named var1; guard_sampler() keeps the message in its column name, reported
# here instead of a non-finite value.
guard_sampler <- function(sampler) {
  force(sampler)
  function(...) {
    out <- tryCatch(sampler(...), error = function(e) e)
    if (!inherits(out, "error")) return(out)
    msg <- substr(gsub("[\r\n]+", " ", conditionMessage(out)), 1L, 400L)
    matrix(
      NA_real_, 1L, 1L, dimnames = list(NULL, paste0("sampler_error: ", msg))
    )
  }
}

report_failed_chains <- function(chains) {
  msg <- vapply(seq_along(chains), function(i) {
    z <- chains[[i]]
    cn <- colnames(z)
    if (length(cn) == 1L && startsWith(cn, "sampler_error: ")) {
      return(sub("^sampler_error: ", "", cn))
    }
    single_cell <- length(z) == 1L && !is.finite(suppressWarnings(
      as.numeric(z[1L])
    ))
    if (!is.numeric(z) || single_cell) {
      first <- if (length(z)) as.character(z[1L]) else "empty result"
      return(paste0("returned no sample matrix (", first, ")"))
    }
    ""
  }, character(1))
  if (!any(nzchar(msg))) return(invisible(NULL))
  stop(
    "Sampler failed in chain(s) ",
    paste(which(nzchar(msg)), collapse = ", "), ": ",
    paste(unique(msg[nzchar(msg)]), collapse = " | "),
    call. = FALSE
  )
}

# One row per chain and sample column holding NA/NaN/Inf, with the first
# offending draw, so a failure points at the parameter or sampler statistic.
nonfinite_sample_table <- function(chains) {
  rows <- lapply(seq_along(chains), function(i) {
    nf <- !is.finite(chains[[i]])
    cols <- which(colSums(nf) > 0L)
    if (!length(cols)) return(NULL)
    data.frame(
      chain = i,
      column = colnames(chains[[i]])[cols],
      n_nonfinite = unname(colSums(nf[, cols, drop = FALSE])),
      first_draw = vapply(cols, function(j) which(nf[, j])[1L], integer(1)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(
      chain = integer(0), column = character(0),
      n_nonfinite = numeric(0), first_draw = integer(0)
    )
  }
  rownames(out) <- NULL
  out
}

# Coerce BAMLSS posterior samples to a matrix
posterior_sample_matrix <- function(model) {
  do.call(rbind, extract_sample_chains(model))
}

# Match a BAMLSS parameter name to a posterior sample column
match_sample_parameter <- function(sample_names, parameter_name) {
  j <- match(parameter_name, sample_names)
  if (!is.na(j)) return(j)
  hit <- which(endsWith(sample_names, parameter_name))
  if (length(hit) == 1L) return(hit)
  NA_integer_
}

# Compute split R-hat for one posterior parameter
split_rhat <- function(chains, parameter_name) {
  if (length(chains) < 2L) return(NA_real_)
  z <- lapply(chains, function(mat) {
    j <- match_sample_parameter(colnames(mat), parameter_name)
    if (is.na(j)) return(numeric(0))
    as.numeric(mat[, j])
  })
  if (any(lengths(z) < 8L) ||
      any(vapply(z, function(x) any(!is.finite(x)), logical(1)))) {
    return(NA_real_)
  }
  n <- min(lengths(z))
  n <- n - (n %% 2L)
  if (n < 8L) return(NA_real_)
  
  h <- n / 2L
  split_chains <- unlist(lapply(z, function(x) {
    x <- tail(x, n)
    list(x[seq_len(h)], x[(h + 1L):n])
  }), recursive = FALSE)
  
  means <- vapply(split_chains, mean, numeric(1))
  vars <- vapply(split_chains, stats::var, numeric(1))
  W <- mean(vars)
  B <- h * stats::var(means)
  tol <- sqrt(.Machine$double.eps)
  
  if (!is.finite(W) || !is.finite(B)) return(NA_real_)
  if (W <= tol) return(NA_real_)
  
  var_plus <- ((h - 1) / h) * W + B / h
  sqrt(max(var_plus / W, 0))
}

# Metadata for the parametric coefficients of a processed frame
build_parametric_metadata <- function(frame) {
  template_names <- names(parameter_template(frame))
  rows <- list()
  pos <- 0L
  
  for (p_name in names(frame$x)) {
    mm <- frame$x[[p_name]]$model.matrix
    if (!is.null(mm)) {
      mm <- as.matrix(mm)
      cn <- colnames(mm)
      if (is.null(cn)) cn <- paste0("X", seq_len(ncol(mm)))
      for (j in seq_along(cn)) {
        leaf <- paste0(p_name, ".p.", cn[j])
        idx <- match(leaf, template_names)
        if (is.na(idx)) {
          hit <- which(endsWith(template_names, leaf))
          if (length(hit) == 1L) idx <- hit
        }
        if (is.na(idx)) next
        pos <- pos + 1L
        rows[[pos]] <- data.frame(
          Parameter = p_name,
          Effect_Type = "structural-fixed",
          Covariate = cn[j],
          Contrast = cn[j],
          Term = cn[j],
          Parameter_Name = template_names[idx],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  meta <- frame$candidate_metadata
  if (!is.null(meta) && length(meta)) {
    for (block in meta) {
      if (identical(block$component, "nonlinear") ||
          !length(block$beta_names)) {
        next
      }
      labels <- block$design_colnames
      if (is.null(labels) || length(labels) != length(block$beta_names)) {
        labels <- paste0("coefficient_", seq_along(block$beta_names))
      }
      for (j in seq_along(block$beta_names)) {
        lab <- labels[j]
        display <- if (length(block$beta_names) == 1L &&
                       (identical(lab, block$variable) ||
                        lab %in% c("1", "(Intercept)"))) {
          block$variable
        } else if (startsWith(lab, block$variable)) {
          lab
        } else {
          paste0(block$variable, ": ", lab)
        }
        pos <- pos + 1L
        rows[[pos]] <- data.frame(
          Parameter = block$parameter,
          Effect_Type = "dynamic-linear",
          Covariate = block$variable,
          Contrast = lab,
          Term = display,
          Parameter_Name = block$beta_names[j],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out <- out[!duplicated(out[c("Parameter", "Parameter_Name")]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Extract posterior summaries for parametric effects
extract_parametric_summaries <- function(
    model, include_intercept = TRUE, strict = TRUE) {
  meta <- attr(model, "parametric_metadata")
  if (is.null(meta) || !nrow(meta)) {
    stop(
      "Model has no parametric metadata; refit with fit_spatial_model().",
      call. = FALSE
    )
  }
  if (!include_intercept) {
    meta <- meta[!grepl("Intercept", meta$Term, fixed = TRUE), , drop = FALSE]
  }
  if (!nrow(meta)) return(meta)
  
  chains <- extract_sample_chains(model)
  mat <- do.call(rbind, chains)
  idx <- vapply(
    meta$Parameter_Name, function(nm) match_sample_parameter(colnames(mat), nm),
    integer(1)
  )
  if (anyNA(idx)) {
    missing <- meta$Parameter_Name[is.na(idx)]
    msg <- paste0(
      "Parametric coefficients missing from the posterior samples: ",
      paste(utils::head(missing, 10L), collapse = ", "),
      if (length(missing) > 10L) " ..." else "", "."
    )
    if (isTRUE(strict)) {
      stop(msg, call. = FALSE)
    } else {
      warning(msg, call. = FALSE)
    }
  }
  
  summ <- lapply(seq_len(nrow(meta)), function(i) {
    if (is.na(idx[i])) {
      return(data.frame(
        Estimate = NA_real_, Sd = NA_real_, Lower_95 = NA_real_,
        Median = NA_real_,
        Upper_95 = NA_real_, ESS = NA_real_, Rhat = NA_real_
      ))
    }
    pname <- meta$Parameter_Name[i]
    z <- as.numeric(mat[, idx[i]])
    q <- stats::quantile(z, c(.025, .5, .975), names = FALSE, type = 8)
    ess_chain <- vapply(chains, function(ch) {
      j <- match_sample_parameter(colnames(ch), pname)
      if (is.na(j)) return(NA_real_)
      effective_sample_size(ch[, j])
    }, numeric(1))
    ess_total <- if (all(is.finite(ess_chain))) sum(ess_chain) else NA_real_
    data.frame(
      Estimate = mean(z),
      Sd = if (length(z) > 1L) stats::sd(z) else NA_real_,
      Lower_95 = q[1L], Median = q[2L], Upper_95 = q[3L],
      ESS = ess_total,
      Rhat = split_rhat(chains, pname)
    )
  })
  cbind(
    meta,
    do.call(rbind, summ),
    Chains = length(chains),
    Posterior_Draws = nrow(mat),
    stringsAsFactors = FALSE
  )
}

# Assess the post-selection posterior sampler
assess_posterior_sampler <- function(
    model, reference_min_chains = 2L, reference_min_draws = 400L,
    ess_reference = 100,
    split_drift_reference = 0.50, rhat_reference = 1.05) {
  
  cf <- extract_parametric_summaries(
    model, include_intercept = TRUE, strict = TRUE
  )
  chains <- extract_sample_chains(model)
  draws_per_chain <- vapply(chains, nrow, integer(1))
  n_chains <- length(chains)
  
  movement <- lapply(cf$Parameter_Name, function(pname) {
    z <- lapply(chains, function(ch) {
      j <- match_sample_parameter(colnames(ch), pname)
      if (is.na(j)) return(numeric(0))
      as.numeric(ch[, j])
    })
    
    unique_counts <- vapply(z, function(x) {
      length(unique(x[is.finite(x)]))
    }, integer(1))
    
    chain_sd <- vapply(z, function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 2L) return(NA_real_)
      stats::sd(x)
    }, numeric(1))
    
    list(
      min_unique = min(unique_counts),
      max_unique = max(unique_counts),
      min_sd = suppressWarnings(
        min(chain_sd[is.finite(chain_sd)], na.rm = TRUE)
      ),
      frozen = any(unique_counts <= 1L)
    )
  })
  
  cf$Minimum_Unique_Draws <- vapply(movement, `[[`, integer(1), "min_unique")
  cf$Maximum_Unique_Draws <- vapply(movement, `[[`, integer(1), "max_unique")
  cf$Frozen_Chain <- vapply(movement, `[[`, logical(1), "frozen")
  
  split_drift_one <- function(z) {
    z <- as.numeric(z)
    z <- z[is.finite(z)]
    n <- length(z)
    if (n < 20L) return(NA_real_)
    h <- floor(n / 2L)
    a <- z[seq_len(h)]
    b <- z[(h + 1L):n]
    den <- stats::sd(z)
    if (!is.finite(den) || den <= sqrt(.Machine$double.eps)) {
      return(NA_real_)
    }
    abs(mean(a) - mean(b)) / den
  }
  
  drift <- vapply(cf$Parameter_Name, function(pname) {
    d <- vapply(chains, function(ch) {
      j <- match_sample_parameter(colnames(ch), pname)
      if (is.na(j)) return(NA_real_)
      split_drift_one(ch[, j])
    }, numeric(1))
    if (any(!is.finite(d))) NA_real_ else max(d)
  }, numeric(1))
  
  ess <- cf$ESS
  rhat <- cf$Rhat
  finite_ess <- ess[is.finite(ess)]
  finite_drift <- drift[is.finite(drift)]
  finite_rhat <- rhat[is.finite(rhat)]
  min_ess_observed <- if (length(finite_ess)) min(finite_ess) else NA_real_
  max_drift_observed <- if (length(finite_drift)) max(finite_drift) else Inf
  max_rhat_observed <- if (n_chains >= 2L && length(finite_rhat)) {
    max(finite_rhat)
  } else {
    NA_real_
  }
  
  rhat_check <- n_chains >= as.integer(reference_min_chains) &&
    length(rhat) > 0L &&
    all(is.finite(rhat)) &&
    max_rhat_observed <= rhat_reference
  
  parameter_movement <- !any(cf$Frozen_Chain)
  
  diagnostic_indicators <- c(
    chain_count = n_chains >= as.integer(reference_min_chains),
    retained_draws_per_chain = length(draws_per_chain) > 0L &&
      all(draws_per_chain >= as.integer(reference_min_draws)),
    effective_sample_size = length(ess) > 0L && all(is.finite(ess)) &&
      min_ess_observed >= ess_reference,
    split_stability = length(drift) > 0L && all(is.finite(drift)) &&
      max_drift_observed <= split_drift_reference,
    split_rhat = rhat_check,
    parameter_movement = parameter_movement
  )
  
  coefficient_diagnostics <- cf
  coefficient_diagnostics$Split_Drift <- drift
  coefficient_diagnostics$ESS_Below_Reference <-
    is.finite(coefficient_diagnostics$ESS) &
    coefficient_diagnostics$ESS < ess_reference
  coefficient_diagnostics$Rhat_Above_Reference <-
    is.finite(coefficient_diagnostics$Rhat) &
    coefficient_diagnostics$Rhat > rhat_reference
  coefficient_diagnostics$Split_Drift_Above_Reference <-
    is.finite(coefficient_diagnostics$Split_Drift) &
    coefficient_diagnostics$Split_Drift > split_drift_reference
  
  coefficient_diagnostics$Diagnostics_Outside_Reference <- vapply(
    seq_len(nrow(coefficient_diagnostics)), function(i) {
      diagnostics <- character(0)
      if (isTRUE(coefficient_diagnostics$Frozen_Chain[i])) {
        diagnostics <- c(diagnostics, "frozen chain")
      }
      if (isTRUE(coefficient_diagnostics$ESS_Below_Reference[i])) {
        diagnostics <- c(diagnostics, "ESS below reference")
      }
      if (isTRUE(coefficient_diagnostics$Rhat_Above_Reference[i])) {
        diagnostics <- c(diagnostics, "R-hat above reference")
      }
      if (isTRUE(coefficient_diagnostics$Split_Drift_Above_Reference[i])) {
        diagnostics <- c(diagnostics, "split drift above reference")
      }
      paste(diagnostics, collapse = "; ")
    }, character(1)
  )
  
  list(
    diagnostics_within_reference = all(diagnostic_indicators),
    diagnostic_indicators = diagnostic_indicators,
    diagnostics_outside_reference =
      names(diagnostic_indicators)[!diagnostic_indicators],
    n_chains = n_chains,
    chain_count_reference = as.integer(reference_min_chains),
    draws_per_chain = draws_per_chain,
    min_draws_per_chain = min(draws_per_chain),
    draws_per_chain_reference = as.integer(reference_min_draws),
    minimum_ess_observed = min_ess_observed,
    ess_reference = ess_reference,
    max_split_drift_observed = max_drift_observed,
    split_drift_reference = split_drift_reference,
    max_rhat_observed = max_rhat_observed,
    rhat_reference = rhat_reference,
    coefficient_diagnostics = coefficient_diagnostics
  )
}

# Posterior chains in saved segments, so a killed run resumes from the
# last full segment of each chain instead of restarting
sample_chains_resumable <- function(
    run_one, start, n_draws, n_burn, thin, segment_iter, chains, cores,
    seed, stem, key) {
  seg_iter <- max(1L, segment_iter %/% thin) * thin
  chain_file <- function(i) sprintf("%s_chain%d.rds", stem, i)
  fit_file <- paste0(stem, "_fit.rds")
  load_fit <- function() {
    fit <- if (file.exists(fit_file)) {
      tryCatch(readRDS(fit_file), error = function(e) NULL)
    }
    if (is.null(fit) || !identical(attr(fit, "resume_key"), key)) NULL else fit
  }
  
  load_chain <- function(i) {
    obj <- if (file.exists(chain_file(i))) {
      tryCatch(readRDS(chain_file(i)), error = function(e) NULL)
    }
    if (is.null(obj) || !identical(obj$key, key)) NULL else obj
  }
  save_atomic <- function(obj, file) {
    tmp <- paste0(file, ".tmp")
    saveRDS(obj, tmp)
    file.rename(tmp, file)
  }
  is_done <- function(i) {
    obj <- load_chain(i)
    !is.null(obj) && NROW(obj$draws) >= n_draws
  }
  
  # The chain continues from the last row of the raw sampler output, which
  # is in the sampler's own parameter names. Summary columns are dropped.
  next_state <- function(raw) {
    if (is.null(raw)) stop("Sampler output was not captured.", call. = FALSE)
    raw <- as.matrix(raw)
    cols <- !grepl("^(logLik|logPost|DIC|pd|alpha)$|\\.edf$", colnames(raw))
    state <- raw[nrow(raw), cols]
    if (!length(state) || any(!is.finite(state))) {
      stop("Last sampler state has no finite parameters.", call. = FALSE)
    }
    state
  }
  
  run_chain <- function(i) {
    obj <- load_chain(i)
    if (is.null(obj)) {
      obj <- list(key = key, burned = 0L, draws = NULL, state = start)
    }
    while (obj$burned < n_burn || NROW(obj$draws) < n_draws) {
      burning <- obj$burned < n_burn
      n_it <- if (burning) {
        min(seg_iter, ceiling((n_burn - obj$burned) / thin) * thin)
      } else {
        min(seg_iter, (n_draws - NROW(obj$draws)) * thin)
      }
      if (!is.null(seed)) {
        set.seed(seed + 1000L * i + obj$burned + NROW(obj$draws))
      }
      res <- run_one(obj$state, n_it)
      z <- extract_sample_chains(res$fit)[[1L]]
      if (!nrow(z)) stop("A segment returned no draws.", call. = FALSE)
      state <- next_state(res$raw)
      if (i == 1L && is.null(load_fit())) {
        res$fit$samples <- NULL
        attr(res$fit, "resume_key") <- key
        save_atomic(res$fit, fit_file)
      }
      if (burning) {
        obj$burned <- obj$burned + n_it
      } else {
        obj$draws <- rbind(obj$draws, z)
      }
      obj$state <- state
      save_atomic(obj, chain_file(i))
    }
    invisible(TRUE)
  }
  
  todo <- which(!vapply(seq_len(chains), is_done, logical(1)))
  
  # Two short runs, the second started from the first one's last state,
  # so a name or format mismatch stops the job in minutes, not hours
  if (length(todo) && is.null(load_chain(todo[1L]))) {
    first <- run_one(start, 2L * thin)
    extract_sample_chains(first$fit)
    second <- run_one(next_state(first$raw), 2L * thin)
    extract_sample_chains(second$fit)
    message("  Resume check passed.")
  }
  
  if (length(todo) > 1L && cores > 1L) {
    parallel::mclapply(
      todo,
      function(i) {
        out <- try(run_chain(i), silent = TRUE)
        if (inherits(out, "try-error")) message("Chain ", i, ": ", out)
        out
      },
      mc.cores = min(cores, length(todo)), mc.preschedule = FALSE
    )
    todo <- which(!vapply(seq_len(chains), is_done, logical(1)))
    if (length(todo)) {
      message(
        "  Chain(s) ", paste(todo, collapse = ", "),
        " unfinished; resuming them one at a time."
      )
    }
  }
  for (i in todo) run_chain(i)
  
  fit <- load_fit()
  if (is.null(fit)) {
    fit <- run_one(start, 2L * thin)$fit
    fit$samples <- NULL
    attr(fit, "resume_key") <- key
    save_atomic(fit, fit_file)
  }
  list(
    fit = fit,
    samples = lapply(seq_len(chains), function(i) {
      coda::mcmc(utils::head(load_chain(i)$draws, n_draws))
    })
  )
}

# Fit selected Bayesian spatial model
fit_spatial_model <- function(
    label,
    formula_list,
    family_name,
    data,
    offset = NULL,
    formula_base = NULL,
    candidate_inclusion = NULL,
    candidate_covariates = NULL,
    block_selected_linear = FALSE,
    strict_sampler_validation = FALSE,
    c_linear = 2.0,
    c0 = 316.227766,
    cplus = 316.227766,
    a_lambda = 1,
    b_lambda = 1,
    start = NULL,
    optimizer = NULL,
    sampler = NULL,
    n_iter = 6000,
    n_burn = 2000,
    thin = 10,
    chains = 4L,
    cores = NULL,
    chain_seed = NULL,
    reference_min_chains = 2L,
    reference_min_draws = 400L,
    ess_reference = 100,
    split_drift_reference = 0.50,
    rhat_reference = 1.05,
    chain_checkpoint = NULL,
    segment_iter = 2000L) {
  
  fam_obj <- get_model_family(family_name)
  
  transform_fun <- NULL
  posterior_frame <- NULL
  source_selected_frame <- NULL
  
  if (!is.null(candidate_inclusion)) {
    
    if (is.null(formula_base)) {
      stop(
        "`formula_base` is required when `candidate_inclusion` is supplied.",
        call. = FALSE
      )
    }
    
    candidate_inclusion <- as.matrix(candidate_inclusion)
    
    if (!is.logical(candidate_inclusion)) {
      storage.mode(candidate_inclusion) <- "logical"
    }
    
    expected_rows <- if (!is.null(candidate_covariates)) {
      candidate_covariates
    } else {
      rownames(candidate_inclusion)
    }
    
    if (
      is.null(rownames(candidate_inclusion)) ||
      is.null(colnames(candidate_inclusion)) ||
      !identical(rownames(candidate_inclusion), expected_rows) ||
      !identical(colnames(candidate_inclusion), fam_obj$names)
    ) {
      stop(
        "Selected allocation matrix has incompatible row/column names.",
        call. = FALSE
      )
    }
    
    if (isTRUE(block_selected_linear)) {
      
      formula_list <- build_selected_posterior_formula(
        inclusion = candidate_inclusion,
        formula_base = formula_base,
        data = data,
        fam_obj = fam_obj,
        a_lambda = a_lambda,
        b_lambda = b_lambda
      )
      
      transform_fun <- make_selected_block_model_transform(
        candidate_inclusion = candidate_inclusion,
        formula_base = formula_base,
        data = data,
        offset = offset,
        c_linear = c_linear,
        c0 = c0,
        cplus = cplus,
        a_lambda = a_lambda,
        b_lambda = b_lambda
      )
      
      source_runtime <- new_runtime_state(
        cache_frames = FALSE,
        cache_blocks = TRUE
      )
      
      source_selected_frame <- get_model_frame(
        inclusion = candidate_inclusion,
        formula_base = formula_base,
        data = data,
        fam_obj = fam_obj,
        offset = offset,
        runtime = source_runtime,
        c_linear = c_linear,
        a_lambda = a_lambda,
        b_lambda = b_lambda,
        c0 = c0,
        cplus = cplus
      )
      
      parsed_block_formula <- parse_distributional_formula(
        formula_list,
        fam_obj
      )
      
      posterior_frame <- bamlss::bamlss.frame(
        parsed_block_formula,
        data = data,
        family = fam_obj,
        offset = offset
      )
      
      posterior_frame <- transform_fun(posterior_frame)
      
    } else {
      
      transform_fun <- make_selected_model_transform(
        candidate_inclusion = candidate_inclusion,
        formula_base = formula_base,
        data = data,
        offset = offset,
        c_linear = c_linear,
        c0 = c0,
        cplus = cplus,
        a_lambda = a_lambda,
        b_lambda = b_lambda
      )
    }
  }
  
  if (identical(optimizer, FALSE) && is.null(start)) {
    stop("optimizer = FALSE needs a finite `start` vector.", call. = FALSE)
  }
  
  if (!is.null(start)) {
    start <- flatten_parameters(start)
    if (!is.numeric(start) || is.null(names(start)) ||
        any(!nzchar(names(start))) || any(!is.finite(start))) {
      stop(
        "Post-selection sampler start must be a finite named numeric vector.",
        call. = FALSE
      )
    }
  }
  
  if (isTRUE(block_selected_linear) && !is.null(start)) {
    if (is.null(source_selected_frame) || is.null(posterior_frame)) {
      stop("Blocked posterior frames were not initialized.", call. = FALSE)
    }
    
    start <- project_start_to_blocked_posterior(
      start = start,
      source_frame = source_selected_frame,
      target_frame = posterior_frame,
      fam_obj = fam_obj,
      a_lambda = a_lambda,
      b_lambda = b_lambda
    )
    
    blocking_validation <- attr(start, "blocking_validation")
    
    message(
      sprintf(
        paste0(
          "  Blocked posterior equivalence: ",
          "max|delta eta|=%.3g; ",
          "delta structural prior=%.3g; ",
          "delta candidate prior=%.3g; ",
          "delta logLik=%.3g."
        ),
        blocking_validation$maximum_eta_difference,
        blocking_validation$structural_prior_difference,
        blocking_validation$candidate_prior_difference,
        blocking_validation$loglik_difference
      )
    )
  }
  
  if (identical(optimizer, FALSE) && is.null(sampler)) {
    sampler <- bamlss::sam_GMCMC
  }
  if (is.function(sampler)) sampler <- guard_sampler(sampler)
  chains <- as.integer(chains)
  if (!is.finite(chains) || chains < 1L) {
    stop("`chains` must be a positive integer.", call. = FALSE)
  }
  if (is.null(cores)) cores <- chains
  cores <- as.integer(cores)
  if (!is.finite(cores) || cores < 1L) {
    stop("`cores` must be a positive integer or NULL.", call. = FALSE)
  }
  cores <- min(cores, chains)
  if (!is.null(chain_seed)) {
    chain_seed <- as.integer(chain_seed)[1L]
    if (!is.finite(chain_seed)) {
      stop("`chain_seed` must be NULL or a finite integer.", call. = FALSE)
    }
  }
  
  parsed_formula_list <- parse_distributional_formula(formula_list, fam_obj)
  
  chain_samples <- list()
  fit <- NULL
  completed_chains <- 0L
  batch_id <- 0L
  
  resumable <- !is.null(chain_checkpoint)
  if (resumable) {
    if (is.null(start)) {
      stop("Resumable sampling needs a named `start` vector.", call. = FALSE)
    }
    segment_iter <- as.integer(segment_iter)[1L]
    if (!is.finite(segment_iter) || segment_iter < thin) {
      stop("`segment_iter` must be at least `thin`.", call. = FALSE)
    }
    resume_key <- paste(
      label, n_iter, n_burn, thin, chains, chain_seed,
      paste(names(start), collapse = ","), signif(sum(start), 8),
      paste(deparse(formula_list), collapse = ""),
      format(getOption("rbst.structural_sd", "bamlss default")),
      format(getOption("rbst.mrf_theta", "bamlss default")),
      sep = "|"
    )
    resumed <- sample_chains_resumable(
      run_one = function(st, n_it) {
        raw <- NULL
        keep_raw <- function(...) {
          raw <<- sampler(...)
          raw
        }
        fit <- bamlss::bamlss(
          parsed_formula_list, family = fam_obj, data = data,
          offset = offset, transform = transform_fun,
          start = to_sampler_start(st),
          optimizer = optimizer, sampler = keep_raw, n.iter = n_it + 1L,
          burnin = 1L, thin = thin, cores = 1L, combine = FALSE,
          results = FALSE
        )
        list(fit = fit, raw = raw)
      },
      start = start, n_draws = (n_iter - n_burn) %/% thin,
      n_burn = n_burn, thin = thin, segment_iter = segment_iter,
      chains = chains, cores = cores, seed = chain_seed,
      stem = chain_checkpoint, key = resume_key
    )
    fit <- resumed$fit
    chain_samples <- resumed$samples
  }
  
  run_chains <- function(k) {
    bamlss::bamlss(
      parsed_formula_list,
      family = fam_obj,
      data = data,
      offset = offset,
      transform = transform_fun,
      start = to_sampler_start(start),
      optimizer = optimizer,
      sampler = sampler,
      n.iter = n_iter,
      burnin = n_burn,
      thin = thin,
      cores = k,
      combine = FALSE,
      results = FALSE
    )
  }
  
  # A forked chain killed (usually for memory) returns nothing; the chains of
  # the batch then run one at a time here, which needs less memory and shows
  # the real error.
  sample_batch <- function(k) {
    invisible(gc(verbose = FALSE))
    out <- tryCatch(
      {
        fit_k <- run_chains(k)
        list(fit = fit_k, samples = extract_sample_chains(fit_k))
      },
      error = function(e) e
    )
    if (!inherits(out, "error")) return(out)
    if (k < 2L) stop(conditionMessage(out), call. = FALSE)
    message(
      "  Parallel chains failed (", conditionMessage(out),
      "); running them one at a time."
    )
    parts <- lapply(seq_len(k), function(i) {
      if (!is.null(chain_seed)) {
        set.seed(chain_seed + completed_chains + i - 1L)
      }
      fit_i <- run_chains(1L)
      list(fit = fit_i, samples = extract_sample_chains(fit_i))
    })
    list(
      fit = parts[[1L]]$fit,
      samples = unlist(lapply(parts, `[[`, "samples"), recursive = FALSE)
    )
  }
  
  while (!resumable && completed_chains < chains) {
    batch_id <- batch_id + 1L
    batch_size <- min(cores, chains - completed_chains)
    if (!is.null(chain_seed)) set.seed(chain_seed + completed_chains)
    
    batch <- sample_batch(batch_size)
    batch_fit <- batch$fit
    batch_samples <- batch$samples
    if (length(batch_samples) != batch_size) {
      stop(
        "Posterior sampler returned ", length(batch_samples),
        " chain(s) for a requested batch of ", batch_size, ".",
        call. = FALSE
      )
    }
    
    chain_samples <- c(chain_samples, lapply(batch_samples, coda::mcmc))
    if (is.null(fit)) {
      fit <- batch_fit
      fit$samples <- NULL
    }
    
    completed_chains <- completed_chains + batch_size
    rm(batch, batch_fit, batch_samples)
    invisible(gc(verbose = FALSE))
  }
  
  if (length(chain_samples) != chains) {
    stop("Posterior chain assembly failed.", call. = FALSE)
  }
  fit$samples <- coda::mcmc.list(chain_samples)
  
  # logLik, DIC and pd from all assembled draws, not from the template fit
  # (one segment or one batch of chains)
  ms <- tryCatch({
    smp <- lapply(fit$samples, function(z) {
      z <- as.matrix(z)
      coda::mcmc(z[, !colnames(z) %in% c("DIC", "pd"), drop = FALSE])
    })
    bamlss::samplestats(
      coda::mcmc.list(smp), x = fit$x, y = fit$y, family = fit$family
    )
  }, error = function(e) {
    warning("Falha no recalculo de logLik/DIC/pd: ", conditionMessage(e),
            call. = FALSE)
    NULL
  })
  if (is.list(ms) && length(ms)) {
    if (is.null(fit$model.stats)) fit$model.stats <- list()
    fit$model.stats$sampler[names(ms)] <- ms
  } else if (!is.null(fit$model.stats$sampler)) {
    # stale values come from one segment only, so they are not reported
    st <- fit$model.stats$sampler
    for (v in intersect(c("DIC", "pd", "logLik"), names(st))) st[[v]] <- NA
    fit$model.stats$sampler <- st
  }
  
  parametric_meta <- NULL
  if (!is.null(candidate_inclusion)) {
    
    parametric_frame <- if (isTRUE(block_selected_linear)) {
      if (is.null(posterior_frame)) {
        stop("Blocked posterior reporting frame is unavailable.", call. = FALSE)
      }
      posterior_frame
    } else {
      rt_report <- new_runtime_state(cache_frames = FALSE, cache_blocks = TRUE)
      get_model_frame(
        candidate_inclusion, formula_base, data, fam_obj, offset, rt_report,
        c_linear = c_linear, a_lambda = a_lambda, b_lambda = b_lambda,
        c0 = c0, cplus = cplus
      )
    }
    
    parametric_meta <- build_parametric_metadata(parametric_frame)
    if (!nrow(parametric_meta)) {
      warning(
        "Selected model contains no parametric coefficients.", call. = FALSE
      )
    }
    
    smat_check <- posterior_sample_matrix(fit)
    idx_check <- vapply(
      parametric_meta$Parameter_Name,
      function(nm) match_sample_parameter(colnames(smat_check), nm),
      integer(1)
    )
    if (anyNA(idx_check)) {
      stop(
        "Posterior samples are missing parametric coefficient(s): ",
        paste(
          parametric_meta$Parameter_Name[is.na(idx_check)], collapse = ", "
        ),
        call. = FALSE
      )
    }
  }
  
  attr(fit, "parametric_metadata") <- parametric_meta
  
  sampler_diagnostics <- NULL
  if (identical(optimizer, FALSE) && !is.null(parametric_meta)) {
    sampler_diagnostics <- assess_posterior_sampler(
      fit,
      reference_min_chains = reference_min_chains,
      reference_min_draws = reference_min_draws,
      ess_reference = ess_reference,
      split_drift_reference = split_drift_reference,
      rhat_reference = rhat_reference
    )
    attr(fit, "sampler_diagnostics") <- sampler_diagnostics
    
    frozen_detected <- (
      !is.null(sampler_diagnostics) &&
        !is.null(sampler_diagnostics$coefficient_diagnostics) &&
        any(sampler_diagnostics$coefficient_diagnostics$Frozen_Chain %in% TRUE)
    )
    
    if (isTRUE(strict_sampler_validation) && frozen_detected) {
      stop(
        "Posterior sampling produced at least one frozen coefficient chain. ",
        "Goodness-of-fit, information criteria, posterior summaries, and ",
        "predictive diagnostics must not be computed from this fit.",
        call. = FALSE
      )
    }
    
    if (!isTRUE(sampler_diagnostics$diagnostics_within_reference)) {
      warning(
        "Posterior diagnostics outside their reference values: ",
        paste(
          sampler_diagnostics$diagnostics_outside_reference, collapse = ", "
        ),
        ". Summaries are kept; read them with the ESS, R-hat and drift checks.",
        call. = FALSE
      )
    }
  }
  
  attr(fit, "sampler_settings") <- list(
    label = label,
    candidate_inclusion = candidate_inclusion,
    block_selected_linear = isTRUE(block_selected_linear),
    c_linear = c_linear,
    c0 = c0,
    cplus = cplus,
    a_lambda = a_lambda,
    b_lambda = b_lambda,
    reference_min_chains = reference_min_chains,
    reference_min_draws = reference_min_draws,
    ess_reference = ess_reference,
    split_drift_reference = split_drift_reference,
    rhat_reference = rhat_reference,
    requested_chains = chains,
    concurrent_chains = cores,
    chain_seed = chain_seed,
    diagnostics_within_reference = if (is.null(sampler_diagnostics)) NA else {
      isTRUE(sampler_diagnostics$diagnostics_within_reference)
    },
    optimizer = if (identical(optimizer, FALSE)) {
      "disabled"
    } else {
      "enabled/default"
    },
    sampler = if (identical(optimizer, FALSE)) {
      "sam_GMCMC from retained model-search state"
    } else {
      "default/requested"
    }
  )
  fit
}

# Create model summary table
create_model_summary_table <- function(
    model, include_intercept = TRUE, base_family = "serif") {
  df_fixed <- extract_parametric_summaries(
    model, include_intercept = include_intercept, strict = TRUE
  )
  if (!nrow(df_fixed)) {
    return(
      gt::gt(data.frame(Message = "No parametric coefficients available."))
    )
  }
  
  cols <- c(
    "Parameter", "Effect_Type", "Term", "Estimate", "Sd", "Lower_95",
    "Median", "Upper_95", "ESS", "Rhat"
  )
  shown <- df_fixed[, cols, drop = FALSE]
  shown$Term <- stringr::str_remove(shown$Term, "\\.{3}\\d+")
  shown$Term <- stringr::str_remove(shown$Term, "^b\\.")
  shown$Term <- stringr::str_remove(shown$Term, "TRUE$")
  shown$Term <- format_term_label(shown$Term)
  shown$Effect_Type <- ifelse(
    shown$Effect_Type == "dynamic-linear", "Selected", "Structural"
  )
  
  gt::gt(shown, groupname_col = "Parameter") |>
    gt::tab_header(
      title = "Conditional Post-Selection Coefficient Summary",
      subtitle = paste("Distribution family:", model$family$family)
    ) |>
    gt::fmt_number(
      columns = c(Estimate, Sd, Lower_95, Median, Upper_95, ESS, Rhat),
      decimals = 4
    ) |>
    gt::cols_label(
      Effect_Type = "Source", Term = "Covariate", Estimate = "Mean", Sd = "SD",
      Lower_95 = "2.5% CrI", Median = "Median", Upper_95 = "97.5% CrI",
      ESS = "ESS", Rhat = "R-hat"
    ) |>
    gt::tab_options(
      row_group.font.weight = "bold", column_labels.font.weight = "bold",
      table.font.names = base_family
    )
}

# Predictive PMF on 0..z_max (no renormalisation); reports the tail mass
predictive_pmf_matrix <- function(
    fam, preds, y_obs, tail_tol = 1e-10, zmax_cap = 100000L) {
  n <- length(y_obs)
  y_top <- max(y_obs)
  z_max <- max(1L, y_top)
  zmax_cap <- min(as.integer(zmax_cap), as.integer(max(1000L, y_top * 3L)))
  
  while (z_max < zmax_cap) {
    cdf_z <- family_cdf(family = fam, q = rep(z_max, n), parameters = preds)
    if (any(!is.finite(cdf_z)) || all(1 - cdf_z <= tail_tol)) break
    z_max <- as.integer(min(zmax_cap, max(z_max + 1L, ceiling(z_max * 2))))
  }
  
  cdf_z <- family_cdf(family = fam, q = rep(z_max, n), parameters = preds)
  support <- 0:z_max
  pmf <- matrix(0, nrow = n, ncol = length(support))
  for (j in seq_along(support)) {
    val <- as.numeric(fam$d(rep(support[j], n), preds, log = FALSE))
    val[!is.finite(val) | val < 0] <- 0
    pmf[, j] <- val
  }
  
  list(
    pmf = pmf, support = support,
    tail = ifelse(is.finite(cdf_z), pmax(0, 1 - cdf_z), NA_real_),
    gap = max(abs(rowSums(pmf) - cdf_z), na.rm = TRUE)
  )
}

# Compute predictive scoring rules (plug-in at the posterior-mean parameters)
compute_predictive_scores <- function(
    model, data, y_col, tail_tol = 1e-10, zmax_cap = 100000L) {
  out <- data.frame(
    QS = NA_real_, LS = NA_real_, SPS = NA_real_, CRPS = NA_real_,
    RPS_K = NA_real_, TailMass = NA_real_, PmfCdfGap = NA_real_
  )
  y_obs <- as.integer(data[[y_col]])
  n <- length(y_obs)
  if (!n || any(!is.finite(y_obs)) || any(y_obs < 0L)) return(out)
  
  preds <- predict_bamlss_parameters(
    model = model, newdata = NULL, drop = FALSE
  )
  if (any(vapply(preds, length, integer(1)) != n)) {
    stop(
      "In-sample predictions are not aligned with the evaluation data.",
      call. = FALSE
    )
  }
  
  fam <- model$family
  logp <- tryCatch(
    as.numeric(fam$d(y_obs, preds, log = TRUE)),
    error = function(e) rep(NA_real_, n)
  )
  out$LS <- if (any(is.finite(logp))) mean(logp[is.finite(logp)]) else NA_real_
  if (is.null(fam$p)) return(out)
  
  pm <- predictive_pmf_matrix(fam, preds, y_obs, tail_tol, zmax_cap)
  pmf <- pm$pmf
  p_obs <- pmf[cbind(seq_len(n), y_obs + 1L)]
  sum_pmf2 <- rowSums(pmf^2)
  out$QS <- mean(2 * p_obs - sum_pmf2)
  out$SPS <- mean(p_obs / sqrt(pmax(sum_pmf2, .Machine$double.xmin)))
  
  # CRPS only if the omitted tail is negligible; RPS_K truncates at max(y)
  # (counting measure on 0..K-1) and is always finite.
  cdf <- t(apply(pmf, 1L, cumsum))
  sq <- (cdf - outer(y_obs, pm$support, "<=")) ^ 2
  k_top <- max(y_obs)
  out$RPS_K <- mean(rowSums(sq[, pm$support < k_top, drop = FALSE]))
  out$TailMass <- max(pm$tail)
  out$PmfCdfGap <- pm$gap
  if (is.finite(out$TailMass) && out$TailMass <= tail_tol) {
    out$CRPS <- mean(rowSums(sq))
  }
  out
}

# Exceedance Brier scores and central-interval score (plug-in)
compute_threshold_scores <- function(
    model, data, y_col, thresholds = NULL, alpha = 0.1,
    tail_tol = 1e-10, zmax_cap = 100000L) {
  y_obs <- as.integer(data[[y_col]])
  n <- length(y_obs)
  if (is.null(thresholds)) thresholds <- default_exceedance_thresholds(y_obs)
  preds <- predict_bamlss_parameters(
    model = model, newdata = NULL, drop = FALSE
  )
  pm <- predictive_pmf_matrix(model$family, preds, y_obs, tail_tol, zmax_cap)
  cdf <- t(apply(pm$pmf, 1L, cumsum))
  
  thresholds <- thresholds[thresholds < max(pm$support)]
  exceed <- do.call(rbind, lapply(thresholds, function(k) {
    p_hat <- pmin(pmax(1 - cdf[, k + 1L], 0), 1)
    hit <- as.numeric(y_obs > k)
    data.frame(
      Threshold = k, ObsFreq = mean(hit), PredFreq = mean(p_hat),
      Brier = mean((p_hat - hit)^2)
    )
  }))
  
  lo <- rowSums(cdf < alpha / 2)
  up <- rowSums(cdf < 1 - alpha / 2)
  undefined <- up > max(pm$support)
  width <- up - lo
  is_loss <- width + (2 / alpha) * (lo - y_obs) * (y_obs < lo) +
    (2 / alpha) * (y_obs - up) * (y_obs > up)
  interval <- data.frame(
    Alpha = alpha, N_Undefined = sum(undefined),
    IntervalScore = if (any(undefined)) NA_real_ else mean(is_loss),
    Coverage = if (any(undefined)) {
      NA_real_
    } else {
      mean(y_obs >= lo & y_obs <= up)
    },
    MeanWidth = if (any(undefined)) NA_real_ else mean(width)
  )
  list(exceedance = exceed, interval = interval)
}

# Watson uniformity test with Monte Carlo calibration
watson_uniformity_test <- function(u, n_sim = 1999L, seed = 20260914L) {
  u <- sort(as.numeric(u[is.finite(u)]))
  n <- length(u)
  if (n < 4L) return(list(statistic = NA_real_, p.value = NA_real_))
  stat_fun <- function(z) {
    z <- sort(z)
    nn <- length(z)
    a <- z - (2 * seq_len(nn) - 1) / (2 * nn)
    sum((a - mean(a))^2) + 1 / (12 * nn)
  }
  observed <- stat_fun(u)
  n_sim <- as.integer(n_sim)
  if (n_sim < 1L) return(list(statistic = observed, p.value = NA_real_))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  simulated <- replicate(n_sim, stat_fun(stats::runif(n)))
  list(
    statistic = observed,
    p.value = (1 + sum(simulated >= observed)) / (n_sim + 1)
  )
}

# Zamanzade TB2 uniformity test with Monte Carlo calibration
tb2_uniformity_test <- function(
    u, m = NULL, n_sim = 1999L, seed = 20260915L) {
  
  u <- sort(as.numeric(u[is.finite(u)]))
  u <- u[u >= 0 & u <= 1]
  n <- length(u)
  
  if (n < 4L) {
    return(list(statistic = NA_real_, p.value = NA_real_, m = NA_integer_))
  }
  
  if (is.null(m)) {
    m <- max(1L, floor(sqrt(n)))
  }
  m <- as.integer(m)[1L]
  
  if (!is.finite(m) || m < 1L || m >= n) {
    stop("`m` must be an integer in 1, ..., n - 1.", call. = FALSE)
  }
  
  stat_fun <- function(z) {
    z <- sort(as.numeric(z))
    nn <- length(z)
    
    if (
      nn != n ||
      any(!is.finite(z)) ||
      any(z < 0) ||
      any(z > 1)
    ) {
      return(NA_real_)
    }
    
    i <- seq_len(nn)
    lo_rank <- pmax(i - m, 0L)
    hi_rank <- pmin(i + m, nn + 1L)
    
    z_ext <- c(0, z, 1)
    lower <- z_ext[lo_rank + 1L]
    upper <- z_ext[hi_rank + 1L]
    spacing <- upper - lower
    
    cdf_lower <- ifelse(lo_rank == 0L, 0, lo_rank / nn)
    cdf_upper <- ifelse(hi_rank == nn + 1L, 1, hi_rank / nn)
    cdf_increment <- cdf_upper - cdf_lower
    
    valid <- is.finite(spacing) &
      is.finite(cdf_increment) &
      spacing > 0 &
      cdf_increment > 0
    
    if (!all(valid)) {
      return(NA_real_)
    }
    
    weight <- cdf_increment /
      sum(cdf_increment)
    
    sum(log(spacing / cdf_increment) * weight)
  }
  
  observed <- stat_fun(u)
  
  if (!is.finite(observed)) {
    return(list(statistic = NA_real_, p.value = NA_real_, m = m))
  }
  
  n_sim <- as.integer(n_sim)[1L]
  
  if (!is.finite(n_sim) || n_sim < 1L) {
    return(list(statistic = observed, p.value = NA_real_, m = m))
  }
  
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  
  set.seed(as.integer(seed))
  
  simulated <- replicate(n_sim, stat_fun(stats::runif(n)))
  
  simulated <- simulated[is.finite(simulated)]
  
  p_value <- if (length(simulated)) {
    (1 + sum(simulated <= observed)) /
      (length(simulated) + 1)
  } else {
    NA_real_
  }
  
  list(statistic = observed, p.value = p_value, m = m)
}

# Compute model diagnostics
compute_all_metrics <- function(
    model_obj, label, data, y_col, space_col, listw_obj = NULL,
    residual_seed = 20260913L, watson_mc = 1999L,
    tb2_mc = 1999L, tb2_m = NULL) {
  
  rqr <- tryCatch(
    quantile_residuals(
      model = model_obj, data = data, y_col = y_col, seed = residual_seed
    ),
    error = function(e) {
      warning(
        "Residual diagnostics failed for `", label, "`: ", conditionMessage(e),
        call. = FALSE
      )
      rep(NA_real_, nrow(data))
    }
  )
  
  finite_rqr <- is.finite(rqr)
  n_rqr <- sum(finite_rqr)
  moran_i <- NA_real_
  moran_p <- NA_real_
  
  if (!is.null(listw_obj) && n_rqr > 0L) {
    df_agg <- data.frame(space = data[[space_col]], resid = rqr) |>
      dplyr::group_by(space) |>
      dplyr::summarise(
        mean_resid = if (any(is.finite(resid))) {
          mean(resid[is.finite(resid)])
        } else {
          NA_real_
        },
        .groups = "drop"
      ) |>
      dplyr::arrange(match(space, attributes(listw_obj)$region.id))
    mt <- tryCatch(
      spdep::moran.test(
        df_agg$mean_resid, listw_obj,
        zero.policy = TRUE, na.action = stats::na.omit
      ),
      error = function(e) NULL
    )
    if (!is.null(mt)) {
      moran_i <- unname(mt$estimate[[1L]])
      moran_p <- mt$p.value
    }
  }
  
  u <- stats::pnorm(rqr[finite_rqr])
  u <- u[is.finite(u)]
  
  extract_htest <- function(expr) {
    ans <- tryCatch(expr, error = function(e) NULL)
    if (is.null(ans)) return(c(statistic = NA_real_, p.value = NA_real_))
    stat <- if (!is.null(ans$statistic)) {
      as.numeric(ans$statistic)[1L]
    } else {
      NA_real_
    }
    pval <- if (!is.null(ans$p.value)) as.numeric(ans$p.value)[1L] else NA_real_
    c(statistic = stat, p.value = pval)
  }
  
  ks <- if (length(u) >= 4L) {
    extract_htest(suppressWarnings(stats::ks.test(u, "punif")))
  } else {
    c(statistic = NA_real_, p.value = NA_real_)
  }
  ad <- if (length(u) >= 4L) {
    extract_htest(goftest::ad.test(u, null = "punif"))
  } else {
    c(statistic = NA_real_, p.value = NA_real_)
  }
  ns <- if (length(u) >= 4L) {
    extract_htest(ddst::ddst.uniform.test(u, compute.p = TRUE))
  } else {
    c(statistic = NA_real_, p.value = NA_real_)
  }
  watson <- if (length(u) >= 4L) {
    watson_uniformity_test(u, n_sim = watson_mc, seed = residual_seed + 1L)
  } else {
    list(statistic = NA_real_, p.value = NA_real_)
  }
  tb2 <- if (length(u) >= 4L) {
    tb2_uniformity_test(u, m = tb2_m, n_sim = tb2_mc, seed = residual_seed + 2L)
  } else {
    list(statistic = NA_real_, p.value = NA_real_, m = NA_integer_)
  }
  
  scores <- tryCatch(
    compute_predictive_scores(model_obj, data, y_col),
    error = function(e) data.frame(
      QS = NA_real_, LS = NA_real_, SPS = NA_real_, CRPS = NA_real_,
      RPS_K = NA_real_, TailMass = NA_real_, PmfCdfGap = NA_real_
    )
  )
  
  data.frame(
    Model = label,
    WAIC = tryCatch(
      bamlss::WAIC(model_obj)$WAIC2, error = function(e) NA_real_
    ),
    DIC = tryCatch(bamlss::DIC(model_obj)$DIC, error = function(e) NA_real_),
    LogLik = tryCatch(
      as.numeric(stats::logLik(model_obj)), error = function(e) NA_real_
    ),
    RQR_Finite_Pct = if (length(rqr)) 100 * n_rqr / length(rqr) else NA_real_,
    Extreme_RQR_Pct = if (n_rqr > 0L) {
      100 * sum(abs(rqr[finite_rqr]) > 3) / n_rqr
    } else {
      NA_real_
    },
    Moran_I = moran_i,
    Moran_P = moran_p,
    KS = unname(ks["statistic"]),
    KS_P = unname(ks["p.value"]),
    AD = unname(ad["statistic"]),
    AD_P = unname(ad["p.value"]),
    Watson = unname(watson$statistic),
    Watson_P = unname(watson$p.value),
    NS = unname(ns["statistic"]),
    NS_P = unname(ns["p.value"]),
    TB2 = unname(tb2$statistic),
    TB2_P = unname(tb2$p.value),
    TB2_m = unname(tb2$m),
    QS = scores$QS[1L],
    LS = scores$LS[1L],
    SPS = scores$SPS[1L],
    CRPS = scores$CRPS[1L],
    RPS_K = scores$RPS_K[1L],
    TailMass = scores$TailMass[1L],
    PmfCdfGap = scores$PmfCdfGap[1L],
    stringsAsFactors = FALSE
  )
}

# Compute expected values E[Y]
compute_expected_value <- function(
    model,
    newdata = NULL) {
  
  preds <- predict_bamlss_parameters(
    model = model,
    newdata = newdata,
    drop = FALSE
  )
  
  fam <- model$family$family
  
  mu <- if (!is.null(preds$mu)) {
    as.numeric(preds$mu)
  } else if (!is.null(preds$lambda)) {
    as.numeric(preds$lambda)
  } else {
    stop("Predictions contain neither `mu` nor `lambda`.", call. = FALSE)
  }
  
  if (identical(fam, "RSB")) {
    s <- as.numeric(preds$s)
    
    return(ifelse(is.finite(s) & s > 0, Inf, mu))
  }
  
  if (fam %in% c("ZINBI", "ZIPIG")) {
    return((1 - as.numeric(preds$nu)) * mu)
  }
  
  if (identical(fam, "ZANBI")) {
    sigma <- pmax(as.numeric(preds$sigma), .Machine$double.xmin)
    
    ppos <- stats::pnbinom(0, mu = mu, size = 1 / sigma, lower.tail = FALSE)
    
    return((1 - as.numeric(preds$nu)) * mu / pmax(ppos, .Machine$double.xmin))
  }
  
  if (identical(fam, "ZOINB")) {
    sigma <- pmax(as.numeric(preds$sigma), .Machine$double.xmin)
    
    p1 <- stats::dnbinom(1, mu = mu, size = 1 / sigma)
    
    pgt1 <- stats::pnbinom(1, mu = mu, size = 1 / sigma, lower.tail = FALSE)
    
    return(
      (1 - as.numeric(preds$nu)) *
        (
          as.numeric(preds$tau) + (1 - as.numeric(preds$tau)) *
            ((mu - p1) / pmax(pgt1, .Machine$double.xmin))
        )
    )
  }
  
  mu
}

# Compute location / regular-process intensity
compute_location_intensity <- function(
    model,
    newdata = NULL) {
  
  preds <- predict_bamlss_parameters(
    model = model,
    newdata = newdata,
    drop = FALSE
  )
  
  if (!is.null(preds$lambda)) {
    return(as.numeric(preds$lambda))
  }
  
  if (!is.null(preds$mu)) {
    return(as.numeric(preds$mu))
  }
  
  stop("Location parameter missing from BAMLSS predictions.", call. = FALSE)
}

# Plot response distribution
plot_response_distribution <- function(
    y, title = "Response Distribution", base_family = "serif") {
  cols <- spatial_palette()
  y <- as.numeric(y)
  y <- y[is.finite(y) & y >= 0]
  if (!length(y)) {
    stop("`y` contains no finite non-negative observations.", call. = FALSE)
  }
  
  p1 <- ggplot2::ggplot(data.frame(y = y), ggplot2::aes(x = y)) +
    ggplot2::geom_histogram(bins = 30, fill = cols["primary"]) +
    ggplot2::labs(
      title = "(a) Raw counts", x = expression(y[st]), y = "Frequency"
    ) +
    theme_bamlss_spatial(base_family = base_family)
  
  p2 <- ggplot2::ggplot(data.frame(y = y[y > 0]), ggplot2::aes(x = y)) +
    ggplot2::geom_histogram(bins = 30, fill = cols["primary"]) +
    ggplot2::labs(
      title = expression(paste("(b) ", y[st] > 0)), x = expression(y[st]),
      y = NULL) +
    theme_bamlss_spatial(base_family = base_family)
  
  p3 <- ggplot2::ggplot(data.frame(y = y[y > 1]), ggplot2::aes(x = y)) +
    ggplot2::geom_histogram(bins = 30, fill = cols["primary"]) +
    ggplot2::labs(
      title = expression(paste("(c) ", y[st] > 1)), x = expression(y[st]),
      y = NULL) +
    theme_bamlss_spatial(base_family = base_family)
  
  p4 <- ggplot2::ggplot(
    data.frame(y = y[y > 1 & y <= 50]), ggplot2::aes(x = y)
  ) +
    ggplot2::geom_histogram(bins = 25, fill = cols["primary"]) +
    ggplot2::labs(
      title = expression(paste("(d) ", 1 < y[st], " and ", y[st] <= 50)),
      x = expression(y[st]), y = NULL
    ) +
    theme_bamlss_spatial(base_family = base_family)
  
  p5 <- ggplot2::ggplot(data.frame(y = log1p(y)), ggplot2::aes(x = y)) +
    ggplot2::geom_histogram(bins = 30, fill = cols["primary"]) +
    ggplot2::labs(
      title = "(e) Log transformation", x = expression(log(1 + y[st])),
      y = NULL
    ) +
    theme_bamlss_spatial(base_family = base_family)
  
  probs <- seq(0, 1, by = 0.01)
  qd <- data.frame(
    p = probs, q = stats::quantile(y, probs, names = FALSE, type = 8)
  )
  p6 <- ggplot2::ggplot(qd, ggplot2::aes(x = p, y = q)) +
    ggplot2::geom_step(color = cols["primary"]) +
    ggplot2::labs(
      title = "(f) Empirical quantiles", x = "Probability", y = NULL
    ) +
    theme_bamlss_spatial(base_family = base_family)
  
  patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 6) +
    patchwork::plot_annotation(title = title) &
    ggplot2::theme(text = ggplot2::element_text(family = base_family))
}

# Plot temporal aggregates
plot_temporal_trends <- function(
    data, time_col, y_col, time_unit = "year", base_family = "serif") {
  cols <- spatial_palette()
  df_time <- data |>
    dplyr::mutate(
      agg_block = lubridate::floor_date(.data[[time_col]], unit = time_unit),
      month = lubridate::month(.data[[time_col]], label = TRUE)
    )
  y_agg <- df_time |>
    dplyr::group_by(agg_block) |>
    dplyr::summarise(
      Total = sum(.data[[y_col]], na.rm = TRUE), .groups = "drop"
    )
  p_agg <- ggplot2::ggplot(y_agg, ggplot2::aes(x = agg_block, y = Total)) +
    ggplot2::geom_col(fill = cols["primary"]) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::comma(Total)),
      vjust = -0.5, size = 3.5, fontface = "bold", family = base_family
    ) +
    ggplot2::labs(
      x = NULL, y = "Total",
      title = sprintf("Totals by %s", stringr::str_to_title(time_unit))
    ) +
    theme_bamlss_spatial(base_family = base_family)
  y_mon <- df_time |>
    dplyr::group_by(month) |>
    dplyr::summarise(
      Total = sum(.data[[y_col]], na.rm = TRUE), .groups = "drop"
    )
  p_month <- ggplot2::ggplot(y_mon, ggplot2::aes(x = month, y = Total)) +
    ggplot2::geom_col(fill = cols["secondary"]) +
    ggplot2::labs(
      x = NULL, y = "Total", title = "Seasonality by calendar month"
    ) +
    theme_bamlss_spatial(base_family = base_family)
  y_line <- df_time |>
    dplyr::group_by(!!rlang::sym(time_col)) |>
    dplyr::summarise(
      Total = sum(.data[[y_col]], na.rm = TRUE), .groups = "drop"
    )
  p_time <- ggplot2::ggplot(
    y_line, ggplot2::aes(x = .data[[time_col]], y = Total)
  ) +
    ggplot2::geom_col(width = 20, fill = cols["neutral"], alpha = 0.9) +
    ggplot2::geom_smooth(
      method = "loess", span = 0.2,
      color = cols["highlight"], fill = "gray70", alpha = 0.3
    ) +
    ggplot2::labs(
      x = NULL, y = "Total", title = "Monthly totals and LOESS trend"
    ) +
    theme_bamlss_spatial(base_family = base_family)
  ((p_agg + p_month) / p_time + patchwork::plot_layout(heights = c(1, 1.5))) &
    ggplot2::theme(text = ggplot2::element_text(family = base_family))
}

# Plot spatial outcome
plot_spatial_risk <- function(
    data, shape, space_col, y_col, shape_id, label_col = shape_id,
    aggregation = c("sum", "mean"), base_family = "serif") {
  
  aggregation <- match.arg(aggregation)
  s_agg <- data |>
    dplyr::group_by(!!rlang::sym(space_col)) |>
    dplyr::summarise(
      Total = if (aggregation == "mean") {
        z <- .data[[y_col]]
        if (all(!is.finite(z))) NA_real_ else mean(z[is.finite(z)])
      } else {
        sum(.data[[y_col]], na.rm = TRUE)
      },
      .groups = "drop"
    )
  
  map_data <-
    shape |> dplyr::left_join(s_agg, by = stats::setNames(space_col, shape_id))
  label_data <- sf_label_points(map_data)
  
  legend_title <- if (aggregation == "mean") "Mean" else "Total"
  finite_values <- map_data$Total[is.finite(map_data$Total)]
  signed <-
    length(finite_values) && any(finite_values < 0) && any(finite_values > 0)
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = map_data, ggplot2::aes(fill = Total), color = "gray80",
      linewidth = 0.1
    )
  
  p <- if (signed) {
    p + ggplot2::scale_fill_gradient2(
      low = "#2980b9", mid = "white", high = "#c0392b",
      midpoint = 0, name = legend_title, na.value = "transparent"
    )
  } else {
    p + ggplot2::scale_fill_gradient(
      low = "white", high = "#b10026", name = legend_title,
      na.value = "transparent"
    )
  }
  
  p +
    ggplot2::geom_sf_text(
      data = label_data, ggplot2::aes(label = !!rlang::sym(label_col)),
      size = 2.5, check_overlap = TRUE, family = base_family
    ) +
    ggplot2::theme_void() +
    ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
    ggspatial::annotation_north_arrow(
      location = "br", which_north = "true",
      style = ggspatial::north_arrow_fancy_orienteering
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      text = ggplot2::element_text(family = base_family)
    )
}

# Plot randomized quantile residual diagnostics
plot_rqr_diagnostics <- function(
    model, data, y_col, title = "Diagnostics", seed = 20260913L,
    base_family = "serif") {
  cols <- spatial_palette()
  res <- tryCatch(
    quantile_residuals(model = model, data = data, y_col = y_col, seed = seed),
    error = function(e) {
      warning(
        "Randomized quantile residuals failed: ", conditionMessage(e),
        call. = FALSE
      )
      numeric(0)
    }
  )
  res <- res[is.finite(res)]
  
  if (length(res) < 10L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text", x = 0, y = 0,
          label = "Too few finite randomized quantile residuals.",
          family = base_family
        ) +
        ggplot2::labs(title = title) +
        ggplot2::theme_void()
    )
  }
  
  n <- length(res)
  res_sorted <- sort(res)
  theo_q <- stats::qnorm(stats::ppoints(n))
  
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(
      ".Random.seed", envir = .GlobalEnv
    )
  }, add = TRUE)
  
  set.seed(as.integer(seed) + 1L)
  sim_mat <- apply(matrix(stats::rnorm(n * 500L), ncol = 500L), 2L, sort)
  
  df_qq <- data.frame(
    Observed = res_sorted,
    Theoretical = theo_q,
    Lower = apply(sim_mat, 1L, stats::quantile, probs = 0.025, names = FALSE),
    Upper = apply(sim_mat, 1L, stats::quantile, probs = 0.975, names = FALSE)
  )
  
  p1 <- ggplot2::ggplot(df_qq, ggplot2::aes(x = Theoretical, y = Observed)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = Lower, ymax = Upper), fill = "grey90", alpha = 0.5
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = Lower), color = "black", linewidth = 0.6
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = Upper), color = "black", linewidth = 0.6
    ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1, color = "black", linetype = "dashed",
      linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 1.5, alpha = 0.8, color = "black") +
    ggplot2::labs(
      title = paste0(title, ": QQ plot"),
      x = expression(paste("Quantiles of ", N(0, 1))),
      y = expression(hat(r)[st])
    ) +
    theme_bamlss_spatial(base_family = base_family)
  
  p2 <- ggplot2::ggplot(
    df_qq, ggplot2::aes(x = Theoretical, y = Observed - Theoretical)
  ) +
    ggplot2::geom_hline(
      yintercept = 0, color = cols["highlight"], linewidth = 1
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = Lower - Theoretical, ymax = Upper - Theoretical),
      fill = "grey85", alpha = 0.5
    ) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6) +
    ggplot2::geom_smooth(
      method = "gam", formula = y ~ s(x, bs = "cs"), color = cols["primary"],
      se = FALSE, linewidth = 0.5
    ) +
    ggplot2::labs(
      title = "Worm plot",
      x = expression(paste("Quantiles of ", N(0, 1))),
      y = expression(hat(r)[(i)] - z[(i)])
    ) +
    ggplot2::coord_cartesian(ylim = c(-2, 2)) +
    theme_bamlss_spatial(base_family = base_family)
  
  df_res <- data.frame(Index = seq_len(n), Residual = res)
  p3 <- ggplot2::ggplot(df_res, ggplot2::aes(x = Index, y = Residual)) +
    ggplot2::geom_hline(
      yintercept = c(-3, 3), color = cols["highlight"], linetype = "dotted",
      linewidth = 1
    ) +
    ggplot2::geom_point(alpha = 0.4, size = 0.5) +
    ggplot2::geom_point(
      data = df_res[abs(df_res$Residual) > 3, , drop = FALSE],
      color = cols["highlight"], size = 1.2
    ) +
    ggplot2::labs(
      title = "Residual index plot", x = "Observation index",
      y = expression(hat(r)[st])
    ) +
    theme_bamlss_spatial(base_family = base_family)
  
  (p1 + p2 + p3) &
    ggplot2::theme(text = ggplot2::element_text(family = base_family))
}

# Plot randomized PIT diagnostics
get_pit_plots <- function(
    model, data, y_col, label, seed = 20260913L, base_family = "serif") {
  cols <- spatial_palette()
  rqr <- tryCatch(
    quantile_residuals(model = model, data = data, y_col = y_col, seed = seed),
    error = function(e) {
      warning(
        "Randomized PIT calculation failed: ", conditionMessage(e),
        call. = FALSE
      )
      numeric(0)
    }
  )
  rqr <- rqr[is.finite(rqr)]
  
  if (length(rqr) < 10L) {
    dummy <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0.5, y = 0.5, label = "PIT unavailable", size = 3,
        family = base_family
      ) +
      ggplot2::labs(title = label) +
      ggplot2::theme_void()
    return(list(hist = dummy, ecdf = dummy))
  }
  
  pit <- stats::pnorm(rqr)
  p_hist <- ggplot2::ggplot(data.frame(pit = pit), ggplot2::aes(x = pit)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)), bins = 15,
      fill = cols["primary"], color = "white", alpha = 0.8
    ) +
    ggplot2::geom_hline(
      yintercept = 1, linetype = "dashed", color = cols["highlight"]
    ) +
    ggplot2::labs(title = label, x = expression(hat(u)[st]), y = "Density") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      text = ggplot2::element_text(family = base_family)
    )
  
  n <- length(pit)
  x <- sort(pit)
  y <- seq_len(n) / n
  k <- 1.36 / sqrt(n)
  df_ecdf <- data.frame(
    x = x, y = y, lower = pmax(0, x - k), upper = pmin(1, x + k)
  )
  
  p_ecdf <- ggplot2::ggplot(df_ecdf, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5
    ) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, linetype = "dashed", color = cols["highlight"]
    ) +
    ggplot2::geom_line(color = cols["primary"], linewidth = 0.8) +
    ggplot2::labs(
      title = label, x = expression(hat(u)[st]),
      y = expression(hat(F)[n](u))
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      text = ggplot2::element_text(family = base_family)
    ) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1))
  
  list(hist = p_hist, ecdf = p_ecdf)
}

# Plot parametric effects
plot_fixed_effects <- function(model, base_family = "serif") {
  cols <- spatial_palette()
  df_fixed <- extract_parametric_summaries(
    model, include_intercept = FALSE, strict = TRUE
  )
  if (!nrow(df_fixed)) {
    return(
      list(data = df_fixed, plot = ggplot2::ggplot() + ggplot2::theme_void())
    )
  }
  
  df_fixed$Term <- stringr::str_remove(df_fixed$Term, "\\.{3}\\d+")
  df_fixed$Term <- stringr::str_remove(df_fixed$Term, "^b\\.")
  df_fixed$Term <- stringr::str_remove(df_fixed$Term, "TRUE$")
  df_fixed$Display <- format_term_label(df_fixed$Term)
  
  plot_df <- df_fixed
  plot_df$Parameter <- ifelse(
    plot_df$Parameter == "s", "zeta", plot_df$Parameter
  )
  g <- ggplot2::ggplot(
    plot_df, ggplot2::aes(x = Estimate, y = stats::reorder(Display, Estimate))
  ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = Lower_95, xmax = Upper_95), orientation = "y",
      width = 0.2, color = cols["primary"], linewidth = 0.9
    ) +
    ggplot2::geom_point(size = 2.3, color = cols["highlight"]) +
    ggplot2::facet_wrap(
      ~Parameter, scales = "free_x", labeller = ggplot2::label_parsed
    ) +
    ggplot2::labs(
      title = "Conditional post-selection parametric effects",
      x = expression(paste("Posterior coefficient ", beta)), y = NULL
    ) +
    theme_bamlss_spatial(base_family = base_family)
  list(data = df_fixed, plot = g)
}

# Extract smooth effect
extract_smooth_effect <- function(
    model, dataset, term_full, var_name, parameter = "mu") {
  nd_source <- prepare_bamlss_newdata(
    model = model, data = dataset, parameter = parameter,
    require_complete = FALSE
  )
  complete_rows <- attr(nd_source, "complete_rows")
  
  if (!any(complete_rows)) {
    stop(
      "No complete reference row is available for parameter `", parameter, "`.",
      call. = FALSE
    )
  }
  if (!var_name %in% names(nd_source)) {
    stop(
      "Smooth covariate `", var_name,
      "` is not present among the predictors for parameter `", parameter, "`.",
      call. = FALSE
    )
  }
  
  x <- dataset[[var_name]]
  if (!is.numeric(x)) {
    stop("Smooth covariate `", var_name, "` must be numeric.", call. = FALSE)
  }
  
  x <- x[is.finite(x)]
  if (length(unique(x)) < 2L) {
    stop("Insufficient finite variation in `", var_name, "`.", call. = FALSE)
  }
  
  grid_vals <- seq(min(x), max(x), length.out = 100L)
  reference_row <- which(complete_rows)[1L]
  nd <- nd_source[rep(reference_row, length(grid_vals)), , drop = FALSE]
  nd[[var_name]] <- grid_vals
  
  summarise_draws <- function(z) {
    c(
      mean = mean(z),
      lower = stats::quantile(z, 0.025, names = FALSE),
      upper = stats::quantile(z, 0.975, names = FALSE)
    )
  }
  preds <- stats::predict(
    model, newdata = nd, model = parameter, term = term_full,
    match.names = TRUE, FUN = summarise_draws, intercept = FALSE, drop = FALSE
  )
  pred_mat <- as.data.frame(preds)
  mean_col <- grep("mean", names(pred_mat), ignore.case = TRUE)
  lower_col <- grep("lower", names(pred_mat), ignore.case = TRUE)
  upper_col <- grep("upper", names(pred_mat), ignore.case = TRUE)
  
  if (length(mean_col) != 1L ||
      length(lower_col) != 1L ||
      length(upper_col) != 1L) {
    stop(
      "Unexpected prediction structure for smooth term `", term_full,
      "` in parameter `", parameter, "`.", call. = FALSE
    )
  }
  
  data.frame(
    Value = grid_vals, eff = pred_mat[[mean_col]],
    lower = pred_mat[[lower_col]], upper = pred_mat[[upper_col]],
    Term = term_full, Parameter = parameter, stringsAsFactors = FALSE)
}

# Plot smooth effects
plot_smooth_effects <- function(
    model, dataset, smooth_terms_dict, parameter = "mu", base_family = "serif"
) {
  errors <- character(0)
  pieces <- lapply(names(smooth_terms_dict), function(term) {
    tryCatch(
      extract_smooth_effect(
        model, dataset, term, smooth_terms_dict[[term]], parameter
      ),
      error = function(e) {
        errors <<- c(errors, paste0(term, ": ", conditionMessage(e)))
        NULL
      }
    )
  })
  df_all <- dplyr::bind_rows(pieces)
  if (!nrow(df_all)) {
    stop(
      "No requested smooth effect could be extracted. ",
      paste(errors, collapse = " | "), call. = FALSE
    )
  }
  
  df_all$Term <- format_term_label(df_all$Term)
  ggplot2::ggplot(df_all, ggplot2::aes(x = Value, y = eff)) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper), fill = "#2980b9", alpha = 0.3
    ) +
    ggplot2::geom_line(color = "#2980b9", linewidth = 1) +
    ggplot2::facet_wrap(~Term, scales = "free", ncol = 3) +
    ggplot2::labs(
      x = "Covariate value",
      y = bquote(
        "Partial effect" ~ hat(f)(x) ~ "on" ~ hat(eta)[.(as.name(parameter))]
      )
    ) +
    theme_bamlss_spatial(base_family = base_family)
}

# Plot spatio-temporal model-implied fatalities
plot_expected_spatiotemporal_risk <- function(
    model, data, shape, space_col, time_col, shape_id, label_col = shape_id,
    time_unit = "year", functional = c("mean", "location"),
    base_family = "serif", scale_max = NULL) {
  
  functional <- match.arg(functional)
  fam <- model$family$family
  data <- as.data.frame(data)
  
  data$posterior_functional <- if (identical(functional, "mean")) {
    compute_expected_value(model)
  } else {
    compute_location_intensity(model)
  }
  
  if (length(data$posterior_functional) != nrow(data)) {
    stop(
      "Model predictions are not aligned with the plotting data.", call. = FALSE
    )
  }
  
  df_agg <- data |>
    dplyr::mutate(
      time_agg = lubridate::floor_date(.data[[time_col]], unit = time_unit)
    ) |>
    dplyr::group_by(!!rlang::sym(space_col), time_agg) |>
    dplyr::summarise(
      total = sum(posterior_functional, na.rm = TRUE), .groups = "drop"
    ) |>
    dplyr::mutate(time_agg = format(
      time_agg, if (identical(time_unit, "year")) "%Y" else "%Y-%m"
    ))
  
  map_data <- shape |>
    dplyr::left_join(df_agg, by = stats::setNames(space_col, shape_id))
  
  label_data <- sf_label_points(map_data)
  
  title_label <- if (identical(functional, "mean")) {
    sprintf("%s expected fatalities", stringr::str_to_title(time_unit))
  } else {
    sprintf("%s location intensity", stringr::str_to_title(time_unit))
  }
  legend_label <- "Location\nintensity"
  if (identical(functional, "mean")) legend_label <- "Expected\nfatalities"
  
  # An externally supplied upper limit gives the same colour meaning to
  # expected-fatality maps from distinct families. Do not pool location
  # intensities (RSB) with expected fatality counts.
  if (!is.null(scale_max)) {
    scale_max <- suppressWarnings(as.numeric(scale_max))
    if (length(scale_max) != 1L || !is.finite(scale_max) || scale_max <= 0) {
      stop("`scale_max` must be one positive finite number.", call. = FALSE)
    }
    if (any(df_agg$total > scale_max, na.rm = TRUE)) {
      stop("`scale_max` is below the largest mapped value.", call. = FALSE)
    }
  }
  fill_limits <- if (is.null(scale_max)) c(0, NA_real_) else c(0, scale_max)

  # Display-only encoding: white indicates fewer than one expected fatality
  # in a district-year. Do not alter the underlying fitted expectations.
  # Location-intensity maps (RSB) retain their unthresholded predictions.
  white_below <- if (identical(functional, "mean")) 1 else 0
  map_data$display_total <- map_data$total
  map_data$display_total[
    !is.na(map_data$display_total) & map_data$display_total < white_below
  ] <- 0

  ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = map_data, ggplot2::aes(fill = display_total), color = "gray80",
      linewidth = 0.1
    ) +
    ggplot2::scale_fill_gradientn(
      # Pure white for zero/near-zero, then the sequential YlOrRd palette.
      colors = c("#FFFFFF", "#ffffcc", "#fed976", "#fd8d3c", "#f03b20",
        "#800026"),
      trans = scales::pseudo_log_trans(sigma = 1),
      limits = fill_limits,
      breaks = function(lim) {
        b <- c(0, 1, 10, 50, 250, 1000, 5000, 10000, 25000,
               50000, 100000)
        b[b >= lim[1L] & b <= lim[2L]]
      },
      labels = scales::label_comma(),
      name = legend_label,
      na.value = "#E5E7EB",
      guide = ggplot2::guide_colorbar(
        barheight = ggplot2::unit(2.5, "in"),
        barwidth = ggplot2::unit(0.15, "in"),
        ticks.linewidth = 1,
        frame.colour = "black"
      )
    ) +
    ggplot2::geom_sf_text(
      data = label_data, ggplot2::aes(label = !!rlang::sym(label_col)),
      size = 2.0, color = "black", check_overlap = TRUE, family = base_family
    ) +
    ggplot2::facet_wrap(~time_agg, ncol = 4) +
    ggplot2::labs(
      title = title_label,
      subtitle = paste("Model:", fam),
      caption = if (identical(functional, "mean"))
        "White = fewer than 1 expected fatality per district-year; grey = missing data"
      else NULL
    ) +
    ggplot2::theme_void() +
    ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
    ggspatial::annotation_north_arrow(
      location = "br", which_north = "true",
      style = ggspatial::north_arrow_fancy_orienteering
    ) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      text = ggplot2::element_text(family = base_family),
      legend.position = "right"
    )
}

# Readable label for a BAMLSS sample column
trace_label <- function(x) {
  par <- sub("^([A-Za-z]+)\\.[ps]\\..*$", "\\1", x)
  par[par == x] <- ""
  term <- sub("^[A-Za-z]+\\.[ps]\\.", "", x)
  term <- sub("^s\\((.*)\\)\\.(.*)$", "\\1 [\\2]", term)
  term <- sub("^\\(Intercept\\)$", "Intercept", term)
  term <- sub("^year_factor", "Year ", term)
  term <- sub("^lin\\(~?(.*)\\)$", "\\1", term)
  term <- vapply(term, function(z) {
    m <- regmatches(z, regexec("^(.*) \\[(.*)\\]$", z))[[1L]]
    if (length(m) == 3L) paste0(format_term_label(m[2L]), " [", m[3L], "]")
    else format_term_label(z)
  }, character(1), USE.NAMES = FALSE)
  trimws(paste0(ifelse(nzchar(par), paste0(par, ": "), ""), term))
}

# Trace plots of all retained parameters, one page per `per_page` panels
plot_trace_pages <- function(
    model, per_page = 12L, ncol = 3L, base_family = "serif") {
  chains <- extract_sample_chains(model)
  nm <- colnames(chains[[1L]])
  skip <- grepl("\\.(accepted|alpha|edf)$", nm) |
    nm %in% c("logLik", "logPost", "DIC", "pd", "alpha")
  keep <- nm[!skip]
  cols <- grDevices::hcl.colors(length(chains), "Dark 3")
  one <- function(v) {
    ess <- sum(vapply(chains, function(ch) {
      effective_sample_size(ch[, v])
    }, numeric(1)))
    rh <- split_rhat(chains, v)
    key <- sprintf("%s\nR-hat %s, ESS %s", trace_label(v),
                   ifelse(is.finite(rh), sprintf("%.2f", rh), "NA"),
                   ifelse(is.finite(ess), sprintf("%.0f", ess), "NA"))
    data.frame(
      panel = key, rhat = rh,
      draw = unlist(lapply(chains, function(ch) seq_len(nrow(ch)))),
      chain = factor(rep(seq_along(chains), vapply(chains, nrow, 1L))),
      value = unlist(lapply(chains, function(ch) ch[, v]))
    )
  }
  info <- lapply(keep, one)
  worst <- vapply(info, function(d) {
    max(d$rhat[1L], 0, na.rm = TRUE)
  }, numeric(1))
  info <- info[order(-worst)]
  pages <- split(seq_along(info), ceiling(seq_along(info) / per_page))
  lapply(seq_along(pages), function(k) {
    d <- do.call(rbind, info[pages[[k]]])
    d$panel <- factor(d$panel, levels = unique(d$panel))
    ggplot2::ggplot(d, ggplot2::aes(draw, value, colour = chain)) +
      ggplot2::geom_line(linewidth = 0.25, alpha = 0.85) +
      ggplot2::facet_wrap(~panel, ncol = ncol, scales = "free_y") +
      ggplot2::scale_colour_manual(values = cols, name = "Chain") +
      ggplot2::labs(
        title = paste("Trace plots:", model$family$family),
        subtitle = bquote(
          "Parameters ordered by decreasing" ~ hat(R) ~
            "(page" ~ .(k) ~ "of" ~ .(length(pages)) * ")"
        ),
        x = "Saved draw", y = "Posterior value"
      ) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        nrow = 1L, override.aes = list(linewidth = 1.2)
      )) +
      ggplot2::theme_bw(base_family = base_family) +
      ggplot2::theme(
        legend.position = "bottom",
        strip.text = ggplot2::element_text(size = 7.5),
        panel.grid.minor = ggplot2::element_blank()
      )
  })
}

# Write the trace plot pages of one model to a multi-page PDF
save_trace_pdf <- function(model, file, width = 11, height = 8.5) {
  pages <- plot_trace_pages(model)
  grDevices::cairo_pdf(file, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (p in pages) print(p)
  invisible(file)
}

# Format term labels for output
format_term_label <- function(x) {
  x <- as.character(x)
  x <- stringr::str_remove(x, "^s\\(")
  x <- stringr::str_remove(x, "\\)$")
  x <- stringr::str_remove_all(x, "(_z|_s|_log1p_z|_log1p|_km_km2|_km|TRUE)$")
  
  x <- stringr::str_replace_all(x, "_lag1", " (Lag 1)")
  x <- stringr::str_replace_all(x, "_within", " (Within)")
  x <- stringr::str_replace_all(x, "_between", " (Between)")
  x <- stringr::str_replace_all(x, "^dist_", "Distance to ")
  x <- stringr::str_replace_all(x, "^distance_", "Distance to ")
  x <- stringr::str_replace_all(x, "bord", "Border ")
  
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_to_title(x)
  
  x <- stringr::str_replace_all(x, "\\bLng\\b", "LNG")
  x <- stringr::str_replace_all(x, "\\bNtl\\b", "NTL")
  x <- stringr::str_replace_all(x, "\\bNdvi\\b", "NDVI")
  x <- stringr::str_replace_all(x, "\\bLoo\\b", "(LOO)")
  
  return(trimws(x))
}

# Extract prediction variables
bamlss_prediction_variables <- function(model, parameter = NULL) {
  fobj <- model$formula
  if (is.null(fobj)) {
    stop("The fitted BAMLSS object does not contain a formula.", call. = FALSE)
  }
  
  fobj <- as.list(fobj)
  if (!is.null(parameter)) {
    if (is.null(names(fobj)) || !parameter %in% names(fobj)) {
      stop(
        "Distributional parameter `", parameter,
        "` is not present in the fitted BAMLSS formula.", call. = FALSE
      )
    }
    fobj <- fobj[parameter]
  }
  
  formulae <- list()
  collect_formulae <- function(x) {
    if (inherits(x, "formula")) {
      formulae[[length(formulae) + 1L]] <<- x
      return(invisible(NULL))
    }
    if (!is.list(x)) return(invisible(NULL))
    if (!is.null(x$fake.formula) && inherits(x$fake.formula, "formula")) {
      formulae[[length(formulae) + 1L]] <<- x$fake.formula
      return(invisible(NULL))
    }
    if (!is.null(x$formula) && inherits(x$formula, "formula")) {
      formulae[[length(formulae) + 1L]] <<- x$formula
      return(invisible(NULL))
    }
    for (z in x) collect_formulae(z)
    invisible(NULL)
  }
  
  for (z in fobj) collect_formulae(z)
  if (!length(formulae)) {
    stop(
      "No usable formula could be recovered from the fitted BAMLSS object.",
      call. = FALSE
    )
  }
  
  all_variables <- unique(unlist(lapply(formulae, all.vars), use.names = FALSE))
  response_variables <- unique(unlist(lapply(formulae, function(f) {
    tt <- stats::terms(f)
    if (isTRUE(attr(tt, "response") == 1L) && length(f) >= 3L) {
      all.vars(f[[2L]])
    } else {
      character(0)
    }
  }), use.names = FALSE))
  
  setdiff(all_variables, response_variables)
}

# Prepare BAMLSS newdata
prepare_bamlss_newdata <- function(
    model, data, parameter = NULL, require_complete = TRUE) {
  data <- as.data.frame(data)
  if (!nrow(data)) stop("Prediction data contain no rows.", call. = FALSE)
  
  vars <- bamlss_prediction_variables(model, parameter = parameter)
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars)) {
    stop(
      "Prediction data are missing model variable(s): ",
      paste(missing_vars, collapse = ", "), ".", call. = FALSE
    )
  }
  
  if (length(vars)) {
    nd <- data[, vars, drop = FALSE]
  } else {
    nd <- data.frame(.bamlss_prediction_row = seq_len(nrow(data)))
  }
  
  complete <- stats::complete.cases(nd)
  numeric_vars <- names(nd)[vapply(nd, is.numeric, logical(1))]
  for (v in numeric_vars) complete <- complete & is.finite(nd[[v]])
  
  if (isTRUE(require_complete) && any(!complete)) {
    bad <- vapply(
      nd,
      function(x) if (is.numeric(x)) {
        sum(is.na(x) | !is.finite(x))
      } else {
        sum(is.na(x))
      },
      integer(1)
    )
    bad <- bad[bad > 0L]
    stop(
      "Prediction data contain missing/non-finite model covariates: ",
      paste(sprintf("%s=%d", names(bad), bad), collapse = ", "),
      ". Prediction rows will not be silently dropped.", call. = FALSE
    )
  }
  
  attr(nd, "complete_rows") <- complete
  attr(nd, "prediction_variables") <- vars
  nd
}

# Predict BAMLSS distributional parameters
predict_bamlss_parameters <- function(
    model,
    newdata = NULL,
    drop = FALSE) {
  
  if (is.null(newdata)) {
    return(stats::predict(model, type = "parameter", drop = drop))
  }
  
  nd <- prepare_bamlss_newdata(
    model = model,
    data = newdata,
    require_complete = TRUE
  )
  
  stats::predict(model, newdata = nd, type = "parameter", drop = drop)
}

# PIG cumulative distribution function
pig_cdf <- function(q, mu, sigma) {
  n <- max(length(q), length(mu), length(sigma))
  q <- rep_len(as.numeric(q), n)
  mu <- rep_len(as.numeric(mu), n)
  sigma <- rep_len(as.numeric(sigma), n)
  
  if (any(!is.finite(mu)) || any(mu <= 0)) {
    stop("PIG `mu` must contain finite positive values.", call. = FALSE)
  }
  if (any(!is.finite(sigma)) || any(sigma <= 0)) {
    stop("PIG `sigma` must contain finite positive values.", call. = FALSE)
  }
  
  out <- numeric(n)
  below_support <- q < 0
  out[below_support] <- 0
  
  use_poisson <- (!below_support & sigma <= 1e-4)
  use_pig <- (!below_support & sigma > 1e-4)
  
  if (any(use_poisson)) {
    out[use_poisson] <- stats::ppois(
      q = q[use_poisson], lambda = mu[use_poisson], lower.tail = TRUE,
      log.p = FALSE
    )
  }
  
  if (any(use_pig)) {
    out[use_pig] <- gamlss.dist::pPIG(
      q = q[use_pig], mu = mu[use_pig], sigma = sigma[use_pig],
      lower.tail = TRUE, log.p = FALSE
    )
  }
  
  pmin(pmax(out, 0), 1)
}

# Distribution cumulative distribution function
family_cdf <- function(family, q, parameters) {
  fam_name <- as.character(family$family)[1L]
  
  if (identical(fam_name, "PIG")) {
    return(pig_cdf(q = q, mu = parameters$mu, sigma = parameters$sigma))
  }
  
  if (identical(fam_name, "ZIPIG")) {
    n <- max(
      length(q), length(parameters$mu), length(parameters$sigma),
      length(parameters$nu)
    )
    qq <- rep_len(as.numeric(q), n)
    nu <- rep_len(as.numeric(parameters$nu), n)
    
    if (any(!is.finite(nu)) || any(nu < 0) || any(nu > 1)) {
      stop("ZIPIG `nu` must lie in [0, 1].", call. = FALSE)
    }
    
    base_cdf <- pig_cdf(q = qq, mu = parameters$mu, sigma = parameters$sigma)
    out <- nu + (1 - nu) * base_cdf
    out[qq < 0] <- 0
    return(pmin(pmax(out, 0), 1))
  }
  
  if (is.null(family$p) || !is.function(family$p)) {
    stop("Family `", fam_name, "` does not provide a CDF.", call. = FALSE)
  }
  
  out <- as.numeric(family$p(q, parameters))
  if (length(out) != length(q)) {
    stop(
      "Family CDF returned ", length(out), " values for ", length(q),
      " requested quantiles.", call. = FALSE
    )
  }
  
  out
}

# Extract point on surface for spatial labels
sf_label_points <- function(x) {
  if (!inherits(x, "sf")) stop("`x` must be an sf object.", call. = FALSE)
  
  original_crs <- sf::st_crs(x)
  geometry_type <- unique(
    as.character(sf::st_geometry_type(x, by_geometry = TRUE))
  )
  polygonal <- geometry_type %in% c("POLYGON", "MULTIPOLYGON")
  
  if (!any(polygonal)) return(x)
  
  if (sf::st_is_longlat(x)) {
    ll <- sf::st_transform(x, 4326)
    bb <- sf::st_bbox(ll)
    
    lon0 <- as.numeric((bb["xmin"] + bb["xmax"]) / 2)
    lat0 <- as.numeric((bb["ymin"] + bb["ymax"]) / 2)
    
    zone <- floor((lon0 + 180) / 6) + 1L
    epsg <- if (lat0 >= 0) 32600L + zone else 32700L + zone
    
    projected <- sf::st_transform(x, epsg)
    sf::st_geometry(projected) <-
      sf::st_point_on_surface(sf::st_geometry(projected))
    
    return(sf::st_transform(projected, original_crs))
  }
  
  out <- x
  sf::st_geometry(out) <- sf::st_point_on_surface(sf::st_geometry(out))
  
  out
}

# Compute randomized quantile residuals reproducibly
quantile_residuals <- function(model, data, y_col, seed = 20260913L) {
  if (missing(data) || is.null(data)) {
    stop(
      "`data` must be supplied for response-alignment validation.",
      call. = FALSE
    )
  }
  if (missing(y_col) ||
      is.null(y_col) ||
      length(y_col) != 1L ||
      !nzchar(y_col)) {
    stop("`y_col` must identify exactly one response column.", call. = FALSE)
  }
  
  data <- as.data.frame(data)
  if (!y_col %in% names(data)) {
    stop(
      "Response variable `", y_col, "` is not present in the evaluation data.",
      call. = FALSE
    )
  }
  
  y <- suppressWarnings(as.numeric(data[[y_col]]))
  if (length(y) != nrow(data) ||
      any(!is.finite(y)) ||
      any(y < 0) ||
      any(abs(y - round(y)) > 1e-8)) {
    stop(
      "Quantile residuals need finite non-negative integer counts.",
      call. = FALSE
    )
  }
  
  y_model <- model$y
  if (is.data.frame(y_model)) y_model <- y_model[[1L]]
  if (is.matrix(y_model)) y_model <- y_model[, 1L]
  
  y_model <- suppressWarnings(as.numeric(y_model))
  
  if (!length(y_model)) {
    stop(
      "The fitted BAMLSS object does not retain its response vector.",
      call. = FALSE
    )
  }
  if (length(y_model) != length(y) ||
      !isTRUE(all.equal(y, y_model, tolerance = 0, check.attributes = FALSE))) {
    stop(
      "The response in `data` is not aligned with the fitted model response.",
      call. = FALSE
    )
  }
  
  fam <- model$family
  pars <- predict_bamlss_parameters(model = model, newdata = NULL, drop = FALSE)
  
  missing_parameters <- setdiff(model$family$names, names(pars))
  if (length(missing_parameters)) {
    stop(
      "BAMLSS prediction is missing distributional parameter(s).", call. = FALSE
    )
  }
  
  prediction_lengths <- vapply(pars, length, integer(1))
  if (any(prediction_lengths != length(y))) {
    stop(
      "BAMLSS parameter predictions are not aligned with the fitted sample.",
      call. = FALSE
    )
  }
  
  lower <- tryCatch(
    family_cdf(family = fam, q = y - 1, parameters = pars),
    error = function(e) stop(
      "Failed to evaluate lower predictive CDF for family `", fam$family, "`: ",
      conditionMessage(e), call. = FALSE
    )
  )
  
  upper <- tryCatch(
    family_cdf(family = fam, q = y, parameters = pars),
    error = function(e) stop(
      "Failed to evaluate upper predictive CDF for family `", fam$family, "`: ",
      conditionMessage(e), call. = FALSE
    )
  )
  
  n <- length(y)
  lower[y == 0] <- 0
  
  if (any(!is.finite(lower)) || any(!is.finite(upper))) {
    stop("Predictive CDF produced non-finite probabilities.", call. = FALSE)
  }
  
  tolerance <- 1e-10
  lower <- pmin(pmax(lower, 0), 1)
  upper <- pmin(pmax(upper, 0), 1)
  
  bad_order <- upper + tolerance < lower
  if (any(bad_order)) {
    stop(
      "Predictive CDF is non-monotone for ", sum(bad_order), " observation(s).",
      call. = FALSE
    )
  }
  upper <- pmax(upper, lower)
  
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(
      ".Random.seed", envir = .GlobalEnv
    )
  }, add = TRUE)
  
  set.seed(as.integer(seed))
  u <- stats::runif(n = n, min = lower, max = upper)
  
  probability_floor <- 1e-12
  u <- pmin(pmax(u, probability_floor), 1 - probability_floor)
  residuals <- stats::qnorm(u)
  
  if (length(residuals) != nrow(data) || any(!is.finite(residuals))) {
    stop(
      "Quantile residuals are not finite or not aligned with the data.",
      call. = FALSE
    )
  }
  
  residuals
}
# Rank-based convergence diagnostics (Vehtari et al., 2021) ------------------

# Split every chain (columns) in two halves; draws are in rows
split_chains_matrix <- function(x) {
  n <- nrow(x) - nrow(x) %% 2L
  h <- n / 2L
  cbind(x[seq_len(h), , drop = FALSE], x[h + seq_len(h), , drop = FALSE])
}

rank_normalise <- function(x) {
  r <- rank(x, ties.method = "average")
  matrix(stats::qnorm((r - 3 / 8) / (length(x) + 1 / 4)), nrow = nrow(x))
}

rhat_core <- function(m) {
  n <- nrow(m)
  w <- mean(apply(m, 2L, stats::var))
  b <- n * stats::var(colMeans(m))
  if (!is.finite(w) || !is.finite(b) || w <= 0) return(NA_real_)
  sqrt(((n - 1) / n * w + b / n) / w)
}

# Effective sample size from multi-chain autocovariances (Geyer pairs)
ess_core <- function(m) {
  n <- nrow(m)
  k <- ncol(m)
  w <- mean(apply(m, 2L, stats::var))
  var_plus <- (n - 1) / n * w + stats::var(colMeans(m))
  if (!is.finite(var_plus) || var_plus <= 0) return(NA_real_)
  acov <- rowMeans(vapply(seq_len(k), function(j) {
    a <- stats::acf(
      m[, j], type = "covariance", lag.max = n - 1L, plot = FALSE
    )
    a$acf[, 1L, 1L] * n / (n - 1)
  }, numeric(n)))
  rho <- 1 - (w - acov) / var_plus
  rho[1L] <- 1
  pairs <- rho[2L * seq_len(floor(n / 2)) - 1L] +
    rho[2L * seq_len(floor(n / 2))]
  neg <- which(pairs < 0)
  if (length(neg)) pairs <- pairs[seq_len(neg[1L] - 1L)]
  tau <- -1 + 2 * sum(cummin(pairs))
  n * k / max(tau, 1 / log10(n * k))
}

# R-hat (rank-normalised, folded), bulk and tail ESS; x is draws by chains
rank_diagnostics <- function(x) {
  x <- as.matrix(x)
  out <- c(Rhat = NA_real_, ESS_bulk = NA_real_, ESS_tail = NA_real_)
  if (nrow(x) < 8L || ncol(x) < 2L || any(!is.finite(x))) return(out)
  s <- split_chains_matrix(x)
  r <- c(
    rhat_core(rank_normalise(s)),
    rhat_core(rank_normalise(abs(s - stats::median(s))))
  )
  q <- stats::quantile(s, c(0.05, 0.95), names = FALSE)
  tails <- c(ess_core(1 * (s <= q[1L])), ess_core(1 * (s <= q[2L])))
  if (any(is.finite(r))) out["Rhat"] <- max(r, na.rm = TRUE)
  out["ESS_bulk"] <- ess_core(rank_normalise(s))
  if (any(is.finite(tails))) out["ESS_tail"] <- min(tails, na.rm = TRUE)
  out
}

# Draw-wise functionals of the posterior, chain by chain ---------------------

# Draws by observations as an n x S matrix, whatever the orientation returned
orient_draws <- function(a, n) {
  a <- as.matrix(a)
  if (nrow(a) == n) return(a)
  if (ncol(a) == n) return(t(a))
  stop(
    "Predicted draws do not match the number of observations.", call. = FALSE
  )
}

posterior_functional_draws <- function(model, n, k_vals, nsamps = 500L) {
  fam <- model$family
  smp <- model$samples
  chain_list <- if (inherits(smp, "mcmc.list")) {
    lapply(seq_along(smp), function(i) smp[[i]])
  } else {
    list(smp)
  }
  lapply(chain_list, function(ch) {
    m <- model
    m$samples <- ch
    pr <- stats::predict(
      m, type = "parameter", drop = FALSE, nsamps = nsamps,
      FUN = function(x) x
    )
    pr <- lapply(pr, orient_draws, n = n)
    idx <- unique(round(seq(1, ncol(pr[[1L]]), length.out = nsamps)))
    pr <- lapply(pr, function(a) a[, idx, drop = FALSE])
    loc <- pr[[if ("mu" %in% names(pr)) "mu" else "lambda"]]
    scal <- t(vapply(seq_along(idx), function(s) {
      par_s <- lapply(pr, function(a) a[, s])
      c(
        mean_location = mean(loc[, s]),
        vapply(k_vals, function(k) {
          1 - mean(family_cdf(fam, rep(k, n), par_s))
        }, numeric(1))
      )
    }, numeric(length(k_vals) + 1L)))
    colnames(scal) <- c("mean_location", paste0("P(Y>", k_vals, ")"))
    list(loc = loc, scal = scal)
  })
}

# Thresholds for exceedance events, taken from the observed counts
default_exceedance_thresholds <- function(y_obs) {
  pos <- y_obs[y_obs > 0L]
  unique(c(
    0L, if (length(pos)) as.integer(stats::quantile(pos, c(.5, .9), type = 1))
  ))
}

# Functional-level R-hat/ESS and ESS growth with the number of draws
compute_functional_diagnostics <- function(
    model, data, y_col, nsamps = 500L, rhat_cut = 1.01) {
  y_obs <- as.integer(data[[y_col]])
  n <- length(y_obs)
  k_vals <- default_exceedance_thresholds(y_obs)
  dr <- posterior_functional_draws(model, n, k_vals, nsamps)
  S <- min(vapply(dr, function(z) nrow(z$scal), integer(1)))
  cols <- colnames(dr[[1L]]$scal)
  stack_f <- function(f) {
    vapply(dr, function(z) z$scal[seq_len(S), f], numeric(S))
  }
  
  scalar <- do.call(rbind, lapply(cols, function(f) {
    data.frame(
      Functional = f, Draws_per_chain = S, Chains = length(dr),
      t(rank_diagnostics(stack_f(f))), check.names = FALSE
    )
  }))
  curve <- do.call(rbind, lapply(cols, function(f) {
    x <- stack_f(f)
    do.call(rbind, lapply(c(0.25, 0.5, 0.75, 1), function(a) {
      m <- max(16L, floor(a * S))
      data.frame(
        Functional = f, Draws_per_chain = m,
        t(rank_diagnostics(x[seq_len(m), , drop = FALSE])),
        check.names = FALSE
      )
    }))
  }))
  
  per_obs <- t(vapply(seq_len(n), function(i) {
    rank_diagnostics(
      vapply(dr, function(z) z$loc[i, seq_len(S)], numeric(S))
    )
  }, numeric(3)))
  colnames(per_obs) <- c("Rhat", "ESS_bulk", "ESS_tail")
  fin <- function(v) v[is.finite(v)]
  location <- data.frame(
    N = n, N_Undefined = sum(!is.finite(per_obs[, "Rhat"])),
    Max_Rhat = if (length(fin(per_obs[, "Rhat"]))) {
      max(fin(per_obs[, "Rhat"]))
    } else {
      NA_real_
    },
    Share_Rhat_Above_Cut = mean(fin(per_obs[, "Rhat"]) > rhat_cut),
    Min_ESS_bulk = min(fin(per_obs[, "ESS_bulk"]), Inf),
    Min_ESS_tail = min(fin(per_obs[, "ESS_tail"]), Inf)
  )
  list(scalar = scalar, curve = curve, location = location)
}

# Result tables and PIT figures --------------------------------

# Simultaneous ECDF envelope for n uniform draws (Saeilynoja et al., 2022)
ecdf_envelope <- function(
    n, level = 0.95, n_grid = min(n, 300L), n_sim = 4000L,
    seed = 20260913L) {
  z <- seq_len(n_grid) / (n_grid + 1)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit(
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv),
    add = TRUE
  )
  set.seed(seed)
  p_min <- vapply(seq_len(n_sim), function(b) {
    x <- findInterval(z, sort(stats::runif(n)))
    min(2 * pmin(
      stats::pbinom(x, n, z),
      stats::pbinom(x - 1L, n, z, lower.tail = FALSE)
    ))
  }, numeric(1))
  gam <- min(1, stats::quantile(p_min, 1 - level, names = FALSE))
  list(
    z = z, lower = stats::qbinom(gam / 2, n, z) / n,
    upper = stats::qbinom(1 - gam / 2, n, z) / n, gamma = gam
  )
}

# Four panels: histograms (20, 50 bins), ECDF and ECDF difference with bands
plot_pit_panels <- function(
    u, label = NULL, level = 0.95, bins = c(20L, 50L),
    base_family = "serif") {
  u <- u[is.finite(u)]
  n <- length(u)
  fill_col <- unname(spatial_palette()["primary"])
  base_theme <- ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      text = ggplot2::element_text(family = base_family),
      plot.title = ggplot2::element_text(size = 10, hjust = 0.5)
    )
  
  hist_panel <- function(b, tag) {
    br <- seq(0, 1, length.out = b + 1L)
    cnt <- tabulate(
      pmin(findInterval(u, br, rightmost.closed = TRUE), b), nbins = b
    )
    df <- data.frame(mid = (br[-1L] + br[-length(br)]) / 2, n = cnt)
    band <- stats::qbinom(c((1 - level) / 2, 1 - (1 - level) / 2), n, 1 / b)
    ggplot2::ggplot(df, ggplot2::aes(mid, n)) +
      ggplot2::annotate(
        "rect", xmin = 0, xmax = 1, ymin = band[1L], ymax = band[2L],
        fill = "grey85"
      ) +
      ggplot2::geom_col(
        width = 1 / b, fill = fill_col, colour = "grey30", linewidth = 0.15,
        alpha = 0.85
      ) +
      ggplot2::geom_hline(
        yintercept = n / b, colour = "grey35", linewidth = 0.3
      ) +
      ggplot2::labs(title = tag, x = expression(hat(u)[st]), y = "Count") +
      base_theme
  }
  
  env <- ecdf_envelope(n, level)
  df <- data.frame(
    z = env$z, Fn = stats::ecdf(u)(env$z), lo = env$lower, hi = env$upper
  )
  p_ecdf <- ggplot2::ggplot(df, ggplot2::aes(z)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lo, ymax = hi), fill = "grey85"
    ) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, colour = "grey35", linewidth = 0.3
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = Fn), colour = fill_col, linewidth = 0.5
    ) +
    ggplot2::labs(
      title = "(c) ECDF plot", x = expression(hat(u)[st]),
      y = expression(hat(F)[n](u))
    ) +
    base_theme
  p_diff <- ggplot2::ggplot(df, ggplot2::aes(z)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lo - z, ymax = hi - z), fill = "grey85"
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.3) +
    ggplot2::geom_line(
      ggplot2::aes(y = Fn - z), colour = fill_col, linewidth = 0.5
    ) +
    ggplot2::labs(
      title = "(d) ECDF difference plot", x = expression(hat(u)[st]),
      y = expression(hat(F)[n](u) - u)
    ) +
    base_theme
  
  panels <- list(
    hist_panel(bins[1L], sprintf("(a) Histogram, %d bins", bins[1L])),
    hist_panel(bins[2L], sprintf("(b) Histogram, %d bins", bins[2L])),
    p_ecdf, p_diff
  )
  out <- patchwork::wrap_plots(panels, nrow = 1L)
  if (!is.null(label)) out <- out + patchwork::plot_annotation(title = label)
  out
}

# ECDF-difference panels of all families in one figure
plot_ecdf_difference_grid <- function(
    u_list, level = 0.95, ncol = 4L, base_family = "serif") {
  df <- do.call(rbind, lapply(names(u_list), function(nm) {
    u <- u_list[[nm]][is.finite(u_list[[nm]])]
    env <- ecdf_envelope(length(u), level)
    data.frame(
      Family = nm, z = env$z, d = stats::ecdf(u)(env$z) - env$z,
      lo = env$lower - env$z, hi = env$upper - env$z
    )
  }))
  df$Family <- factor(df$Family, levels = names(u_list))
  ggplot2::ggplot(df, ggplot2::aes(z)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "grey85") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.3) +
    ggplot2::geom_line(
      ggplot2::aes(y = d), colour = spatial_palette()[["primary"]],
      linewidth = 0.5
    ) +
    ggplot2::facet_wrap(~ Family, ncol = ncol) +
    ggplot2::labs(
      x = expression(hat(u)[st]), y = expression(hat(F)[n](u) - u)
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      text = ggplot2::element_text(family = base_family),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA)
    )
}

# Formatting helpers -----------------------------------------------------------

tex_escape <- function(x) gsub("([&%$#_{}])", "\\\\\\1", as.character(x))

fmt_col <- function(x, digits = 2, best = c("none", "min", "max")) {
  best <- match.arg(best)
  s <- ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "--")
  if (best != "none" && any(is.finite(x))) {
    target <- if (best == "min") min(x[is.finite(x)]) else max(x[is.finite(x)])
    hit <- is.finite(x) & x == target
    s[hit] <- paste0("\\textbf{", s[hit], "}")
  }
  s
}

fmt_p <- function(p) {
  ifelse(
    is.finite(p),
    ifelse(p < 0.001, "$<$0.001", formatC(p, digits = 3, format = "f")), "--"
  )
}

col_or_na <- function(d, nm) {
  if (!is.null(d) && nm %in% names(d)) d[[nm]] else rep(NA_real_, nrow(d))
}

param_symbol <- function(x) {
  map <- c(
    mu = "$\\mu$", sigma = "$\\sigma$", nu = "$\\nu$", tau = "$\\tau$",
    lambda = "$\\lambda$", s = "$\\zeta$"
  )
  out <- unname(map[as.character(x)])
  ifelse(is.na(out), as.character(x), out)
}

param_group <- function(x) {
  map <- c(
    mu = "mu", lambda = "mu", sigma = "sigma", nu = "nu", s = "nu",
    tau = "tau"
  )
  unname(map[as.character(x)])
}

pretty_cov <- function(x) {
  lab <- format_term_label(x)
  lab <- gsub("\\bDag\\b", "DAG", lab)
  tex_escape(lab)
}

# Table environments; the caption goes above the table
latex_table <- function(
    df, caption, label, align = NULL, note = NULL, pre_header = NULL,
    long = FALSE, size = "footnotesize") {
  df[] <- lapply(df, as.character)
  if (is.null(align)) align <- paste0("l", strrep("r", ncol(df) - 1L))
  header <- paste0(paste(names(df), collapse = " & "), " \\\\")
  body <- paste0(apply(df, 1L, paste, collapse = " & "), " \\\\")
  note_line <- if (!is.null(note)) {
    paste0("\\par\\smallskip{\\scriptsize ", note, "}")
  }
  if (long) {
    return(c(
      paste0("{\\", size),
      paste0("\\begin{longtable}{", align, "}"),
      paste0("\\caption{", caption, "}\\label{", label, "}\\\\"),
      "\\toprule", pre_header, header, "\\midrule", "\\endfirsthead",
      "\\toprule", pre_header, header, "\\midrule", "\\endhead",
      "\\bottomrule", "\\endfoot", body, "\\end{longtable}",
      note_line, "}"
    ))
  }
  c(
    "\\begin{table}[!htbp]", "\\centering", paste0("\\", size),
    paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"),
    "\\begin{adjustbox}{max width=\\textwidth}",
    paste0("\\begin{tabular}{", align, "}"), "\\toprule", pre_header, header,
    "\\midrule", body, "\\bottomrule", "\\end{tabular}", "\\end{adjustbox}",
    note_line, "\\end{table}"
  )
}

read_csv_if <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

# Table builders ---------------------------------------------------------------

build_family_comparison <- function(fit, interval) {
  fit <- fit[order(fit$DIC), , drop = FALSE]
  fam <- fit[["Distribution Family"]]
  j <- if (is.null(interval)) rep(NA_integer_, nrow(fit)) else {
    match(fam, interval$Family)
  }
  pick <- function(nm) if (all(is.na(j))) rep(NA_real_, nrow(fit)) else {
    interval[[nm]][j]
  }
  waic <- col_or_na(fit, "WAIC")
  dic <- col_or_na(fit, "DIC")
  data.frame(
    Family = fam,
    WAIC = fmt_col(waic, 1), `$\\Delta$WAIC` = fmt_col(waic - min(waic), 1),
    DIC = fmt_col(dic, 1), `$\\Delta$DIC` = fmt_col(dic - min(dic), 1),
    LS = fmt_col(col_or_na(fit, "Logarithmic Score"), 3, "max"),
    QS = fmt_col(col_or_na(fit, "Quadratic Score"), 3, "max"),
    SPS = fmt_col(col_or_na(fit, "Spherical Score"), 3, "max"),
    `RPS$_K$` = fmt_col(col_or_na(fit, "RPS (K = max y)"), 3, "min"),
    CRPS = fmt_col(col_or_na(fit, "CRPS"), 3, "min"),
    `IS$_{90}$` = fmt_col(pick("IntervalScore"), 2, "min"),
    `Cov.$_{90}$ (\\%)` = fmt_col(100 * pick("Coverage"), 1),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

build_calibration_table <- function(fit) {
  fit <- fit[order(fit$DIC), , drop = FALSE]
  data.frame(
    Family = fit[["Distribution Family"]],
    KS = fmt_p(col_or_na(fit, "KS p-value")),
    AD = fmt_p(col_or_na(fit, "Anderson-Darling p-value")),
    Watson = fmt_p(col_or_na(fit, "Watson p-value")),
    Neyman = fmt_p(col_or_na(fit, "Neyman smooth p-value")),
    `TB$_2$` = fmt_p(col_or_na(fit, "TB2 p-value")),
    `Moran's $I$` = fmt_col(col_or_na(fit, "Moran's I"), 3),
    `$p$ (Moran)` = fmt_p(col_or_na(fit, "Moran's p-value")),
    `$|r|>3$ (\\%)` = fmt_col(col_or_na(fit, "Extreme RQR (|r|>3) (%)"), 2),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

build_sampler_table <- function(fam_order, samp, fscalar, floc) {
  rows <- lapply(fam_order, function(f) {
    s <- if (is.null(samp)) NULL else samp[samp$Family == f, , drop = FALSE]
    sc <- if (is.null(fscalar)) NULL else {
      fscalar[fscalar$Family == f & grepl("^P\\(Y", fscalar$Functional), ,
              drop = FALSE]
    }
    lc <- if (is.null(floc)) NULL else floc[floc$Family == f, , drop = FALSE]
    g <- function(d, nm, fun) {
      if (is.null(d) || !nrow(d) || !(nm %in% names(d))) return(NA_real_)
      v <- d[[nm]][is.finite(d[[nm]])]
      if (length(v)) fun(v) else NA_real_
    }
    data.frame(
      Family = f,
      `Draws/chain` = fmt_col(g(s, "Minimum_draws_per_chain", min), 0),
      `max $\\hat R$` = fmt_col(g(s, "Maximum_split_Rhat", max), 3),
      `min ESS` = fmt_col(g(s, "Minimum_ESS_observed", min), 0),
      `max $\\hat R$ ` = fmt_col(g(lc, "Max_Rhat", max), 3),
      `min bulk` = fmt_col(g(lc, "Min_ESS_bulk", min), 0),
      `min tail` = fmt_col(g(lc, "Min_ESS_tail", min), 0),
      `max $\\hat R$  ` = fmt_col(g(sc, "Rhat", max), 3),
      `min bulk ` = fmt_col(g(sc, "ESS_bulk", min), 0),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

build_selection_table <- function(mix, pip, forms, thr = 0.5) {
  fam <- mix$Family
  count <- function(f, grp, nl = FALSE) {
    d <- pip[pip$Family == f & is.finite(pip$PIP_Association) &
               pip$PIP_Association > thr & param_group(pip$Parameter) == grp, ]
    if (!nl) return(nrow(d))
    k <- paste(d$Family, d$Parameter, d$Covariate)
    sum(forms$Selected_Form[match(k, paste(forms$Family, forms$Parameter,
                                           forms$Covariate))] == "nonlinear",
        na.rm = TRUE)
  }
  data.frame(
    Family = fam,
    `Birth` = fmt_col(100 * mix$Birth_accept, 1),
    `Death` = fmt_col(100 * mix$Death_accept, 1),
    `Swap` = fmt_col(100 * mix$Swap_accept, 1),
    `Switches` = fmt_col(mix$Postburn_model_switches, 0),
    `ESS (size)` = fmt_col(mix$ESS_model_size, 1),
    `MCSE assoc.` = fmt_col(mix$Max_association_PIP_MCSE, 3),
    `MCSE NL` = fmt_col(mix$Max_nonlinear_PIP_MCSE, 3),
    `Decided` = ifelse(mix$Selection_threshold_stable, "Yes", "No"),
    `$\\mu$` = vapply(fam, count, numeric(1), grp = "mu"),
    `$\\sigma$` = vapply(fam, count, numeric(1), grp = "sigma"),
    `$\\nu$` = vapply(fam, count, numeric(1), grp = "nu"),
    `$\\tau$` = vapply(fam, count, numeric(1), grp = "tau"),
    `NL` = vapply(fam, function(f) {
      sum(vapply(c("mu", "sigma", "nu", "tau"), function(g) {
        count(f, g, nl = TRUE)
      }, numeric(1)))
    }, numeric(1)),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

build_audit_table <- function(audit, agree = NULL, fam_order = NULL) {
  if (is.null(fam_order)) fam_order <- unique(audit$Family)
  rows <- lapply(fam_order, function(f) {
    a <- audit[audit$Family == f, , drop = FALSE]
    r <- if (is.null(agree)) NULL else agree[agree$Family == f, , drop = FALSE]
    has_rep <- !is.null(r) && nrow(r) > 0L
    data.frame(
      Family = f,
      Included = sum(a$Decision == "Included"),
      Excluded = sum(a$Decision == "Excluded"),
      Inconclusive = sum(a$Decision == "Inconclusive"),
      `Median model` = sum(a$Median_Model),
      BFDR = sum(a$BFDR_Selected),
      Replicates = if (has_rep) as.character(r$Replicates[1L]) else "1",
      `Max. $\\hat R$` = if (has_rep) {
        fmt_col(max(r$Split_Rhat, na.rm = TRUE), 3)
      } else "--",
      Discordant = if (has_rep) sum(!r$Concordant) else "--",
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

build_pip_matrix <- function(pip, forms, grp, fam_order, support = 0.30,
                             thr = 0.5) {
  d <- pip[param_group(pip$Parameter) == grp, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  covs <- unique(d$Covariate)
  keep <- vapply(covs, function(v) {
    any(d$PIP_Association[d$Covariate == v] >= support, na.rm = TRUE)
  }, logical(1))
  covs <- covs[keep]
  if (!length(covs)) return(NULL)
  top <- vapply(covs, function(v) {
    max(d$PIP_Association[d$Covariate == v], na.rm = TRUE)
  }, numeric(1))
  covs <- covs[order(-top)]
  out <- data.frame(Covariate = pretty_cov(covs), stringsAsFactors = FALSE)
  for (f in fam_order) {
    out[[f]] <- vapply(covs, function(v) {
      r <- d[d$Family == f & d$Covariate == v, , drop = FALSE]
      if (!nrow(r) || !is.finite(r$PIP_Association[1L])) return("--")
      s <- formatC(r$PIP_Association[1L], digits = 2, format = "f")
      if (r$PIP_Association[1L] > thr) s <- paste0("\\textbf{", s, "}")
      fm <- forms$Selected_Form[
        forms$Family == f & forms$Parameter == r$Parameter[1L] &
          forms$Covariate == v]
      if (length(fm) && identical(fm[1L], "nonlinear")) {
        s <- paste0(s, "$^{\\mathrm{NL}}$")
      }
      s
    }, character(1))
  }
  out
}

build_effect_table <- function(fam, pip, forms, coefs, evid, thr = 0.5) {
  d <- pip[pip$Family == fam & is.finite(pip$PIP_Association) &
             pip$PIP_Association > thr, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  key <- function(x) paste(x$Family, x$Parameter, x$Covariate, sep = "|")
  d <- d[order(match(param_group(d$Parameter), c("mu", "sigma", "nu", "tau")),
               -d$PIP_Association), , drop = FALSE]
  i <- match(key(d), key(coefs))
  e <- match(key(d), key(evid))
  fm <- forms$Selected_Form[match(key(d), key(forms))]
  ci <- function(lo, hi) {
    ifelse(is.finite(lo) & is.finite(hi), paste0(
      "[", formatC(lo, digits = 2, format = "f"), ", ",
      formatC(hi, digits = 2, format = "f"), "]"), "--")
  }
  data.frame(
    Parameter = param_symbol(d$Parameter),
    Covariate = pretty_cov(d$Covariate),
    Form = ifelse(is.na(fm), "--", fm),
    PIP = fmt_col(d$PIP_Association, 2),
    Mean = fmt_col(coefs$Estimate[i], 3),
    `95\\% CrI` = ci(coefs$Lower_95[i], coefs$Upper_95[i]),
    `$\\hat R$` = fmt_col(coefs$Rhat[i], 3),
    ESS = fmt_col(coefs$ESS[i], 0),
    `Model-averaged mean` = fmt_col(evid$Model_Averaged_Mean[e], 3),
    `95\\% CrI ` = ci(evid$Model_Averaged_Q025[e], evid$Model_Averaged_Q975[e]),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

build_exceedance_table <- function(exc, fam_order) {
  ks <- sort(unique(exc$Threshold))
  obs <- vapply(ks, function(k) exc$ObsFreq[exc$Threshold == k][1L], numeric(1))
  rows <- lapply(fam_order, function(f) {
    cells <- lapply(ks, function(k) {
      r <- exc[exc$Family == f & exc$Threshold == k, , drop = FALSE]
      if (!nrow(r)) return(c("--", "--"))
      c(fmt_col(100 * r$PredFreq, 2), fmt_col(r$Brier, 4))
    })
    out <- data.frame(Family = f, stringsAsFactors = FALSE)
    for (m in seq_along(ks)) {
      out[[paste0("pred", m)]] <- cells[[m]][1L]
      out[[paste0("brier", m)]] <- cells[[m]][2L]
    }
    out
  })
  tab <- do.call(rbind, rows)
  names(tab) <- c("Family", rep(c("Pred. (\\%)", "Brier"), length(ks)))
  head <- paste0(
    "& ", paste0(
      "\\multicolumn{2}{c}{$Y>", ks, "$ (obs. ",
      formatC(100 * obs, digits = 2, format = "f"), "\\%)}",
      collapse = " & "), " \\\\"
  )
  rule <- paste0(
    "\\cmidrule(lr){", 2 * seq_along(ks), "-", 2 * seq_along(ks) + 1, "}",
    collapse = ""
  )
  list(table = tab, pre_header = c(head, rule))
}

build_totals_table <- function(dir_mod, fam_order) {
  rows <- lapply(fam_order, function(f) {
    a <- read_csv_if(file.path(dir_mod, f, "annual_observed_vs_fitted.csv"))
    if (is.null(a)) {
      return(data.frame(Family = f, Fitted = "--", Ratio = "--", MAE = "--",
                        RMSE = "--", stringsAsFactors = FALSE))
    }
    data.frame(
      Family = f, Fitted = fmt_col(sum(a$fitted_fatalities), 0),
      Ratio = fmt_col(sum(a$observed_fatalities) / sum(a$fitted_fatalities), 2),
      MAE = fmt_col(mean(abs(a$prediction_error)), 2),
      RMSE = fmt_col(sqrt(mean(a$prediction_error^2)), 2),
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  names(tab) <- c(
    "Family", "Fitted total", "Obs./fitted", "MAE", "RMSE"
  )
  tab
}

plot_pip_heatmap <- function(pip, forms, fam_order, support = 0.30) {
  d <- pip[is.finite(pip$PIP_Association), , drop = FALSE]
  d$Group <- factor(
    param_group(d$Parameter), c("mu", "sigma", "nu", "tau"),
    c("mu~'/'~lambda", "sigma", "nu~'/'~zeta", "tau")
  )
  keep <- tapply(d$PIP_Association, d$Covariate, max) >= support
  d <- d[d$Covariate %in% names(keep)[keep], , drop = FALSE]
  d$Label <- format_term_label(d$Covariate)
  d$Family <- factor(d$Family, fam_order)
  d$Selected <- d$PIP_Association > 0.5
  ggplot2::ggplot(d, ggplot2::aes(Family, Label, fill = PIP_Association)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(Selected, "*", "")), size = 4
    ) +
    ggplot2::facet_wrap(
      ~ Group, nrow = 1L, scales = "free_y", labeller = ggplot2::label_parsed
    ) +
    ggplot2::scale_fill_gradient(
      low = "white", high = spatial_palette()[["highlight"]],
      limits = c(0, 1), name = "PIP"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "serif"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}

# Writes every table (.tex) and figure that can be built from the CSV outputs
write_result_tables <- function(dir_mod, out_dir = file.path(dir_mod, "tables"),
                                main_family = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  save_tex <- function(lines, name) {
    writeLines(lines, file.path(out_dir, paste0(name, ".tex")))
    message("  wrote ", name, ".tex")
  }
  attempt <- function(name, expr) {
    tryCatch(expr, error = function(e) {
      message("  [tables] ", name, " skipped: ", conditionMessage(e))
    })
  }
  rd <- function(nm) read_csv_if(file.path(dir_mod, nm))
  fit <- rd("model_fit_diagnostics.csv")
  if (is.null(fit)) stop("model_fit_diagnostics.csv not found.", call. = FALSE)
  fam_order <- fit[["Distribution Family"]][order(fit$DIC)]
  if (is.null(main_family)) main_family <- fam_order[1L]
  pip <- rd("pip_mcse.csv")
  forms <- rd("selected_effect_forms.csv")
  coefs <- rd("coefficient_diagnostics.csv")
  evid <- rd("effect_evidence.csv")
  mix <- rd("model_space_mixing_summary.csv")
  samp <- rd("sampler_diagnostics.csv")
  fsc <- rd("functional_diagnostics_scalar.csv")
  floc <- rd("functional_diagnostics_location.csv")
  intv <- rd("predictive_interval_scores.csv")
  exc <- rd("predictive_exceedance_scores.csv")
  
  attempt("family comparison", save_tex(latex_table(
    build_family_comparison(fit, intv),
    paste(
      "Information criteria and in-sample predictive scores by family,",
      "ordered by DIC. Scores use the posterior-mean parameters; bold marks",
      "the best value. RPS$_K$ is the ranked probability score truncated at",
      "the largest observed count. CRPS is reported only when the omitted",
      "tail mass is at most $10^{-10}$ (it is infinite for RSB)."
    ),
    "tab:family_comparison",
    note = paste(
      "IS$_{90}$ and Cov.$_{90}$ are the interval score and coverage of the",
      "central 90\\% predictive interval; they are undefined when the upper",
      "quantile lies beyond the evaluated support."
    )
  ), "tab_family_comparison"))
  
  attempt("calibration", save_tex(latex_table(
    build_calibration_table(fit),
    paste(
      "Uniformity tests ($p$-values) of the randomized PIT values and",
      "residual spatial autocorrelation (Moran's $I$ on district means of",
      "the quantile residuals), by family."
    ),
    "tab:calibration"
  ), "tab_calibration"))
  
  attempt("sampler diagnostics", {
    tab <- build_sampler_table(fam_order, samp, fsc, floc)
    save_tex(latex_table(
      tab,
      paste(
        "Convergence of the post-selection chains. Coefficients are",
        "summarised by split $\\hat R$ and total ESS. The location",
        "parameter of each observation and the exceedance probabilities",
        "$\\mathbb{P}(Y>k)$ are summarised by rank-normalised $\\hat R$ and",
        "bulk and tail ESS."
      ),
      "tab:sampler",
      pre_header = c(
        paste0("& & \\multicolumn{2}{c}{Coefficients} & ",
               "\\multicolumn{3}{c}{Location, by observation} & ",
               "\\multicolumn{2}{c}{$\\mathbb{P}(Y>k)$} \\\\"),
        "\\cmidrule(lr){3-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-9}"
      )
    ), "tab_sampler")
  })
  
  attempt("selection and mixing", {
    save_tex(latex_table(
      build_selection_table(mix, pip, forms),
      paste(
        "Model-space exploration by family. Acceptance rates (\\%) of birth,",
        "death and swap moves, post-burn-in model switches, effective sample",
        "size of the model size, largest Monte Carlo standard errors of the",
        "inclusion probabilities, whether every inclusion decision is",
        "conclusive given the Monte Carlo error (Decided), and number of",
        "terms selected by parameter (NL, nonlinear forms)."
      ),
      "tab:selection_mixing"
    ), "tab_selection_mixing")
  })
  
  attempt("allocation audit", {
    au <- rd("selection_audit.csv")
    if (!is.null(au)) {
      save_tex(latex_table(
        build_audit_table(au, rd("selection_replicates.csv"), fam_order),
        paste(
          "Inclusion decisions by family over all components. A component is",
          "included (excluded) when the Monte Carlo interval of its inclusion",
          "probability lies above (below) 0.5, and inconclusive otherwise.",
          "Median model and BFDR give the numbers of components selected by",
          "the median probability rule and by the Bayesian false discovery",
          "rate rule at level 0.10. The last three columns refer to pooled",
          "independent searches."
        ),
        "tab:allocation_audit"
      ), "tab_allocation_audit")
    }
  })
  
  for (grp in c("mu", "sigma", "nu", "tau")) {
    attempt(paste("PIP", grp), {
      m <- build_pip_matrix(pip, forms, grp, fam_order)
      if (!is.null(m)) {
        save_tex(latex_table(
          m,
          paste0(
            "Posterior inclusion probabilities of the associations in the ",
            switch(grp, mu = "location ($\\mu$; $\\lambda$ for RSB)",
                   sigma = "dispersion ($\\sigma$)",
                   nu = "zero-inflation or shape ($\\nu$; $\\zeta$ for RSB)",
                   tau = "one-inflation ($\\tau$)"),
            " parameter. Bold marks the median probability model ",
            "(PIP$>0.5$) and NL a nonlinear form. Only covariates with PIP ",
            "of at least 0.30 in some family are shown."
          ),
          paste0("tab:pip_", grp), long = TRUE
        ), paste0("tab_pip_", grp))
      }
    })
  }
  
  for (f in fam_order) {
    attempt(paste("effects", f), {
      t <- build_effect_table(f, pip, forms, coefs, evid)
      if (!is.null(t)) {
        save_tex(latex_table(
          t,
          paste0(
            f, ": selected associations. Mean, 95\\% credible interval, ",
            "$\\hat R$ and ESS come from the post-selection chains; the ",
            "model-averaged summary comes from the retained states of the ",
            "model search and includes the draws in which the term is absent."
          ),
          paste0("tab:effects_", f), long = TRUE,
          note = paste(
            "Coefficients are on the link scale of each parameter;",
            "nonlinear terms are shown in the smooth-effect figures."
          )
        ), paste0("tab_effects_", f))
      }
    })
  }
  
  attempt("tail calibration", {
    e <- build_exceedance_table(exc, fam_order)
    save_tex(latex_table(
      e$table,
      paste(
        "Predicted share of district-months with more than $k$ fatalities and",
        "Brier score of the exceedance forecast, by family. The observed",
        "share is in the column heading."
      ),
      "tab:exceedance", pre_header = e$pre_header
    ), "tab_exceedance")
  })
  
  attempt("annual totals", {
    save_tex(latex_table(
      build_totals_table(dir_mod, fam_order),
      paste(
        "Agreement between observed and fitted annual fatalities by district.",
        "Fitted totals use the conditional mean, which is infinite for RSB",
        "and therefore not shown."
      ),
      "tab:totals"
    ), "tab_totals")
  })
  
  attempt("PIP heatmap", {
    ggplot2::ggsave(
      file.path(out_dir, "Figure_PIP_Heatmap.pdf"),
      plot_pip_heatmap(pip, forms, fam_order), width = 13, height = 7,
      device = grDevices::cairo_pdf
    )
  })
  attempt("ESS growth", {
    cv <- rd("functional_diagnostics_curve.csv")
    if (!is.null(cv)) {
      ggplot2::ggsave(
        file.path(out_dir, "Figure_ESS_Growth.pdf"),
        ggplot2::ggplot(cv, ggplot2::aes(
          Draws_per_chain, ESS_bulk, colour = Functional)) +
          ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
          ggplot2::facet_wrap(~ Family, scales = "free_y") +
          ggplot2::labs(x = "Draws per chain", y = "Bulk ESS") +
          ggplot2::theme_bw(base_size = 10),
        width = 12, height = 7, device = grDevices::cairo_pdf
      )
    }
  })
  invisible(out_dir)
}

#------------------------------------------------------------------------------
# Audit of the covariate allocation
#------------------------------------------------------------------------------

# Batch-means mean, MCSE and ESS of the columns of a 0/1 matrix
batch_means_stats <- function(X, n_batch = NULL) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  na <- rep(NA_real_, p)
  if (n < 20L) return(list(mean = colMeans(X), mcse = na, ess = na))
  if (is.null(n_batch)) n_batch <- max(2L, floor(sqrt(n)))
  m <- floor(n / n_batch)
  idx <- rep(seq_len(n_batch), each = m)
  bm <- rowsum(X[seq_len(m * n_batch), , drop = FALSE], idx) / m
  pm <- colMeans(X)
  mcse <- apply(bm, 2L, stats::sd) / sqrt(n_batch)
  ess <- ifelse(mcse > 0, pmin(n, pm * (1 - pm) / mcse^2), NA_real_)
  list(mean = pm, mcse = mcse, ess = ess)
}

# Stack the indicator trace of a search in a draws x components matrix, and
# describe its columns
search_indicator_matrix <- function(search) {
  tr <- search$trace_inclusion
  first <- tr[[1L]]
  X <- do.call(rbind, lapply(tr, function(z) as.numeric(z)))
  desc <- data.frame(
    Key = rep(rownames(first), times = ncol(first)),
    Parameter = rep(colnames(first), each = nrow(first)),
    stringsAsFactors = FALSE
  )
  list(X = X, index = desc)
}

# Attempts that could have changed a component, per cell
attempt_counts <- function(search, state) {
  d <- search$diagnostics
  pick <- if (state == "in") {
    c("candidate_death_attempted", "candidate_swap_out_attempted")
  } else {
    c("candidate_birth_attempted", "candidate_swap_in_attempted")
  }
  mats <- Filter(Negate(is.null), d[pick])
  if (!length(mats)) return(NULL)
  Reduce(`+`, lapply(mats, function(z) as.matrix(z) * 1))
}

# Allocation audit of a search, one row per component and parameter.
# A component is Included (Excluded) when the Monte Carlo interval of its
# posterior inclusion probability lies above (below) the cut, and
# Inconclusive otherwise. A chain that never moves has no Monte Carlo
# error; it is called Included (Excluded) only if the sampler tried to
# change it at least `min_attempts` times.
selection_audit <- function(
    search, threshold = 0.5, z = 1.96, alpha = 0.10, min_attempts = 50L,
    n_batch = NULL) {
  ind <- search_indicator_matrix(search)
  st <- batch_means_stats(ind$X, n_batch)
  out <- ind$index
  out$PIP <- unname(st$mean)
  out$MCSE <- unname(st$mcse)
  out$ESS <- unname(st$ess)
  
  adm <- search$admissibility
  ok <- rep(TRUE, nrow(out))
  if (!is.null(adm)) {
    ok <- mapply(
      function(k, p) isTRUE(adm[k, p]), out$Key, out$Parameter
    )
  }
  out <- out[ok, , drop = FALSE]
  
  cell <- function(m, k, p) {
    if (is.null(m)) return(rep(NA_real_, length(k)))
    as.numeric(m[cbind(match(k, rownames(m)), match(p, colnames(m)))])
  }
  n_in <- attempt_counts(search, "in")
  n_out <- attempt_counts(search, "out")
  out$Change_Attempts <- ifelse(
    out$PIP >= 0.5, cell(n_in, out$Key, out$Parameter),
    cell(n_out, out$Key, out$Parameter)
  )
  
  fixed <- out$PIP <= 0 | out$PIP >= 1
  lower <- out$PIP - z * out$MCSE
  upper <- out$PIP + z * out$MCSE
  dec <- rep("Inconclusive", nrow(out))
  dec[!fixed & is.finite(lower) & lower > threshold] <- "Included"
  dec[!fixed & is.finite(upper) & upper < threshold] <- "Excluded"
  tested <- fixed & is.finite(out$Change_Attempts) &
    out$Change_Attempts >= min_attempts
  dec[tested & out$PIP >= 1] <- "Included"
  dec[tested & out$PIP <= 0] <- "Excluded"
  out$Fixed_Chain <- fixed
  out$Decision <- dec
  out$Median_Model <- out$PIP > threshold
  
  # Bayesian FDR rule: largest set whose mean posterior null probability
  # does not exceed alpha
  ord <- order(out$PIP, decreasing = TRUE)
  fdr <- cumsum(1 - out$PIP[ord]) / seq_along(ord)
  take <- which(fdr <= alpha)
  cut <- if (length(take)) max(take) else 0L
  out$BFDR_Selected <- FALSE
  if (cut > 0L) out$BFDR_Selected[ord[seq_len(cut)]] <- TRUE
  out$Expected_FDR <- NA_real_
  out$Expected_FDR[ord] <- fdr
  
  out$Level <- ifelse(startsWith(out$Key, "NL::"), "nonlinear", "association")
  out$Variable <- sub("^(NL|L)::", "", out$Key)
  rownames(out) <- NULL
  out
}

# Agreement of independent searches of one family
selection_replicates <- function(searches, threshold = 0.5, z = 1.96) {
  R <- length(searches)
  if (R < 2L) return(data.frame())
  mats <- lapply(searches, search_indicator_matrix)
  idx <- mats[[1L]]$index
  n <- min(vapply(mats, function(m) nrow(m$X), numeric(1)))
  h <- floor(n / 2L)
  if (h < 10L) return(data.frame())
  
  chains <- unlist(lapply(mats, function(m) list(
    m$X[seq_len(h), , drop = FALSE],
    m$X[(h + 1L):(2L * h), , drop = FALSE]
  )), recursive = FALSE)
  rhat <- vapply(seq_len(ncol(chains[[1L]])), function(j) {
    M <- vapply(chains, function(ch) ch[, j], numeric(h))
    W <- mean(apply(M, 2L, stats::var))
    B <- h * stats::var(colMeans(M))
    if (W <= 0) return(if (B <= 0) 1 else Inf)
    sqrt(((h - 1) / h * W + B / h) / W)
  }, numeric(1))
  
  stats_by_rep <- lapply(mats, function(m) batch_means_stats(m$X))
  pip <- vapply(stats_by_rep, function(s) s$mean, numeric(ncol(chains[[1L]])))
  mc <- vapply(stats_by_rep, function(s) s$mcse, numeric(ncol(chains[[1L]])))
  if (is.null(dim(pip))) {
    pip <- matrix(pip, ncol = R)
    mc <- matrix(mc, ncol = R)
  }
  spread <- apply(pip, 1L, function(v) diff(range(v)))
  se <- sqrt(rowSums(mc^2, na.rm = TRUE))
  se[se <= 0] <- NA_real_
  decision <- pip > threshold
  
  out <- idx
  out$Replicates <- R
  out$PIP_Min <- apply(pip, 1L, min)
  out$PIP_Max <- apply(pip, 1L, max)
  out$PIP_Spread <- spread
  out$Spread_Over_MCSE <- spread / se
  out$Split_Rhat <- rhat
  out$Concordant <- rowSums(decision) %in% c(0L, R)
  out$Level <- ifelse(startsWith(out$Key, "NL::"), "nonlinear", "association")
  out$Variable <- sub("^(NL|L)::", "", out$Key)
  adm <- searches[[1L]]$admissibility
  if (!is.null(adm)) {
    ok <- mapply(function(k, p) isTRUE(adm[k, p]), out$Key, out$Parameter)
    out <- out[ok, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

# Merge independent searches of one family into a single search object
pool_search_fits <- function(searches, max_rhat = 1.1, max_spread = 3) {
  searches <- Filter(Negate(is.null), searches)
  R <- length(searches)
  if (R == 1L) return(searches[[1L]])
  
  base <- searches[[1L]]
  tr <- unlist(lapply(searches, `[[`, "trace_inclusion"), recursive = FALSE)
  base$trace_inclusion <- tr
  base$trace_logpost <- unlist(lapply(searches, `[[`, "trace_logpost"))
  base$trace_model_size <- unlist(lapply(searches, `[[`, "trace_model_size"))
  base$retained_states <- unlist(
    lapply(searches, `[[`, "retained_states"), recursive = FALSE
  )
  reps <- unlist(
    lapply(searches, `[[`, "representative_states"), recursive = FALSE
  )
  base$representative_states <- reps[!duplicated(names(reps))]
  
  comp <- Reduce("+", lapply(tr, function(z) z * 1)) / length(tr)
  adm <- base$admissibility
  if (!is.null(adm)) comp[!adm] <- NA_real_
  pips <- aggregate_component_pips(comp)
  base$component_inclusion_probabilities <- comp
  base$posterior_inclusion_probabilities <- pips$association
  base$posterior_nonlinear_probabilities <- pips$nonlinear
  base$posterior_nonlinear_given_association <-
    pips$nonlinear_given_association
  
  agree <- selection_replicates(searches)
  audit <- selection_audit(base)
  d <- base$diagnostics
  ind <- aggregate_trace_indicators(tr)
  for (lev in c("association", "nonlinear")) {
    dg <- indicator_trace_diagnostics(ind[[lev]])
    key <- if (lev == "association") "association_pip_mcse" else
      "nonlinear_pip_mcse"
    d[[key]] <- dg$mcse
    if (lev == "nonlinear") d$nonlinear_ess <- dg$ess
  }
  d$pooled_replicates <- R
  d$kernel_operational <- all(vapply(
    searches, function(s) isTRUE(s$diagnostics$kernel_operational), logical(1)
  ))
  reasons <- character(0)
  if (!d$kernel_operational) reasons <- c(reasons, "kernel_operational")
  if (nrow(agree) && any(agree$Split_Rhat > max_rhat, na.rm = TRUE)) {
    reasons <- c(reasons, "replicate_split_rhat")
  }
  if (nrow(agree) &&
      any(agree$Spread_Over_MCSE > max_spread, na.rm = TRUE)) {
    reasons <- c(reasons, "replicate_pip_spread")
  }
  if (nrow(agree) && any(!agree$Concordant)) {
    reasons <- c(reasons, "replicate_discordant_decision")
  }
  if (any(audit$Decision == "Inconclusive")) {
    reasons <- c(reasons, "inconclusive_components")
  }
  d$mixing_within_reference <- length(reasons) == 0L
  d$mixing_diagnostics_outside_reference <- reasons
  d$replicate_agreement <- agree
  base$diagnostics <- d
  base
}

# Selected components per family in a long audit table: how many families
# include each component
selection_concordance <- function(audit_all) {
  if (!nrow(audit_all)) return(data.frame())
  key <- paste(audit_all$Level, audit_all$Variable, audit_all$Parameter)
  split_by <- split(seq_len(nrow(audit_all)), key)
  do.call(rbind, lapply(split_by, function(i) {
    a <- audit_all[i, , drop = FALSE]
    data.frame(
      Level = a$Level[1L], Variable = a$Variable[1L],
      Parameter = a$Parameter[1L], Families = nrow(a),
      Included = sum(a$Decision == "Included"),
      Excluded = sum(a$Decision == "Excluded"),
      Inconclusive = sum(a$Decision == "Inconclusive"),
      Median_Model = sum(a$Median_Model),
      BFDR_Selected = sum(a$BFDR_Selected),
      Mean_PIP = mean(a$PIP), Min_PIP = min(a$PIP), Max_PIP = max(a$PIP),
      stringsAsFactors = FALSE
    )
  }))
}
