### Racing Diffusion Model

### Single accumulator

# pigt, digt, rwaldt Copyright (C) 2013  Trisha Van Zandt distributed with:
# Logan, Van Zandt, Verbruggen, and Wagenmakers (2014).  On the ability to
# inhibit thought and action: General and special theories of an act of control.
# Psychological Review.

# Comments and changes added by Andrew Heathcote and Gabriel Tillman

rWald <- function(n, B, v, A)
  # random deviate function for single acumulator
{
  rwaldt <- function(n, k, l, tiny = 1e-6) {
    # random sample of n from a Wald (or Inverse Gaussian)
    # k = criterion, l = rate, assumes sigma=1 Browninan motion
    # about same speed as statmod rinvgauss
    
    rlevy <- function(n = 1, m = 0, c = 1) {
      if (any(c < 0))
        stop("c must be positive")
      c / qnorm(1 - runif(n) / 2) ^ 2 + m
    }
    
    flag     <- l > tiny
    x        <- rep(NA, times = n)
    x[!flag] <- rlevy(sum(!flag), 0, k[!flag] ^ 2)
    mu       <- k / l
    lambda   <- k ^ 2
    y        <- rnorm(sum(flag)) ^ 2
    mu.0     <- mu[flag]
    lambda.0 <- lambda[flag]
    x.0      <- mu.0 + mu.0 ^ 2 * y / (2 * lambda.0) -
      sqrt(4 * mu.0 * lambda.0 * y + mu.0 ^ 2 * y ^ 2) * mu.0 /
      (2 * lambda.0)
    
    z             <- runif(length(x.0))
    test          <- mu.0 / (mu.0 + x.0)
    x.0[z > test] <- mu.0[z > test] ^ 2 / x.0[z > test]
    x[flag]       <- x.0
    x[x < 0]      <- max(x)
    x
  }
  
  # Act as if negative v never terminates, cluge to do single accumulator
  # case by passing negative v
  if (length(v) != n)
    v <- rep(v, length.out = n)
  if (length(B) != n)
    B <- rep(B, length.out = n)
  if (length(A) != n)
    A <- rep(A, length.out = n)
  
  # Kluge to return -Inf for negative rates, so can implment one accumulator case
  out      <- numeric(n)
  ok       <- v > 0
  nok      <- sum(ok)
  bs       <- B[ok] + runif(nok, 0, A[ok])
  out[ok]  <- rwaldt(nok, k = bs, l = v[ok])
  out[!ok] <- Inf
  out
}


dWald <- function(t, v, B, A)
  # density for single accumulator
{
  digt <- function(t, k = 1, l = 1, a = .1, tiny = 1e-10) {
    # pdf of inverse gaussian at t with k +/- a/2 uniform variability
    # returns digt.0 if a<1e-10
    
    digt.0 <- function(t, k = 1, l = 1) {
      # pdf of inverse gaussian at t with no k variability
      # much faster than statmod's dinvgauss funciton
      
      lambda <- k ^ 2
      l0     <- l == 0
      e      <- numeric(length(t))
      if (any(!l0)) {
        mu     <- k[!l0] / l[!l0]
        e[!l0] <-
          -(lambda[!l0] / (2 * t[!l0])) * (t[!l0] ^ 2 / mu ^ 2 - 2 * t[!l0] / mu  + 1)
      }
      if (any(l0))
        e[l0]   <- -.5 * lambda[l0] / t[l0]
      x         <- exp(e + .5 * log(lambda) - .5 * log(2 * t ^ 3 * pi))
      x[t <= 0] <- 0
      x
    }
    
    options(warn = -1)
    if (length(k) != length(t))
      k <- rep(k, length.out = length(t))
    if (length(l) != length(t))
      l <- rep(l, length.out = length(t))
    if (length(a) != length(t))
      a <- rep(a, length.out = length(t))
    
    tpos <- t <= 0
    
    atiny    <- a <= tiny & !tpos
    a[atiny] <- 0
    
    ltiny        <- (l <= tiny) & !atiny & !tpos
    notltiny     <- (l > tiny) & !atiny & !tpos
    l[l <= tiny] <- 0
    x            <- numeric(length(t))
    
    # No threshold variability
    if (any(atiny))
      x[atiny] <- digt.0(t = t[atiny], k = k[atiny], l = l[atiny])
    
    # Threshold variability
    if (any(!atiny)) {
      if (any(notltiny)) {
        # rate non-zero
        
        sqr.t   <- sqrt(t[notltiny])
        
        term.1a <-
          -(a[notltiny] - k[notltiny] + t[notltiny] * l[notltiny]) ^ 2 / (2 * t[notltiny])
        term.1b <-
          -(a[notltiny] + k[notltiny] - t[notltiny] * l[notltiny]) ^ 2 / (2 * t[notltiny])
        term.1  <-
          (exp(term.1a) - exp(term.1b)) / sqrt(2 * pi * t[notltiny])
        
        term.2a <- log(.5) + log(l[notltiny])
        term.2b <-
          2 * pnorm((-k[notltiny] + a[notltiny]) / sqr.t + sqr.t * l[notltiny]) -
          1
        term.2c <-
          2 * pnorm((k[notltiny] + a[notltiny]) / sqr.t - sqr.t * l[notltiny]) -
          1
        term.2d <- term.2b + term.2c
        term.2  <- exp(term.2a) * term.2d
        
        term.3  <- term.1 + term.2
        term.4  <- log(term.3) - log(2) - log(a[notltiny])
        x[notltiny] <- exp(term.4)
      }
      
      if (any(ltiny)) {
        # rate zero
        log.t    <- log(t[ltiny])
        term.1   <- -.5 * (log(2) + log(pi) + log.t)
        term.2   <- (k[ltiny] - a[ltiny]) ^ 2 / (2 * t[ltiny])
        term.3   <- (k[ltiny] + a[ltiny]) ^ 2 / (2 * t[ltiny])
        term.4   <- (exp(-term.2) - exp(-term.3))
        term.5   <- term.1 + log(term.4) - log(2) - log(a[ltiny])
        x[ltiny] <- exp(term.5)
      }
      
    }
    
    x[x < 0 | is.nan(x)] <- 0
    x
  }
  
  out      <- numeric(length(t))
  ok       <- v > 0
  out[ok]  <- digt(t[ok],
                   k = B[ok] + A[ok] / 2,
                   l = v[ok],
                   a = A[ok] / 2)
  out[!ok] <- 0
  out
}


# # Check density and normalization

# # Without Starting Point Variability
# par(mfrow=c(1,2))
# n=1e6; A=0; v=2; B=1
# integrate(dWald,lower=0,upper=Inf,A=A,v=v,B=B)
# sim <- rWald(n=n,A=A,v=v,B=B)
# bad <- sim>10; mean(bad) # if this is large the check isnt valid
# dns <- density(sim[!bad])
# x <- dns$x[dns$x>0]
# d <- dWald(x,A=A,v=v,B=B)
# plot(dns)
# lines(x,d,col="red")

# # With Starting Point Variability
# n=1e6; A=1; v=2; B=1
# integrate(dWald,lower=0,upper=Inf,A=A,v=v,B=B)
# sim <- rWald(n=n,A=A,v=v,B=B)
# bad <- sim>10; mean(bad) # if this is large the check isnt valid
# dns <- density(sim[!bad])
# x <- dns$x[dns$x>0]
# d <- dWald(x,A=A,v=v,B=B)
# plot(dns)
# lines(x,d,col="red")

pWald <- function(t, v, B, A)
  # cumulative density for single accumulator
{
  pigt <- function(t, k = 1, l = 1, a = .1, tiny = 1e-10) {
    # cdf of inverse gaussian at t with k +/- a/2 uniform variability
    # returns pigt.0 if a<=0
    
    pigt.0 <- function(t, k = 1, l = 1) {
      # cdf of inverse gaussian at t with no k variability
      # much faster than statmod's pinvgauss funciton
      
      mu       <- k / l
      lambda   <- k ^ 2
      e        <- exp(log(2 * lambda) - log(mu))
      add      <- sqrt(lambda / t) * (1 + t / mu)
      sub      <- sqrt(lambda / t) * (1 - t / mu)
      p.1      <- 1 - pnorm(add)
      p.2      <- 1 - pnorm(sub)
      x        <- exp(e + log(p.1)) + p.2
      
      x[t < 0] <- 0
      x
    }
    
    options(warn = -1)
    if (length(k) != length(t))
      k  <- rep(k, length.out = length(t))
    if (length(l) != length(t))
      l  <- rep(l, length.out = length(t))
    if (length(a) != length(t))
      a  <- rep(a, length.out = length(t))
    
    tpos <- t <= 0
    
    atiny    <- a <= tiny & !tpos
    a[atiny] <- 0
    
    ltiny        <- (l <= tiny) & !atiny & !tpos
    notltiny     <- (l > tiny) & !atiny & !tpos
    l[l <= tiny] <- 0
    
    x <- numeric(length(t))
    
    # No threshold variability
    if (any(atiny))
      x[atiny] <- pigt.0(t[atiny], k[atiny], l[atiny])
    
    # Threshold variability
    if (any(!atiny)) {
      if (any(notltiny)) {
        # rate non-zero
        
        log.t   <- log(t[notltiny])
        sqr.t   <- sqrt(t[notltiny])
        
        term.1a <- .5 * log.t - .5 * log(2 * pi)
        term.1b <-
          exp(-((k[notltiny] - a[notltiny] - t[notltiny] * l[notltiny]) ^ 2 / t[notltiny]) /
                2)
        term.1c <-
          exp(-((k[notltiny] + a[notltiny] - t[notltiny] * l[notltiny]) ^ 2 / t[notltiny]) /
                2)
        term.1  <- exp(term.1a) * (term.1b - term.1c)
        
        term.2a <- exp(2 * l[notltiny] * (k[notltiny] - a[notltiny]) +
                         log(pnorm(-(
                           k[notltiny] - a[notltiny] + t[notltiny] * l[notltiny]
                         ) / sqr.t)))
        term.2b <- exp(2 * l[notltiny] * (k[notltiny] + a[notltiny]) +
                         log(pnorm(-(
                           k[notltiny] + a[notltiny] + t[notltiny] * l[notltiny]
                         ) / sqr.t)))
        term.2  <- a[notltiny] + (term.2b - term.2a) / (2 * l[notltiny])
        
        term.4a <-
          2 * pnorm((k[notltiny] + a[notltiny]) / sqr.t - sqr.t * l[notltiny]) -
          1
        term.4b <-
          2 * pnorm((k[notltiny] - a[notltiny]) / sqr.t - sqr.t * l[notltiny]) -
          1
        term.4c <-
          .5 * (t[notltiny] * l[notltiny] - a[notltiny] - k[notltiny] + .5 / l[notltiny])
        term.4d <-
          .5 * (k[notltiny] - a[notltiny] - t[notltiny] * l[notltiny] - .5 / l[notltiny])
        term.4  <- term.4c * term.4a + term.4d * term.4b
        
        x[notltiny] <- (term.4 + term.2 + term.1) / (2 * a[notltiny])
      }
      
      if (any(ltiny)) {
        # rate zero
        sqr.t   <- sqrt(t[ltiny])
        log.t   <- log(t[ltiny])
        term.5a <- 2 * pnorm((k[ltiny] + a[ltiny]) / sqr.t) - 1
        term.5b <- 2 * pnorm(-(k[ltiny] - a[ltiny]) / sqr.t) - 1
        term.5  <-
          (-(k[ltiny] + a[ltiny]) * term.5a - (k[ltiny] - a[ltiny]) * term.5b) /
          (2 * a[ltiny])
        
        term.6a <-
          -.5 * (k[ltiny] + a[ltiny]) ^ 2 / t[ltiny] - .5 * log(2) - .5 * log(pi) + .5 *
          log.t - log(a[ltiny])
        term.6b <-
          -.5 * (k[ltiny] - a[ltiny]) ^ 2 / t[ltiny] - .5 * log(2) - .5 * log(pi) + .5 *
          log.t - log(a[ltiny])
        term.6 <- 1 + exp(term.6b) - exp(term.6a)
        
        x[ltiny] <- term.5 + term.6
      }
      
    }
    
    x[x < 0 | is.nan(x)] <- 0
    x
  }
  
  out <- numeric(length(t))
  ok <- v > 0
  out[ok] <- pigt(t[ok],
                  k = B[ok] + A[ok] / 2,
                  l = v[ok],
                  a = A[ok] / 2)
  out[!ok] <- 0
  out
  
}


# # Check cumulative density
# par(mfrow=c(1,2))
# 
# n=1e6; A=0; v=2; B=1
# sim <- rWald(n=n,A=A,v=v,B=B)
# probs=1:990/1000
# qs <- quantile(sim,probs=probs)
# cd <- pWald(qs,A=A,v=v,B=B)
# plot(qs,probs,type="l")
# lines(qs,cd,col="red")
# 
# n=1e6; A=1; v=2; B=1
# sim <- rWald(n=n,A=A,v=v,B=B)
# probs=1:990/1000
# qs <- quantile(sim,probs=probs)
# cd <- pWald(qs,A=A,v=v,B=B)
# plot(qs,probs,type="l")
# lines(qs,cd,col="red")

### Race model

rWaldRace <- function(n, v, B, A, t0, gf = 0, return.ttf = FALSE)
  # random function for Wald race.
{
  B[B < 0] <- 0 # Protection for negatives
  A[A < 0] <- 0
  n_v      <- ifelse(is.null(dim(v)), length(v), dim(v)[1])
  ttf      <- matrix(t0 + rWald(n * n_v, B = B, v = v, A = A), nrow = n_v)
  if (return.ttf)
    return(ttf)
  resp <- apply(ttf, 2, which.min)
  out  <-
    data.frame(RT = ttf[cbind(resp, 1:n)], R = apply(ttf, 2, which.min))
  
  if (gf[1] > 0) {
    is.gf <- as.logical(rbinom(dim(out)[1], 1, gf))
    out$RT[is.gf] <- NA
    out$R[is.gf] <- 1
  }
  
  out
}


n1Wald <- function(dt, v, B, A, t0, gf = 0)
  # Generates defective PDF for responses on node=1, dt (decison time) is a vector of times
{
  B[B < 0] <- 0 # Protection for negatives
  A[A < 0] <- 0
  n_acc <- ifelse(is.null(dim(v)), length(v), dim(v)[1])
  if (is.null(dim(dt)))
    dt <- matrix(rep(dt, each = n_acc), nrow = n_acc)
  dt <- dt - t0
  
  is.go <- !is.na(dt[1, ])
  n.go <- sum(is.go)
  
  if (!is.matrix(v))
    v <- matrix(rep(v, n.go), nrow = n_acc)
  if (!is.matrix(B))
    B <- matrix(rep(B, n.go), nrow = n_acc)
  if (!is.matrix(A))
    A <- matrix(rep(A, n.go), nrow = n_acc)
  
  # Winner
  dt[1, is.go] <-
    (1 - gf[1]) * dWald(dt[1, is.go], A = A[1, ], v = v[1, ], B = B[1, ])
  if (n_acc > 1)
    for (i in 2:n_acc)
      dt[1, is.go] <-
    dt[1, is.go] * (1 - pWald(dt[i, is.go], A = A[i, ], v = v[i, ], B = B[i, ]))
  
  dt[1, !is.go] <- gf[1]
  
  dt[1, ]
}




#### Alternative Code

# Code here was written by Gabriel Tillman. Simplified
# the density and probability function, while allowing for the specification
# of the diffusion coefficient (i.e., d). The diffusion coeffient parameter
# is available for the dwald_alt function but not the pwald_alt function. 
# Still working on math for that, when start point variability is included.

dwald_alt <- function(t, v, b, A, d) {
  if(A==0) A <- 1e-5
  s = (1 / sqrt(t)) / d
  alpha = (b - A - (t * v)) / (t * s)
  beta  = (b - (t * v)) / (t * s)
  1 / A * (-v * pnorm(alpha) + s * dnorm(alpha) +
             v * pnorm(beta) - s * dnorm(beta))
}

pwald_alt <- function(t, v, b, A) {
  require(pracma)
  if(A==0) A <- 1e-5 
  B = b - A
  (
    A / 2 + 0.3989422804014327 * (exp(1) ^ (-((B - t * v) ^ 2 / (2 * t))) -
        exp(1) ^ (-((A + B - t * v) ^ 2 / (2 * t)))) * 
      t ^ 0.5 + 0.5 * (B - 0.5 / v - t * v) *
      (-1 + erfc((-(B / sqrt(t)) + sqrt(t) * v) / sqrt(2))) +
      0.5 * (-A - B + 0.5 / v + t * v) * (-1 + erfc((-((A + B) / sqrt(t)) + 
                                                       sqrt(t) * v) /
        sqrt(2))) + ((-(1 / 2)) * exp(1) ^ (2 * B * v) * 
                       erfc((B + t * v) / (sqrt(2) * sqrt(t))) +
                       (1 / 2) * exp(1) ^ (2 * (A + B) * v) * 
                       erfc((A + B + t * v) / (sqrt(2) * sqrt(t)))) / 
      (2 * v)) / A
}


# Check dwald_alt integrates to 1
# With Starting Point Variability
# n=1e6; A=0; v=2; B=1
# integrate(dwald_alt,lower=0,upper=Inf,A=A,v=v,b=B+A,d=1)
# sim <- rWald(n=n,A=A,v=v,B=B)
# bad <- sim>10; mean(bad) # if this is large the check isnt valid
# dns <- density(sim[!bad])
# x <- dns$x[dns$x>0]
# d <- dWald(x,A=A,v=v,B=B)
# plot(dns)
# lines(x,d,col="red")


# Check pwald_alt density 
# n=1e6; A=0; v=2; B=1
# sim <- rWald(n=n,A=A,v=v,B=B)
# probs=1:990/1000
# qs <- quantile(sim,probs=probs)
# cd <- pwald_alt(qs,A=A,v=v,b=B+A)
# plot(qs,probs,type="l")
# lines(qs,cd,col="red")
