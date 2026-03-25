scDblFinder_simple <- function(sce) {

  if (ncol(sce) < 100)
    warning("scDblFinder might not work well with very low numbers of cells.")

  k <- scDblFinder:::.defaultKnnKs(NULL, ncol(sce))

  # Feature selection
  sel_features <- row.names(sce)
  if (length(sel_features) > 1352)
    sel_features <- scDblFinder::selFeatures(sce[sel_features, ], NULL, nfeatures = 1352, propMarkers = 0)
  sce <- sce[sel_features, ]

  # Generate artificial doublets
  artificialDoublets <- min(25000, max(1500, ceiling(ncol(sce) * 0.8)))
  message("Creating ~", artificialDoublets, " artificial doublets...")

  ad <- scDblFinder::getArtificialDoublets(
    SingleCellExperiment::counts(sce),
    n = artificialDoublets,
    clusters = NULL,
    propRandom = 0
  )

  ado <- ad$origins
  ad <- ad$counts

  # Setup classification matrix
  src <- factor(rep(1:2, c(ncol(sce), ncol(ad))), labels = c("real", "artificial"))
  ctype <- factor(rep(c(1L, 2L), c(ncol(sce), ncol(ad))), labels = c("real", "doublet"))
  inclInTrain <- rep(c(TRUE, TRUE), c(ncol(sce), ncol(ad)))

  # Combine counts
  ad <- convert_matrix_type(ad)
  e <- cbind(SingleCellExperiment::counts(sce), ad[row.names(sce), ])
  matrix_dir <- tempfile()
  e <- write_matrix_dir(e, matrix_dir)

  lsizes <- Matrix::colSums(e)
  nfeatures <- Matrix::colSums(e > 0L)
  nAbove2 <- Matrix::colSums(e > 2L)

  # Calculate scores
  message("Running cxds2...")
  cxds_score <- cxds2(e, whichDbls = which(ctype == "doublet"))

  # PCA
  pca <- .defaultProcessing(e, dims = 20)

  # Evaluate kNN
  origins <- as.factor(c(rep(NA, ncol(sce)), as.character(ado)))
  d <- scDblFinder:::.evaluateKNN(
    pca,
    ctype,
    origins,
    expected = NULL,
    k = k)$d

  d$lsizes <- lsizes
  d$nfeatures <- nfeatures
  d$nAbove2 <- nAbove2
  d$src <- src
  d$cxds_score <- cxds_score
  d$include.in.training <- inclInTrain
  d$cluster <- NULL

  # Score doublets
  d <- scDblFinder:::.scDblscore(d, scoreType = "xgb", addVals = pca[, 1:19, drop = FALSE],
                                 threshold = TRUE, dbr = NULL, dbr.sd = NULL,
                                 nrounds = 0.25, dbr.per1k = 0.008, max_depth = 4,
                                 iter = 3, features = NULL, verbose = TRUE,
                                 metric = "logloss", filterUnidentifiable = TRUE,
                                 unident.th = 0.2, BPPARAM = SerialParam(progressbar = TRUE))

  # Return annotated SCE
  rowData(sce)$scDblFinder.selected <- TRUE
  scDblFinder:::.scDblAddCD(sce, d)
}

.defaultProcessing <- function(e, dims = NULL, doNorm = NULL) {
  if (is.null(doNorm))
    doNorm <- ncol(e) <= 50000
  if (doNorm) {
    tryCatch({
      e <- normalizeCounts(e)
    }, error = function(er) {
      warning("Error in calculating norm factors:", er)
    })
  }
  if (is.null(dims))
    dims <- 20

  # built-in support for BPCells
  pca <- Seurat::RunPCA(e, npcs = 20, verbose = FALSE)
  pca <- pca@cell.embeddings
  colnames(pca) <- paste0('PC', seq_len(dims))

  pca
}


createDoublets <- function(x, dbl.idx, clusters=NULL, resamp=0.5,
                           halfSize=0.5, adjustSize=FALSE, prefix="dbl."){
  adjustSize <- scDblFinder:::.checkPropArg(as.numeric(adjustSize),FALSE)
  halfSize <- scDblFinder:::.checkPropArg(as.numeric(halfSize),FALSE)
  resamp <- scDblFinder:::.checkPropArg(as.numeric(resamp),FALSE)
  if(is(x,"SingleCellExperiment")) x <- counts(x)
  if(adjustSize>1 || adjustSize<0)
    stop("`adjustSize` should be a logical or a number between 0 and 1.")
  if(halfSize>1 || halfSize<0)
    stop("`adjustSize` should be a logical or a number between 0 and 1.")
  wAd <- sample.int(nrow(dbl.idx), size=round(adjustSize*nrow(dbl.idx)))
  wNad <- setdiff(seq_len(nrow(dbl.idx)),wAd)

  xa <- x[,dbl.idx[wNad,1],drop=FALSE]
  xb <- x[,dbl.idx[wNad,2],drop=FALSE]
  # these are now small and BPCells doesn't support '+'
  xa <- as(xa, 'dgCMatrix')
  xb <- as(xb, 'dgCMatrix')
  x1 <- xa + xb

  if(length(wAd)>1){
    if(is.null(clusters)) stop("If `adjustSize=TRUE`, clusters must be given.")
    dbl.idx <- as.data.frame(dbl.idx[wAd,,drop=FALSE])
    ls <- Matrix::colSums(x)
    csz <- vapply(split(ls,clusters), FUN=median, FUN.VALUE=numeric(1))
    dbl.idx$ls.ratio <- ls[dbl.idx[,1]]/(ls[dbl.idx[,1]]+ls[dbl.idx[,2]])
    ls1 <- csz[as.character(clusters[dbl.idx[,1]])]
    ls2 <- csz[as.character(clusters[dbl.idx[,2]])]
    dbl.idx$factor <- (dbl.idx$ls.ratio+ls1/(ls1+ls2))/2
    dbl.idx$factor[dbl.idx$factor>0.8] <- 0.8
    dbl.idx$factor[dbl.idx$factor<0.2] <- 0.2
    dbl.idx$ls <- (ls[dbl.idx[,1]]+ls[dbl.idx[,2]])
    x2 <- x[,dbl.idx[,1]]*dbl.idx$factor+x[,dbl.idx[,2]]*(1-dbl.idx$factor)
    x2 <- tryCatch(x2 %*% diag(dbl.idx$ls/Matrix::colSums(x2)),
                   error=function(e) t(t(x2)/Matrix::colSums(x2)))
    x1 <- cbind(x1,x2)
    rm(x2)
  }
  x <- x1
  rm(x1)
  if(halfSize>0){
    wAd <- sample.int(nrow(dbl.idx), size=ceiling(halfSize*nrow(dbl.idx)))
    if(length(wAd)>0)    x[,wAd] <- x[,wAd]/2
  }
  if(resamp>0){
    if(resamp!=halfSize) wAd <- sample.int(ncol(x), ceiling(resamp*ncol(x)))
    if(length(wAd)>0)
      x[,wAd] <- matrix(as.integer(rpois(nrow(x)*length(wAd),
                                         as.numeric(as.matrix(x[,wAd])))),
                        nrow=nrow(x))
  }else{
    x <- round(x)
  }
  x <- as(x,"CsparseMatrix")
  colnames(x) <- paste0( prefix, seq_len(ncol(x)) )
  x
}

cxds2 <- function (x, whichDbls = c(), ntop = 500, binThresh = NULL) {

  if (is.null(binThresh)) {
    if (is(x, "sparseMatrix")) {
      pNonZero <- length(x@x)/length(x)

    } else {
      len_x <- dim(x)[1] * dim(x)[2]
      num_non_zero <- sum(colSums(x > 0))
      pNonZero <- num_non_zero/len_x
    }

    if (pNonZero > 0.5) {
      pNonZero <- rowSums(x > 0)/ncol(x)
      x <- x[head(order(pNonZero), ntop), ]

      if (is(x, "sparseMatrix")) {
        binThresh <- max(1L, as.numeric(quantile(x@x, mean(pNonZero) * 0.5)))

      } else if (is(x, "IterableMatrix")) {
        binThresh <- max(1L, median(BPCells::colQuantiles(x, probs = 0.5)))
      } else {
        binThresh <- max(1L, median(x))
      }
    }
    else {
      binThresh <- 1L
    }
  }
  Bp <- x <- x >= binThresh
  ps <- rowMeans(x)
  if (nrow(x) > ntop) {
    hvg <- order(ps * (1 - ps), decreasing = TRUE)[seq_len(ntop)]
    Bp <- x <- x[hvg, ]
    ps <- ps[hvg]
  }
  if (length(whichDbls) > 0)
    Bp <- Bp[, -whichDbls]

  prb <- outer(ps, 1 - ps)
  prb <- prb + t(prb)

  Bp <- write_matrix_dir(Bp, dir = tempfile())
  Bpt <- 1 - Matrix::t(Bp)
  Bpt <- transpose_storage_order(Bpt)
  obs <- Bp %*% Bpt

  tstart <- Sys.time()
  obs <- as.matrix(obs)
  print(Sys.time()-tstart)

  obs <- obs + Matrix::t(obs)
  S <- suppressWarnings({
    stats::pbinom(as.matrix(obs) - 1, prob = prb, size = ncol(Bp),
                  lower.tail = FALSE, log.p = TRUE)
  })
  if (all(S == 0)) {
    return(rep(0L, ncol(x)))
  }
  if (any(w <- is.infinite(S))) {
    smin <- min(S[!is.infinite(S)])
    S[S < smin] <- smin
  }
  S <- as(S, 'IterableMatrix')
  Sx <- S %*% x

  pb <- txtProgressBar(min = 0, max = ncol(x), style = 3)

  # Convert subset to matrix before summing
  s <- apply_by_col(x, function(val, row, col) {
    setTxtProgressBar(pb, col)
    sum(as.matrix(Sx[row, col]))
  }) |> unlist()

  close(pb)

  s <- -s
  s <- s - min(s)
  s/max(s)
}

# Overwrite the original function in the package's namespace
# assignInNamespace("createDoublets", createDoublets, ns="scDblFinder")
