# Disable Bioconductor version checks for offline cluster nodes
options(BiocManager.check_repositories = FALSE)

options(warn = -1)

perf_wrapper <- function(perf_names, tst_scores, ytest){
  perf_df <- lapply(perf_names, function(perf_name){
    as.data.frame(t(sapply(tst_scores, function(preds){
      if(perf_name=="mxe"){preds <- pmax(pmin(preds, 1 - 1e-15), 1e-15)}  
      rocr_pred <- prediction(preds, as.numeric(as.character(ytest)))
      if(perf_name %in% c("acc", "f", "err")){
        curr_perf <- performance(rocr_pred, perf_name)  
        return(curr_perf@y.values[[1]][which.min(abs(curr_perf@x.values[[1]]-0.5))])
      }else{
        curr_perf <- performance(rocr_pred, perf_name)
        return(as.numeric(curr_perf@y.values))
      }
    })))
  })
  names(perf_df) <- perf_names
  perf_df <- do.call(rbind, perf_df)
  return(perf_df)
}

# Calculate confusion matrix elements and derived metrics
confusion_matrix_wrapper <- function(tst_scores, ytest, threshold = 0.5) {
  # Convert test labels to numeric
  ytest_numeric <- as.numeric(as.character(ytest))
  
  # Calculate confusion matrix elements for each method
  confusion_metrics <- sapply(tst_scores, function(preds) {
    # Convert predictions to binary using threshold
    pred_binary <- as.numeric(preds >= threshold)
    
    # Calculate confusion matrix elements
    tp <- sum(pred_binary == 1 & ytest_numeric == 1)
    fp <- sum(pred_binary == 1 & ytest_numeric == 0)
    tn <- sum(pred_binary == 0 & ytest_numeric == 0)
    fn <- sum(pred_binary == 0 & ytest_numeric == 1)
    
    # Calculate derived metrics
    # Matthews Correlation Coefficient
    mcc_denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc <- if(mcc_denom == 0) 0 else (tp * tn - fp * fn) / mcc_denom
    
    # Precision and Recall
    precision <- if(tp + fp == 0) 0 else tp / (tp + fp)
    recall <- if(tp + fn == 0) 0 else tp / (tp + fn)
    
    # Specificity
    specificity <- if(tn + fp == 0) 0 else tn / (tn + fp)
    
    # Balanced Accuracy
    sensitivity <- recall  # Same as recall
    balanced_acc <- (sensitivity + specificity) / 2
    
    return(c(
      tp = tp, fp = fp, tn = tn, fn = fn,
      mcc = mcc, precision = precision, recall = recall, 
      specificity = specificity, balanced_acc = balanced_acc
    ))
  })
  
  # Convert to data frame with proper structure
  # confusion_metrics should be a matrix with 9 rows (metrics) and n columns (methods)
  # After transpose: n rows (methods) and 9 columns (metrics)
  confusion_transposed <- t(confusion_metrics)
  confusion_df <- as.data.frame(confusion_transposed)
  
  # The rownames should be the metric names, but we transposed so they become column names
  expected_colnames <- c("tp", "fp", "tn", "fn", "mcc", "precision", "recall", "specificity", "balanced_acc")
  if(ncol(confusion_df) == length(expected_colnames)) {
    colnames(confusion_df) <- expected_colnames
  }
  
  # Now transpose back to get metrics as rows and methods as columns (like perf_df)
  confusion_df <- as.data.frame(t(confusion_df))
  
  # Remove any problematic row names
  if(!is.null(rownames(confusion_df))) {
    # Ensure row names are valid
    rownames(confusion_df) <- make.names(rownames(confusion_df), unique = TRUE)
  }
  
  return(confusion_df)
}



####  Helpers related to inverse gamma distribution
mv2ab <- function(m, v){
  a <- 2 + m^2/v
  b <- m * (a-1)
  return(list(alpha=a, beta=b))
}

ab2mv <- function(a, b){
  m <- b / (a-1)
  v <- b^2 / ((a-1)^2*(a-2))
  return(list(mean=m, var=v))
}


####  Take subset of dataset
reduceSize <- function(dat, y, N){
  reduced_ctrl <- sample(which(y==0)) 
  reduced_case <- sample(which(y==1))
  #identical(sort(reduced_ctrl), which(y==0))
  
  reduced_ctrl <- reduced_ctrl[1:N] 
  reduced_case <- reduced_case[1:N]
  
  reduced_indices <- c(reduced_ctrl, reduced_case)
  dat <- dat[, reduced_indices]
  y <- y[reduced_indices]
  return(list(dat=dat, y=y))
}


splitBatch <- function(condition, N_batch){
  # split samples into case / control groups
  case_ind <- which(condition==1)
  ctrl_ind <- which(condition==0)
  
  # split each condition group into N_batch batches
  batches_ind_case <- split(case_ind, sample(N_batch,length(case_ind),replace=TRUE))
  batches_ind_ctrl <- split(ctrl_ind, sample(N_batch,length(ctrl_ind),replace=TRUE))
  #print((sum(sapply(batches_ind_case,length))==length(case_ind)) & (sum(sapply(batches_ind_ctrl,length))==length(ctrl_ind)))
  
  # combine case / control samples in each batch
  batches_ind <- list()
  for(i in 1:N_batch){
    batches_ind[[i]] <- sort(c(batches_ind_case[[i]], batches_ind_ctrl[[i]]))
  }
  return(batches_ind)
}


subsetBatch <- function(condition, N_sample_size, N_batch){
  # split samples into case / control groups
  case_ind <- which(condition==1)
  ctrl_ind <- which(condition==0)
  
  # number of controls and cases to take
  N_ctrl <- N_sample_size / 2
  N_case <- N_sample_size / 2
  if(N_ctrl*N_batch > length(ctrl_ind) | N_case*N_batch > length(case_ind)){stop("Not enough samples to subset!")}
  
  # split each condition group into N_batch batches
  batches_ind_case <- split(case_ind, sample(N_batch,length(case_ind),replace=TRUE))
  while(any(sapply(batches_ind_case, length) < N_case)){batches_ind_case <- split(case_ind, sample(N_batch,length(case_ind),replace=TRUE))}
  batches_ind_case <- lapply(batches_ind_case, function(x){x[1:N_case]})
  
  batches_ind_ctrl <- split(ctrl_ind, sample(N_batch,length(ctrl_ind),replace=TRUE))
  while(any(sapply(batches_ind_ctrl, length) < N_ctrl)){batches_ind_ctrl <- split(ctrl_ind, sample(N_batch,length(ctrl_ind),replace=TRUE))}
  batches_ind_ctrl <- lapply(batches_ind_ctrl, function(x){x[1:N_ctrl]})
  
  # combine case / control samples in each batch
  batches_ind <- list()
  for(i in 1:N_batch){
    batches_ind[[i]] <- sort(c(batches_ind_case[[i]], batches_ind_ctrl[[i]]))
  }
  return(batches_ind)
}


####  Simulate batch effect based on ComBat assumption
simBatch <- function(dat, condition, batches_ind, batch, hyper_pars){
  n_batches <- sapply(batches_ind, length) # number of samples in each batch
  n_genes <- nrow(dat)
  
  ## Organize hyper batch parameters
  batch_par <- list()
  for(i in 1:length(n_batches)){
    batch_par[[i]] <- sapply(hyper_pars, function(item){item[i]}) # mean, sd of gaussian; alpha, beta for InvGamma
  }
    
  ## Simulate batch parameters from hyper-pars
  gamma <- delta2 <- list()
  for(i in 1:length(n_batches)){
    gamma[[i]] <- rnorm(n_genes, mean=batch_par[[i]]["hyper_mu"], sd=batch_par[[i]]["hyper_sd"])
            delta2[[i]] <- MCMCpack::rinvgamma(n_genes, shape=batch_par[[i]]["hyper_alpha"], scale=batch_par[[i]]["hyper_beta"])
  }
    
  ## Simulate batch effect
  # fit linear model to data with no batch parameters, calculate residual variance
  X <- model.matrix(~Condition, data=data.frame(Condition=condition))
  beta <- solve(t(X) %*% X) %*% t(X) %*% t(dat)
  resid <- dat - t(X %*% beta)
  #range(apply(resid,1,mean)); range(apply(resid,1,var))
  
  # spike-in batch variance: multiply by condition adjusted data with delta
  resid_varbatch <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat), dimnames=dimnames(dat))
  for(j in 1:length(n_batches)){
    curr_resid <- resid[, batches_ind[[j]]]
    spikein_var <- lapply(1:n_batches[j], function(col_ind){curr_resid[, col_ind] * sqrt(delta2[[j]])})
    resid_varbatch[, batches_ind[[j]]] <- do.call(cbind, spikein_var) 
  }
  #sapply(1:5, function(k){mean(apply(resid[, batches_ind[[k]]],1,var))})
  #sapply(1:5, function(k){mean(apply(resid_varbatch[, batches_ind[[k]]],1,var))})
  
  # construct mean batch parameter design matrix using gamma
  X_batch <- model.matrix(~-1+Batch, data=data.frame(Batch=factor(batch)))
  gamma_vec <- do.call(rbind, gamma)  #apply(gamma_vec,1,mean)
  
  # new data with added batch effect
  new_dat <- t(cbind(X, X_batch) %*% rbind(beta, gamma_vec)) + resid_varbatch
  if(!identical(rownames(new_dat), rownames(dat))){stop("BUG in simBatch function!")
  }else{colnames(new_dat) <- colnames(dat)}
  
  res <- list(new_dat=new_dat, batch_par=batch_par)
  return(res)
}


####  Gene-wise normalize datasets (z-score scaling)
normalizeData <- function(dat){
  dat_norm <- t(apply(dat, 1, scale, center=TRUE, scale=TRUE))
  dimnames(dat_norm) <- dimnames(dat)
  return(dat_norm)
}


#### Train pipeline: fit learner on training data only (no test data to prevent leakage)
trainPipe <- function(train_set, train_label, lfit=learner_fit){
  pred_res <- lfit(trn_set=train_set, y_trn=train_label)
  return(pred_res)
}


LogLossBinary <- function(actual, predicted, eps=1e-15) {
  if(class(actual)=="factor"){actual <- as.numeric(as.character(actual))}
  predicted = pmin(pmax(predicted, eps), 1-eps)
  - (sum(actual * log(predicted) + (1 - actual) * log(1 - predicted))) / length(actual)
}

AccuracyBinary <- function(actual, predicted) {
  if(class(actual)=="factor"){actual <- as.numeric(as.character(actual))}
  predicted_bi <- as.numeric(predicted >= 0.5)
  sum(predicted_bi==actual) / length(actual)
}


predSVM <- function(trn_set, y_trn){
  # Strip names to prevent metadata-alignment errors in tune()
  colnames(trn_set) <- NULL
  rownames(trn_set) <- NULL
  y_trn <- as.vector(y_trn)
  names(y_trn) <- NULL
  
  y_trn_factor <- as.factor(y_trn)
  tab <- table(y_trn_factor)
  min_samples <- min(tab)
  
  if (min_samples < 2) {
    stop(sprintf("SVM training failed: minimum 2 samples per class required (found classes: %s)", paste(names(tab), "=", as.numeric(tab), collapse=", ")))
  }
  
  # Ensure cross-validation folds do not exceed minimum samples in any class
  n_folds <- min(4, min_samples)
  tune_ctrl <- e1071::tune.control(sampling="cross", cross=n_folds)
  
  if (any(is.na(trn_set)) || any(is.infinite(trn_set))) {
    stop(sprintf("SVM training failed: trn_set contains NAs or Infs"))
  }
  
  if (any(is.na(y_trn))) {
    stop(sprintf("SVM training failed: y_trn contains NAs"))
  }

  # Hard cast to raw types to strip all hidden attributes/classes
  X_train <- matrix(as.numeric(as.matrix(t(trn_set))), nrow=ncol(trn_set))
  y_train <- as.factor(as.character(as.vector(y_trn)))
  
  cat(sprintf("DEBUG SVM: X_train=%d x %d, y_train=%d, levels=[%s], X_class=%s, y_class=%s\n", 
              nrow(X_train), ncol(X_train), length(y_train), 
              paste(levels(y_train), collapse=", "),
              class(X_train), class(y_train)))

  obj <- tryCatch({
    e1071::tune(e1071::svm, train.x=X_train, train.y=y_train,
                type="C-classification",
                kernel="linear",
                tunecontrol=tune_ctrl,
                ranges=list(cost=exp(seq(from=-10, to=10, by=2))))
  }, error = function(e) {
    cat(sprintf("WARNING: SVM tuning failed (%s). Falling back to default cost=1.\n", e$message))
    return(NULL)
  })
  
  if (is.null(obj) || is.null(obj$best.parameters)) {
    best_cost <- 1
  } else {
    best_cost <- obj$best.parameters[1, "cost"]
  }
  
  mod_svm <- e1071::svm(x=X_train, y=y_train,
                 type="C-classification", kernel="linear",
                 cost=best_cost, probability=TRUE)
  
  pred_train_svm_raw <- predict(mod_svm, X_train, probability=TRUE)
  probs <- attr(pred_train_svm_raw, "probabilities")
  
  if (is.null(probs)) {
    stop("SVM prediction failed: probabilities not found in output attributes.")
  }
  
  # Ensure we extract the probability for the correct class (labeled "1")
  if ("1" %in% colnames(probs)) {
    pred_train_svm <- probs[, "1", drop=TRUE]
  } else if (ncol(probs) >= 2) {
    # Fallback to the second column if "1" is missing
    cat("WARNING: SVM probability column '1' not found. Available columns:", paste(colnames(probs), collapse=", "), "\n")
    pred_train_svm <- probs[, 2, drop=TRUE]
  } else {
    stop(sprintf("SVM prediction failed: unexpected probability matrix dimensions %d x %d", nrow(probs), ncol(probs)))
  }
  
  res <- list(mod=mod_svm, pred_trn_prob=as.vector(pred_train_svm), pred_tst_prob=NULL,
              pred_trn_class=NULL, pred_tst_class=NULL)
  return(res)
}

predLasso_pp <- function(trn_set, y_trn, ...){
  mod <- glmnet::cv.glmnet(x=t(trn_set), y=as.numeric(as.character(y_trn)), family="binomial", ...)
  pred_trn_prob <- as.vector(predict(mod, newx=t(trn_set), s="lambda.1se", type="response"))
  pred_trn_class <- as.vector(predict(mod, newx=t(trn_set), s="lambda.1se", type="class"))
  return(list(mod=mod, pred_trn_prob=pred_trn_prob, pred_tst_prob=NULL,
              pred_trn_class=pred_trn_class, pred_tst_class=NULL))
}

predRF_pp <- function(trn_set, y_trn){
  trn_set <- as.matrix(trn_set); trn_transposed <- as.data.frame(t(trn_set))
  colnames(trn_transposed) <- paste0("V", 1:ncol(trn_transposed))
  data <- data.frame(y=as.factor(y_trn), trn_transposed)
  mod <- ranger::ranger(y ~ ., data = data, write.forest=TRUE)
  mod_prob <- ranger::ranger(y ~ ., data = data, write.forest=TRUE, probability=TRUE)
  
  pred_trn_prob <- predict(mod_prob, data = trn_transposed)$predictions[, "1"]
  pred_trn_class <- predict(mod, data = trn_transposed)$predictions
  
  return(list(mod = mod_prob, pred_trn_prob=pred_trn_prob, pred_tst_prob=NULL,
              pred_trn_class=pred_trn_class, pred_tst_class=NULL))
}

predRF_fs_pp <- function(trn_set, y_trn){
  data <- data.frame(y=as.factor(y_trn), t(trn_set))
  tout <- colttests(as.matrix(data[,-1]), data[,1], tstatOnly = TRUE)
  data <- data[, c("y",rownames(tout)[order(abs(tout[,1]), decreasing = T)[1:20]])]
  mod <- ranger(y ~ ., data = data, write.forest = TRUE, probability = T)
  pred_trn_prob <- predict(mod, data=data[,-1])$predictions[,2]
  return(list(mod=mod, pred_trn_prob=pred_trn_prob, pred_tst_prob=NULL,
              pred_trn_class=NULL, pred_tst_class=NULL))
}

predNnet_pp <- function(trn_set, y_trn){
  if(is.null(trn_set) || is.null(y_trn)) {
    stop("Training set or labels are NULL")
  }
  
  trn_set <- as.matrix(trn_set); trn_transposed <- as.data.frame(t(trn_set))
  colnames(trn_transposed) <- paste0("V", 1:ncol(trn_transposed))
  data <- data.frame(y=as.factor(y_trn), trn_transposed)

  
  n_samples <- nrow(data)
  network_size <- min(10, max(2, n_samples %/% 3))
  
  mod <- nnet::nnet(y ~ ., data = data, size = network_size, MaxNWts = 50000, 
              decay = 0.01, linout = FALSE, trace = FALSE, maxit = 200)
  
  pred_trn_prob <- as.vector(predict(mod, newdata = trn_transposed))
  pred_trn_class <- as.vector(predict(mod, newdata = trn_transposed, type="class"))
  
  return(list(mod=mod, pred_trn_prob=pred_trn_prob, pred_tst_prob=NULL,
              pred_trn_class=pred_trn_class, pred_tst_class=NULL))
}

predLogistic_pp <- function(trn_set, y_trn){
  trn_set <- as.matrix(trn_set); trn_transposed <- as.data.frame(t(trn_set))
  colnames(trn_transposed) <- paste0("V", 1:ncol(trn_transposed))
  data <- data.frame(y=as.factor(y_trn), trn_transposed)
  
  mod <- glm(y ~ ., data = data, family = binomial())
  
  pred_trn_prob <- as.vector(predict(mod, newdata = trn_transposed, type="response"))
  pred_trn_class <- as.vector(ifelse(pred_trn_prob >= 0.5, "1", "0"))
  
  return(list(mod=mod, pred_trn_prob=pred_trn_prob, pred_tst_prob=NULL,
              pred_trn_class=pred_trn_class, pred_tst_class=NULL))
}

predElnet_pp <- function(trn_set, y_trn){
  mod <- glmnet::cv.glmnet(x=t(trn_set), y=as.numeric(as.character(y_trn)), family="binomial", alpha=0.5)
  pred_trn_prob <- as.vector(predict(mod, newx=t(trn_set), s="lambda.1se", type="response"))
  pred_trn_class <- as.vector(predict(mod, newx=t(trn_set), s="lambda.1se", type="class"))
  return(list(mod=mod, pred_trn_prob=pred_trn_prob, pred_tst_prob=NULL,
              pred_trn_class=pred_trn_class, pred_tst_class=NULL))
}

predRDA_pp <- function(trn_set, y_trn){
  # Regularized Discriminant Analysis (RDA) using klaR package
  # Transpose data: R (features x samples) -> klaR expects (samples x features)
  X_train <- t(trn_set)
  y_train <- as.factor(y_trn)
  
  # Ensure data is in matrix format
  if (!is.matrix(X_train)) {
    X_train <- as.matrix(X_train)
  }
  
  # RDA hyperparameters as recommended in classifier_suggestion.md
  # gamma = 0.3: mostly shared covariance (0 = LDA, 1 = QDA)
  # lambda = 0.6: heavy diagonal shrinkage (0 = full covariance, 1 = diagonal)
  mod <- klaR::rda(
    x = X_train,
    grouping = y_train,
    gamma = 0.3,    # covariance shrinkage
    lambda = 0.6    # diagonal shrinkage
  )
  
  # Generate training predictions
  pred_obj <- predict(mod, X_train)
  
  # For binary classification, extract probability of positive class
  if (nlevels(y_train) == 2) {
    # Probability of positive class (second column, class "1")
    pred_trn_prob <- pred_obj$posterior[, "1"]
  } else {
    # Multiclass: pick max probability
    pred_trn_prob <- apply(pred_obj$posterior, 1, max)
  }
  
  # Class predictions
  pred_trn_class <- as.character(pred_obj$class)
  
  return(list(mod=mod, pred_trn_prob=pred_trn_prob, pred_tst_prob=NULL,
              pred_trn_class=pred_trn_class, pred_tst_class=NULL))
}

predKNN_pp <- function(trn_set, y_trn){
  n_samples <- ncol(trn_set)
  k_values <- c(3, 5, 7, 9, 11)
  k_values <- k_values[k_values < n_samples]
  
  if(length(k_values) == 0) {
    k_opt <- min(3, n_samples - 1)
  } else {
    cv_accuracy <- sapply(k_values, function(k) {
      if(n_samples < 10) {
        correct <- 0
        for(i in 1:n_samples) {
          train_idx <- setdiff(1:n_samples, i)
          trn_cv <- t(trn_set[, train_idx])
          tst_cv <- t(trn_set[, i, drop=FALSE])
          rownames(trn_cv) <- NULL
          rownames(tst_cv) <- NULL
          pred <- class::knn(train = trn_cv, test = tst_cv, cl = y_trn[train_idx], k = k)
          if(as.character(pred) == as.character(y_trn[i])) correct <- correct + 1
        }
        return(correct / n_samples)
      } else {
        folds <- cut(seq(1, n_samples), breaks = 5, labels = FALSE)
        correct <- 0
        total <- 0
        for(fold in 1:5) {
          test_idx <- which(folds == fold)
          train_idx <- which(folds != fold)
          if(length(test_idx) > 0 && length(train_idx) > 0) {
            trn_fold <- t(trn_set[, train_idx])
            tst_fold <- t(trn_set[, test_idx])
            rownames(trn_fold) <- NULL
            rownames(tst_fold) <- NULL
            pred <- class::knn(train = trn_fold, test = tst_fold, cl = y_trn[train_idx], k = k)
            correct <- correct + sum(as.character(pred) == as.character(y_trn[test_idx]))
            total <- total + length(test_idx)
          }
        }
        return(if(total > 0) correct / total else 0)
      }
    })
    k_opt <- k_values[which.max(cv_accuracy)]
  }
  
  pred_trn_class <- character(n_samples)
  for(i in 1:n_samples) {
    train_idx <- setdiff(1:n_samples, i)
    trn_subset <- t(trn_set[, train_idx])
    tst_subset <- t(trn_set[, i, drop=FALSE])
    rownames(trn_subset) <- NULL
    rownames(tst_subset) <- NULL
    pred_trn_class[i] <- as.character(class::knn(train = trn_subset, test = tst_subset, 
                                         cl = y_trn[train_idx], k = k_opt))
  }
  
  pred_trn_prob <- as.numeric(pred_trn_class == "1")
  mod <- list(k = k_opt, train_data = t(trn_set), train_labels = y_trn)
  
  return(list(mod=mod, pred_trn_prob=pred_trn_prob, pred_tst_prob=NULL,
              pred_trn_class=pred_trn_class, pred_tst_class=NULL))
}

predXGBoost_pp <- function(trn_set, y_trn){
  y_numeric <- as.numeric(as.factor(y_trn)) - 1
  train_matrix <- xgboost::xgb.DMatrix(data = t(trn_set), label = y_numeric)
  
  params <- list(
    objective = "binary:logistic",
    eval_metric = "logloss",
    max_depth = 3,
    eta = 0.1,
    subsample = 0.8,
    colsample_bytree = 0.8,
    verbose = 0
  )
  
  n_samples <- ncol(trn_set)
  nrounds <- if(n_samples < 50) 50 else if(n_samples < 200) 100 else 200
  
  if(n_samples >= 20) {
    cv_result <- xgboost::xgb.cv(
      params = params,
      data = train_matrix,
      nrounds = nrounds,
      nfold = min(5, n_samples %/% 4),
      early_stopping_rounds = 10,
      verbose = 0,
      showsd = FALSE
    )
    if(!is.null(cv_result$best_iteration) && length(cv_result$best_iteration) > 0 && !is.na(cv_result$best_iteration)) {
      best_nrounds <- cv_result$best_iteration
    } else {
      if(!is.null(cv_result$evaluation_log) && nrow(cv_result$evaluation_log) > 0) {
        best_nrounds <- which.min(cv_result$evaluation_log$test_logloss_mean)
      } else {
        best_nrounds <- min(50, nrounds)
      }
    }
  } else {
    best_nrounds <- min(30, nrounds)
  }
  
  mod <- xgboost::xgb.train(params = params, data = train_matrix, nrounds = best_nrounds, verbose = 0)
  
  pred_trn_prob <- predict(mod, train_matrix)
  pred_trn_class <- as.character(as.numeric(pred_trn_prob >= 0.5))
  
  return(list(mod=mod, pred_trn_prob=pred_trn_prob, pred_tst_prob=NULL,
              pred_trn_class=pred_trn_class, pred_tst_class=NULL))
}

predWrapper <- function(mod, tst_set, function_name){
  if(function_name=='logistic'){
    newdata <- as.data.frame(t(tst_set))
    colnames(newdata) <- paste0("V", 1:ncol(newdata))
    res <- as.vector(predict(mod, newdata = newdata, type="response"))
  }else if(function_name=='lasso'){
    res <- as.vector(predict(mod, newx=t(tst_set), s="lambda.1se", type="response"))
  }else if(function_name=='elnet' || function_name=='elasticnet'){
    res <- as.vector(predict(mod, newx=t(tst_set), s="lambda.1se", type="response"))
  }else if(function_name=='svm'){
    res <- predict(mod, t(tst_set), probability=TRUE)
    res <- attr(res, "probabilities")[,"1"]
  }else if(function_name=='rf'){
    newdata <- as.data.frame(t(tst_set))
    colnames(newdata) <- paste0("V", 1:ncol(newdata))
    res <- predict(mod, data = newdata)$predictions[, "1"]
  }else if(function_name=='nnet' || function_name=='nn'){
    newdata <- as.data.frame(t(tst_set))
    colnames(newdata) <- paste0("V", 1:ncol(newdata))
    res <- as.vector(predict(mod, newdata = newdata))
  }else if(function_name=='knn'){
    tst_for_knn <- t(tst_set)
    rownames(tst_for_knn) <- NULL
    res <- class::knn(train = mod$train_data, test = tst_for_knn, cl = mod$train_labels, k = mod$k)
    res <- as.numeric(as.character(res) == "1")
  }else if(function_name=='xgboost'){
    test_matrix <- xgboost::xgb.DMatrix(data = t(tst_set))
    res <- predict(mod, test_matrix)
  }else if(function_name=='rda'){
    # Regularized Discriminant Analysis prediction
    tst_transposed <- t(tst_set)
    pred_obj <- predict(mod, tst_transposed)
    # For binary classification, extract probability of positive class
    if (ncol(pred_obj$posterior) == 2) {
      res <- pred_obj$posterior[, "1"]
    } else {
      res <- apply(pred_obj$posterior, 1, max)
    }
    # Handle NA predictions by replacing with 0.5 (neutral probability)
    if (any(is.na(res))) {
      cat("[WARNING] RDA predictions contain", sum(is.na(res)), "NA values, replacing with 0.5\n")
      res[is.na(res)] <- 0.5
    }
  }
  
  # Final check for NA values in predictions
  if (any(is.na(res))) {
    cat("[WARNING] Predictions contain", sum(is.na(res)), "NA values after", function_name, "prediction, replacing with 0.5\n")
    res[is.na(res)] <- 0.5
  }
  
  return(res)
}

# Named dispatch table for classifiers — defined after all pred* functions so
# R can resolve the function objects when the list is evaluated.
CLASSIFIER_FUNCTIONS <- list(
  rda        = predRDA_pp,
  logistic   = predLogistic_pp,
  lasso      = predLasso_pp,
  elnet      = predElnet_pp,
  elasticnet = predElnet_pp,
  svm        = predSVM,
  rf         = predRF_pp,
  nnet       = predNnet_pp,
  nn         = predNnet_pp,
  knn        = predKNN_pp,
  xgboost    = predXGBoost_pp,
  rf_fs      = predRF_fs_pp
)

getPredFunctions <- function(learner_type) {
  fn <- CLASSIFIER_FUNCTIONS[[learner_type]]
  if (is.null(fn)) stop(sprintf("Classifier not supported: '%s'", learner_type))
  fn
}

quiet_melt <- function(...) {
  suppressMessages(melt(...))
}

load_results <- function(study_type = c("simulation", "real_4studies", "real_6studies"), 
                        study_name = NULL, metric = NULL, base_dir = "/scripts/evaluations/robustifying") {
  
  study_type <- match.arg(study_type)
  
  # Define possible result directories in order of preference
  possible_dirs <- switch(study_type,
    "simulation" = c(
      file.path(base_dir, "results")
    ),
    "real_4studies" = c(
      file.path(base_dir, "results", "real_4studies"),  # Unified structure (preferred)
      file.path(base_dir, "results"),                   # Copied files
      file.path(base_dir, "results_real_4studies")      # Original structure
    ),
    "real_6studies" = c(
      file.path(base_dir, "results", "real_6studies"),  # Unified structure (preferred)
      file.path(base_dir, "results"),                   # Copied files  
      file.path(base_dir, "results_real_6studies")      # Original structure
    )
  )
  
  # For real data, construct filename
  if (study_type != "simulation" && !is.null(study_name) && !is.null(metric)) {
    filename <- sprintf("test%s_%s.csv", study_name, metric)
    cat("DEBUG: Looking for specific file:", filename, "\n")
    
    # Try each directory until we find the file
    for (dir in possible_dirs) {
      filepath <- file.path(dir, filename)
      cat("DEBUG: Checking filepath:", filepath, "\n")
      cat("DEBUG: File exists:", file.exists(filepath), "\n")
      if (file.exists(filepath)) {
        cat("DEBUG: Found file at:", filepath, "\n")
        return(list(dir = dir, file = filepath))
      }
    }
    
    stop(sprintf("Could not find %s in any of the expected directories: %s", 
                filename, paste(possible_dirs, collapse = ", ")))
  }
  
  # For simulation data or when just getting directory
  for (dir in possible_dirs) {
    if (dir.exists(dir)) {
      csv_files <- list.files(dir, pattern = "\\.csv$")
      
      # Check if this directory contains the right type of files
      if (study_type == "simulation") {
        # Simulation files have patterns like "lasso_auc_batchN20_m0_v1.csv"
        sim_pattern_files <- grep("_batchN[0-9]+_m[0-9]+_v[0-9]+\\.csv$", csv_files)
        if (length(sim_pattern_files) >= 1) {
          dir <- sub("/$", "", dir)
          return(list(dir = dir, file = NULL))
        }
      } else {
        # Real data files have patterns like "testGSE37250_SA_auc.csv"
        real_pattern_files <- grep("^test[A-Za-z0-9_]+_[a-z]+\\.csv$", csv_files)
        if (length(real_pattern_files) >= 10) {
          dir <- sub("/$", "", dir)
          return(list(dir = dir, file = NULL))
        }
      }
    }
  }
  
  stop(sprintf("Could not find results directory for %s in: %s", 
              study_type, paste(possible_dirs, collapse = ", ")))
}
