med_iqr <- function(x, get_n = FALSE){
  
  med_x <- round(median(x, na.rm = T), 1)
  iqr_x <- round(quantile(x, probs = c(0.25, 0.75), na.rm = T), 1)
  
  res <- "-"
  
  if(!is.na(med_x)){
    res <- paste0(med_x, " (", iqr_x[1], " - ", iqr_x[2],")")
  }
  
  if(get_n){
    n <- sum(!is.na(x))  
    res <- paste0(res, " n=", n)
  }
  
  return(res)
  
}

mean_sd <- function(x, get_n = FALSE){
  
  mean_x <- round(mean(x, na.rm = T), 2)
  sd_x <- round(sd(x, na.rm = T), 2)
  
  res <- "-"
  
  if(!is.na(mean_x)){
    res <- paste0(mean_x, " ± ", sd_x)
  }
  
  if(get_n){
    n <- sum(!is.na(x))  
    res <- paste0(res, " n=", n)
  }
  return(res)
  
}

number_percent <- function(df, gr1, gr2){
  
  df %>% group_by(!!sym(gr1), !!sym(gr2)) %>%
    summarise(n = n()) %>%
    mutate(perc = round(n/sum(n) * 100, 2)) %>%
    mutate(comb =  paste0(n, " (", perc,"%)")) %>%
    tidyr::pivot_wider(id_cols = c(gr1, gr2), names_from = gr1, values_from = "comb") %>%
    mutate(across(.cols = c(-!!sym(gr2)),.fns = function(x) replace(x, which(is.na(x)), "0 (0.00%)")))
  
}

number_percent_one_factor <- function(df, gr1){
  
  df %>% group_by(!!sym(gr1)) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(perc = round(n/sum(n) * 100, 2)) %>%
    mutate(comb =  paste0(n, " (", perc,"%)")) %>%
    tidyr::pivot_wider(id_cols = c(gr1), names_from = gr1, values_from = "comb") #%>%
  # mutate(across(.cols = c(-!!sym(gr1)),.fns = function(x) replace(x, which(is.na(x)), "0 (0.00%)")))
  
}

mut_num_perc <- function(x, char_sel){
  
  abs <- sum(x == char_sel)
  rel <- round(abs/length(x)* 100, digits = 2)
  
  return(paste0(abs, " (", rel, "%)"))
}

wilcox_test_table <- function(df, group, feature, as_kable = FALSE, feature.levels = NULL){
  
  if(!is.null(feature.levels)){
    df <- df %>% mutate(!!sym(feature) := factor(!!sym(feature), levels = use.levels))
  }
  
  wt_table <- df %>% group_by(!!sym(group), .drop = FALSE) %>% dplyr::count(name = "Number of Subjects") %>% 
    tidyr::pivot_wider(names_from = !!sym(group), values_from = `Number of Subjects`) %>%
    replace(is.na(.), 0)
  
  uniq_groups <- df %>% pull(!!sym(group)) %>% unique()
  
  wt_res <- wilcox.test(x = df %>% dplyr::filter(!!sym(group) %in% uniq_groups[1]) %>% pull(!!sym(feature)),
                        y = df %>% dplyr::filter(!!sym(group) %in% uniq_groups[2]) %>% pull(!!sym(feature)))
  
  kable_table <- df %>% group_by(!!sym(group), .drop = FALSE) %>%
    summarise(val =  med_iqr(!!sym(feature))) %>%
    tidyr::pivot_wider(names_from = group, values_from = "val")
  
  kable_table <- cbind(kable_table, "P-value" = wt_res$p.value)

  if(as_kable){
    return(knitr::kable(kable_table, align = "c") %>%
             kableExtra::kable_styling(full_width = F, position = "left")%>%
             kableExtra::column_spec(1, bold = T) %>%
             kableExtra::collapse_rows(columns = 3:5, valign = "middle"))
  }else{
    
    return(kable_table)
  }
  
  
}

welch_test_table <- function(df, group, feature, as_kable = FALSE, feature.levels = NULL){
  
  if(!is.null(feature.levels)){
    df <- df %>% mutate(!!sym(feature) := factor(!!sym(feature), levels = use.levels))
  }
  
  wt_table <- df %>% group_by(!!sym(group), .drop = FALSE) %>% dplyr::count(name = "Number of Subjects") %>% 
    tidyr::pivot_wider(names_from = !!sym(group), values_from = `Number of Subjects`) %>%
    replace(is.na(.), 0)
  
  uniq_groups <- df %>% pull(!!sym(group)) %>% unique()
  
  wt_res <- t.test(x = df %>% dplyr::filter(!!sym(group) %in% uniq_groups[1]) %>% pull(!!sym(feature)),
                   y = df %>% dplyr::filter(!!sym(group) %in% uniq_groups[2]) %>% pull(!!sym(feature)),
                   paired = FALSE, var.equal = F, alternative= "two.sided")
  
  kable_table <- df %>% group_by(!!sym(group), .drop = FALSE) %>%
    summarise(val =  mean_sd(!!sym(feature))) %>%
    tidyr::pivot_wider(names_from = group, values_from = "val")
  
  kable_table <- cbind(kable_table, "P-value" = wt_res$p.value)
  
  # kable_table <- kable_table %>% mutate(
  #                               `P-value-short` = dplyr::case_when(
  #                                                 `P-value` >= 0.01 ~ as.character(round(`P-value`, 2)),
  #                                                 `P-value` < 0.0001 ~ "< 0.0001",
  #                                                 `P-value` < 0.001 ~ "< 0.001",
  #                                                 `P-value` < 0.01 ~ "< 0.01")) 
  
  if(as_kable){
    return(knitr::kable(kable_table, align = "c") %>%
             kableExtra::kable_styling(full_width = F, position = "left")%>%
             kableExtra::column_spec(1, bold = T) %>%
             kableExtra::collapse_rows(columns = 3:5, valign = "middle"))
  }else{
    
    return(kable_table)
  }
  
  
}

fisher_test_table <- function(df, gr_row, gr_col, as_kable = FALSE, use.levels = NULL){
  
  if(!is.null(use.levels)){
    df <- df %>% mutate(!!sym(gr_row) := factor(!!sym(gr_row), levels = use.levels))
  }
  
  ft_table <- df %>% group_by(!!sym(gr_row), !!sym(gr_col), .drop = FALSE) %>% dplyr::count(name = "Number of Subjects") %>% 
    tidyr::pivot_wider(names_from = !!sym(gr_row), values_from = `Number of Subjects`) %>%
    replace(is.na(.), 0)
  
  ft_res <- fisher.test(ft_table[,-1])
  
  kable_table <- df %>% group_by(!!sym(gr_row), !!sym(gr_col), .drop = FALSE) %>% dplyr::count(name = "n") %>% 
    group_by(!!sym(gr_col), .drop = FALSE) %>% mutate(val = paste0(n, " (", round(n/sum(n) * 100, 2), "%)" )) %>% dplyr::select(-n) %>%
    tidyr::pivot_wider(names_from = gr_col, values_from = "val") %>%
    replace(is.na(.), "0 (0.00%)")
  
  kable_table <- cbind(kable_table, "Odds-Ratio" = round(ft_res$estimate, 2), "P-value" = ft_res$p.value)
  
  kable_table <- kable_table %>% mutate(`P-value` = dplyr::case_when(
    `P-value` >= 0.01 ~ as.character(round(`P-value`, 2)),
    `P-value` < 0.0001 ~ "< 0.0001",
    `P-value` < 0.001 ~ "< 0.001",
    `P-value` < 0.01 ~ "< 0.01")) 
  
  if(as_kable){
    return(knitr::kable(kable_table, align = "c") %>%
             kableExtra::kable_styling(full_width = F, position = "left")%>%
             kableExtra::column_spec(1, bold = T) %>%
             kableExtra::collapse_rows(columns = 3:5, valign = "middle"))
  }else{
    
    return(kable_table)
  }
  
  
}

chisq_test_table <- function(df, gr_row, gr_col, as_kable = FALSE){
  
  chi_table <- df %>% group_by(!!sym(gr_row), !!sym(gr_col)) %>% dplyr::count(name = "Number of Subjects") %>%
    tidyr::pivot_wider(names_from = !!sym(gr_col), values_from = `Number of Subjects`) %>%
    replace(is.na(.), 0)
  
  chisq_res <- chisq.test(chi_table[,-1])
  
  kable_table <- df %>% group_by(!!sym(gr_row), !!sym(gr_col), .drop = FALSE) %>% dplyr::count(name = "n") %>% 
    group_by(!!sym(gr_col), .drop = FALSE) %>% mutate(val = paste0(n, " (", round(n/sum(n) * 100, 2), "%)" )) %>% dplyr::select(-n) %>%
    tidyr::pivot_wider(names_from = gr_col, values_from = "val") %>%
    replace(is.na(.), "0 (0.00%)")
  
  if(as_kable){
    return(knitr::kable(cbind(kable_table, "Odds-Ratio" =  "-", "P-value" = round(chisq_res$p.value, 2)), align = "c") %>%
             kableExtra::kable_styling(full_width = F, position = "left")%>%
             kableExtra::column_spec(1, bold = T) %>%
             kableExtra::collapse_rows(columns = dim(kable_table)[2]:(dim(kable_table)[2] + 2), valign = "middle"))
  }else{
    
    return(cbind(kable_table, "Odds-Ratio" =  "-", "P-value" = round(chisq_res$p.value, 2)))
  }
}

em.gmm.mu.sigma.pi <- function(x, K, epsilon = 1e-10, do.plot = F, iterations = 100){
  library("gtools")
  # Dimensionen der Daten festlegen
  N <- dim(x)[1]
  J <- dim(x)[2]
  
  # Initialisierung der mu, sigma, pi
  mu <- matrix(nr=K,nc=J,runif(K*J))  
  pi <- runif(K)
  pi <- pi/sum(pi)
  sigma <- matrix(nr=K,nc=J,rgamma(K*J,shape=0.5,rate=1))  
  
  log.ll <- c()
  g.nom <- c()
  denom <- c()
  
  gamma <- t(sapply(1:N,function(i){sample(x=c(rep(0,K-1),1),size=K)}))
  #gamma <- rdirichlet(K,rep(1,K))
  for(l in 1:iterations){
    
    # M-Schritt aus gammas Parameter bestimmen
    for(k in 1:K){
      N.k <- sum(gamma[,k])
      for(j in 1:J){
        mu[k,j] <- sum(gamma[,k] * x[,j])
        mu[k,j] <- mu[k,j]/N.k
        
        sigma[k,j] <- sum(gamma[,k] * (x[,j] - mu[k,j])^2)
        sigma[k,j] <- sigma[k,j]/N.k
      }
      pi[k] <- N.k/N
    }
    
    # E-Schritt: Berechnen der gammas
    g.nom <- matrix(nr=N,nc=K)
    for(k in 1:K){
      tmp <- matrix(nc=J,nr=N)
      for(j in 1:J){
        tmp[,j] <- dnorm(x[,j],mean=mu[k,j],sd=sqrt(sigma[k,j]))
      } 
      t <- apply(tmp,1,prod)
      g.nom[,k] <- t * pi[k]
    }
    denom <- rowSums(g.nom)
    
    for(k in 1:K){
      gamma[,k] <- g.nom[,k]/denom
    }
    
    log.ll[l] <- sum(log(denom))
    
    if(is.infinite(log.ll[l]) || is.na(log.ll[l])){
      gamma <- t(sapply(1:N,function(i){sample(x=c(rep(0,K-1),1),size=K)}))
      l <- 1
      k <- 1
      log.ll <- c()
      next;
    }
    
    if(length(log.ll)>1){
      
      if(is.na(abs(log.ll[l] - log.ll[l-1]))){
        gamma <- t(sapply(1:N,function(i){sample(x=c(rep(0,K-1),1),size=K)}))
        l <- 1
        k <- 1
        log.ll <- c()
        next;
      }
      
      if(abs(log.ll[l] - log.ll[l-1]) < epsilon){
        #print(paste("Break nach ",l," Iterationen!",sep=""))
        break
      }
    }
  }
  if(do.plot){
    plot(log.ll,type="l",ylab="logLikelihood",xlab="Iteration",main="Likelihood des EM bis Abbruch")
  }
  return(list(mu=mu,sigma=sigma,pi=pi,log.ll=log.ll, gamma = gamma));
}

max.em.gmm <- function(x, K=2, epsilon = 1e-2, do.plot = T, iterations = 100, runs = 100, plot.file = ""){
  
  em.list <- list()
  max.run <- c()
  max.ll <- -Inf
  for(i in 1:runs){
    em.list[[i]] <- em.gmm.mu.sigma.pi(x=x,K=K,epsilon=epsilon,do.plot=do.plot,iterations=iterations)
    if(max(em.list[[i]]$log.ll) > max.ll){
      max.run <- i
      max.ll <- max(em.list[[i]]$log.ll)
    }
  }
  
  return(em.list[[max.run]])
}

# calculation of f1 score
f1 <- function (data, lev = NULL, model = NULL) {
  
  precision <- posPredValue(data$pred, data$obs, positive = "pass")
  recall  <- sensitivity(data$pred, data$obs, postive = "pass")
  f1_val <- (2 * precision * recall) / (precision + recall)
  names(f1_val) <- c("F1")
  
  return(f1_val)
}

plot_ROC_5cv <- function(roc_merged_df){
  
  ggplot2::ggplot(data = roc_merged_df) +
    ggplot2::geom_line(ggplot2::aes(x = fpr, 
                                    y = tpr, 
                                    colour = CV_fold,
                                    linetype = CV_fold,
                                    size = CV_fold,
                                    alpha = CV_fold)) + 
    ggplot2::scale_alpha_manual(values = c(rep(0.8, 5), 1)) + 
    ggplot2::scale_size_manual(values = c(rep(1, 5), 1.5)) + 
    ggplot2::scale_linetype_manual(values = c(8, 6:3, 1)) +
    ggplot2::scale_color_manual(values = c("#DF8F44FF", "#00A1D5FF","#B24745FF", 
                                           "#79AF97FF", "#6A6599FF", "#000000")) +
    ggplot2::scale_x_continuous(expand = expansion(add = c(0.01, 0.01))) + 
    ggplot2::scale_y_continuous(expand = expansion(add = c(0.01, 0.01))) + 
    ggplot2::geom_abline(slope = 1, linetype = "dashed") +
    ggplot2::xlab("False positive rate") +
    ggplot2::ylab("True positive rate") + 
    ggplot2::theme_bw() + 
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = "bottom",
                   axis.text = ggplot2::element_text(colour = "black", size = 12),
                   strip.background = ggplot2::element_rect(fill = "white"),
                   axis.title = ggplot2::element_text(face = "bold", size = 12),
                   strip.text = ggplot2::element_text(face = "bold", size = 12)) +
    ggplot2::guides(col = ggplot2::guide_legend(nrow = 2))
}

# Training and validation of SVM model 
#-------------------------------------
# - 5-fold cross validation for model performance
# - 10-fold cross validation for parameter tuning within each cross validation

svm_cross_validation <- function(data_set, svm_method = "svmLinear", number_cv = 5){
  
  cv_index_mat <- matrix(nrow = number_cv, data = sample(seq_along(data_set[, 1])))
  
  roc_plt_list <- list()
  auc_list <- list()
  
  svmModels <- list()
  svmFitProbs_list <- list()
  svmFitClasses_list <- list()
  svmFitPredictions_list <- list()
  overall_perf <- c()
  
  approx_tpr <- list()
  mean_fpr <- seq(0, 1, length.out = 100)
  
  for (i in 1:nrow(cv_index_mat)) {
    
    train_df <- data_set[as.vector(cv_index_mat[-i, ]), ]
    test_df <- data_set[as.vector(cv_index_mat[i, ]), ]
    
    ctrl <- trainControl(method = "cv",
                         savePred = T,
                         classProb = T)
    
    svmFit <- c()
    
    if(svm_method == "svmLinear"){
      
      svmFit <- train(
        y ~ .,
        data = train_df,
        method = "svmLinear",
        trControl = ctrl,
        preProcess = c("center", "scale"),
        tuneGrid = expand.grid(C = seq(0, 2, length = 20)))
      
    }else if(svm_method == "svmRadial"){
      
      svmFit <- train(y ~ ., 
                      data = train_df,
                      method = "svmRadial", 
                      trControl = ctrl, 
                      preProcess = c("center","scale"), 
                      tuneLength = 10)
      
    }else if(svm_method == "svmPoly"){
      
      svmFit <- train(y ~ ., 
                      data = train_df, 
                      method = "svmPoly", 
                      trControl = ctrl, 
                      preProcess = c("center","scale"), 
                      tuneLength = 4)
      
    }else{
      
      stop("Please choose for svm_method between \"svmLinear\", \"svmRadial\", and \"svmPoly\"")
      
    }
    
    svmModels[[i]] <- svmFit
    
    svmFitProbs_list[[i]] <- predict(svmFit, test_df, type = "prob")
    svmFitClasses_list[[i]] <- predict(svmFit, test_df, type = "raw")
    svmFitPredictions_list[[i]] <- cbind(test_df$y, svmFitClasses_list[[i]], svmFitProbs_list[[i]])
    
    names(svmFitPredictions_list[[i]]) <- c("Observed", "Predicted", "probAC", "probnoAC")
    
    auc_list[[i]] <- svmFitPredictions_list[[i]] %>% roc_auc(truth = Observed, probAC)
    
    roc_plt_list[[i]] <- svmFitPredictions_list[[i]] %>%
      roc_curve(truth = Observed, probAC) %>%
      dplyr::mutate(CV_fold = paste0(i, ". run (auROC = ", format(auc_list[[i]]$.estimate, digits = 2), ")"),
                    fpr = 1 - specificity) %>%
      dplyr::rename(tpr = sensitivity) %>%
      dplyr::select(-specificity) %>% 
      dplyr::arrange(-.threshold)  %>% # rm threshold column for later merge with mean ROC coordinates
      dplyr::relocate(fpr, tpr, .threshold)
    
    overall_perf[i] <- caret::confusionMatrix(svmFitClasses_list[[i]], test_df$y)$overall["Accuracy"]
    
  }
  
  overall_auc <- do.call("rbind", auc_list) %>% 
    summarize(.estimate = mean(.estimate)) %>% 
    pull(.estimate)
  
  average_ROC <- genscore::avgrocs(lapply(roc_plt_list, function(tab) 
    tab %>% 
      dplyr::select(tpr, fpr) %>% 
      as.matrix()), 12, 8
  )
  
  # Add [0, 0] starting point to avg ROC
  average_ROC <- rbind(c(0,0), average_ROC)
  
  colnames(average_ROC) <- c("tpr", "fpr")
  roc_mean_df <- average_ROC %>% as_tibble() %>% 
    mutate(CV_fold = paste0("approx. mean ROC (mean auROC = ",
                            round(overall_auc, digits = 2),")"))
  
  rm(average_ROC)
  
  roc_single_df <- do.call("rbind", roc_plt_list) %>% dplyr::select(fpr, tpr, CV_fold)
  roc_merged_single_mean <- rbind(roc_single_df, roc_mean_df) %>% 
    dplyr::mutate(CV_fold = as.factor(CV_fold))
  
  svmFitPredictions_df <- data.frame(do.call("rbind", svmFitPredictions_list),
                                     run = rep(
                                       seq_along(svmFitPredictions_list),
                                       each = nrow(svmFitPredictions_list[[1]])
                                     ))
  
  perf_plot <- plot_ROC_5cv(roc_merged_df = roc_merged_single_mean)
  
  return(
    list(roc_plt_list= roc_plt_list,
         auc_list = auc_list,
         svmFitPredictions_df = svmFitPredictions_df,
         svmFitProbs_list = svmFitProbs_list,
         svmFitClasses_list = svmFitClasses_list,
         svmFitPredictions_list = svmFitPredictions_list,
         overall_perf = overall_perf,
         approx_tpr = approx_tpr,
         mean_fpr = mean_fpr,
         roc_merged_single_mean = roc_merged_single_mean,
         perf_plot = perf_plot
    )
  )
}