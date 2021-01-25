compute_path_lengths <- function(.data){
  distances(.data) %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column(var='name2') %>% 
    pivot_longer(cols=-name2,values_to='PL') %>% 
    dplyr::filter(PL != Inf, name2 != name, PL != 1) %>%  # remove all infys and nodes with itself and take out adjacent nodes
    dplyr::filter(PL == min(PL))
}

compute_similarity <- function(.data, method, top_n=80, these_names){
  this_method <- rlang::ensym(method)
  sim_df <- similarity(.data, method=method) %>% 
    as.matrix() %>% 
    as.data.frame()
  
  to_character <- these_names
  # to_character <- as.character(1:length(V(.data)))
  colnames(sim_df) <- to_character
  rownames(sim_df) <- to_character
  
  sim_df <- sim_df %>% 
    rownames_to_column(var='name2') %>% 
    pivot_longer(cols=-name2) %>% 
    filter(name != name2) %>% 
    arrange(desc(value)) %>% 
    rename((!!this_method):= value) %>% 
    slice(1:top_n)  # keep top 20 since pairs are duplicated
  
  return(sim_df)
}

make_standardizer <- function(train_data, ...){
  require(caret)
  # fit preprocesser, assuming the only methods being passed will be c('center', 'scale') or 'range'
  kwargs <- list(...)
  std_fit <- preProcess(train_data, ...)  # pass the corresponding arguments
  this_method <- kwargs$method
  
  train_data <- c()  # avoid multiple copies of data being saved
  
  if(length(this_method)==1){  # if statement to save the transformation
    unfit_func <- function(new_data){
      num_data <- dplyr::select(new_data, where(is.numeric))  # TODO confirm this next line works
      sweep(sweep(num_data, MARGIN = 2, STATS=std_fit$ranges[1, ], FUN='+'),
            MARGIN=2, STATS=(std_fit$ranges[2, ] - std_fit$ranges[1, ]), FUN='*')  # inverse range func
    }
    things_to_save <- as.data.frame(std_fit$ranges)
    rownames(things_to_save) <- c('min', 'max')
  }else{
    unfit_func <- function(new_data){
      num_data <- dplyr::select(new_data, where(is.numeric))
      sweep(sweep(num_data, MARGIN=2, STATS=std_fit$std, FUN='*'), MARGIN=2, STATS=std_fit$mean, FUN='+')  # inverse z score
    }
    things_to_save <- cbind.data.frame(std_fit$mean, std_fit$std)
    colnames(things_to_save) <- c('mean', 'std')
  }
  
  standardize <- function(new_data){  # function to use the processor fit 
    predict(std_fit, newdata=new_data)
  }
  
  return(list(standardize=standardize, unfit_func=unfit_func, df=things_to_save))
}



make_predictor <- function(pred_func){
  function(model, new_data, type, response, thresh=.5){
    probs <- pred_func(model, new_data)  # return probs
    
    if(type=='prob'){
      return(probs)
    }
    
    preds <- if_else(probs>thresh, 1, 0)
    
    return(factor(preds, levels = c(0,1), labels=levels(new_data[, response])))
  }
}

make_roc <- function(model, this_data, pred_func, response, method, fit_thresh=T, ...){
  # kwargs are for the prediction func
  require(pROC)
  require(ggthemes)
  probs <- pred_func(model, this_data, type='prob')
  
  fit <- roc(this_data[, response], probs)  # fit roc curve
  
  metrics <- data.frame(sen = fit$sensitivities, spec = fit$specificities, thresh = fit$thresholds) %>% 
    mutate(metrics = sen+ spec)
  
  if(fit_thresh){
    opt_thresh <- metrics$thresh[which(metrics$metrics == max(metrics$metrics))]
  } else{
    opt_thresh=.5
  }
  this_auc <- fit$auc
  nice_preds <- pred_func(model, this_data, type='pred', response, thresh=opt_thresh)
  
  cm <- confusionMatrix(this_data[, response], nice_preds, ...)  # beware this may not be ordered correctly
  
  this_plot <- metrics %>% 
    pivot_longer(., cols=c('sen', 'spec')) %>%   # combine sensitivity and specificity
    ggplot(., aes(x=thresh)) + 
    geom_line(aes(y=value, color=name)) +
    geom_vline(xintercept = opt_thresh, linetype='dashed', color = 'limegreen') +
    geom_text(aes(x=opt_thresh, y=.5, label=paste('Optimal Threshold:', round(opt_thresh,3)))) +
    geom_text(aes(x=opt_thresh, y=.2, label=paste('AUC:', round(this_auc, 3)))) +
    theme_minimal() + 
    labs(x='Threshold', y='', color='Metric', title=paste('Method:', method)) +
    scale_color_colorblind(labels=c('Sensitivity', 'Specificity'))  # relabel metrics
  
  return(list(cm=cm, plot=this_plot))
}