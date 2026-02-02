### PLS-DA model

# x: matrix or data.frame or list (if multi-block).
# y: matrix or data.frame
# ncomp: If NULL, estima las maximas. Si NULL, automatico.
# train: Percentage to take as training set

plsda = function(x, y,
                 ncomp = NULL, cvFolds = 5, rep = 10, perm = 20,
                 scaling = c("none", "center", "standard", "softBlock", "hardBlock")[3],
                 scalingY = c("none", "center", "standard")[2],
                 algo = c("nipals","mbpls")[1],
                 train = 1, alpha = 0.05, parallel = FALSE){

  if(!inherits(x, "list")){
    if(!all(sapply(x, is.numeric))) {
      warning('Categorical variables automatically converted to dummy variables and default filters for NAs and CV applied.\n If the user does not want to use any default filtering, consider using Preparing function before executing pca to convert categorical variables into dummies.')
      P = Preparing(x, includeFactors = T, CVfilter = 0.01, excludeNA = 0.2)
      x = P$x
    } else {
      P = Preparing(x, CVfilter = 0.00001, excludeNA = 1) #Evitar nearZeroVariance variables pero no aplicar ningun filtro extra
      x = P$x
    }

    desSummary = data.frame(c("Samples","Variables","Excluded Near-Zero Variables","Missing Values"),
                            c(nrow(x),ncol(x),length(which(P$removed$Variables$Problem == 'Low CV')),
                              paste0(length(which(is.na(x))), " (",round((length(which(is.na(x))) / (nrow(x) * ncol(x))) * 100, 2), "%)")), stringsAsFactors = FALSE)
    colnames(desSummary) = NULL

  } else{
    b_names = names(x)
    P = lapply(x, function(y){
      if(!all(sapply(y, is.numeric))) {
        warning('Categorical variables automatically converted to dummy variables and default filters for NAs and CV applied.\n If the user does not want to use any default filtering, consider using Preparing function before executing pca to convert categorical variables into dummies.')
        Preparing(y, includeFactors = T, CVfilter = 0.01, excludeNA = 0.2)
      } else {
        Preparing(y, CVfilter = 0.00001, excludeNA = 1) #Evitar nearZeroVariance variables pero no aplicar ningun filtro extra
      }
    })
    x = setNames(lapply(b_names, function(y) P[[y]]$x), b_names)

    desSummary = setNames(lapply(b_names, function(y){
      des = data.frame(c("Samples","Variables","Excluded Near-Zero Variables","Missing Values"),
                       c(nrow(x[[y]]),ncol(x[[y]]),length(which(P[[y]]$removed$Variables$Problem == 'Low CV')),
                         paste0(length(which(is.na(x[[y]]))), " (",round((length(which(is.na(x[[y]]))) / (nrow(x[[y]]) * ncol(x[[y]]))) * 100, 2), "%)")), stringsAsFactors = FALSE)
      colnames(des) = NULL
      return(des)
    }),b_names)

  }

  if(!inherits(x, "list") && scaling %in% c('softBlock','hardBlock')) blocks = rep(1, times = ncol(x)) else blocks = NULL

  if (class(x)[1]=='list') {
    if (!scaling %in% c("softBlock", "hardBlock")) warning('When blocks are considered we recommend using either softBlock or hardBlock scaling')
    if (length(unique(sapply(x, nrow)))!=1) return(stop('All blocks must have the same observations'))
    blocks = rep(seq_along(x), times = sapply(x, ncol) )
    x = do.call(cbind, x)
  }

  X = x
  Y = y

  #Convert into dummy categorical response
  if(!all(sapply(y, is.numeric))) {
    Py = suppressWarnings(Preparing(y, includeFactors = T, CVfilter = 0.01, excludeNA = 0.2))
    y = Py$x
  } else {
    warning('We recommed using PLS instead of PLS-DA model')
    Py = suppressWarnings(Preparing(y, CVfilter = 0.00001, excludeNA = 1)) #Evitar nearZeroVariance variables pero no aplicar ningun filtro extra
    y = Py$x
  }

  ## Divide data into train and test taking into account group ratios

  if (train < 1) {
    ntrain = round(nrow(x)*train, 0)
    group_indices = split(1:nrow(x),Y)
    train_indices = lapply(group_indices, function(i){
      ntrain_group = round(length(i)*train,0)
      sample(i,ntrain_group)
    } )
    rowtrain = unlist(train_indices)
    xTest = x[-rowtrain, ,drop=FALSE]
    x = x[rowtrain,, drop=FALSE]
    yTest = y[-rowtrain,, drop=FALSE ]
    y = y[rowtrain,,drop=FALSE]
    YTest = Y[-rowtrain,,drop=FALSE]
    Y = Y[rowtrain,,drop=FALSE]
  } else {
    xTest = yTest = NULL
    YTest = NULL
  }
  X2 = x
  Y2 = y

  #Calcular los parametros optimos

  if(is.null(ncomp)){

    res_cross = plsda_cross(X = x, Y = y, Y2=Y, scaling = scaling, blocks = blocks, scalingY = scalingY, ncomp = ncomp, folds = cvFolds, rep = rep, parallel = parallel, algo = algo)

    r2_array = simplify2array(lapply(res_cross, function(res) res$R2))
    mean_r2_per_component = apply(r2_array, 1, mean)

    q2_array = simplify2array(lapply(res_cross, function(res) res$Q2))
    mean_q2_per_component = apply(q2_array, 1, mean)

    rmse_array = simplify2array(lapply(res_cross, function(res) res$RMSE))
    mean_rmse_per_component = apply(rmse_array, 1, mean)

    ber_array = simplify2array(lapply(res_cross, function(res) res$BER))
    mean_ber_per_component = apply(ber_array, 1, mean)

    # Improvement of 0.01 in R2 and Q2 positive increment
    r2_diff = diff(mean_r2_per_component)
    q2_diff = diff(mean_q2_per_component)

    ncomp = kneedle::kneedle(1:nrow(r2_array), mean_ber_per_component, decreasing = T)[1] #BER    La libreria esta basada en el paper Finding a Kneedle in a Haystack: Detecting Knee Points in System Behavior

  } else{
    res_cross = pls_cross_plsda(X = x, Y = y, Y2=Y, scaling = scaling, blocks = blocks, scalingY = scalingY, ncomp = ncomp, folds = cvFolds, rep = rep, parallel = parallel, algo = algo)

    q2_array = simplify2array(lapply(res_cross, function(res) res$Q2))
    mean_q2_per_component = apply(q2_array, 1, mean)

    rmse_array = simplify2array(lapply(res_cross, function(res) res$RMSE))
    mean_rmse_per_component = apply(rmse_array, 1, mean)

  }

  ## Scaling X

  if (scaling != 'none') {
    escalado = Scaling(x, scaling= scaling, blocks = blocks)
    x = escalado$x
    centrado = escalado$centering
    escalado = escalado$scaling
  } else{
    escalado = rep(1, times = ncol(x))
    centrado = rep(0, times = ncol(x))
  }

  ## Scaling Y

  if (scalingY != 'none') {
    escaladoY = Scaling(y, scaling= scalingY, blocks = NULL)
    y = escaladoY$x
    centradoY = escaladoY$centering
    escaladoY = escaladoY$scaling
  } else{
    escaladoY = rep(1, times = ncol(y))
    centradoY = rep(0, times = ncol(y))
  }

  ## Modelo PLS

  mypls = nipals_pls(x,y,ncomp = ncomp)

  coef_cv = array(unlist(lapply(res_cross, `[[`, "coef_cv")),  dim = c(dim(res_cross[[1]]$coef_cv)[1:3], rep * dim(res_cross[[1]]$coef_cv)[4]))
  Jack = p.jack(mypls$coefficients, coef_cv, ncomp, alpha)

  resum_tab = data.frame(
    variable   = rownames(mypls$coefficients),
    reshape2::melt(mypls$coefficients, varnames = c("variable", "response"), value.name = "coefficient")[, -1],
    pValJK = as.vector(Jack$pval),
    LCI_JK = as.vector(Jack$LCI_coef),
    UCI_JK = as.vector(Jack$UCI_coef)
  )
  resum_tab = split(resum_tab, resum_tab$response)

  resum_coef = lapply(resum_tab, function(df) {
    out = df[, c("variable", "coefficient", "pValJK", "LCI_JK", "UCI_JK")]
    rownames(out) = out$variable
    out$variable = NULL
    return(out)
  })

  ypred = as.matrix(x) %*% mypls$coefficients
  SCEY = sum(colSums(ypred**2))

  SCEYa =  sapply(1:ncomp, function(j) sum(tcrossprod(mypls$scores[, j], mypls$loadingsY[, j])^2))
  vip = sqrt(ncol(x) * rowSums(sweep(mypls$weights^2, 2,SCEYa, "*")) / SCEY)

  R2 = R2Y = f1score = MCC = Etot = BER = numeric(ncomp)

  if (is.null(xTest)) {
    for (i in 1:ncomp) {
      xpred = tcrossprod(mypls$scores[,1:i,drop=FALSE], mypls$loadings[,1:i,drop=FALSE])
      ypred = tcrossprod(mypls$scores[,1:i,drop=FALSE],mypls$loadingsY[,1:i,drop=FALSE])
      #R2X
      SCR = sum(colSums((xpred-x)**2))
      SCT = sum(colSums(x**2))
      R2[i] = 1- (SCR/SCT)
      #R2Y
      SCRY = sum(colSums((ypred-y)**2))
      SCTY = sum(colSums(y**2))
      R2Y[i] = 1- (SCRY/SCTY)
      #F1-score
      real_class = apply(y, 1, which.max)
      pred_class = apply(ypred, 1, which.max)

      conf = table(factor(real_class, levels=1:length(table(real_class))), factor(pred_class, levels=1:length(table(real_class))))
      #Define  TP, FN,FP,TN per class
      TP = diag(conf)
      FN = rowSums(conf) - TP
      FP = colSums(conf) - TP
      TN = sum(conf) - TP - FN - FP

      #F1
      f1score[i] = f1s(TP,FP,FN, 0.000001)
      MCC[i] = MCCs(TP,FP,FN,TN, 0.000001)

      #Error
      Etot[i] = mean(real_class != pred_class)
      BER[i] = BER(TP, FP, FN, TN, epsilon = 0.00001)

    }

    r2_diff = diff(R2)
    r2y_diff = diff(R2Y)

    q2_diff = diff(mean_q2_per_component[1:ncomp])

    resum = data.frame("comp" = factor(1:ncomp),
                       "R2X" = c(R2[1], r2_diff),
                       "cumR2X" = R2,
                       "R2Y" = c(R2Y[1], r2y_diff),
                       "cumR2Y" = R2Y,
                       "Q2" = c(mean_q2_per_component[1], q2_diff),
                       "cumQ2" = mean_q2_per_component[1:ncomp],
                       "F1score" = f1score,
                       "MCC" = MCC,
                       "Etot" = Etot,
                       "BER" = BER)

  } else{

    xTest = sweep(xTest, 2, centrado, FUN = "-")
    xTest = sweep(xTest, 2, escalado, FUN = "/")

    yTest = sweep(yTest, 2, centradoY, FUN = "-")
    yTest = sweep(yTest, 2, escaladoY, FUN = "/")

    R2test = R2Ytest = f1scoretest = MCCtest = BERtest = Etottest = numeric(ncomp)

    for (i in 1:ncomp) {

      weightsStar =  try(suppressWarnings(mypls$weights[,1:i,drop=FALSE]%*%solve(crossprod(mypls$loadings[,1:i,drop=FALSE], mypls$weights[,1:i,drop=FALSE]))),silent = TRUE)
      if(inherits(weightsStar,'try-error')) weightsStar = mypls$weights[,1:i,drop=FALSE]%*%corpcor::pseudoinverse(crossprod(mypls$loadings[,1:i,drop=FALSE], mypls$weights[,1:i,drop=FALSE]))
      coefficients = tcrossprod(weightsStar, mypls$loadingsY[,1:i,drop=FALSE])
      Ttest = as.matrix(xTest) %*% weightsStar

      xpred = tcrossprod(mypls$scores[,1:i,drop=FALSE], mypls$loadings[,1:i,drop=FALSE])
      ypred = tcrossprod(mypls$scores[,1:i,drop=FALSE],mypls$loadingsY[,1:i,drop=FALSE])
      #R2X
      SCR = sum(colSums((xpred-x)**2))
      SCT = sum(colSums(x**2))
      R2[i] = 1- (SCR/SCT)
      #R2Y
      SCRY = sum(colSums((ypred-y)**2))
      SCTY = sum(colSums(y**2))
      R2Y[i] = 1- (SCRY/SCTY)

      #F1score
      real_class = apply(y, 1, which.max)
      pred_class = apply(ypred, 1, which.max)

      conf = table(factor(real_class, levels=1:length(table(real_class))), factor(pred_class, levels=1:length(table(real_class))))
      #Define  TP, FN,FP,TN per class
      TP = diag(conf)
      FN = rowSums(conf) - TP
      FP = colSums(conf) - TP
      TN = sum(conf) - TP - FN - FP

      #F1
      f1score[i] = f1s(TP,FP,FN, 0.000001)
      MCC[i] = MCCs(TP,FP,FN,TN, 0.000001)

      #Error
      Etot[i] = mean(real_class != pred_class)
      BER[i] = BER(TP, FP, FN, TN, epsilon = 0.00001)


      xpred = tcrossprod(Ttest[,1:i,drop=FALSE], mypls$loadings[,1:i,drop=FALSE])
      ypred = as.matrix(xTest) %*% coefficients
      #R2Xtest
      SCR = sum(colSums((xpred-xTest)**2))
      SCT = sum(colSums(xTest**2))
      R2test[i] = 1- (SCR/SCT)
      #R2Ytest
      SCRY = sum(colSums((ypred-yTest)**2))
      SCTY = sum(colSums(yTest**2))
      R2Ytest[i] = 1- (SCRY/SCTY)
      #F1scoretest
      real_class = apply(yTest, 1, which.max)
      pred_class = apply(ypred, 1, which.max)
      conf = table(factor(real_class, levels=1:length(table(real_class))), factor(pred_class, levels=1:length(table(real_class))))
      #Define  TP, FN,FP,TN per class
      TP = diag(conf)
      FN = rowSums(conf) - TP
      FP = colSums(conf) - TP
      TN = sum(conf) - TP - FN - FP
      #F1scoretest
      f1scoretest[i] = f1s(TP,FP,FN, 0.000001)
      MCCtest[i] = MCCs(TP,FP,FN,TN, 0.000001)
      #Errortest
      Etottest[i] = mean(real_class != pred_class)
      BERtest[i] = BER(TP, FP, FN, TN, epsilon = 0.00001)

    }
    r2test_diff = diff(R2test)
    r2ytest_diff = diff(R2Ytest)

    r2_diff = diff(R2)
    r2y_diff = diff(R2Y)

    q2_diff = diff(mean_q2_per_component[1:ncomp])

    resum = data.frame("comp" = factor(1:ncomp),
                       "R2X" = c(R2[1], r2_diff),
                       "cumR2X" = R2,
                       "R2Y" = c(R2Y[1], r2y_diff),
                       "cumR2Y" = R2Y,
                       "R2Xtest" = c(R2test[1], r2test_diff),
                       "cumR2Xtest" = R2test,
                       "R2Ytest" = c(R2Ytest[1], r2ytest_diff),
                       "cumR2Ytest" = R2Ytest,
                       "Q2" = c(mean_q2_per_component[1], q2_diff),
                       "cumQ2" = mean_q2_per_component[1:ncomp],
                       "RMSE" = mean_rmse_per_component[1:ncomp],
                       "F1score" = f1score,
                       "F1scoretest" = f1scoretest,
                       "MCC" = MCC,
                       "MCCtest" = MCCtest,
                       "Etot" = Etot,
                       "Etottest" = Etottest,
                       "BER" = BER,
                       "BERtest" = BERtest)

  }

  res_val = plsValidation(X2,Y2,scaling, scalingY,ncomp,cvFolds,perm,algo, parallel)

  return(list("scoresX" = mypls$scores,
              "loadingsX" = mypls$loadings,
              "scoresY" = mypls$scoresY,
              "loadingsY" = mypls$loadingsY,
              "weights" = mypls$weights,
              "weightStar" = mypls$weightsStar,
              "coefficients" = mypls$coefficients,
              "vip" = vip,
              "X" = x, # escalada y sin variables con baja variabilidad
              "Y" = y, # escalada
              "scaling" = list("center" = centrado, "scale" = escalado, #para que cuadre con los scorePLot de cuando hicimos PCA
                               "centerY" = centradoY, "scaleY" = escaladoY), # completar con block??
              "summary" = resum, ## completar si se desea dar + info de lo que se ha hecho
              "coefficients_summary" = resum_coef,
              "ncomp" = ncomp,
              "test" = list("yTest" = yTest, "xTest" = xTest),
              "validation" = as.data.frame(res_val),
              "cv_results" = res_cross,
              "PreproSummary" = desSummary,
              "input" = list("scalingType" = scaling,"scalingTypeY" = scalingY, "X" = X, "Y" = Y,
                             'algo' = algo, "alpha" = alpha, "perm" = perm, "blocks" = blocks)))

}






### PLS predictions

# x: PLS object

plsPredict = function(x, new = NULL, plot = TRUE) {
  ## proyectar nuevas observaciones antes de predecir

  x$explVar = data.frame("comp" = factor(1:ncomp),
                         "percVar" = round(100*x$summary$R2X,4),
                         "cumPercVar" = round(100*x$summary$cumR2X,4))

  scorePlot(x, newObs = new)

  new = as.data.frame(new)
  new = new[,colnames(x$X)]

  new = sweep(new, 2, x$scaling$center, FUN = "-")
  new = sweep(new, 2, x$scaling$scale, FUN = "/")

  prediction = as.matrix(new) %*% x$coefficients

  return(prediction)

}




### PLS error metrics

# x: PLS object

plsError = function(x) {

  ## Metrics:   RMSE,
  # plots???

  # Train


  # Repeated CV --> especialmente para cuando no haya test
  # Calcular error medio por fold, por run y global


  # Test


}



### PLS plots

# x: Object return by pls() function

plsdaPlot = function(x,
                   type = c("ncomp","R2vsQ2", "loadings", "scoresX",  "loadingsX", "scoresY", "loadingsY",
                            "weights", "linearity", "overfitting", "R2", "corr", "biplot", "coef"),
                   comp = NULL,
                   col = c('main', 'complete', 'cblindfriendly','sunshine','hot','warm','grass','oficial')[1],
                   colBy = NULL,
                   shape = NULL,
                   shapeBy = NULL,
                   ellipses = NULL,  # only for scores
                   selVars = NULL,
                   labels = NULL,
                   labelTop = NULL, #percentage expressed in 0-1 range
                   repel = TRUE, #not implemented for scores
                   newObs = NULL #data.frame
) {

  if(!(type %in% c("ncomp","R2vsQ2", "loadings", "scoresX",  "loadingsX", "scoresY", "loadingsY", "weights", "linearity", "overfitting", "R2", "corr", "biplot", "coef"))) return(stop('Please use one of: R2vsQ2, loadings, scoresX, loadingsX, scoresY, loadingsY, weights, linearity, overffiting, R2, correl, biplot.'))

  eigenvalues = apply(x$scoresX, 2, function(t) crossprod(t)/(nrow(x$X)-1) )

  x$explVar = data.frame("comp" = factor(1:x$ncomp),
                         "percVar" = round(100*x$summary$R2X,4),
                         "cumPercVar" = round(100*x$summary$cumR2X,4),
                         "eigenVal" = eigenvalues)

  ## ncomp

  if (type == "ncomp") {

    ## Plot ropls
    mean_per_comp = function(cv_results, metric, comp = NULL) {
      arr = simplify2array(lapply(cv_results, function(res) res[[metric]]))
      m = apply(arr, 1, mean)
      if (!is.null(comp)) m = m[1:comp]
      return(m)
    }

    mean_r2_per_component = mean_per_comp(x$cv_results, 'R2', comp)
    mean_q2_per_component = mean_per_comp(x$cv_results, 'Q2', comp)

    if(is.null(comp)) comp = length(mean_r2_per_component)

    df = data.frame(Component = rep(1:comp, 2),Metric = rep(c("R2", "Q2"), each = comp),Value = c(mean_r2_per_component, mean_q2_per_component))

    if ((length(col) == 1 && !(col %in% c('main', 'complete', 'cblindfriendly', 'sunshine','hot','warm','grass','oficial'))) || length(col)>1){
      if (length(col)!= 3) return(stop('Either provide 3 colors or use one of the default color palettes'))
      color_palette = col
      custom_colors = setNames(col[-2], unique(df$Metric))
    } else{
      color_palette = colorbiostat(3, palette = col)
      custom_colors = setNames(color_palette[-2], unique(df$Metric))
    }

    ncomp_ropls = {
     idx = which(!(diff(mean_r2_per_component) >= 0.01 & diff(mean_q2_per_component) > 0))
    if (length(idx) == 0) length(mean_r2_per_component) else idx[1]}

    ggp1 = ggplot(df, aes(x = Component, y = Value, color = Metric)) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      scale_color_manual(values = custom_colors) +
      labs(
        title = "R² and Q²",
        x = "Component",
        y = "Metric"
      ) +
      theme_minimal(base_size = 10 ) +
      scale_x_continuous(breaks = df$Component) +
      geom_vline(xintercept = ncomp_ropls, alpha = 0.5)

    plot_metric = function(values, title, ylab, decreasing = FALSE, base_size = 10, color = 'skyblue'){
      comp = length(values)
      df = data.frame(Component = 1:comp, Value = values)
      ncomp = kneedle::kneedle(1:comp, values, decreasing = decreasing)[1]
      ggp = ggplot(df, aes(Component, Value)) +
        geom_line(linewidth = 1, color = color) +
        geom_point(size = 3, color = color) +
        geom_vline(xintercept = ncomp, alpha = 0.5) +
        scale_x_continuous(breaks = 1:comp) +
        labs(title = title, x = "Component", y = ylab) +
        theme_minimal(base_size = base_size)
      return(ggp)
    }

    ggp2 = plot_metric(mean_per_comp(x$cv_results, "F1score", comp),title = "F1-score",ylab = "F1-score",color = color_palette[2])
    ggp3 = plot_metric(mean_per_comp(x$cv_results, "MCC", comp),title = "Matthews Correlation Coefficient",ylab = "MCC",color = color_palette[2] )
    ggp4 = plot_metric(mean_per_comp(x$cv_results, "Etot", comp),title = "Total Error",ylab = "Etot", decreasing = T,color = color_palette[2] )
    ggp5 = plot_metric(mean_per_comp(x$cv_results, "BER", comp),title = "Balanced Error Rate",ylab = "BER", decreasing = T,color = color_palette[2] )
    ggp6 = plot_metric(mean_per_comp(x$cv_results, "RMSE", comp),title = "Root Mean Squared Error",ylab = "RMSE", decreasing = T,color = color_palette[2] )

    ggp = patchwork::wrap_plots(ggp1, ggp6, ggp2, ggp3, ggp4, ggp5, nrow = 2) +
      patchwork::plot_annotation(
        title = "Results per Component",
        theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5))
      )

  }

  ## R2vsQ2

  if (type == "R2vsQ2") {

    r2_array = simplify2array(lapply(x$cv_results, function(res) res$R2))
    mean_r2_per_component = apply(r2_array, 1, mean)
    q2_array = simplify2array(lapply(x$cv_results, function(res) res$Q2))
    mean_q2_per_component = apply(q2_array, 1, mean)

    if(!is.null(comp)){
      mean_r2_per_component = mean_r2_per_component[1:comp]; mean_q2_per_component = mean_q2_per_component[1:comp]
    }  else {
      comp = length(mean_r2_per_component)
    }

    df = data.frame(Component = rep(1:comp, 2),Metric = rep(c("R2", "Q2"), each = comp),Value = c(mean_r2_per_component, mean_q2_per_component))

    #if(col=='main') col = "skyblue3"

    if ((length(col) == 1 && !(col %in% c('main', 'complete', 'cblindfriendly', 'sunshine','hot','warm','grass','oficial'))) || length(col)>1){
      if (length(col)!= 2) return(stop('Either provide 2 colors or use one of the default color palettes'))
      custom_colors = setNames(col, unique(df$Metric))
    } else{
      color_palette = colorbiostat(3, palette = col)
      custom_colors = setNames(color_palette[-2], unique(df$Metric))
    }

    ggp = ggplot(df, aes(x = Component, y = Value, color = Metric)) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      scale_color_manual(values = custom_colors) +
      labs(
        title = "R² and Q² per PLS Component",
        x = "Component",
        y = "Value"
      ) +
      theme_minimal(base_size = 10 ) +
      scale_x_continuous(breaks = df$Component)
  }

  ### Scores

  if (type == "scoresX") {

    if(is.null(shape)) shape = 18
    if(is.null(labels)) labels = FALSE
    if(is.null(ellipses)) ellipses = TRUE

    x$scores = x$scoresX
    colBy = x$input$Y

    ggp = scorePlot(x,
                    comp = comp,
                    col = col,
                    colBy = colBy,
                    shape = shape,
                    shapeBy = shapeBy,
                    ellipses = ellipses,
                    labels = labels,
                    newObs = newObs)

  }

  if (type == "scoresY") {

    if(is.null(shape)) shape = 18
    if(is.null(labels)) labels = FALSE
    if(is.null(ellipses)) ellipses = FALSE

    x$scores = x$scoresY
    colBy = x$input$Y

    ggp = scorePlot(x,
                    comp = comp,
                    col = col,
                    colBy = colBy,
                    shape = shape,
                    shapeBy = shapeBy,
                    ellipses = ellipses,
                    labels = labels,
                    newObs = newObs)

  }

  ### Loadings

  if(type == 'loadings'){
    if(is.null(shape)) shape = 'arrow'
    if(is.null(labels)) labels = TRUE

    ggp = loadingPlotPLS(x,
                         comp = comp,
                         col = col,
                         colBy = colBy,
                         shape = shape,
                         selVars = selVars,
                         labels = labels,
                         labelTop = labelTop,
                         repel = repel)
  }

  if (type == "loadingsX") {

    if(is.null(shape)) shape = 'arrow'
    if(is.null(labels)) labels = TRUE

    x$loadings = x$loadingsX

    ggp = loadingPlot(x,
                      comp = comp,
                      col = col,
                      colBy = colBy,
                      shape = shape,
                      selVars = selVars,
                      labels = labels,
                      labelTop = labelTop,
                      repel = repel)

  }

  if (type == "loadingsY") {

    if(is.null(shape)) shape = 'arrow'
    if(is.null(labels)) labels = TRUE

    x$loadings = x$loadingsY

    x$explVar = data.frame("comp" = factor(1:x$ncomp),
                           "percVar" = round(100*x$summary$R2Y,4),
                           "cumPercVar" = round(100*x$summary$cumR2Y,4))

    ggp = loadingPlot(x,
                      comp = comp,
                      col = col,
                      colBy = colBy,
                      shape = shape,
                      selVars = selVars,
                      labels = labels,
                      labelTop = labelTop,
                      repel = repel)

  }

  ### Weights

  if (type == 'weights'){

    if(is.null(shape)) shape = 'arrow'
    if(is.null(labels)) labels = TRUE

    ggp = weightsPlot(x,
                      comp = comp,
                      col = col,
                      colBy = colBy,
                      shape = shape,
                      selVars = selVars,
                      labels = labels,
                      labelTop = labelTop,
                      repel = repel)
  }

  ### Correlation

  if (type == "corr") {

    if(is.null(shape)) shape = 'arrow'
    if(is.null(labels)) labels = TRUE

    x$loadings = x$loadingsX

    ggp = corrPlot(x,
                   comp = comp,
                   col = col,
                   colBy = colBy,
                   shape = shape,
                   selVars = selVars,
                   labels = labels)

  }

  ### Biplot

  if (type == "biplot"){

    if(is.null(shape)) shape = c(18,'arrow')
    if(is.null(labels)) labels = c(FALSE, TRUE)
    if(is.null(ellipses)) ellipses = FALSE

    colBy = x$input$Y

    ggp = biPlotPLS(x,
                    comp = comp,
                    col = col,
                    colBy = colBy,
                    shape = shape,
                    shapeBy = shapeBy,
                    selVars = selVars,
                    labels = labels,
                    labelTop = labelTop,
                    ellipses = ellipses,
                    repel = repel,
                    newObs = newObs)

  }

  ### Coef

  if (type == "coef"){

    if(inherits(x$coefficients_summary, "list")){

      plots = list()

      for (i in 1:length(x$coefficients_summary)) {

        coefficients_summary = x$coefficients_summary[[i]]
        coefficients_summary$Variable = rownames(coefficients_summary)
        coefficients_summary$Significant = coefficients_summary$pValJK < 0.05

        coefficients_summary$Variable = factor(coefficients_summary$Variable, levels = coefficients_summary$Variable)

        # Create plot
        ggp = ggplot(coefficients_summary, aes(x = Variable, y = coefficient)) +
          geom_col(fill = "skyblue", width = 0.6) +  # Bars
          geom_errorbar(aes(ymin = LCI_JK, ymax = UCI_JK), width = 0.2, color = "gray40") +  # Error bars
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
          geom_text(
            data = coefficients_summary[coefficients_summary$Significant == TRUE, ],
            aes(x = Variable, y = UCI_JK + 0.05, label = "*"),  # Asterisk above bar
            size = 6, color = "black"
          ) +
          labs(
            x = "Variable",
            y = "Coefficients",
            subtitle = paste("Class: ", names(x$coefficients_summary)[i])
          ) +
          theme_minimal(base_size = 12)

        plots[[i]] = ggp

      }

      ggp = patchwork::wrap_plots(plots, ncol = 2) +
        patchwork::plot_annotation(
          title = "Coefficients with Jack-Knifed Confidence Intervals",
          theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))


    } else {

      coefficients_summary = x$coefficients_summary
      coefficients_summary$Variable = rownames(coefficients_summary)
      coefficients_summary$Coefficient = coefficients_summary[,1]
      coefficients_summary$Significant = coefficients_summary$pValJK < 0.05

      coefficients_summary$Variable = factor(coefficients_summary$Variable, levels = coefficients_summary$Variable)

      # Create plot
      ggp = ggplot(coefficients_summary, aes(x = Variable, y = Coefficient)) +
        geom_col(fill = "skyblue", width = 0.6) +  # Bars
        geom_errorbar(aes(ymin = LCI_JK, ymax = UCI_JK), width = 0.2, color = "gray40") +  # Error bars
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
        geom_text(
          data = coefficients_summary[coefficients_summary$Significant == TRUE, ],
          aes(x = Variable, y = UCI_JK + 0.05, label = "*"),  # Asterisk above bar
          size = 6, color = "black"
        ) +
        labs(
          x = "Variable",
          y = "Coefficients",
          title = "Coefficients with Jack-Knifed Confidence Intervals"
        ) +
        theme_minimal(base_size = 12)

    }

  }

  ### R2

  if (type == 'R2'){

    x$loadings = x$loadingsX
    x$scores = x$scoresX
    ggp = R2varcomp(x, col)

  }

  ### Validation

  if (type == 'overfitting'){

    ggp = validationPlot(x,
                         col = col)

  }

  ### Linearity

  if (type =='linearity'){

    if(is.null(comp)) comp = 1
    scores_df = data.frame(t = x$scoresX[,comp,drop=FALSE], u = x$scoresY[,comp,drop=FALSE])
    colnames(scores_df) = c("t", "u")

    cor_coef = cor(scores_df$t, scores_df$u)

    ggp = ggplot(scores_df, aes(x = t, y = u)) +
      geom_hline(yintercept = 0, linetype = "dashed") +  # Horizontal dashed line
      geom_vline(xintercept = 0, linetype = "dashed") +  # Vertical dashed line
      coord_fixed(ratio = 1) +  # Maintain 1:1 aspect ratio and x-limits
      theme_minimal() +
      theme(legend.position = "bottom") +
      labs(title = "Scores u/t",
           x = paste0('t', comp, ' (', round(x$explVar[comp,'percVar'],2), '%)'),
           y = paste0('u', comp, ' (', round(100*x$summary[comp,'R2Y'],2), '%)'))

    ggp = ggp + geom_point(shape = 16, size = 2, color = 'skyblue3') +
      annotate("text",
               x = -Inf, y = Inf,
               label = paste0("r = ", round(cor_coef, 2)),
               hjust = -0.1, vjust = 1.1,
               size = 4, color = "black")

  }

  return(ggp)
}




