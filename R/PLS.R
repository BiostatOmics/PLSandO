
### PLS model

# x: matrix or data.frame or list (if multi-block). Ver esto porque prefiera que metamos la opcion de multiblock como algoritmo
# y: matrix or data.frame
# ncomp: If NULL, estima las maximas. Si NA, automatico.
# train: Percentage to take as training set

pls = function(x, y,
               ncomp = NULL, cvFolds = 10, rep = 10, perm = 20,
               scaling = c("none", "center", "standard", "softBlock", "hardBlock")[3],
               scalingY = c("none", "center", "standard")[2],
               algo = c("nipals","mbpls")[1],
               train = 1, alpha = 0.05, parallel = FALSE){

  if(algo=='mbpls'){

    if(any(is.na(x))) return(stop('mbpls can not be applied when missing values are present in the data. Consider using nipals algorithm instead, which supports missing values, together with a suitable block-scaling approach'))

    if (!inherits(x,'list')) return(stop('To perform multiblock PLS, please provide the block information as a list where each element corresponds to one block'))
    if (length(unique(sapply(x, nrow)))!=1) return(stop('All blocks must have the same observations'))

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

    if(!all(sapply(y, is.numeric))) {
      if(ncol(y)==1) warning('We recommed using PLS-DA instead of PLS model')
      warning('Categorical variables automatically converted to dummy variables and default filters for NAs and CV applied.\n If the user does not want to use any default filtering, consider using Preparing function before executing pca to convert categorical variables into dummies.')
      Py = Preparing(y, includeFactors = T, CVfilter = 0.01, excludeNA = 0.2)
      y = Py$x
    } else {
      Py = Preparing(y, CVfilter = 0.00001, excludeNA = 1) #Evitar nearZeroVariance variables pero no aplicar ningun filtro extra
      y = Py$x
    }

    ## Divide data into train and test

    if (train < 1) {
      ntrain = round(nrow(x[[1]])*train, 0)
      rowtrain = sample(1:nrow(x[[1]]), ntrain)
      xTest = setNames(lapply(b_names, function(y) x[[y]][-rowtrain, ,drop=FALSE]), b_names)
      x = setNames(lapply(b_names, function(y) x[[y]][rowtrain, ,drop=FALSE]), b_names)
      yTest = y[-rowtrain,, drop=FALSE ]
      y = y[rowtrain,,drop=FALSE]
    } else {
      xTest = yTest = NULL
    }
    X2 = x
    Y2 = y

    #Calcular los parametros optimos

    if(is.null(ncomp)){

      res_cross = pls_cross(X = x,Y = y,scaling = scaling, blocks = NULL, scalingY = scalingY, ncomp = ncomp, folds = cvFolds, rep = rep, algo = algo, parallel = parallel)

      r2_array = simplify2array(lapply(res_cross, function(res) res$R2))
      mean_r2_per_component = apply(r2_array, 1, mean)

      q2_array = simplify2array(lapply(res_cross, function(res) res$Q2))
      mean_q2_per_component = apply(q2_array, 1, mean)

      rmse_array = simplify2array(lapply(res_cross, function(res) res$RMSE))
      mean_rmse_per_component = apply(rmse_array, 1, mean)

      # Improvement of 0.01 in R2 and Q2 positive increment
      r2_diff = diff(mean_r2_per_component)
      q2_diff = diff(mean_q2_per_component)

      ncomp = {
        idx = which(!(diff(mean_r2_per_component) >= 0.01 & diff(mean_q2_per_component) > 0))
        if (length(idx) == 0) length(mean_r2_per_component) else idx[1]}

    } else{
      res_cross = pls_cross(X = x,Y = y,scaling = scaling, blocks = NULL, scalingY = scalingY, ncomp = ncomp, folds = cvFolds, rep = rep, algo = algo, parallel = parallel)

      q2_array = simplify2array(lapply(res_cross, function(res) res$Q2))
      mean_q2_per_component = apply(q2_array, 1, mean)

      rmse_array = simplify2array(lapply(res_cross, function(res) res$RMSE))
      mean_rmse_per_component = apply(rmse_array, 1, mean)

    }

    ## Scaling X

    if (scaling != 'none') {
      escalado = lapply(x, function(y) Scaling(y, scaling = scaling, blocks = NULL))
      x = lapply(escalado, `[[`, "x")
      centrado = lapply(escalado, `[[`, "centering")
      escalado = lapply(escalado, `[[`, "scaling")
    } else{
      escalado = lapply(b_names, function(b) rep(1, times = ncol(x[[b]])))
      centrado = lapply(b_names, function(b) rep(0, times = ncol(x[[b]])))
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

    mypls = nipals_mbpls(x,y,ncomp = ncomp)
    block_coef = lapply(mypls$block_coefficients, function(b) apply(b[,,1:ncomp, drop = FALSE], c(1,2), sum))

    block_coef_cv = setNames(lapply(seq_along(b_names), function(b) array(unlist(lapply(res_cross, function(r) r$block_coef_cv[[b]])), dim = c(dim(res_cross[[1]]$block_coef_cv[[1]])[1:3], dim(res_cross[[1]]$block_coef_cv[[1]])[4] * rep)) ), b_names)
    blockJack = setNames(lapply(seq_along(b_names), function(b) p.jack(block_coef[[b]], block_coef_cv[[b]], ncomp, alpha) ), b_names)

    if(ncol(y)==1){

      block_tables = lapply(b_names, function(b) {
        coef_mat = data.frame( "coefficient" = block_coef[[b]],
                               "pValJK" = blockJack[[b]]$pval,
                               "LCI_JK" = blockJack[[b]]$LCI_coef,
                               "UCI_JK" = blockJack[[b]]$UCI_coef,
                               "Block" = b)
      })
      resum_coef = do.call(rbind, block_tables)

    } else{ # TO DO: Revise that it is correct with multiple response

      block_tables = lapply(b_names, function(b) {

        melted = reshape2::melt(block_coef[[b]], varnames = c("variable", "response"), value.name = "coefficient")
        data.frame(
          variable   = melted$variable,
          response   = melted$response,
          coefficient = melted$coefficient,
          pValJK = as.vector(blockJack[[b]]$pval),
          LCI_JK = as.vector(blockJack[[b]]$LCI_coef),
          UCI_JK = as.vector(blockJack[[b]]$UCI_coef),
          Block = b
        ) })

      resum_tab = do.call(rbind, block_tables) # rownames(resum_tab) <- NULL
      resum_tab = split(resum_tab, resum_tab$response)

      resum_coef = lapply(resum_tab, function(df) {
        out = df[, c("variable", "coefficient", "pValJK", "LCI_JK", "UCI_JK", "Block")]
        rownames(out) = out$variable
        out$variable = NULL
        return(out)
      })
    }

    ypred = do.call('+',lapply(b_names, function(b) as.matrix(x[[b]])%*%block_coef[[b]]))
    SCEY = sum(colSums(ypred**2))

    p = sum(sapply(x, ncol))
    W = matrix(0,p,ncomp)
    for(j in 1:ncomp){
      temp = unlist(lapply(seq_along(b_names), function(b){
        temp = mypls$block_weight[[b]][,j]*mypls$weights[b,j]
        names(temp) = paste(b_names[b], names(temp), sep = "_")
        temp} ))

      W[,j] = temp
    }
    rownames(W) = names(temp)

    SCEYa =  sapply(1:ncomp, function(a) sum(tcrossprod(mypls$super_scores[, a], mypls$loadingsY[, a])^2))
    vip = sqrt(p * rowSums(sweep(W^2, 2,SCEYa, "*")) / SCEY)

    R2 = R2Y = numeric(ncomp)

    if (is.null(xTest)) {
      for (i in 1:ncomp) {
        #R2X
        xpred = setNames(lapply(b_names, function(b) tcrossprod(mypls$super_scores[,1:i,drop=FALSE], mypls$block_loadings[[b]][,1:i,drop=FALSE])),b_names)
        SCR = sum(unlist(lapply(b_names, function(b) colSums((xpred[[b]]-x[[b]])**2))))
        SCT = sum(unlist(lapply(b_names, function(b) colSums(x[[b]]**2))))
        R2[i] = 1- (SCR/SCT)
        #R2Y
        block_coef = lapply(mypls$block_coefficients, function(b) apply(b[,,1:i, drop = FALSE], c(1,2), sum))
        ypred = do.call('+',lapply(b_names, function(b) as.matrix(x[[b]])%*%block_coef[[b]]))

        SCRY = sum(colSums((ypred-y)**2))
        SCTY = sum(colSums(y**2))
        R2Y[i] = 1- (SCRY/SCTY)
      }

      r2_diff = c(R2[1], diff(R2))
      r2y_diff = c(R2Y[1],diff(R2Y))

      q2_diff = diff(mean_q2_per_component[1:ncomp])

      resum = data.frame("comp" = factor(1:ncomp),
                         "R2X" = r2_diff,
                         "cumR2X" = R2,
                         "R2Y" = r2y_diff,
                         "cumR2Y" = R2Y,
                         "Q2" = c(mean_q2_per_component[1], q2_diff),
                         "cumQ2" = mean_q2_per_component[1:ncomp],
                         "RMSE" = mean_rmse_per_component[1:ncomp])

    } else{

      xTest = setNames(lapply(seq_along(x), function(b) {
        tmp = sweep(xTest[[b]], 2, centrado[[b]], "-")
        return(sweep(tmp, 2, escalado[[b]], "/"))
      }),b_names)

      yTest = sweep(yTest, 2, centradoY, FUN = "-")
      yTest = sweep(yTest, 2, escaladoY, FUN = "/")

      R2test = R2Ytest = numeric(ncomp)

      for (i in 1:ncomp) {

        #R2X
        xpred = setNames(lapply(b_names, function(b) tcrossprod(mypls$super_scores[,1:i,drop=FALSE], mypls$block_loadings[[b]][,1:i,drop=FALSE])),b_names)
        SCR = sum(unlist(lapply(b_names, function(b) colSums((xpred[[b]]-x[[b]])**2))))
        SCT = sum(unlist(lapply(b_names, function(b) colSums(x[[b]]**2))))
        R2[i] = 1- (SCR/SCT)
        #R2Y
        block_coef = lapply(mypls$block_coefficients, function(b) apply(b[,,1:i, drop = FALSE], c(1,2), sum))
        ypred = do.call('+',lapply(b_names, function(b) as.matrix(x[[b]])%*%block_coef[[b]]))

        SCRY = sum(colSums((ypred-y)**2))
        SCTY = sum(colSums(y**2))
        R2Y[i] = 1- (SCRY/SCTY)

        #R2Xtest
        if(i == 1){
          t_b = setNames(lapply(b_names, function(b) (xTest[[b]]%*%mypls$block_weight[[b]][,1,drop=FALSE])/sqrt(ncol(xTest[[b]])) ), b_names)
          T_scores = do.call(cbind, t_b)
          super_scores = T_scores %*% mypls$weights[,1,drop=FALSE] / drop(crossprod(mypls$weights[,1,drop=FALSE]))
          xpredT = setNames(lapply(b_names, function(b) tcrossprod(super_scores, mypls$block_loadings[[b]][,1,drop=FALSE])),b_names)
          xpredTest = xpredT
        } else{
          X_def = setNames(lapply(b_names, function(b) xTest[[b]] - xpredT[[b]]), b_names)
          t_b = setNames(lapply(b_names, function(b) (X_def[[b]]%*%mypls$block_weight[[b]][,i,drop=FALSE])/sqrt(ncol(xTest[[b]])) ), b_names)
          T_scores = do.call(cbind, t_b)
          super_scores = T_scores %*% mypls$weights[,i,drop=FALSE] / drop(crossprod(mypls$weights[,i,drop=FALSE]))
          xpredT = setNames(lapply(b_names, function(b) tcrossprod(super_scores, mypls$block_loadings[[b]][,i,drop=FALSE])),b_names)
          xpredTest = setNames(lapply(b_names, function(b) xpredTest[[b]] + xpredT[[b]] ),b_names)
        }

        SCR = sum(unlist(lapply(b_names, function(b) colSums((xpredTest[[b]]-xTest[[b]])**2))))
        SCT = sum(unlist(lapply(b_names, function(b) colSums(xTest[[b]]**2))))
        R2test[i] = 1- (SCR/SCT)
        #R2Ytest
        ypred = do.call('+',lapply(b_names, function(b) as.matrix(xTest[[b]])%*%block_coef[[b]]))
        SCRY = sum(colSums((ypred-yTest)**2))
        SCTY = sum(colSums(yTest**2))
        R2Ytest[i] = 1- (SCRY/SCTY)
      }
      r2test_diff = diff(R2test)
      r2ytest_diff = diff(R2Ytest)

      r2_diff = c(R2[1], diff(R2))
      r2y_diff = c(R2Y[1],diff(R2Y))

      q2_diff = diff(mean_q2_per_component[1:ncomp])

      resum = data.frame("comp" = factor(1:ncomp),
                         "R2X" = r2_diff,
                         "cumR2X" = R2,
                         "R2Y" = r2y_diff,
                         "cumR2Y" = R2Y,
                         "R2Xtest" = c(R2test[1], r2test_diff),
                         "cumR2Xtest" = R2test,
                         "R2Ytest" = c(R2Ytest[1], r2ytest_diff),
                         "cumR2Ytest" = R2Ytest,
                         "Q2" = c(mean_q2_per_component[1], q2_diff),
                         "cumQ2" = mean_q2_per_component[1:ncomp],
                         "RMSE" = mean_rmse_per_component[1:ncomp])

    }

    bip = sqrt(length(b_names)*rowSums( (mypls$weights^2) * matrix(r2y_diff, nrow = length(b_names), ncol = ncomp, byrow = TRUE) )/mean_r2_per_component[ncomp])
    names(bip) = b_names

    res_val = plsValidation(X2,Y2,scaling, scalingY,ncomp,cvFolds,perm,algo,parallel)

    explVarBlock = setNames(lapply(b_names, function(y){
      expl = data.frame("comp" = factor(1:ncomp),
                        "percVar" = round(c(mypls$explvarB[[y]][1] ,diff(as.numeric(mypls$explvarB[[y]]))),4),
                        "cumPercVar" = t(round(mypls$explvarB[[y]],4)))
    }),b_names)

    desSummary = setNames(lapply(b_names, function(y){
      des = data.frame(c("Samples","Variables","Excluded Near-Zero Variables","Missing Values"),
                       c(nrow(x[[y]]),ncol(x[[y]]),length(which(P[[y]]$removed$Variables$Problem == 'Low CV')),
                         paste0(length(which(is.na(x[[y]]))), " (",round((length(which(is.na(x[[y]]))) / (nrow(x[[y]]) * ncol(x[[y]]))) * 100, 2), "%)")), stringsAsFactors = FALSE)
      colnames(des) = NULL
      return(des)
    }),b_names)

    return(list("SscoresX" = mypls$super_scores, #super_scores
                "BscoresX" = mypls$block_scores,
                "BloadingsX" = mypls$block_loadings,
                "scoresY" = mypls$scoresY,
                "loadingsY" = mypls$loadingsY,
                "weights" = mypls$weights,
                "Bweights" = mypls$block_weight,
                "coefficients" = mypls$block_coefficients,
                "BexplVar" = explVarBlock,
                "vip" = vip,
                "bip" = bip,
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
                "input" = list("scalingType" = scaling,"scalingTypeY" = scalingY, "X" = X2, "Y" = Y2 ,
                               'algo' = algo, "alpha" = alpha, "perm" = perm, "blocks" = NULL)))

    ## ACABA EL MBPLS

  } else{

    ## EMPIEZA EL PLS

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

    if (inherits(x,'list')) {
      if (!scaling %in% c("softBlock", "hardBlock")) warning('When blocks are considered we recommend using either softBlock or hardBlock scaling')
      if (length(unique(sapply(x, nrow)))!=1) return(stop('All blocks must have the same observations'))
      blocks = rep(seq_along(x), times = sapply(x, ncol) )
      x = do.call(cbind, x)
    }

    if(!all(sapply(y, is.numeric))) {
      if(ncol(y)==1) warning('We recommed using PLS-DA instead of PLS model')
      warning('Categorical variables automatically converted to dummy variables and default filters for NAs and CV applied.\n If the user does not want to use any default filtering, consider using Preparing function before executing pca to convert categorical variables into dummies.')
      Py = Preparing(y, includeFactors = T, CVfilter = 0.01, excludeNA = 0.2)
      y = Py$x
    } else {
      Py = Preparing(y, CVfilter = 0.00001, excludeNA = 1) #Evitar nearZeroVariance variables pero no aplicar ningun filtro extra
      y = Py$x
    }

    ## Divide data into train and test

    if (train < 1) {
      ntrain = round(nrow(x)*train, 0)
      rowtrain = sample(1:nrow(x), ntrain)
      xTest = x[-rowtrain, ,drop=FALSE]
      x = x[rowtrain,, drop=FALSE]
      yTest = y[-rowtrain,, drop=FALSE ]
      y = y[rowtrain,,drop=FALSE]
    } else {
      xTest = yTest = NULL
    }
    X2 = x
    Y2 = y

    #Calcular los parametros optimos

    if(is.null(ncomp)){

      res_cross = pls_cross(X = x,Y = y,scaling = scaling, blocks = blocks, scalingY = scalingY, ncomp = ncomp, folds = cvFolds, rep = rep, algo = algo, parallel = parallel)

      r2_array = simplify2array(lapply(res_cross, function(res) res$R2))
      mean_r2_per_component = apply(r2_array, 1, mean)

      q2_array = simplify2array(lapply(res_cross, function(res) res$Q2))
      mean_q2_per_component = apply(q2_array, 1, mean)

      rmse_array = simplify2array(lapply(res_cross, function(res) res$RMSE))
      mean_rmse_per_component = apply(rmse_array, 1, mean)

      # Improvement of 0.01 in R2 and Q2 positive increment
      r2_diff = diff(mean_r2_per_component)
      q2_diff = diff(mean_q2_per_component)

      ncomp = {
        idx = which(!(diff(mean_r2_per_component) >= 0.01 & diff(mean_q2_per_component) > 0))
        if (length(idx) == 0) length(mean_r2_per_component) else idx[1]}

    } else{
      res_cross = pls_cross(X = x,Y = y,scaling = scaling, blocks = blocks, scalingY = scalingY, ncomp = ncomp, folds = cvFolds, rep = rep, algo = algo, parallel = parallel)

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

    if(ncol(y)==1){
      resum_coef = data.frame("coefficient" = mypls$coefficients,
                       "pValJK" = Jack$pval,
                       "LCI_JK" = Jack$LCI_coef,
                       "UCI_JK" = Jack$UCI_coef)
      names(resum_coef)[1] = "coefficient"
    } else{
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
    }

    ypred = as.matrix(x) %*% mypls$coefficients
    SCEY = sum(colSums(ypred**2))

    SCEYa =  sapply(1:ncomp, function(j) sum(tcrossprod(mypls$scores[, j], mypls$loadingsY[, j])^2))
    vip = sqrt(ncol(x) * rowSums(sweep(mypls$weights^2, 2,SCEYa, "*")) / SCEY)

    R2 = R2Y = numeric(ncomp)

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
                         "RMSE" = mean_rmse_per_component[1:ncomp])

    } else{

      xTest = sweep(xTest, 2, centrado, FUN = "-")
      xTest = sweep(xTest, 2, escalado, FUN = "/")

      yTest = sweep(yTest, 2, centradoY, FUN = "-")
      yTest = sweep(yTest, 2, escaladoY, FUN = "/")

      R2test = R2Ytest = numeric(ncomp)

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
                         "RMSE" = mean_rmse_per_component[1:ncomp])

    }

    res_val = plsValidation(X2,Y2,scaling, scalingY,ncomp,cvFolds,perm,algo,parallel)

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
                "input" = list("scalingType" = scaling,"scalingTypeY" = scalingY, "X" = X2, "Y" = Y2 ,
                               'algo' = algo, "alpha" = alpha, "perm" = perm, "blocks" = blocks)))
  }

}




### PLS predictions

# x: PLS object

plsPredict = function(x, new = NULL, plot = TRUE) {
  ## proyectar nuevas observaciones antes de predecir

  if(x$input$algo=='nipals'){

    x$explVar = data.frame("comp" = factor(1:ncomp),
                           "percVar" = round(100*x$summary$R2X,4),
                           "cumPercVar" = round(100*x$summary$cumR2X,4))

    if(plot) scorePlot(x, newObs = new)

    new = as.data.frame(new)
    new = new[,colnames(x$X)]

    new = sweep(new, 2, x$scaling$center, FUN = "-")
    new = sweep(new, 2, x$scaling$scale, FUN = "/")

    prediction = as.matrix(new) %*% x$coefficients

  } else{

    if(!inherits(new, "list")) return(stop('Please provide the new data structured in blocks'))

    x$explVar = data.frame("comp" = factor(1:ncomp),
                           "percVar" = round(100*x$summary$R2X,4),
                           "cumPercVar" = round(100*x$summary$cumR2X,4))

    if(plot) scorePlotmb(x, newObs = new)

    new = lapply(new, as.matrix)
    b_names = names(x$X)
    new = lapply(b_names, function(b) new[[b]][,colnames(x$X[[b]])])

    centrado = x$scaling$center
    escalado = x$scaling$scaling

    new = setNames(lapply(seq_along(x), function(b) {
      tmp = sweep(xTest[[b]], 2, centrado[[b]], "-")
      return(sweep(tmp, 2, escalado[[b]], "/"))
    }),b_names)

    block_coef = lapply(mypls$coefficients, function(b) apply(b[,,1:i, drop = FALSE], c(1,2), sum))
    prediction = do.call('+',lapply(b_names, function(b) as.matrix(new[[b]])%*%block_coef[[b]]))

  }

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

plsPlot = function(x,
                   type = c("R2vsQ2", "loadings", "scoresX",  "loadingsX", "scoresY", "loadingsY",
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

  if(!(type %in% c("R2vsQ2", "loadings", "scoresX",  "loadingsX", "scoresY", "loadingsY", "weights", "linearity", "overfitting", "R2", "corr", "biplot", "coef"))) return(stop('Please use one of: R2vsQ2, loadings, scoresX, loadingsX, scoresY, loadingsY, weights, linearity, overfitting, R2, correl, biplot.'))

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

  ### Validation

  if (type == 'overfitting'){

    ggp = validationPlot(x,
                         col = col)

  }

  ### Coef

  #TO DO: Quizas estaria bien indicar de que bloque provienen

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
            subtitle = paste("Response variable: ", names(x$coefficients_summary)[i])
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


  if(x$input$algo == 'nipals'){

    eigenvalues = apply(x$scoresX, 2, function(t) crossprod(t)/(nrow(x$X)-1) )

    x$explVar = data.frame("comp" = factor(1:x$ncomp),
                           "percVar" = round(100*x$summary$R2X,4),
                           "cumPercVar" = round(100*x$summary$cumR2X,4),
                           "eigenVal" = eigenvalues)

    ### Scores

    if (type == "scoresX") {

      if(is.null(shape)) shape = 18
      if(is.null(labels)) labels = FALSE
      if(is.null(ellipses)) ellipses = TRUE

      x$scores = x$scoresX

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
      if(is.null(ellipses)) ellipses = TRUE

      x$scores = x$scoresY

      x$explVar = data.frame("comp" = factor(1:x$ncomp),
                             "percVar" = round(100*x$summary$R2Y,4),
                             "cumPercVar" = round(100*x$summary$cumR2Y,4))

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

    ### R2

    if (type == 'R2'){

      x$loadings = x$loadingsX
      x$scores = x$scoresX
      ggp = R2varcomp(x, col)

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

  } else{

    x$explVar = data.frame("comp" = factor(1:x$ncomp),
                           "percVar" = round(100*x$summary$R2X,4),
                           "cumPercVar" = round(100*x$summary$cumR2X,4))

    ## Scores

    if (type == "scoresX") {


    if(is.null(shape)) shape = 18
    if(is.null(labels)) labels = FALSE
    if(is.null(ellipses)) ellipses = TRUE

    x$scores = x$SscoresX

    ggp = scorePlotmb(x,
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
      if(is.null(ellipses)) ellipses = TRUE

      x$explVar = data.frame("comp" = factor(1:x$ncomp),
                             "percVar" = round(100*x$summary$R2Y,4),
                             "cumPercVar" = round(100*x$summary$cumR2Y,4))

      x$scores = x$scoresY

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

    ## Loadings

    if (type == "loadings") {

      if(is.null(shape)) shape = 'arrow'
      if(is.null(labels)) labels = TRUE

      ggp = loadingPlotPLSmb(x,
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

      ggp = loadingPlotmb(x,
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

      ggp = weightsPlotmb(x,
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

    ### biPlot

    if (type == 'biPlot'){

      #TO DO: Verificar que esta funcionando

      if(is.null(shape)) shape = 'arrow'
      if(is.null(labels)) labels = TRUE

      ggp = biPlotmb(x,
                     comp = comp,
                     col = col,
                     colBy = colBy,
                     shape = shape,
                     selVars = selVars,
                     labels = labels,
                     labelTop = labelTop,
                     repel = repel)
    }

    ### R2

    if (type == 'R2'){

      ggp = R2varcompmb(x, col)

    }

    ### Linearity

    if (type =='linearity'){

      if(is.null(comp)) comp = 1
      scores_df = data.frame(t = x$SscoresX[,comp,drop=FALSE], u = x$scoresY[,comp,drop=FALSE])
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

  }

  return(ggp)
}



### PLS outliers

# x: Object return by pca() function
# K: Number of componentes to be used for outlier detection. If NULL, all components in the object will be used.

plsOutliers = function(x, K = NULL,
                       method = c("T2", "RSS", "both")[3],
                       conf = 99) {

  if (is.null(K)) K = x$ncomp

  if (method == "T2") {

    out = T2plot(x = x, K = K, conf = conf)
    return(out)

  }

  if (method == "RSS") {

    out = RSSplot(x = x, K = K, conf = conf)
    return(out)

  }

  if (method == "both") {

    par(mfrow = c(1,2))
    outS = T2plot(x = x, K = K, conf = conf)
    outM = RSSplot(x = x, K = K, conf = conf)
    par(mfrow = c(1,1))

    return(c(outS, outM))

  }



}


# x:  pca object
# outliers:  output of pcaOutliers

outlierContrib = function(x, outliers, labelSize = 1) {

  # Severe outliers (T2)
  severos = outliers$SevereOutliers
  num = length(severos)

  if (num == 0) {
    cat("There are no severe outliers in the data or they were not computed.\n")
    miscontrT2 = NULL
  } else {
    cat (paste0("There are ", num, " severe outliers in the data.\n"))
    miscontrT2 = contribT2(X = x$X, scores = x$scores, loadings = x$loadings,
                           eigenval = x$explVar$eigenVal, observ = severos,
                           cutoff = 2)
    colnames(miscontrT2) = severos
    if (num > 10) {
      cat (paste0("Only the 10 most severe outliers will be plotted.\n"))
      num = 10
      severos = rownames(sort(outliers$T2, decreasing = TRUE)[1:10])
    }

    for (i in severos) {
      barplot(sort(miscontrT2[,i], decreasing = TRUE),
              las=2, cex.names = labelSize, border = NA, col = "red3",
              main = paste0("Severe outlier: ", i))
    }

  }


  # Moderate outliers (RSS)
  modera = outliers$ModerateOutliers
  num = length(modera)

  if (num == 0) {
    cat("There are no moderate outliers in the data or they were not computed.\n")
    miscontrRSS = NULL
  } else {
    cat (paste0("There are ", num, " moderate outliers in the data.\n"))
    miscontrRSS = ContriSCR(E = outliers$E, SCR = outliers$RSS)
    miscontrRSS = miscontrRSS[modera,,drop = FALSE]
    if (num > 10) {
      cat (paste0("Only the 10 most moderate outliers will be plotted.\n"))
      num = 10
      modera = rownames(sort(outliers$RSS, decreasing = TRUE)[1:10])
    }
    for (i in modera) {
      barplot(sort(miscontrRSS[i,], decreasing = TRUE),
              las=2, cex.names = labelSize, border = NA, col = "blue4",
              main = paste0("Moderate outlier: ", i))
    }
  }

  return(list(contribT2 = t(miscontrT2), contribRSS = miscontrRSS))

}

### NIPALS pseudo-imputation

#x :  pca object
#scaled: If we want to return the values in the ranges of the original values or not

imputeNipals = function(x, scaled = FALSE){
  X = x$X
  Xest = x$scores  %*% t(x$loadings)
  X[is.na(X)] = Xest[is.na(X)]
  if(!scaled){
    X = sweep(X, 2, x$scaling$scale, FUN = "*")
    X = sweep(X, 2, x$scaling$center, FUN = "+")
  }
  return(X)

}

#### Variable selection

plsVarSel = function(x,
                     type = c('Jack','Perm','VIP','sMC')[1],
                     cvFolds = 5, rep = 10, threshold = NULL, parallel = FALSE){

  if (!(type%in%c('Jack','Perm','VIP','sMC'))) return(stop('Please use one of: Jack, Perm, VIP, sMC.'))

  if(type == 'Jack'){
    #Mirar para los casos multibloque
    if (is.null(threshold)) threshold = 0.05

    rJK = JK(x, cvFolds = cvFolds, rep = rep, threshold = threshold, parallel = parallel)

    if(ncol(x$Y)==1){

      if(x$input$algo=='nipals'){
        coefficients_summary = data.frame("Coefficient" = x$coefficients,
                                          "pValJK" = rJK$pval,
                                          "LCI_JK" = rJK$LCI_coef,
                                          "UCI_JK" = rJK$UCI_coef)
        vars = rJK$pval
        rownames(vars) = rownames(x$coefficients)
        vars = vars[vars<threshold,,drop=FALSE]
      } else{
        block_tables = lapply(b_names, function(b) {
          coef_mat = data.frame( "coefficient" = block_coef[[b]],
                                 "pValJK" = blockJack[[b]]$pval,
                                 "LCI_JK" = blockJack[[b]]$LCI_coef,
                                 "UCI_JK" = blockJack[[b]]$UCI_coef,
                                 "Block" = b)
        })
        coefficients_summary = do.call(rbind, block_tables)

        vars = coefficients_summary[,'pValJK',drop=FALSE]
        vars = vars[vars<threshold,,drop=FALSE]
      }

      names(coefficients_summary)[1] = "Coefficient"
      coefficients_summary$Variable = rownames(coefficients_summary)
      coefficients_summary$Significant = coefficients_summary$pValJK < 0.05

      coefficients_summary$Variable = factor(coefficients_summary$Variable, levels = coefficients_summary$Variable)

      # Create plot
      ggplot(coefficients_summary, aes(x = Variable, y = Coefficient)) +
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


    } else{

      resum_tab = data.frame(
        variable   = rownames(x$coefficients),
        reshape2::melt(x$coefficients, varnames = c("variable", "response"), value.name = "coefficient")[, -1],
        pValJK = as.vector(rJK$pval),
        LCI_JK = as.vector(rJK$LCI_coef),
        UCI_JK = as.vector(rJK$UCI_coef))

      resum_tab = split(resum_tab, resum_tab$response)

      vars = lapply(resum_tab, function(df) {
        out = df[, c("variable", "pValJK")]
        rownames(out) = out$variable
        out$variable = NULL
        return(out) })

      resum_coef = lapply(resum_tab, function(df) {
        out = df[, c("variable", "coefficient", "pValJK", "LCI_JK", "UCI_JK")]
        rownames(out) = out$variable
        out$variable = NULL
        return(out) })

      plots = list()

      for (i in 1:length(resum_coef)) {

        coefficients_summary = resum_coef[[i]]
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
            subtitle = paste("Response variable: ", names(x$coefficients_summary)[i])
          ) +
          theme_minimal(base_size = 12)

        plots[[i]] = ggp

      }

      patchwork::wrap_plots(plots, ncol = 2) +
        patchwork::plot_annotation(
          title = "Coefficients with Jack-Knifed Confidence Intervals",
          theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0)))

    }

    #Anadir los plots de los coeficientes de regresion
  }

  if(type == 'Perm'){

    #Mirar para los casos multibloque
    if (is.null(threshold)) threshold = 0.05
    rPerm = perm(x, R = rep, threshold = threshold)

    if(ncol(x$Y)==1){

      if(x$input$algo=='nipals'){
        coefficients_summary = data.frame("Coefficient" = x$coefficients,
                                          "pValPerm" = as.vector(rPerm$pval),
                                          "LCI_Perm" = rPerm$LCI_coef,
                                          "UCI_Perm" = rPerm$UCI_coef)
      } else{
        coefficients_summary = data.frame("Coefficient" = x$coefficients_summary[,'coefficient',drop=FALSE],
                                          "pValPerm" = as.vector(rPerm$pval),
                                          "LCI_Perm" = rPerm$LCI_coef,
                                          "UCI_Perm" = rPerm$UCI_coef)
      }

      names(coefficients_summary)[1] = "Coefficient"
      coefficients_summary$Variable = rownames(coefficients_summary)
      coefficients_summary$Significant = coefficients_summary$pValPerm < 0.05

      coefficients_summary$Variable = factor(coefficients_summary$Variable, levels = coefficients_summary$Variable)

      # Create plot
      ggplot(coefficients_summary, aes(x = Variable, y = Coefficient)) +
        geom_col(fill = "skyblue", width = 0.6) +  # Bars
        geom_errorbar(aes(ymin = LCI_Perm, ymax = UCI_Perm), width = 0.2, color = "gray40") +  # Error bars
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
        geom_text(
          data = coefficients_summary[coefficients_summary$Significant == TRUE, ],
          aes(x = Variable, y = UCI_Perm + 0.05, label = "*"),  # Asterisk above bar
          size = 6, color = "black"
        ) +
        labs(
          x = "Variable",
          y = "Coefficients",
          title = "Coefficients with Confidence Intervals obtained from permutation testing",
          subtitle = 'IC defines the range of non-significance'
        ) +
        theme_minimal(base_size = 12)

      vars = coefficients_summary[,'pValPerm',drop=FALSE]
      vars = vars[vars<threshold,,drop=FALSE]
    } else{

      resum_tab = data.frame(
        variable   = rownames(x$coefficients),
        reshape2::melt(x$coefficients, varnames = c("variable", "response"), value.name = "coefficient")[, -1],
        pValPerm = as.vector(rPerm$pval),
        LCI_Perm = as.vector(rPerm$LCI_coef),
        UCI_Perm = as.vector(rPerm$UCI_coef))

      resum_tab = split(resum_tab, resum_tab$response)

      vars = lapply(resum_tab, function(df) {
        out = df[, c("variable", "pValPerm")]
        rownames(out) = out$variable
        out$variable = NULL
        return(out) })

      resum_coef = lapply(resum_tab, function(df) {
        out = df[, c("variable", "coefficient", "pValPerm", "LCI_Perm", "UCI_Perm")]
        rownames(out) = out$variable
        out$variable = NULL
        return(out) })

      plots = list()

      for (i in 1:length(resum_coef)) {

        coefficients_summary = resum_coef[[i]]
        coefficients_summary$Variable = rownames(coefficients_summary)
        coefficients_summary$Significant = coefficients_summary$pValPerm < 0.05

        coefficients_summary$Variable = factor(coefficients_summary$Variable, levels = coefficients_summary$Variable)

        # Create plot
        ggp = ggplot(coefficients_summary, aes(x = Variable, y = coefficient)) +
          geom_col(fill = "skyblue", width = 0.6) +  # Bars
          geom_errorbar(aes(ymin = LCI_Perm, ymax = UCI_Perm), width = 0.2, color = "gray40") +  # Error bars
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
          geom_text(
            data = coefficients_summary[coefficients_summary$Significant == TRUE, ],
            aes(x = Variable, y = UCI_Perm + 0.05, label = "*"),  # Asterisk above bar
            size = 6, color = "black"
          ) +
          labs(
            x = "Variable",
            y = "Coefficients",
            subtitle = paste("Response variable: ", names(x$coefficients_summary)[i])
          ) +
          theme_minimal(base_size = 12)

        plots[[i]] = ggp

      }

      patchwork::wrap_plots(plots, ncol = 2) +
        patchwork::plot_annotation(
          title = "Coefficients with Confidence Intervals obtained from permutation testing",
          subtitle = 'IC defines the range of non-significance',
          theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0)))


    }
  }

  if(type == 'VIP'){
    if (is.null(threshold)) threshold = 1
    vars = as.data.frame(x$vip[x$vip>threshold])
    colnames(vars) = 'VIP'
  }

  if(type == 'sMC'){

    if(x$input$algo == 'mbpls') return(stop('Not available for multiblock models'))
    if (is.null(threshold)) threshold = 0.05

    if(ncol(x$Y)==1){
      rsMC = SMC(x$X, x$coefficients)
      fcrit = qf(1-threshold,1,nrow(x$X)-2)
      vars = as.data.frame(rsMC[rsMC>fcrit])
      colnames(vars) = 'sMC'
    } else{

      vars = list()
      for (i in 1:ncol(x$Y)) {
        rsMC = SMC(x$X, x$coefficients[,i,drop=FALSE])
        fcrit = qf(1-threshold,1,nrow(x$X)-2)
        bvars = as.data.frame(rsMC[rsMC>fcrit])
        colnames(bvars) = 'sMC'
        vars[[i]] = bvars
      }
      names(vars) = colnames(x$Y)
    }

  }

  return(vars)
}





