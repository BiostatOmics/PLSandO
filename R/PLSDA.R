### PLS-DA model

#' Partial Least Squares Discriminant Analysis (PLS-DA) and Multi-block PLS-DA regression models
#'
#' @description Performs Partial Least Squares Discriminant Analysis and multi-block PLS-DA regression models using the NIPALS algorithm.
#' Includes automatic preprocessing, filtering of low-variance variables, automatic cross-validation for optimal
#' component selection and additional functionalities for model validation.
#'
#' @param x Matrix, data.frame, or list of blocks (for MBPLSDA algorithm).
#' @param y Matrix or data.frame containing the response variables.
#' @param ncomp Number of components. If \code{NULL}, it is estimated via cross-validation.
#' @param cvFolds Number of cross-validation folds. Default is 10.
#' @param rep Number of cross-validation repetitions. Default is 10.
#' @param perm Number of permutations for model validation. Default is 20.
#' @param scaling Scaling method for X: "none", "center", "standard" (default), "softBlock", or "hardBlock".
#' @param scalingY Scaling method for Y: "none", "center" (default), or "standard".
#' @param algo Algorithm to use: "nipals" (standard PLS, default) or "mbpls" (Multi-block PLS-DA).
#' @param train Training set percentage (0 to 1). Default is 1 (all data).
#' @param alpha Significance level for Jack-knife confidence intervals. Default is 0.05.
#' @param parallel Logical. If \code{TRUE}, cross-validation runs in parallel.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{scoresX} or \code{Sscores}: Matrix of scores for the X block (Super-scores in MB-PLS).
#'   \item \code{loadingsX} or \code{BloadingsX}: Matrix or list of loadings for the X block (Block-loadings in MB-PLS).
#'   \item \code{scoresY}: Matrix of scores for the Y block.
#'   \item \code{loadingsY}: Matrix of loadings for the Y block.
#'   \item \code{weights} or \code{Bweights}: Matrix or list of weights (Block-weights in MB-PLS).
#'   \item \code{Sweights}: Matrix of Super-weights (only for MB-PLS).
#'   \item \code{BscoresX}: List of block scores (only for MB-PLS).
#'   \item \code{coefficients}: Matrix or list (in multiresponse scenarios or multi-block) of regression coefficients for the model.
#'   \item \code{coefficients_summary}: Data frame or list with coefficients, Jack-knife p-values, and confidence intervals.
#'   \item \code{vip}: Vector of Variable Importance in Projection (VIP) values.
#'   \item \code{bip}: Vector of Block Importance in Projection (BIP) values (only for MB-PLS).
#'   \item \code{X}: The preprocessed and scaled data matrix (or list of matrices) used for model fitting.
#'   \item \code{Y}: The preprocessed and scaled Y data matrix.
#'   \item \code{scaling}: Data frame or list (only in MB-PLS) with centering and scaling parameters of each variable for both X and Y blocks.
#'   \item \code{summary}: Data frame with R2X, R2Y, Q2, RMSE, F1score, Mathews Correlation Coefficient (MCC), Balanced Error Rate (BER) and Total Error (Etot) metrics per component, including test set results if applicable.
#'   \item \code{BexplVar}: List of explained variance per block (only for MB-PLS).
#'   \item \code{validation}: Data frame with R2, Q2, PRESS, RMSE and Similarity obtained in permutation testing.
#'   \item \code{cv_results}: Detailed results from the cross-validation procedure.
#'   \item \code{PreproSummary}: Summary of data dimensions, missing values, and excluded low-variance variables.
#'   \item \code{ncomp}: The final number of components used in the model (optimized via cross-validation).
#'   \item \code{test}: List containing the processed \code{xTest} and \code{yTest} data if a training percentage was specified.
#'   \item \code{input}: List of input parameters used (\code{scalingType}, \code{algo}, \code{alpha}, original X and Y, etc.).
#' }
#' @export

plsda = function(x, y,
                 ncomp = NULL, cvFolds = 10, rep = 10, perm = 20,
                 scaling = c("none", "center", "standard", "softBlock", "hardBlock")[3],
                 scalingY = c("none", "center", "standard")[2],
                 algo = c("nipals","mbplsda")[1],
                 train = 1, alpha = 0.05, parallel = FALSE){

  if(algo=='mbplsda'){

    if(any(is.na(x))) return(stop('mbplsda can not be applied when missing values are present in the data. Consider using nipals algorithm instead, which supports missing values, together with a suitable block-scaling approach'))

    if (!inherits(x,'list')) return(stop('To perform multiblock PLSDA, please provide the block information as a list where each element corresponds to one block'))
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

    X = x
    Y = y

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

      res_cross = plsda_cross(X = x, Y = y, Y2=Y, scaling = scaling, blocks = NULL, scalingY = scalingY, ncomp = ncomp, folds = cvFolds, rep = rep, parallel = parallel, algo = algo)

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

      ncomp = if(length(mean_ber_per_component)< 3) which.min(mean_ber_per_component) else kneedle::kneedle(1:length(mean_ber_per_component), mean_ber_per_component, decreasing = T)[1] #BER    La libreria esta basada en el paper Finding a Kneedle in a Haystack: Detecting Knee Points in System Behavior

    } else{
      res_cross = plsda_cross(X = x, Y = y, Y2=Y, scaling = scaling, blocks = NULL, scalingY = scalingY, ncomp = ncomp, folds = cvFolds, rep = rep, parallel = parallel, algo = algo)

      q2_array = simplify2array(lapply(res_cross, function(res) res$Q2))
      mean_q2_per_component = if(inherits(q2_array,'numeric')) mean(q2_array) else apply(q2_array, 1, mean)

      rmse_array = simplify2array(lapply(res_cross, function(res) res$RMSE))
      mean_rmse_per_component = if(inherits(q2_array,'numeric')) mean(rmse_array) else apply(rmse_array, 1, mean)

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

    } else{

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

    R2 = R2Y = f1score = MCC = Etot = BER = numeric(ncomp)

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
                         "F1score" = f1score,
                         "MCC" = MCC,
                         "Etot" = Etot,
                         "BER" = BER)

    } else{

      xTest = setNames(lapply(seq_along(x), function(b) {
        tmp = sweep(xTest[[b]], 2, centrado[[b]], "-")
        return(sweep(tmp, 2, escalado[[b]], "/"))
      }),b_names)

      yTest = sweep(yTest, 2, centradoY, FUN = "-")
      yTest = sweep(yTest, 2, escaladoY, FUN = "/")

      R2test = R2Ytest = f1scoretest = MCCtest = BERtest = Etottest = numeric(ncomp)

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

    bip = sqrt(length(b_names)*rowSums( (mypls$weights^2) * matrix(r2y_diff, nrow = length(b_names), ncol = ncomp, byrow = TRUE) )/R2Y[ncomp])
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
                "Sweights" = mypls$weights,
                "Bweights" = mypls$block_weight,
                "coefficients" = mypls$block_coefficients,
                "BexplVar" = explVarBlock,
                "vip" = vip,
                "bip" = bip,
                "X" = x, # escalada y sin variables con baja variabilidad
                "Y" = y, # escalada
                "scaling" = list("center" = centrado, "scale" = escalado, #para que cuadre con los scorePLot de cuando hicimos PCA
                                 "centerY" = centradoY, "scaleY" = escaladoY),
                "summary" = resum, ## completar si se desea dar + info de lo que se ha hecho
                "coefficients_summary" = resum_coef,
                "ncomp" = ncomp,
                "test" = list("yTest" = yTest, "xTest" = xTest),
                "validation" = as.data.frame(res_val),
                "cv_results" = res_cross,
                "PreproSummary" = desSummary,
                "input" = list("scalingType" = scaling,"scalingTypeY" = scalingY, "X" = X, "Y" = Y,
                               'algo' = algo, "alpha" = alpha, "perm" = perm, "blocks" = NULL, "model" = 'plsda')))

    ## ACABA EL MBPLSDA

  } else{

    #EMPIEZA EL PLSDA
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
      missNATest = any(is.na(xTest))
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

    missNA = any(is.na(x))

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

      ncomp = if(length(mean_ber_per_component)< 3) which.min(mean_ber_per_component) else kneedle::kneedle(1:length(mean_ber_per_component), mean_ber_per_component, decreasing = T)[1] #BER    La libreria esta basada en el paper Finding a Kneedle in a Haystack: Detecting Knee Points in System Behavior

    } else{
      res_cross = plsda_cross(X = x, Y = y, Y2=Y, scaling = scaling, blocks = blocks, scalingY = scalingY, ncomp = ncomp, folds = cvFolds, rep = rep, parallel = parallel, algo = algo)

      q2_array = simplify2array(lapply(res_cross, function(res) res$Q2))
      mean_q2_per_component = if(inherits(q2_array,'numeric')) mean(q2_array) else apply(q2_array, 1, mean)

      rmse_array = simplify2array(lapply(res_cross, function(res) res$RMSE))
      mean_rmse_per_component = if(inherits(q2_array,'numeric')) mean(rmse_array) else apply(rmse_array, 1, mean)
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

    if(missNA) x = impute_nipals(X = x, mypls = mypls)
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

        if(missNATest) xTest = impute_nipals(xTest, yTest, inObs = FALSE, ncomp = i)
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
                               'algo' = algo, "alpha" = alpha, "perm" = perm, "blocks" = blocks, "model" ='plsda')))

  }
}


### PLSDA predictions

### PLSDA predictions

#' Prediction for PLS-DA Models
#'
#' @description Uses previously fitted PLSDA or MB-PLSDA model to obtain the predictions of new observations.
#'
#' @param x A plsda object returned by the \code{plsda()} function.
#' @param new A matrix, data.frame, or list of blocks (for MB-PLSDA or for PLSDA with block-scaling) containing
#'   the new observations to be predicted.
#' @param plot Logical. If \code{TRUE}, the function generates a score plot
#'   showing the projection of the new observations.
#'
#' @return A matrix containing the predicted values for the response variable(s) Y.
#' @export

plsdaPredict = function(x, new = NULL, plot = TRUE) {
  ## proyectar nuevas observaciones antes de predecir

  if(x$input$algo=='nipals'){

    x$explVar = data.frame("comp" = factor(1:x$ncomp),
                           "percVar" = round(100*x$summary$R2X,4),
                           "cumPercVar" = round(100*x$summary$cumR2X,4))
    x$scores = x$scoresX

    if(plot) print(scorePlot(x, newObs = new))

    new = as.data.frame(new)
    new = new[,colnames(x$X)]

    new = sweep(new, 2, x$scaling$center, FUN = "-")
    new = sweep(new, 2, x$scaling$scale, FUN = "/")

    if(any(rowSums(is.na(new))>0)){
      scores = project.obs.nipals.pls(x, new, 1:x$ncomp)
      prediction = tcrossprod(scores, x$loadingsY)
    } else{
      prediction = as.matrix(new) %*% x$coefficients
    }
    prediction = sweep(prediction, 2, x$scaling$scaleY, FUN = "*")
    plsprediction = sweep(prediction, 2, x$scaling$centerY, FUN = "+")
    prediction = as.data.frame(sub(".*_","", colnames(plsprediction)[apply(plsprediction, 1, which.max)]))
    colnames(prediction) = unique(sub("_.*","", colnames(plsprediction)[apply(plsprediction, 1, which.max)]))

  } else{

    if(!inherits(new, "list")) return(stop('Please provide the new data structured in blocks'))

    x$explVar = data.frame("comp" = factor(1:x$ncomp),
                           "percVar" = round(100*x$summary$R2X,4),
                           "cumPercVar" = round(100*x$summary$cumR2X,4))
    x$Bscores = x$BscoresX; x$Sscores = x$SscoresX

    if(plot) scorePlotmb(x, newObs = new)

    new = lapply(new, as.matrix)
    b_names = names(x$X)
    new = setNames(lapply(b_names, function(b) new[[b]][,colnames(x$X[[b]])]),b_names)

    centrado = x$scaling$center
    escalado = x$scaling$scale

    new = setNames(lapply(seq_along(new), function(b) {
      tmp = sweep(new[[b]], 2, centrado[[b]], "-")
      return(sweep(tmp, 2, escalado[[b]], "/"))
    }),b_names)

    if(any(unlist(lapply(new, function(x) any(rowSums(is.na(x))>0))))){
      return(stop('NAs are not allowed in new observations. Consider using PLS with block-scaling rather than MBPLS'))
    } else{
      block_coef = lapply(x$coefficients, function(b) apply(b[,,1:x$ncomp, drop = FALSE], c(1,2), sum))
      prediction = do.call('+',lapply(b_names, function(b) as.matrix(new[[b]])%*%block_coef[[b]]))
    }

    prediction = sweep(prediction, 2, x$scaling$scaleY, FUN = "*")
    plsprediction = sweep(prediction, 2, x$scaling$centerY, FUN = "+")
    prediction = as.data.frame(sub(".*_","", colnames(plsprediction)[apply(plsprediction, 1, which.max)]))
    colnames(prediction) = unique(sub("_.*","", colnames(plsprediction)[apply(plsprediction, 1, which.max)]))

  }

  return(list('prediction' = prediction, 'plsprediction' = plsprediction))

}

### PLSDA plots

#' Visualize PLSDA Results
#'
#' @description Generates diagnostic and results plots for objects created with the \code{plsda()} function.
#'
#' @param x A plsda object returned by the \code{plsda()} function.
#' @param type Character. The type of visualization to generate:
#' \itemize{
#'   \item \code{"ncomp"}: Line plots with the suggestion of the optimal number of components based on R2vsQ2 criterion, RMSE, F1-score, MCC, Etot and BER.
#'   \item \code{"R2vsQ2"}: Line plot of R2 and Q2 metrics per component.
#'   \item \code{"loadings"}: Joint visualization of X and Y loadings (for standard PLSDA).
#'   \item \code{"scoresX"}: Scatter plot of observation projections for X blocks.
#'   \item \code{"loadingsX"}: Scatter plot of variable contributions for X blocks.
#'   \item \code{"scoresY"}: Scatter plot of observation projections for Y block.
#'   \item \code{"loadingsY"}: Scatter plot of variable contributions for Y block.
#'   \item \code{"weights"}: Scatter plot of PLS weights.
#'   \item \code{"linearity"}: Scatter plot of X vs Y scores to check the inner relation.
#'   \item \code{"overfitting"}: Permutation test results to check for model overfitting.
#'   \item \code{"R2"}: Variance explained per variable/component.
#'   \item \code{"corr"}: Correlation circle plot (variables vs components).
#'   \item \code{"biplot"}: Simultaneous visualization of scores and loadings.
#'   \item \code{"coef"}: Bar plot of regression coefficients with Jack-knife confidence intervals.
#' }
#' @param comp Numeric vector. Components to plot (e.g., \code{c(1,2)}). By default, components 1 and 2 are used.
#' @param col Character. A predefined color palette name (e.g., "main", "oficial", ...) or a vector of custom colors.
#' @param colBy Variable used to color the points.
#' \itemize{
#'   \item For \strong{Scores}: A column name from the dataset or an external vector.
#'   \item For \strong{Loadings/Correlation}: Can be one of: \code{"contrib"} (to color by variable contribution)
#'   or \code{"cos2"} (to color by the quality of representation).
#' }
#' By default, no variable is used.
#' @param shape Numeric or character. The shape of the points (numeric) in the Score plots or 'arrow', 'point' for loading plots.
#' @param shapeBy Variable used to change point shapes. Can be a column name from the original dataset or an external vector/factor. Must be categorical.
#' @param ellipses Logical. If \code{TRUE}, draws 95\% confidence ellipses for groups defined in \code{colBy}.
#' @param selVars Numeric. The number of top variables to display in loading, correlation, or biplots, selected by their importance ("contrib" or "cos2"). Useful for decluttering plots with many variables.
#' @param labels Logical or Character vector. If \code{TRUE}, uses row names as labels. If Character, uses the provided vector as labels.
#' @param labelTop Numeric (0 to 1). Percentage of variables to label based on their importance.
#' @param repel Logical. If \code{TRUE}, uses \code{ggrepel} to prevent label overlap.
#' @param newObs Optional. A data frame (or list for MB-PLSDA) containing new observations to project onto the existing model.
#'
#' @return Depending on the model type:
#' \itemize{
#'   \item For \strong{Standard PLSDA}: Returns a \code{ggplot} object.
#'   \item For \strong{Multi-block PLSDA}:
#'     \itemize{
#'       \item Prints individual plots for each data block.
#'       \item Returns (or prints) the global "Super-score/Super-weights" plot.
#'     }
#' }
#' @export

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
      ncomp = if(comp < 3) if(decreasing) which.min(values) else which.max(values) else kneedle::kneedle(1:comp, values, decreasing = decreasing)[1]
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

  ### Validation

  if (type == 'overfitting'){

    ggp = validationPlot(x,
                         col = col)

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

      if (is.null(newObs)){
        colBy = x$input$Y
      } else{
        colBy = as.data.frame(c(x$input$Y[,1], as.factor(rep('newObs',nrow(Xnew)))))
        colnames(colBy) = 'colBy'
      }

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
      if (is.null(newObs)){
        colBy = x$input$Y
      } else{
        colBy = as.data.frame(c(x$input$Y[,1], as.factor(rep('newObs',nrow(Xnew)))))
        colnames(colBy) = 'colBy'
      }

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

      if (is.null(newObs)){
        colBy = x$input$Y
      } else{
        colBy = as.data.frame(c(x$input$Y[,1], as.factor(rep('newObs',nrow(Xnew)))))
        colnames(colBy) = 'colBy'
      }

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
      x$Bscores = x$BscoresX
      x$Bloadings = x$BloadingsX

      if (is.null(newObs)){
        colBy = x$input$Y
      } else{
        colBy = as.data.frame(c(x$input$Y[,1], as.factor(rep('newObs',nrow(Xnew)))))
        colnames(colBy) = 'colBy'
      }

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
      if (is.null(newObs)){
        colBy = x$input$Y
      } else{
        colBy = as.data.frame(c(x$input$Y[,1], as.factor(rep('newObs',nrow(Xnew)))))
        colnames(colBy) = 'colBy'
      }

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

    if (type == 'corr'){
      return(stop('In MBPLSDA corrPlot is not provided'))
    }

    ### biPlot

    if (type == 'biplot'){

      if(is.null(shape)) shape = c(18,'arrow')
      if(is.null(labels)) labels = c(FALSE, TRUE)
      if(is.null(ellipses)) ellipses = FALSE

      if (is.null(newObs)){
        colBy = x$input$Y
      } else{
        colBy = as.data.frame(c(x$input$Y[,1], as.factor(rep('newObs',nrow(Xnew)))))
        colnames(colBy) = 'colBy'
      }

      ggp = biPlotPLSmb(x,
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

#' Variable Selection for PLSDA models
#'
#' @description Performs variable selection using different metrics such as Jack-knife resampling,
#' Permutation tests, Variable Importance in Projection (VIP), or Significance Multivariate Correlation (sMC).
#'
#' @param x A pls object returned by the \code{plsda()} function.
#' @param type Selection method: "Jack" (default), "Perm", "VIP", or "sMC".
#' @param cvFolds Number of folds for cross-validation (used in Jack-knife). Default is 5.
#' @param rep Number of repetitions or permutations. Default is 10.
#' @param threshold Numerical value for selection (e.g., p-value < 0.05 or VIP > 1).
#'   If \code{NULL}, default statistical thresholds are applied (p-value = 0.05 and VIP = 1).
#' @param parallel Logical. If \code{TRUE}, computations are performed in parallel.
#'
#' @return Data frame or list containing the selected variables and their respective scores/p-values. Results are plotted
#' for Jack and Perm variable selection.
#'
#' @export

plsdaVarSel = function(x,
                     type = c('Jack','Perm','VIP','sMC')[1],
                     cvFolds = 5, rep = 10, threshold = NULL, parallel = FALSE){

  return(plsVarSel(x,type,cvFolds,rep,threshold,parallel))
}

#' PLSDA Outlier Detection
#'
#' @description Identifies atypical observations using Hotelling's T2 distance
#' and Residual Sum of Squares (RSS).
#'
#' @param x An object returned by the \code{plsda()} function.
#' @param ncomp Integer. Number of components to use for outlier calculation. If NULL, all components in the plsda object will be used. By default, NULL.
#' @param method Character. Detection method: "T2" (Hotellings-T2), "RSS" (Residuals Sum of Squares), or "both".
#' @param conf Numeric. Confidence level for the critical limit (e.g., 95, 99).
#'
#' @return A list of identified moderate and severe outliers.
#' @export

plsdaOutliers = function(x, ncomp = NULL,
                       method = c("T2", "RSS", "both")[3],
                       conf = 99){
  return(plsOutliers(x,ncomp,method,conf))

}

#' Outlier Contribution Analysis
#'
#' @description Analyzes which variables contribute most to an observation being
#' identified as an outlier.
#'
#' @param x An object returned by the \code{plsda()} function.
#' @param outliers Output from the \code{plsdaOutliers()} function.
#' @param labelSize Numeric. Font size for axis labels in the contribution plots.
#'
#' @return Bar plots of contributions and a list of numerical contribution values.
#' @export

plsdaOutlierContrib = function(x, outliers, labelSize = 1){
  return(plsOutlierContrib(x,outliers,labelSize ))
}
