# Multilevel approach

#' Multilevel Decomposition for Multivariate Analysis
#'
#' @description
#' Performs a multilevel decomposition of a data matrix into within-subject and
#' between-subject variations. This is particularly useful for crossover or
#' repeated measures designs where the individual "signature" needs to be
#' separated from the treatment effect.
#'
#' @param x A matrix, data.frame
#' @param y Matrix or data.frame containing the response variables.
#' @param design A named vector identifying the groups/subjects for decomposition.
#' @param ncomp  Number of components. If \code{NULL}, it is estimated based on design groups.
#' @param algo Character. Algorithm to use: "nipals", or "svd" (Single Value Decomposition).
#' @param scaling Character. Scaling method: "none", "center" (default), "standard", "softBlock", or "hardBlock".
#' @param scalingY Scaling method for Y: "none", "center" (default), or "standard".
#' @param cvFolds Number of cross-validation folds. Default is 10.
#' @param rep Number of cross-validation repetitions. Default is 10.
#' @param perm Number of permutations for model validation. Default is 20.
#' @param train Training set percentage (0 to 1). Default is 1 (all data).
#' @param alpha Significance level for Jack-knife confidence intervals. Default is 0.05.
#' @param parallel Logical. If \code{TRUE}, cross-validation runs in parallel.
#' @param method Character. The multivariate method to use: "pca", "pls", or "plsda".
#'
#' @return A list containing the decomposed matrices (Xm, Xw, Xb) and the
#' resulting multivariate model. On pca both the models of the within and between parts are returned.
#' @export

multilevel = function(x, y = NULL, design =NULL,
               ncomp = NULL,
               algo = c("nipals", "svd")[1],
               scaling = c("none", "center", "standard")[2],
               scalingY = c("none", "center", "standard")[2],
               cvFolds = 10, rep = 10, perm = 20,
               train = 1, alpha = 0.05, parallel = FALSE,
               method = c('pca','pls','plsda')[1])
{
  x = as.data.frame(x)
  if(is.null(ncomp)) ncomp = min(min(nrow(x)-1, ncol(x)), max(length(table(design))-1, ceiling(nrow(x)/length(table(design)))))

  if(!all(sapply(x, is.numeric))) {
    warning('Categorical variables automatically converted to dummy variables and default filters for NAs and CV applied.\n If the user does not want to use any default filtering, consider using Preparing function before executing pca to convert categorical variables into dummies.')
    P = Preparing(x, includeFactors = T, CVfilter = 0.01, excludeNA = 0.2)
    x = P$x
  } else {
    P = Preparing(x, CVfilter = 0.00001, excludeNA = 1) #Evitar nearZeroVariance variables pero no aplicar ningun filtro extra
    x = P$x
  }

  X = x

  if (scaling != 'none') {
    escalado = Scaling(x, scaling= 'center', blocks = NULL)
    x = escalado$x
  }

  #Calculate the offset
  Xm = X - x

  #Calculate the within subject variation (differences within groups of subjects) (el de interes) (total variation due to the treatment)
  design = as.data.frame(design)
  if(all(!is.na(suppressWarnings(as.numeric(rownames(design)))))) {
    rownames(design) = paste0('r', rownames(design))
  }
  if(all(order(rownames(design))!= order(rownames(X)))) rownames(design) = rownames(X)
  design = design[rownames(x),,drop=FALSE]
  colnames(design) = c('group')

  Xw = x

  for (g in unique(design$group)){
    idx <- g == design$group
    Xw[idx,] = sweep(Xw[idx,,drop=FALSE],2, colMeans(Xw[idx,,drop=FALSE]))
  }

  #Calculate the between subject variation (differences across groups of subjects) (within treatment variation)
  Xb = x - Xw

  if (method=='pca'){
    pca_xw = pca(Xw, ncomp = ncomp, algo ="nipals", scaling = scaling)
    pca_xb = pca(Xb, ncomp = ncomp, algo ="nipals", scaling = scaling)

    return(list('Within' = pca_xw, 'Between' = pca_xb, 'Xm' = Xm,'Xw'= Xw,'Xb' = Xb, 'input' = list('method' = method,'algo' = algo, 'scaling' = scaling, 'design' = design, 'model' = 'multi')))
  } else if (method=='pls'){
    Y = y

    if (scalingY != 'none') {
      escalado = Scaling(y, scaling= 'center', blocks = NULL)
      Y = escalado$x
    }

    # #Calculate the offset
    # Ym = Y - y (este y era el escalado)
    # Yw = y
    #
    # for (g in unique(design$group)){
    #   idx <- g == design$group
    #   Yw[idx,] = sweep(Yw[idx,,drop=FALSE],2, colMeans(Yw[idx,,drop=FALSE]))
    # }
    #
    # #Calculate the between subject variation (differences across groups of subjects)
    # Yb = y - Yw

    mypls = pls(Xw, Y, ncomp = ncomp, scaling = scaling, scalingY = scalingY, cvFolds = cvFolds, rep = rep, perm = perm, train = train, alpha = alpha, parallel = parallel)

    return(list('model' = mypls, 'Xm' = Xm, 'Xw' = Xw, 'Xb' = Xb, 'input' = list('method' = method,'algo' = algo, 'scaling' = scaling, 'scalingY'= scalingY, 'design' = design, 'model' = 'multi')))
  } else{

    Py = suppressWarnings(Preparing(y, includeFactors = T, CVfilter = 0.01, excludeNA = 0.2))
    Y = Py$x

    mypls = plsda(Xw, Y, ncomp = ncomp, scaling = scaling, scalingY = scalingY, cvFolds = cvFolds, rep = rep, perm = perm, train = train, alpha = alpha, parallel = parallel)

    return(list('model' = mypls, 'Xm' = Xm, 'Xw' = Xw, 'Xb' = Xb, 'input' = list('method' = method,'algo' = algo, 'scaling' = scaling, 'scalingY'= scalingY, 'design' = design, 'model' = 'multi')))
  }
}

#' Plotting Multilevel Multivariate Models
#'
#' @description
#' Visualizes the results of a multilevel analysis.
#'
#' @param x An object returned by the \code{multilevel} function.
#' @param type Character. The type of plot (e.g., "scores", "loadings", "biplot", "scree").
#' @param comp Numeric vector. Components to plot (e.g., \code{c(1,2)}). By default, components 1 and 2 are used.
#' @param col Character. A predefined color palette name (e.g., "main", "oficial", ...)  or a vector of custom colors.
#' @param colBy Variable used to color the points.
#' \itemize{
#'   \item For \strong{Scores}: A column name from the dataset or an external vector.
#'   \item For \strong{Loadings/Correlation}: Can be one of: \code{"contrib"} (to color by variable contribution)
#'   or \code{"cos2"} (to color by the quality of representation).
#'   \item For \strong{biPlot}: A column name from the dataset or an external vector to color the observations.
#' }
#' By default, no variable is used.
#' @param shape Numeric or character. The shape of the points (numeric) in the Score plots or 'arrow', 'point' for loading plots.
#' @param shapeBy Variable used to change point shapes. Can be a column name from the original dataset or an external vector/factor. Must be categorical.
#' @param ellipses Logical. If \code{TRUE}, draws 95\% confidence ellipses for groups defined in \code{colBy}.
#' @param selVars Numeric. The number of top variables to display in loading, correlation, weights or biplots, selected by their importance ("contrib" or "cos2"). Useful for decluttering plots with many variables.
#' @param labels Logical or Character vector. If \code{TRUE}, uses row names as labels. If Character, uses the provided vector as labels.
#' @param labelTop Numeric (0 to 1). Percentage of variables to label based on their importance (contribution or cos2).
#' @param repel Logical. If \code{TRUE}, uses \code{ggrepel} to prevent label overlap (recommended for loadings and biplots).
#' @param newObs Optional. A data frame (or list of data frames for MB-PCA) containing new observations to project onto the existing model.
#' @param newdesign A named vector identifying the groups/subjects for the new observations.
#'
#' @return A ggplot object (or patchwork of plots).
#' @export

multilevelPlot = function(x,
                     type = c("scree", "scores", "loadings","R2vsQ2", "scoresX",  "loadingsX", "scoresY", "loadingsY",
                              "weights", "linearity", "overfitting", "R2", "R2X","R2Y", "corr", "corrX", "corrY", "biplot", "coef", "ncomp"),
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
                     newObs = NULL, #data.frame
                     newdesign = NULL
) {

  if(x$input$method == 'pca'){
    x$Between$Xm = x$Within$Xm =x$Xm
    x$Between$multi = x$Within$multi =  'multi'
    x$Between$newdesign = x$Within$newdesign = newdesign

    if(type%in% c("ncomp", "R2vsQ2", "scoresX",  "loadingsX", "scoresY", "loadingsY", "weights", "linearity", "overfitting", "coef","corrY","R2Y")) return(stop('Plot not available for PCA models, please select one of scree, scores, loadings, R2, corr, biplot'))

    if(type=='R2X') type = "R2"
    if(type=='corrX') type = 'corr'

    x$input$design$group = as.character(x$input$design$group)
    colBy = if(is.null(colBy)) if(type=='scores' | type == 'biplot') x$input$design else colBy else colBy
    ellipses = if(is.null(ellipses)) ellipses = FALSE else ellipses

    ggp1 = pcaPlot(x$Between, type = type, comp = comp, col = col, colBy = colBy, shape = shape, shapeBy = shapeBy,
                   ellipses = ellipses, selVars = selVars, labels = labels, labelTop = labelTop, repel = repel, newObs = newObs)

    title = ggp1$labels$title

    ggp1 = ggp1 + labs(title = paste(title, "Between subject variation"))

    ggp2 = pcaPlot(x$Within, type = type, comp = comp, col = col, colBy = colBy, shape = shape, shapeBy = shapeBy,
                   ellipses = FALSE, selVars = selVars, labels = labels, labelTop = labelTop, repel = repel, newObs = newObs)
    ggp2 = ggp2 + labs(title = paste(title, "Within subject variation"))

    ggp = (ggp1 | ggp2) +
      patchwork::plot_layout(widths = c(1, 1)) +
      patchwork::plot_layout(guides = 'collect') &
      theme(legend.position = 'bottom')

  }

  if(x$input$method == 'pls'){

    if(type=='R2') type = "R2X"
    if(type=='corr') type = 'corrX'

    x$model$Xm = x$Xm
    x$model$multi = 'multi'
    x$model$newdesign = newdesign

    if(type%in% c("scree", "scores","loadings", "ncomp")) return(stop('Plot not available for PLS models, please select one of R2vsQ2, scoresX, loadingsX, scoresY, loadingsY, weights, linearity, overfitting, coef, R2X, R2Y, corrX, corrY, biplot'))
    x$input$design$group = as.character(x$input$design$group)
    colBy = if(is.null(colBy)) if(type=='scoresX' | type=='scoresY' | type == 'biplot') x$input$design
    ellipses = if(is.null(ellipses)) ellipses = FALSE

    ggp1 = plsPlot(x$model, type = type, comp = comp, col = col, colBy = colBy, shape = shape, shapeBy = shapeBy,
                   ellipses = ellipses, selVars = selVars, labels = labels, labelTop = labelTop, repel = repel, newObs = newObs)

    title = ggp1$labels$title

    ggp1 = ggp1 + labs(title = paste(title, "Between subject variation"))

    ggp2 = plsPlot(x$Within, type = type, comp = comp, col = col, colBy = colBy, shape = shape, shapeBy = shapeBy,
                   ellipses = FALSE, selVars = selVars, labels = labels, labelTop = labelTop, repel = repel, newObs = newObs)
    ggp2 = ggp2 + labs(title = paste(title, "Within subject variation"))

    ggp = (ggp1 | ggp2) +
      patchwork::plot_layout(widths = c(1, 1)) +
      patchwork::plot_layout(guides = 'collect') &
      theme(legend.position = 'bottom')
  }

  if(x$input$method == 'plsda'){

    if(type=='R2') type = "R2X"
    if(type=='corr') type = 'corrX'

    x$model$Xm = x$Xm
    x$model$multi = 'multi'
    x$model$newdesign = newdesign

    if(type%in% c("scree", "scores","loadings")) return(stop('Plot not available for PLSDA models, please select one of ncomp, R2vsQ2, scoresX, loadingsX, scoresY, loadingsY, weights, linearity, overfitting, coef, R2X, R2Y, corrX, corrY, biplot'))
    x$input$design$group = as.character(x$input$design$group)
    shapeBy = if(is.null(colBy)) if(type=='scoresX' | type == 'scoresY' | type == 'biplot') x$input$design

    ggp = plsdaPlot(x$model, type = type, comp = comp, col = col, colBy = colBy, shape = shape, shapeBy = shapeBy,
                   ellipses = ellipses, selVars = selVars, labels = labels, labelTop = labelTop, repel = repel, newObs = newObs)

    title = ggp$labels$title

    ggp = ggp + labs(title = paste(title, "Within subject variation"))

  }

  return(ggp)


}

