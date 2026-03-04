## PCA analysis

#' Principal Component Analysis (PCA) and Multi-block PCA
#'
#' @description Performs principal component analysis using various algorithms
#' (SVD, NIPALS) with support for multi-block data structures. Includes
#' automatic preprocessing and filtering of low-variance variables.
#'
#' @param x A matrix, data.frame, or a list of blocks (for MB-PCA).
#' @param ncomp Integer. Number of components to extract. If NULL, the maximum number of principal components allowed by the data is extracted (min(number_observations-1,number_variables)).
#' @param algo Character. Algorithm to use: "nipals", "svd" (Single Value Decomposition), or "mbpca" (multi-block PCA).
#' @param scaling Character. Scaling method: "none", "center", "standard" (default), "softBlock", or "hardBlock".
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{scores} or \code{Sscores}: Matrix of scores for each component (Super-scores in MB-PCA).
#'   \item \code{loadings} or \code{Sweights}: Matrix of variable loadings for each component (Super-weights in MB-PCA).
#'   \item \code{Bscores}: List of block scores (only for MB-PCA).
#'   \item \code{Bloadings}: List of block loadings (only for MB-PCA).
#'   \item \code{X}: The preprocessed and scaled data matrix (or list of matrices) used for model fitting.
#'   \item \code{scaling}: Data frame or list (only in MB-PCA) with centering and scaling parameters for each variable.
#'   \item \code{explVar}: Data frame with eigenvalues, percentage of variance, and cumulative variance per component.
#'   \item \code{BexplVar}: List of explained variance per block (only for MB-PCA).
#'   \item \code{PreproSummary}: Summary of samples, variables, excluded low-CV variables, and missing values.
#'   \item \code{ncomp}: Number of principal components extracted.
#'   \item \code{input}: A list with input parameters (\code{scalingType}, original \code{X}, and \code{algo}).
#' }
#' @export
#'

pca = function(x,
               ncomp = NULL,
               algo = c("nipals", "svd", "mbpca")[1],
               scaling = c("none", "center", "standard", "softBlock", "hardBlock")[3])
  {

  if(algo == 'mbpca'){

    if (!inherits(x, "list")) return(stop('To perform multiblock PCA, please provide the block information as a list where each element corresponds to one block'))
    if (length(unique(sapply(x, nrow)))!=1) return(stop('All blocks must have the same observations'))

    b_names = names(x)

    P = lapply(x, function(y){
      if(!all(sapply(y, is.numeric))) {
        cat('Warning: Categorical variables automatically converted to dummy variables and default filters for NAs and CV applied.\n If the user does not want to use any default filtering, consider using Preparing function before executing pca to convert categorical variables into dummies.\n')
        Preparing(y, includeFactors = T, CVfilter = 0.01, excludeNA = 0.2)
      } else {
        Preparing(y, CVfilter = 0.00001, excludeNA = 1) #Evitar nearZeroVariance variables pero no aplicar ningun filtro extra
      }
    })
    x = setNames(lapply(b_names, function(y) P[[y]]$x), b_names)

    X = x

    if(is.null(ncomp)) ncomp = min(sapply(X, function(x) min(nrow(x)-1, ncol(x))))

    if (scaling != 'none') {
      blocks = NULL
      escalado = setNames(lapply(b_names, function(y){
        if(scaling%in%c('softBlock', 'hardBlock')) blocks = rep(1, ncol(x[[y]]))
        return(Scaling(x[[y]], scaling = scaling, blocks = blocks))
      } ),b_names)

      x = setNames(lapply(b_names, function(y) escalado[[y]]$x), b_names)
      centrado = setNames(lapply(b_names, function(y) escalado[[y]]$centering), b_names)
      escalado = setNames(lapply(b_names, function(y) escalado[[y]]$scaling), b_names)
    } else{
      escalado = setNames(lapply(x,function(y) rep(1, times = ncol(y))),b_names)
      centrado = setNames(lapply(x,function(y) rep(0, times = ncol(y))),b_names)
      x = setNames(lapply(b_names, function(y) as.matrix(x[[y]])), b_names)
    }

    scal_param = setNames(lapply(b_names, function(y){
      center = centrado[[y]]
      scale = escalado[[y]]
      return(data.frame("center" = center , "scale" = scale))
    }),b_names)

    #Habria que discutir si metemos el caso de valores faltantes a MBPCA
    mypca = nipals_mbpca(x, ncomp = ncomp)

    totalVar = sum(sapply(x, function(y) sum(apply(y, 2, function(z) var(z,na.rm = TRUE)))))

    explVar = data.frame("comp" = factor(1:ncomp),
                         "eigenVal" = mypca$eigen,
                         "percVar" = round(100*mypca$eigen/totalVar,4),
                         "cumPercVar" = round(100*cumsum(mypca$eigen/totalVar),4))
    explVarBlock = setNames(lapply(b_names, function(y){
      expl = data.frame("comp" = factor(1:ncomp),
                        "percVar" = round(c(mypca$explvarB[[y]][1] ,diff(as.numeric(mypca$explvarB[[y]]))),4),
                        "cumPercVar" = t(round(mypca$explvarB[[y]],4)))
    }),b_names)

    desSummary = setNames(lapply(b_names, function(y){
      des = data.frame(c("Samples","Variables","Excluded Near-Zero Variables","Missing Values"),
                       c(nrow(x[[y]]),ncol(x[[y]]),length(which(P[[y]]$removed$Variables$Problem == 'Low CV')),
                         paste0(length(which(is.na(x[[y]]))), " (",round((length(which(is.na(x[[y]]))) / (nrow(x[[y]]) * ncol(x[[y]]))) * 100, 2), "%)")), stringsAsFactors = FALSE)
      colnames(des) = NULL
      return(des)
    }),b_names)

    return(list("Bscores" = mypca$block_scores,
                "Bloadings" = mypca$block_loadings,
                "Sscores" = mypca$super_scores,
                "Sweights" = mypca$super_loadings,
                "X" = x, ### con escalado y sin variables con baja variabilidad
                "scaling" = scal_param,
                "explVar" = explVar, ### Quizas habria que cambiar este nombre a summary para que cuadre con lo del PLS
                "BexplVar" = explVarBlock,
                "PreproSummary" = desSummary, ## completar si se desea dar + info de lo que se ha hecho
                "ncomp" = ncomp,
                "input" = list("scalingType" = scaling, "X" = X, 'algo' = algo, "model" = 'pca')))


  } else{

    if (!inherits(x, "list")){
      if(!all(sapply(x, is.numeric))) {
        cat('Warning: Categorical variables automatically converted to dummy variables and default filters for NAs and CV applied.\n If the user does not want to use any default filtering, consider using Preparing function before executing pca to convert categorical variables into dummies.\n')
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
          cat('Warning: Categorical variables automatically converted to dummy variables and default filters for NAs and CV applied.\n If the user does not want to use any default filtering, consider using Preparing function before executing pca to convert categorical variables into dummies.\n')
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

    X = x

    if(!inherits(x, "list") && scaling %in% c('softBlock','hardBlock')) blocks = rep(1, times = ncol(x)) else blocks = NULL

    if (inherits(x, "list")) {
      if (!scaling %in% c("softBlock", "hardBlock")) cat('Warning: When blocks are considered we recommend using either softBlock or hardBlock scaling\n')
      if (length(unique(sapply(x, nrow)))!=1) return(stop('All blocks must have the same observations'))
      blocks = rep(seq_along(x), times = sapply(x, ncol) )
      x = do.call(cbind, x)
    }

    if (is.null(ncomp)) ncomp = min(ncol(x), nrow(x)-1)

    if (scaling != 'none') {
      escalado = Scaling(x, scaling= scaling, blocks = blocks)
      x = escalado$x
      centrado = escalado$centering
      escalado = escalado$scaling
    } else{
      escalado = rep(1, times = ncol(x))
      centrado = rep(0, times = ncol(x))
    }

    totalVar = sum(apply(x, 2, function(x) var(x,na.rm = TRUE)))

    if(any(is.na(x)) && algo == 'svd') return(stop('Missing values detected. Please remove or impute them, or use the nipals algorithm instead'))

    mypca = if(algo == 'svd') svd_pca(x, ncomp = ncomp) else nipals_pca(x, ncomp = ncomp)

    explVar = data.frame("comp" = factor(1:length(mypca$eigen)),
                         "eigenVal" = mypca$eigen,
                         "percVar" = round(100*mypca$eigen/totalVar,4),
                         "cumPercVar" = round(100*cumsum(mypca$eigen/totalVar),4))

    return(list("scores" = mypca$scores,
                "loadings" = mypca$loadings,
                "X" = x, ### con escalado y sin variables con baja variabilidad
                "scaling" = data.frame("center" = centrado , "scale" = escalado),
                "explVar" = explVar,
                "PreproSummary" = desSummary, ## completar si se desea dar + info de lo que se ha hecho
                "ncomp" = ncomp,
                "input" = list("scalingType" = scaling, "X" = X, 'algo' = algo, "model" = 'pca')))
  }

}


### PCA plots

#' Visualize PCA Results
#'
#' @description Generates diagnostic and results plots for objects created with the \code{pca()} function.
#'
#' @param x An object returned by the \code{pca()} function.
#' @param type Character. The type of visualization to generate:
#' \itemize{
#'   \item \code{"scree"}: Bar plot of explained variance per component.
#'   \item \code{"scores"}: Scatter plot of observation projections.
#'   \item \code{"loadings"}: Scatter plot of variable contributions.
#'   \item \code{"corr"}: Correlation circle plot (variables vs components).
#'   \item \code{"biplot"}: Simultaneous visualization of scores and loadings.
#'   \item \code{"R2"}: Variance explained per variable/component.
#' }
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
#' @param size Numeric. The size of the points in the Score plots.
#' @param ellipses Logical. If \code{TRUE}, draws 95\% confidence ellipses for groups defined in \code{colBy}.
#' @param selVars Numeric. The number of top variables to display in loading, correlation, or biplots, selected by their importance ("contrib" or "cos2"). Useful for decluttering plots with many variables.
#' @param labels Logical or Character vector. If \code{TRUE}, uses row names as labels. If Character, uses the provided vector as labels.
#' @param labelTop Numeric (0 to 1). Percentage of variables/observations to label based on their importance (contribution or cos2).
#' @param repel Logical. If \code{TRUE}, uses \code{ggrepel} to prevent label overlap (recommended for loadings and biplots).
#' @param newObs Optional. A data frame (or list of data frames for MB-PCA) containing new observations to project onto the existing model.
#'
#' @return Depending on the model type:
#' \itemize{
#'   \item For \strong{Standard PCA}: Returns a \code{ggplot} object.
#'   \item For \strong{Multi-block PCA}:
#'     \itemize{
#'       \item Prints individual plots for each data block.
#'       \item Returns (or prints) the global "Super-score/Super-weights" plot.
#'     }
#' }
#' @export

pcaPlot = function(x,
                   type = c("scree", "scores", "loadings", "corr", "biplot", "R2"),
                   comp = NULL,
                   col = c('main', 'complete', 'cblindfriendly','sunshine','hot','warm','grass','oficial')[1],
                   colBy = NULL,
                   shape = NULL,
                   shapeBy = NULL,
                   size = NULL,
                   ellipses = NULL,  # only for scores
                   selVars = NULL,
                   labels = NULL,
                   labelTop = NULL, #percentage expressed in 0-1 range
                   repel = TRUE,
                   newObs = NULL #data.frame or list of data.frames in case of mbpca
                   ) {

  if(!(type %in% c("scree", "scores", "loadings", "corr", "biplot", "R2"))) return(stop('Please use one of: scree, scores, loadings, corr, biplot, R2.'))

  if (type == "scree") {

    if(col=='main') col = "skyblue3"
    if(is.null(comp)) comp = 1:min(nrow(x$explVar),if(x$input$algo!='mbpca') 20 else 10)

    df = x$explVar[comp,c("comp", "percVar")]
    etiq = paste0(round(df$percVar, 1), "%")

    ggp = ggplot(df, aes(comp, percVar, group = 1))
    ggp = ggp + geom_bar(stat = "identity", fill = col, color = col)
    ggp = ggp + geom_line(color = "grey50", linewidth = 0.5)
    ggp = ggp + labs(title = "Scree plot", x = "Comp", y = "% Explained Variance")
    ggp = ggp + geom_text(label = etiq, vjust = -0.3, hjust = 0.4)

    media = 100 * (1/nrow(x$explVar))
    ggp = ggp + geom_hline(yintercept=media, linetype=2, color="red3")

    if(x$input$algo=='mbpca'){

      explvar = unlist(lapply(x$BexplVar, function(x) x$percVar))
      ncomp = as.factor(rep(1:x$ncomp, length(x$X)))
      b_names = as.factor(rep(names(x$X),each = x$ncomp))
      df_block = data.frame(explvar, ncomp, b_names)
      etiq = paste0(round(as.numeric(df_block$explvar), 1), "%")

      ggp_blocks = ggplot(df_block, aes(ncomp,explvar))+
        geom_bar(stat = "identity", color= col, fill= col)+
        geom_line(aes(group = b_names), color = "grey50", size = 0.5) +
        geom_text(aes(label = etiq), vjust = -0.3, hjust = 0.4) +
        geom_hline(yintercept = media, linetype = 2, color = "red3") +
        facet_wrap(~b_names, ncol=length(comp), strip.position = 'top')+
        labs(
          title = "Explained Variance per Block",
          x = "Comp",
          y = "% Explained Variance"
        )

      ggp = ggp + labs(title = "Globally Explained Variance", x = "Comp", y = "% Explained Variance")

      ggp = ggp + ggp_blocks + patchwork::plot_layout(ncol = 2) +
        patchwork::plot_annotation(title = "Scree plot", theme = theme(plot.title = element_text(size = 16, face = "bold")))

    }

  }

  ### Scores

  if (type == "scores") {

    if(is.null(shape)) shape = 18
    if(is.null(size)) size = 3 # redundant I think as scorePlot child function defines "size = 3" as default value
    if(is.null(labels)) labels = FALSE
    if(is.null(ellipses)) ellipses = TRUE
    if(!is.null(colBy) & length(colBy)>1) colBy = as.data.frame(colBy)
    if(!is.null(shapeBy) & length(shapeBy)>1) shapeBy = as.data.frame(shapeBy)

    if(x$input$algo!='mbpca'){
      ggp = scorePlot(x,
                      comp = comp,
                      col = col,
                      colBy = colBy,
                      shape = shape,
                      shapeBy = shapeBy,
                      size = size,
                      ellipses = ellipses,
                      labels = labels,
                      labelTop = labelTop,
                      repel = repel,
                      newObs = newObs)
    } else{
      ggp = scorePlotmb(x,
                        comp = comp,
                        col = col,
                        colBy = colBy,
                        shape = shape,
                        shapeBy = shapeBy,
                        ellipses = ellipses,
                        labels = labels,
                        labelTop = labelTop,
                        repel = repel,
                        newObs = newObs)
    }
  }

  ### Loadings

  if (type == "loadings") {

    if(is.null(shape)) shape = 'arrow'
    if(is.null(labels)) labels = TRUE

    if(x$input$algo!='mbpca'){
      ggp = loadingPlot(x,
                comp = comp,
                col = col,
                colBy = colBy,
                shape = shape,
                selVars = selVars,
                labels = labels,
                labelTop = labelTop,
                repel = repel,
                newObs = newObs)
    } else{
      if(!is.null(newObs)) return(stop('In MBPCA projection of new variables is not provided.'))
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

  }

  ### Correlation

  if (type == "corr") {

    if(is.null(shape)) shape = 'arrow'
    if(is.null(labels)) labels = TRUE

    if(x$input$algo!='mbpca'){
      ggp = corrPlot(x,
                     comp = comp,
                     col = col,
                     colBy = colBy,
                     shape = shape,
                     selVars = selVars,
                     labels = labels,
                     labelTop = labelTop,
                     repel = repel,
                     newObs = newObs)
    } else return(stop('In MBPCA corrPlot is not provided'))

  }

  ### Biplot

  if (type == "biplot"){

    if(is.null(shape)) shape = c(18,'arrow')
    if(is.null(labels)) labels = c(FALSE, TRUE)
    if(is.null(ellipses)) ellipses = FALSE

    if(!is.null(colBy) & length(colBy)>1) colBy = as.data.frame(colBy)
    if(!is.null(shapeBy) & length(shapeBy)>1) shapeBy = as.data.frame(shapeBy)

    if(x$input$algo!='mbpca'){
      ggp = biPlot(x,
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
    } else{
      ggp = biPlotmb(x,
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

  }

  ### R2

  if (type == 'R2'){

    if(x$input$algo!='mbpca'){
      ggp = R2varcomp(x, col)
    } else ggp = R2varcompmb(x, col)

  }

  return(ggp)
}

### PCA outliers

#' PCA Outlier Detection
#'
#' @description Identifies atypical observations using Hotelling's T2 distance
#' and Residual Sum of Squares (RSS).
#'
#' @param x An object returned by the \code{pca()} function.
#' @param ncomp Integer. Number of components to use for outlier calculation. If NULL, all components in the pca object will be used. By default, NULL.
#' @param method Character. Detection method: "T2" (Hotellings-T2), "RSS" (Residuals Sum of Squares), or "both".
#' @param conf Numeric. Confidence level for the critical limit (e.g., 95, 99).
#'
#' @return A list of identified moderate and severe outliers.
#' @export

pcaOutliers = function(x, ncomp = NULL,
                       method = c("T2", "RSS", "both")[3],
                       conf = 99) {

  if (is.null(ncomp)) ncomp = x$ncomp

  if (method == "T2") {

    out = T2plot(x = x, K = ncomp, conf = conf)
    return(out)

  }

  if (method == "RSS") {

    out = RSSplot(x = x, K = ncomp, conf = conf)
    return(out)

  }

  if (method == "both") {

    par(mfrow = c(1,2))
    outS = T2plot(x = x, K = ncomp, conf = conf)
    outM = RSSplot(x = x, K = ncomp, conf = conf)
    par(mfrow = c(1,1))

    return(c(outS, outM))

  }

}

#' Outlier Contribution Analysis
#'
#' @description Analyzes which variables contribute most to an observation being
#' identified as an outlier.
#'
#' @param x An object returned by the \code{pca()} function.
#' @param outliers Output from the \code{pcaOutliers()} function.
#' @param labelSize Numeric. Font size for axis labels in the contribution plots.
#' @param specificObs Character or vector of characters. Names of the specific observations to plot.
#'
#' @return Bar plots of contributions and a list of numerical contribution values.
#' @export

outlierContrib = function(x, outliers, labelSize = 1, specificObs = NULL) {

  # Severe outliers (T2)
  if (!is.null(specificObs)) {
    severos = specificObs
    num = length(intersect(severos, outliers$SevereOutliers))
  } else {
    severos = outliers$SevereOutliers
    num = length(severos)
  }

  if (num == 0) {
    cat("There are no severe outliers in the data or they were not computed.\n")
    miscontrT2 = NULL
  } else {
    cat (paste0("There are ", num, " severe outliers in the data provided.\n"))
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
              las=2, cex.names = labelSize, border = NA, col = "lightcoral",
              main = paste0("Severe outlier: ", i))
    }

  }


  # Moderate outliers (RSS)
  if (!is.null(specificObs)) {
    modera = specificObs
    num = length(intersect(modera, outliers$ModerateOutliers))
  } else {
    modera = outliers$ModerateOutliers
    num = length(modera)
  }

  if (num == 0) {
    cat("There are no moderate outliers in the data or they were not computed.\n")
    miscontrRSS = NULL
  } else {
    cat (paste0("There are ", num, " moderate outliers in the data provided.\n"))
    miscontrRSS = ContriSCR(E = outliers$E, SCR = outliers$RSS)
    miscontrRSS = miscontrRSS[modera,,drop = FALSE]
    if (num > 10) {
      cat (paste0("Only the 10 most moderate outliers will be plotted.\n"))
      num = 10
      modera = rownames(sort(outliers$RSS, decreasing = TRUE)[1:10])
    }
    for (i in modera) {
      barplot(sort(miscontrRSS[i,], decreasing = TRUE),
              las=2, cex.names = labelSize, border = NA, col = "skyblue",
              main = paste0("Moderate outlier: ", i))
    }
  }

  return(list(contribT2 = if(!is.null(miscontrT2)) t(miscontrT2) else NULL, contribRSS = miscontrRSS))

}

### NIPALS pseudo-imputation

#' Missing Value Imputation via NIPALS
#'
#' @description Reconstructs missing values (NAs) in the original matrix using the underlying structure captured by the PCA/PLS model.
#'
#' @param x An object returned by the \code{pca()} function.
#' @param scaled Logical. If TRUE, returns data in the model's scale (centered/scaled).
#' If FALSE, returns data in the original units. By default, FALSE.
#'
#' @return A matrix with original values and imputed NAs.
#' @export

imputeNipals = function(x, scaled = FALSE){
  X = x$X
  if(x$input$model == 'pca'){
    Xest = tcrossprod(x$scores, x$loadings)
  } else{
    Xest = tcrossprod(x$scoresX, x$loadingsX)
  }
  X[is.na(X)] = Xest[is.na(X)]

  if(!scaled){
    X = sweep(X, 2, x$scaling$scale, FUN = "*")
    X = sweep(X, 2, x$scaling$center, FUN = "+")
  }
  return(X)

}

