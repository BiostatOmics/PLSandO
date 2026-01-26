

# Multilevel approach

# We will want to separate the sum of squares of the offset, the between and the within variance

#First the data has to be mean centered. The data becomes a sum of the mean-centered data and the offset

multilevel = function(x, y = NULL, design =NULL,
               ncomp = NULL,
               algo = c("nipals", "svd")[1],
               scaling = c("none", "center", "standard")[2],
               scalingY = c("none", "center", "standard")[3],
               method = c('pca','pls','plsda')[1])
{

  if(!all(sapply(x, is.numeric))) {
    warning('Categorical variables automatically converted to dummy variables and default filters for NAs and CV applied.\n If the user does not want to use any default filtering, consider using Preparing function before executing pca to convert categorical variables into dummies.')
    P = Preparing(x, includeFactors = T, CVfilter = 0.01, excludeNA = 0.2)
    x = P$x
  } else {
    P = Preparing(x, CVfilter = 0.00001, excludeNA = 1) #Evitar nearZeroVariance variables pero no aplicar ningun filtro extra
    x = P$x
  }

  X = x

  if (is.null(ncomp)) ncomp = min(ncol(x), nrow(x)-1)

  if (scaling != 'none') {
    escalado = Scaling(x, scaling= 'center', blocks = NULL)
    x = escalado$x
  }

  #Calculate the offset

  Xm = X - x

  #Calculate the within subject variation (differences within groups of subjects) (el de interes) (total variation due to the treatment)

  design = as.data.frame(design)
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

  #To be done: Comprobar que este sea el sitio correcto de escalar los datos y que no tenga que ser al principio

  if (scaling == 'standard'){
    escalado = Scaling(Xw, scaling= scaling, blocks = NULL)
    Xw = escalado$x
    centrado_xw = escalado$centering
    escalado_xw = escalado$scaling

    escalado = Scaling(Xb, scaling= scaling, blocks = NULL)
    Xb = escalado$x
    centrado_xb = escalado$centering
    escalado_xb = escalado$scaling
  }

  if (method=='pca'){
    pca_xw = pca(Xw, ncomp = ncomp, algo ="nipals", scaling = "none")
    pca_xb = pca(Xb, ncomp = ncomp, algo ="nipals", scaling = "none")

    return(list('Within' = pca_xw, 'Between' = pca_xb, 'arguments' = list('method' = method)))
  } else if (method=='pls'){
    Y = y

    if (scalingY != 'none') {
      escalado = Scaling(y, scaling= 'center', blocks = NULL)
      y = escalado$x
    }

    #Calculate the offset

    Ym = Y - y

    Yw = y

    for (g in unique(design$group)){
      idx <- g == design$group
      Yw[idx,] = sweep(Yw[idx,,drop=FALSE],2, colMeans(Yw[idx,,drop=FALSE]))
    }

    #Calculate the between subject variation (differences across groups of subjects)
    Yb = y - Yw

    if (scalingY == 'standard'){
      escalado = Scaling(Yw, scaling= scaling, blocks = NULL)
      Yw = escalado$x
      centrado_yw = escalado$centering
      escalado_yw = escalado$scaling

      escalado = Scaling(Yb, scaling= scaling, blocks = NULL)
      Yb = escalado$x
      centrado_yb = escalado$centering
      escalado_yb = escalado$scaling
    }

    mypls = pls(Xw, Yw, ncomp = ncomp, scaling = "none", scalingY = 'none' )

    return(mypls)
  } else{

    #PLS-DA


  }

}


multilevelPlot = function(x,
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

  if(x$arguments$method){



  }


}
