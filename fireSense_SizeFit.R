# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense_SizeFit",
  description = "Fit statistical models describing the empirical distribution of fire sizes.",
  keywords = c("fire size distribution", "tapered Pareto", "optimization", "fireSense", "statistical model"),
  authors=c(person("Jean", "Marchal", email="jean.d.marchal@gmail.com", role=c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("1.2.0.9000"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_SizeFit.Rmd"),
  reqdPkgs = list("DEoptimR", "magrittr", "numDeriv", "PtProcess"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", default, min, max, "parameter description")),
    defineParameter("formula", "list", NULL, desc = 'a (named) list with two components called "beta" and "theta" of class "formula" : a symbolic description of the model to be fitted (beta and theta parameters of the tapered Pareto).'),
    defineParameter("a", "numeric", NULL, desc = 'range parameter a of the tapered Pareto. The random variable x take values on the interval a <= x < Inf. Values outside of this range are removed with a warning.'),
    defineParameter("link", "list", list(beta = "log", theta = "identity"), desc = 'a (named) list of two link functions for the beta parameter of the tapered Pareto. see ?family. It should contain two components named "beta" and "theta".'),
    defineParameter(name = "start", class = "numeric", default = NULL,
                    desc = 'optional (named) list of starting values for the parameters to be estimated. Should contain two components named "beta" and "theta".'),
    defineParameter(name = "lb", class = "numeric", default = NULL, desc = 'optional lower bounds for the parameters to be estimated. Should contain two components named "beta" and "theta".'),
    defineParameter(name = "ub", class = "numeric", default = NULL, desc = 'optional upper bounds for the parameters to be estimated. Should contain two components named "beta" and "theta".'),
    defineParameter(name = "nlminb.control", class = "numeric", default = list(iter.max = 5000L, eval.max=5000L),
                    desc = "optional list of control parameters to be passed to the nlminb optmizer. See ?nlminb"),
    defineParameter(name = "trace", class = "numeric", default = 0,
                    desc = "non-negative integer. If > 0, tracing information on the progress of the optimization is produced every trace iteration. Defaults to 0 which indicates no trace information is to be printed."),
    defineParameter(name = "data", class = "character", default = NULL,
      desc = "optional. A character vector indicating the names of objects present in the sim environment, in which
              to look for variables with which to predict. Objects can be data.frames. If omitted, or if variables
              are not found in the data objects, variables are searched in the sim environment."),
    defineParameter(name = "mapping", class = "character", default = NULL,
      desc = "optional. Named character vector to map variable names in the formula to those in the data objects.
              Names of unmapped variables are used directly to look for variables in data objects or in the sim environment."),
    defineParameter(name = "initialRunTime", class = "numeric", default = NA, desc = "optional. Simulation time at which to start this module. If omitted, start at start(sim)."),
    defineParameter(name = "intervalRunModule", class = "numeric", default = NA, desc = "optional. Interval in simulation time units between two module runs.")
  ),
  inputObjects = data.frame(objectName="dataFireSense_SizeFit",
                            objectClass="data.frame",
                            sourceURL="",
                            other=NA_character_,
                            stringsAsFactors=FALSE),
  outputObjects = data.frame(
    objectName = "fireSense_SizeFit",
    objectClass = "fireSense_SizeFit",
    other = NA_character_,
    stringsAsFactors = FALSE
  )
))

## Toolbox: set of functions used internally by the module
  ## Compute the order of magnitude
  oom <- function(x){
    10^(ceiling(log10(abs(x))))
  }

  ## Function to pass to the optimizer
  objFun <- function(params, scalMx, mmBeta, mmTheta, ntBeta, nt, lnBeta, lnTheta, y, envData){

    ## Parameters scaling
      params <- drop(params %*% scalMx)

    beta <- drop(mmBeta %*% params[1:ntBeta])
    theta <- drop(mmTheta %*% params[(ntBeta + 1L):nt])

    if(length(beta) == 1L) beta <- rep_len(beta, length(theta)) ## Recycled if needed
    if(length(theta) == 1L) theta <- rep_len(theta, length(beta)) ## Recycled if needed

    ## link implementation
      beta <- lnBeta$linkinv(beta)
      theta <- lnTheta$linkinv(theta)

    if(any(beta <= 0L) || anyNA(beta) || any(is.infinite(beta)) || any(theta <= 0L) || anyNA(theta) || any(is.infinite(theta))){
      return(1e20)
    } else {
      return(eval(nll))
    }
  }

## event types
#   - type `init` is required for initialiazation

doEvent.fireSense_SizeFit = function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    ### check for more detailed object dependencies:
    ### (use `checkObject` or similar)

    # do stuff for this event
    sim <- sim$fireSense_SizeFitInit(sim)

    # schedule future event(s)
    #sim <- scheduleEvent(sim, p(sim)$.plotInitialTime, "fireSense_SizeFit", "plot")
    #sim <- scheduleEvent(sim, p(sim)$.saveInitialTime, "fireSense_SizeFit", "save")
  } else if (eventType == "run") {

    sim <- sim$fireSense_SizeFitRun(sim)

  } else if (eventType == "save") {
    # ! ----- EDIT BELOW ----- ! #
    # do stuff for this event

    # e.g., call your custom functions/methods here
    # you can define your own methods below this `doEvent` function

    # schedule future event(s)

    # e.g.,
    # sim <- scheduleEvent(sim, time(sim) + increment, "fireSense_SizeFit", "save")

    # ! ----- STOP EDITING ----- ! #
  } else {
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  return(invisible(sim))
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

fireSense_SizeFitInit <- function(sim) {
  
  sim <- scheduleEvent(sim, eventTime = if (is.na(p(sim)$initialRunTime)) start(sim) else p(sim)$initialRunTime, "fireSense_SizeFit", "run")
  sim

}

fireSense_SizeFitRun <- function(sim) {

  envData <- new.env(parent = envir(sim))
  on.exit(rm(envData))
  list2env(as.list(envir(sim)), envir = envData)

  ## Handle data arg
  if (!is.null(p(sim)$data))
    lapply(p(sim)$data, function(x, envData) if (is.list(sim[[x]])) list2env(sim[[x]], envir = envData), envData = envData)

  ## Check formula for beta
  if (is.empty.model(p(sim)$formula$beta))
    stop("fireSense_SizeFit> The specified formula (beta) is empty")
  
  ## Check formula for theta
  if (is.empty.model(p(sim)$formula$theta))
    stop("fireSense_SizeFit> The specified formula (theta) is empty")

  termsBeta <- terms.formula(formulaBeta <- p(sim)$formula$beta)
  termsTheta <- terms.formula(formulaTheta <- p(sim)$formula$theta)
    
  if (attr(termsBeta, "response"))
    y <- yBeta <- as.character(attr(termsBeta, "variables")[[2L]])
  else
    stop("fireSense_SizeFit> Incomplete formula (beta), the LHS is missing.")
  
  if (attr(termsTheta, "response"))
    yTheta <- as.character(attr(termsTheta, "variables")[[2L]])
  else
    stop("fireSense_SizeFit> Incomplete formula (theta), the LHS is missing.")
  
  if (!identical(yBeta, yTheta))
    stop("fireSense_SizeFit> the response variables for beta and theta must be identical.")

  
  ## Mapping variables names to data
  if (!is.null(p(sim)$mapping)) {
    
    for (i in 1:length(p(sim)$mapping)) {
      
      attr(termsBeta, "term.labels") <- gsub(pattern = names(p(sim)$mapping[i]),
                                             replacement = p(sim)$mapping[i], x = attr(termsBeta, "term.labels"))
      attr(termsTheta, "term.labels") <- gsub(pattern = names(p(sim)$mapping[i]),
                                              replacement = p(sim)$mapping[i], x = attr(termsTheta, "term.labels"))
    }
    
  }
  
  formulaBeta <- reformulate(attr(termsBeta, "term.labels"), y, attr(termsBeta, "intercept"))
  formulaTheta <- reformulate(attr(termsTheta, "term.labels"), y, attr(termsTheta, "intercept"))
  
    
  ## Coerce lnBeta to a link-glm object
  if (is.character(p(sim)$link$beta)) {
    lnBeta <- make.link(p(sim)$link$beta)
  } else if (is(p(sim)$link$beta, "link-glm")) {
    ## Do nothing
  } else {
    ## Try to coerce to link-glm class
    lnBeta <- make.link(p(sim)$link$beta)
  }

  ## Coerce lnTheta to a link-glm object
  if (is.character(p(sim)$link$theta)) {
    lnTheta <- make.link(p(sim)$link$theta)
  } else if (is(p(sim)$link$theta, "link-glm")) {
    ## Do nothing
  } else {
    lnTheta <- make.link(p(sim)$link$theta)
  }

  ## No tracing if trace < 0
  trace <- if(p(sim)$trace < 0) 0 else p(sim)$trace

  ## If there are rows in the dataset where y < a, remove them
  rm <- envData[[y]] < p(sim)$a

  if (all(rm)){
    
    stop("fireSense_SizeFit> All rows contain values outside of the range a <= x < Inf.")
    
  } else if (any(rm)) {
    
    lapply(unique(c(all.vars(termsBeta), all.vars(termsTheta))), function(x) assign(x = x, value = envData[[x]][!rm], envir = envData))
    warning(paste("fireSense_SizeFit> Removed", sum(rm), "rows containing values outside of the range a <= x < Inf."), immediate. = TRUE)
    
  }

  ## Number of terms
  ntBeta <- length(labels(termsBeta)) + attr(termsBeta, "intercept")
  ntTheta <- length(labels(termsTheta)) + attr(termsTheta, "intercept")
  nt <- ntBeta + ntTheta

  ## Define the scaling matrices. This is used later in the optimization process
  ## to rescale parameter values between 0 and 1, i.e. put all variables on the same scale.
  scalMx <- matrix(0L, nt, nt)
  diag(scalMx) <- 1L

  ## Design matrices
  mmBeta <- model.matrix(termsBeta, envData)
  mmTheta <- model.matrix(termsTheta, envData)

  ## Define parameter bounds automatically if they are not supplied by user
  ## First defined the bounds for DEoptim, the first optimizer
    ## Beta

      switch(lnBeta$name,
             log = {
               DEoptimUpperBound <-
                 if (is.null(p(sim)$ub$beta)) {
                   ## Automatically estimate an upper boundary for each parameter
                   (lm(update(formulaBeta, log(.) ~ .),
                       y = FALSE,
                       model = FALSE,
                       data = envData) %>%
                      coef %>%
                      abs) * 1.1
                 } else {
                   ## User-defined bounds
                   rep_len(p(sim)$ub$beta, ntBeta) ## Recycled if necessary
                 }

               DEoptimLowerBound <-
                 if (is.null(p(sim)$lb$beta)) {
                   ## Automatically estimate a lower boundary for each parameter
                   -DEoptimUpperBound
                 } else {
                   ## User-defined bounds
                   rep_len(p(sim)$lb$beta, ntBeta) ## Recycled if necessary
                 }
             }, identity = {
               DEoptimUpperBound <-
                 if (is.null(p(sim)$ub$beta)) {
                   ## Automatically estimate an upper boundary for each parameter
                   (lm(formulaBeta,
                       y = FALSE,
                       model = FALSE,
                       data = envData) %>%
                      coef %>%
                      abs) * 1.1
                 } else {
                   ## User-defined bounds
                   rep_len(p(sim)$ub$beta, ntBeta) ## Recycled if necessary
                 }

               DEoptimLowerBound <-
                 if(is.null(p(sim)$lb$beta)){
                   ## Automatically estimate an lower boundary for each parameter
                   rep_len(1e-30, ntBeta) ## Enforce non-negativity, recycled if necessary
                 } else {
                   ## User-defined bounds
                   rep_len(p(sim)$lb$beta, ntBeta) ## Recycled if necessary
                 }
             }, stop(paste("fireSense_SizeFit> Link function", p(sim)$link$beta, "(beta) is not supported.")))

  ## Theta
  switch(lnTheta$name,
         log = {
           DEoptimUpperBound <- c(DEoptimUpperBound,
                                  if (is.null(p(sim)$ub$theta)) {
                                    ## Automatically estimate an upper boundary for each parameter
                                    (lm(update(formulaTheta, log(.) ~ .),
                                        y = FALSE,
                                        model = FALSE,
                                        data = envData) %>%
                                       coef %>%
                                       abs) * 1.1
                                  } else {
                                    ## User-defined bounds
                                    rep_len(p(sim)$ub$theta, ntTheta) ## Recycled if necessary
                                  })

           DEoptimLowerBound <- c(DEoptimLowerBound,
                                  if(is.null(p(sim)$lb$theta)){
                                    ## Automatically estimate a lower boundary for each parameter
                                    -DEoptimUpperBound[(ntBeta + 1L):nt]
                                  } else {
                                    ## User-defined bounds
                                    rep_len(p(sim)$lb$theta, ntTheta) ## Recycled if necessary
                                  })
         }, identity = {
           DEoptimUpperBound <- c(DEoptimUpperBound,
                                  if (is.null(p(sim)$ub$theta)) {
                                    ## Automatically estimate an upper boundary for each parameter
                                    (lm(formulaTheta,
                                        y = FALSE,
                                        model = FALSE,
                                        data = envData) %>%
                                       coef %>%
                                       abs) * 1.1
                                  } else {
                                    ## User-defined bounds
                                    rep_len(p(sim)$ub$theta, ntTheta) ## Recycled if necessary
                                  })

           DEoptimLowerBound <- c(DEoptimLowerBound,
                                  if (is.null(p(sim)$lb$theta)) {
                                    ## Automatically estimate an lower boundary for each parameter
                                    rep_len(1e-30, ntTheta) ## Enforce non-negativity, recycled if necessary
                                  } else {
                                    ## User-defined bounds
                                    rep_len(p(sim)$lb$theta, ntTheta) ## Recycled if necessary
                                  })
         }, stop(paste("fireSense_SizeFit> Link function", p(sim)$link$theta, "(theta) is not supported.")))

  ## Then, define lower and upper bounds for the second optimizer (nlminb)
  ## Beta
    nlminbUpperBound <-
      if(is.null(p(sim)$ub$beta)){
        rep_len(Inf, ntBeta)
      } else {
        ## User-defined lower bounds for parameters to be estimated
        DEoptimUpperBound[1:ntBeta]
      }
  
    nlminbLowerBound <-
      if(is.null(p(sim)$lb$beta)){
        switch(lnBeta$name,
               log = -Inf,            ## log-link, default: -Inf for terms and 0 for breakpoints/knots
               identity = 1e-30) %>%  ## identity link, default: enforce non-negativity
          rep_len(ntBeta)
      } else {
        ## User-defined lower bounds for parameters to be estimated
        DEoptimLowerBound[1:ntBeta]
      }

  ## Theta
    nlminbUpperBound <- c(nlminbUpperBound,
      if(is.null(p(sim)$ub$theta)){
        rep_len(Inf, ntTheta)
      } else {
        ## User-defined lower bounds for parameters to be estimated
        DEoptimUpperBound[(ntBeta + 1L):nt]
      })
  
    nlminbLowerBound <- c(nlminbLowerBound,
      if(is.null(p(sim)$lb$theta)){
        switch(lnTheta$name,
               log = -Inf,            ## log-link, default: -Inf for terms and 0 for breakpoints/knots
               identity = 1e-30) %>%  ## identity link, default: enforce non-negativity
          rep_len(ntTheta)
      } else {
        ## User-defined lower bounds for parameters to be estimated
        DEoptimLowerBound[(ntBeta + 1L):nt]
      })

  ## Define the log-likelihood function (objective function)
  sim$nll <- parse(text = paste0("-sum(dtappareto(envData[[\"", y, "\"]], lambda=beta, theta=theta, a=", p(sim)$a, ", log=TRUE))"))

  ## If starting values are not supplied
  if (is.null(p(sim)$start)) {
    ## First optimizer, get rough estimates of the parameter values
    ## Use these estimates to compute the order of magnitude of these parameters
    #     opDE <- DEoptim::DEoptim(objFun,lower=DEoptimLowerBound,upper=DEoptimUpperBound, control = DEoptim::DEoptim.control(itermax = 500L, trace = trace))

    JDE <- list(iter = 0L)
    i <- 0L
    while(JDE$iter == 0L && i < 30){
      i <- i + 1L
      JDE.call <- quote(JDEoptim(fn = objFun, lower = DEoptimLowerBound, upper = DEoptimUpperBound, trace = if(trace > 0) TRUE else FALSE, triter = trace))
      JDE.call[names(formals(objFun)[-1])] <- parse(text = formalArgs(objFun)[-1])
      JDE <- suppressWarnings(eval(JDE.call))
    }

    ## Update scaling matrix
    diag(scalMx) <- oom(JDE$par)


    ## Second optimization with nlminb()
    ## Brute-force to make models converge & select the best fit (according to AICc criterion)
    svList <- c(lapply(1:500,function(i)pmin(pmax(rnorm(length(JDE$par),0L,2L)/10 + unname(JDE$par/oom(JDE$par)), nlminbLowerBound), nlminbUpperBound)),
                list(unname(JDE$par/oom(JDE$par))))

    nlminb.call <- quote(nlminb(start=sv, objective = objFun, lower = nlminbLowerBound, upper = nlminbUpperBound,
                                control=c(p(sim)$nlminb.control, list(trace = trace))))
    nlminb.call[names(formals(objFun)[-1L])] <- parse(text = formalArgs(objFun)[-1L])

    out <- lapply(svList,function(sv){
      op <- eval(nlminb.call)

      i <- 1L

      ## When there is no convergence and restart is possible, run nlminb() again
      while(as.integer(gsub("[\\(\\)]", "", regmatches(op$message, gregexpr("\\(.*?\\)", op$message))[[1L]])) %in% 7:14 & i < 3L){
        i <- i + 1L
        op <- eval(nlminb.call)
      }
      op
    })

    ## Select best minimum amongst all trials
    op <- out[[which.min(sapply(out,"[[","objective"))]]

    ## If starting values are supplied
  } else if (is.list(p(sim)$start)) {

    nlminb.call <- quote(nlminb(start=sv, objective = objFun, lower = nlminbLowerBound, upper = nlminbUpperBound,
                                control=c(p(sim)$nlminb.control, list(trace = trace))))
    nlminb.call[names(formals(objFun)[-1L])] <- parse(text = formalArgs(objFun)[-1L])

    ## List of vectors of user-defined starting values
    out <- lapply(p(sim)$start, function(sv) {
      op <- eval(nlminb.call)

      i <- 1L

      ## When there is no convergence and restart is possible, run nlminb() again
      while(as.integer(gsub("[\\(\\)]", "", regmatches(op$message, gregexpr("\\(.*?\\)", op$message))[[1L]])) %in% 7:14 & i < 3L) {
        i <- i + 1L
        op <- eval(nlminb.call)
      }
      op
    })

    ## Select best minimum amongst all trials
    op <- out[[which.min(sapply(out,"[[","objective"))]]
  } else if (is.vector(p(sim)$start)) {
    ## Vector of user-defined starting values
    op <- nlminb(p(sim)$start,
                 objective = objFun,
                 lower = nlminbLowerBound,
                 upper = nlminbUpperBound,
                 control = p(sim)$nlminb.control)
  }

  ## Compute the standard errors around the estimates
  hess.call <- quote(numDeriv::hessian(func = objFun, x = op$par))
  hess.call[names(formals(objFun)[-1L])] <- parse(text = formalArgs(objFun)[-1L])
  hess <- eval(hess.call)
  se <- try(drop(sqrt(diag(solve(hess))) %*% scalMx), silent = TRUE)

  ## Negative values in the Hessian matrix suggest that the algorithm did not converge
  if(anyNA(se)) warning("fireSense_SizeFit> nlminb: algorithm did not converge", noBreaks. = TRUE)

  ## Parameters scaling: Revert back estimated coefficients to their original scale
  op$par <- drop(op$par %*% scalMx)

  sim$fireSense_SizeFitted <- list(formula = p(sim)$formula,
                                   linkFunBeta = lnBeta,
                                   linkFunTheta = lnTheta,
                                   coefBeta = setNames(op$par[1:ntBeta], colnames(mmBeta)),
                                   coefTheta = setNames(op$par[(ntBeta + 1L):nt], colnames(mmTheta)),
                                   seBeta = setNames(se[1:ntBeta], colnames(mmBeta)),
                                   seTheta = setNames(se[(ntBeta + 1L):nt], colnames(mmTheta)))

  class(sim$fireSense_SizeFitted) <- "fireSense_SizeFit"
  
  if (!is.na(p(sim)$intervalRunModule))
    sim <- scheduleEvent(sim, time(sim) + p(sim)$intervalRunModule, "fireSense_SizePredict", "run")
  
  sim
}
