# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense_SizeFit",
  description = "Fit statistical models describing the empirical distribution of fire sizes. A tapered Pareto distribution is assumed.",
  keywords = c("fire size distribution", "tapered Pareto", "optimization", "fireSense", "statistical model"),
  authors=c(person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut", "cre"))),
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
    defineParameter("formula", "list", NULL, 
      desc = 'a (named) list with two components called "beta" and "theta" of class "formula" :
              a symbolic description of the model to be fitted (beta and theta being parameters
              of the tapered Pareto distribution).'),
    defineParameter("a", "numeric", NULL, 
      desc = 'range parameter a of the tapered Pareto. The random variable x take values on the
              interval a <= x < Inf. Values outside of this range are ignored with a warning.'),
    defineParameter("link", "list", list(beta = "log", theta = "identity"), 
      desc = 'a (named) list with two components called "beta" and "theta" specifying model 
              links for the "beta" and "theta" parameters of the tapered Pareto. These can be 
              character strings or objects of class "link-glm". For more info see ?family.'),
    defineParameter(name = "start", class = "list", default = NULL,
      desc = 'optional (named) list with two components called "beta" and "theta" specifying 
              starting values for the parameters to be estimated. Those are passed to nlminb
              and can be numeric vectors, or lists of numeric vectors.'),
    defineParameter(name = "lb", class = "numeric", default = NULL,
      desc = 'optional (named) list with two components called "beta" and "theta" specifying
              lower bounds for the parameters to be estimated. These should be numeric values.'),
    defineParameter(name = "ub", class = "numeric", default = NULL, 
      desc = 'optional (named) list with two components called "beta" and "theta" specifying
              upper bounds for the parameters to be estimated. These should be numeric values.'),
    defineParameter(name = "nlminb.control", class = "numeric", default = list(iter.max = 5e3L, eval.max = 5e3L),
      desc = "optional list of control parameters to be passed to the nlminb optmizer. See ?nlminb"),
    defineParameter(name = "trace", class = "numeric", default = 0,
      desc = "non-negative integer. If > 0, tracing information on the progress of the 
              optimization is produced every trace iteration. Defaults to 0 which indicates no
              trace information should be printed."),
    defineParameter(name = "data", class = "character", default = NULL,
      desc = "optional. A character vector indicating the names of objects present in the 
              simList environment, in which to look for variables with which to predict. Objects
              should be data.frames. If omitted, or if variables are not found in data objects,
              variables are searched in the simList environment."),
    defineParameter(name = "initialRunTime", class = "numeric", default = NA,
      desc = "optional. Simulation time at which to start this module. If omitted, start at start(simList)."),
    defineParameter(name = "intervalRunModule", class = "numeric", default = NA, 
      desc = "optional. Interval in simulation time units between two runs of this module.")
  ),
  inputObjects = data.frame(objectName = "dataFireSense_SizeFit",
                            objectClass = "data.frame",
                            sourceURL = "",
                            other = NA_character_,
                            stringsAsFactors = FALSE),
  outputObjects = data.frame(
    objectName = "fireSense_SizeFit",
    objectClass = "fireSense_SizeFit",
    other = NA_character_,
    stringsAsFactors = FALSE
  )
))


## event types
#   - type `init` is required for initialiazation

doEvent.fireSense_SizeFit = function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {

    sim <- sim$fireSense_SizeFitInit(sim)

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

  ## Toolbox: set of functions used internally by the module
    ## Compute the order of magnitude
      oom <- function(x) 10 ^ (ceiling(log10(abs(x))))
    
    ## Function to pass to the optimizer
      objfun <- function(params, scalMx, mmBeta, mmTheta, nBeta, n, lnBeta, lnTheta, y, envData){
        
        ## Parameters scaling
        params <- drop(params %*% scalMx)
        
        beta <- drop(mmBeta %*% params[1:nBeta])
        theta <- drop(mmTheta %*% params[(nBeta + 1L):n])
        
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
  
    ## Nlminb wrapper
      objNlminb <- function(start, objective, lower, upper, control) {
        
        nlminb.call <- quote(nlminb(start = start, objective = objective, lower = lower, upper = upper, control = control))
        nlminb.call[names(formals(objective)[-1L])] <- parse(text = formalArgs(objective)[-1L])
        
        o <- eval(nlminb.call)
        
        i <- 1L
        
        ## When there is no convergence and restart is possible, run nlminb() again
        while(as.integer(gsub("[\\(\\)]", "", regmatches(o$message, gregexpr("\\(.*?\\)", o$message))[[1L]])) %in% 7:14 & i < 3L){
          i <- i + 1L
          o <- eval(nlminb.call)
        }            
        o
      }  
      
      
  envData <- new.env(parent = envir(sim))
  on.exit(rm(envData))
  list2env(as.list(envir(sim)), envir = envData)

  ## Handle data arg
  if (!is.null(p(sim)$data))
    lapply(p(sim)$data, function(x, envData) if (is.list(sim[[x]])) list2env(sim[[x]], envir = envData), envData = envData)

  ## Check formula for beta
  if (is.empty.model(p(sim)$formula$b))
    stop("fireSense_SizeFit> The formula (beta) describes an empty model.")
  
  ## Check formula for theta
  if (is.empty.model(p(sim)$formula$t))
    stop("fireSense_SizeFit> The formula (theta) describes an empty model.")

  termsBeta <- terms.formula(formulaBeta <- p(sim)$formula$b)
  termsTheta <- terms.formula(formulaTheta <- p(sim)$formula$t)
    
  if (attr(termsBeta, "response")) y <- yBeta <- as.character(formulaBeta[[2L]])
  else stop("fireSense_SizeFit> Incomplete formula (beta), the LHS is missing.")
  
  if (attr(termsTheta, "response")) yTheta <- as.character(formulaTheta[[2L]])
  else stop("fireSense_SizeFit> Incomplete formula (theta), the LHS is missing.")
  
  if (!identical(yBeta, yTheta))
    stop("fireSense_SizeFit> the response variables for beta and theta must be identical.")
    
  ## Coerce lnBeta to a link-glm object
  if (is.character(p(sim)$link$b)) {
    lnBeta <- make.link(p(sim)$link$b)
  } else if (is(p(sim)$link$b, "link-glm")) {
    ## Do nothing
  } else lnBeta <- make.link(p(sim)$link$b) ## Try to coerce to link-glm class

  ## Coerce lnTheta to a link-glm object
  if (is.character(p(sim)$link$t)) {
    lnTheta <- make.link(p(sim)$link$t)
  } else if (is(p(sim)$link$t, "link-glm")) {
    ## Do nothing
  } else lnTheta <- make.link(p(sim)$link$t)

  ## No tracing if trace < 0
  trace <- if(p(sim)$trace < 0) 0 else p(sim)$trace

  ## If there are rows in the dataset where y < a, remove them
  rm <- envData[[y]] < p(sim)$a

  if (all(rm)){
    
    stop("fireSense_SizeFit> All rows contain values outside of the range a <= x < Inf.")
    
  } else if (any(rm)) {
    
    lapply(unique(c(all.vars(termsBeta), all.vars(termsTheta))), function(x) assign(x = x, value = envData[[x]][!rm], envir = envData))
    warning(paste("fireSense_SizeFit> Ignored", sum(rm), "rows containing values outside of the range a <= x < Inf."), immediate. = TRUE)
    
  }

  ## Number of terms
  nBeta <- length(labels(termsBeta)) + attr(termsBeta, "intercept")
  nTheta <- length(labels(termsTheta)) + attr(termsTheta, "intercept")
  n <- nBeta + nTheta

  ## Define the scaling matrices. This is used later in the optimization process
  ## to rescale parameter values between 0 and 1, i.e. put all variables on the same scale.
  scalMx <- matrix(0, n, n)
  diag(scalMx) <- 1

  ## Design matrices
  mmBeta <- model.matrix(termsBeta, envData)
  mmTheta <- model.matrix(termsTheta, envData)

  ## Define parameter bounds automatically if they are not supplied by user
  ## First defined the bounds for DEoptim, the first optimizer
    ## Beta
      switch(lnBeta$name,
             log = {
               DEoptimUB <-
                 if (is.null(p(sim)$ub$b)) {
                   ## Automatically estimate an upper boundary for each parameter
                   (lm(update(formulaBeta, log(.) ~ .),
                       y = FALSE,
                       model = FALSE,
                       data = envData) %>%
                      coef %>%
                      abs) * 1.1
                 } else rep_len(p(sim)$ub$b, nBeta) ## User-defined bounds (recycled if necessary)

               DEoptimLB <-
                 if (is.null(p(sim)$lb$b))
                   -DEoptimUB ## Automatically estimate a lower boundary for each parameter
                 else
                   rep_len(p(sim)$lb$b, nBeta) ## User-defined bounds (recycled if necessary)
             }, identity = {
               DEoptimUB <-
                 if (is.null(p(sim)$ub$b)) {
                   ## Automatically estimate an upper boundary for each parameter
                   (lm(formulaBeta,
                       y = FALSE,
                       model = FALSE,
                       data = envData) %>%
                      coef %>%
                      abs) * 1.1
                 } else rep_len(p(sim)$ub$b, nBeta) ## User-defined bounds (recycled if necessary)

               DEoptimLB <-
                 if (is.null(p(sim)$lb$b))
                   rep_len(1e-30, nBeta) ## Enforce non-negativity (reecycled if necessary)
                 else
                   rep_len(p(sim)$lb$b, nBeta) ## User-defined bounds (recycled if necessary)
             }, stop(paste("fireSense_SizeFit> Link function", p(sim)$link$b, "(beta) is not supported.")))

  ## Theta
    switch(lnTheta$name,
           log = {
             DEoptimUB <- 
               c(DEoptimUB,
                 if (is.null(p(sim)$ub$t)) {
                   ## Automatically estimate an upper boundary for each parameter
                   (lm(update(formulaTheta, log(.) ~ .),
                       y = FALSE,
                       model = FALSE,
                       data = envData) %>%
                      coef %>%
                      abs) * 1.1
                 } else rep_len(p(sim)$ub$t, nTheta)
               ) ## User-defined bounds (recycled if necessary)
  
             DEoptimLB <-
               c(DEoptimLB,
                 if (is.null(p(sim)$lb$t)) -DEoptimUB[(nBeta + 1L):n] ## Automatically estimate a lower boundary for each parameter
                 else rep_len(p(sim)$lb$t, nTheta) ## User-defined bounds (recycled if necessary)
               )
           }, identity = {
             DEoptimUB <-
               c(DEoptimUB,
                 if (is.null(p(sim)$ub$t)) {
                   ## Automatically estimate an upper boundary for each parameter
                   (lm(formulaTheta,
                       y = FALSE,
                       model = FALSE,
                       data = envData) %>%
                      coef %>%
                      abs) * 1.1
                 } else rep_len(p(sim)$ub$t, nTheta) ## User-defined bounds (recycled if necessary)
               )
  
             DEoptimLB <- 
               c(DEoptimLB,
                 if (is.null(p(sim)$lb$t)) rep_len(1e-30, nTheta) ## Enforce non-negativity (recycled if necessary)
                 else rep_len(p(sim)$lb$t, nTheta) ## User-defined bounds (recycled if necessary)
               )
           }, stop(paste("fireSense_SizeFit> Link function", p(sim)$link$t, "(theta) is not supported.")))

  ## Then, define lower and upper bounds for the second optimizer (nlminb)
    ## Beta
      nlminbUB <-
        if(is.null(p(sim)$ub$b)) rep_len(Inf, nBeta)
        else DEoptimUB[1:nBeta] ## User-defined bounds
    
      nlminbLB <-
        if (is.null(p(sim)$lb$b)) {
          
          switch(lnBeta$name,
                 log = -Inf,            ## log-link, default: -Inf for terms and 0 for breakpoints/knots
                 identity = 1e-30) %>%  ## identity link, default: enforce non-negativity
            rep_len(nBeta)
          
        } else DEoptimLB[1:nBeta] ## User-defined bounds
  
    ## Theta
      nlminbUB <- c(nlminbUB,
        if (is.null(p(sim)$ub$t)) rep_len(Inf, nTheta)
        else DEoptimUB[(nBeta + 1L):n] ## User-defined bounds
      )
    
      nlminbLB <- c(nlminbLB,
        if (is.null(p(sim)$lb$t)) {
          
          switch(lnTheta$name,
                 log = -Inf,            ## log-link, default: -Inf for terms and 0 for breakpoints/knots
                 identity = 1e-30) %>%  ## identity link, default: enforce non-negativity
            rep_len(nTheta)
        
        } else DEoptimLB[(nBeta + 1L):n] ## User-defined bounds
      )

  ## Define the log-likelihood function (objective function)
  sim$nll <- parse(text = paste0("-sum(dtappareto(envData[[\"", y, "\"]], lambda=beta, theta=theta, a=", p(sim)$a, ", log=TRUE))"))

  ## If starting values are not supplied
  if (is.null(p(sim)$start)) {
    ## First optimizer, get rough estimates of the parameter values
    ## Use these estimates to compute the order of magnitude of these parameters
    #     opDE <- DEoptim::DEoptim(objfun,lower=DEoptimLB,upper=DEoptimUB, control = DEoptim::DEoptim.control(itermax = 500L, trace = trace))

      JDE <- list(iter = 0L)
      i <- 0L
      while(JDE$iter == 0L && i < 30){
        i <- i + 1L
        JDE.call <- quote(JDEoptim(fn = objfun, lower = DEoptimLB, upper = DEoptimUB, trace = if(trace > 0) TRUE else FALSE, triter = trace, maxiter = 10))
        JDE.call[names(formals(objfun)[-1])] <- parse(text = formalArgs(objfun)[-1])
        JDE <- suppressWarnings(eval(JDE.call))
      }
  
      ## Update scaling matrix
      diag(scalMx) <- oom(JDE$par)


    ## Second optimization with nlminb()
    ## Brute-force to make models converge & select the best fit (according to AICc criterion)
      svList <- c(lapply(1:500,function(i)pmin(pmax(rnorm(length(JDE$par),0L,2L)/10 + unname(JDE$par/oom(JDE$par)), nlminbLB), nlminbUB)),
                  list(unname(JDE$par/oom(JDE$par))))
      
      out <- lapply(svList, objNlminb, objective = objfun, lower = nlminbLB, upper = nlminbUB, control = c(p(sim)$nlminb.control, list(trace = trace)))
      
      ## Select best minimum amongst all trials
      out <- out[[which.min(sapply(out, "[[", "objective"))]]

  } else if (is.list(p(sim)$start)) { ## If starting values are supplied as a list of vectors of starting values

    out <- lapply(p(sim)$start, objNlminb, objective = objfun, lower = nlminbLB, upper = nlminbUB, control = c(p(sim)$nlminb.control, list(trace = trace)))

    ## Select best minimum amongst all trials
    out <- out[[which.min(sapply(out, "[[", "objective"))]]
    
  } else if (is.vector(p(sim)$start)) { ## If starting values are supplied as a vector of starting values
      
    out <- objNlminb(p(sim)$start, objfun, nlminbLB, nlminbUB, c(p(sim)$nlminb.control, list(trace = trace)))

  }

  ## Compute the standard errors around the estimates
  hess.call <- quote(numDeriv::hessian(func = objfun, x = o$par))
  hess.call[names(formals(objfun)[-1L])] <- parse(text = formalArgs(objfun)[-1L])
  hess <- eval(hess.call)
  se <- try(drop(sqrt(diag(solve(hess))) %*% scalMx), silent = TRUE)

  ## Negative values in the Hessian matrix suggest that the algorithm did not converge
  if(anyNA(se)) warning("fireSense_SizeFit> nlminb: algorithm did not converge", noBreaks. = TRUE)

  ## Parameters scaling: Revert back estimated coefficients to their original scale
  o$par <- drop(o$par %*% scalMx)

  sim$fireSense_SizeFitted <- 
    list(formula = p(sim)$formula,
         link = list(beta = lnBeta, theta = lnTheta),
         coef = list(beta = setNames(o$par[1:nBeta], colnames(mmBeta)),
                     theta = setNames(o$par[(nBeta + 1L):n], colnames(mmTheta))),
         se = list(beta = setNames(se[1:nBeta], colnames(mmBeta)),
                   theta = setNames(se[(nBeta + 1L):n], colnames(mmTheta))))
  
  class(sim$fireSense_SizeFitted) <- "fireSense_SizeFit"
  
  if (!is.na(p(sim)$intervalRunModule))
    sim <- scheduleEvent(sim, time(sim) + p(sim)$intervalRunModule, "fireSense_SizePredict", "run")
  
  sim
  
}
