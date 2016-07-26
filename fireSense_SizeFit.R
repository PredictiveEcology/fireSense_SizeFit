# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense_SizeFit",
  description = "Fit statistical models describing the empirical distribution of fire sizes.",
  keywords = c("fire size distribution", "tapered Pareto", "optimization", "fireSense", "statistical model"),
  authors=c(person("Jean", "Marchal", email="jean.d.marchal@gmail.com", role=c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("0.0.0.9000"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_SizeFit.Rmd"),
  reqdPkgs = list("data.table", "DEoptimR", "magrittr", "numDeriv", "PtProcess"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", default, min, max, "parameter description")),
    defineParameter("formula", "list", NULL, desc = 'a (named) list with two components called "beta" and "theta" of class "formula" : a symbolic description of the model to be fitted (beta and theta parameters of the tapered Pareto).'),
    defineParameter("a", "numeric", NULL, desc = 'range parameter a of the tapered Pareto. The random variable x take values on the interval a <= x < Inf. Values outside of this range are removed with a warning.'),
    defineParameter("link", "list", list(beta = "log", theta = "identity"), desc = 'a (named) list of two link functions for the beta parameter of the tapered Pareto. see ?family. It should contain two components named "beta" and "theta".'),
    defineParameter(name = "start", class = "numeric", default = NULL, 
                    desc = 'optional (named) list of starting value(s) for the parameter(s) to be estimated. Should contain two components named "beta" and "theta".'),
    defineParameter(name = "lb", class = "numeric", default = NULL, desc = 'optional lower bound(s) for the parameter(s) to be estimated. Should contain two components named "beta" and "theta".'),
    defineParameter(name = "ub", class = "numeric", default = NULL, desc = 'optional upper bound(s) for the parameter(s) to be estimated. Should contain two components named "beta" and "theta".'),
    defineParameter(name = "nlminb.control", class = "numeric", default = list(iter.max = 5000L, eval.max=5000L),
                    desc = "optional list of control parameters to be passed to the nlminb optmizer. See ?nlminb"),
    defineParameter(name = "trace", class = "numeric", default = 0,
                    desc = "non-negative integer. If > 0, tracing information on the progress of the optimization is produced every trace iteration. Defaults to 0 which indicates no trace information is to be printed.")    # defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    # defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    # defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    # defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur")
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
  objFun <- function(params, scalMx, mmBeta, mmTheta, ntBeta, nt, lnBeta, lnTheta, y, sim){
    
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
    #sim <- scheduleEvent(sim, params(sim)$fireSense_SizeFit$.plotInitialTime, "fireSense_SizeFit", "plot")
    #sim <- scheduleEvent(sim, params(sim)$fireSense_SizeFit$.saveInitialTime, "fireSense_SizeFit", "save")
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

  ## Init/Fit function
    fireSense_SizeFitInit <- function(sim) {
      ## Check formula for beta
      if (is.empty.model(params(sim)$fireSense_SizeFit$formula$beta))
        stop("fireSense_SizeFit> The specified formula (beta) is empty")
      
      ## Check formula for theta
      if (is.empty.model(params(sim)$fireSense_SizeFit$formula$theta))
        stop("fireSense_SizeFit> The specified formula (theta) is empty")
      
      ## Coerce lnBeta to a link-glm object
      if (is.character(params(sim)$fireSense_SizeFit$link$beta)) {
        lnBeta <- make.link(params(sim)$fireSense_SizeFit$link$beta) 
      } else if (is(params(sim)$fireSense_SizeFit$link$beta, "link-glm")) {
        ## Do nothing
      } else {
        ## Try to coerce to link-glm class
        lnBeta <- make.link(params(sim)$fireSense_SizeFit$link$beta)
      }
      
      ## Coerce lnTheta to a link-glm object
      if (is.character(params(sim)$fireSense_SizeFit$link$theta)) {
        lnTheta <- make.link(params(sim)$fireSense_SizeFit$link$theta)
      } else if (is(params(sim)$fireSense_SizeFit$link$theta, "link-glm")) {
        ## Do nothing
      } else {
        lnTheta <- make.link(params(sim)$fireSense_SizeFit$link$theta)
      }
      
      ## No tracing if trace < 0
      trace <- if(params(sim)$fireSense_SizeFit$trace < 0) 0 else params(sim)$fireSense_SizeFit$trace
      
      
      tfBeta <- terms.formula(params(sim)$fireSense_SizeFit$formula$beta)
      
      if (attr(tfBeta, "response")) {
        y <- yBeta <- as.character(attr(tfBeta, "variables")[[2L]])
      } else {
        stop("fireSense_SizeFit> Incomplete formula (beta), the LHS is missing")
      }
      
      tfTheta <- terms.formula(params(sim)$fireSense_SizeFit$formula$theta)
      
      if (attr(tfTheta, "response")) {
        yTheta <- as.character(attr(tfTheta, "variables")[[2L]])
      } else {
        stop("fireSense_SizeFit> Incomplete formula (theta), the LHS is missing")
      }
      
      if (!identical(yBeta, yTheta)) {
        stop("fireSense_SizeFit> the response variables for beta and theta must be identical")
      } else {
        yExpr <- parse(text = y)
      }
      
      ## Coerce data to data.table if not already one
      sim$fireSense_SizeFit$data <- if (!is(sim$dataFireSense_SizeFit, "data.table")) {
        data.table(sim$dataFireSense_SizeFit)
      } else {
        sim$dataFireSense_SizeFit
      }
      
      ## If there are rows in the dataset where y < a, remove them
      rm <- sim$fireSense_SizeFit$data[, y] < params(sim)$fireSense_SizeFit$a
      
      if (all(rm)){
        stop("fireSense_SizeFit> All rows contain values outside of the range a <= x < Inf.")
      } else if (any(rm)) {
        sim$fireSense_SizeFit$data <- sim$fireSense_SizeFit$data[!rm, ]
        warning("fireSense_SizeFit> Removed", sum(rm), "rows containing values outside of the range a <= x < Inf")  
      }
      
      ## Number of terms
      ntBeta <- length(labels(tfBeta)) + attr(tfBeta, "intercept")
      ntTheta <- length(labels(tfTheta)) + attr(tfTheta, "intercept")
      nt <- ntBeta + ntTheta
      
      ## Define the scaling matrices. This is used later in the optimization process
      ## to rescale parameter values between 0 and 1, i.e. put all variables on the same scale.
      scalMx <- matrix(0L, nt, nt)
      diag(scalMx) <- 1L
      
      ## Design matrices
      mmBeta <- model.matrix(tfBeta, sim$fireSense_SizeFit$data)
      mmTheta <- model.matrix(tfTheta, sim$fireSense_SizeFit$data)  
      
      ## Define parameter bounds automatically if they are not supplied by user
      ## First defined the bounds for DEoptim, the first optimizer
      ## Beta
      switch(lnBeta$name, 
             log = {
               DEoptimUpperBound <- 
                 if (is.null(params(sim)$fireSense_SizeFit$ub$theta)) {
                   ## Automatically estimate an upper boundary for each parameter       
                   (lm(update(params(sim)$fireSense_SizeFit$formula$theta, log(.) ~ .),
                       y = FALSE,
                       model = FALSE,
                       data = sim$fireSense_SizeFit$data) %>%
                      coef %>%
                      abs) * 1.1
                 } else {
                   ## User-defined bounds
                   rep_len(params(sim)$fireSense_SizeFit$ub$theta, ntBeta) ## Recycled if necessary
                 }
               
               DEoptimLowerBound <- 
                 if (is.null(params(sim)$fireSense_SizeFit$lb$theta)) {
                   ## Automatically estimate a lower boundary for each parameter 
                   -DEoptimUpperBound
                 } else {
                   ## User-defined bounds
                   rep_len(params(sim)$fireSense_SizeFit$lb$theta, ntBeta) ## Recycled if necessary
                 }       
             }, identity = {
               DEoptimUpperBound <-
                 if (is.null(params(sim)$fireSense_SizeFit$ub$theta)) {
                   ## Automatically estimate an upper boundary for each parameter
                   (lm(params(sim)$fireSense_SizeFit$formula$theta,
                       y = FALSE,
                       model = FALSE,
                       data = sim$fireSense_SizeFit$data) %>%
                      coef %>%
                      abs) * 1.1
                 } else {
                   ## User-defined bounds
                   rep_len(params(sim)$fireSense_SizeFit$ub$theta, ntBeta) ## Recycled if necessary
                 } 
               
               DEoptimLowerBound <-
                 if(is.null(params(sim)$fireSense_SizeFit$lb$theta)){
                   ## Automatically estimate an lower boundary for each parameter 
                   rep_len(1e-30, ntBeta) ## Enforce non-negativity, recycled if necessary
                 } else {
                   ## User-defined bounds
                   rep_len(params(sim)$fireSense_SizeFit$lb$theta, ntBeta) ## Recycled if necessary
                 }     
             }, stop(paste("fireSense_SizeFit> Link function", params(sim)$fireSense_SizeFit$link$theta, "(theta) is not supported.")))
      
      ## Theta
      switch(lnTheta$name, 
             log = {
               DEoptimUpperBound <- c(DEoptimUpperBound,
                                      if (is.null(params(sim)$fireSense_SizeFit$ub$theta)) {
                                        ## Automatically estimate an upper boundary for each parameter       
                                        (lm(update(params(sim)$fireSense_SizeFit$formula$theta, log(.) ~ .),
                                            y = FALSE,
                                            model = FALSE,
                                            data = sim$fireSense_SizeFit$data) %>%
                                           coef %>%
                                           abs) * 1.1
                                      } else {
                                        ## User-defined bounds
                                        rep_len(params(sim)$fireSense_SizeFit$ub$theta, ntTheta) ## Recycled if necessary
                                      })
               
               DEoptimLowerBound <- c(DEoptimLowerBound,
                                      if(is.null(params(sim)$fireSense_SizeFit$lb$theta)){
                                        ## Automatically estimate a lower boundary for each parameter 
                                        -DEoptimUpperBound[(ntBeta + 1L):nt]
                                      } else {
                                        ## User-defined bounds
                                        rep_len(params(sim)$fireSense_SizeFit$lb$theta, ntTheta) ## Recycled if necessary
                                      })
             }, identity = {
               DEoptimUpperBound <- c(DEoptimUpperBound,
                                      if (is.null(params(sim)$fireSense_SizeFit$ub$theta)) {
                                        ## Automatically estimate an upper boundary for each parameter
                                        (lm(params(sim)$fireSense_SizeFit$formula$theta,
                                            y = FALSE,
                                            model = FALSE,
                                            data = sim$fireSense_SizeFit$data) %>%
                                           coef %>%
                                           abs) * 1.1
                                      } else {
                                        ## User-defined bounds
                                        rep_len(params(sim)$fireSense_SizeFit$ub$theta, ntTheta) ## Recycled if necessary
                                      })
               
               DEoptimLowerBound <- c(DEoptimLowerBound,
                                      if (is.null(params(sim)$fireSense_SizeFit$lb$theta)) {
                                        ## Automatically estimate an lower boundary for each parameter 
                                        rep_len(1e-30, ntTheta) ## Enforce non-negativity, recycled if necessary
                                      } else {
                                        ## User-defined bounds
                                        rep_len(params(sim)$fireSense_SizeFit$lb$theta, ntTheta) ## Recycled if necessary
                                      })
             }, stop(paste("fireSense_SizeFit> Link function", params(sim)$fireSense_SizeFit$link$theta, "(theta) is not supported.")))
      
      ## Then, define lower and upper bounds for the second optimizer (nlminb)
      ## Beta
      nlminbUpperBound <-
        if(is.null(params(sim)$fireSense_SizeFit$ub$beta)){
          rep_len(Inf, ntBeta)
        } else {
          ## User-defined lower bounds for parameters to be estimated
          DEoptimUpperBound[1:ntBeta]
        }
      
      nlminbLowerBound <- 
        if(is.null(params(sim)$fireSense_SizeFit$lb$beta)){
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
                            if(is.null(params(sim)$fireSense_SizeFit$ub$theta)){
                              rep_len(Inf, ntTheta)
                            } else {
                              ## User-defined lower bounds for parameters to be estimated
                              DEoptimUpperBound[(ntBeta + 1L):nt]
                            })
      
      nlminbLowerBound <- c(nlminbLowerBound,
                            if(is.null(params(sim)$fireSense_SizeFit$lb$theta)){
                              switch(lnTheta$name,
                                     log = -Inf,            ## log-link, default: -Inf for terms and 0 for breakpoints/knots
                                     identity = 1e-30) %>%  ## identity link, default: enforce non-negativity
                                rep_len(ntTheta)
                            } else {
                              ## User-defined lower bounds for parameters to be estimated
                              DEoptimLowerBound[(ntBeta + 1L):nt]
                            })
      
      ## Define the log-likelihood function (objective function)
      sim$nll <- parse(text = paste0("-sum(dtappareto(sim$fireSense_SizeFit$data[, ", y, "], lambda=beta, theta=theta, a=", params(sim)$fireSense_SizeFit$a, ", log=TRUE))"))
      
      ## If starting values are not supplied
      if (is.null(params(sim)$fireSense_SizeFit$start)) {
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
                                    control=c(params(sim)$fireSense_SizeFit$nlminb.control, list(trace = trace))))
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
      } else if (is.list(params(sim)$fireSense_SizeFit$start)) {
        
        nlminb.call <- quote(nlminb(start=sv, objective = objFun, lower = nlminbLowerBound, upper = nlminbUpperBound, 
                                    control=c(params(sim)$fireSense_SizeFit$nlminb.control, list(trace = trace))))
        nlminb.call[names(formals(objFun)[-1L])] <- parse(text = formalArgs(objFun)[-1L])
        
        ## List of vectors of user-defined starting values
        out <- lapply(params(sim)$fireSense_SizeFit$start, function(sv) {
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
      } else if (is.vector(params(sim)$fireSense_SizeFit$start)) {
        ## Vector of user-defined starting values
        op <- nlminb(params(sim)$fireSense_SizeFit$start, 
                     objective = objFun,
                     lower = nlminbLowerBound,
                     upper = nlminbUpperBound,
                     control = params(sim)$fireSense_SizeFit$nlminb.control)
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
      
      fit <- list(formula = params(sim)$fireSense_SizeFit$formula,
                  linkFunBeta = lnBeta,
                  linkFunTheta = lnTheta,
                  coefBeta = setNames(op$par[1:ntBeta], colnames(mmBeta)),
                  coefTheta = setNames(op$par[(ntBeta + 1L):nt], colnames(mmTheta)),
                  seBeta = setNames(se[1:ntBeta], colnames(mmBeta)),
                  seTheta = setNames(se[(ntBeta + 1L):nt], colnames(mmTheta)))
      fit$fitted.values$beta <- fit$linkFunBeta$linkinv(drop(mmBeta %*% fit$coefBeta))
      fit$fitted.values$theta <- fit$linkFunTheta$linkinv(drop(mmTheta %*% fit$coefTheta))
      sim$fireSense_SizeFit <- fit
      class(sim$fireSense_SizeFit) <- "fireSense_SizeFit"
      invisible(sim)
    }
