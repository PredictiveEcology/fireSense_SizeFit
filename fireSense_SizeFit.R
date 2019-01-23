# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense_SizeFit",
  description = "Fit statistical models of the fire size distribution. A tapered
                 Pareto distribution is assumed. This distribution has three
                 parameters: a, beta and theta. a is the lower truncation point 
                 and is assumed to be known a priori, beta controls the rate of
                 frequency decrease as the fire size increases, and theta
                 governs the location of the exponential taper. This module can 
                 be used to assess the relation between two parameters of the 
                 tapered Pareto distribution, beta and theta, and environmental 
                 controls of the fire size distribution.",
  keywords = c("fire size distribution", "tapered Pareto", "optimization", "fireSense", "statistical model"),
  authors=person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut", "cre")),
  childModules = character(),
  version = list(SpaDES.core = "0.1.0", fireSense_SizeFit = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_SizeFit.Rmd"),
  reqdPkgs = list("DEoptim", "magrittr", "numDeriv", "parallel", "PtProcess"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", default, min, max, "parameter description")),
    defineParameter("formula", "list", list(beta = NA, theta = NA), 
                    desc = "a named list with two elements, 'beta' and 'theta',
                            describing the model to be fitted. 'beta' and 
                            'theta' should be formulas (see `?formula`)"),
    defineParameter("a", "numeric", NULL, 
                    desc = "lower truncation point a of the tapered Pareto. The
                            random variable x take values on the interval 
                            a <= x < Inf. Values outside of this range are 
                            ignored with a warning."),
    defineParameter("link", "list", list(beta = "log", theta = "log"), 
                    desc = "a named list with two elements, 'beta' and 'theta',
                            specifying link functions for beta and theta
                            parameters of the tapered Pareto. These can be
                            character strings or objects of class link-glm. For
                            more additional details see `?family`."),
    defineParameter(name = "data", class = "character", 
                    default = "dataFireSense_SizeFit",
                    desc = "a character vector indicating the names of objects 
                            in the `simList` environment in which to look for
                            variables present in the model formula. `data`
                            objects should be data.frames. If variables are not
                            found in `data` objects, they are searched in the
                            `simList` environment."),
    defineParameter(name = "start", class = "list", default = NULL,
                    desc = "optional named list with two elements, 'beta' and 
                            'theta', specifying starting values for the 
                            coefficients to be estimated. Those are passed to
                            `nlminb` and can be a single vector, or a list of
                            vectors. In the latter case, only the best solution,
                            that is, the one which minimizes the most the
                            objective function, is kept."),
    defineParameter(name = "lb", class = "numeric", default = NULL,
                    desc = "optional named list with two elements, 'beta' 
                            and 'theta', specifying lower bounds for the 
                            coefficients to be estimated. These must be
                            finite and will be recycled if necessary to match
                            `length(coefficients)`."),
    defineParameter(name = "ub", class = "numeric", default = NULL, 
                    desc = "optional named list with two elements, 'beta' and
                            'theta', specifying numeric vectors of upper bounds
                            for the coefficients to be estimated. These 
                            must be finite and will be recycled if necessary to
                            match `length(coefficients)`."),
    defineParameter(name = "itermax", class = "integer", default = 2000,
                    desc = "integer defining the maximum number of iterations 
                            allowed (DEoptim optimizer). Default is 2000."),
    defineParameter(name = "nTrials", class = "integer", default = 500, 
                    desc = "if start is not supplied, nTrials defines 
                            the number of trials, or searches, to be performed
                            by the nlminb optimizer in order to find the best
                            solution. Default is 500."),
    defineParameter(name = "nCores", class = "integer", default = 1, 
                    desc = "non-negative integer. Defines the number of logical
                            cores to be used for parallel computation. The
                            default value is 1, which disables parallel 
                            computing."),
    defineParameter(name = "trace", class = "numeric", default = 0,
                    desc = "non-negative integer. If > 0, tracing information on
                            the progress of the optimization are printed every
                            `trace` iteration. If parallel computing is enable, 
                            nlminb trace logs are written to disk. Log files are
                            prefixed with 'fireSense_SizeFit_trace'
                            followed by the subprocess pid. Default is 0, which
                            turns off tracing."),    
    defineParameter(name = "nlminb.control", class = "numeric",
                    default = list(iter.max = 5e3L, eval.max = 5e3L),
                    desc = "optional list of control parameters to be passed to
                            the `nlminb` optimizer. See `?nlminb`."),
    defineParameter(name = ".runInitialTime", class = "numeric", default = start(sim),
                    desc = "when to start this module? By default, the start 
                            time of the simulation."),
    defineParameter(name = ".runInterval", class = "numeric", default = NA, 
                    desc = "optional. Interval between two runs of this module,
                            expressed in units of simulation time."),
    defineParameter(name = ".saveInitialTime", class = "numeric", default = NA, 
                    desc = "optional. When to start saving output to a file."),
    defineParameter(name = ".saveInitialTime", class = "numeric", default = NA, 
                    desc = "optional. Interval between save events."),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = expectsInput(
    objectName = "dataFireSense_SizeFit",
    objectClass = "data.frame",
    sourceURL = NA_character_,
    desc = "One or more data.frames in which to look for variables present in the model formula."
  ),
  outputObjects = createsOutput(
    objectName = "fireSense_SizeFit",
    objectClass = "fireSense_SizeFit",
    desc = "A fitted model object of class fireSense_SizeFit."
  )
))


## event types
#   - type `init` is required for initialiazation

doEvent.fireSense_SizeFit = function(sim, eventTime, eventType, debug = FALSE) 
{
  switch(
    eventType,
    init = { sim <- sizeFitInit(sim) },
    run = { sim <- sizeFitRun(sim) },
    save = { sim <- sizeFitSave(sim) },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  
  invisible(sim)
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

sizeFitInit <- function(sim) 
{
  moduleName <- current(sim)$moduleName
  
  # Checking parameters
  stopifnot(P(sim)$trace >= 0)
  stopifnot(P(sim)$nCores >= 1)
  stopifnot(P(sim)$itermax >= 1)
  stopifnot(P(sim)$nTrials >= 1)
  if (!is(P(sim)$formula$beta, "formula")) stop(moduleName, "> The supplied object for the 'formula' parameter (beta) is not of class formula.")
  if (!is(P(sim)$formula$theta, "formula")) stop(moduleName, "> The supplied object for the 'formula' parameter (theta) is not of class formula.")
  if (is.null(P(sim)$a)) stop(moduleName, "> Parameter 'a' is missing.")
  stopifnot(P(sim)$a > 0)
  
  sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, moduleName, "run")
  
  if (!is.na(P(sim)$.saveInitialTime))
    sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, moduleName, "save", .last())
  
  invisible(sim)
}

sizeFitRun <- function(sim)
{
  moduleName <- current(sim)$moduleName
  currentTime <- time(sim, timeunit(sim))
  endTime <- end(sim, timeunit(sim))
  
  ## Toolbox: set of functions used internally by the module
    ## Compute the order of magnitude
      oom <- function(x) 10 ^ (ceiling(log10(abs(x))))
    
    ## Function to pass to the optimizer
      objfun <- function(params, sm, mmB, mmT, nB, n, linkinvB, linkinvT, y, envData) 
      {
        ## Parameters scaling
        params <- drop(params %*% sm)
        
        beta <- drop(mmB %*% params[1:nB])
        theta <- drop(mmT %*% params[(nB + 1L):n])
        
        if(length(beta) == 1L) beta <- rep_len(beta, length(theta)) ## Recycled if needed
        if(length(theta) == 1L) theta <- rep_len(theta, length(beta)) ## Recycled if needed
        
        ## link implementation
        beta <- linkinvB(beta)
        theta <- linkinvT(theta)
        
        if(any(beta <= 0L) || anyNA(beta) || any(is.infinite(beta)) || any(theta <= 0L) || anyNA(theta) || any(is.infinite(theta))) return(1e20)
        else return(eval(nll))
      }
  
    ## Nlminb wrapper
      objNlminb <- function(start, objective, lower, upper, control)
      {
        nlminb.call <- quote(nlminb(start = start, objective = objective, lower = lower, upper = upper, control = control))
        nlminb.call[names(formals(objective)[-1L])] <- parse(text = formalArgs(objective)[-1L])
        
        op <- options(warn = -1)
        o <- eval(nlminb.call)
        
        i <- 1L
        
        ## When there is no convergence and restart is possible, run nlminb() again
        while(as.integer(gsub("[\\(\\)]", "", regmatches(o$message, gregexpr("\\(.*?\\)", o$message))[[1L]])) %in% 7:14 & i < 3L)
        {
          i <- i + 1L
          o <- eval(nlminb.call)
        }
        
        options(op)
        o
      }  
      
  # Create a container to hold the data
  envData <- new.env(parent = envir(sim))
  on.exit(rm(envData))

  # Load inputs in the data container
  list2env(as.list(envir(sim)), envir = envData)
  
  for (x in P(sim)$data) 
  {
    if (!is.null(sim[[x]])) 
    {
      if (is.data.frame(sim[[x]]))
      {
        list2env(sim[[x]], envir = envData)
      }
      else stop(moduleName, "> '", x, "' is not a data.frame.")
    }
  }

  ## Check formula for beta
  if (is.empty.model(P(sim)$formula$b))
    stop(moduleName, "> The formula (beta) describes an empty model.")
  
  ## Check formula for theta
  if (is.empty.model(P(sim)$formula$t))
    stop(moduleName, "> The formula (theta) describes an empty model.")

  termsBeta <- terms.formula(formulaBeta <- P(sim)$formula$b)
  termsTheta <- terms.formula(formulaTheta <- P(sim)$formula$t)
    
  if (attr(termsBeta, "response")) y <- yBeta <- as.character(formulaBeta[[2L]])
  else stop(moduleName, "> Incomplete formula (beta), the LHS is missing.")
  
  if (attr(termsTheta, "response")) yTheta <- as.character(formulaTheta[[2L]])
  else stop(moduleName, "> Incomplete formula (theta), the LHS is missing.")
  
  if (!identical(yBeta, yTheta))
    stop(moduleName, "> The response variable for beta and theta must be identical.")
    
  allxy = unique(sort(c(all.vars(termsBeta, "variables"),
                 all.vars(termsTheta, "variables"))))
  
  missing <- !allxy %in% ls(envData, all.names = TRUE)
  
  if (s <- sum(missing))
    stop(moduleName, "> '", allxy[missing][1L], "'",
         if (s > 1) paste0(" (and ", s-1L, " other", if (s>2) "s", ")"),
         " not found in data objects nor in the simList environment.")
  

  ## Coerce lnB to a link-glm object
  lnB <- 
    if (is.character(P(sim)$link$b))
    {
      make.link(P(sim)$link$b)
    }
    else if (is(P(sim)$link$b, "link-glm"))
    {
      P(sim)$link$b
    } 
    else lnB <- make.link(P(sim)$link$b) ## Try to coerce to link-glm class
  
  ## Coerce lnT to a link-glm object
  lnT <- 
    if (is.character(P(sim)$link$t))
    {
      make.link(P(sim)$link$t)
    } 
    else if (is(P(sim)$link$t, "link-glm"))
    {
      P(sim)$link$t
    } 
    else lnT <- make.link(P(sim)$link$t)

  linkinvB <- lnB$linkinv
  linkinvT <- lnT$linkinv
  
  ## No tracing if trace < 0
  trace <- P(sim)$trace

  ## If there are rows in the dataset where y < a, remove them
  rm <- envData[[y]] < P(sim)$a

  if (all(rm)) 
  { 
    stop(moduleName, "> All x values are outside of the range a <= x < Inf.")
  } 
  else if (any(rm))
  {
    lapply(unique(c(all.vars(termsBeta), all.vars(termsTheta))), function(x) assign(x = x, value = envData[[x]][!rm], envir = envData))
    warning(moduleName, "> Ignored ", sum(rm), " rows containing values outside of the range a <= x < Inf.", immediate. = TRUE)
  }

  ## Number of terms
  nB <- length(labels(termsBeta)) + attr(termsBeta, "intercept")
  nT <- length(labels(termsTheta)) + attr(termsTheta, "intercept")
  n <- nB + nT

  ## Define the scaling matrices. This is used later in the optimization process
  ## to rescale parameter values between 0 and 1, i.e. put all variables on the same scale.
  sm <- matrix(0, n, n)
  diag(sm) <- 1

  ## Design matrices
  mmB <- model.matrix(termsBeta, envData)
  mmT <- model.matrix(termsTheta, envData)

  ## Define parameter bounds automatically if they are not supplied by user
  ## First defined the bounds for DEoptim, the first optimizer
    DEoptimUB <- c(
      if (is.null(P(sim)$ub$b)) 
      {
        ## Automatically estimate an upper boundary for each parameter
        (suppressWarnings(
          tryCatch(
            glm(
              formulaBeta, ## family gaussian link 
              family = gaussian(link = lnB$name),
              y = FALSE,
              model = FALSE,
              data = envData
            ),
            error = function(e) stop(
              moduleName, "> Automated estimation of upper bounds", 
              " (beta) failed, please set the 'beta' element of ",
              "the 'ub' parameter."
            )
          )
        ) %>% coef %>% abs) * 1.1 -> ub
        
        if (anyNA(ub))
          stop(
            moduleName, "> Automated estimation of upper bounds (beta) failed, ",
            "please set the 'beta' element of the 'ub' parameter."
          )
        else ub
      }
      else rep_len(P(sim)$ub$b, nB), ## User-defined bounds (recycled if necessary)
      
      if (is.null(P(sim)$ub$t))
      {
        ## Automatically estimate an upper boundary for each parameter
        (suppressWarnings(
          tryCatch(
            glm(
              formulaTheta,
              family = gaussian(link = lnT$name),
              y = FALSE,
              model = FALSE,
              data = envData
            ),
            error = function(e) stop(
              moduleName, "> Automated estimation of upper bounds", 
              " (theta) failed, please set the 'theta' element of ",
              "the 'ub' parameter."
            )
          )
        ) %>% coef %>% abs) * 1.1 -> ub
        
        if (anyNA(ub))
          stop(
            moduleName, "> Automated estimation of upper bounds (theta) failed, ",
            "please set the 'theta' element of the 'ub' parameter."
          )
        else ub
      } 
      else rep_len(P(sim)$ub$t, nT) ## User-defined bounds (recycled if necessary)
    )

    ## Beta
    DEoptimLB <- 
      switch(lnB$name,
             log = {
               
               if (is.null(P(sim)$lb$b)) -abs(DEoptimUB[1:nB]) * 3 ## Automatically estimate a lower boundary for each parameter
               else rep_len(P(sim)$lb$b, nB) ## User-defined bounds (recycled if necessary)
               
             }, identity = {
               
               if (is.null(P(sim)$lb$b)) -abs(DEoptimUB[1:nB]) * 3
               else rep_len(P(sim)$lb$b, nB) ## User-defined bounds (recycled if necessary)
               
             }, stop(moduleName, "> Link function ", P(sim)$link$b, 
                     " (beta) is not supported by the process for automated estimation of lower bounds.")
      )

  ## Theta
    DEoptimLB <- c(DEoptimLB,
      switch(lnT$name,
           log = {
             
             if (is.null(P(sim)$lb$t)) -abs(DEoptimUB[(nB + 1L):n]) * 3 ## Automatically estimate a lower boundary for each parameter
             else rep_len(P(sim)$lb$t, nT) ## User-defined bounds (recycled if necessary)
             
           }, identity = {
             
             if (is.null(P(sim)$lb$t)) -abs(DEoptimUB[(nB + 1L):n]) * 3
             else rep_len(P(sim)$lb$t, nT) ## User-defined bounds (recycled if necessary)
             
           }
      )
    )
  
    ## Then, define lower and upper bounds for the second optimizer (nlminb)
    ## Beta
      nlminbUB <-
        if(is.null(P(sim)$ub$b)) rep_len(Inf, nB)
        else DEoptimUB[1:nB] ## User-defined bounds
    
      nlminbLB <-
        if (is.null(P(sim)$lb$b)) rep_len(-Inf, nB)
        else DEoptimLB[1:nB] ## User-defined bounds
  
    ## Theta
      nlminbUB <- c(nlminbUB,
        if (is.null(P(sim)$ub$t)) rep_len(Inf, nT)
        else DEoptimUB[(nB + 1L):n] ## User-defined bounds
      )
    
      nlminbLB <- c(nlminbLB,
        if (is.null(P(sim)$lb$t)) rep_len(-Inf, nT)
        else DEoptimLB[(nB + 1L):n] ## User-defined bounds
      )

  ## Define the log-likelihood function (objective function)
  sim$nll <- parse(text = paste0("-sum(dtappareto(envData[[\"", y, "\"]], lambda=beta, theta=theta, a=", P(sim)$a, ", log=TRUE))"))

  if (P(sim)$nCores > 1) 
  {
    cl <- parallel::makePSOCKcluster(names = P(sim)$nCores)
    on.exit(stopCluster(cl))
    parallel::clusterEvalQ(cl, library("PtProcess"))
  }
  
  ## If starting values are not supplied
  if (is.null(P(sim)$start))
  {
    ## First optimizer, get rough estimates of the parameter values
    ## Use these estimates to compute the order of magnitude of these parameters

      control <- list(itermax = P(sim)$itermax, trace = P(sim)$trace)
      if(P(sim)$nCores > 1) control$cluster <- cl
      
      DEoptimCall <- quote(DEoptim(fn = objfun, lower = DEoptimLB, upper = DEoptimUB, control = do.call("DEoptim.control", control)))
      DEoptimCall[names(formals(objfun)[-1])] <- parse(text = formalArgs(objfun)[-1])
      DEoptimBestMem <- eval(DEoptimCall) %>% `[[` ("optim") %>% `[[` ("bestmem")
      
      ## Update scaling matrix
      diag(sm) <- oom(DEoptimBestMem)

      getRandomStarts <- function(.) pmin(pmax(rnorm(length(DEoptimBestMem),0L,2L)/10 + unname(DEoptimBestMem/oom(DEoptimBestMem)), nlminbLB), nlminbUB)
      start <- c(lapply(1:P(sim)$nTrials, getRandomStarts), list(unname(DEoptimBestMem/oom(DEoptimBestMem))))
  } 
  else 
  {
    start <- if (is.list(P(sim)$start))
    {
      diag(sm) <- lapply(P(sim)$start, oom) %>%
        do.call("rbind", .) %>%
        colMeans
      
      lapply(P(sim)$start, function(x) x / diag(sm))
    } 
    else 
    {
      diag(sm) <- oom(P(sim)$start)
      P(sim)$start / diag(sm)
    }
  }
  
  out <- if (is.list(start)) 
  {
    if (P(sim)$nCores > 1) 
    {
      if (trace) clusterEvalQ(cl, sink(paste0("fireSense_SizeFit_trace.", Sys.getpid())))
      
      out <- clusterApplyLB(cl = cl, x = start, fun = objNlminb, objective = objfun, lower = nlminbLB, upper = nlminbUB, control = c(P(sim)$nlminb.control, list(trace = trace)))
      
      if (trace) clusterEvalQ(cl, sink())
    } 
    else out <- lapply(start, objNlminb, objective = objfun, lower = nlminbLB, upper = nlminbUB, control = c(P(sim)$nlminb.control, list(trace = trace)))
    
    ## Select best minimum amongst all trials
    out[[which.min(sapply(out, "[[", "objective"))]]
  } 
  else objNlminb(start, objfun, nlminbLB, nlminbUB, c(P(sim)$nlminb.control, list(trace = trace)))
  
  ## Compute the standard errors around the estimates
  hess.call <- quote(numDeriv::hessian(func = objfun, x = out$par))
  hess.call[names(formals(objfun)[-1L])] <- parse(text = formalArgs(objfun)[-1L])
  hess <- eval(hess.call)
  se <- suppressWarnings(tryCatch(drop(sqrt(diag(solve(hess))) %*% sm), error = function(e) NA))

  convergence <- TRUE
  
  if (out$convergence) 
  {
    convergence <- FALSE
    convergDiagnostic <- paste0("nlminb optimizer did not converge (", out$message, ")")
    warning(moduleName, "> ", convergDiagnostic, immediate. = TRUE)
  } 
  else if(anyNA(se)) 
  {
    ## Negative values in the Hessian matrix suggest that the algorithm did not converge
    convergence <- FALSE
    convergDiagnostic <- "nlminb optimizer reached relative convergence, saddle point?"
    warning(moduleName, "> ", convergDiagnostic, immediate. = TRUE)
  }
  else 
  {
    convergDiagnostic <- out$message
  }

  ## Parameters scaling: Revert back estimated coefficients to their original scale
  out$par <- drop(out$par %*% sm)

  sim$sizeFitted <- 
    list(formula = P(sim)$formula,
         a = P(sim)$a,
         link = list(beta = lnB, theta = lnT),
         coef = list(beta = setNames(out$par[1:nB], colnames(mmB)),
                     theta = setNames(out$par[(nB + 1L):n], colnames(mmT))),
         se = list(beta = setNames(se[1:nB], colnames(mmB)),
                   theta = setNames(se[(nB + 1L):n], colnames(mmT))),
         LL = -out$objective,
         AIC = 2 * length(out$par) + 2 * out$objective,
         convergence = convergence,
         convergenceDiagnostic = convergDiagnostic
    )
  
  class(sim$sizeFitted) <- "fireSense_SizeFit"
  
  if (!is.na(P(sim)$.runInterval)) # Assumes time only moves forward
    sim <- scheduleEvent(sim, currentTime + P(sim)$.runInterval, moduleName, "run")
  
  invisible(sim)
}


sizeFitSave <- function(sim) 
{
  moduleName <- current(sim)$moduleName
  timeUnit <- timeunit(sim)
  currentTime <- time(sim, timeUnit)
  
  saveRDS(
    sim$sizeFitted, 
    file = file.path(paths(sim)$out, paste0("fireSense_SizeFitted_", timeUnit, currentTime, ".rds"))
  )
  
  if (!is.na(P(sim)$.saveInterval))
    sim <- scheduleEvent(sim, currentTime + P(sim)$.saveInterval, moduleName, "save", .last())
  
  invisible(sim)
}
