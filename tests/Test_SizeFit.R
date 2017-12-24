library(SpaDES)

## TODO: Expliquer que les coefficients du trace sont scalés et également que les bounds doivent être scalées (ou modif code).

modulePath <- normalizePath("..")

# Define simulation parameters
times <- list(start = 1, end = 1, timeunit = "year")
modules <- list("fireSense_SizeFit")
paths <- list(
  modulePath = modulePath
)

# Examples of model formula
formula_1v <- list(beta = formula(size_fire ~ MDC_06),
                   theta = formula(size_fire ~ MDC_57))

formula_6v <- list(beta = formula(size_fire ~ MDC_06 + hw + dt + wt + ot),
                   theta = formula(size_fire ~ MDC_57 + hw + dt + wt + ot))

# Define module parameters
parameters <- list(
  fireSense_SizeFit = list(
    formula = formula_1v,
    data = "dataFireSense_SizeFit",
    a = 1,
    itermax = 200,
    trace = 10 # Print progress every 10 iterations
  )
)

# Define from where and how data will be loaded in the simList environment
inputs <- data.frame(
  objectName = "dataFireSense_SizeFit",
  file = normalizePath("../inputs/dataFireSense_SizeFit.rds"),
  fun = "readRDS",
  package = "base",
  loadTime = 1
)

# Create the simList
sim <- simInit(
  times = times,
  modules = modules,
  params = parameters,
  paths = paths,
  inputs = inputs
)

sim <- spades(sim)
sim$fireSense_SizeFitted

