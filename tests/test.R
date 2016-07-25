library(SpaDES)

mySim <- simInit(
  times = list(start = 1, end = 2, timeunit = "year"),
  modules = list("fireSense_SizeFit"),
  paths = list(modulePath = " # replace with empty string instead"),
  params = list(fireSense_SizeFit = list(
    formula = list(beta = formula(SUP_HA ~ MDC_JUN + HW + DIST + O + WATER),
                   theta = formula(SUP_HA ~ MDC_MJ + HW + DIST + O + WATER)),
    a = 1,
    trace = 5
  )),
  inputs = data.frame(
    files = "Z:/dataFireSense_SizeFit.RData",
    objectName = "dataFireSense_SizeFit",
    functions = "load",
    package = "base",
    loadTime = 0,
    stringsAsFactors = FALSE)
)

spades(mySim, debug = TRUE)

