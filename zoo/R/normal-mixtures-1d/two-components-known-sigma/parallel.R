library(parallel)

# Only works on macOS and Linux
if (.Platform$OS.type != "windows") {
  results <- mclapply(1:10, function(i) {
    Sys.sleep(1)
    return(i * 2)
  }, mc.cores = detectCores())
  print(results)
} else {
  print("mclapply() is not available on Windows.")
}