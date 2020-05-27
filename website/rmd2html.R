library(rmarkdown)

args = commandArgs(TRUE)
rmd_fn = args[1]

print("########################")
sprintf("# RMD: %s #", rmd_fn)
print("########################")

rmarkdown::render(rmd_fn)

