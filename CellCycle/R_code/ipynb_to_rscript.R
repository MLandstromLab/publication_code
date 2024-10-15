library(rmarkdown)
library(knitr)
library(xfun)

args = commandArgs(trailingOnly=TRUE)

input <- args[1]

# First convert to .Rmd...
convert_ipynb(input)

# ..then convert .Rmd into .R
purl(with_ext(input, "Rmd"), output = with_ext(input, "R"))

# (Optional) Remove .Rmd file
#file.remove(with_ext(input, "Rmd"))
