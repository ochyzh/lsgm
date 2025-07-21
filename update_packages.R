#These are the steps to updating packages.
# After you edit the data.R file for data descriptions and add any new datasets and functions.
#install.packages("roxygen2")
library(roxygen2)
roxygen2::roxygenise()
devtools::check()
devtools::build()
devtools::load_all()
