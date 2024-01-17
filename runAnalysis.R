#load the quarto library
library(quarto)



#render all the documents in order
qmdFiles <- list.files(pattern="qmd") 
sapply(qmdFiles,quarto_render)

