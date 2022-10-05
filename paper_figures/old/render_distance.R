library(stringr)

cell.lines <- c(
  "786-0",
  "A549",   
  "ALS-WT",
  "DU145",
  "HEPG2",
  "OVCAR4"
)

for (cell.line in cell.lines) {
  rmarkdown::render('distance.R', output_file=str_c(cell.line, '.html'))
}
