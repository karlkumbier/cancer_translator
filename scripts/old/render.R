files <- c(
  'cell_line_bioactivity.R',
  'cell_line_classification.R',
  'marker_bioactivity.R',
  'marker_classification.R'
)

for (f in files) {
  rmarkdown::render(f)
}
