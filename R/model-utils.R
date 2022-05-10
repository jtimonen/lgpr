# Colorize string
colorize_string <- function(x, col) {
  if (interactive()) {
    x <- paste0(col, x, "\u001b[0m")
  }
  x
}

# Number string
number_string <- function(x) {
  col <- "\u001b[34;1m" # bold blue
  colorize_string(x, col)
}

# Variable name string
variable_string <- function(x) {
  col <- "\u001b[33m" # orange
  colorize_string(x, col)
}

# Parameter name string
parameter_string <- function(x) {
  col <- "\u001b[31m" # red
  colorize_string(x, col)
}

# A small number
small_number <- function() {
  1e-15
}
