require(rmarkdown)
require(utils)

my_markdown_renderer <- function(text) {
  rmd_file_name <- "temp.Rmd"
  yaml_header <- "---
output:
  html_document:
    toc: true
    fig_caption: true
    theme: flatly
---"

  content <- paste0(
    yaml_header,
    "\n",
    "\n",
    "```{r, echo=FALSE}\n",
    text,
    "\n",
    "```\n"
  )

  write(content, rmd_file_name)

  rmarkdown::render(rmd_file_name)
  enc <- utils::URLencode(gsub("Rmd$", "html", rmd_file_name))
  print(enc)
  utils::browseURL(paste0("file://", enc))
}


text <- 'd1 = list(
  a = 42,
  b = "foo",
  c = c("elem1","elem2","elem3"))

  d1'
