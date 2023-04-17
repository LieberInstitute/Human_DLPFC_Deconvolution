styler::style_dir(
  path = here::here("code"),
  filetype = c("R", "Rmd"),
  transformers = biocthis::bioc_style()
)

code_dirs <- list.files(here::here("code"))
my_dirs <- code_dirs[grepl("^0[1-9]",code_dirs)]

purrr::walk(my_dirs,  ~styler::style_dir(
  path = here::here("code", .x),
  filetype = c("R", "Rmd"),
  transformers = biocthis::bioc_style()
)
)
