# Vignettes

The source vignette `netmatchRI.Rmd` is kept here for GitHub-facing
documentation. It is excluded from source-package builds in this local setup
because Pandoc is not available on the current machine.

To re-enable vignette building, install Pandoc, remove the `vignettes/*.Rmd`
entry from `.Rbuildignore`, add `rmarkdown` back to `Suggests`, and add
`VignetteBuilder: knitr` to `DESCRIPTION`.
