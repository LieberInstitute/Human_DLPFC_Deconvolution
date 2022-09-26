#!/usr/bin/env R

#
# Submitting calculations using bash scripts generated with sgejobs.
#
#

#if (!requireNamespace("remotes", quietly = TRUE))
#  install.packages("remotes")
## Install with:
#remotes::install_github('LieberInstitute/sgejobs')
library(sgejobs)
job_single('jhpce_job', create_logdir = FALSE)
