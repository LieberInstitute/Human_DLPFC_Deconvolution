library(getopt)
library(sessioninfo)

# Import command-line parameters
spec <- matrix(
    c(
        c("n_donors", "run_num"),
        c("n", "r"),
        rep("1", 2),
        rep("character", 2),
        rep("Add variable description here", 2)
    ),
    ncol = 5
)
opt <- getopt(spec)

print("Using the following parameters:")
print(opt)

session_info()

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/
