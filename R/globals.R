# Includes global variables to avoid the issue:
# 'no visible binding for global variable [variable name]'
#
utils::globalVariables(c(
   "tmp_n_t"
#    "Y_name",
#    "all.equations.rhs",
#    "DF.Fit.Predict"
# #   "X",
))

# The following is to avoid the NOTE: "unable to verify current time"
Sys.setenv('_R_CHECK_SYSTEM_CLOCK_' = 0)
# Solution from:
# https://stackoverflow.com/questions/63613301/r-cmd-check-note-unable-to-verify-current-time
