##########################################################
# Starting script to the module 'SIR models of epidemics'
##########################################################
# implements the basic SIR model, and plots simulation results

###################################
# FUNCTION DEFINITIONS
###################################

###
# sir_model <- function(t, x, parameters)
# Use: calculates the derivatives for the SIR model
# Input:
#      t: time (not used here, because there is no explicit time dependence)
#      x: vector of the current values of all variables (S, I, R)
#      parameters: list containing all parameter values used in the model
# Output:
#      der: vector of derivatives

# To use the lsoda function, the model function has to be a function of t (time),
# x (the vector of the variables) and model parameters.

sir_model <- function(t, x, parameters) {
    # beta and r are not global variables. To be able to refer to their current values,
    # you have to specify that you want to use the value of beta and r from 'parameters'.
    # Similarly, the variables of the model need to be taken from the vector x. This is done by the 'with' function.

    with(as.list(c(parameters, x)), {
        dS <- -beta * S * I
        dI <- +beta * S * I - r * I
        dR <- r * I # Note: because S+I+R=constant, this equation could actually be omitted,
        # and R at any time point could simply be calculated as N-S-I.
        der <- c(dS, dI, dR)
        return(list(der)) # the output must be returned as a list
    }) # end of 'with'
} # end of function definition


###########################
# MAIN PROGRAM
###########################

### LOAD LIBRARIES
# load R library for ordinary differential equation solvers
# you might need to install this package (`install.packages("deSolve")`)
library(deSolve)

### INITIALIZE PARAMETER SETTINGS

parameters <- list(beta = 1e-3, r = 1e-1) # set the parameters of the model
initial_conditions <- c(S = 499, I = 1, R = 0) # set the initial values
timepoints <- seq(0, 100, 0.1) # set the time points for evaluation

# Calculate and print R_0 on the screen
population_size <- sum(initial_conditions)
r_0 <- with(parameters, {
    beta * population_size / r
})
print(paste("R_0 =", r_0), quote = FALSE)

### SIMULATE THE MODEL

## Use lsoda to solve the differential equations numerically.
## (for convenience we convert the output to a data.frame)

sir_simulation <- as.data.frame(
    lsoda(
        y = initial_conditions,
        times = timepoints,
        func = sir_model,
        parms = parameters
    )
)

### PLOT THE OUTPUT

# If you remove the # before pdf(...) and dev.off(), the output will be written in a pdf file,
# in the working directory. If you don't, a window containing your graph will just pop up.

#pdf("startingscript1.pdf")
# par(cex=1.7)
# Plot S as a function of time, in blue, and add the graphs I and R over time
# in red and dark green respectively. Call help(plot) for further details.

plot(x = timepoints, y = sir_simulation$S,
    type = "l", col = "blue", lwd = 3,
    ylim = c(0, sum(initial_conditions)),
    xlab = "time", ylab = "number of individuals"
)
# Note: if you follow the order of arguments in a function call you can omit the argument names.
lines(timepoints, sir_simulation$I,
    type = "l", col = "red", lwd = 3
)
lines(timepoints, sir_simulation$R,
    type = "l", col = "darkgreen", lwd = 3
)
# Add a legend to the graph
legend(70, 400,
    legend = c("S", "I", "R"),
    col = c("blue", "red", "darkgreen"),
    lty = 1, lwd = 2
)
#dev.off()
