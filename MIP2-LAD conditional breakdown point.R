# This program calculates the breakdown point of the LAD 
# regression of the starsdata. We use the algorithm 
# of GiloniPadbergSIAM called MIP2 on pg.12
library('gurobi')
library("Matrix")
library("slam")
# Create list called model to be used later when solving lp
model <- list()
# x2 is a vector containing all of the x values in the starsdata
x2 <- c(4.37, 4.56, 4.26, 4.56, 4.3, 4.46, 3.84, 4.57, 4.26, 4.37, 3.49, 4.43, 4.48, 4.01, 4.29, 4.42, 4.23, 4.42,4.23, 3.49, 4.29, 4.29, 4.42, 4.49, 4.38, 4.42, 4.29, 4.38, 4.22, 3.48, 4.38, 4.56, 4.45, 3.49, 4.23, 4.62, 4.53, 4.45, 4.53, 4.43, 4.38, 4.45, 4.5, 4.45, 4.55, 4.45, 4.42 )

# Define the variable count to be length of x2 because 
# we will use it a lot later throughout the program.
count <- length(x2)

# Now we will slowly start putting together the pieces of 
# the constraint matrix. We will first set up xnew to 
# become a 47x4 matrix with column of 1s, column x2, negative
# column of 1s, negative column x2.
x1 <- rep(1,count)


# We make xnew a 47(number of rows) by 1 matrix filled with 1s.
xnew <- matrix(x1, count,1)


# Now redefine xnew to be column of 1s followed by a column 
# of stars data(x2), then negative column of 1s and negative 
# column of star data. xnew will now be a 47x4 matrix.
xnew <- cbind(xnew, x2, -x1, -x2)

# Now we extend the the constraint matrix to contain the 
# coefficients of the other variables in the constraint matrix
# in the algorithm. Recall, diag(n) where n is a number
# creates an identity matrix of size nxn. Also, 
# diag( x = 0, count) creates a 47x47 zero matrix
# (since count = 47).
xnew <- cbind(xnew, diag(count), -diag(count), diag(count), -diag(count), diag( x = 0, count), diag(x = 0, count))

# As of now, xnew represents the first complete row of the 
# constraint matrix. Now we will append the other rows beneath 
# it to complete the constraint matrix. 
# See GiloniPadbergSIAM pg. 12 for details of algorithm.
x3 <- cbind(diag(x = 0, nrow = count, ncol = 4), diag(x = 0, count), diag(x=0, count), diag(count), diag(x=0, count), -10000*diag(count), diag(0,count))

x4 <- cbind(diag(x = 0, nrow = count, ncol = 4), diag(x = 0, count), diag(x=0, count), diag(x=0, count), diag(count), diag(x=0, count), 10000*diag(count))

x5 <- cbind(diag(x = 0, nrow = count, ncol = 4), diag(count), diag(count),diag(x = 0, count), diag(x = 0, count), 10000*diag(count), 10000*diag(count) )

x6 <- cbind(diag(x = 0, nrow = count, ncol = 4), diag(x = 0, count), diag(x=0, count), diag(x = 0, count), diag(x=0, count), diag(count), diag(count))

# When you are taking a summation, the coefficients 
# will just be a vector not a matrix.We can easily do this 
# using the rep() function.
x7 <- c(rep(0,4), rep(1, count*2), rep(-1, count*2), rep(0, count*2))

x8 <- c(rep(0,count*2 + 4), rep(1,count*2), rep(0,count*2))
xnew <- rbind(xnew, x3, x4, x5, x6, x7, x8)

# Great, now we have completed the constraint matrix and it
# is called xnew. xnew = constraint matrix

# It's very important to make sparse = TRUE because the 
# matrix has so many 0 values and making it sparse as opposed
# to dense will speed up the calculations exponentially.
model$A <- Matrix(xnew, sparse = TRUE)
# Again, we use rep() here because we're dealing with summations
model$obj <- c(rep(0, count*4 +4), rep(1, count*2))
# modesense is what you want to do the your optimization 
# problem, minimize or maximize. "min" or "max"?
model$modelsense <- "min"
# rhs is the right hand side of the constraint matrix
model$rhs <- c(rep(0,count*3), rep(10000, count), rep(1, count), 0, .5)
#sense is the comparison sign in the constraint matrix
model$sense <- c(rep('=', count), rep('<', count*4 +1), '>')
# vtype is the variable type. Common are B for Binary and 
# C for Continuous.
model$vtype <- c(rep('C', count*4 + 4), rep('B', count *2))


params <- list(OutputFlag=0, IntFeasTol = 1e-9, MIPGAP =0)
result <- gurobi(model, params)
print('Solution:')
print(result$objval)




