#This program will generate the weights for a given matrix of m rows and n columns

# Generate a 100x5 Matrix of random numbers between 0 and 100
X <- replicate(5, runif(100, min = 0, max = 100))
#Make variables that equal the dimensions to use repeatedly in later iterations
# n is rows and m is columns
n <- dim(X)[1]
m <- dim(X)[2]

# We use this because we eventually take a subset of 60% of the data
ell <- .6*n

#The subset is defined by first calculating |x˘i − m˘|, 
# i = 1, . . . , n,
# where xi is the ith row of the scaled predictor matrix
# X˘ (the matrix X where each predictor is scaled 
# to be in the range [0, 1] by subtracting the minimum value 
# and then dividing by the maximum value), m˘ is the vector of 
# coordinatewise medians of X˘ , and ∥ · ∥ is the vector norm.
# Then, the ℓ observations with smallest distances are used 
# as the clean subset, called XS.

#First we define min_c and max_c which find the minimum and
# maximum values, respectively, of each column.
#They should each be vectors with 5 elements
min_c <- apply(X, 2, min)
max_c <- apply(X, 2, max)

#Now I will initialize all of the matrices that I will use
x_new <- matrix(nrow = 100, ncol =5)
x_use <- matrix(nrow = 100, ncol = 1)
h <- matrix(nrow = 100, ncol =1)
w <- matrix(nrow = 100, ncol =1)

#To make the matrix X`: iterate through each column in X, 
# subtract the min of that column(min_c), then divide it by the 
# max of that column(max_c)
for (k in 1:m) {
  x_new[,k] <- (X[,k] -min_c[k])/ max_c[k]
  #Now we subtract the median of each column of x_new to get 
  # (X`-m`). In other words x_new = (X`-m`).
  x_new[,k] <- x_new[,k] - median(x_new[,k])
}

#In order to get Xs we need to take the norm of each row 
# of x_new. Now we will fill up the column vector x_use
# to contain the norms of each row. Since there are 100 rows
# it will be 100X1 matrix. We calculate norm by summing the 
# absolute values of the difference of all the points in the 
# row. temp is the value of the sum of the absolute values 
# of each row so we put that into x_use. 
for (i in 1:n) {
  temp = 0
  for (k in 1:m) {
    temp = temp + abs(x_new[i,k])
  }
  x_use[i,1] <- temp
}
# To get Xs we need to get the rows in the bottom 60% of norms.
# Recall that x_use is a column vector where each row contains
# the value of the norm for its corresponding row in x_new.
# We want the rows of x_new in the bottom 60% of norms.
# We do this by redefining x_use to be X where the rows are 
# in descending order of norms.

# We then define X_s to be the first 60 rows in the newly 
# sorted x_use matrix
x_use <- X[order(x_use), ] 
X_s <- x_use[1:ell, ]

# Now we define h[i] to be the set of leverage values for an 
# observation xi relative to the clean subset
# h[i] = xi(X′S XS)^(−1)xi′. (' means transpose)
# To begin finding h[i], we define guts = (X′S XS)^(−1).
# solve() finds the inverse of a function, and t(x) is the 
# transpose of matrix x.
guts <- solve(t(X_s)%*%(X_s))

# Now we have to multiply guts by row i of xi on the left and 
# row i of xi'on the right. For some reason, you need to use
# %*% to multiply the matrices without any issues. 
for (i in 1:n) {
  h[i,1] <- X[i,]%*%guts%*%as.matrix(X[i,])
}

# Great, now we have successfully generated h[i].
# The weights, w[i] are defined to be w[i] = sqrt(min(h)/h[i])
min_h <- min(h)

for (i in 1:n) {
  w[i,1] <- sqrt(min_h/h[i,1])
}





