# Load necessary libraries
library(Matrix)  # For sparse matrix operations
library(glmnet)  # For Ridge regression (optional)

# Step 1: Create the Laplacian matrix from the adjacency matrix
# A: adjacency matrix (n_features x n_features)
create_degree = function(x_coor,y_coor,dist_thresh,anno_structure){
  D = matrix(,nrow = length(x_coor),ncol = length(y_coor))
  location_coor = data.frame(Xcoor = x_coor,Ycoor = y_coor)
  i = 1
  while(i <= length(x_coor)){
    cnt = 0
    j = 1
    while(j<=length(y_coor)){
      if( i != j){
        dist_locations = sqrt(sum((location_coor[i,] - location_coor[j,])^2))
        if(dist_locations <= dist_thresh && anno_structure[i] == anno_structure[j]){
          cnt = cnt + 1
        }
      }
      j = j + 1
    }
    D[i,i] = cnt
    i = i + 1
    print(i)
  }
  return(D)
}

create_adj = function(x_coor,y_coor,dist_thresh,anno_structure){
  D = matrix(0,nrow = length(x_coor),ncol = length(y_coor))
  location_coor = data.frame(Xcoor = x_coor,Ycoor = y_coor)
  i = 1
  while(i <= length(x_coor)){
    cnt = 0
    j = 1
    while(j<=length(y_coor)){
      if( i != j){
        dist_locations = sqrt(sum((location_coor[i,] - location_coor[j,])^2))
        if(dist_locations <= dist_thresh && anno_structure[i] == anno_structure[j]){
          D[i,j] = 1
        }
      }
      j = j + 1
    }
    i = i + 1
    print(i)
  }
  return(D)
}


create_laplacian <- function(A) {
  D <- diag(rowSums(A))  # Degree matrix
  L = matrix(0,nrow = dim(A)[1], ncol = dim(A)[2])
  i = 1
  while(i <= dim(A)[1]){
    j = 1
    while (j <= dim(A)[2]) {
      if((i == j) && (D[i,i] != 0) ){
        L[i,j] = 1
      }else if((i != j) && (A[i,j] == 1)){
        L[i,j] = -1 / sqrt(D[i,i] * D[j,j])
      }else{
        L[i,j] = 0
      }
      j = j + 1
    }
    i = i + 1
    print(i)
  }
  # Laplacian matrix
  return(L)
}

#test
# Example adjacency matrix (for 4 features)
#A <- matrix(c(0, 1, 1, 0,
#              1, 0, 1, 0,
#              1, 1, 0, 1,
#              0, 0, 1, 0),
#            nrow=4, byrow=TRUE)

# Compute the Laplacian matrix
L <- create_laplacian(A)

# Step 2: Define the loss function
# X: feature matrix (n_samples x n_features)
# y: target variable (n_samples)
# alpha1: L2 regularization coefficient
# alpha2: Graph regularization coefficient
graph_regularized_loss <- function(w, X, y, alpha1, alpha2, L,weights) {
  residuals <- y - X %*% w  # Residual sum of squares
  rss <- sum(residuals^2)
  #print(rss)

  l1_term <- alpha1 * sum(weights * abs(w))  # L1 regularization
  #print(l1_term)

  graph_term <- alpha2 * sum(t(w) %*% L %*% w)  # Graph regularization


  return(rss + l1_term + graph_term)
}

# Step 3: Fit the model using optimization
# Define a function that wraps the loss function for `optim`
fit_graph_regularized_lm <- function(X, y, alpha1, alpha2, L) {
  n_features <- ncol(X)

  # Objective function to minimize
  objective <- function(w) {
    graph_regularized_loss(w, X, y, alpha1, alpha2, L)
  }

  # Initial coefficients (can start with zeros or Ridge regression result)
  w0 <- rep(0, n_features)


  # Use optim() to minimize the objective function
  fit <- optim(w0, objective)

  return(fit$par)  # Return the fitted coefficients
}
