### Chow-Liu Algorithm ###

#install.packages("mlbench")
#install.packages("gRbase")
source("http://bioconductor.org/biocLite.R")
biocLite(c("graph", "RBGL", "Rgraphviz"))
install.packages("gRbase", dep = TRUE)
#install.packages("gRbase", repos="http://R-Forge.R-project.org",dependencies = TRUE)
#install.packages("gRbase_1.8-3.4.tar.gz", repos = NULL, type = "source")
install.packages("gRapHD_0.2.6.tar.gz", repos = NULL, type = "source")
library(gRapHD)
library(mlbench)
library(ggplot2)
library(MASS)

X = UScrime
head(UScrime)
### Variables Meaning ###
#LF = Labor force participation rate
#So = Indicator variable for a Southern State (0 or 1)
#Ed = Mean years of schooling
#M.F = number of males per 1000 females.
#U1 = Unemployment rate of urban males 14-24
#U2 = Unemployment rate of urban males 35-39
#Prob = probability of imprisonment
#Ineq = Income inequality
#NW = Number of non-whites per 1000 people.

#optimal forest/maximum likelihood tree for the data
mf = minForest(data=X, stat="BIC")
plot(mf)
#optimal decomposable graph/triangulated graph that minimizes the BIC.
stepf = stepw(mf,X, stat="BIC")
plot(stepf,numIter=1000)

#Has the Markov Properties (pairwise, local, global)
#Ed and M.F are conditionally independent given LF.
#Ineq and NW are conditionally independent given So.

#Interpret the output. 
#Chosen subset (So, Ed, M.F, U1, U2, Prob)

library(dplyr)
str(X)
Xsub = select(X, So, Ed, M.F, U1, U2, Prob)
#best tree
mf2 = minForest(data=Xsub, stat="BIC")
plot(mf2)
#best decomp. graph
stepf2 = stepw(mf2,Xsub, stat="BIC")
plot(stepf2,numIter=1000)

#M.F & Prob are conditionally independent given its neighbors (U1, U2, Ed, So)
#Also, M.F & So are conditionally independent given Ed.

#Comparing the best trees and the best decomposable graphs, the full version 
#definitely displays more branches and interactions between the variables in many directions.
#However, when having an interest in a specific group of variables,
#the reduced version might be more effective in quickly identifying the relationships.


##Using Cross Validation and Stability Approach to Regularization (StARS)##
install.packages("spcov")
library(gRbase)
library(spcov)
source("l12.supplement.r")
data(PimaIndiansDiabetes)
x <- as.matrix(PimaIndiansDiabetes[,1:8])

#best tree
mf3 = minForest(data=x, stat="BIC")
plot(mf3)
#best decomp. graph
stepf3 = stepw(mf3,x, stat="BIC")
plot(stepf3,numIter=1000)
#High-dimensional GGMs with
#a) cross-validation
#b) stability selection with stars.thresh = 0.05
library(glasso)
library(cvTools)
library(gRain)
library(igraph)
library(huge)

##Cross-validation
K <- 5
s <- cvFolds(n=nrow(x), K=K, type="random")
number.candidates <- 20 
number.edges <- vector(length=number.candidates)
m <- vector(length=number.candidates)
sd <- vector(length=number.candidates)
lambda <-  2 ^ (c(1:number.candidates) / (2*number.candidates)) - 1
for (i in 1:number.candidates) 
{
  print(lambda[i])
  l <- cv(K=K, x=x, s=s, lambda=lambda[i])
  number.edges[i] <- l$number.edges
  m[i] <- l$m
  sd[i] <- l$sd
}
plot(lambda, m, xlab=expression(lambda), ylab=expression(f), main="", type="b", col="blue")

# Model selection using CV (Cross Validation)
minimizer <- which.min(m)
c <- cov(x)
d <- glasso(c, rho=lambda[minimizer], penalize.diagonal = FALSE)
g <- d$wi != 0
diag(g) <- FALSE
g1 <- graph.adjacency(g!=0)
layout(1)
plot(g1, vertex.size = 3, vertex.label = NA, edge.arrow.size = .2)


##model Selection Using Stability Selection (StARS)
lasso <- huge(x, method = "glasso")
stars <- huge.select(lasso, criterion = "stars", stars.thresh=0.05)
plot(stars)
