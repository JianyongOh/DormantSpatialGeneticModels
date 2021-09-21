#Retrospective model for diversity with seeds in space

#Questions

#Cases: NoSeed, SeedNoMove, SeedMove (Assume that our 2 organisms always "start" as plants)

#Initialise 2 organisms on the 1D Torus
n = 2
L <- 50
#choose n s.t. L divisible by n
organism <- matrix(append(c(1,1),rep(1,n)),nrow = 2,byrow=TRUE)
#m == repetitions
m = 1000
#simulaton variables
u <- 0.1
r <- 1
R <- 10
recolrate <- 1
exrate <- 1
moverate <- 1

#run simulation
#ct == time spent active
#cs == time spent dormant
#instructions to run the simulation m times while resetting variables that need to be reset
start_time <- Sys.time()
for (k in 1:17){
  ct <- c()
  cs <- c()
  for (i in 1:m){
    organism <- matrix(append(c(1,k),rep(1,n)),nrow = 2,byrow=TRUE)
    temptime = 0
    T = 0
    S = 0
    t <- 0
    CommonAncestorFound <- FALSE
    while (!CommonAncestorFound){
      InsideBall <- rep(FALSE,n)
      
      Xplant <- positions_P(organism)
      lnum_P <- lnumber_P(Xplant)
      pnum_P <- rep(0,n)
      prob_P <- 0
      if (length(Xplant)>0){
        for (j in 1:n){
          pnum_P[j] <- sum(lnum_P == j)/length(Xplant)
          prob_P <- prob_P + pnum_P[j]*(1-pbinom(0,j,u))
        }
      }
      rate_P <- prob_P*(length(Xplant)/L)*recolrate
      
      Xexchange <- positions_E(organism)
      lnum_E <- lnumber_E(Xexchange)
      pnum_E <- rep(0,n)
      prob_E <- 0
      for (j in 1:n){
        pnum_E[j] <- sum(lnum_E == j)/length(Xexchange)
        prob_E <- prob_E + pnum_E[j]*(1-pbinom(0,j,u))
      }
      rate_E <- prob_E*(length(Xexchange)/L)*exrate
      
      Xseed <- positions_S(organism)
      lnum_S <- lnumber_S(Xseed)
      pnum_S <- rep(0,n)
      prob_S <- 0
      if (length(Xseed)>0){
        for (j in 1:n){
          pnum_S[j] <- sum(lnum_S == j)/length(Xseed)
          prob_S <- prob_S + pnum_S[j]*(1-pbinom(0,j,u))
        }
      }
      rate_S <- prob_S*(length(Xseed)/L)*moverate
      
      rate <- rate_P + rate_E + rate_S
      t <- t + rexp(1, rate)
      #probability that a given event is a extinction recolonisation event
      if (runif(1)<(rate_P/rate)){
        #AnB vector of probabilities that centre of event falls in intersection of exactly j lineages and at least one lineage jumps
        AnB = rep(0,n)
        for (j in 1:n){
          AnB[j] <- pnum_P[j]*(1-pbinom(0,j,u))*length(Xplant)/L
        }
        #AgivenB is vector of probabilities that conditional on one lineage jumping, the centre of event intersects j lineages
        AgivenB <- AnB/rate_P
        lcount <- sample(1:n,size = 1,prob = AgivenB)
        id = ceiling(runif(1,0,length(Xplant[lnum_P==lcount])))
        x <- Xplant[lnum_P==lcount][id]
        
        loc <- ceiling(runif(1, x - r - 1, x + r))
        if (loc > L) {
          loc <- loc - L
        }
        if (loc < 1) {
          loc <- loc + L
        }
        
        for (j in 1:n){
          if (organism[1,j] >= x-r & organism[1,j] <= x+r & organism[2,j] == 1){
            InsideBall[j] <- TRUE
          }
          else if (organism[1,j] == 1 & x == L & organism[2,j] == 1){
            InsideBall[j] <- TRUE
          }
          else if (organism[1,j] == L & x == 1 & organism[2,j] == 1){
            InsideBall[j] <- TRUE
          }
        }
        
        condbinom <- rep(0,lcount)
        for (j in 1:lcount){
          condbinom[j] <- dbinom(j,lcount,u)/(1-dbinom(0,lcount,u))
        }
        lpart <- sample(1:lcount,1,prob = condbinom)
        idx <- ((1:n)[InsideBall])[sample.int(length((1:n)[InsideBall]),size=lpart)]
        organism[1,idx] <- loc
        
        if (lpart >= 2){
          T = T + n*(t - temptime)
          CommonAncestorFound <- TRUE
          ct = append(ct,T)
          cs = append(cs,S)
        }
      }
      
      else if (runif(1) < (rate_E)/(rate_E + rate_S)){
        AnB = rep(0,n)
        for (j in 1:n){
          AnB[j] <- pnum_E[j]*(1-pbinom(0,j,u))*length(Xexchange)/L
        }
        AgivenB <- AnB/rate_E
        lcount <- sample(1:n,size = 1,prob = AgivenB)
        id = ceiling(runif(1,0,length(Xexchange[lnum_E==lcount])))
        x <- Xexchange[lnum_E==lcount][id]
        
        for (j in 1:n){
          if (organism[1,j] >= x-r & organism[1,j] <= x+r){
            InsideBall[j] <- TRUE
          }
          else if (organism[1,j] == 1 & x == L){
            InsideBall[j] <- TRUE
          }
          else if (organism[1,j] == L & x == 1){
            InsideBall[j] <- TRUE
          }
        }
        
        condbinom <- rep(0,lcount)
        for (j in 1:lcount){
          condbinom[j] <- dbinom(j,lcount,u)/(1-dbinom(0,lcount,u))
        }
        lpart <- sample(1:lcount,1,prob = condbinom)
        
        for (j in 0:n){
          if (sum(organism[2,]==1)==j){
            S = S + (n-j)*(t - temptime)
            T = T + j*(t-temptime)
          }
        }
        temptime = t
        
        idx <- ((1:n)[InsideBall])[sample.int(length((1:n)[InsideBall]),size=lpart)]
        organism[2,idx] <- 1 - organism[2,idx]
      }
      
      else{
        AnB = rep(0,n)
        for (j in 1:n){
          AnB[j] <- pnum_S[j]*(1-pbinom(0,j,u))*length(Xseed)/L
        }
        AgivenB <- AnB/rate_S
        lcount <- sample(1:n,size = 1,prob = AgivenB)
        id = ceiling(runif(1,0,length(Xseed[lnum_S==lcount])))
        x <- Xseed[lnum_S==lcount][id]
        
        for (j in 1:n){
          if (organism[1,j] >= x-r & organism[1,j] <= x+r & organism[2,j] == 0){
            InsideBall[j] <- TRUE
          }
          else if (organism[1,j] == 1 & x == L & organism[2,j] == 0){
            InsideBall[j] <- TRUE
          }
          else if (organism[1,j] == L & x == 1 & organism[2,j] == 0){
            InsideBall[j] <- TRUE
          }
        }
        
        condbinom <- rep(0,lcount)
        for (j in 1:lcount){
          condbinom[j] <- dbinom(j,lcount,u)/(1-dbinom(0,lcount,u))
        }
        lpart <- sample(1:lcount,1,prob = condbinom)
        idx <- ((1:n)[InsideBall])[sample.int(length((1:n)[InsideBall]),size=lpart)]
        for (j in 1:lpart){  
          loc <- ceiling(runif(1, x - r - 1, x + r))
          if (loc > L) {
            loc <- loc - L
          }
          if (loc < 1) {
            loc <- loc + L
          }
          organism[1,idx[j]] <- loc
        }
      }
    }
  }
  df = data.frame(T = ct, S = cs)
  for (i in 1:17){
    if (k == i){
      path = paste("C:/Users/ohjia/Desktop/MORSE/URSS/Data0.10//timedata",i,".csv",sep="")
      write.csv(df,path)
    }
  }
  
}
end_time <- Sys.time()
end_time - start_time

#Defined Functions

#positions function to find only relevant ball centers
positions_P <- function(organism){
  X <- c()
  for (j in 1:n){
    if (organism[2,j] == 1){
      X <- union(X,(organism[1,j]-1):(organism[1,j]+1))
    }
  }
  if (length(X) > 0){
    for (k in 1:length(X)){
      if (X[k]>L){
        X[k] <- X[k] - L
      }
      if (X[k]<1){
        X[k] <- X[k] + L
      }
    }
  }
  unique(X)
  return(X)
}
Xplant <- positions_P(organism)

#lnumber function to find number of lineages intersect ball at each relevant ball center
lnumber_P <- function(X){
  if (length(X)>0){
    lnum <- rep(0,length(X))
    for (j in 1:length(X)){
      for (k in 1:n){
        if (organism[1,k] >= X[j]-r & organism[1,k] <= X[j]+r & organism[2,k] == 1){
          lnum[j] <- lnum[j]+1
        }
        else if (organism[1,k] == 1 & X[j] == L & organism [2,k] == 1){
          lnum[j] <- lnum[j]+1
        }
        else if (organism[1,k] == L & X[j] == 1 & organism[2,k] == 1){
          lnum[j] <- lnum[j]+1
        }
      }
    }
  }
  else{
    lnum <- c()
  }
  return(lnum)
}
lnum_P <- lnumber_P(Xplant)

#pnum is the vector of probabilities of number of organisms being contained in the ball given that ball contains at least one
pnum_P <- rep(0,n)
#prob is the probability that a plant jumps given that a relevant ball center was chosen
prob_P <- 0
if (length(Xplant)>0){
  for (j in 1:n){
    pnum_P[j] <- sum(lnum_P == j)/length(Xplant)
    prob_P <- prob_P + pnum_P[j]*(1-pbinom(0,j,u))
  }
}
rate_P <- prob_P*(length(Xplant)/L)*recolrate

#similar functions defined for the exchange event
positions_E <- function(organism){
  X <- c()
  for (j in 1:n){
    X <- union(X,(organism[1,j]-1):(organism[1,j]+1))
  }
  for (k in 1:length(X)){
    if (X[k]>L){
      X[k] <- X[k] - L
    }
    if (X[k]<1){
      X[k] <- X[k] + L
    }
  }
  unique(X)
  return(X)
}
Xexchange <- positions_E(organism)

#similar functions defined for the exchange event
lnumber_E <- function(X){
  if (length(X)>0){
    lnum <- rep(0,length(X))
    for (j in 1:length(X)){
      for (k in 1:n){
        if (organism[1,k] >= X[j]-r & organism[1,k] <= X[j]+r){
          lnum[j] <- lnum[j]+1
        }
        else if (organism[1,k] == 1 & X[j] == L){
          lnum[j] <- lnum[j]+1
        }
        else if (organism[1,k] == L & X[j] == 1){
          lnum[j] <- lnum[j]+1
        }
      }
    }
  }
  else{
    lnum <- c()
  }
  return(lnum)
}
lnum_E <- lnumber_E(Xexchange)

pnum_E <- rep(0,n)
prob_E <- 0
for (j in 1:n){
  pnum_E[j] <- sum(lnum_E == j)/length(Xexchange)
  prob_E <- prob_E + pnum_E[j]*(1-pbinom(0,j,u))
}
rate_E <- prob_E*(length(Xexchange)/L)*exrate

#similar function defined for move event
positions_S <- function(organism){
  X <- c()
  for (j in 1:n){
    if (organism[2,j] == 0){
      X <- union(X,(organism[1,j]-1):(organism[1,j]+1))
    }
  }
  if (length(X) > 0){
    for (k in 1:length(X)){
      if (X[k]>L){
        X[k] <- X[k] - L
      }
      if (X[k]<1){
        X[k] <- X[k] + L
      }
    }
  }
  unique(X)
  return(X)
}
Xseed <- positions_S(organism)

#similar function defined for move event
lnumber_S <- function(X){
  if (length(X)>0){
    lnum <- rep(0,length(X))
    for (j in 1:length(X)){
      for (k in 1:n){
        if (organism[1,k] >= X[j]-r & organism[1,k] <= X[j]+r & organism[2,k] == 0){
          lnum[j] <- lnum[j]+1
        }
        else if (organism[1,k] == 1 & X[j] == L & organism [2,k] == 0){
          lnum[j] <- lnum[j]+1
        }
        else if (organism[1,k] == L & X[j] == 1 & organism[2,k] == 0){
          lnum[j] <- lnum[j]+1
        }
      }
    }
  }
  else{
    lnum <- c()
  }
  return(lnum)
}
lnum_S <- lnumber_S(Xseed)

pnum_S <- rep(0,n)
prob_S <- 0
if (length(Xseed)>0){
  for (j in 1:n){
    pnum_S[j] <- sum(lnum_S == j)/length(Xseed)
    prob_S <- prob_S + pnum_S[j]*(1-pbinom(0,j,u))
  }
}
rate_S <- prob_S*(length(Xseed)/L)*moverate

