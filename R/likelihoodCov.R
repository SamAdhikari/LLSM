 likelihoodCov =
   function (Y, Z, intercept, XX, Beta) 
   {
     llik = 0
     Zdist = as.matrix(dist(Z))
     nn = nrow(Y)
     for (ii in 2:nn) {
       for (jj in 1:(ii - 1)) {
         Xij = XX[ii, jj, ]
         Xji = XX[jj,ii, ]
         dij = Zdist[ii, jj]
         pij = logitInverseCov(intercept, dij, Xij, Beta)
         pji = logitInverseCov(intercept,dij,Xji,Beta)
         if (Y[ii, jj] == 1) {
           llik = llik + (Y[ii, jj]) * log(pij)
         }
         if (Y[ii, jj] == 0) {
           llik = llik + (1 - Y[ii, jj]) * log(1 - pij)
         }
         if (Y[jj, ii] == 1) {
           llik = llik + (Y[jj, ii]) * log(pji)
         }
         if (Y[jj, ii] == 0) {
           llik = llik + (1 - Y[jj, ii]) * log(1 - pji)
         }
       }
     }
     return(llik)
   }