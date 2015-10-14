#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// // [[Rcpp::export]]
//mat procrustes(mat Z00, mat C, mat Z,int dd){
//    mat Znew = C * Z;
//    mat Z00new = Z00.t();
//    mat X = Z00new * Znew;
//    mat U;
//    mat V;
//    arma::vec s;
//    svd(U,s,V,X);
//    mat Ut = U.t();
//    mat VU = V * Ut;
//    mat zfinal = Znew * VU;
//    return zfinal;    
//}

// [[Rcpp::export]]
arma::mat distMat(int nn, int dd, arma::mat ZZ)
{
    arma::mat dMat(nn,nn,arma::fill::zeros);
    int ii,jj,kk;
    double tmp;
    for(ii = 0 ; ii <= (nn-1) ; ii++){
      for(jj = 0 ; jj <= ii ; jj++){
          tmp = 0.0;
    for(kk = 0 ; kk < dd ; kk++){
    //	double VV = ZZ(ii+kk*nn) - ZZ(jj+kk*nn);
        double VV = ZZ(ii,kk) - ZZ(jj,kk);
 		tmp = tmp + VV*VV;
	}
    dMat(ii,jj) = sqrt(tmp);
    dMat(jj,ii) = sqrt(tmp);
      }
    }
    return(dMat);
}

double logitInverse(double x){
    return 1.0/(1.0 + exp(-1.0 * x));
}


// [[Rcpp::export]]
double FullLogLik(arma::mat YY, arma::mat ZZ, 
    double intercept,int nn,int dd){
	arma::mat dMat = distMat(nn,dd,ZZ);
	double total = 0.0;
	double tmp1,tmp2;
	for(int ii = 1; ii < nn; ii++){ //we exclude diagonal elements
		for(int jj = 0; jj < ii; jj++){
            double vv1 = (intercept - dMat(ii,jj)); 
    		double vv2 = (intercept - dMat(jj,ii));
			tmp1 = logitInverse(vv1);
			tmp2 = logitInverse(vv2);
			if(YY(ii,jj) == 1){
				total = total + log(tmp1);
			}else if(YY(ii,jj) == 0){
				total = total + log(1.0 - tmp1);
			}
			if(YY(jj,ii) == 1){
				total = total + log(tmp2);
			}else if(YY(jj,ii) == 0){
				total = total + log(1.0 - tmp2);
			}
		}		
	}
	return total;
}
    
//// [[Rcpp::export]]
//double FullLogLik(arma::vec YY, arma::mat ZZ, 
//    double intercept,int nn,int dd){
  //  arma::mat dMat(nn,nn);
//	arma::mat dMat = distMat(nn,dd,ZZ);
//	double total = 0.0;
//	double tmp1,tmp2;
//	for(int ii = 1; ii < nn; ii++){ //we exclude diagonal elements
//		for(int jj = 0; jj < ii; jj++){
//            double vv1 = (intercept - dMat(ii,jj)); 
//    		double vv2 = (intercept - dMat(jj,ii));
//			tmp1 = logitInverse(vv1);
//			tmp2 = logitInverse(vv2);
//			if(YY(jj*nn + ii) == 1){
//				total = total + log(tmp1);
//			}else if(YY(jj*nn + ii) == 0){
//				total = total + log(1.0 - tmp1);
//			}
//			if(YY(jj + ii*nn) == 1){
//				total = total + log(tmp2);
//			}else if(YY(jj + ii*nn) == 0){
//				total = total + log(1.0 - tmp2);
//			}
//		}		
//	}
//	return total;
//	//delete[] dMat;
  //  }
    
            
// [[Rcpp::export]]
double likelihoodi(int ii,int dd,int nn, arma::mat Yt,
    arma::mat Zt,double intercept)
{
    ii = ii - 1.0;
    double lliki = 0.0;
    arma::mat dMat =  distMat(nn,dd,Zt);
  //  dMat.print("dMat:");
  //  #compute loglikelihood for the entries of Yt except the diagonal
    for(int jj =0; jj < nn; jj++){
        if(jj == ii){
            lliki = lliki + 0.0;
        }
        if(jj != ii){
            double dij = dMat(ii,jj);
            //Rprintf("dij %f",dij,"\n");
            double v1 = intercept - dij;
            double pij = logitInverse(v1);
            //Rprintf("pij %f",pij);
            if(Yt(ii,jj) == 1){
                lliki = lliki + (Yt(ii,jj))*log(pij);
                }
            if(Yt(ii,jj) == 0){
                lliki = lliki + (1-Yt(ii,jj))*log(1-pij);
                }
            if(Yt(jj,ii) == 1){
                lliki = lliki + (Yt(jj,ii)*log(pij));
                }
            if(Yt(jj,ii) == 0){
                lliki = lliki + (1-Yt(jj,ii))*log(1-pij);
                }
        } }
    return lliki;
}    

// [[Rcpp::export]]
arma::mat varZero(arma::mat Phi,int TT, int dd,arma::mat Zvar){
    arma::mat PhiKr = kron(Phi,Phi); //Kronecker product of Phi
    //PhiKr.print("PhiKr:");
    arma::mat II = arma::eye<arma::mat>(dd*dd,dd*dd);
    arma::mat CC = II - PhiKr;
    //CC.print("CC:");
    arma::mat CCinv = inv(CC);
    //CCinv.print("CCinv:");
    arma::vec ZVarvec(dd*dd);
    int ii = 0;
    for(int kk =0;kk < dd;kk++){
        for(int ll=0;ll< dd;ll++){
            ZVarvec(ii) = Zvar(ll,kk);
            ii = ii + 1;
        }
    }
    //ZVarvec.print("ZVarvec:");
    arma::vec Var0vec = CCinv * ZVarvec;
    //Var0vec.print("Var0vec:");
    arma::mat Var0(&Var0vec[0],dd,dd);
    //Var0.print("Var0");
    return Var0;
}

// [[Rcpp::export]]
double logMVNkern(arma::vec Z, arma::vec mu, arma::mat Sigma){
    arma::mat Zmat(Z);
    arma::mat mumat(mu);
    arma::mat Zdiff = Zmat - mumat;
   // Zdiff.print("Zdiff:");
    arma::mat ZdiffTr = Zdiff.t();
    //ZdiffTr.print("ZdiffTr:");
    arma::mat Sigmainv = inv(Sigma);
    //Sigmainv.print("Sigmainv:");
    arma::mat logkern = (ZdiffTr * Sigmainv * Zdiff);
  //logkern.print("logkern:");
    double logkernVec = -0.5* logkern(0,0);
    return logkernVec;    
}

// [[Rcpp::export]]
double logdmvnorm(arma::vec Z, arma::vec mu, arma::mat Sigma,int dd){
    double logkern = logMVNkern(Z,mu,Sigma);
    double logdetSigma = log(det(Sigma));
    double logNorm = -0.5*dd*log(2*M_PI)-(1.0/2.0)*logdetSigma+logkern;
    return logNorm;    
}


// [[Rcpp::export]]
double Zllik(List ZZ,int TT,int dd,arma::vec nn,
             arma::mat Phi,arma::mat ZVar,arma::vec gList,arma::vec posPrev)
{
    double sumn = 0.0;
    arma::vec mu;
    arma::mat Ztminus1;   
    arma::mat var0 = varZero(Phi,TT,dd,ZVar);
    double llik = 0.0;    
    for(int tt =0; tt < TT; tt ++){
        arma::mat Zt = ZZ(tt);
        arma::mat ZtTr = Zt.t();
        //ZtTr.print("ZtTr:");
        if(tt == 0){
            for(int ii = 0; ii < nn[tt]; ii++){
                arma::vec Z1i = ZtTr.col(ii); 
                //  Z1i.print("Z1i:");
                mu = arma::zeros<arma::vec>(2);
                llik = llik + logdmvnorm(Z1i,mu,var0,dd);
                //    llik.print("llik:");
            }
            Ztminus1 = ZtTr;   
        }
        if(tt > 0){
            if( 0 < tt & tt < (TT-1)){
                for(int ii =0; ii < nn[tt]; ii++){
                    arma::mat Zti = ZtTr.col(ii);
                    if((gList(ii+sumn)==1)|(gList(ii+sumn)==2)){
                        double ij = posPrev(ii+sumn);
                        arma::vec Ztminus1i = Ztminus1.col(ij);
                        mu = Phi * Ztminus1i;
                        llik = llik + logdmvnorm(Zti,mu,ZVar,dd);
                        }else{
                            mu = arma::zeros<arma::vec>(2);
                            llik = llik + logdmvnorm(Zti,mu,var0,dd);
                        }
                    }
            Ztminus1 = ZtTr;
            }
            if(tt==(TT-1)){
                for(int ii =0; ii < nn[tt]; ii++){
                    arma::mat Zti = ZtTr.col(ii);
                    if(gList(ii+sumn)==1){
                        double ij = posPrev(ii+sumn);
                        arma::vec Ztminus1i = Ztminus1.col(ij);
                        mu = Phi * Ztminus1i;
                        llik = llik + logdmvnorm(Zti,mu,ZVar,dd);
                    }else{
                        mu = arma::zeros<arma::vec>(2);
                        llik = llik + logdmvnorm(Zti,mu,var0,dd);
                    }
                }
                Ztminus1 = ZtTr;                
            }
        }
        sumn = sumn + nn[tt];
    }
    return llik;
}


//// [[Rcpp::export]]
//double Zllik(List ZZ,int TT,int dd,arma::vec nn,
//    arma::mat Phi, arma::mat ZVar)
//{
//    arma::vec mu;
//    arma::mat Ztminus1;   
//    arma::mat var0 = varZero(Phi,TT,dd,ZVar);
//    double llik = 0.0;    
//    for(int tt =0; tt < TT; tt ++){
//        arma::mat Zt = ZZ(tt);
//        arma::mat ZtTr = Zt.t();
//        //ZtTr.print("ZtTr:");
//        if(tt == 0){
//            for(int ii = 0; ii < nn[tt]; ii++){
//                arma::vec Z1i = ZtTr.col(ii); 
//              //  Z1i.print("Z1i:");
//                mu = arma::zeros<arma::vec>(2);
//                llik = llik + logdmvnorm(Z1i,mu,var0,dd);
//            //    llik.print("llik:");
//                }
//            Ztminus1 = ZtTr;   
//        }
//        if(tt > 0){
//            for(int ii =0; ii < nn[tt]; ii++){
//                arma::mat Zti = ZtTr.col(ii);
//                arma::vec Ztminus1i = Ztminus1.col(ii);
//                mu = Phi * Ztminus1i;
//                llik = llik + logdmvnorm(Zti,mu,ZVar,dd);
//            }
//            Ztminus1 = ZtTr;
//        }
//    }
//    return llik;
//}
//
//// [[Rcpp::export]]
//List Zupdate1(arma::mat Yt,arma::mat Zt,arma::mat ZNext, int TT,
//        double Intercept,int dd,int nn,arma::mat Phi,arma::mat var,
//                   arma::vec llikOld,arma::vec acct,arma::vec tunet)
//{
//        arma::vec Zsmt(dd);
//        arma::vec Znewsm(dd);
//        arma::vec ZsmNext(dd);
//        arma::vec ZsmPrev(dd);
//        arma::vec meanOld1(dd);
//        arma::vec meanOld2(dd);
//        arma::vec meanNew1(dd);
//        arma::vec meanNew2(dd);
//        double prior1;
//        double prior2;
//        double priorNew1;
//        double priorNew2;
//        double llikNew;
//        double logratio;
//        arma::mat var0 = varZero(Phi,TT, dd,var);
//    // print(var0)
//        arma::mat Znewt(&Zt[0],nn,dd);
//        arma::mat Znew = Znewt.t();//matrix to store updated values         
//        arma::mat Ztt = Zt.t();
//        arma::mat ZNextt = ZNext.t();
//         //update for t = 1
//            //update Z_t by row    
//            for(int i=0; i < nn; i ++){
//                Zsmt = Ztt.col(i);
//              //  Zsmt.print("Zsmt: ");
//             //   Zsmt.print();
//                ZsmNext = ZNextt.col(i);
//             //   ZsmNext.print();
//                ZsmPrev = arma::zeros<arma::vec>(dd);
//            //    ZsmPrev.print();
//                meanOld2 = Phi * Zsmt;
//                //prior density at current values
//                prior1 = logdmvnorm(Zsmt,ZsmPrev,var0,dd);
//                prior2 = logdmvnorm(ZsmNext,meanOld2,var,dd);
////                Rprintf("pO1 %f",prior1);
////                Rprintf("pO2 %f",prior2);
//                //propose new vector
//                for(int ww = 0 ; ww < dd ; ww++){ 
//        	         Znewsm(ww) = Zsmt(ww) + tunet(i)*rnorm(1,0.0,1.0)[0];
//                 }	
//            //    Znewsm.print("Znewsm: ");
//                Znew.col(i) = Znewsm;
//         //       Znew.print();
//                //prior density at proposed values
//                meanNew2 = Phi * Znewsm;
//           //     meanNew2.print();
//                priorNew1 =  logdmvnorm(Znewsm,ZsmPrev,var0,dd);
//                priorNew2 = logdmvnorm(ZsmNext,meanNew2,var,dd);
//                llikNew = likelihoodi((i+1),dd,nn,Yt,Znew.t(),Intercept);
//                //compute log acceptance ratio
//                double ll = llikOld(i);
//                logratio = llikNew-ll+priorNew1-prior1+priorNew2 - prior2;
//          //      Rprintf("logratio %f",logratio);
//                if(log(runif(1,0.0,1.0)[0]) < logratio){
//            //            Rprintf("step %f",1.0);
//                        Ztt.col(i) = Znewsm; //move to proposed values
//                        acct(i) = acct(i)+1; //update acceptance
//                        llikOld(i) = llikNew; //update the likelihood matrix
//                    }else{
//                //        Rprintf("step %f",1.0);
//                        Znew.col(i) = Zsmt;
//                        
//                        } //otherwise stay at current state
//            }
//       Zt = Ztt.t();
//       List rslt(3);
//       rslt(0) = Zt;
//       rslt(1) = acct;
//       rslt(2) = llikOld;
//       return(rslt);
//}
//// [[Rcpp::export]]
//List Zupdatet(arma::mat Yt,arma::mat Zt,arma::mat ZNext,arma::mat ZPrev,int TT,
//            double Intercept,int dd,
//            int nn,arma::mat Phi,arma::mat var,
//              arma::vec llikOld,arma::vec acct,arma::vec tunet)
//{
//    
//    arma::vec Zsmt(dd);
//    arma::vec Znewsm(dd);
//    arma::vec ZsmNext(dd);
//    arma::vec ZsmPrev(dd);
//    arma::vec meanOld1(dd);
//    arma::vec meanOld2(dd);
//    arma::vec meanNew1(dd);
//    arma::vec meanNew2(dd);
//    double prior1;
//    double prior2;
//    double priorNew1;
//    double priorNew2;
//    double llikNew;
//    double logratio;
//    arma::mat Znewt(&Zt[0],nn,dd);
//    arma::mat Znew = Znewt.t();//arma::matrix to store updated values         
//    arma::mat Ztt = Zt.t();
//    arma::mat ZPrevt = ZPrev.t();
//    arma::mat ZNextt = ZNext.t();
////    mat Znew(&Zt[0],nn,dd); //matrix to store updated values         
//    for(int j=0; j < nn;j++){
//        Zsmt = Ztt.col(j);
//        ZsmPrev = ZPrevt.col(j);
//        ZsmNext = ZNextt.col(j); 
//        //propose new z_i
//        for(int ww = 0 ; ww < dd ; ww++){ 
//        	         Znewsm(ww) = Zsmt(ww) + tunet(j)*rnorm(1,0.0,1.0)[0];
//                 }	
//        Znew.col(j) = Znewsm;
//        meanOld1 = Phi*ZsmPrev;
//        meanOld2 = Phi*Zsmt;                        
//        //prior density at current value
//        prior1 = logdmvnorm(Zsmt,meanOld1,var,dd);
//        prior2 = logdmvnorm(ZsmNext,meanOld2,var,dd);
//        //priors at new z_i
//        meanNew1 = Phi*ZsmPrev;
//        meanNew2 = Phi*Znewsm;
//        priorNew1 = logdmvnorm(Znewsm,meanNew1,var,dd);
//        priorNew2 = logdmvnorm(ZsmNext,meanNew2,var,dd);            	
//        //compute likelihood
//        llikNew = likelihoodi((j+1),dd,nn,Yt,Znew.t(),Intercept);    	
//        logratio = llikNew-llikOld(j)+priorNew1-prior1+
//            priorNew2 - prior2;	    
//            //logratio
//        if(log(runif(1,0.0,1.0)[0]) < logratio){
//            Ztt.col(j) = Znewsm;
//            acct(j) = acct(j)+1;
//            llikOld(j) = llikNew;
//        }else{
//            Znew.col(j) = Zsmt;
//        }
//    }
//       Zt = Ztt.t();
//       List rslt(3);
//       rslt(0) = Zt;
//       rslt(1) = acct;
//       rslt(2) = llikOld;
//       return(rslt);
//       }
//// [[Rcpp::export]]
//List ZupdateTT(arma::mat Yt,arma::mat Zt, arma::mat ZPrev, int TT,
//        double Intercept,int dd,int nn,arma::mat Phi,arma::mat var,
//        arma::vec llikOld,arma::vec acct,arma::vec tunet)
//{
//    
//    arma::vec Zsmt(dd);
//    arma::vec Znewsm(dd);
//    arma::vec ZsmPrev(dd);
//    arma::vec meanOld1(dd);
//    arma::vec meanNew1(dd);
//    double prior1;
//    double priorNew1;
//    double llikNew;
//    double logratio;
//    arma::mat Znewt(&Zt[0],nn,dd);
//    arma::mat Znew = Znewt.t();//arma::matrix to store updated values         
//    arma::mat Ztt = Zt.t();
//    arma::mat ZPrevt = ZPrev.t();
//    for(int k=0; k<nn;k++){
//                Zsmt = Ztt.col(k);
//                ZsmPrev = ZPrevt.col(k);
//                //propose new value
//                 for(int ww = 0 ; ww < dd ; ww++){ 
//    		         Znewsm(ww) = Zsmt(ww) + tunet(k)*rnorm(1,0.0,1.0)[0];
//                 }			
//              //   Znewsm.print("Znewsm:");
//                Znew.col(k) = Znewsm;
//            //    Znew.col(k).print();
//                //prior at current value
//                meanOld1 = Phi*ZsmPrev;
//                prior1 = logdmvnorm(Zsmt,meanOld1,var,dd);
//                //prior at new values
//                meanNew1 = Phi*ZsmPrev;
//                priorNew1 = logdmvnorm(Znewsm,meanNew1,var,dd);
//                //loglikelihood at new value
//                llikNew = likelihoodi((k+1),dd,nn,Yt,Znew.t(),Intercept);        
//                logratio = llikNew-llikOld(k)+priorNew1-prior1;
//                //compute logratio
//                if(log(runif(1,0.0,1.0)[0]) < logratio){
//                        Ztt.col(k) = Znewsm;
//                        acct(k) = acct(k)+1;
//                        llikOld(k) = llikNew;
//                }else{
//                        Znew.col(k) = Zsmt; }
//            }
//       Zt = Ztt.t();
//       List rslt(3);
//       rslt(0) = Zt;
//       rslt(1) = acct;
//       rslt(2) = llikOld;
//       return(rslt);
//       }

// [[Rcpp::export]]
List Zupdate1(arma::mat Yt,arma::mat Zt,arma::mat ZNext, int TT,
              double Intercept,int dd,int nn,arma::mat Phi,arma::mat var,
              arma::vec llikOld,arma::vec acct,arma::vec tunet, arma::vec gList, 
              arma::vec posNext)
{
    arma::vec Zsmt(dd);
    arma::vec Znewsm(dd);
    arma::vec ZsmNext(dd);
    arma::vec ZsmPrev(dd);
    arma::vec meanOld1(dd);
    arma::vec meanOld2(dd);
    arma::vec meanNew1(dd);
    arma::vec meanNew2(dd);
    double prior1;
    double prior2;
    double priorNew1;
    double priorNew2;
    double llikNew;
    double logratio;
    arma::mat var0 = varZero(Phi,TT, dd,var);
    arma::mat Znewt(&Zt[0],nn,dd);
    arma::mat Znew = Znewt.t();//matrix to store updated values         
    arma::mat Ztt = Zt.t();
    arma::mat ZNextt = ZNext.t();
    //update for t = 1
    //update Z_t by row    
    for(int i=0; i < nn; i ++){
	//Rprintf("I : %d\n",i);
        Zsmt = Ztt.col(i);
        ZsmPrev = arma::zeros<arma::vec>(dd);
        //prior density at current values
        prior1 = logdmvnorm(Zsmt,ZsmPrev,var0,dd);
         //propose new vector
        for(int ww = 0 ; ww < dd ; ww++){ 
            Znewsm(ww) = Zsmt(ww) + tunet(i)*rnorm(1,0.0,1.0)[0];
        }    
        Znew.col(i) = Znewsm;
        //prior density at proposed values
        priorNew1 =  logdmvnorm(Znewsm,ZsmPrev,var0,dd);        
        llikNew = likelihoodi((i+1),dd,nn,Yt,Znew.t(),Intercept);                        
        double ll = llikOld(i);
        if(gList(i) == 1){
            int ij = posNext(i)-1; //to make it compatible with C++ 0 indexing
	//    Rprintf("IJ : %d\n",ij);
            ZsmNext = ZNextt.col(ij);
	//    ZsmNext.print();	            
            meanOld2 = Phi * Zsmt;
            meanNew2 = Phi * Znewsm;            
            prior2 = logdmvnorm(ZsmNext,meanOld2,var,dd);
            priorNew2 = logdmvnorm(ZsmNext,meanNew2,var,dd);
	 //   Rprintf('priorNew2 %f',priorNew2);
            //compute log acceptance ratio            
            logratio = llikNew-ll+priorNew1-prior1+priorNew2-prior2;
        }else(logratio = llikNew-ll+priorNew1-prior1);
        if(log(runif(1,0.0,1.0)[0]) < logratio){
            Ztt.col(i) = Znewsm; //move to proposed values
            acct(i) = acct(i)+1; //update acceptance
            llikOld(i) = llikNew; //update the likelihood matrix
        }else{
            Znew.col(i) = Zsmt;            
        } //otherwise stay at current state
    }
    Zt = Ztt.t();
 //   Zt.print();	
    List rslt(3);
    rslt(0) = Zt;
    rslt(1) = acct;
    rslt(2) = llikOld;
    return(rslt);
}
// [[Rcpp::export]]
List Zupdatet(arma::mat Yt,arma::mat Zt,arma::mat ZNext,arma::mat ZPrev,int TT,
              double Intercept,int dd,
              int nn,arma::mat Phi,arma::mat var,
              arma::vec llikOld,arma::vec acct,arma::vec tunet,arma::vec gList,arma::vec posPrev, arma::vec posNext)
{    
    arma::vec Zsmt(dd);
    arma::vec Znewsm(dd);
    arma::vec ZsmNext(dd);
    arma::vec ZsmPrev(dd);
    arma::vec meanOld1(dd);
    arma::vec meanOld2(dd);
    arma::vec meanNew1(dd);
    arma::vec meanNew2(dd);
    double prior1;
    double prior2;
    double priorNew1;
    double priorNew2;
    double llikNew;
    double logratio;
    arma::mat Znewt(&Zt[0],nn,dd);
    arma::mat Znew = Znewt.t();//matrix to store updated values         
    arma::mat Ztt = Zt.t();
    arma::mat ZPrevt = ZPrev.t();
    arma::mat ZNextt = ZNext.t();
    arma::mat var0 = varZero(Phi,TT, dd,var);
    for(int j=0; j < nn;j++){
        Zsmt = Ztt.col(j);
        for(int ww = 0 ; ww < dd ; ww++){ 
            Znewsm(ww) = Zsmt(ww) + tunet(j)*rnorm(1,0.0,1.0)[0];
        }           
        Znew.col(j) = Znewsm;
        if((gList(j) == 3)|(gList(j)==4)){
            ZsmPrev = arma::zeros<arma::vec>(dd);
            //prior density at current values
            prior1 = logdmvnorm(Zsmt,ZsmPrev,var0,dd);
            //propose new vector
            //prior density at proposed values
            priorNew1 =  logdmvnorm(Znewsm,ZsmPrev,var0,dd); 
        }else{
            int ij = posPrev(j)-1;
            ZsmPrev = ZPrevt.col(ij);                
            //propose new z_i
            for(int ww = 0 ; ww < dd ; ww++){ 
                Znewsm(ww) = Zsmt(ww)+tunet(j)*rnorm(1,0.0,1.0)[0];
            }    
            Znew.col(j) = Znewsm;
            meanOld1 = Phi*ZsmPrev;
            //prior density at current value
            prior1 = logdmvnorm(Zsmt,meanOld1,var,dd);                
            //priors at new z_i
            meanNew1 = Phi*ZsmPrev;        
            priorNew1 = logdmvnorm(Znewsm,meanNew1,var,dd);
        }
         //compute likelihood
        llikNew = likelihoodi((j+1),dd,nn,Yt,Znew.t(),Intercept);    
        if((gList(j) == 1) |(gList(j) == 3)){
            int jk = posNext(j)-1;
            ZsmNext = ZNextt.col(jk);         
            meanOld2 = Phi*Zsmt;    
            prior2 = logdmvnorm(ZsmNext,meanOld2,var,dd);
            meanNew2 = Phi*Znewsm;
            priorNew2 = logdmvnorm(ZsmNext,meanNew2,var,dd);                
            logratio = llikNew-llikOld(j)+priorNew1-prior1+
                priorNew2 - prior2;        
        }else(logratio = llikNew-llikOld(j)+priorNew1-prior1);
        //logratio
        if(log(runif(1,0.0,1.0)[0]) < logratio){
            Ztt.col(j) = Znewsm;
            acct(j) = acct(j)+1;
            llikOld(j) = llikNew;
        }else{
            Znew.col(j) = Zsmt;
        }
    }
    Zt = Ztt.t();
    List rslt(3);
    rslt(0) = Zt;
    rslt(1) = acct;
    rslt(2) = llikOld;
    return(rslt);
}

// [[Rcpp::export]]
List ZupdateTT(arma::mat Yt,arma::mat Zt, arma::mat ZPrev, int TT,
               double Intercept,int dd,int nn,arma::mat Phi,arma::mat var,
               arma::vec llikOld,arma::vec acct,arma::vec tunet,arma::vec gList,arma::vec posPrev)
{    
    arma::vec Zsmt(dd);
    arma::vec Znewsm(dd);
    arma::vec ZsmPrev(dd);
    arma::vec meanOld1(dd);
    arma::vec meanNew1(dd);
    double prior1;
    double priorNew1;
    double llikNew;
    double logratio;
    arma::mat Znewt(&Zt[0],nn,dd);
    arma::mat Znew = Znewt.t();//matrix to store updated values         
    arma::mat Ztt = Zt.t();
    arma::mat ZPrevt = ZPrev.t();
    arma::mat var0 = varZero(Phi,TT, dd,var);
    for(int k=0; k<nn;k++){
        Zsmt = Ztt.col(k);
        //propose new value
        for(int ww = 0 ; ww < dd ; ww++){ 
            Znewsm(ww) = Zsmt(ww)+tunet(k)*rnorm(1,0.0,1.0)[0];
        }    		                
        if(gList(k) == 0){ //if not present at (TT-1)
            ZsmPrev = arma::zeros<arma::vec>(dd);
            //prior density at current values
            prior1 = logdmvnorm(Zsmt,ZsmPrev,var0,dd);
            //propose new vector
            //prior density at proposed values
            priorNew1 =  logdmvnorm(Znewsm,ZsmPrev,var0,dd); 
        }else{
            int ij = posPrev(k)-1;
            ZsmPrev = ZPrevt.col(ij);
            Znew.col(k) = Znewsm;
            //prior at current value
            meanOld1 = Phi*ZsmPrev;
            prior1 = logdmvnorm(Zsmt,meanOld1,var,dd);
            //prior at new values
            meanNew1 = Phi*ZsmPrev;
            priorNew1 = logdmvnorm(Znewsm,meanNew1,var,dd);
        }
        //loglikelihood at new value
        llikNew = likelihoodi((k+1),dd,nn,Yt,Znew.t(),Intercept);        
        logratio = llikNew-llikOld(k)+priorNew1-prior1;
        //compute logratio
        if(log(runif(1,0.0,1.0)[0]) < logratio){
            Ztt.col(k) = Znewsm;
            acct(k) = acct(k)+1;
            llikOld(k) = llikNew;
        }else{
            Znew.col(k) = Zsmt; }
    }
    Zt = Ztt.t();
    List rslt(3);
    rslt(0) = Zt;
    rslt(1) = acct;
    rslt(2) = llikOld;
    return(rslt);
}


// // **   ///////////////////////


/*		Functions for LSM
 */

double logdnormCpp(double x, double mu, double sigma) {
    double d = (x - mu)/sigma;
    double ret = -log(sqrt(2*PI)*sigma) + (-0.5*d);
   return(ret);
}


//Compute logpriors 
double LogpriorBeta(double Beta, double MuBeta,double SigmaBeta)
{
    double sigma = sqrt(SigmaBeta);    
    double val = logdnormCpp(Beta, MuBeta, sigma);
    return val;
}

// [[Rcpp::export]]
List ZupdateLSM(arma::mat Y,arma::mat Z,double Intercept,int dd,
            int nn,arma::mat var,
            arma::vec llikOld,arma::vec acc,arma::vec tune)
{
        arma::vec Zsm(dd);
        arma::vec Znewsm(dd);
        arma::vec ZsmPrev(dd);
        double prior;
        double priorNew;
        double llikNew;
        double logratio;
    // print(var0)
        arma::mat Znew(&Z[0],nn,dd);
        arma::mat Znewt = Znew.t();//matrix to store updated values         
        arma::mat Ztt = Z.t();
         //update for t = 1
            //update Z_t by row    
            for(int i=0; i < nn; i ++){
                Zsm = Ztt.col(i); 
                ZsmPrev = arma::zeros<arma::vec>(dd);
            //    ZsmPrev.print(); 
                //prior density at current values
                prior = logdmvnorm(Zsm,ZsmPrev,var,dd);
                // Rprintf("pO1 %f",prior1);
//                Rprintf("pO2 %f",prior2);
                //propose new vector
                for(int ww = 0 ; ww < dd ; ww++){ 
                     Znewsm(ww) = Zsm(ww) + tune(i)*rnorm(1,0.0,1.0)[0];
                 }    
            //    Znewsm.print("Znewsm: ");
                Znewt.col(i) = Znewsm;
         //       Znew.print();
                //prior density at proposed values
           //     meanNew2.print();
                priorNew =  logdmvnorm(Znewsm,ZsmPrev,var,dd);
                llikNew = likelihoodi((i+1),dd,nn,Y,Znewt.t(),Intercept);
                //compute log acceptance ratio
                double ll = llikOld(i);
                logratio = llikNew-ll+priorNew-prior;
          //      Rprintf("logratio %f",logratio);
                if(log(runif(1,0.0,1.0)[0]) < logratio){
            //            Rprintf("step %f",1.0);
                        Ztt.col(i) = Znewsm; //move to proposed values
                        acc(i) = acc(i)+1; //update acceptance
                        llikOld(i) = llikNew; //update the likelihood matrix
                    }else{
                //        Rprintf("step %f",1.0);
                        Znewt.col(i) = Zsm;
                        
                        } //otherwise stay at current state
            }
       Z = Ztt.t();
       return List::create(
           _["Z"] = Z,
           _["acc"] = acc,
           _["llikOld"] = llikOld
           );
}


List updateInterceptLSM(arma::mat Y,arma::mat Z,int nn,
    int dd, double Intercept, double mu,double sigmasq,
    double tuneInt,double llik,int acc){
    double intnew;
	double lpIntnew;
	double llikIntnew;
    double lpInt;
    lpInt = LogpriorBeta(Intercept,mu,sigmasq);
    intnew = Intercept + tuneInt*rnorm(1,0.0,1.0)[0];
    lpIntnew = LogpriorBeta(intnew, mu,sigmasq);
    llikIntnew = FullLogLik(Y,Z,intnew,nn,dd);
	double logratio = lpIntnew-lpInt+llikIntnew-llik;
	if(log(runif(1,0.0,1.0)[0]) < logratio){
		Intercept = intnew;
		llik = llikIntnew;
		acc = acc + 1;    
	}
    return List::create(
        _["Intercept"] = Intercept,
        _["acc"] = acc,
        _["llik"] =llik
       );
}


// [[Rcpp::export]]
List MCMCcppLSM(arma::mat Y,arma::mat Z,double Intercept,int nn,int dd,int niter,
        double tuneInt, arma::vec tuneZ, int accInt, arma::vec accZ,
        double MuInt, double VarInt, arma::mat VarZ, double A, double B)
{
    arma::vec llikOldZ(nn);
    arma::vec D(dd);
    List Zupdt;
    List Intupdt;
    arma::vec InterceptFinal(niter);
    List ZFinal(niter);
    List ZVarFinal(niter);
    arma::vec Likelihood(niter);        
    for(int iter =0; iter< niter; iter ++)
    {
        for(int ii = 0; ii < nn; ii ++)
        {
            double ll = likelihoodi((ii+1),dd,nn,Y,Z,Intercept);
            llikOldZ(ii) = ll;
        }
        
        Zupdt = ZupdateLSM(Y,Z,Intercept,dd,nn,VarZ,llikOldZ,accZ,tuneZ);        
        arma::mat Zint = Zupdt["Z"];        
        arma::vec accZint = Zupdt["acc"];
        arma::mat Ztr = Z.t();
        arma::mat ZintTr = Zint.t();
        //this is my way around updating acceptances and Zs        
        for(int aa = 0; aa < nn; aa++){
            accZ(aa) = accZint(aa);
            Ztr.col(aa) = ZintTr.col(aa);
        }
        Z = Ztr.t();
        
 //       accZ.print("accZ");
        double llikAll = FullLogLik(Y,Z,Intercept,nn,dd);
  //     Rprintf("llikAll %f", llikAll);                         

        Intupdt = updateInterceptLSM(Y, Z,nn,dd,Intercept,
                MuInt,VarInt,tuneInt,llikAll,accInt);
        Intercept = Intupdt["Intercept"];
        accInt = Intupdt["acc"];
        llikAll = Intupdt["llik"];        
   		for(int dZ=0; dZ < dd; dZ++){
            for(int nZ=0; nZ < nn; nZ++){
                arma::vec Zcol = Z.col(dZ);
				D(dZ) += (Zcol(nZ))*(Zcol(nZ));
		}   }
          
        for(int vv = 0; vv < dd; vv++){	
            D(vv) = D(vv)/2.0 + B;
            double C = A + nn/2.0;
            VarZ(vv,vv) = 1.0/(rgamma(1,C,((1.0)/D(vv)))[0]);
        }
        //STORE UPDATES
        InterceptFinal(iter) = Intercept;
        ZFinal(iter) = Z;
        Likelihood(iter) = llikAll;
        ZVarFinal(iter) = VarZ;	     
    }    
    return List::create(
        Named("Intercept") = InterceptFinal,
        Named("Z") = ZFinal,
        Named("ZVar") = ZVarFinal,
        Named("Likelihood") = Likelihood,
        Named("accInt") = accInt,
        Named("accZ") = accZ
        );
}


