#'@name getIPD
#'@title Reconstruct individual patient data (IPD) from Scanned Kaplan-Meier(K-M) curves
#'@description After the raw dataset is processed using the \code{\link{preprocess}} function, we can use the \code{getIPD()} function to reconstruct the IPD.
#'       Here the total number of events (tot.events) is an optional input; and the treatment arm can be arbitrarily assigned to label the
#'       patients' treatment group (Typically, 0 for the control group, and 1 for the treatment group). \cr\cr
#'       The output is the reconstructed IPD in the form of a three-column table (i.e.,time, patient status, and treatment group ID). \cr\cr
#'       In addition, in order to evaluate the accuracy of our reconstruction process, we will calculate the survival probabilities at each read-in time points
#'       based on the reconstructed IPD, then compare them with the corresponding read-in survival probabilities. The test statistics are also included in the
#'       output.
#'
#'@usage getIPD(prep,armID=1,tot.events=NULL)
#'@param prep the class object returned from the \code{preprocess()} function.
#'@param armID the arbitrary lable used as the group indicator for the reconstructed IPD. Typically 0 for the control group and 1 for the treatment group.
#'@param tot.events the total number of events. This may not be available for some published curves, thus this input is optional.
#'@return \code{getIPD()} returns a list object, including four items as follows. \cr \cr
#'        IPD:  the estimated individual patient in a three-column table (i.e.time, status, and treatment group indicator).  \cr \cr
#'        Points:   the data frame shows estimations of parameters at each read-in time points. \cr \cr
#'        riskmat:  the data frame shows index of read-in points within each time interval, as well as the estimated numbers of censored patients and events within each time interval.\cr \cr
#'        kstest:   the test statistics and p value of Kolmogorov-Smirnov test when comparing the distributions of estimated and read-in K-M curves. The null hypothesis is the read-in and estimated survival probabilities are from the same distribution.\cr\cr
#'        precision:   a list shows the root mean squre error(RMSE), mean absolute error and max absolute error which measure the differences between the estimated and read-in survival probabilities. \cr\cr
#'        endpts:  the number of patients remaining at the end of trial.\cr \cr
#'@importFrom dplyr mutate group_by summarise first last
#'@importFrom stats na.omit ks.test
#'@importFrom survival Surv survfit

#'@export
#'@examples
#'
#' # Radiationdata$radio is a dataset exported from ScanIt software ================
#' radio <- Radiationdata$radio
#'
#' # Load time points when the patients numbers =======
#' # at risk reported (i.e. trisk in month) =========
#' trisk <- Radiationdata$trisk
#'
#' # Load the numbers of patients at risk reported (i.e. nrisk) ======
#' # at the time points (trisk) ========
#' nrisk.radio <- Radiationdata$nrisk.radio
#'
#' # Use the trisk and nrisk as input for preprocess and reconstruction ============
#' pre_radio_1 <- preprocess(dat=Radiationdata$radio, trisk=trisk,
#'              nrisk=nrisk.radio,totalpts=NULL,maxy=100)
#' est_radio_1 <- getIPD(prep=pre_radio_1,armID=0,tot.events=NULL)
#'
#' # Output include reconstructed individual patients data =========================
#' head(est_radio_1$IPD)
#'
#' # When trisk and nrisk were not available, then we must input ====================
#' # the initial number of patients   ===============================================
#' pre_radio_2 <- preprocess(dat=Radiationdata$radio, totalpts=213,maxy=100)
#' est_radio_2 <- getIPD(prep=pre_radio_2,armID=0,tot.events=NULL)
#'
#' # Output include reconstructed individual patients data ==========================
#' head(est_radio_2$IPD)

#'@references Guyot P, Ades AE, Ouwens MJ, Welton NJ. Enhanced secondary analysis of survival data: reconstructing the data from published Kaplan-Meier survival curves. BMC Med Res Methodol.2012; 1:9.


getIPD <- function(prep,armID=1,tot.events=NULL){
  interval <- id <- risk <- censor <- event <- estsurv <- surv <- NULL
  ## extract input data from list, and initial vectors----------------------------------------------------------
  dat <- prep[[1]]
  riskmat <- prep[[2]]
  endpts <- prep[[3]]
  ori_dat <- prep[[4]]

  TT <- dat$time
  SS <- dat$surv
  total <- nrow(dat)

  t.risk<-riskmat[,4]
  lower<-riskmat[,2]
  upper<-riskmat[,3]
  n.risk<-riskmat[,5]

  ninterval <- nrow(riskmat)
  ncensor=rep(0,ninterval)
  lasti<-rep(1,ninterval)

  cen <- rep(0,total)
  nhat <- rep(n.risk[1]+1,total+1)
  d <- rep(0,total)
  KM.hat<-rep(1,total)

  ## function to estimate the censor time, and censor number within interval i:evenly distributed
  est_cen <- function(ncen,i,low,upp,TT){
    k <- upp[i]-low[i]+1
    if (k<0) {stop("Riskmat error: upper<lower")}
    if ((k==0) && (ncen<=0)){cen<-0}
    if ((k==0) && (ncen>0)){ cen<-ncen}

    if (k>0) {
      cen <- rep(0,k)
      if (ncen>0){
        cen.t <- rep(0,ncen)
        for (j in 1:ncen){
          cen.t[j]<- TT[low[i]] +
            j*(TT[upp[(i)]]-TT[low[i]])/(ncen+1)
        }
        m=0
        for (j in low[i]:upp[i]){
          m <- m+1
          if (is.na(TT[j+1])) {
            cen[m]<- sum((TT[j]<=cen.t))
          } else{
            cen[m] <- sum((TT[j]<=cen.t) & (cen.t<TT[j+1]))
          }}
        cen[is.na(cen)] <- 0
      }
    }
    ncen <- sum(cen)
    re <- list(cen,ncen)
    invisible(re)
  } # function




  ## calculate the values for interval 1 to (ninterval-1)------------------------------------
  if (ninterval>1)
  {
    for (i in 1:(ninterval-1)){
      #First approximation of number of censored on interval i
      ncensor[i]<- round(n.risk[i]*SS[lower[i+1]]/SS[lower[i]]- n.risk[i+1])

      #Adjust ncensor until converge to estimation n.hat = n.risk at start of interval (i+1)
      while( ((nhat[lower[i+1]]>n.risk[i+1])&&(ncensor[i]<(n.risk[i]-n.risk[i+1]+1)))
             ||((nhat[lower[i+1]]<n.risk[i+1])&&(ncensor[i]>0)) ){

        estcen <- est_cen (ncen=ncensor[i],i=i,low=lower,upp=upper,TT=TT)
        cen[lower[i]:upper[i]] <- estcen[[1]]
        ncensor[i] <- estcen[[2]]

        #Find number of events and patients at risk for each points
        nhat[lower[i]]<-n.risk[i]
        las<-lasti[i]

        # estimate d[k] for interval i
        for (k in lower[i]:upper[i]){
          if (k==1){
            d[k]<-0
            KM.hat[k]<-1
          }
          else {
            if (KM.hat[las]!=0) {d[k]<-round(nhat[k]*(1-(SS[k]/KM.hat[las])))}
            else {d[k] <- 0}
            KM.hat[k]<-KM.hat[las]*(1-(d[k]/nhat[k]))
          }
          nhat[k+1]<-nhat[k]-d[k]-cen[k]
          if (d[k] != 0) las<-k
          if (nhat[k+1]<0) {nhat[k+1] <- 0}
        }  # end for k loop

        # update ncensor for interval i
        ncensor[i]<- ncensor[i]+(nhat[lower[i+1]]-n.risk[i+1])
      } # end while loop

      ## repare for next interval
      n.risk[i+1]<-nhat[lower[i+1]]
      lasti[(i+1)]<-las
    }  # end for loop
  }  # end if ninterval>1



  ## calculate the value for the last interval-----------------------------------------------------------

  ## for the case have more than one intervals, assume the last interval has the same censor rate as before
  if (ninterval>1)
  { ## the number of events happend in the last interval
    if (is.null(tot.events))
       {leftd <- 0}
    else
       {temp <- sum(d[1:upper[ninterval-1]])
       leftd <- tot.events-temp
       leftd <- ifelse(leftd<=0,0,leftd)
    }
    ## patients number at the end of trial
    mm <- ifelse(is.null(endpts),0,endpts)
    ## the number of censoring happened in the last interval
    ncensor[ninterval] <- min(mean(ncensor[1:(ninterval-1)])*(TT[total]-t.risk[ninterval])
                              /(t.risk[ninterval]-t.risk[ninterval-1]),(n.risk[ninterval]-mm-leftd))
  }

  ## for the case only have one interval (no nrisk and trisk information), assume ncensor=0
  if (ninterval==1)
  {if (is.null(tot.events))
  {ncensor[ninterval] <- 0}  ## since no censor information available, assume ncensor=0
    else {ncensor[ninterval] <- n.risk[ninterval]-tot.events}
  }


  ## estimate the censoring for each point
  estcen <- est_cen (ncen=ncensor[ninterval],i=ninterval,low=lower,upp=upper,TT=TT)
  cen[lower[ninterval]:upper[ninterval]] <- estcen[[1]]
  ncensor[ninterval] <- estcen[[2]]


  ## estimate number of events and patients at risk for each point
  nhat[lower[ninterval]]<-n.risk[ninterval]
  las<-lasti[ninterval]
  for (k in lower[ninterval]:upper[ninterval]){
    if (k==1){
      d[k]<-0
      KM.hat[k]<-1
    }
    else {
      if (KM.hat[las]!=0) {d[k]<-round(nhat[k]*(1-(SS[k]/KM.hat[las])))}
      else {d[k] <- 0}
      KM.hat[k]<-KM.hat[las]*(1-(d[k]/nhat[k]))
    }
    nhat[k+1]<-nhat[k]-d[k]-cen[k]
    if (nhat[k+1]<0) {
      nhat[k+1] <- 0
      cen[k] <- nhat[k]-d[k]
    }

    if (d[k] != 0) las<-k
  }  # end for k loop



  ### summarize the estimation of number of censor numbers, event numbers, and patients numbers -------------------
  ### at risk for each interval -----------------------------------------------------------------------------------
  dat <- dat %>%
    mutate(risk=nhat[1:total],censor=cen,event=d)
  options(dplyr.summarise.inform = FALSE)
  riskmat <- dat %>%
    group_by(interval) %>%
    summarise(lower=first(id),upper=last(id),
              trisk=t.risk[interval[1]],
              nrisk=n.risk[interval[1]],
              nrisk.hat=first(risk),
              censor.hat=sum(censor),
              event.hat=sum(event)
              )

  ## reconstruct individual patients records--------------------------------------------------------------------
  ipd=data.frame(time=numeric(0),status=integer(0),treat=integer(0))
  for (i in 1:total){
    if (d[i]>0) {
      time <- rep(dat$time[i],d[i])
      status <- rep(1,d[i])
      treat <- rep(armID,d[i])
      ipd <- rbind(ipd,cbind(time,status,treat))
    }
    if (cen[i]>0) {
      if (i < total) {time <- rep((dat$time[i]+dat$time[i+1])/2,cen[i])}
      if (i == total) {time <- rep(dat$time[i],cen[i])}
      status <- rep(0,cen[i])
      treat <- rep(armID,cen[i])
      ipd <- rbind(ipd,cbind(time,status,treat))
    }
  }  # end for loop

  if (nhat[total+1]>0)
  {
    time <- rep(dat$time[total],nhat[total+1])
    status <- rep(0,nhat[total+1])
    treat <- rep(armID,nhat[total+1])
    ipd <- rbind(ipd,cbind(time,status,treat))
  }


  ## summary features of the estimation ------------------------------------------------------------
  ## function to find survival rate for giving time from a KM curve defined by IPD
  anypoint <- function(a,x,y){
    if (length(x)!=length(y)){
      stop("The lengths of time and surv vectors should be the same!")
    }
    n <- length(x)
    for (i in 1:n)
      if (a>x[i]) {
        i=i+1
      } else if (a==x[i]){
        re=y[i]
        break
      } else {
        re=y[i-1]
        break
      }
    if (a<x[1]) {re=1}
    if (a>x[n]) {re=y[n]}
    return(re)
  }

  ## find survival estimation from ipd for each time points
  fit <- survfit(Surv(ipd$time,ipd$status)~1,data=ipd)
  dat <- dat %>%
    mutate(estsurv=round(unlist(lapply(time,anypoint,x=fit$time,y=fit$surv)),3)) %>%
    mutate(diff=round((estsurv-surv),3))
  dat <- na.omit(dat)

  ## estimate the precision of the estimation
  ## root mean squre
  rmse <- round(sqrt(sum(dat$diff^2)/(NROW(dat)-1)),3)
  ## mean absolute error
  mean_ae <- round(sum(abs(dat$diff))/NROW(dat),3)
  ## max absolute error
  max_ae <- round(max(dat$diff),3)
  precision <- c(rmse,mean_ae,max_ae)
  names(precision) <- c("RMSE","mean_abserror","max_abserror")
  ## compare the read in and estimated survival rates:Kolmogorov-smirnov test
  t1 <- suppressWarnings(ks.test(dat$estsurv,dat$surv))
  t1 <- c(t1[[1]],t1[[2]])
  names(t1) <- c('D',"p-value")

  ## return
  re <- list(IPD=ipd, Points=dat,riskmat=as.data.frame(riskmat),kstest=t1,precision=precision, endpts=endpts)
  class(re)="getKM"
  invisible(re)
}

