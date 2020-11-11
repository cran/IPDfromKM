#'@name preprocess
#'@title Preprocess the read-in coordinates
#'@description Preprocess the raw coordinates into an appropriate format for reconstruct IPD. Returns include the clean dataset and a table displaying the index of read-in points within each time interval.\cr\cr
#'
#'@details The \code{preprocess()} function process the coordinates dataset extrated from a  published K-M curve using \code{\link{getpoints}} function, or software such as \href{https://www.digitizeit.de}{DigitizeIt} or \href{https://www.amsterchem.com/scanit.html}{ScanIt}.  \cr
#'         In most of published Kaplan-Meier curves, we can also find several numbers of patients at risk under the x-axis. These numbers at risk, and the time
#'         reported them, should be manually input in the form of vectors (nrisk and trisk). However, when these information is not available, we can leave the "trisk" and
#'         "nrisk" parameter as "NULL". In this case, the initial number of patients "totalpts" should be input.  \cr\cr
#'         Sample dataset can be found in \code{\link{Radiationdata}}.
#'@usage preprocess(dat,trisk=NULL,nrisk=NULL,totalpts=NULL,maxy=100)
#'@param dat a two-column dataset with the first column being times, and the second the survival probabilities extracted from a published K-M curve using \code{\link{getpoints}} function, or software such as ScanIt or DigitizeIt.
#'@param trisk a vector containing risk time points (i.e., times points at which the number of patients at risk are reported). This often can be found under the x-axis of a K-M curve. The default value is NULL.
#'@param nrisk a vector containing the numbers of patients at risk reported at the risk time points. This often can be found under the x-axis of a K-M curve. The default value is NULL.
#'@param totalpts the initial number of patients, with a default value of NULL. However, when both trisk and nrisk are NULL, this number is required for the estimation.
#'@param maxy the scale of survival probability. Set maxy=100 when the probabilities are reported in percentages (e.g., 70\%). Set maxy=1 when the probabilities are reported using decimal numbers (e.g, 0.7).

#'@return \code{preprocess()} returns a list object, including four items as follows. \cr
#'
#'        preprocessdat: the two-column(i.e.,time, survival) table after preprocessing \cr \cr
#'        intervalIndex: a table displaying the index of read-in points within each time interval.\cr \cr
#'        endpts: the number of patients remaining at the end of the trial.\cr \cr
#'        inputdat: the read-in dataset.\cr \cr
#'
#'@importFrom dplyr filter mutate arrange slice group_by lag lead desc select distinct summarise row_number
#'@importFrom stats na.omit quantile time wilcox.test isoreg
#'@importFrom dplyr "%>%"
#'
#'@export
#'@examples
#'
#'
#' # Radiationdata$radio is a dataset exported from ScanIt software ================
#' radio <- Radiationdata$radio
#'
#' # Load time points when the patients number =======
#' # at risk reported (i.e. trisk in month) ======
#' trisk <- Radiationdata$trisk
#'
#' # Load the numbers of patients at risk reported (i.e. nrisk) =======
#' # at the time points (trisk) ======
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
#'
#'@references Guyot P, Ades AE, Ouwens MJ, Welton NJ. Enhanced secondary analysis of survival data:
#'reconstructing the data from published Kaplan-Meier survival curves. BMC Med Res Methodol.2012; 1:9.




preprocess <- function(dat,trisk=NULL,nrisk=NULL,totalpts=NULL,maxy=100){
  surv <- survspan <- is_outliers <- takeaway <- n <- interval <- id <- NULL
  ## Formalize the dataset to appropriate values and orders
  dat <- as.data.frame(dat)
  if (NCOL(dat)!=2) {
    stop('The dataset should have exactly two columns.')}
  if (NROW(dat)<5){
    stop('Not enough read in points')
  }
  names(dat) <- c("time","surv")
  ## rescale to survival rate to within [0,1] in case using percentage scale
  if (maxy==100) dat$surv <- dat$surv/100
  inputdat <- dat
  dat <- na.omit(dat)
  dat <- dat[order(dat$time),]

  ## remove outliers to avoid the influence of ridiculous clicks when tracing the curves
  ## detect outlier points using Tukey's fences
  is_outlier <- function(x,k=3.0){
    Q <- quantile(x, probs = c(0.25, 0.75), na.rm = T)
    iqr <- diff(Q)
    (x <= (Q[1]-k*iqr)) | (x >= (Q[2]+k*iqr))
  }

  dat <- dat %>%
      mutate(survspan=abs(surv-lag(surv,default=surv[1]))) %>%
      mutate(is_outliers=is_outlier(x=survspan,k=0.5)) %>%
      mutate(takeaway= is_outliers
                       & (lag(is_outliers,default=F)==F)
                       & lead(is_outliers,default=T))

  subdat <- dat %>% filter(!takeaway) %>% select(time,surv)

  ## As time increasing, make sure survival rates were non-increasing
  for (i in 2:NROW(subdat)){
     if (subdat$surv[i]>subdat$surv[i-1]) {
        subdat$surv[i] <- subdat$surv[i-1]
      }
    }
    subdat <- subdat %>% distinct()
    subdat <- as.data.frame(subdat)

  ## make sure one time point has at most two survival reads
  subdat <- subdat %>%
      group_by(time) %>%
      arrange(desc(surv),.by_group=T) %>%
      slice(c(1,n())) %>%
      distinct()
  subdat <- as.data.frame(subdat)


  ## formalize dataset to make sure that the start point was (0,1), all survival
  ## rates within [0,1], and time always positive and monotonic increasing.
  subdat <- subdat %>% filter(surv>=0&surv<=1&time>=0) %>% arrange(time)
  if (subdat[1,1]!=0|subdat[1,2]!=1) {
    subdat <-  rbind(c(time=0,surv=1),subdat)}
  subdat <- subdat %>% mutate(id=row_number())

  ## formate nrsik, trisk and find the number of intervals as nint
  if (!(is.null(nrisk)) & !(is.null(trisk))){
    ## trim nrisk to have at most one zero at the end, then match the length of nrisk and trisk
    nint <- length(nrisk)
    while ((nrisk[nint]==0) && (nrisk[nint-1]==0)) nint <- nint-1
    nrisk <- nrisk[1:nint]
    trisk <- trisk[1:nint]
  } else if (is.null(totalpts)) {
      stop ("If there is no nrisk vector availabel, at least need to know the total patients' at the beginning.")
  } else {
        nrisk <- totalpts
        trisk <- 0
        nint <- 1
        }


  ## add one column to dataset, indicate which interval each point is belonged to
  locate_interval <- function(x,trisk){
    nint <- length(trisk)
    interval <- 1
    i <- 1
    while (i<nint) {
      i=i+1
      if (x >=trisk[i]) {interval=i}
    }
   return(interval)
  }

  subdat<- subdat %>%
    mutate(interval=unlist(lapply(time,locate_interval,trisk=trisk)))

  ## find out the riskmat table
  riskmat <- subdat %>%
    group_by(interval) %>%
    summarise(lower=first(id),
              upper=last(id),
              t.risk=trisk[interval[1]],
              n.risk=nrisk[interval[1]])

  ## record the endpts: when the curve end at the horizontal segment and there are patients alive at the end of the trial,
  ## record the number of patients left in the endpts (useful for IPD estimaiton); otherwise, null
  if (max(riskmat$t.risk) < max(trisk)) {
    endpts=min(nrisk)
  } else {endpts=NULL}


  ## output the preKM object
  riskmat=as.data.frame(riskmat)
  re=list(preprocessdat=subdat,intervalIndex=riskmat,endpts=endpts,inputdat=inputdat)
  class(re)="preKM"
  invisible(re)
}

