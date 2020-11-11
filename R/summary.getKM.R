#'@name summary.getKM
#'@title Print the summary of the IPD estimation
#'@description   Generate descriptive summary for objects returned by other functions\cr\cr
#' @param object the object returned by other functions.
#' @param ... ignored arguments
#' @details \code{summary()} prints the objects returned by other functions.
#'
#'@importFrom gridExtra grid.table grid.arrange tableGrob
#'@importFrom gridExtra ttheme_default
#'@importFrom grid textGrob gpar

#'@examples
#'
#' # Radiationdata$radio is a dataset exported from ScanIt software ================
#' radio <- Radiationdata$radio
#'
#' # Load time points when the patients number =========
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
#' summary(est_radio_1)
#'
#' # When trisk and nrisk were not available, then we must input ====================
#' # the initial number of patients ===============================================
#' pre_radio_2 <- preprocess(dat=Radiationdata$radio, totalpts=213,maxy=100)
#' est_radio_2 <- getIPD(prep=pre_radio_2,armID=0,tot.events=NULL)
#'
#' # Output include reconstructed individual patients data ==========================
#' head(est_radio_2$IPD)
#' summary(est_radio_2)
#'
#'
#'@references Guyot P, Ades AE, Ouwens MJ, Welton NJ. Enhanced secondary analysis of survival data: reconstructing the data from published Kaplan-Meier survival curves. BMC Med Res Methodol.2012; 1:9.


#'@export


summary.getKM <- function (object,...){
  obj <- object
  cat("\n")
  cat("The function read in ",nrow(obj$Points), "points from the K-M curve, and",
      nrow(obj$riskmat),"numbers of patients at risk.","\n")
  cat("Thus the read-in points are divided into ",
      nrow(obj$riskmat), "time intervals.","\n")
  cat("\n")
  print(as.data.frame(obj$riskmat),row.names=F)
  cat("\n")
  cat("The root-mean-square error between estimated and read-in survival probabilities is",
      format(obj$precision[[1]],digits=4,format='f'),". \n")
  cat("The mean absolute error between estimated and read-in survival probabilities is",
      format(obj$precision[[2]],digits=4,format='f'),". \n")
  cat("The max absolute error between estimated and read-in survival probabilities is",
      format(obj$precision[[3]],digits=4,format='f'),". \n")
  cat("\n")
  cat("The Kolmogorov-Smirnov test:","\n")
  cat("Test statistics D=", format(obj$kstest[[1]],digits=4,format='f'),
      "      p-value=",format(obj$kstest[[2]],digits=4,format='f'),"\n")
  cat("Null hypothesis: distributions of the read-in and estimated survival probabilities are the same.\n")
}






