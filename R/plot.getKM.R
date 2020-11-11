#'@name plot.getKM
#'@title Graph and compare the K-M curve from reconstructed IPD with the read-in coordinates
#'@description  Graph the survival curve based on the reconstructed IPD, and compare it with the input coordinates.
#'        The output includes three graphs: (1) The estimated K-M curve versus read-in; (2) The estimated numbers of patients
#'        at risk versus reported; and (3) The estimated survival probabilities minus read-in survival probabilities over time.
#'
#'
#' @param x the object returned by other functions.
#' @param ... ignored arguments
#'@return NULL


#'@importFrom survival Surv survfit
#'@importFrom ggplot2 ggplot aes geom_line xlim ylim labs theme geom_point geom_hline
#'@importFrom ggplot2 scale_color_manual element_rect geom_ribbon
#'@importFrom gridExtra grid.arrange
#'@importFrom stats na.omit quantile



#'@examples
#'
#' # Radiationdata$radio is a dataset exported from ScanIt software ================
#' radio <- Radiationdata$radio
#'
#' # Load time points when the patients number ========
#' # at risk reported (i.e. trisk in month) =========
#' trisk <- Radiationdata$trisk
#'
#' # Load the numbers of patients at risk reported (i.e. nrisk) ========
#' # at the time points (trisk) =============
#' nrisk.radio <- Radiationdata$nrisk.radio
#'
#' ##### Use the trisk and nrisk as input ==========
#' pre_radio_1 <- preprocess(dat=Radiationdata$radio, trisk=trisk,nrisk=nrisk.radio,maxy=100)
#' est_radio_1 <- getIPD(prep=pre_radio_1,armID=0,tot.events=NULL)
#' # Output include reconstructed individual patients data
#' head(est_radio_1$IPD)
#' # Plot
#' plot(est_radio_1)
#'
#' ##### When trisk and nrisk were not available, then we must input ========
#' ##### the initial number of patients                              ========
#' pre_radio_2 <- preprocess(dat=Radiationdata$radio, totalpts=213,maxy=100)
#' est_radio_2 <- getIPD(prep=pre_radio_2,armID=0,tot.events=NULL)
#' # Output include reconstructed individual patients data
#' head(est_radio_2$IPD)
#' # Plot
#' plot (est_radio_2)
#'
#'@references Guyot P, Ades AE, Ouwens MJ, Welton NJ. Enhanced secondary analysis of survival data: reconstructing the data from published Kaplan-Meier survival curves. BMC Med Res Methodol.2012; 1:9.

#'@export

plot.getKM <- function (x,...)

{
  surv <- estsurv <- lwr <- upr <- From <- n <- type <- NULL

  obj <- x
  ipd <- obj$IPD
  dat <- obj$Points
  npoints <- NROW(dat)

  ## begin graph

  ### graph 1: Compare KM curves
  xrange <- max(dat$time)*1.10
  rmse <- sqrt(sum(dat$diff^2)/(npoints-1))
  dat <- dat %>%
    mutate(lwr=surv-5*rmse) %>%
    mutate(upr=surv+5*rmse)

  p1 <- ggplot(data=dat,aes(x=time))+
    geom_line(aes(y=surv,color="black"),size=0.5)+
    xlim(0,xrange)+ylim(0,1)+
    labs(title="Compare KM curves",
         subtitle="Shadow area covers 99% CI of the estimation.",
         x="Time",
         y="Survival Probability")+
    geom_point(aes(y=estsurv,color="red"),size=0.5)+
    scale_color_manual(labels = c("Estimated","Read-in"), values = c("black", "red"))+
    geom_ribbon(data=dat,aes(x=time,ymin=lwr,ymax=upr),fill="red",alpha=0.3,linetype=2)+
    theme(
      axis.line.x = element_line(color="darkslategray4", size = .5),
      axis.line.y = element_line(color="darkslategray4", size = .5),
      panel.background=element_blank(),
      axis.text.x=element_text(colour="darkslategray4"))+
    theme(legend.position=c(0.7,0.8),
          legend.text = element_text(size=10),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.title = element_blank())





  ## graph2: comparing estimated risk and true risk
  riskmat <- as.data.frame(obj$riskmat)
  yrange <- riskmat$nrisk[1]*1.10
  if (NROW(riskmat)>1)
  {
    gtab1 <- data.frame(time=dat$time,n=dat$risk,type=rep("Estimated",npoints))
    gtab2 <- data.frame(time=riskmat$trisk,n=riskmat$nrisk,type=rep("Reported",NROW(riskmat)))
    zz <- rbind(gtab1,gtab2)
    p2 <- ggplot(zz, aes(x=time, y=n, color=type)) +
      xlim(0,xrange)+ylim(0,yrange)+
      scale_color_manual(values = c("black", "red")) +
      geom_point(data=zz[zz$type=="Reported",],size=2) +
      geom_line(data=zz[zz$type=="Estimated", ],size = 1) +
      labs(title="Compare numbers of patients at risk",
           x="Time",
           y="Numbers at risk")+
      theme(
           axis.line.x = element_line(color="darkslategray4", size = .5),
           axis.line.y = element_line(color="darkslategray4", size = .5),
           panel.background=element_blank(),
           axis.text.x=element_text(colour="darkslategray4"))+
      theme(legend.position=c(0.7,0.8),
            legend.text = element_text(size=10),
            legend.key = element_rect(colour = NA, fill = NA),
            legend.title = element_blank())


  }

  ## graph3: comparing survival probabilites
  p3 <- ggplot(dat, aes(x=time, y=diff)) +
    geom_point(color="black") + ylim(-0.10,0.10)+ xlim(0,xrange)+
    geom_hline(yintercept=0, color = "red",size=0.5)+
    labs(title="Difference between estimated and read-in suvival probabilities",
         x="Time",
         y="Estimated - Readins")+
    theme(
      axis.line.x = element_line(color="darkslategray4", size = .5),
      axis.line.y = element_line(color="darkslategray4", size = .5),
      panel.background=element_blank(),
      axis.text.x=element_text(colour="darkslategray4"))+
    theme(legend.position='none')



  if (NROW(riskmat)>1) {
    grid.arrange(p1,p2,p3,nrow=2,
                 layout_matrix=rbind(c(1,2),c(3,3)),
                 widths=c(2.7,2.7),heights=c(2.5,1.5))
  } else {
    grid.arrange(p1,p3,nrow=2,ncol=1,
                 layout_matrix=rbind(c(1),c(3)),
                 widths=c(4),heights=c(2.5,2.5))
  }



}

