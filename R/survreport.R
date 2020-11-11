#'@name survreport
#'@title Survival analysis on the reconstructed IPD
#'@description Graph the Kaplan-Meier curves and the cumulative hazard curves for the
#'             reconstructed IPD (from the output of \code{\link{getIPD}} function). Also report
#'             the survival times with confidence intervals for a given vector of survival probabilities,
#'             as well as the landmark survival probabilities of interest.(for example, if set interval=6,
#'             the survival probability will be reported at every six months)
#'
#'@usage survreport(ipd1,ipd2=NULL,arms=1,interval=6,s=c(0.75,0.5,0.25),showplots=TRUE)
#'@param ipd1 a three-column (i.e., time, status, and treatment indicator)table of IPD for treatment 1.
#'@param ipd2 a three-column (i.e., time, status, and treatment indicator)table of IPD for treatment 2.
#'@param arms number of treatment group. Can be either 1 or 2.
#'@param interval  length of the time interval for which the landmark survival probabilities are of interest. The default is at every 6 months.
#'@param s a vector with survival probabilities for which the corresponding survival times are reported. e.g., s=0.5 means that the median survival time is desired.
#'@param showplots indicate if the survival plots are displayed or not in the plot window
#'@return \code{survreport()} returns a list object. \cr
#'
#'@export
#'@importFrom survival Surv survfit coxph
#'@importFrom ggplot2 ggplot aes geom_step labs theme element_line element_blank scale_color_manual
#'@importFrom ggplot2 element_text geom_segment expand_limits element_rect annotate
#'@importFrom gridExtra ttheme_default arrangeGrob
#'@importFrom grid textGrob unit
#'@importFrom stats na.omit quantile
#'
#'@examples
#' ### Get data from the sample dataset=======================
#' radio <- Radiationdata$radio
#' radioplus <- Radiationdata$radioplus
#' trisk <- Radiationdata$trisk
#' nrisk_radio <- Radiationdata$nrisk.radio
#' nrisk_radioplus <- Radiationdata$nrisk.radioplus
#' ### Estimate the IPD for the Radiotherapy treatment group ====================
#' pre_radio <- preprocess(dat=radio, trisk=trisk,nrisk=nrisk_radio,maxy=100)
#' est_radio <- getIPD(prep=pre_radio,armID=0,tot.events=NULL)
#' ### Estimate the IPD for the Radiotherapy plus treatment group ====================
#' pre_radio_plus <- preprocess(dat=radioplus, trisk=trisk,nrisk=nrisk_radioplus,maxy=100)
#' est_radio_plus <- getIPD(prep=pre_radio_plus,armID=1,tot.events=NULL)
#' ### survival report for one arm ===================
#' surv1 <- survreport(ipd1=est_radio$IPD,arms=1,interval=6,s=c(0.75,0.5,0.25),showplots=FALSE)
#' print(surv1)
#' surv1 <- survreport(ipd1=est_radio$IPD,arms=1,interval=10,s=seq(0,1,0.2),showplots=TRUE)
#' print(surv1)
#' ### survival report for two arms ===================
#' surv2 <- survreport(ipd1=est_radio$IPD,ipd2=est_radio_plus$IPD,arms=2,
#'                     interval=8,s=c(0.75,0.5,0.25),showplots=TRUE)
#' print(surv2)
#'
#'@references Guyot P, Ades AE, Ouwens MJ, Welton NJ. Enhanced secondary analysis of survival data: reconstructing the data from published Kaplan-Meier survival curves. BMC Med Res Methodol.2012; 1:9.


survreport<- function(ipd1,ipd2=NULL,arms=1,interval=6, s=c(0.75,0.5,0.25),showplots=TRUE){
  surv <- upperCI <- lowerCI <- cumhaz <- cuml <- cumh <- Treatment <- NULL
  ## check the input data set
  if (is.null(ipd1)) {stop('Need to input the first IPD dataset:')}
  if ((arms==2) & (is.null(ipd2))) {stop('Need to input the second IPD dataset if want to analyze two arms:')}

  if (sum(colnames(ipd1)==c("time","status","treat"))<3)
  {stop('The dataset should have exactly three columns named as "time","status" and "treat".')}


  if  ((!is.null(ipd2)) & (sum(colnames(ipd1)==c("time","status","treat"))<3))
  {stop('The dataset should have exactly three columns named as "time","status" and "treat".')}

  ### for survival probability vector s
  s <- s[s!=0]
  s <- s[s!=1]

  ### Peform KM analysis and Print the estimated quantities for the first arm---------------
  fit1 <- survfit(Surv(time,status)~1,data=ipd1)
  # survival time for certain survival probabilities of interest
  qt1 <- data.frame(q=s,qt = quantile(fit1,probs=(1-s)))
  names(qt1) <- c("surv", "time","0.95LCI","0.95UCI")
  qt1 <- qt1[order(-qt1$surv),]
  new.names <- NULL
  for (i in 1:length(s))
  {
    if (is.na(qt1$time[i])) {qt1[i,3] <- qt1[i,4] <-NA}
    new.names <- c(new.names,paste("survprob = ",qt1$surv[i],sep=""))
  }
  rownames(qt1) <- new.names
  qt1=round(qt1,4)

  # survival probablities at some certain time points
  nn <- floor(ipd1$time[nrow(ipd1)]/interval)
  sur <- low <- upp <- se <- new.names <- NULL
  for (i in 1:nn)
  {
    sur <- c(sur,summary(fit1,times=(i*interval))$surv)
    se <- c(se,summary(fit1,times=(i*interval))$std.err)
    low=c(low,summary(fit1,times=(i*interval))$lower)
    upp <- c(upp,summary(fit1,times=(i*interval))$upper)
    new.names <- c(new.names,paste("time = ",i*interval,sep=""))
  }
  survtable1=data.frame(survival=sur,se=se,lower=low,upper=upp)
  survtable1=round(survtable1,4)
  colnames(survtable1)=c("Surv","SE","0.95LCI","0.95UCI")
  rownames(survtable1)=new.names

  group1 <- list(survtime=qt1[,2:4],survprob=survtable1)
  group2 <- NULL
  ### Peform KM analysis and Print the estimated quantities for the second arm---------------
  if (!is.null(ipd2)){
    fit2 <- survfit(Surv(time,status)~1,data=ipd2)

    # survival time for certain survival probabilities of interest
    qt2 <- data.frame(q = s,qt = quantile(fit2,probs=(1-s)))
    names(qt2) <- c("surv", "time","0.95LCI","0.95UCI")
    qt2 <- qt2[order(-qt2$surv),]
    new.names <- NULL
    for (i in 1:length(s))
    {
      if (is.na(qt2$time[i])) {qt2[i,3] <- qt2[i,4] <-NA}
      new.names <- c(new.names,paste("survprob = ",qt2$surv[i],sep=""))
    }
    rownames(qt2) <- new.names
    qt2=round(qt2,4)

    # survival probablities at some certain time points
    nn <- floor(ipd2$time[nrow(ipd2)]/interval)
    sur <- low <- upp <- se <- new.names <- NULL
    for (i in 1:nn)
    {
      sur <- c(sur,summary(fit2,times=(i*interval))$surv)
      se <- c(se,summary(fit2,times=(i*interval))$std.err)
      low=c(low,summary(fit2,times=(i*interval))$lower)
      upp <- c(upp,summary(fit2,times=(i*interval))$upper)
      new.names <- c(new.names,paste("time = ",i*interval,sep=""))
    }
    survtable2=data.frame(survival=sur,se=se,lower=low,upper=upp)
    survtable2=round(survtable2,4)
    colnames(survtable2)=c("Surv","SE","0.95LCI","0.95UCI")
    rownames(survtable2)=new.names
    group2 <- list(survtime=qt2[,2:4],survprob=survtable2)
  }


  ### graph KM for one curve
  if (arms==1 && showplots==TRUE){
    df <- data.frame(time=fit1$time,surv=fit1$surv,cumhaz=fit1$cumhaz,
                     upperCI=fit1$upper,lowerCI=fit1$lower,
                     cuml=fit1$cumhaz-1.96*fit1$std.chaz,cumh=fit1$cumhaz+1.96*fit1$std.chaz)
    ## survival function p1

    p1 <- ggplot(df,aes(x = time)) +
      ylim(0,1)+
      geom_step(aes(y = surv),size = 1,color="orange4") +
      #geom_step(aes(y = upperCI),colour = "orange4", linetype = "dashed", size = 0.5) +
      #geom_step(aes(y = lowerCI),colour = "orange4", linetype = "dashed", size = 0.5) +
      geom_ribbon(aes(x=time,ymin=lowerCI,ymax=upperCI),fill="orange4",alpha=0.3,linetype=2)+
      labs(x="time",y="Survival probability",subtitle="Survival function") +
      theme(
        axis.line.x = element_line(color="darkslategray4", size = .5),
        axis.line.y = element_line(color="darkslategray4", size = .5),
        panel.background=element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        axis.text.x=element_text(colour="darkslategray4"))
    qt11 <- qt1[!is.na(qt1$time),]
    qt12 <- qt1[is.na(qt1$time),]

    if (nrow(qt11)>0)
      p1 <- p1+geom_segment(data = qt11, aes(x = time, y = surv, xend = time, yend = 0), lty = 2) +
      geom_segment(data = qt11, aes(x = 0, y =surv, xend = time, yend = surv), lty = 2)

    if (nrow(qt12)>0)
      p1 <- p1+geom_segment(data = qt12, aes(x = 0, y =surv, xend = ipd1$time[nrow(ipd1)],
                                             yend = surv), lty = 2)

    ## cumulative hazard function p2
    p2 <- ggplot(df,aes(x = time)) +
      geom_step(aes(y = cumhaz),size = 1,color="orange4") +
      geom_ribbon(aes(x=time,ymin=cuml,ymax=cumh),fill="orange4",alpha=0.3,linetype=2)+
      labs(x="time",y="Cumulative hazard",subtitle="Cumulative risk") +
      theme(
        axis.line.x = element_line(color="darkslategray4", size = .5),
        axis.line.y = element_line(color="darkslategray4", size = .5),
        panel.background=element_blank(),
        axis.text.x=element_text(colour="darkslategray4"))

    ## survival probabilities and survival time estimations  p3 p4
    p3 <- gridExtra::tableGrob(survtable1,theme = ttheme_default(base_size = 8))
    p4 <- gridExtra::tableGrob(qt1[,2:4],theme = ttheme_default(base_size = 8))
    ## display the graph
    gridExtra::grid.arrange(p1,p2,
                            arrangeGrob(p3,top=textGrob("Survival probabilities estimation",y=unit(0.5,"npc"))),
                            arrangeGrob(p4,top=textGrob("Survival time estimation",y=unit(-1.5,"npc"))),
                            layout_matrix=rbind(c(1,2),c(4,3)),
                            widths=c(2.7,2.7),heights=c(4.0,5.0))

  }


  if (arms==2 && showplots==TRUE){
    dfm <- data.frame(time=c(fit1$time,fit2$time),surv=c(fit1$surv,fit2$surv),
                      cumhaz=c(fit1$cumhaz,fit2$cumhaz),
                      upperCI=c(fit1$upper,fit2$upper),
                      lowerCI=c(fit1$lower,fit2$lower),
                      cuml=c(fit1$cumhaz-1.96*fit1$std.chaz,fit2$cumhaz-1.96*fit2$std.chaz),
                      cumh=c(fit1$cumhaz+1.96*fit1$std.chaz,fit2$cumhaz+1.96*fit2$std.chaz),
                      Treatment=c(rep(1,NROW(fit1$time)),rep(2,NROW(fit2$time))))

    ipdm <- rbind(ipd1,ipd2)
    fitm <- coxph(Surv(time,status)~treat, data=ipdm)
    ## survival function p5
    p5 <- ggplot(dfm,aes(x = time,y=surv,col=factor(Treatment))) +
      ylim(0,1)+
      geom_step(aes(y = surv),size = 1) +
      geom_ribbon(data=dfm[which(dfm$Treatment==1),],aes(x=time,ymin=lowerCI,ymax=upperCI),
                  fill="darkslategray4",alpha=0.2,linetype = 0)+
      geom_ribbon(data=dfm[which(dfm$Treatment==2),],aes(x=time,ymin=lowerCI,ymax=upperCI),
                  fill="orange4",alpha=0.2,linetype = 0)+
      scale_color_manual(values = c("darkslategray4", "orange4")) +
      annotate(geom = "text", x = 0.1, y = 0.20,
               label = paste("logHR = ",round(fitm$coefficients,3)," (SE = ",round(fitm$var[[1]]^0.5,3),")",sep=""), hjust = 0) +
      annotate(geom = "text", x = 0.1, y = 0.12,
               label = paste("Logrank test p-value = ",round(summary(fitm)$sctest[[3]],3),sep=""), hjust = 0) +
      labs(x="time",y="Survival probability",subtitle="Survival function",col="Treatment") +
      theme(
        axis.line.x = element_line(color="darkslategray4", size = .5),
        axis.line.y = element_line(color="darkslategray4", size = .5),
        panel.background=element_blank(),
        axis.text.x=element_text(colour="darkslategray4"))+
      theme(legend.position=c(0.75,0.9),
            legend.text = element_text(size=10),
            legend.key = element_rect(colour = NA, fill = NA),
            legend.title = element_text(size=10))

    ## cumulative hazard function p6
    p6 <- ggplot(dfm,aes(x = time,y=cumhaz,col=factor(Treatment))) +
      geom_step(aes(y = cumhaz),size = 1) +
      geom_ribbon(data=dfm[which(dfm$Treatment==1),],aes(x=time,ymin=cuml,ymax=cumh),
                  fill="darkslategray4",alpha=0.2,linetype = 0)+
      geom_ribbon(data=dfm[which(dfm$Treatment==2),],aes(x=time,ymin=cuml,ymax=cumh),
                  fill="orange4",alpha=0.2,linetype = 0)+
      scale_color_manual(values = c("darkslategray4", "orange4")) +
      labs(x="time",y="Cumulative hazard",subtitle="Cumulative risk",col="Treatment") +
      theme(
        axis.line.x = element_line(color="darkslategray4", size = .5),
        axis.line.y = element_line(color="darkslategray4", size = .5),
        panel.background=element_blank(),
        axis.text.x=element_text(colour="darkslategray4"))+
      theme(legend.position=c(0.2,0.75),
            legend.text = element_text(size=10),
            legend.key = element_rect(colour = NA, fill = NA),
            legend.title = element_text(size=10))

    ## display the graph
    gridExtra::grid.arrange(p5,p6,nrow=2,ncol=1,
                            widths=8,heights=c(4.5,4.5))

  }
  re <- list(arm1=group1,arm2=group2)
  class(re)="survKM"
  invisible(re)

}



