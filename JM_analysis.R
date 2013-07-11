#####################################################################################################################################################################################################
## input data
#####################################################################################################################################################################################################
m <- read.csv("../../02 Processed Data/00 For Review/relapse_analysis_set.csv")

#####################################################################################################################################################################################################
## packages required
#####################################################################################################################################################################################################
library(survival)
library(nlme)
library(JM)
library(lattice)

#####################################################################################################################################################################################################
## run joint models for each assay test (around 62)
## longitudinal model with random intercept and slope
## joint model fit using spline-PH-aGH method
## Wald test performed on two coefficients (value and slope) for event process
#####################################################################################################################################################################################################
pops <- unique(m$Assay_test)
out <- NULL; ops <- NULL; asa <- NULL
for(i in pops) {
    dat <- m[m$Assay_test==i&!is.na(m$log2_val),]
    dat$event <- 1-dat$DURCREMF
    dat$months <- dat$days_from_comrem1/30; dat$timetoevent <- dat$DURCREM/30
    dat$months[dat$months>dat$timetoevent] <- dat$timetoevent[dat$months>dat$timetoevent]
    fitLME <- try(lme(log2_val~months*TRT,random=~months|participantId,data=dat),silent=T)
    if(class(fitLME)!="try-error") {
        pp <- predict(fitLME); dat$prdv <- as.numeric(pp)
        fitSURV <- coxph(Surv(timetoevent,event)~TRT,data=dat[!duplicated(dat$participantId),],x=T)
        dForm <- list(fixed=~1+TRT,indFixed=c(2,4),random=~1,indRandom=2)
        fit.JM <- try(jointModel(fitLME,fitSURV,timeVar="months",method="spline-PH-aGH",parameterization="both",derivForm=dForm),silent=T)
        if(class(fit.JM)!="try-error") {
            md <- try(summary(fit.JM),silent=T)
            if(class(md)!="try-error") {
                ops <- c(ops,i)
                out <- rbind(out,md$'CoefTable-Event')
                mda <- anova(fit.JM,proccess="Event"); asa <- rbind(asa,mda$aovTab.T["Assoct(all)",])
            }
        }
    }
}

#####################################################################################################################################################################################################
## summarize test statistics with FDR correction
#####################################################################################################################################################################################################
rownames(asa) <- ops; asa$adjp <- p.adjust(asa$Pr,method="BH")
itc <- out[rownames(out)=="Assoct",]; rownames(itc) <- 1:nrow(itc); itc <- as.data.frame(itc); itc <- cbind(Assay_test=ops,itc)
slp <- out[rownames(out)=="Assoct.s",]; rownames(slp) <- 1:nrow(slp); slp <- as.data.frame(slp); slp <- cbind(Assay_test=ops,slp)

#####################################################################################################################################################################################################
## figures to show individual profile and its estimated trajectory for selected assay tests
#####################################################################################################################################################################################################
pdf("../../04 Analysis Output/01 Figures/sig_after_comrem1.pdf",width=8.5,height=11)
## ANCA
i <- "ANCA"
dat <- m[m$Assay_test==i&!is.na(m$log2_val),]
dat$event <- 1-dat$DURCREMF; dat$months <- dat$days_from_comrem1/30; dat$timetoevent <- dat$DURCREM/30
dat$months[dat$months>dat$timetoevent] <- dat$timetoevent[dat$months>dat$timetoevent]
print(xyplot(log2_val~months|as.factor(event)*TRT,group=participantId,data=dat,type="l",ylab="log2(value)",xlab="months from COMREM1",main=i,layout=c(1,4)),split=c(1,1,2,1),more=T)
fitLME <- try(lme(log2_val~months*TRT,random=~1+months|participantId,data=dat),silent=T)
pp <- predict(fitLME); dat$prdv <- as.numeric(pp)
print(xyplot(prdv~months|as.factor(event)*TRT,group=participantId,data=dat,type="l",xlab="months from COMREM1",ylab="predicted value",main=i,layout=c(1,4)),split=c(2,1,2,1))

## Activated Memory B cells (CD80, CD86)
i <- "Activated Memory B cells (CD80, CD86)"
dat <- m[m$Assay_test==i&!is.na(m$log2_val),]
dat$event <- 1-dat$DURCREMF; dat$months <- dat$days_from_comrem1/30; dat$timetoevent <- dat$DURCREM/30
dat$months[dat$months>dat$timetoevent] <- dat$timetoevent[dat$months>dat$timetoevent]
print(xyplot(log2_val~months|as.factor(event)*TRT,group=participantId,data=dat,type="l",ylab="log2(value)",xlab="months from COMREM1",main=i,layout=c(1,4)),split=c(1,1,2,1),more=T)
fitLME <- try(lme(log2_val~months*TRT,random=~1+months|participantId,data=dat),silent=T)
pp <- predict(fitLME); dat$prdv <- as.numeric(pp)
print(xyplot(prdv~months|as.factor(event)*TRT,group=participantId,data=dat,type="l",xlab="months from COMREM1",ylab="predicted value",main=i,layout=c(1,4)),split=c(2,1,2,1))

## Memory B cells (IgD- IgM-)
i <- "Memory B cells (IgD- IgM-)"
dat <- m[m$Assay_test==i&!is.na(m$log2_val),]
dat$event <- 1-dat$DURCREMF; dat$months <- dat$days_from_comrem1/30; dat$timetoevent <- dat$DURCREM/30
dat$months[dat$months>dat$timetoevent] <- dat$timetoevent[dat$months>dat$timetoevent]
print(xyplot(log2_val~months|as.factor(event)*TRT,group=participantId,data=dat,type="l",ylab="log2(value)",xlab="months from COMREM1",main=i,layout=c(1,4)),split=c(1,1,2,1),more=T)
fitLME <- try(lme(log2_val~months*TRT,random=~1+months|participantId,data=dat),silent=T)
pp <- predict(fitLME); dat$prdv <- as.numeric(pp)
print(xyplot(prdv~months|as.factor(event)*TRT,group=participantId,data=dat,type="l",xlab="months from COMREM1",ylab="predicted value",main=i,layout=c(1,4)),split=c(2,1,2,1))

## Memory B cells (IgD- IgM+)
i <- "Memory B cells (IgD- IgM+)"
dat <- m[m$Assay_test==i&!is.na(m$log2_val),]
dat$event <- 1-dat$DURCREMF; dat$months <- dat$days_from_comrem1/30; dat$timetoevent <- dat$DURCREM/30
dat$months[dat$months>dat$timetoevent] <- dat$timetoevent[dat$months>dat$timetoevent]
print(xyplot(log2_val~months|as.factor(event)*TRT,group=participantId,data=dat,type="l",ylab="log2(value)",xlab="months from COMREM1",main=i,layout=c(1,4)),split=c(1,1,2,1),more=T)
fitLME <- try(lme(log2_val~months*TRT,random=~1+months|participantId,data=dat),silent=T)
pp <- predict(fitLME); dat$prdv <- as.numeric(pp)
print(xyplot(prdv~months|as.factor(event)*TRT,group=participantId,data=dat,type="l",xlab="months from COMREM1",ylab="predicted value",main=i,layout=c(1,4)),split=c(2,1,2,1))
dev.off()

#####################################################################################################################################################################################################
## compute conditional probability of surviving later time for an arbitrary trajectory of marker
#####################################################################################################################################################################################################
set.seed(123456)
i <- "ANCA"
dat <- m[m$Assay_test==i&!is.na(m$log2_val),]
dat$event <- 1-dat$DURCREMF; dat$months <- dat$days_from_comrem1/30; dat$timetoevent <- dat$DURCREM/30
dat$months[dat$months>dat$timetoevent] <- dat$timetoevent[dat$months>dat$timetoevent]
fitLME <- try(lme(log2_val~months*TRT,random=~1+months|participantId,data=dat),silent=T)
fitSURV <- coxph(Surv(timetoevent,event)~TRT,data=dat[!duplicated(dat$participantId),],x=T)
dForm <- list(fixed=~1+TRT,indFixed=c(2,4),random=~1,indRandom=2)
fit.JM <- try(jointModel(fitLME,fitSURV,timeVar="months",method="spline-PH-aGH",parameterization="both",derivForm=dForm),silent=T)

## build a specified trajectory
trj <- NULL
v0 <- median(dat$log2_val[dat$TRT=="CYC"&!is.na(dat$log2_val)&dat$months==0])
trj <- rbind(trj,data.frame(id="CYC_nochg",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nochg",log2_val=v0,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nochg",log2_val=v0,months=6,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_upup",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_upup",log2_val=v0+.5,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_upup",log2_val=v0+1,months=6,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_noup",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_noup",log2_val=v0,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_noup",log2_val=v0+1,months=6,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nodn",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nodn",log2_val=v0,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nodn",log2_val=v0-1,months=6,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_dndn",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_dndn",log2_val=v0-.5,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_dndn",log2_val=v0-1,months=6,TRT="CYC"))
v0 <- median(dat$log2_val[dat$TRT=="RTX"&!is.na(dat$log2_val)&dat$months==0])
trj <- rbind(trj,data.frame(id="RTX_nochg",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nochg",log2_val=v0,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nochg",log2_val=v0,months=6,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_upup",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_upup",log2_val=v0+.5,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_upup",log2_val=v0+1,months=6,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_noup",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_noup",log2_val=v0,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_noup",log2_val=v0+1,months=6,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nodn",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nodn",log2_val=v0,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nodn",log2_val=v0-1,months=6,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_dndn",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_dndn",log2_val=v0-.5,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_dndn",log2_val=v0-1,months=6,TRT="RTX"))
md <- survfitJM(fit.JM,trj,idVar="id",); xrg=c(6,14)

pdf("../../04 Analysis Output/01 Figures/surv_pred_ANCA.pdf",width=6,height=6)
## CYC
plot(1,.5,type="n",xlab="Months",ylab="Predicted survival",main="CYC",xlim=xrg,ylim=c(0,1))
tmp <- md$summaries[[1]]; lines(tmp[,1],tmp[,2],col="black",lwd=2); lines(tmp[,1],tmp[,4],col="black",lty="dashed"); lines(tmp[,1],tmp[,5],col="black",lty="dashed")
tmp <- md$summaries[[2]]; lines(tmp[,1],tmp[,2],col="red",lwd=2); lines(tmp[,1],tmp[,4],col="orange",lty="dashed"); lines(tmp[,1],tmp[,5],col="orange",lty="dashed")
tmp <- md$summaries[[3]]; lines(tmp[,1],tmp[,2],col="orange",lwd=2); lines(tmp[,1],tmp[,4],col="red",lty="dashed"); lines(tmp[,1],tmp[,5],col="red",lty="dashed")
tmp <- md$summaries[[4]]; lines(tmp[,1],tmp[,2],col="green",lwd=2); lines(tmp[,1],tmp[,4],col="green",lty="dashed"); lines(tmp[,1],tmp[,5],col="green",lty="dashed")
tmp <- md$summaries[[5]]; lines(tmp[,1],tmp[,2],col="blue",lwd=2); lines(tmp[,1],tmp[,4],col="blue",lty="dashed"); lines(tmp[,1],tmp[,5],col="blue",lty="dashed")
legend("bottomleft",legend=c("no change","constantly up","stable and up","stable and down","constantly down"),col=c("black","red","orange","green","blue"),lwd=2,lty="solid",bty="n")
## RTX
plot(1,.5,type="n",xlab="Months",ylab="Predicted survival",main="RTX",xlim=xrg,ylim=c(0,1))
tmp <- md$summaries[[6]]; lines(tmp[,1],tmp[,2],col="black",lwd=2); lines(tmp[,1],tmp[,4],col="black",lty="dashed"); lines(tmp[,1],tmp[,5],col="black",lty="dashed")
tmp <- md$summaries[[7]]; lines(tmp[,1],tmp[,2],col="red",lwd=2); lines(tmp[,1],tmp[,4],col="orange",lty="dashed"); lines(tmp[,1],tmp[,5],col="orange",lty="dashed")
tmp <- md$summaries[[8]]; lines(tmp[,1],tmp[,2],col="orange",lwd=2); lines(tmp[,1],tmp[,4],col="red",lty="dashed"); lines(tmp[,1],tmp[,5],col="red",lty="dashed")
tmp <- md$summaries[[9]]; lines(tmp[,1],tmp[,2],col="green",lwd=2); lines(tmp[,1],tmp[,4],col="green",lty="dashed"); lines(tmp[,1],tmp[,5],col="green",lty="dashed")
tmp <- md$summaries[[10]]; lines(tmp[,1],tmp[,2],col="blue",lwd=2); lines(tmp[,1],tmp[,4],col="blue",lty="dashed"); lines(tmp[,1],tmp[,5],col="blue",lty="dashed")
legend("bottomleft",legend=c("no change","constantly up","stable and up","stable and down","constantly down"),col=c("black","red","orange","green","blue"),lwd=2,lty="solid",bty="n")
dev.off()


i <- "Activated Memory B cells (CD80, CD86)"
dat <- m[m$Assay_test==i&!is.na(m$log2_val),]
dat$event <- 1-dat$DURCREMF; dat$months <- dat$days_from_comrem1/30; dat$timetoevent <- dat$DURCREM/30
dat$months[dat$months>dat$timetoevent] <- dat$timetoevent[dat$months>dat$timetoevent]
fitLME <- try(lme(log2_val~months*TRT,random=~1+months|participantId,data=dat),silent=T)
fitSURV <- coxph(Surv(timetoevent,event)~TRT,data=dat[!duplicated(dat$participantId),],x=T)
dForm <- list(fixed=~1+TRT,indFixed=c(2,4),random=~1,indRandom=2)
fit.JM <- try(jointModel(fitLME,fitSURV,timeVar="months",method="spline-PH-aGH",parameterization="both",derivForm=dForm),silent=T)

## build a specified trajectory
trj <- NULL
v0 <- median(dat$log2_val[dat$TRT=="CYC"&!is.na(dat$log2_val)&dat$months==0])
trj <- rbind(trj,data.frame(id="CYC_nochg",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nochg",log2_val=v0,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nochg",log2_val=v0,months=6,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_upup",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_upup",log2_val=v0+.5,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_upup",log2_val=v0+1,months=6,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_noup",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_noup",log2_val=v0,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_noup",log2_val=v0+1,months=6,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nodn",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nodn",log2_val=v0,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nodn",log2_val=v0-1,months=6,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_dndn",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_dndn",log2_val=v0-.5,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_dndn",log2_val=v0-1,months=6,TRT="CYC"))
v0 <- median(dat$log2_val[dat$TRT=="RTX"&!is.na(dat$log2_val)&dat$months==0])
trj <- rbind(trj,data.frame(id="RTX_nochg",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nochg",log2_val=v0,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nochg",log2_val=v0,months=6,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_upup",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_upup",log2_val=v0+.5,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_upup",log2_val=v0+1,months=6,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_noup",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_noup",log2_val=v0,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_noup",log2_val=v0+1,months=6,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nodn",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nodn",log2_val=v0,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nodn",log2_val=v0-1,months=6,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_dndn",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_dndn",log2_val=v0-.5,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_dndn",log2_val=v0-1,months=6,TRT="RTX"))
md <- survfitJM(fit.JM,trj,idVar="id",); xrg=c(6,14)

pdf("../../04 Analysis Output/01 Figures/surv_pred_ActMemBCD80CD86.pdf",width=6,height=6)
## CYC
plot(1,.5,type="n",xlab="Months",ylab="Predicted survival",main="CYC",xlim=xrg,ylim=c(0,1))
tmp <- md$summaries[[1]]; lines(tmp[,1],tmp[,2],col="black",lwd=2); lines(tmp[,1],tmp[,4],col="black",lty="dashed"); lines(tmp[,1],tmp[,5],col="black",lty="dashed")
tmp <- md$summaries[[2]]; lines(tmp[,1],tmp[,2],col="red",lwd=2); lines(tmp[,1],tmp[,4],col="orange",lty="dashed"); lines(tmp[,1],tmp[,5],col="orange",lty="dashed")
tmp <- md$summaries[[3]]; lines(tmp[,1],tmp[,2],col="orange",lwd=2); lines(tmp[,1],tmp[,4],col="red",lty="dashed"); lines(tmp[,1],tmp[,5],col="red",lty="dashed")
tmp <- md$summaries[[4]]; lines(tmp[,1],tmp[,2],col="green",lwd=2); lines(tmp[,1],tmp[,4],col="green",lty="dashed"); lines(tmp[,1],tmp[,5],col="green",lty="dashed")
tmp <- md$summaries[[5]]; lines(tmp[,1],tmp[,2],col="blue",lwd=2); lines(tmp[,1],tmp[,4],col="blue",lty="dashed"); lines(tmp[,1],tmp[,5],col="blue",lty="dashed")
legend("bottomleft",legend=c("no change","constantly up","stable and up","stable and down","constantly down"),col=c("black","red","orange","green","blue"),lwd=2,lty="solid",bty="n")
## RTX
plot(1,.5,type="n",xlab="Months",ylab="Predicted survival",main="RTX",xlim=xrg,ylim=c(0,1))
tmp <- md$summaries[[6]]; lines(tmp[,1],tmp[,2],col="black",lwd=2); lines(tmp[,1],tmp[,4],col="black",lty="dashed"); lines(tmp[,1],tmp[,5],col="black",lty="dashed")
tmp <- md$summaries[[7]]; lines(tmp[,1],tmp[,2],col="red",lwd=2); lines(tmp[,1],tmp[,4],col="orange",lty="dashed"); lines(tmp[,1],tmp[,5],col="orange",lty="dashed")
tmp <- md$summaries[[8]]; lines(tmp[,1],tmp[,2],col="orange",lwd=2); lines(tmp[,1],tmp[,4],col="red",lty="dashed"); lines(tmp[,1],tmp[,5],col="red",lty="dashed")
tmp <- md$summaries[[9]]; lines(tmp[,1],tmp[,2],col="green",lwd=2); lines(tmp[,1],tmp[,4],col="green",lty="dashed"); lines(tmp[,1],tmp[,5],col="green",lty="dashed")
tmp <- md$summaries[[10]]; lines(tmp[,1],tmp[,2],col="blue",lwd=2); lines(tmp[,1],tmp[,4],col="blue",lty="dashed"); lines(tmp[,1],tmp[,5],col="blue",lty="dashed")
legend("bottomleft",legend=c("no change","constantly up","stable and up","stable and down","constantly down"),col=c("black","red","orange","green","blue"),lwd=2,lty="solid",bty="n")
dev.off()


i <- "Memory B cells (IgD- IgM-)"
dat <- m[m$Assay_test==i&!is.na(m$log2_val),]
dat$event <- 1-dat$DURCREMF; dat$months <- dat$days_from_comrem1/30; dat$timetoevent <- dat$DURCREM/30
dat$months[dat$months>dat$timetoevent] <- dat$timetoevent[dat$months>dat$timetoevent]
fitLME <- try(lme(log2_val~months*TRT,random=~1+months|participantId,data=dat),silent=T)
fitSURV <- coxph(Surv(timetoevent,event)~TRT,data=dat[!duplicated(dat$participantId),],x=T)
dForm <- list(fixed=~1+TRT,indFixed=c(2,4),random=~1,indRandom=2)
fit.JM <- try(jointModel(fitLME,fitSURV,timeVar="months",method="spline-PH-aGH",parameterization="both",derivForm=dForm),silent=T)

## build a specified trajectory
trj <- NULL
v0 <- median(dat$log2_val[dat$TRT=="CYC"&!is.na(dat$log2_val)&dat$months==0])
trj <- rbind(trj,data.frame(id="CYC_nochg",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nochg",log2_val=v0,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nochg",log2_val=v0,months=6,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_upup",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_upup",log2_val=v0+.5,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_upup",log2_val=v0+1,months=6,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_noup",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_noup",log2_val=v0,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_noup",log2_val=v0+1,months=6,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nodn",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nodn",log2_val=v0,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nodn",log2_val=v0-1,months=6,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_dndn",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_dndn",log2_val=v0-.5,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_dndn",log2_val=v0-1,months=6,TRT="CYC"))
v0 <- median(dat$log2_val[dat$TRT=="RTX"&!is.na(dat$log2_val)&dat$months==0])
trj <- rbind(trj,data.frame(id="RTX_nochg",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nochg",log2_val=v0,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nochg",log2_val=v0,months=6,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_upup",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_upup",log2_val=v0+.5,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_upup",log2_val=v0+1,months=6,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_noup",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_noup",log2_val=v0,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_noup",log2_val=v0+1,months=6,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nodn",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nodn",log2_val=v0,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nodn",log2_val=v0-1,months=6,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_dndn",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_dndn",log2_val=v0-.5,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_dndn",log2_val=v0-1,months=6,TRT="RTX"))
md <- survfitJM(fit.JM,trj,idVar="id",); xrg=c(6,14)

pdf("../../04 Analysis Output/01 Figures/surv_pred_MemBIgDnIgMn.pdf",width=6,height=6)
## CYC
plot(1,.5,type="n",xlab="Months",ylab="Predicted survival",main="CYC",xlim=xrg,ylim=c(0,1))
tmp <- md$summaries[[1]]; lines(tmp[,1],tmp[,2],col="black",lwd=2); lines(tmp[,1],tmp[,4],col="black",lty="dashed"); lines(tmp[,1],tmp[,5],col="black",lty="dashed")
tmp <- md$summaries[[2]]; lines(tmp[,1],tmp[,2],col="red",lwd=2); lines(tmp[,1],tmp[,4],col="orange",lty="dashed"); lines(tmp[,1],tmp[,5],col="orange",lty="dashed")
tmp <- md$summaries[[3]]; lines(tmp[,1],tmp[,2],col="orange",lwd=2); lines(tmp[,1],tmp[,4],col="red",lty="dashed"); lines(tmp[,1],tmp[,5],col="red",lty="dashed")
tmp <- md$summaries[[4]]; lines(tmp[,1],tmp[,2],col="green",lwd=2); lines(tmp[,1],tmp[,4],col="green",lty="dashed"); lines(tmp[,1],tmp[,5],col="green",lty="dashed")
tmp <- md$summaries[[5]]; lines(tmp[,1],tmp[,2],col="blue",lwd=2); lines(tmp[,1],tmp[,4],col="blue",lty="dashed"); lines(tmp[,1],tmp[,5],col="blue",lty="dashed")
legend("bottomleft",legend=c("no change","constantly up","stable and up","stable and down","constantly down"),col=c("black","red","orange","green","blue"),lwd=2,lty="solid",bty="n")
## RTX
plot(1,.5,type="n",xlab="Months",ylab="Predicted survival",main="RTX",xlim=xrg,ylim=c(0,1))
tmp <- md$summaries[[6]]; lines(tmp[,1],tmp[,2],col="black",lwd=2); lines(tmp[,1],tmp[,4],col="black",lty="dashed"); lines(tmp[,1],tmp[,5],col="black",lty="dashed")
tmp <- md$summaries[[7]]; lines(tmp[,1],tmp[,2],col="red",lwd=2); lines(tmp[,1],tmp[,4],col="orange",lty="dashed"); lines(tmp[,1],tmp[,5],col="orange",lty="dashed")
tmp <- md$summaries[[8]]; lines(tmp[,1],tmp[,2],col="orange",lwd=2); lines(tmp[,1],tmp[,4],col="red",lty="dashed"); lines(tmp[,1],tmp[,5],col="red",lty="dashed")
tmp <- md$summaries[[9]]; lines(tmp[,1],tmp[,2],col="green",lwd=2); lines(tmp[,1],tmp[,4],col="green",lty="dashed"); lines(tmp[,1],tmp[,5],col="green",lty="dashed")
tmp <- md$summaries[[10]]; lines(tmp[,1],tmp[,2],col="blue",lwd=2); lines(tmp[,1],tmp[,4],col="blue",lty="dashed"); lines(tmp[,1],tmp[,5],col="blue",lty="dashed")
legend("bottomleft",legend=c("no change","constantly up","stable and up","stable and down","constantly down"),col=c("black","red","orange","green","blue"),lwd=2,lty="solid",bty="n")
dev.off()


i <- "Memory B cells (IgD- IgM+)"
dat <- m[m$Assay_test==i&!is.na(m$log2_val),]
dat$event <- 1-dat$DURCREMF; dat$months <- dat$days_from_comrem1/30; dat$timetoevent <- dat$DURCREM/30
dat$months[dat$months>dat$timetoevent] <- dat$timetoevent[dat$months>dat$timetoevent]
fitLME <- try(lme(log2_val~months*TRT,random=~1+months|participantId,data=dat),silent=T)
fitSURV <- coxph(Surv(timetoevent,event)~TRT,data=dat[!duplicated(dat$participantId),],x=T)
dForm <- list(fixed=~1+TRT,indFixed=c(2,4),random=~1,indRandom=2)
fit.JM <- try(jointModel(fitLME,fitSURV,timeVar="months",method="spline-PH-aGH",parameterization="both",derivForm=dForm),silent=T)

## build a specified trajectory
trj <- NULL
v0 <- median(dat$log2_val[dat$TRT=="CYC"&!is.na(dat$log2_val)&dat$months==0])
trj <- rbind(trj,data.frame(id="CYC_nochg",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nochg",log2_val=v0,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nochg",log2_val=v0,months=6,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_upup",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_upup",log2_val=v0+.5,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_upup",log2_val=v0+1,months=6,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_noup",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_noup",log2_val=v0,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_noup",log2_val=v0+1,months=6,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nodn",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nodn",log2_val=v0,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_nodn",log2_val=v0-1,months=6,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_dndn",log2_val=v0,months=0,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_dndn",log2_val=v0-.5,months=3,TRT="CYC"))
trj <- rbind(trj,data.frame(id="CYC_dndn",log2_val=v0-1,months=6,TRT="CYC"))
v0 <- median(dat$log2_val[dat$TRT=="RTX"&!is.na(dat$log2_val)&dat$months==0])
trj <- rbind(trj,data.frame(id="RTX_nochg",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nochg",log2_val=v0,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nochg",log2_val=v0,months=6,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_upup",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_upup",log2_val=v0+.5,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_upup",log2_val=v0+1,months=6,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_noup",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_noup",log2_val=v0,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_noup",log2_val=v0+1,months=6,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nodn",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nodn",log2_val=v0,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_nodn",log2_val=v0-1,months=6,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_dndn",log2_val=v0,months=0,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_dndn",log2_val=v0-.5,months=3,TRT="RTX"))
trj <- rbind(trj,data.frame(id="RTX_dndn",log2_val=v0-1,months=6,TRT="RTX"))
md <- survfitJM(fit.JM,trj,idVar="id",); xrg=c(6,14)

pdf("../../04 Analysis Output/01 Figures/surv_pred_MemBIgDnIgMp.pdf",width=6,height=6)
## CYC
plot(1,.5,type="n",xlab="Months",ylab="Predicted survival",main="CYC",xlim=xrg,ylim=c(0,1))
tmp <- md$summaries[[1]]; lines(tmp[,1],tmp[,2],col="black",lwd=2); lines(tmp[,1],tmp[,4],col="black",lty="dashed"); lines(tmp[,1],tmp[,5],col="black",lty="dashed")
tmp <- md$summaries[[2]]; lines(tmp[,1],tmp[,2],col="red",lwd=2); lines(tmp[,1],tmp[,4],col="orange",lty="dashed"); lines(tmp[,1],tmp[,5],col="orange",lty="dashed")
tmp <- md$summaries[[3]]; lines(tmp[,1],tmp[,2],col="orange",lwd=2); lines(tmp[,1],tmp[,4],col="red",lty="dashed"); lines(tmp[,1],tmp[,5],col="red",lty="dashed")
tmp <- md$summaries[[4]]; lines(tmp[,1],tmp[,2],col="green",lwd=2); lines(tmp[,1],tmp[,4],col="green",lty="dashed"); lines(tmp[,1],tmp[,5],col="green",lty="dashed")
tmp <- md$summaries[[5]]; lines(tmp[,1],tmp[,2],col="blue",lwd=2); lines(tmp[,1],tmp[,4],col="blue",lty="dashed"); lines(tmp[,1],tmp[,5],col="blue",lty="dashed")
legend("bottomleft",legend=c("no change","constantly up","stable and up","stable and down","constantly down"),col=c("black","red","orange","green","blue"),lwd=2,lty="solid",bty="n")
## RTX
plot(1,.5,type="n",xlab="Months",ylab="Predicted survival",main="RTX",xlim=xrg,ylim=c(0,1))
tmp <- md$summaries[[6]]; lines(tmp[,1],tmp[,2],col="black",lwd=2); lines(tmp[,1],tmp[,4],col="black",lty="dashed"); lines(tmp[,1],tmp[,5],col="black",lty="dashed")
tmp <- md$summaries[[7]]; lines(tmp[,1],tmp[,2],col="red",lwd=2); lines(tmp[,1],tmp[,4],col="orange",lty="dashed"); lines(tmp[,1],tmp[,5],col="orange",lty="dashed")
tmp <- md$summaries[[8]]; lines(tmp[,1],tmp[,2],col="orange",lwd=2); lines(tmp[,1],tmp[,4],col="red",lty="dashed"); lines(tmp[,1],tmp[,5],col="red",lty="dashed")
tmp <- md$summaries[[9]]; lines(tmp[,1],tmp[,2],col="green",lwd=2); lines(tmp[,1],tmp[,4],col="green",lty="dashed"); lines(tmp[,1],tmp[,5],col="green",lty="dashed")
tmp <- md$summaries[[10]]; lines(tmp[,1],tmp[,2],col="blue",lwd=2); lines(tmp[,1],tmp[,4],col="blue",lty="dashed"); lines(tmp[,1],tmp[,5],col="blue",lty="dashed")
legend("bottomleft",legend=c("no change","constantly up","stable and up","stable and down","constantly down"),col=c("black","red","orange","green","blue"),lwd=2,lty="solid",bty="n")
dev.off()
