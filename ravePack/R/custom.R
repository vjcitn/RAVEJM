
custom_km_plot <- function(md,off) {
    par(mar=c(11,16,4,1))
    plot(md$time,md$surv,type="n",xlim=c(0,449),ylim=c(0,1),xlab=NA,ylab=NA,font.axis=2,font.lab=2,cex.axis=1.3)
    nstr=length(md$strata); stt=0; #lntp=c("solid","dash","dotdash","twodash")
    grp=names(md$strata); grp=sapply(strsplit(grp,"="),function(x) x[2])
    for(i in 1:nstr) {
        nln=md$strata[i]; shtg=grp[i]
        xd=md$time[(stt+1):(nln+stt)]; cns=md$n.censor[(stt+1):(nln+stt)]
        yd=md$surv[(stt+1):(nln+stt)]; nrk=md$n.risk[(stt+1):(nln+stt)]
        stt=stt+nln
        lines(xd[rep(1:2,nln-1)+rep(0:(nln-2),each=2)],yd[rep(1:(nln-1),each=2)],lty=i,lwd=2)
        points(xd[cns>0],yd[cns>0],pch=20+i,cex=1.5,lwd=2)
        lnr0=md$n[i]
        lnr1=ifelse(sum(xd<=100)==0,md$n[i],min(nrk[xd<=100]))
        lnr2=ifelse(sum(xd<=200)==0,md$n[i],min(nrk[xd<=200]))
        lnr3=ifelse(sum(xd<=300)==0,md$n[i],min(nrk[xd<=300]))
        lnr4=ifelse(sum(xd<=400)==0,md$n[i],min(nrk[xd<=400]))
        shtg=gsub("Relapsing Disease","Relapsing",shtg)
        shtg=gsub("New Diagnosis","New",shtg)
        shtg=gsub("No Renal Disease at Baseline","No Renal",shtg)
        shtg=gsub("Renal Disease at Baseline","Renal",shtg)
        shtg=gsub("Granulomatosis with Polyangiitis \\(GPA\\)","GPA",shtg)
        shtg=gsub("Microscopic Polyangiitis \\(MPA\\)","MPA",shtg)
        mtext(shtg,at=-290+off*10,adj=0,line=1+1.2*(i+off),font=2,side=1)
        mtext(c(lnr0,lnr1,lnr2,lnr3,lnr4),at=c(0,100,200,300,400),side=1,line=1+1.2*(i+off),font=2)
    }
mtext("Number at risk",side=1,line=2.2+1.2*(i+off),font=2, cex=1.2)
    legend("bottomleft",lty=1:nstr,lwd=2,pt.cex=1.5,cex=1.7,pch=21:(20+nstr),legend=paste(grp," (n=",md$n,")",sep=""),bg="white",bty="n")
    mtext("Time from Complete Remission to Relapse (days)",side=1,line=2.7,font=2, cex=1.2)
    mtext("Probability of Remaining in Complete Remission",side=2,line=3,font=2, cex=1.2)
}
