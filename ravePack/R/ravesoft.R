setClass("raveSet", contains="ExpressionSet")
#> m[1,]
#    participantId sequenceNum vstnum     combdt days_from_rand
#1 ITN021AI_001002       20008     V8 2005-06-28            179
#  days_from_comrem1 TRT                       AAVTYPE ANCASTAT NEWDXYN RENFLBL
#1                 0 RTX Wegener's Granulomatosis (WG)      PR3      No      No
#  TSCORE PRED0 GLCDOSE CREMFL CREMFLDT relapse_class TMCREM DURCREM DURCREMF
#1      0   Yes       0                                  181     358        1
#         clin_status value Assay_group     Assay_test log2_val
#1 complete remission  2.16        Flow B cells (CD19) 1.111031
#
#setClass("raveSet", contains="ExpressionSet")

#
#new_defn,USUBJID,TRT,value,VISDT,BL_date,days_from_BL,elisa,TSCORE,GLCDOSE,GLCCUM,xv,anal_svfl,anal_lmfl,flare_mod,censor_day,log2_value,flare_status,anal_flare
#B cells (CD19),1002,RTX,2.16,28-Jun-05,30-Dec-04,179,29.3,0,0,70,V8,,,Normal,536,1.111031312,Normal,
#


###################################################
### code chunk number 4: domff
###################################################
 makeFlowFrame = function(subj, NMARKERS=3) {
  dates = names(subj)
  rs = sapply(subj, nrow)
  allok = which(rs == NMARKERS)
  ok = min(allok)
  fulln = as.character(subj[[ok]]$Assay_test)      
  ncol = length(dates)
  nrow = NMARKERS
  targ = matrix(NA, nr=nrow, nc=ncol)
  subj = subj[allok]
  rownames(targ) = fulln
  for (i in 1:length(subj)) {
     curdat = subj[[ i ]]
     dat = curdat$log2_val
     names(dat) = as.character(curdat$Assay_test)
     stopifnot(all(names(dat) == fulln))
     targ[ names(dat), i ] = dat
     }
  cn = paste("s", as.character(subj[[1]]$participantId[1]), "d", dates, sep="_")
  colnames(targ) = cn
  targ
 }


grabSspec = function(rs, vbl, idvar="id") {
  ids = pData(rs)[[idvar]]
  dat = pData(rs)[[vbl]]
  sdat = sapply(split(dat, ids), "[", 1)
  sdat
}


getz = function(ind, raveset, depvar="B cells (CD19)", timevar="days", predtime = seq(0,200,25),
    idvar = "id", trtvar = "tx", fudge=.01) {
  rl200 = ravef[ , which(ravef$days < 200)]
 dep = exprs(raveset)[ depvar, ]
 tim = pData(raveset)[[timevar]]
 id = pData(raveset)[[idvar]]
 tx = pData(raveset)[[trtvar]]
 basef = na.omit(data.frame(dep, time=tim, id=id, tx=tx))
 names(basef) = c(depvar, timevar, idvar, trtvar)
 sod = split(basef, basef$id)[[ind]]
 thecall = match.call()
 Y = sod[[depvar]]
 T = sod[[timevar]]
 if (length(Y)>5) {
   mod = lm(Y~poly(T,2),x=TRUE,y=TRUE)
   modc = summary(mod)$coef
   z3 = modc[,1]/(modc[,2]+fudge)
#   } else if (length(Y)>2) {
##   mod = lm(Y~T,x=TRUE,y=TRUE)
#   modc = summary(mod)$coef
#   z3 = c(modc[,1]/(modc[,2]+fudge),0)
   } else {
   mod = NULL
   z3 = c(mean(Y)/(sd(Y)+fudge),0,0)
   }
 if (!is.null(mod)) preds = predict(mod, newdata=list(T=predtime))
 else preds = rep(z3[1], length(predtime))
 ans = list(z3=z3, preds=preds, indat=sod, mod=mod, depvar=depvar, idvar=idvar, timevar=timevar, call=thecall, predtime=predtime)
 class(ans) = "qmod"
 ans
}
print.qmod = function(x, ...) {
  cat("qmod instance with z-scores\n")
  print(x$z3)
}
plot.qmod = function(x, ...) {
  Y = x$indat[[x$depvar]]
  T = x$indat[[x$timevar]]
  rng = range(c(Y,x$preds))
  plot(Y~T, xlab=x$timevar, ylab=x$depvar, ylim=rng, ...)
  lines(x$preds~x$predtime)
  invisible(NULL)
}
getcox = function(vbl="B cells (CD19)") {
 require(multicore)
 allz = mclapply(1:length(unique(ravef$id)), function(i) getz(i, 
   ravef[, ravef$days < 200], depvar=vbl))
 pc = try(prcomp(nn <- na.omit(t(sapply(allz, function(x)x$z3)))))
 ann = attributes(nn)
 kp = 1:length(timeon)
 if (inherits(pc, "try-error")) coxmod=NA
 else coxmod=coxph(Surv(timeon[kp], event[kp]=="Disease Flare")~pc$x[,1]+pc$x[,2]+pc$x[,3])
 ans = list(vbl=vbl, coxmod=coxmod, allz=allz)
 class(ans) = "ravecox"
 ans
}
print.ravecox = function(x, ...) {
 cat(paste("ravecox instance for vbl", x$vbl, " with coxmod summary:\n"))
 print(summary(x$coxmod))
}
