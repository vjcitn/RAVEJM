
 setClass("JMset", contains="list", representation(ids="character"))
 setMethod("show", "JMset", function(object) cat("JMset with", length(object), "elements\n"))

JMset = function(list) {
 new("JMset", list, ids=names(list))
}

setMethod("[", "JMset", function(x, i, j, ..., drop=FALSE) {
 new("JMset", as(x, "list")[i], ids=x@ids[i])
})

setGeneric("filter", function(x, f) standardGeneric("filter"))
setMethod("filter", c("JMset", "function"), function(x, f) {
 xf = lapply(x, f)
 ok = sapply(xf, function(z) !is.null(z))
 nids = x@ids[ok]
 new("JMset", xf[ok], ids=nids)
})

setAs("JMset", "data.frame", function (from, to, strict=TRUE) {
  as.data.frame(do.call(rbind, from))
})

addStartEnd = function(x, timevbl="days_from_comrem1",
  outcomevbl="clin_status", eventcode="relapse") {
  stopifnot(is(x, "data.frame"))
  nrec = nrow(x)
  nf = x[-nrec,]
  nf$start = x[[timevbl]][-nrec]
  nf$end = x[[timevbl]][-1]
  nf$event = 1*(x[[outcomevbl]][-1] == eventcode)
  nf[[outcomevbl]] = x[[outcomevbl]][-1]
  nf
}


#MSET = JMset(msplit)

.plotOne = function(x, id)  {
  with(x, {
   if (FALSE) {
      plot(days_from_rand, log2_val, col = as.numeric(Assay_test), pch = 19, main=id)
      legend(0,4,legend=levels(Assay_test), pch=19, col=1:length(levels(Assay_test)))
      x$remissionStarts = rep(na.omit(with(x, unique(days_from_rand[days_from_comrem1==0])))[1], nrow(x))
      suppressWarnings({relapseStarts = na.omit(with(x, min(days_from_rand[clin_status == "relapse"])))})
      if (length(remissionStarts)>1) warning(paste("multiple remission start dates", id))
      abline(v=remissionStarts, col="gray", lwd=3)
      abline(v=relapseStarts, col="orange", lwd=3)
      timerange = range(days_from_rand)
      remissionTimeScale = timerange - remissionStarts[1]
#      axis(3, at=c(-25,0,25)+remissionStarts, labels=c(-25,0,25))
#      mtext("time from compl. remiss.", 3)
      invisible(NULL)
      }
      x$remissionStarts = rep(na.omit(with(x, unique(days_from_rand[days_from_comrem1==0])))[1], nrow(x))
      suppressWarnings({x$relapseStarts = na.omit(with(x, min(days_from_rand[clin_status == "relapse"])))})
      p = ggplot(x, aes(x=days_from_rand, y=log2_val, colour=Assay_test)) + geom_point() +
          geom_vline(xintercept=x$remissionStarts, colour="green")
      if (is.finite(x$relapseStarts))
           return(p + geom_vline(xintercept=relapseStarts, colour="orange"))
      return(p)
      })
}

setGeneric("plotOne", function(set, ind) standardGeneric("plotOne")) #c("JMset", "numeric"), 
setMethod("plotOne", c("JMset", "numeric"), function(set, ind) {
  .plotOne(set[[ind]], set@ids[ind] )
})

collectFour = function(x, inds) {
  stopifnot(length(inds) ==4)
  myl = x[inds]
  nl = lapply(1:4, function(ind) {
      curx = myl[[ind]]
      minrem = min(curx$days_from_comrem1, na.rm=TRUE)
      curx$remissionStarts = rep(na.omit(with(curx, unique(days_from_rand[days_from_comrem1==minrem])))[1], nrow(curx))
      suppressWarnings({curx$relapseStarts = rep(na.omit(with(curx, min(days_from_rand[clin_status == "relapse"]))), nrow(curx))})
      curx$ids = rep(x@ids[inds[ind]], nrow(curx))
      curx$row = ifelse(ind <= 2, 1, 2)
      curx$col = ifelse(ind %% 2, 1, 2)
      curx
  })
  do.call(rbind, nl)
}

.old.plotFour = function(x, inds) {
      cl = collectFour(x, inds)
      remis = data.frame(t(sapply(split(cl, as.character(cl$participantId)), function(x) c(remis=x$remissionStarts[1], row=x$row[1], col=x$col[1]))))
      p = ggplot(cl, aes(x=days_from_rand, y=log2_val, colour=Assay_test)) + geom_point() +
          geom_vline(xintercept=remis, colour="green", data=remis) + facet_grid(row~col)
      return(p)
}

#remrefdata = sapply(split(MM$remissionStarts, MM$participantId), "[", 1)
#relrefdata = sapply(split(MM$relapseStarts, MM$participantId), "[", 1)
#xyplot(log2_val~days_from_rand|participantId, data=MM, col=as.numeric(MM$Assay_test), pch=19, 
#  groups=MM$Assay_Test, auto.key=TRUE, panel=
# function(...) {panel.xyplot(...); panel.abline(v=remrefdata[packet.number()],col="green"); 
#               panel.abline(v=relrefdata[packet.number()], col="red")})

xyplFour = function(df) {
 remrefdata = sapply(split(df$remissionStarts, as.character(df$participantId)), "[", 1)
 relrefdata = sapply(split(df$relapseStarts, as.character(df$participantId)), "[", 1)
 xyplot(log2_val~days_from_rand|participantId, data=df, col=as.numeric(df$Assay_test), pch=19, 
  groups=df$Assay_Test, auto.key=TRUE, panel=
 function(...) {panel.xyplot(...); panel.abline(v=remrefdata[packet.number()],col="green"); 
               panel.abline(v=relrefdata[packet.number()], col="red")})
}

.plotFour = function(x, inds) {
 xyplFour(x[inds,])
}

setGeneric("plotFour", function(x, inds) standardGeneric("plotFour"))
setMethod("plotFour", c("JMset", "numeric"), function(x, inds) {
  d = collectFour(x, inds)
  xyplFour(d)
})
