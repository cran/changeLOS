"progdismodel" <- function(data, plot=FALSE, file1, file2, file3, max.time, 
                           lwd=2, cex=1.2,
                           lty1=1, lty2=2, lty3=3,
                           color1=4, color2=1, color3=2) {
## --------------------------------------------------------------------------------
## Title: R-function progdismodel
## ---------------------------------------------------------------------------------
## Author: Matthias Wangler
##         mw@imbi.uni-freiburg.de
##         Institute of Med. Biometry and Med. Computer Science
##         Stefan-Meier-Strasse 26, D-79104 Freiburg,
##         http://www.imbi.uni-freiburg.de
## ---------------------------------------------------------------------------------
## Description: 
## ---------------------------------------------------------------------------------
## Required Packages: -
## ---------------------------------------------------------------------------------
## Usage: progdismodel(data, file1, file2, file3, max.time,
##                     lwd=2, cex=1.2, lty1=1, lty2=2, lty3=3,
##                     color1=4, color2=1, color3=2)
##
## data: data.frame of the form data.frame( id, j.01, j.02, j.03, j.12, j.13, cens):
##
## id:      id (patient id, admision id, ...)
## j.01:    observed time for jump from "0" to "1"
## j.02:    observed time for jump from "0" to "2"
## j.03:    observed time for jump from "0" to "3"
## j.12:    observed time for jump from "1" to "2"
## j.13:    observed time for jump from "1" to "3"
## cens:    observed time for censoring
## ---------------------------------------------------------------------------------
## Value: 
## ---------------------------------------------------------------------------------
## License: GPL 2
##----------------------------------------------------------------------------------
## History:    19.09.2005, Matthias Wangler
##                         first version
## ---------------------------------------------------------------------------------
    

  ## prepare the data
  ## --> one row for each transition
  ob <- prepare.los.data(x = data)
  
  ## describe the model

  ## 6x6 matrix to describe the possible transitions
  tra <- rbind(c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE),
               c(FALSE,TRUE,FALSE,FALSE,TRUE,TRUE),
               c(FALSE,FALSE,TRUE,FALSE,FALSE,FALSE),
               c(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE),
               c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE),
               c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE))
  
  my.model <- msmodel(c("0", "1", "2", "3", "4", "5"), tra, cens.name = "cens")

  
  ## prepare the data for the 'progressive disability model':

  ob$to <- factor(ob$to, levels=my.model$state.names)

  ## set 'to' = 4 if 'from' = 1 and 'to' = 2
  ob[,3][ob[,2]==1 & ob[,3]==2] <- 4

  ## set 'to' = 5 if 'from' = 1 and 'to' = 3
  ob[,3][ob[,2]==1 & ob[,3]==3] <- 5

  ## compute the Aalen-Johansen estimator for the matrices of transition posibilities
  my.trans <- trans(my.model, ob)
  my.aj <- aj(my.trans, s=0, t=max(my.trans$times))

  p00 <- my.aj$matrices[1,1,]
  p01 <- my.aj$matrices[1,2,]
  p02 <- my.aj$matrices[1,3,]
  p03 <- my.aj$matrices[1,4,]
  p04 <- my.aj$matrices[1,5,]
  p05 <- my.aj$matrices[1,6,]
  
  my.times.par <- my.aj$times[ (p00+p02+p03) > 0 & (p01+p04+p05) > 0 & (p03+p05) > 0 ]


  if( !missing(max.time) ) {
    my.times.par <-  my.times.par[my.times.par <= max.time]
  }
  

  p00 <- p00[my.aj$times %in% my.times.par]
  p01 <- p01[my.aj$times %in% my.times.par]
  p02 <- p02[my.aj$times %in% my.times.par]
  p03 <- p03[my.aj$times %in% my.times.par]
  p04 <- p04[my.aj$times %in% my.times.par]
  p05 <- p05[my.aj$times %in% my.times.par]

  prob.death <- (p03 + p05)   ## death

  prob.death.rfa <- p03 / (p00+p02+p03) ## death, given risk factor absent

  prob.death.rfp <- p05 / (p01+p04+p05) ## death, given risk factor present


  my.ar <- prob.death.rfp - prob.death.rfa ## attributable risk
  
  my.par <- (prob.death - prob.death.rfa) / prob.death ## population attributable risk

  names(my.par) <- NULL
  names(my.ar) <- NULL
  names(prob.death) <- NULL
  names( prob.death.rfa) <- NULL
  names( prob.death.rfp) <- NULL

  if( plot == TRUE ) {
    if( missing(file1) ) {
      x11()
    }
    else {
      postscript(file=file1)
    }          
    
    plot(y=prob.death, x=my.times.par, type="l",
         xlab="Time t", ylab="mortality",
         main="Mortality",
         lwd=lwd,
         lty=lty1,
         col=color1,
         ylim=c(min(prob.death, prob.death.rfa, prob.death.rfp), max(prob.death, prob.death.rfa, prob.death.rfp) ) )
    
    lines(y=prob.death.rfa, x=my.times.par, type="l", lwd=lwd, lty=lty2, col=color2)
    lines(y=prob.death.rfp, x=my.times.par, type="l", lwd=lwd, lty=lty3, col=color3)
    
    legend(x=max(my.times.par)/2,
           y=max(prob.death, prob.death.rfa, prob.death.rfp)/4,
           legend=c("P( death )","P( death | risk factor absent)","P( death | risk factor present)") ,
           bty="o", pch=c(15,15,15), cex=cex,
           lty=c(lty1,lty2,lty3), col=c(color1,color2,color3))
    
    
    if( !missing(file1) ) {
      dev.off()
    }
    
    
    if( missing(file2) ) {
      x11()
    }
    else {
      postscript(file=file2)
    }
    
    plot(y=my.ar, x=my.times.par, type="l",
         xlab="Time t", ylab="AR",
         main="Attributable Mortality",
         lwd=lwd,
         col=color1)
    
    lines(y=rep(0,length(my.times.par)), x=my.times.par, type="l",lty=2)
    
    if( !missing(file2) ) {
      dev.off()
    }

    if( missing(file3) ) {
      x11()
    }
    else {
      postscript(file=file3)
    }
    
    plot(y=my.par, x=my.times.par, type="l",
         xlab="Time t", ylab="PAR",
         main="Population Attributable Mortality",
         lwd=lwd,
         col=color1)
    
    lines(y=rep(0,length(my.times.par)), x=my.times.par, type="l",lty=2)
    
    if( !missing(file3) ) {
      dev.off()
    }
  }

  nt <- length(my.times.par)

  a <- prob.death[nt]
  b <- prob.death.rfa[nt]
  c <- prob.death.rfp[nt]
  d <- my.ar[nt]
  e <- my.par[nt]

  s <- matrix( c(a,b,c,d,e), ncol=1)

  rownames(s) <- c("P(death)", "P(death | risk factor absent)",
                   "P(death | risk factor present)", "Attriburable Mortality",
                   "Population Attrributable Mortality")
    
  colnames(s) <- ""
    
  ## return the result
  res <- list(trans=my.trans,
              aj=my.aj,
              times.par=my.times.par,
              PAR=my.par,
              AR=my.ar,
              death=prob.death,
              death.given.rfa=prob.death.rfa,
              death.given.rfp=prob.death.rfp,
              tab=s)


  return(res)
}
