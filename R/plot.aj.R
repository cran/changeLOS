"plot.aj" <-
function(x,from,to,xlab=expression(paste(Time, " ", italic(t))),
                    ylab= eval(substitute(expression(paste("Estimate of ", P[{a}][{b}], "(", italic(s), ",", italic(t), ")")),
                               list(a=from[1],b=to[1],s=x$start))),
                    xlim = c(x$start,max(x$times)), ylim=c(0,1),
                    lab=c(10,10,7),
                    txt=eval(substitute(expression(paste(hat(P)[{a}][{b}], "(",italic( s), ",", italic(t), ")")),
                             list(a=from[1],b=to[1],s=x$start))),
                    x.txt=(xlim[2]+xlim[1])/2,y.txt=ylim[2]*0.9,
                    col=1,...){
## ----------------------------------------------------------------------------
## Title: plot.aj
## ----------------------------------------------------------------------------
## Author: Matthias Wangler, mw@imbi.uni-freiburg.de
## Institute of Med. Biometry and Med. Computer Science
## Stefan-Meier-Strasse 26, D-79104 Freiburg,
## http://www.imbi.uni-freiburg.de
## ----------------------------------------------------------------------------
## Description: plot function for an object of class 'aj'
## ----------------------------------------------------------------------------
## Required Packages: 
## ----------------------------------------------------------------------------
## Usage: plot.aj(x, from, to, xlab, ylab, xlim, ylim, lab, txt, x.txt, y.tyt, col)
##
## x:           an object of the class 'aj'
## from:        a character vector naming the states 'from'
## to:          a character vector naming the states 'to'
## xlab:        a title for the x axis
## ylab:        a title for the y axis
## xlim:        the x limits (min,max) of the plot
## ylim:        the y limits (min,max) of the plot
## lab:         A numerical vector of the form 'c(x, y, len)' which modifies
##              the way that axes are annotated.  The values of 'x' and 'y'
##              give the (approximate) number of tickmarks on the x and y
##              axes and 'len' specifies the label size.
## txt:         one or more character strings or expressions specifying a
##              text to be written.
## x.txt,y.txt: the x and y co-ordinates to be used to position the text
## col:         the color of the lines
## ----------------------------------------------------------------------------
## Value: a matrix with ncol = 1 + length(from)
##
## column 1:                vector of the time points (the x coordinates)
## column 2 to column ncol: vector of the y coordinates
## ----------------------------------------------------------------------------
## Notes: -
## ----------------------------------------------------------------------------
## Example: > data(los.data)
##          > my.observ <- prepare.los.data(x=los.data)
##          > my.model <- msmodel(c("0","1","2","3"),cens.name="cens")
##          > my.trans <- trans(model=my.model,observ=my.observ)
##          > my.aj <- aj(my.trans, s=0, t=80)  
##          > plot(my.aj,c("0","0","0","0"),c("0","1","2","3"))  
## ----------------------------------------------------------------------------
## License: GPL 2
##-----------------------------------------------------------------------------
## History: 11.08.2004, Matthias Wangler
##                      the first version
## ----------------------------------------------------------------------------
  if (!inherits(x, "aj")) {
    stop("Argument 'x' must be an object of class \"aj\".")
  }

  if( !is.character(from) ) {
    stop("Argument 'from' must be a character string.")
  }

  if( !is.character(to) ) {
    stop("Argument 'to' must be a character string.")
  }

  if(length(from) != length(to)) {
    stop("The number of states 'from' must be equal to the number of states 'to'.")
  }  
    
  state.names <- dimnames(x$matrix)[[1]]

  if( !(all(from %in% state.names)) ) {
    stop("Argument 'from' must be valid statenames.")
  }

  if( !all((to %in% state.names)) ) {
    stop("Argument 'from' must be valid statenames.")
  }

  states <- 1:length(state.names)
  
  ## save default, for resetting
  def.par <- par(no.readonly = TRUE)
  on.exit(par(def.par))

  column <- 1
  row <- 1
  v <- 1:length(from)
  
  if( length(from) > 1 ) {
    column <- 2

    row <- length(from)%/%2

    if( (length(from)%%2) > 0 ) {
      row <- row + 1

      v <- c(v,length(from)+1)
    }
  }
    
  layout(matrix(v,row,column, byrow=TRUE))
  
  par(lab=lab)

  v <- c(x$times)

  clb <- c("Time")
  
  for( i in 1:length(from)) {    
    if( length(from) > 1 ) {
      ylab <- eval(substitute(expression(paste("Estimate of ", P[{a}][{b}], "(", italic(s), ",", italic(t), ")")),
                              list(a=from[i],b=to[i],s=x$start)))
      txt <- eval(substitute(expression(paste(hat(P)[{a}][{b}], "(",italic( s), ",", italic(t), ")")),
                             list(a=from[i],b=to[i],s=x$start)))
      op <- par(mar=c(4, 5, 2, 1))
      plot(x=x$times,y=x$matrices[states[state.names==from[i]],states[state.names==to[i]],],
           type ="s", lty=1, lwd=2, xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab, col=col, ...)
      par(op)
    }
    else {      
      plot(x=x$times,y=x$matrices[states[state.names==from[i]],states[state.names==to[i]],],
           type ="s", lty=1, lwd=2, xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab, col=col, ...)
    }
    
    v <- c(v,x$matrices[states[state.names==from[i]],states[state.names==to[i]],])

    clb <- c(clb,paste(from[i],to[i]))
           
    text(x.txt, y.txt, txt, ...)
  }

   m <- matrix( v, nrow=length(x$times), ncol=length(from)+1,byrow=FALSE)

   dimnames(m) <- list(NULL,clb)

  return(m)
}

