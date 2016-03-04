
inputf1=`echo $1`
# inputf2=`echo $2`
cat $inputf1 | grep -v nan | ruby -e 'ARGF.each{|x| print x if x.split("\s").size==5 }' | ruby -e 'ARGF.each{|x| print x if x.split("\s")[1..3].join("") !~ /S/ }'> t
# cat $inputf2 | grep -v nan | ruby -e 'ARGF.each{|x| print x if x.split("\s").size==5 }' | ruby -e 'ARGF.each{|x| print x if x.split("\s")[1..3].join("") !~ /S/ }'> s

R --quiet --vanilla > log.txt 2>&1 <<EOF 

  t <-read.table("t", sep=" ")
#   s <-read.table("s", sep=" ")
  xsize <- nrow(t)
#   if( nrow(s) > nrow(t) ) xsize <- nrow(s)
  ysize <- ncol(t)
#   if( ncol(s) > ncol(t) ) ysize <- ncol(s)

  png( "omegaplot.png" )
  plot( (1:nrow(t)), sort( t$V3 ), col="blue", pch=1, xlim=c(0, xsize), ylim=c(0, ysize), xlab="", ylab="")
  par(new=T)
  plot( (1:nrow(t)), sort( t$V4 ), col="red", pch=1, xlim=c(0, xsize), ylim=c(0, ysize), xlab="", ylab="")
  par(new=T)
  plot( (1:nrow(t)), sort( t$V5 ), col="black", pch=1, xlim=c(0, xsize), ylim=c(0, ysize), xlab="", ylab="" )
#   par(new=T)
#   plot( (1:nrow(s)), sort( s$V3 ), col="blue", pch=4, xlim=c(0, xsize), ylim=c(0, ysize), xlab="", ylab="")
#   par(new=T)
#   plot( (1:nrow(s)), sort( s$V4 ), col="red", pch=4, xlim=c(0, xsize), ylim=c(0, ysize), xlab="", ylab="")
#   par(new=T)
#   plot( (1:nrow(s)), sort( s$V5 ), col="black", pch=4, xlim=c(0, xsize), ylim=c(0, ysize), xlab="Sorted Gene Pairs", ylab="omega (dN/dS)")
  abline(h=1.0)
  dev.off()

EOF
