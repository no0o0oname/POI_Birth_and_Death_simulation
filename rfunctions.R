initialize = function(density.per.area=1,species=1,U=1,dim=2,poisson=F) {
 if (dim == 1) {
   if (poisson) {
    n=rpois(1,U*density.per.area)
   }
   else {
    n=U*density.per.area
   }
   x=matrix(0,nrow=n,ncol=2)
   x[,1] = species
   x[,2] = runif(n,min=0,max=U)
 }
 else if (dim == 2) {
   if (poisson) {
    n=rpois(1,U*U*density.per.area)
   }
   else {
    n=U*U*density.per.area
   }
   x=matrix(0,nrow=n,ncol=3)
   x[,1] = species
   x[,2] = runif(n,min=0,max=U)
   x[,3] = runif(n,min=0,max=U)
 }
 else {
  stop(sprintf("dimension should be 1 or 2 (was %d).\n",dim));
 }
 x
}

write.file = function(x,filename) {
 write.table(x,filename,col.names=F,row.names=F)
}

read.singletime.file = function(x,filename) {
 a=as.numeric(read.table(filename,header=F))
 list(time=0, event.num=0, species=a[,1], coord=a[,2:ncol(a)])
}

read.timeseries.file = function(filename) {
 a=read.table(filename,header=F,sep="\n",colClasses="character")
}

get.timeseries.point.count = function(a,dim,splist=NULL) {
 if (dim==1) {
  K=2
 }
 else if (dim==2) {
  K=3
 }
 else {
  stop(sprintf("dimension should be 1 or 2 (was %d).\n",dim));
 }
 if ((length(splist)) > 0) {
  L=length(splist)
  y=data.frame(time=rep(0,nrow(a)), event.num=rep(0,nrow(a)), count=matrix(0,nr=nrow(a),nc=L))
  for (j in 1:L) {
   colnames(y)[2+j] = sprintf("count.s%d",splist[j])
  }
  for (i in 1:nrow(a)) {
   b=as.numeric(strsplit(a[i,], " ")[[1]])
   y$time[i]=b[1]
   y$event.num[i]=b[2]
   species=b[seq(3,length(b),K)]
   for (j in 1:L) {
    y[i,2+j] = sum(species == splist[j])
   }
  }
 }
 else {
  y=data.frame(time=rep(0,nrow(a)), event.num=rep(0,nrow(a)), count=rep(0,nrow(a)))
  for (i in 1:nrow(a)) {
   b=as.numeric(strsplit(a[i,], " ")[[1]])
   y$time[i]=b[1]
   y$event.num[i]=b[2]
   y$count[i]=(length(b)-2)/K
  }
 }
 y
}

get.point.coordinates = function(a,ind,dim) {
 b=as.numeric(strsplit(a[ind,], " ")[[1]])
 n=length(b)
 if (dim==1) {
  species=b[seq(3,n,2)]
  coord=b[seq(4,n,2)]
 }
 else if (dim==2) {
  species=b[seq(3,n,3)]
  coord=cbind(b[seq(4,n,3)], b[seq(5,n,3)])
 }
 else {
  stop(sprintf("dimension should be 1 or 2 (was %d).\n",dim));
 }
 list(time=b[1], event.num=b[2], species=species, coord=coord)
}

plot.neighborhood = function(radius=1) {
 a=seq(0,2*pi,.1)
 x = radius * cos(a)
 y = radius * sin(a)
 p=locator(1,type="p",pch=3)
 while (!is.null(p)) {
  polygon(p$x+x,p$y+y,density=10)
  p=locator(1,type="p",pch=3)
 }
}

plot.points = function(x,species,U,spch=NULL,scol=NULL,scex=NULL,xlab="",ylab="",xlim=c(0,U),ylim=c(0,U),...) {
 sp=19
 if (!is.null(spch)) {
  if (max(species) > length(spch)) {
   stop(sprintf("spch vector does not cover all species (max species %d, spch vector length %d).\n",max(species),length(spch)))
  }
  sp = spch[species]
 }
 sc=species
 if (!is.null(scol)) {
  if (max(species) > length(scol)) {
   stop(sprintf("scol vector does not cover all species (max species %d, scol vector length %d).\n",max(species),length(scol)))
  }
  sc = scol[species]
 }
 sx=1
 if (!is.null(scex)) {
  if (max(species) > length(scex)) {
   stop(sprintf("scex vector does not cover all species (max species %d, scex vector length %d).\n",max(species),length(scex)))
  }
  sx = scex[species]
 }
 plot(x[,1], x[,2], pch=sp, col=sc, cex=sx, xlab=xlab,ylab=ylab, xlim=xlim, ylim=ylim, ...)
}

movie = function(a,U,pause.interval=0.1,spch=NULL,scol=NULL,scex=NULL,xlab="",ylab="",xlim=c(0,U),ylim=c(0,U), ...) {
 for (i in 1:nrow(a)) {
  b=as.numeric(strsplit(a[i,], " ")[[1]])
  n=length(b)
  x=b[seq(4,n,3)]
  y=b[seq(5,n,3)]
  species=b[seq(3,n,3)]

  sp=19
  if (!is.null(spch)) {
   if (max(species) > length(spch)) {
    stop(sprintf("spch vector does not cover all species (max species %d, spch vector length %d).\n",max(species),length(spch)))
   }
   sp = spch[species]
  }
  sc=species
  if (!is.null(scol)) {
   if (max(species) > length(scol)) {
    stop(sprintf("scol vector does not cover all species (max species %d, scol vector length %d).\n",max(species),length(scol)))
   }
   sc = scol[species]
  }
  sx=1
  if (!is.null(scex)) {
   if (max(species) > length(scex)) {
    stop(sprintf("scex vector does not cover all species (max species %d, scex vector length %d).\n",max(species),length(scex)))
   }
   sx = scex[species]
  }
  plot(x,y,pch=sp,col=sc,cex=sx,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab, ...)
  grid()
  title(sprintf("time %.2f %3d events",b[1],b[2]))
  Sys.sleep(pause.interval)
 }
}

movie2file = function(a,U,filename,spch=NULL,scol=NULL,scex=NULL,xlab="",ylab="",xlim=c(0,U),ylim=c(0,U), ...) {
 for (i in 1:nrow(a)) {
  outfile=sprintf("%s%06d.png",filename,i)
  png(outfile)
  b=as.numeric(strsplit(a[i,], " ")[[1]])
  n=length(b)
  x=b[seq(4,n,3)]
  y=b[seq(5,n,3)]
  species=b[seq(3,n,3)]

  sp=19
  if (!is.null(spch)) {
   if (max(species) > length(spch)) {
    stop(sprintf("spch vector does not cover all species (max species %d, spch vector length %d).\n",max(species),length(spch)))
   }
   sp = spch[species]
  }
  sc=species
  if (!is.null(scol)) {
   if (max(species) > length(scol)) {
    stop(sprintf("scol vector does not cover all species (max species %d, scol vector length %d).\n",max(species),length(scol)))
   }
   sc = scol[species]
  }
  sx=1
  if (!is.null(scex)) {
   if (max(species) > length(scex)) {
    stop(sprintf("scex vector does not cover all species (max species %d, scex vector length %d).\n",max(species),length(scex)))
   }
   sx = scex[species]
  }

  plot(x,y,pch=sp,col=sc,cex=sx,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab, ...)
  grid()
  title(sprintf("time %.2f %3d events",b[1],b[2]))
  dev.off()
 }
}

spm2 = function(x,U,grid.resolution=1,xlim=NULL, cumulant=F) {
 if (is.null(dim(x))) {
  dim=1
 }
 else {
  dim=ncol(x)
 }
 if (is.null(xlim)) {
  xlim=c(0,U/2);
 }

 if (dim == 1) {
  # count array
  gn=ceiling(U/grid.resolution)
  gcount=tabulate(floor(x/grid.resolution) + 1, nbins=gn)

  # spatial correlation via 1D FFT
  y=Re(fft(abs(fft(gcount))^2, inverse=T))/length(gcount)
  hd=1:(floor(xlim[2]/grid.resolution))
  hc = y[hd+1]
  z=U*grid.resolution
  spm2 = hc/z
  if (cumulant) {
   de=length(x)/U 
   spm2 = spm2 - de*de;
  }
 }
 else if (dim == 2) {
  # count matrix
  gn=ceiling(U/grid.resolution)
  gcount=matrix(tabulate(floor(x[,1]/grid.resolution)*gn + floor(x[,2]/grid.resolution) + 1, nbins=gn*gn), nr=gn, nc=gn)

  # spatial correlation via 2D FFT
  y=Re(fft(abs(fft(gcount))^2, inverse=T))/length(gcount)
  hd=1:(floor(xlim[2]/grid.resolution))
  hc=rep(0,length(hd))
  cc=rep(0,length(hd))
  for (i in 1:length(hd)) {
   for (j in 1:length(hd)) {
    d = round(sqrt((i-1)^2 + (j-1)^2))
    if ((d >= min(hd)) & (d <= max(hd))) {
     hc[d] = hc[d] + y[i,j]
     cc[d] = cc[d] + 1
    }
   }
  }
  z=U*U*grid.resolution*grid.resolution
  spm2 = hc/cc/z
  if (cumulant) {
   de=nrow(x)/(U*U) 
   spm2 = spm2 - de*de;
  }
 }
 else {
  stop(sprintf("ERROR: dimension should be 1 or 2 (was %d).\n",dim));
 }

 data.frame(dist=hd*grid.resolution,spm2=spm2)
}

spx = function(x,y,U,grid.resolution=1,xlim=NULL, cumulant=F) {
 if (is.null(dim(x))) {
  dim=1
 }
 else {
  dim=ncol(x)
 }
 if (is.null(xlim)) {
  xlim=c(0,U/2);
 }

 if (dim == 1) {
  # count array
  gn=ceiling(U/grid.resolution)
  gcount1=tabulate(floor(x/grid.resolution) + 1, nbins=gn)
  gcount2=tabulate(floor(y/grid.resolution) + 1, nbins=gn)

  # spatial correlation via 1D FFT
  v=convolve(gcount1,gcount2)
  hd=0:(floor(xlim[2]/grid.resolution))
  hc = v[hd+1]
  if (xlim[2] < U/2) {
   i2=1:(length(hd)-1)
  }
  else {
   i2=1:(length(hd)-2)
  }
  hc[i2+1] = hc[i2+1] + v[gn-i2+1]
  z=U*grid.resolution
  spx = hc/z
  if (cumulant) {
   spx = spx - length(x)/U * length(y)/U
  }
 }
 else if (dim == 2) {
  # count matrix
  gn=ceiling(U/grid.resolution)
  gcount1=matrix(tabulate(floor(x[,1]/grid.resolution)*gn + floor(x[,2]/grid.resolution) + 1, nbins=gn*gn), nr=gn, nc=gn)
  gcount2=matrix(tabulate(floor(y[,1]/grid.resolution)*gn + floor(y[,2]/grid.resolution) + 1, nbins=gn*gn), nr=gn, nc=gn)
   
  # spatial correlation via 2D FFT
  v=convolve(gcount1,gcount2)
  hd=0:(floor(xlim[2]/grid.resolution))
  hc=rep(0,length(hd))
  cc=rep(0,length(hd))
  for (i in 1:length(hd)) {
   for (j in 1:length(hd)) {
    d = round(sqrt((i-1)^2 + (j-1)^2))
    if ((d >= min(hd)) & (d <= max(hd))) {
     hc[d+1] = hc[d+1] + v[i,j]
     cc[d+1] = cc[d+1] + 1
    }
   }
  }
  if (xlim[2] < U/2) {
   n2=length(hd)-1
  }
  else {
   n2=length(hd)-2
  }
  for (i in seq(1,n2)) {
   for (j in seq(gn-n2+1,gn)) {
    d = round(sqrt((i-1)^2 + (j-1 - gn)^2))
    if ((d >= min(hd)) & (d <= max(hd))) {
     hc[d+1] = hc[d+1] + v[i,j]
     cc[d+1] = cc[d+1] + 1
    }
   }
  }

  for (i in seq(gn-n2+1,gn)) {
   for (j in seq(1,n2)) {
    d = round(sqrt((i-1 - gn)^2 + (j-1)^2))
    if ((d >= min(hd)) & (d <= max(hd))) {
     hc[d+1] = hc[d+1] + v[i,j]
     cc[d+1] = cc[d+1] + 1
    }
   }
  }

  for (i in seq(gn-n2+1,gn)) {
   for (j in seq(gn-n2+1,gn)) {
    d = round(sqrt((i-1 - gn)^2 + (j-1 - gn)^2))
    if ((d >= min(hd)) & (d <= max(hd))) {
     hc[d+1] = hc[d+1] + v[i,j]
     cc[d+1] = cc[d+1] + 1
    }
   }
  }

  z=U*U*grid.resolution*grid.resolution
  spx = hc/cc/z
  if (cumulant) {
   de1 = nrow(x)/(U*U)
   de2 = nrow(y)/(U*U)
   spx = spx - de1*de2
  }
 }
 else {
  stop(sprintf("ERROR: dimension should be 1 or 2 (was %d).\n",dim));
 }

 data.frame(dist=hd*grid.resolution,spx=spx)
}

plot.event.count = function(x,events=NULL,pch=NULL,xlab="Time",ylab="Number of events",...) {
 # event count in outputfile is cumulative, so calculate here count within time interval
 if (is.null(events)) {
  i=3:ncol(x)
 }
 else {
  i=2+events
 }
 b=x[2:nrow(x),i] - x[1:(nrow(x)-1),i]
 b0=x[2:nrow(x),2] - x[1:(nrow(x)-1),2]
 if (is.null(pch)) {
  pch=sprintf("%d",i-2)
 }
 matplot(x$time[2:nrow(x)],b,pch=pch,xlab=xlab,ylab=ylab, ...)
}

plot.species.density = function(x,U,dim,species=NULL,pch=NULL,xlab="Time",ylab="Density", ...) {
 area=U
 for (i in 2:dim) {
  area = area*U
 }
 if (is.null(species)) {
  i=4:ncol(x)
 }
 else {
  i=species+3
 }
 if (is.null(pch)) {
  pch=sprintf("%d",i-3)
 }
 matplot(x$time,x[,i]/area,pch=pch,xlab=xlab,ylab=ylab, ...)
}

plot.spm2 = function(coord,species,U,grid.resolution=1,selected.species=NULL,pch=NULL,col=NULL,cex=NULL,xlim=NULL,ylim=NULL,xlab="Distance",ylab=NULL, cumulant=F, ...) {
 if (is.null(selected.species)) {
  selected.species=min(species):max(species)
 }
 s=vector("list",length(selected.species))
 if (is.null(xlim)) {
  xlim=c(0,U/2);
 }

 if (is.null(dim(coord))) {
  dimension=1
 }
 else {
  dimension=ncol(coord)
 }
 
 minval = 0
 maxval = 0
 for (si in 1:length(selected.species)) {  
  i=which(species == selected.species[si])
  if (dimension == 1) {
   s[[si]]=spm2(coord[i], U=U, grid.resolution=grid.resolution,xlim=xlim, cumulant=cumulant)
  }
  else {
   s[[si]]=spm2(coord[i,,drop=F], U=U, grid.resolution=grid.resolution,xlim=xlim, cumulant=cumulant)
  }
  mm = max(s[[si]]$spm2)
  if (mm > maxval) {
   maxval=mm
  }
  mm = min(s[[si]]$spm2)
  if (mm < minval) {
   minval=mm
  }
 }

 if (is.null(pch)) {
  pch = sprintf("%d",selected.species)
 }
 else {
   if (length(pch) == 1) {
    pch=rep(pch,length(selected.species))
   }
   else if (length(selected.species) != length(pch)) {
    stop(sprintf("pch vector length (%d) does not match with selected.species (%d).\n",length(pch), length(selected.species)));
   }
 }

 if (is.null(col)) {
  col = selected.species
 }
 else {
   if (length(col) == 1) {
    col = rep(col,length(selected.species))
   }
   else if (length(selected.species) != length(col)) {
    stop(sprintf("col vector length (%d) does not match with selected.species (%d).\n",length(col), length(selected.species)));
   }
 }
 if (is.null(cex)) {
  cex = rep(1,length(selected.species))
 }
 else {
   if (length(cex) == 1) {
    cex = rep(cex,length(selected.species))
   }
   else if (length(selected.species) != length(cex)) {
    stop(sprintf("cex vector length (%d) does not match with selected.species (%d).\n",length(cex), length(selected.species)));
   }
 }

 if (is.null(ylim)) {
  ylim = c(minval,maxval)
 }
 if (is.null(ylab)) {
  if (cumulant) {
   ylab = "2nd spatial cumulant"
  }
  else {
   ylab = "2nd spatial moment"
  }
 }
 
 plot(s[[1]]$dist,s[[1]]$spm2,pch=pch[1],col=col[1],cex=cex[1], xlim=xlim, ylim=ylim,xlab=xlab,ylab=ylab,...)
 if (length(selected.species) > 1) {
  for (si in 2:length(selected.species)) {
   points(s[[si]]$dist,s[[si]]$spm2,pch=pch[si],col=col[si],cex=cex[si])
  }
 }
}

plot.spx = function(coord1,species1,coord2,species2,U,grid.resolution=1,selected.species=NULL,pch=NULL,col=NULL,cex=NULL,xlim=NULL,ylim=NULL,xlab="Distance",ylab="Spatial crosscorrelation", cumulant=F, ...) {
 if (is.null(selected.species)) {
  selected.species=min(c(species1,species2)):max(c(species1,species2))
 }
 s=vector("list",length(selected.species))
 if (is.null(xlim)) {
  xlim=c(0,U/2);
 }

 if (is.null(dim(coord1))) {
  dimension=1
 }
 else {
  dimension=ncol(coord1)
 }

 minval = 0
 maxval = 0
 for (si in 1:length(selected.species)) {  
  i1=which(species1 == selected.species[si])
  i2=which(species2 == selected.species[si])
  s[[si]] = data.frame(dist=0,spx=0)
  if (length(i1)>0 && length(i2)>0) {
   if (dimension == 1) {
    s[[si]]=spx(coord1[i1], coord2[i2], U=U, grid.resolution=grid.resolution,xlim=xlim,cumulant=cumulant)
   }
   else {
    s[[si]]=spx(coord1[i1,,drop=F],coord2[i2,,drop=F], U=U, grid.resolution=grid.resolution,xlim=xlim,cumulant=cumulant)
   }
   mm = max(s[[si]]$spx)
   if (mm > maxval) {
    maxval=mm
   }
   mm = min(s[[si]]$spx)
   if (mm < minval) {
    minval=mm
   }
  }
 }

 if (is.null(pch)) {
  pch = sprintf("%d",selected.species)
 }
 else {
   if (length(pch) == 1) {
    pch=rep(pch,length(selected.species))
   }
   else if (length(selected.species) != length(pch)) {
    stop(sprintf("pch vector length (%d) does not match with selected.species (%d).\n",length(pch), length(selected.species)));
   }
 }

 if (is.null(col)) {
  col = selected.species
 }
 else {
   if (length(col) == 1) {
    col = rep(col,length(selected.species))
   }
   else if (length(selected.species) != length(col)) {
    stop(sprintf("col vector length (%d) does not match with selected.species (%d).\n",length(col), length(selected.species)));
   }
 }
 if (is.null(cex)) {
  cex = rep(1,length(selected.species))
 }
 else {
   if (length(cex) == 1) {
    cex = rep(cex,length(selected.species))
   }
   else if (length(selected.species) != length(cex)) {
    stop(sprintf("cex vector length (%d) does not match with selected.species (%d).\n",length(cex), length(selected.species)));
   }
 }

 if (is.null(ylim)) {
  ylim = c(minval,maxval)
 }
 
 plot(s[[1]]$dist,s[[1]]$spx,pch=pch[1],col=col[1],cex=cex[1], xlim=xlim, ylim=ylim,xlab=xlab,ylab=ylab,...)
 if (length(selected.species) > 1) {
  for (si in 2:length(selected.species)) {
   points(s[[si]]$dist,s[[si]]$spx,pch=pch[si],col=col[si],cex=cex[si])
  }
 }
}










initializeN = function(seed=1,density.per.area=1,species=1,U=1,dim=2,poisson=F) {
 if (dim == 1) {
   if (poisson) {
    n=rpois(1,U*density.per.area)
   }
   else {
    n=U*density.per.area
   }
   x=matrix(0,nrow=n,ncol=2)
   x[,1] = species
   set.seed(seed,kind = "Mersenne-Twister", normal.kind = "Inversion")
   x[,2] = runif(n,min=0,max=U)
 }
 else if (dim == 2) {
   if (poisson) {
    n=rpois(1,U*U*density.per.area)
   }
   else {
    n=U*U*density.per.area
   }
   x=matrix(0,nrow=n,ncol=3)
   x[,1] = species
   set.seed(seed,kind = "Mersenne-Twister", normal.kind = "Inversion")
   x[,2] = runif(n,min=0,max=U)
   set.seed(seed+1,kind = "Mersenne-Twister", normal.kind = "Inversion")
   x[,3] = runif(n,min=0,max=U)
 }
 else {
  stop(sprintf("dimension should be 1 or 2 (was %d).\n",dim));
 }
 x
}


get.coords.one.species.singletime = function(a=table,sp) {
  coord=a[,2:ncol(a)]
  species=a[,1]
  cond = which(species == sp)

if(ncol(a)-1==2){
  x.coord = coord[cond,][,1]
  y.coord = coord[cond,][,2]
  coordResult = cbind(x.coord,y.coord)
} 
else if (ncol(a)-1==1)
{
  x.coord = coord[cond,][,1]
  coordResult = cbind(x.coord)
} 
else 
{  stop(sprintf("dimension should be 1 or 2 (was %d).\n",dim));
}
}




