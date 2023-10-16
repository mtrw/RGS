setwd("/data/gpfs/projects/punim1869/users/mrabanuswall/workspace/RGS")

source("0_llprIdent_setup.R")

require(microbenchmark)

microbenchmark(
  l_ply(unique(fai$chr),function(c){
    l_ply(seq(from=1,to=fai[chr==c]$length-increment,by=increment), function(s){ #COMMENT OUT FOR PARALLEL PROCESSING OF ALIGNMENTS
    #mclapply(mc.cores=8,seq(from=1,to=fai[chr==c]$length-increment,by=increment), function(s){ #UNCOMMENT FOR PARALLEL, SET mc.cores.
      # dev c <- fai[1]$chr; s <- seq(from=1,to=fai[chr==c]$length-increment,by=increment)[2]
      rg <- c(s,s+windowsize)
      
      fname <- paste0("data/alignments/",c,"_",rg[1],"-",rg[2],"_bw",bw,"_dbw",dbw,"_winsize",windowsize,"_ncalcpts",qgap,".dots")
      ce("Searching for ",fname)
      if(!file.exists(fname) | file.size(fname) == 0){
        ce("\tFile missing or miniscule! (",fname,"); Making alignment ...")
        get_lastz_dotplot(
          subjectFile  = refFname,
          subjectSeq   = c,
          subjectRange = rg,
          args = paste(lastzArgs,"--format=general:name1,start1,end1,strand1,length1,name2,start2,end2,strand2,length2,id%,blastid%,cov%"),
          selfAlignment = T,
          lastz_binary=lastzBin,
          save_dots_to_file=fname, #give it a filename to save to if you want.
          plot=F
        )
      } else {
        ce(fname," already exists, skipping.")
      }
    })
  }) %>% invisible
)
























l_ply(unique(fai$chr),function(c){
  #l_ply(Sys.glob(paste0("plots/dots_files/*",c,"*")), function(fname){
  l_ply(seq(from=1,to=fai[chr==c]$length,by=increment)), function(s){
    # dev c <- fai[1]$chr; s <- 1800001
    #s <- sub(".*_(\\d*)-\\d*.*","\\1",fname) %>% as.integer
    rg <- c(s,s+windowsize)
    
    fname <- paste0("data/alignments/",c,"_",rg[1],"-",rg[2],"_bw",bw,"_dbw",dbw,"_winsize",windowsize,"_ncalcpts",qgap,".dots")
    ce("Working on ",fname)
    if(!file.exists(fname) | file.size(fname) == 0){
      stop("\tFile missing or miniscule! (",fname,")")
    }
    ce("\tReading alignments")
    if(nrow(t)==0){
      ce("No alignments")
      return(NULL)
    }
    ce("\tFiltering alignments")
    
    t[,gp:=rep(1:(.N/3),each=3)]
    setkey(t,gp)
    t[,fromdiag:=abs(x[1]-y[1]),by="gp"]
    t[,length:=max(abs(x[1]-x[2]),abs(y[1]-y[2])),by="gp"]
    
    tfplot <- t[fromdiag != 0 & fromdiag <= 1000000 & length >= ll,.(x,y,fromdiag)]
    tfplot[,col:=ifelse(fromdiag <= bw,"black","red")]
    
    tf <- t[fromdiag != 0 & fromdiag <= bw & length >= ll,.(x,y)]
    
    if(nrow(tf)==0){
      ce("No alignments after filtering")
      return(NULL)
    }
    
    ce("\tCasting to overlap-detection-friendly format")
    a <- data.table(chr="C",start=tf[seq(1,.N,by=3)]$x,end=tf[seq(2,.N,by=3)]$x) #aligned regions in format for intersection by foverlaps
    a[,length := end-start]
    
    #query points
    s <- seq(rg[1],rg[2],by=qgap)
    q <- data.table(chr="C",start=s,end=s)
    setkey(a,chr,start,end) #overlap based on 'subject' position
    setkey(q,chr,start,end)
    
    ce("\tPerforming overlap")
    ol <- foverlaps(q,a,type="any",nomatch=0L)
    if(nrow(ol)==0){
      ce("No query point overlaps")
      return(NULL)
    }
    ol[,weight:=length]
    
    ce("\tCalculating density")
    d <- density(ol$i.start,weights=ol$weight,bw=dbw,from=rg[1],to=rg[2],n=round(diff(rg)/qgap)) %>% suppressWarnings()
    
    density_range <- range(d$y)
    weights_per_density <- sum(ol$weight)/sum(d$y)
    
    cutoff_scaled <- rg[1]+diff(rg)*(cutoff/(density_range[2]-density_range[1]))
    cutoff_scaled_under <- cutoff_scaled-diff(rg)*0.01
    
    if(doplot==T){
      ce("\tBegin plotting")
      plot(tfplot[col=="black"]$x,tfplot[col=="black"]$y,col=tfplot[col=="black"]$col,type="l",ylab="Dotplot (after band and length filters)",xlab=NA,main=fname,xlim=rg,ylim=rg)
      lines(tfplot[col=="red"]$x,tfplot[col=="red"]$y,col=tfplot[col=="red"]$col)
      lines(d$x,d$y %>% scale_between(rg[1],rg[2]),col="#4564ed88")
      abline(h=c(cutoff_scaled),col="#3dad8faa")
      text(x=left(rg),y=top(rg),adj=c(0,0.5),labels=paste0("Density range: ",density_range[1],"-",density_range[2],"\nCutoff: ",cutoff)) #,"\nWeights/density: ",weights_per_density))
    }
    ce("\tExtracting regions of interest")
    #extract regions
    min_runlength <- 0
    rled <- (d$y > cutoff) %>% rle
    rled$lengths
    rled$values
    rled$values[rled$lengths < min_runlength & rled$values==FALSE] <- TRUE
    rled <- rle(inverse.rle(rled))
    rled_csum <- cumsum(rled$lengths)
    regions <- data.table(
      chr=c,
      state=rled$values,
      start=d$x[c(1,rled_csum[1:(length(rled$values)-1)]+1)] %>% as.integer,
      end=d$x[rled_csum] %>% as.integer
    )[state==TRUE][,state:=NULL][,ord:=1:.N][,start:=start-qgap][,end:=end+qgap][]
    
    #overlap regions with alignments and expand to maximum span of their start/end positions
    if(nrow(regions)==0){
      ce("No regions passed filters")
      if(doplot==TRUE){
        gplots::textplot(data.table(message="No regions passed filter"))
      }
      return(NULL)
    }
    
    seedregions <- copy(regions)
    a <- data.table(chr="C",xstart=tf[seq(1,.N,by=3)]$x,xend=tf[seq(2,.N,by=3)]$x,ystart=tf[seq(1,.N,by=3)]$y,yend=tf[seq(2,.N,by=3)]$y) #aligned regions in format for intersection by foverlaps
    a[,alnid:=1:.N]
    
    q <- regions[,.(chr="C",regchr=chr,xstart=start,xend=end,regid=ord)]
    setkey(q,chr,xstart,xend)
    setkey(a,chr,xstart,xend)
    ol_ra <- foverlaps(a,q,type="any",nomatch=0L)
    ol_ra <- a[,.(alnid,alnxstart=xstart,alnxend=xend,alnystart=ystart,alnyend=yend)][ol_ra,on="alnid"][,chr:=NULL][]
    
    #overlap regions with alignments and expand to maximum span of their start/end positions
    if(nrow(ol_ra)==0){
      #the seed region falls between two alignments but is too skinny to overlap either.
      ce("No alignment / region seed overlap, figuring it out ...")
      ol_ra <- q[,{
        addbefore <- a[xstart <= beforeme][xstart==max(xstart)]
        addafter <-  a[xend >= afterme][xend==min(xend)]
        rbindlist(list(
          addbefore,
          addafter
        ))
      },by=.(regid,regchr,beforeme=xstart,afterme=xend)][,.(alnid,regid,regchr,alnxstart=xstart,alnxend=xend,alnystart=ystart,alnyend=yend)]
    }
    
    regions <- ol_ra[,.(start=min(c(alnxstart,alnystart)),end=max(c(alnxend,alnyend)) ),by=.(ord=regid,chr=regchr)]
    
    if(doplot==TRUE){
      ce("\tMelting region data to plot")
      pregions <- melt(regions,measure.vars=c("start","end"),id.vars=c("chr","ord"))[order(ord,variable),idx:=rep(1:(.N/2),each=2)][,.SD[1:3],by="idx"]
      lines(pregions$value,rep(cutoff_scaled,nrow(pregions)),col="#8c0606",lwd=2)
      pseedregions <- melt(seedregions,measure.vars=c("start","end"),id.vars=c("chr","ord"))[order(ord,variable),idx:=rep(1:(.N/2),each=2)][,.SD[1:3],by="idx"]
      lines(pseedregions$value,rep(cutoff_scaled_under,nrow(pregions)),col="#1080e8",lwd=2)
    }
    
    #list genes in there also, purely for plotting at this point
    if(doplot==TRUE){
      ce("\tPlotting gene overlaps with array regions")
      setkey(gff1,seqname,start,end)
      setkey(regions,chr,start,end)
      genes <- foverlaps(regions,gff1,type="any")[feature=="mRNA",.(chr,start,end,strand,id,description,region_start=i.start,region_end=i.end,attribute)]
      abline(
        v=genes[,pmean(start,end)],
        col="#bd040422"
      )
      wait()
    }
    
    #check genes are repeated
    tfc <- tf[!is.na(x)][,t := rep(c("start","end"),times=.N/2)][,idx:=rep(1:(.N/2),each=2)][,.(idx,x,chr=c,t)]
    tfc <- dcast(tfc,idx+chr~t,value.var = "x")
    setkey(tfc,chr,start,end)
    setkey(genes,chr,start,end)

    genes_pass <- foverlaps(genes,tfc,type="any")[!is.na(idx),.N,by="id"]$id %>% unique #only if they overlap a (non-diagonal) alignment
    genes <- genes[id %in% genes_pass]
    if(doplot==TRUE){
      abline(
        v=genes[,pmean(start,end)],
        col="#bd0404",
        lwd=0.5
      )

    }
    
    #output
    print(regions)
    out_regions <<- bind_rows(out_regions,regions)
    out_genes <<- bind_rows(out_genes,genes)
    
    if(doplot==TRUE){
      if(nrow(genes)>0){
        gplots::textplot(genes[,.(chr,start,end,strand,id,description)])
      } else {
        gplots::textplot(data.table(message="No genes recorded"))
      }
    }
    return(NULL)
  })
  return(NULL)
}) %>% invisible


































#merge regions
out_regions[,ord:=NULL]
setorder(out_regions,chr,start)
out_regions_tmp <- copy(out_regions)
setorder(out_regions_tmp,chr,start)

r <- 1
while(r==1){
  print("################ NEXT ROUND ##################")
  r <- 0
  i <- 1
  while( i < nrow(out_regions_tmp) ){
    #print(i)
    if(out_regions_tmp[i]$chr != out_regions_tmp[i+1]$chr){
      i <- i + 1
      next
    }
    if(out_regions_tmp[i]$end >= out_regions_tmp[i+1]$start){
      #wait()
      r <- 1
      print(paste0("MERGE: ",i))
      print(out_regions_tmp[i:(i+1)])
      
      out_regions_tmp[i:(i+1)] <- data.table(
        chr=c(out_regions_tmp[i]$chr,NA),
        start=c(out_regions_tmp[i]$start,NA),
        end=c(out_regions_tmp[i]$end,NA)
      )
      
      out_regions_tmp <- out_regions_tmp[!is.na(chr)]
      
      if(i == nrow(out_regions_tmp)){
        ce("Merge on last row, next! (last)")
        next
      } else {
        i <- i + 1
      }
    }
    i <- i+1
  }
}

out_regions_tmp -> out_regions
out_regions[,regionid:=1:.N]
# saveRDS(out_regions,"data/rds/out_regions.rds")
out_regions <- readRDS("data/rds/out_regions.rds")
