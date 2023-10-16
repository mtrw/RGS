

l_ply(unique(fai$chr),function(c){
  l_ply(seq(from=1,to=fai[chr==c]$length-increment,by=increment), function(s){ #COMMENT OUT FOR PARALLEL PROCESSING OF ALIGNMENTS
  #mclapply(mc.cores=8,seq(from=1,to=fai[chr==c]$length-increment,by=increment), function(s){ #UNCOMMENT FOR PARALLEL, SET mc.cores.
    # dev c <- fai[1]$chr; s <- seq(from=1,to=fai[chr==c]$length-increment,by=increment)[2]
    rg <- c(s,s+windowsize)
    # dev rg <- c(s,s+2000000)
    
    fname <- paste0("data/alignments/1",c,"_",rg[1],"-",rg[2],"_bw",bw,"_dbw",dbw,"_winsize",windowsize,"_ncalcpts",qgap,".dots")
    ce("Searching for ",fname)
    if(!file.exists(fname) | file.size(fname) == 0){
      ce("\tFile missing or miniscule! (",fname,"); Making alignment ...")
      get_lastz_dotplot(
        subjectFile  = refFname,
        subjectSeq   = c,
        subjectRange = rg,
        args = paste(lastzArgs,"--format=general:name1,start1,end1,strand1,length1,name2,start2,end2,strand2,length2,id%,blastid%"),
        selfAlignment = T,
        lastz_binary=lastzBin,
        save_dots_to_file=fname,
        plot=F
      )
    } else {
      ce(fname," already exists, skipping.")
    }
  })
}) %>% invisible



