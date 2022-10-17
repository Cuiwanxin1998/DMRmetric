#' Generate DMRs data that can be used for metric calculations
#'
#' This function generatea DMR data that includes chromosome, starting and end position, DMR length, and probe distribution
#'
#' @param DMR the raw DMR that need convert
#' @param chr the column index of chromosome in the raw DMR and The format should be character type 1,2,3...,X,Y
#' @param start the column index of start in the raw DMR
#' @param end the column index of end in the raw DMR
#' @param minProbeLength The minimum number of probes the DMR should have
#' @export
#' @examples
#' DMR_file = system.file("extdata","raw_DMR.csv",package = 'DMRmetric')
#' DMR <- read.csv(DMR_file)
#' DMR$seqnames = as.character(DMR$seqnames)
#' convert_DMR <- preinput(DMR, chr=3, start=4, end=5, minProbeLength=3)

preinput <-
  function(DMR,chr,start,end,minProbeLength=3){
    cat(sprintf("[%s] Generate DMRs data that can be used for metric calculations #\n", Sys.time()))
    cat("Obtain the probes included in the DMR...\n")
    illumina450k <- DMRmetric::illumina450k_hg19
    if(typeof(DMR[,chr]) != 'character'){
      stop("the format of chromosome should be character  type '1', '2', '3'..., 'X', 'Y'")
    }
    probe <- apply(DMR,1,function(data){
      cni <- data[chr]
      St <- data[start]
      Ed <- data[end]
      tmp <- illumina450k[illumina450k$CHR==as.character(cni),]
      index <- tmp[as.integer(St) <= tmp$MAPINFO & as.integer(Ed) >= tmp$MAPINFO,]
      probe <- paste(rownames(index), collapse=', ')
      length <- length(rownames(index))
      return(cbind(probe,length))
    })
    probe = t(probe)
    probe_length = as.integer(probe[,2])
    probe = probe[,1]

    res <- as.data.frame(cbind(chr=DMR[,chr],start=DMR[,start],end=DMR[,end],length=DMR[,end]-DMR[,start],probe=probe, probe_length=probe_length))
    colnames(res) <- c("chr",'start','end','length','probe',"probeLength")
    res$chr <- as.character(res$chr)
    res$start <- as.integer(res$start)
    res$end <- as.integer(res$end)
    res$length <- as.integer(res$length)
    res$probeLength <- as.integer(res$probeLength)
    cat(sprintf("remove the number of probe in DMRs < %i ...\n",minProbeLength))
    res <- res[res$probeLength >= minProbeLength,]
    cat("Output converted DMR...\n")
    cat(sprintf("[%s] Done\n", Sys.time()))
    res
  }

