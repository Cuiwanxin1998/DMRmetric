#' Calculation of differential methylation values for DMRs
#'
#' This function calculates the differential methylation value of the DMR and the number of probes and CpGs included and the length of the DMR
#'
#' @param DMR the converted DMR that converted by function preinput
#' @param beta the beta matrix with case group and control group(The case group is in the front, the group group is in the back)
#' @param human_ref Human genome reference strand from different cancer tissue or different tissue
#' @param Case_name the colnames of case group
#' @param Control_name the colnames of control group
#' @export
#' @examples
#' DMRfile_list = read.table(system.file("extdata","DMRfile.Array.list.txt",package = 'DMRmetric'))
#' beta_file = system.file("extdata","Array_beta.RDS",package = 'DMRmetric')
#' Array_beta = readRDS(beta_file)
#' for(i in 1:nrow(DMRfile_list)){DMRfile_list[i,1] = eval(parse(text=DMRfile_list[i,1]))}
#' case = colnames(Array_beta)[1:5]
#' control = colnames(Array_beta)[6:10]
#' metric.DMRn = calculate_DMRn(DMRfile_list = DMRfile_list,beta = Array_beta, human_ref = 'Tissue', Case_name=case, Control_name=control)
#'
calculate_DMRn <-
  function(DMRfile_list, beta, human_ref=c('Cancer',"Tissue"), Case_name, Control_name){
    human_ref <- match.arg(human_ref)
    DMRn = as.data.frame(matrix(0,10,0))
    cat(sprintf("[%s] # Calculate the number of cpg and the difference methylation #\n", Sys.time()))
    if(human_ref == 'Cancer'){
      ref = DMRmetric::RefCancer
    }else if(human_ref == 'Tissue'){
      ref = DMRmetric::RefTissue
    }else{
      stop("The human reference must be Cancer or Tissue...")
    }
    cat("Calculate the difference methylation value for each probe...\n")
    prb_diff <- apply(beta[,Case_name], 1, mean) - apply(beta[,Control_name], 1, mean)
    cat("Calculate the difference methylation value for each DMR...\n")
    list_res = list()
    for(k in 1:nrow(DMRfile_list)){
      res <- as.data.frame(matrix(NA, 0, 8))
      cat("process the file", paste0(DMRfile_list[k,1],","),"belong to method", DMRfile_list[k,2],'...\n')
      DMR <- read.table(DMRfile_list[k,1],header=T,sep='\t')
      for(i in 1:nrow(DMR)){
        match <- intersect(strsplit(DMR[i,5], split=', ')[[1]], rownames(beta))
        DMR_beta <- beta[match, ]
        probe <- rownames(DMR_beta)
        manchr <- ref[probe,]
        numprbs <- length(probe)
        manchr$diff <- prb_diff[probe]
        manchr <- manchr[order(manchr$CHR,manchr$MAPINFO),]

        chr <- DMR[i,1]
        start <- DMR[i,2]
        end <- DMR[i,3]
        DMRlength<- DMR[i,4]
        molecular <- sum(manchr[,3]*manchr[,5])
        Denominator <- sum(manchr[,4])
        diff <- molecular/Denominator
        diff1 <- mean(prb_diff[probe])
        res=rbind(res,c(chr,start,end,numprbs,Denominator,DMRlength,diff,diff1))
      }
      colnames(res) <- c("chr","Start", "End", "prb","cpg", "length", "Diff1",'Diff2')
      res$Start <- as.integer(res$Start)
      res$End <- as.integer(res$End)
      res$prb <- as.integer(res$prb)
      res$cpg <- as.integer(res$cpg)
      res$length <- as.integer(res$length)
      res$Diff1 <- as.numeric(res$Diff1)

      write.table(res,paste0(DMRfile_list[k,2],'.diff.res'),row.names = F,col.names=T,quote=F,sep='\t')
      res = res[!is.na(res$Diff1),]
      list_res = append(list_res,list(res))
    }

    interval <- as.data.frame(matrix(NA,10,3))
    colnames(interval) <- c("Tprb","Tcpg", "Tlength")
    t = 0
    for(i in 1:nrow(DMRfile_list)){
      tmp_res = list_res[[i]]
      t = max(t,max(tmp_res[,7]))
    }
    count = 0
    while(t < 1){
      t = t*10
      count = count + 1
    }

    cat("Calculate the total number of cpg, probe and DMR length for each interval...\n")
    for(i in 1:nrow(DMRfile_list)){
      res = list_res[[i]]
      res[,7] = abs(res[,7])
      for(n in c(1:10)){
        t1 <- (n-1)*(0.1^count)
        t2 <- n*(0.1^count)
        tmp <- res[which(res[,7] > t1 & res[,7] <= t2),]
        interval[n,1] <- sum(tmp[,4])
        interval[n,2] <- sum(tmp[,5])
        interval[n,3] <- sum(tmp[,6])

      }
      w <- seq(0.1^count/2, (9*0.1^count+10*0.1^count)/2, 0.1^count)
      DMRn_tmp = c()
      for(m in c(1:10)){
        interval$weight <- seq(0.1^count/2, (9*0.1^count+10*0.1^count)/2, 0.1^count)
        tmp <- interval[m:dim(interval)[1], ]
        tmp <- tmp[tmp$Tprb!=0,]
        if(dim(tmp)[1]==0){
          DMRn_tmp <- append(DMRn_tmp,0)
          next;
        }
        Qn_tmp<-sum(tmp$weight*tmp$Tprb^3*tmp$Tcpg/tmp$Tlength)/sum(w[m:10])
        Qn_tmp <- log2(Qn_tmp+1)
        DMRn_tmp <- append(DMRn_tmp,Qn_tmp)
      }
      DMRn = cbind(DMRn,DMRn_tmp)
    }

    colnames(DMRn) = DMRfile_list[,2]
    rownames(DMRn) = as.character(seq(0,1*0.1^(count-1),0.1^count)[1:10])
    cat(sprintf("[%s] Done\n", Sys.time()))
    DMRn
}


