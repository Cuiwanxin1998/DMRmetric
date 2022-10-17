#' calculate the metric Qn and Ql under different threshold for a DMR set.
#'
#' This function will calculate the methyation difference and count the CpG number and the length in each DMR.
#'
#' @param DMRfile_list a file containing formated DMR files of different methods on the same dataset(for each DMR file, the format is : chr	start_pos	end_pos	...).
#' @param merged_file  a merged methylation file of group A and group B.
#' @param method 1 or 2 (used in calculating the methylation difference, e.g. 1 ---count based; 2 ---methylation level based;   ) .
#' @export
#' @examples
#' DMRfile_list = read.table(system.file('extdata','DMRfile.BSseq.list.txt',package='DMRmetric'))
#' for(i in 1:nrow(DMRfile_list)){DMRfile_list[i,1] = eval(parse(text=DMRfile_list[i,1]))}
#' merged_example_file = readRDS(system.file('extdata','BSseq.merge.RDS',package='DMRmetric'))
#' metric.QnQl = calculate_QnQl(DMRfile_list,merged_example_file,method="1")
#'
calculate_QnQl <- function(DMRfile_list, merged_file, method=c("1","2")){
  cat(sprintf("[%s] # Calculate Qn and Ql under different threshold for a DMR set #\n", Sys.time()))

  if(ncol(DMRfile_list)>2){
    stop("the DMR file list must only have two column, the first column is DMR file path, the second column is the name of DMR method...")
  }

  method <- as.integer(match.arg(method))
  if(method == 1){
    cat("calculating the methylation difference based on count...\n")
  }else if(method == 2){
    cat("calculating the methylation difference based on methylation level...\n")
  }else{
    stop("the method argument must be 1 or 2...")
  }
  Qn = as.data.frame(matrix(0,9,0))
  Ql = as.data.frame(matrix(0,9,0))
  for(k in 1:nrow(DMRfile_list)){
    DMR = read.table(DMRfile_list[k,1])
    cat("process the file", paste0(DMRfile_list[k,1],","),"belong to method", DMRfile_list[k,2],'...\n')
    chrom = names(table(DMR[,1]))
    if(sum(chrom %in% names(table(merged_file[,1]))) == 0){
      stop("chromosome of merge file is not match the DMR file...")
    }
    DMR_difmeth = as.data.frame(matrix(0,0,8))
    total_count = 0
    total_length = 0
    sum_ratio = 0
    minum_cpg=5
    maxium_length=10000
    count_interval = rep(0,11)
    length_interval = rep(0,11)
    Qn_tmp = c()
    Ql_tmp = c()
    for(cni in chrom){
      DMR.temp = DMR[DMR[,1] == cni,]
      merge.temp = merged_file[merged_file[,1] == cni,]

      for(i in 1:nrow(DMR.temp)){
        cpg_interval = rep(0,11)
        names(cpg_interval) = seq(0,1,0.1)
        start = DMR.temp[i,2]
        end = DMR.temp[i,3]
        merge.cpg = merge.temp[start <= merge.temp[,2] & merge.temp[,2] <= end,]
        cpg.num = nrow(merge.cpg)
        sample_A = sum(merge.cpg[,3])
        sample_B = sum(merge.cpg[,7])
        valid.A =  which(merge.cpg[,3] != 0)
        valid.B =  which(merge.cpg[,7] != 0)

        DMR_difmeth.1.meth.A = sum(merge.cpg[valid.A,5])
        DMR_difmeth.1.unmeth.A = sum(merge.cpg[valid.A,6])
        DMR_difmeth.2.A = sum(merge.cpg[valid.A,4] * merge.cpg[valid.A,3])
        DMR_difmeth.1.meth.B = sum(merge.cpg[valid.B,9])
        DMR_difmeth.1.unmeth.B = sum(merge.cpg[valid.B,10])
        DMR_difmeth.2.B = sum(merge.cpg[valid.B,8]* merge.cpg[valid.B,7])
        valid.A.B = intersect(valid.A,valid.B)

        pos.dif = abs(merge.cpg[valid.A.B,4]-merge.cpg[valid.A.B,8])
        pos.dif =as.integer(pos.dif*10)
        cpg = table(pos.dif/10)
        cpg_interval[names(cpg)] = cpg_interval[names(cpg)] + cpg
        cpg_interval[10] = cpg_interval[10] + cpg_interval[11]

        total_count_rank = 0
        sum_ratio = 0
        for(t in seq(0,1,0.1)){
          sum_ratio = sum_ratio + (t + 0.05)
          total_count_rank = total_count_rank + (t+0.05)*cpg_interval[as.character(t)]
        }
        DMR_rank=total_count_rank/sum_ratio
        if(sample_A > 0 & sample_B >0){
          DMR_difmeth.1 = abs(
            DMR_difmeth.1.meth.A / (DMR_difmeth.1.meth.A + DMR_difmeth.1.unmeth.A) - DMR_difmeth.1.meth.B /
              (DMR_difmeth.1.meth.B + DMR_difmeth.1.unmeth.B)
          )
          DMR_difmeth.2 = abs(DMR_difmeth.2.A/sample_A - DMR_difmeth.2.B/sample_B)
        }else{
          DMR_difmeth.1 = 0
          DMR_difmeth.2 = 0
        }
        DMR_difmeth = rbind(
          DMR_difmeth,
          cbind(
            chr = cni,
            start = start,
            end = end,
            length = end - start + 1,
            cpg = cpg.num,
            rank = DMR_rank,
            difmeth1 = DMR_difmeth.1,
            difmeth2 = DMR_difmeth.2
          )
        )
      }

    }
    cat("Output the file", paste0(DMRfile_list[k,1],".Meth_diff,"),"store in the path", paste0(getwd(),"...\n"))
    write.table(DMR_difmeth,paste0(DMRfile_list[k,1],".Meth_diff"),row.names = F,col.names = F,quote=F,sep='\t')
    DMR_difmeth = read.table(paste0(DMRfile_list[k,1],".Meth_diff"))
    cat("select the DMR that has maximum DMR length 10000bp and minumum cpg number is 5...\n")
    DMR_difmeth = DMR_difmeth[which(DMR_difmeth[,4] <= maxium_length & DMR_difmeth[,5] >= minum_cpg),]
    c = as.integer(DMR_difmeth[,6+method]*10)
    for(i in 1:length(c)){
      count_interval[c[i]+1] = count_interval[c[i]+1] + DMR_difmeth[i,5]
      length_interval[c[i]+1] = length_interval[c[i]+1] + DMR_difmeth[i,4]
    }
    count_interval[10] = count_interval[10] + count_interval[11]
    length_interval[10] = length_interval[10] + length_interval[11]
    cat("calculate Qn and Qf of the method",paste0(DMRfile_list[k,2],"...\n"))
    for(i in seq(0.9,0.1,-0.1)){
      sum_ratio = sum_ratio + i + 0.05
      total_count = total_count + (i + 0.05)*count_interval[i*10 + 1]
      total_length = total_length + (i + 0.05)*length_interval[i*10 + 1]

      if(total_count > 0){
        Qn_tmp = append(total_count /sum_ratio,Qn_tmp)
        Ql_tmp = append(total_length / sum_ratio,Ql_tmp)
      }else{
        Qn_tmp = append(0,Qn_tmp)
        Ql_tmp = append(0,Ql_tmp)
      }
    }
    Qn = cbind(Qn,Qn_tmp)
    Ql = cbind(Ql,Ql_tmp)
  }
  colnames(Qn) = DMRfile_list[,2]
  rownames(Qn) = seq(0.1,0.9,0.1)
  colnames(Ql) = DMRfile_list[,2]
  rownames(Ql) = seq(0.1,0.9,0.1)
  cat(sprintf("[%s] Done\n", Sys.time()))
  return(list(Qn=Qn,Ql=Ql))
}


