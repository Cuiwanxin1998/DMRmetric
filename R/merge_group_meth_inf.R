#' Merge the methylation report files of samples in each group into one methylation report file.
#'
#' 	It provides two ways of merging for each cytosine: 1. count based: merging the numbers of methylated and unmethylated cytosines from each sample; methylation level based: calculating the average methylation level based on the methylation levels from samples.
#'
#' @param list_methylation_report_file  A file containing the methylation report files of samples with group id.
#' @param meth   the column id to specify the meth counts.
#' @param unmeth the column id to specify the unmeth counts.
#' @export


merge_group_meth_inf <-function(list_methylation_report_file,meth,unmeth){
  cat(sprintf("[%s] Merge the methylation report files of samples in each group into one methylation report file #\n", Sys.time()))
  merge = as.data.frame(matrix(0,0,10))
  Group_name = names(table(list_methylation_report_file[,2]))
  cat("Groups:",paste0(Group_name[1],' and ', Group_name[2]),'\n')
  A.index = which(list_methylation_report_file[,2] == Group_name[1])
  B.index = which(list_methylation_report_file[,2] == Group_name[2])
  cat("process Group",Group_name[1],"......\n")
  for(i in 1:length(A.index)){
    cat("process file",list_methylation_report_file[A.index[i],1],"......\n")
    methy_report = read.table(list_methylation_report_file[A.index[i],1],sep='\t',header=T)
    methy_report = methy_report[which((methy_report[,meth] + methy_report[,unmeth]) > 0),]
    rownames(methy_report) = paste0(methy_report[,1],'_',methy_report[,2])
    if(i == 1){
      merge_A = as.data.frame(cbind(Chr=methy_report[,1],pos=methy_report[,2],
                                    Sample_num_A=1,
                                    avg_meth_level_A=methy_report[,meth]/(methy_report[,meth]+methy_report[,unmeth]),
                                    meth_A=methy_report[,meth],
                                    unmeth_A=methy_report[,unmeth]))
      rownames(merge_A) = paste0(methy_report[,1],'_',methy_report[,2])
    }else{
      merge.index = match(rownames(methy_report),rownames(merge_A))
      no_exist = which(is.na(merge.index))
      if(length(no_exist) != 0){
        addition = as.data.frame(cbind(Chr=methy_report[no_exist,1],pos=methy_report[no_exist,2],
                                       Sample_num_A=0,
                                       avg_meth_level_A=0,
                                       meth_A=0,
                                       unmeth_A=0))
        rownames(addition) = paste0(methy_report[no_exist,1],'_',methy_report[no_exist,2])
        merge_A = rbind(merge_A, addition)
      }

      merge.index = intersect(rownames(methy_report),rownames(merge_A))
      merge_A[merge.index,3] = as.integer(merge_A[merge.index,3])+ 1
      merge_A[merge.index,4]= as.numeric(merge_A[merge.index,4]) + methy_report[merge.index,meth]/(methy_report[merge.index,meth]+methy_report[merge.index,unmeth])
      merge_A[merge.index,5] = as.integer(merge_A[merge.index,5]) + methy_report[merge.index,meth]
      merge_A[merge.index,6]= as.integer(merge_A[merge.index,6]) + methy_report[merge.index,unmeth]
    }
  }

  cat("process Group",Group_name[2],"......\n")
  for(i in 1:length(B.index)){
    cat("process file",list_methylation_report_file[B.index[i],1],"......\n")
    methy_report = read.table(list_methylation_report_file[B.index[i],1],sep='\t',header=T)
    methy_report = methy_report[which((methy_report[,meth] + methy_report[,unmeth]) > 0),]
    rownames(methy_report) = paste0(methy_report[,1],'_',methy_report[,2])
    if(i == 1){
      merge_B = as.data.frame(cbind(Chr=methy_report[,1],pos=methy_report[,2],
                                    Sample_num_B=1,
                                    avg_meth_level_B=methy_report[,meth]/(methy_report[,meth]+methy_report[,unmeth]),
                                    meth_B=methy_report[,meth],
                                    unmeth_B=methy_report[,unmeth]))
      rownames(merge_B) = paste0(methy_report[,1],'_',methy_report[,2])
    }else{
      merge.index = match(rownames(methy_report),rownames(merge_B))
      no_exist = which(is.na(merge.index))
      if(length(no_exist) != 0){
        addition = as.data.frame(cbind(Chr=methy_report[no_exist,1],pos=methy_report[no_exist,2],
                                       Sample_num_B=0,
                                       avg_meth_level_B=0,
                                       meth_B=0,
                                       unmeth_B=0))
        rownames(addition) = paste0(methy_report[no_exist,1],'_',methy_report[no_exist,2])
        merge_B = rbind(merge_B, addition)
      }

      merge.index = intersect(rownames(methy_report),rownames(merge_B))
      merge_B[merge.index,3] = as.integer(merge_B[merge.index,3])+ 1
      merge_B[merge.index,4]= as.numeric(merge_B[merge.index,4]) + methy_report[merge.index,meth]/(methy_report[merge.index,meth]+methy_report[merge.index,unmeth])
      merge_B[merge.index,5] = as.integer(merge_B[merge.index,5]) + methy_report[merge.index,meth]
      merge_B[merge.index,6]= as.integer(merge_B[merge.index,6]) + methy_report[merge.index,unmeth]
    }
  }
  match_A = match(rownames(merge_A),rownames(merge_B))
  match_B = match(rownames(merge_B),rownames(merge_A))
  all = intersect(rownames(merge_A),rownames(merge_B))
  merge = as.data.frame(cbind(chr=merge_A[all,1],pos=merge_A[all,2],Sample_numA=merge_A[all,3],
                              avg_meth_levelA=merge_A[all,4],meth_A = merge_A[all,5],unmeth_A=merge_A[all,6],
                              Sample_numB=merge_B[all,3],
                              avg_meth_levelB=merge_B[all,4],meth_B=merge_B[all,5],unmeth_B=merge_B[all,6]))
  only_A = which(is.na(match_A))
  only_B = which(is.na(match_B))
  addition_A = as.data.frame(
    cbind(
      chr = merge_A[only_A, 1], pos = merge_A[only_A, 2],
      Sample_numA = merge_A[only_A, 3], avg_meth_levelA = merge_A[only_A, 4],
      meth_A = merge_A[only_A, 5],unmeth_A = merge_A[only_A,6],
      Sample_numB=0,avg_meth_levelB = NA, meth_B = NA,unmeth_B = NA
    )
  )
  addition_B = as.data.frame(
    cbind(
      chr = merge_B[only_B, 1], pos = merge_B[only_B, 2],
      Sample_numA = 0, avg_meth_levelA = NA,
      meth_A =NA,unmeth_A = NA,
      Sample_numB=merge_B[only_B, 3],avg_meth_levelB = merge_B[only_B, 4],
      meth_B = merge_B[only_B, 5],unmeth_B = merge_B[only_B,6]
    )
  )
  cat("merge Group",Group_name[1],"and",Group_name[2],"......\n")
  merge=rbind(merge,addition_A,addition_B)
  colnames(merge) = c("Chr",'pos',paste0("Sample_num_",Group_name[1]),
                      paste0("ave_meth_level_",Group_name[1]),
                      paste0("meth_",Group_name[1]),
                      paste0("unmeth",Group_name[1]),
                      paste0("Sample_num_",Group_name[2]),
                      paste0("ave_meth_level_",Group_name[2]),
                      paste0("meth_",Group_name[2]),
                      paste0("unmeth",Group_name[2]))
  merge$pos = as.integer(merge$pos)
  merge=merge[order(merge$Chr,merge$pos),]
  write.table(merge,paste0("merge_",Group_name[1],"_",Group_name[2],'.txt'),quote=F,row.names=F,col.names=T,sep='\t')
}


