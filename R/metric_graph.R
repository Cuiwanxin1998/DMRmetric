#' Draw a line chart for DMR metric.
#'
#' This function will draw a line chart for Qn, Ql and DMRn.
#'
#' @param metric The text file containing the metric DMRn file path and the name of DMR detection method.
#' @param threshold The methylation difference threshold t, default 0-0.9.
#' @param Title The title of rank graph, default NULL.
#' @export
#' @import ggplot2
#' @examples
#' library("ggplot2")
#' metric_graph(metric.DMRn,type="DMRn",Title='GroupA-GroupB')
metric_graph <- function(metric, type = c("Qn","Ql","DMRn"),Title=NULL){
  cat(sprintf("[%s] # Draw a line chart for metric. #\n", Sys.time()))
  type <- match.arg(type)
  if(type == 'Qn'){
    cat("The metirc Qn will be draw...")
  }else if(type == 'Ql'){
    cat("The metric Ql will be draw...")

  }else if(type == 'DMRn'){
    cat("The metric DMRn will be draw...")
  }else{
    stop("The type of metric must be Qn, Qf or DMRn, please check your input...")
  }
  group = colnames(metric)
  threshold = rownames(metric)
  cat("There are", length(group),'methods in the metric...\n')
  DMR_data = as.data.frame(matrix(0,0,3))

  if(is.null(Title)){
    Title='Title'
    cat("The title is NULL...\n")
  }
  for(i in 1:ncol(metric)){
    DMR_data = rbind(DMR_data,cbind(metric=metric[,i],Group=group[i],threshold=threshold))

  }
  DMR_data$metric = as.numeric(DMR_data$metric)
  if(type == "Qn"){
    p = ggplot(data = DMR_data, aes(x = threshold, y = metric)) + geom_line(aes(group =
                                                                                  Group, color = Group)) + geom_point(aes(color = Group, shape = Group)) +
      scale_shape_manual(values = 0:nrow(metric)) +
      theme(panel.background = element_blank(), axis.line = element_line(), axis.text.y=element_text()) +
      ggtitle(Title) + theme(plot.title = element_text(hjust = 0.5)) +
      xlab("methylation difference threshold t") + ylab(expression(italic('Qn')))
    cat("The plot is stored in the path:",getwd(),"as pdf file...\n")
    cat(sprintf("[%s] Done\n", Sys.time()))
    ggsave(p,file='Qn_graph.pdf')
  }
  if(type == "Ql"){
    p = ggplot(data = DMR_data, aes(x = threshold, y = metric)) + geom_line(aes(group =
                                                                                  Group, color = Group)) + geom_point(aes(color = Group, shape = Group)) +
      scale_shape_manual(values = 0:nrow(metric)) +
      theme(panel.background = element_blank(), axis.line = element_line(), axis.text.y=element_text()) +
      ggtitle(Title) + theme(plot.title = element_text(hjust = 0.5)) +
      xlab("methylation difference threshold t") + ylab(expression(italic('Ql')))
    cat("The plot is stored in the path:",getwd(),"as pdf file...\n")
    cat(sprintf("[%s] Done\n", Sys.time()))
    ggsave(p,file='Ql_graph.pdf')
  }
  if(type == 'DMRn'){
    p = ggplot(data = DMR_data, aes(x = threshold, y = metric)) + geom_line(aes(group =
                                                                                  Group, color = Group)) + geom_point(aes(color = Group, shape = Group)) +
      scale_shape_manual(values = 0:nrow(metric)) +
      theme(panel.background = element_blank(), axis.line = element_line(), axis.text.y=element_text()) +
      ggtitle(Title) + theme(plot.title = element_text(hjust = 0.5)) +
      xlab("methylation difference threshold t") + ylab(expression(paste('log'[2],"(",italic('DMRn'),")")))
    cat("The plot is stored in the path:",getwd(),"as pdf file...\n")
    cat(sprintf("[%s] Done\n", Sys.time()))
    ggsave(p,file='DMRn_graph.pdf')
  }

}

