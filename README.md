# DMRmetric

This is an R package to evaluate those DMRs(differential methylation regions) obtained using 450K methylation array data, contains four main functions including preinput(), merge_group_meth_inf(), compute_QnQl(), compute_DMRn() and metric_graph().


The following Running the test sections will guide the user how to use the R package.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

This R package is based on R-4.1.3

```
You need install R-4.1.3
```

### Installing

Please enter R CMD, user can input the follwing command to install our R package.

```
install.packages("devtools")
install.packages("ggplot2")
library("devtools")
install_github("Cuiwanxin1998/DMRmetric ")
```

## Running the tests

### calculate_QnQl

This function evaluate the DMR predicted from the BS-seq data by calculating the scores for the metric Qn and Ql. Two input file are required: a DMR list file and a merged BS-seq data file. The DMR list file contains information for the storage path of DMR predicted by different methods and the name of each method.

```
DMRfile_list = read.table(system.file('extdata','DMRfile.BSseq.list',package='DMRmetric'))
for(i in 1:nrow(DMRfile_list)){DMRfile_list[i,1] = eval(parse(text=DMRfile_list[i,1]))}
merged_example_file = readRDS(system.file('extdata','BSseq.merge.RDS',package='DMRmetric'))'
metric.QnQl = calculate_QnQl(DMRfile_list,merged_example_file,method="1")
```



### preinput

This function generatea DMR data that includes chromosome, starting and end position, DMR length, and probe distribution.

```
library("DMRmetric")
DMR_file = system.file("extdata","raw_DMR.csv",package = 'DMRmetric')
DMR <- read.csv(DMR_file)
DMR$seqnames = as.character(DMR$seqnames)
convert_DMR <- preinput(DMR, chr=3, start=4, end=5, minProbeLength=3)
```

### compute_DMRn

This function evaluate the DMR predicted from the methylation array data by calculating the scores for
the metric DMRn.

```
DMRfile_list = system.file("extdata","DMRfile.Array.list.txt",package = 'DMRmetric')
beta_file = system.file("extdata","Array_beta.RDS",package = 'DMRmetric')
Array_beta = readRDS(beta_file)
for(i in 1:nrow(DMRfile_list)){DMRfile_list[i,1] = eval(parse(text=DMRfile_list[i,1]))}
case = colnames(Array_beta)[1:5] 
control = colnames(Array_beta)[6:10]
metric.DMRn = calculate_DMRn(DMRfile_list = DMRfile_list,beta = Array_beta, human_ref = 'Tissue', Case_name=case, Control_name=control)
```

### metric_graph

This function will draw a line chart for Qn, Ql and DMRn.

```
library("ggplot2")
metric_graph(metric.DMRn,type="DMRn",Title='GroupA-GroupB')
```


## Built With

* [R](https://www.r-project.org/) - R is a free software environment for statistical computing and graphics.

## Authors

* **wx Cui** - *Thesis code writing work* [Central South University](https://cse.csu.edu.cn/)


