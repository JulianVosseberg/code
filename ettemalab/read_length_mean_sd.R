#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
reads<-read.csv(file=args[1], sep="", header=FALSE)
av_length<-mean(rep(reads$V2,reads$V1))
print(av_length)
st_dev<-sd(rep(reads$V2,reads$V1))
print(st_dev)
