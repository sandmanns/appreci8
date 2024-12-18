library(Biostrings)
library(stringr)
args = commandArgs()
args = args[length(args)]
args = strsplit(args, ",")[[1]]
dir2<-paste(args[1],"/",sep="")
dir1<-paste(dir2,"/documents/",sep="")
samples<-read.table(paste(dir2,"SampleNames.txt",sep=""),header=F,stringsAsFactors = F)

bam_readcount<-as.character(args[2])
reference<-as.character(args[3])
dir_bam<-as.character(args[4])

start_pos<-1
end_pos<-length(samples[,1])


for(n in start_pos:end_pos){
  message("Sample ",samples[n,1]," (",n," out of ",length(samples[,1]),")")
  input<-read.table(paste(dir1,samples[n,1],".results_raw.txt",sep=""),
                    header=T,stringsAsFactors =F,quote = "",sep="\t")
  input<-data.frame(input[,c(2,3,4,5,1)],stringsAsFactors = F)
  input1<-input[nchar(input[,3])==1&nchar(input[,4])==1,]
  input2<-input[nchar(input[,3])>1|nchar(input[,4])>1,]
  results1<-data.frame()
  
    #A>C
    input1ac<-input1[input1[,3]=="A"&input1[,4]=="C",]

    if(length(input1ac[,1])>0){
      results1ac<-data.frame(input1ac[,c(1:5)],Nr_Ref=NA,Nr_Alt=NA,DP=NA,VAF=NA,BQ_REF=NA,BQ_ALT=NA,
                             Nr_Ref_fwd=NA,Nr_Alt_fwd=NA,DP_fwd=NA,VAF_fwd=NA,
                             Nr_Ref_rev=NA,Nr_Alt_rev=NA,DP_rev=NA,VAF_rev=NA)
      
      write.table(cbind(input1ac[,1],input1ac[,2],input1ac[,2]),paste(dir1,"temp",start_pos,"_",end_pos,".txt",sep=""),
                  quote=F,row.names = F,sep="\t")
      
        message("SNPs: ",length(input1ac[,1]))
        system(paste("tail -n +2 ",dir1,"temp",start_pos,"_",end_pos,".txt > ",dir1,"temp",start_pos,"_",end_pos,"_2.txt",sep=""),intern = F)
        result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                 dir_bam,"/",samples[n,1],".bam -l ",dir1,"temp",start_pos,"_",end_pos,"_2.txt -w 0",sep=""),intern=T)
        result_raw2<-str_split_fixed(result_raw,pattern = "\t",n=Inf)
        message("Evaluation")
        line<-match(results1ac[,2],result_raw2[,2])
        #line<-line[!is.na(line)]
        results1ac[,8]<-as.numeric(result_raw2[line,4])
        results1ac[,6]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,2])
        results1ac[,10]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,4])
        results1ac[,12]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,6])
        results1ac[,16]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,7])
        results1ac[,7]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,2])
        results1ac[,11]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,4])
        results1ac[,13]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,6])
        results1ac[,17]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,7])
        results1ac[,9]<-results1ac[,7]/results1ac[,8]
        for(o in 6:10){
          results1ac[,14]<-rowSums(data.frame(results1ac[,14],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,6])),na.rm=T)
          results1ac[,18]<-rowSums(data.frame(results1ac[,18],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,7])),na.rm=T)   
        }
       results1ac[,15]<-results1ac[,13]/results1ac[,14]
       results1ac[,19]<-results1ac[,17]/results1ac[,18]
       results1<-results1ac
    }
    
    #A>G
    input1ac<-input1[input1[,3]=="A"&input1[,4]=="G",]
    
    if(length(input1ac[,1])>0){
      results1ac<-data.frame(input1ac[,c(1:5)],Nr_Ref=NA,Nr_Alt=NA,DP=NA,VAF=NA,BQ_REF=NA,BQ_ALT=NA,
                             Nr_Ref_fwd=NA,Nr_Alt_fwd=NA,DP_fwd=NA,VAF_fwd=NA,
                             Nr_Ref_rev=NA,Nr_Alt_rev=NA,DP_rev=NA,VAF_rev=NA)
      
      write.table(cbind(input1ac[,1],input1ac[,2],input1ac[,2]),paste(dir1,"temp",start_pos,"_",end_pos,".txt",sep=""),
                  quote=F,row.names = F,sep="\t")
      
      message("SNPs: ",length(input1ac[,1]))
        system(paste("tail -n +2 ",dir1,"temp",start_pos,"_",end_pos,".txt > ",dir1,"temp",start_pos,"_",end_pos,"_2.txt",sep=""),intern = F)
        result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                 dir_bam,"/",samples[n,1],".bam -l ",dir1,"temp",start_pos,"_",end_pos,"_2.txt -w 0",sep=""),intern=T)
        result_raw2<-str_split_fixed(result_raw,pattern = "\t",n=Inf)
        message("Evaluation")
        if(length(result_raw2)>0){
          line<-match(results1ac[,2],result_raw2[,2])
          #line<-line[!is.na(line)]
          results1ac[,8]<-as.numeric(result_raw2[line,4])
          results1ac[,6]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,2])
          results1ac[,10]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,4])
          results1ac[,12]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,6])
          results1ac[,16]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,7])
          results1ac[,7]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,2])
          results1ac[,11]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,4])
          results1ac[,13]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,6])
          results1ac[,17]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,7])
          results1ac[,9]<-results1ac[,7]/results1ac[,8]
          for(o in 6:10){
            results1ac[,14]<-rowSums(data.frame(results1ac[,14],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,6])),na.rm=T)
            results1ac[,18]<-rowSums(data.frame(results1ac[,18],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,7])),na.rm=T)   
          }
            results1ac[,15]<-results1ac[,13]/results1ac[,14]
            results1ac[,19]<-results1ac[,17]/results1ac[,18]
        }
        if(length(results1)>0){
          results1<-rbind(results1,results1ac)  
        }
        if(length(results1)==0){
          results1<-results1ac
        }
    }

    #A>T
    input1ac<-input1[input1[,3]=="A"&input1[,4]=="T",]
    
    if(length(input1ac[,1])>0){
      results1ac<-data.frame(input1ac[,c(1:5)],Nr_Ref=NA,Nr_Alt=NA,DP=NA,VAF=NA,BQ_REF=NA,BQ_ALT=NA,
                             Nr_Ref_fwd=NA,Nr_Alt_fwd=NA,DP_fwd=NA,VAF_fwd=NA,
                             Nr_Ref_rev=NA,Nr_Alt_rev=NA,DP_rev=NA,VAF_rev=NA)
      
      write.table(cbind(input1ac[,1],input1ac[,2],input1ac[,2]),paste(dir1,"temp",start_pos,"_",end_pos,".txt",sep=""),
                  quote=F,row.names = F,sep="\t")
      
      message("SNPs: ",length(input1ac[,1]))
        system(paste("tail -n +2 ",dir1,"temp",start_pos,"_",end_pos,".txt > ",dir1,"temp",start_pos,"_",end_pos,"_2.txt",sep=""),intern = F)
        result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                 dir_bam,"/",samples[n,1],".bam -l ",dir1,"temp",start_pos,"_",end_pos,"_2.txt -w 0",sep=""),intern=T)
        result_raw2<-str_split_fixed(result_raw,pattern = "\t",n=Inf)
        message("Evaluation")
        if(length(result_raw2)>0){
          line<-match(results1ac[,2],result_raw2[,2])
          #line<-line[!is.na(line)]
          results1ac[,8]<-as.numeric(result_raw2[line,4])
          results1ac[,6]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,2])
          results1ac[,10]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,4])
          results1ac[,12]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,6])
          results1ac[,16]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,7])
          results1ac[,7]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,2])
          results1ac[,11]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,4])
          results1ac[,13]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,6])
          results1ac[,17]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,7])
          results1ac[,9]<-results1ac[,7]/results1ac[,8]
          for(o in 6:10){
            results1ac[,14]<-rowSums(data.frame(results1ac[,14],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,6])),na.rm=T)
            results1ac[,18]<-rowSums(data.frame(results1ac[,18],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,7])),na.rm=T)   
          }
       results1ac[,15]<-results1ac[,13]/results1ac[,14]
       results1ac[,19]<-results1ac[,17]/results1ac[,18]
        }
        if(length(results1)>0){
          results1<-rbind(results1,results1ac)  
        }
        if(length(results1)==0){
          results1<-results1ac
        }
    }
    
    #C>A
    input1ac<-input1[input1[,3]=="C"&input1[,4]=="A",]
    
    if(length(input1ac[,1])>0){
      results1ac<-data.frame(input1ac[,c(1:5)],Nr_Ref=NA,Nr_Alt=NA,DP=NA,VAF=NA,BQ_REF=NA,BQ_ALT=NA,
                             Nr_Ref_fwd=NA,Nr_Alt_fwd=NA,DP_fwd=NA,VAF_fwd=NA,
                             Nr_Ref_rev=NA,Nr_Alt_rev=NA,DP_rev=NA,VAF_rev=NA)
      
      write.table(cbind(input1ac[,1],input1ac[,2],input1ac[,2]),paste(dir1,"temp",start_pos,"_",end_pos,".txt",sep=""),
                  quote=F,row.names = F,sep="\t")
      
      message("SNPs: ",length(input1ac[,1]))
        system(paste("tail -n +2 ",dir1,"temp",start_pos,"_",end_pos,".txt > ",dir1,"temp",start_pos,"_",end_pos,"_2.txt",sep=""),intern = F)
        result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                 dir_bam,"/",samples[n,1],".bam -l ",dir1,"temp",start_pos,"_",end_pos,"_2.txt -w 0",sep=""),intern=T)
        result_raw2<-str_split_fixed(result_raw,pattern = "\t",n=Inf)
        message("Evaluation")
        if(length(result_raw2)>0){
          line<-match(results1ac[,2],result_raw2[,2])
          #line<-line[!is.na(line)]
          results1ac[,8]<-as.numeric(result_raw2[line,4])
          results1ac[,6]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,2])
          results1ac[,10]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,4])
          results1ac[,12]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,6])
          results1ac[,16]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,7])
          results1ac[,7]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,2])
          results1ac[,11]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,4])
          results1ac[,13]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,6])
          results1ac[,17]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,7])
          results1ac[,9]<-results1ac[,7]/results1ac[,8]
          for(o in 6:10){
            results1ac[,14]<-rowSums(data.frame(results1ac[,14],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,6])),na.rm=T)
            results1ac[,18]<-rowSums(data.frame(results1ac[,18],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,7])),na.rm=T)   
          }
       results1ac[,15]<-results1ac[,13]/results1ac[,14]
       results1ac[,19]<-results1ac[,17]/results1ac[,18]
        }
        if(length(results1)>0){
          results1<-rbind(results1,results1ac)  
        }
        if(length(results1)==0){
          results1<-results1ac
        }
    }
    
    #C>G
    input1ac<-input1[input1[,3]=="C"&input1[,4]=="G",]
    
    if(length(input1ac[,1])>0){
      results1ac<-data.frame(input1ac[,c(1:5)],Nr_Ref=NA,Nr_Alt=NA,DP=NA,VAF=NA,BQ_REF=NA,BQ_ALT=NA,
                             Nr_Ref_fwd=NA,Nr_Alt_fwd=NA,DP_fwd=NA,VAF_fwd=NA,
                             Nr_Ref_rev=NA,Nr_Alt_rev=NA,DP_rev=NA,VAF_rev=NA)
      
      write.table(cbind(input1ac[,1],input1ac[,2],input1ac[,2]),paste(dir1,"temp",start_pos,"_",end_pos,".txt",sep=""),
                  quote=F,row.names = F,sep="\t")
      
      message("SNPs: ",length(input1ac[,1]))
        system(paste("tail -n +2 ",dir1,"temp",start_pos,"_",end_pos,".txt > ",dir1,"temp",start_pos,"_",end_pos,"_2.txt",sep=""),intern = F)
        result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                 dir_bam,"/",samples[n,1],".bam -l ",dir1,"temp",start_pos,"_",end_pos,"_2.txt -w 0",sep=""),intern=T)
        result_raw2<-str_split_fixed(result_raw,pattern = "\t",n=Inf)
        message("Evaluation")
        if(length(result_raw2)>0){
          line<-match(results1ac[,2],result_raw2[,2])
          #line<-line[!is.na(line)]
          results1ac[,8]<-as.numeric(result_raw2[line,4])
          results1ac[,6]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,2])
          results1ac[,10]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,4])
          results1ac[,12]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,6])
          results1ac[,16]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,7])
          results1ac[,7]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,2])
          results1ac[,11]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,4])
          results1ac[,13]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,6])
          results1ac[,17]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,7])
          results1ac[,9]<-results1ac[,7]/results1ac[,8]
          for(o in 6:10){
            results1ac[,14]<-rowSums(data.frame(results1ac[,14],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,6])),na.rm=T)
            results1ac[,18]<-rowSums(data.frame(results1ac[,18],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,7])),na.rm=T)   
          }
       results1ac[,15]<-results1ac[,13]/results1ac[,14]
       results1ac[,19]<-results1ac[,17]/results1ac[,18]
        }
        if(length(results1)>0){
          results1<-rbind(results1,results1ac)  
        }
        if(length(results1)==0){
          results1<-results1ac
        }
    }
    
    #C>T
    input1ac<-input1[input1[,3]=="C"&input1[,4]=="T",]
    
    if(length(input1ac[,1])>0){
      results1ac<-data.frame(input1ac[,c(1:5)],Nr_Ref=NA,Nr_Alt=NA,DP=NA,VAF=NA,BQ_REF=NA,BQ_ALT=NA,
                             Nr_Ref_fwd=NA,Nr_Alt_fwd=NA,DP_fwd=NA,VAF_fwd=NA,
                             Nr_Ref_rev=NA,Nr_Alt_rev=NA,DP_rev=NA,VAF_rev=NA)
      
      write.table(cbind(input1ac[,1],input1ac[,2],input1ac[,2]),paste(dir1,"temp",start_pos,"_",end_pos,".txt",sep=""),
                  quote=F,row.names = F,sep="\t")
      
      message("SNPs: ",length(input1ac[,1]))
        system(paste("tail -n +2 ",dir1,"temp",start_pos,"_",end_pos,".txt > ",dir1,"temp",start_pos,"_",end_pos,"_2.txt",sep=""),intern = F)
        result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                 dir_bam,"/",samples[n,1],".bam -l ",dir1,"temp",start_pos,"_",end_pos,"_2.txt -w 0",sep=""),intern=T)
        result_raw2<-str_split_fixed(result_raw,pattern = "\t",n=Inf)
        message("Evaluation")
        if(length(result_raw2)>0){
        line<-match(results1ac[,2],result_raw2[,2])
        #line<-line[!is.na(line)]
       results1ac[,8]<-as.numeric(result_raw2[line,4])
       results1ac[,6]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,2])
       results1ac[,10]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,4])
       results1ac[,12]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,6])
       results1ac[,16]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,7])
       results1ac[,7]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,2])
       results1ac[,11]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,4])
       results1ac[,13]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,6])
       results1ac[,17]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,7])
       results1ac[,9]<-results1ac[,7]/results1ac[,8]
        for(o in 6:10){
           results1ac[,14]<-rowSums(data.frame(results1ac[,14],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,6])),na.rm=T)
           results1ac[,18]<-rowSums(data.frame(results1ac[,18],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,7])),na.rm=T)   
        }
       results1ac[,15]<-results1ac[,13]/results1ac[,14]
       results1ac[,19]<-results1ac[,17]/results1ac[,18]
        }
        if(length(results1)>0){
          results1<-rbind(results1,results1ac)  
        }
        if(length(results1)==0){
          results1<-results1ac
        }
    }
    
    #G>A
    input1ac<-input1[input1[,3]=="G"&input1[,4]=="A",]
    
    if(length(input1ac[,1])>0){
      results1ac<-data.frame(input1ac[,c(1:5)],Nr_Ref=NA,Nr_Alt=NA,DP=NA,VAF=NA,BQ_REF=NA,BQ_ALT=NA,
                             Nr_Ref_fwd=NA,Nr_Alt_fwd=NA,DP_fwd=NA,VAF_fwd=NA,
                             Nr_Ref_rev=NA,Nr_Alt_rev=NA,DP_rev=NA,VAF_rev=NA)
      
      write.table(cbind(input1ac[,1],input1ac[,2],input1ac[,2]),paste(dir1,"temp",start_pos,"_",end_pos,".txt",sep=""),
                  quote=F,row.names = F,sep="\t")
      
      message("SNPs: ",length(input1ac[,1]))
        system(paste("tail -n +2 ",dir1,"temp",start_pos,"_",end_pos,".txt > ",dir1,"temp",start_pos,"_",end_pos,"_2.txt",sep=""),intern = F)
        result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                 dir_bam,"/",samples[n,1],".bam -l ",dir1,"temp",start_pos,"_",end_pos,"_2.txt -w 0",sep=""),intern=T)
        result_raw2<-str_split_fixed(result_raw,pattern = "\t",n=Inf)
        message("Evaluation")
        if(length(result_raw2)>0){
          line<-match(results1ac[,2],result_raw2[,2])
          #line<-line[!is.na(line)]
          results1ac[,8]<-as.numeric(result_raw2[line,4])
          results1ac[,6]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,2])
          results1ac[,10]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,4])
          results1ac[,12]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,6])
          results1ac[,16]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,7])
          results1ac[,7]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,2])
          results1ac[,11]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,4])
          results1ac[,13]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,6])
          results1ac[,17]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,7])
          results1ac[,9]<-results1ac[,7]/results1ac[,8]
          for(o in 6:10){
            results1ac[,14]<-rowSums(data.frame(results1ac[,14],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,6])),na.rm=T)
            results1ac[,18]<-rowSums(data.frame(results1ac[,18],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,7])),na.rm=T)   
          }
       results1ac[,15]<-results1ac[,13]/results1ac[,14]
       results1ac[,19]<-results1ac[,17]/results1ac[,18]
        }
        if(length(results1)>0){
          results1<-rbind(results1,results1ac)  
        }
        if(length(results1)==0){
          results1<-results1ac
        }
    }
    
    #G>C
    input1ac<-input1[input1[,3]=="G"&input1[,4]=="C",]
    
    if(length(input1ac[,1])>0){
      results1ac<-data.frame(input1ac[,c(1:5)],Nr_Ref=NA,Nr_Alt=NA,DP=NA,VAF=NA,BQ_REF=NA,BQ_ALT=NA,
                             Nr_Ref_fwd=NA,Nr_Alt_fwd=NA,DP_fwd=NA,VAF_fwd=NA,
                             Nr_Ref_rev=NA,Nr_Alt_rev=NA,DP_rev=NA,VAF_rev=NA)
      
      write.table(cbind(input1ac[,1],input1ac[,2],input1ac[,2]),paste(dir1,"temp",start_pos,"_",end_pos,".txt",sep=""),
                  quote=F,row.names = F,sep="\t")
      
      message("SNPs: ",length(input1ac[,1]))
        system(paste("tail -n +2 ",dir1,"temp",start_pos,"_",end_pos,".txt > ",dir1,"temp",start_pos,"_",end_pos,"_2.txt",sep=""),intern = F)
        result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                 dir_bam,"/",samples[n,1],".bam -l ",dir1,"temp",start_pos,"_",end_pos,"_2.txt -w 0",sep=""),intern=T)
        result_raw2<-str_split_fixed(result_raw,pattern = "\t",n=Inf)
        message("Evaluation")
        if(length(result_raw2)>0){
          line<-match(results1ac[,2],result_raw2[,2])
          #line<-line[!is.na(line)]
          results1ac[,8]<-as.numeric(result_raw2[line,4])
          results1ac[,6]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,2])
          results1ac[,10]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,4])
          results1ac[,12]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,6])
          results1ac[,16]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,7])
          results1ac[,7]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,2])
          results1ac[,11]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,4])
          results1ac[,13]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,6])
          results1ac[,17]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,7])
          results1ac[,9]<-results1ac[,7]/results1ac[,8]
          for(o in 6:10){
            results1ac[,14]<-rowSums(data.frame(results1ac[,14],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,6])),na.rm=T)
            results1ac[,18]<-rowSums(data.frame(results1ac[,18],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,7])),na.rm=T)   
          }
       results1ac[,15]<-results1ac[,13]/results1ac[,14]
       results1ac[,19]<-results1ac[,17]/results1ac[,18]
        }
        if(length(results1)>0){
          results1<-rbind(results1,results1ac)  
        }
        if(length(results1)==0){
          results1<-results1ac
        }
    }

    #G>T
    input1ac<-input1[input1[,3]=="G"&input1[,4]=="T",]
    
    if(length(input1ac[,1])>0){
      results1ac<-data.frame(input1ac[,c(1:5)],Nr_Ref=NA,Nr_Alt=NA,DP=NA,VAF=NA,BQ_REF=NA,BQ_ALT=NA,
                             Nr_Ref_fwd=NA,Nr_Alt_fwd=NA,DP_fwd=NA,VAF_fwd=NA,
                             Nr_Ref_rev=NA,Nr_Alt_rev=NA,DP_rev=NA,VAF_rev=NA)
      
      write.table(cbind(input1ac[,1],input1ac[,2],input1ac[,2]),paste(dir1,"temp",start_pos,"_",end_pos,".txt",sep=""),
                  quote=F,row.names = F,sep="\t")
      
      message("SNPs: ",length(input1ac[,1]))
        system(paste("tail -n +2 ",dir1,"temp",start_pos,"_",end_pos,".txt > ",dir1,"temp",start_pos,"_",end_pos,"_2.txt",sep=""),intern = F)
        result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                 dir_bam,"/",samples[n,1],".bam -l ",dir1,"temp",start_pos,"_",end_pos,"_2.txt -w 0",sep=""),intern=T)
        result_raw2<-str_split_fixed(result_raw,pattern = "\t",n=Inf)
        message("Evaluation")
        if(length(result_raw2)>0){
          line<-match(results1ac[,2],result_raw2[,2])
          #line<-line[!is.na(line)]
          results1ac[,8]<-as.numeric(result_raw2[line,4])
          results1ac[,6]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,2])
          results1ac[,10]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,4])
          results1ac[,12]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,6])
          results1ac[,16]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,7])
          results1ac[,7]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,2])
          results1ac[,11]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,4])
          results1ac[,13]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,6])
          results1ac[,17]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,7])
          results1ac[,9]<-results1ac[,7]/results1ac[,8]
          for(o in 6:10){
            results1ac[,14]<-rowSums(data.frame(results1ac[,14],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,6])),na.rm=T)
            results1ac[,18]<-rowSums(data.frame(results1ac[,18],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,7])),na.rm=T)   
          }
       results1ac[,15]<-results1ac[,13]/results1ac[,14]
       results1ac[,19]<-results1ac[,17]/results1ac[,18]
        }
        if(length(results1)>0){
          results1<-rbind(results1,results1ac)  
        }
        if(length(results1)==0){
          results1<-results1ac
        }
    }
    
    #T>A
    input1ac<-input1[input1[,3]=="T"&input1[,4]=="A",]
    
    if(length(input1ac[,1])>0){
      results1ac<-data.frame(input1ac[,c(1:5)],Nr_Ref=NA,Nr_Alt=NA,DP=NA,VAF=NA,BQ_REF=NA,BQ_ALT=NA,
                             Nr_Ref_fwd=NA,Nr_Alt_fwd=NA,DP_fwd=NA,VAF_fwd=NA,
                             Nr_Ref_rev=NA,Nr_Alt_rev=NA,DP_rev=NA,VAF_rev=NA)
      
      write.table(cbind(input1ac[,1],input1ac[,2],input1ac[,2]),paste(dir1,"temp",start_pos,"_",end_pos,".txt",sep=""),
                  quote=F,row.names = F,sep="\t")
      
      message("SNPs: ",length(input1ac[,1]))
        system(paste("tail -n +2 ",dir1,"temp",start_pos,"_",end_pos,".txt > ",dir1,"temp",start_pos,"_",end_pos,"_2.txt",sep=""),intern = F)
        result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                 dir_bam,"/",samples[n,1],".bam -l ",dir1,"temp",start_pos,"_",end_pos,"_2.txt -w 0",sep=""),intern=T)
        result_raw2<-str_split_fixed(result_raw,pattern = "\t",n=Inf)
        message("Evaluation")
        if(length(result_raw2)>0){
          line<-match(results1ac[,2],result_raw2[,2])
          #line<-line[!is.na(line)]
          results1ac[,8]<-as.numeric(result_raw2[line,4])
          results1ac[,6]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,2])
          results1ac[,10]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,4])
          results1ac[,12]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,6])
          results1ac[,16]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,7])
          results1ac[,7]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,2])
          results1ac[,11]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,4])
          results1ac[,13]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,6])
          results1ac[,17]<-as.numeric(str_split_fixed(result_raw2[line,6],":",n=Inf)[,7])
          results1ac[,9]<-results1ac[,7]/results1ac[,8]
          for(o in 6:10){
            results1ac[,14]<-rowSums(data.frame(results1ac[,14],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,6])),na.rm=T)
            results1ac[,18]<-rowSums(data.frame(results1ac[,18],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,7])),na.rm=T)   
          }
       results1ac[,15]<-results1ac[,13]/results1ac[,14]
       results1ac[,19]<-results1ac[,17]/results1ac[,18]
        }
        if(length(results1)>0){
          results1<-rbind(results1,results1ac)  
        }
        if(length(results1)==0){
          results1<-results1ac
        } 
    }
    
    #T>C
    input1ac<-input1[input1[,3]=="T"&input1[,4]=="C",]
    
    if(length(input1ac[,1])>0){
      results1ac<-data.frame(input1ac[,c(1:5)],Nr_Ref=NA,Nr_Alt=NA,DP=NA,VAF=NA,BQ_REF=NA,BQ_ALT=NA,
                             Nr_Ref_fwd=NA,Nr_Alt_fwd=NA,DP_fwd=NA,VAF_fwd=NA,
                             Nr_Ref_rev=NA,Nr_Alt_rev=NA,DP_rev=NA,VAF_rev=NA)
      
      write.table(cbind(input1ac[,1],input1ac[,2],input1ac[,2]),paste(dir1,"temp",start_pos,"_",end_pos,".txt",sep=""),
                  quote=F,row.names = F,sep="\t")
      
      message("SNPs: ",length(input1ac[,1]))
        system(paste("tail -n +2 ",dir1,"temp",start_pos,"_",end_pos,".txt > ",dir1,"temp",start_pos,"_",end_pos,"_2.txt",sep=""),intern = F)
        result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                 dir_bam,"/",samples[n,1],".bam -l ",dir1,"temp",start_pos,"_",end_pos,"_2.txt -w 0",sep=""),intern=T)
        result_raw2<-str_split_fixed(result_raw,pattern = "\t",n=Inf)
        message("Evaluation")
        if(length(result_raw2)>0){
          line<-match(results1ac[,2],result_raw2[,2])
          #line<-line[!is.na(line)]
          results1ac[,8]<-as.numeric(result_raw2[line,4])
          results1ac[,6]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,2])
          results1ac[,10]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,4])
          results1ac[,12]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,6])
          results1ac[,16]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,7])
          results1ac[,7]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,2])
          results1ac[,11]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,4])
          results1ac[,13]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,6])
          results1ac[,17]<-as.numeric(str_split_fixed(result_raw2[line,7],":",n=Inf)[,7])
          results1ac[,9]<-results1ac[,7]/results1ac[,8]
          for(o in 6:10){
            results1ac[,14]<-rowSums(data.frame(results1ac[,14],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,6])),na.rm=T)
            results1ac[,18]<-rowSums(data.frame(results1ac[,18],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,7])),na.rm=T)   
          }
       results1ac[,15]<-results1ac[,13]/results1ac[,14]
       results1ac[,19]<-results1ac[,17]/results1ac[,18]
        }
        if(length(results1)>0){
          results1<-rbind(results1,results1ac)  
        }
        if(length(results1)==0){
          results1<-results1ac
        }
    }
    
    #T>G
    input1ac<-input1[input1[,3]=="T"&input1[,4]=="G",]
    
    if(length(input1ac[,1])>0){
      results1ac<-data.frame(input1ac[,c(1:5)],Nr_Ref=NA,Nr_Alt=NA,DP=NA,VAF=NA,BQ_REF=NA,BQ_ALT=NA,
                             Nr_Ref_fwd=NA,Nr_Alt_fwd=NA,DP_fwd=NA,VAF_fwd=NA,
                             Nr_Ref_rev=NA,Nr_Alt_rev=NA,DP_rev=NA,VAF_rev=NA)
      
      write.table(cbind(input1ac[,1],input1ac[,2],input1ac[,2]),paste(dir1,"temp",start_pos,"_",end_pos,".txt",sep=""),
                  quote=F,row.names = F,sep="\t")
      
      message("SNPs: ",length(input1ac[,1]))
        system(paste("tail -n +2 ",dir1,"temp",start_pos,"_",end_pos,".txt > ",dir1,"temp",start_pos,"_",end_pos,"_2.txt",sep=""),intern = F)
        result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                 dir_bam,"/",samples[n,1],".bam -l ",dir1,"temp",start_pos,"_",end_pos,"_2.txt -w 0",sep=""),intern=T)
        result_raw2<-str_split_fixed(result_raw,pattern = "\t",n=Inf)
        message("Evaluation")
        if(length(result_raw2)>0){
          line<-match(results1ac[,2],result_raw2[,2])
          #line<-line[!is.na(line)]
          results1ac[,8]<-as.numeric(result_raw2[line,4])
          results1ac[,6]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,2])
          results1ac[,10]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,4])
          results1ac[,12]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,6])
          results1ac[,16]<-as.numeric(str_split_fixed(result_raw2[line,9],":",n=Inf)[,7])
          results1ac[,7]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,2])
          results1ac[,11]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,4])
          results1ac[,13]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,6])
          results1ac[,17]<-as.numeric(str_split_fixed(result_raw2[line,8],":",n=Inf)[,7])
          results1ac[,9]<-results1ac[,7]/results1ac[,8]
          for(o in 6:10){
            results1ac[,14]<-rowSums(data.frame(results1ac[,14],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,6])),na.rm=T)
            results1ac[,18]<-rowSums(data.frame(results1ac[,18],as.numeric(str_split_fixed(result_raw2[line,o],":",n=Inf)[,7])),na.rm=T)   
          }
       results1ac[,15]<-results1ac[,13]/results1ac[,14]
       results1ac[,19]<-results1ac[,17]/results1ac[,18]
        }
        if(length(results1)>0){
          results1<-rbind(results1,results1ac)  
        }
        if(length(results1)==0){
          results1<-results1ac
        }
    }
    
    write.table(results1,paste(dir1,"temp_SNVs.txt",sep=""),
                quote=F,row.names = F,sep="\t")
    
    #Rest1
    if(length(input2[,1])>0){
        input2.1a<-input2[nchar(input2[,4])==1,]
        if(length(input2.1a[,1])>0){
            message("Dels:",length(input2.1a[,1]))
            results2.1a<-data.frame(input2.1a[,c(1:5)],Nr_Ref=NA,Nr_Alt=NA,DP=NA,VAF=NA,BQ_REF=NA,BQ_ALT=NA,
                                   Nr_Ref_fwd=NA,Nr_Alt_fwd=NA,DP_fwd=NA,VAF_fwd=NA,
                                   Nr_Ref_rev=NA,Nr_Alt_rev=NA,DP_rev=NA,VAF_rev=NA)
            write.table(cbind(input2.1a[,1],input2.1a[,2]+1,input2.1a[,2]+1),paste(dir1,"temp",start_pos,"_",end_pos,".txt",sep=""),
                        quote=F,row.names = F,sep="\t")
            system(paste("tail -n +2 ",dir1,"temp",start_pos,"_",end_pos,".txt > ",dir1,"temp",start_pos,"_",end_pos,"_2.txt",sep=""),intern = F)
            result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                     dir_bam,"/",samples[n,1],".bam -l ",dir1,"temp",start_pos,"_",end_pos,"_2.txt -w 0",sep=""),intern=T)
            result_raw2<-str_split_fixed(result_raw,pattern = "\t",n=Inf)
            counter<-1
            message("Evaluation")
            for(i in 1:length(result_raw2[,1])){
                if((i%%1000)==0){
                    message("i=",i," out of ",length(result_raw2[,1]))
                }
                if(counter<=length(input2.1a[,1])&&
                   (result_raw2[i,1]!=input2.1a[counter,1]||result_raw2[i,2]!=(input2.1a[counter,2]+1))){
                    counter<-counter+1
                }
                if(counter<=length(input2.1a[,1])&&
                   (result_raw2[i,1]==input2.1a[counter,1]&&result_raw2[i,2]==(input2.1a[counter,2]+1))){
                    results2.1a[counter,8]<-as.numeric(result_raw2[i,4])
                    if(substr(results2.1a[counter,3],2,2)=="A"){
                        results2.1a[counter,6]<-as.numeric(strsplit(result_raw2[i,6],split=":")[[1]][2])
                        results2.1a[counter,10]<-as.numeric(strsplit(result_raw2[i,6],split=":")[[1]][4])
                        results2.1a[counter,12]<-as.numeric(strsplit(result_raw2[i,6],split=":")[[1]][6])
                        results2.1a[counter,16]<-as.numeric(strsplit(result_raw2[i,6],split=":")[[1]][7])
                    }
                    if(substr(results2.1a[counter,3],2,2)=="C"){
                        results2.1a[counter,6]<-as.numeric(strsplit(result_raw2[i,7],split=":")[[1]][2])
                        results2.1a[counter,10]<-as.numeric(strsplit(result_raw2[i,7],split=":")[[1]][4])
                        results2.1a[counter,12]<-as.numeric(strsplit(result_raw2[i,7],split=":")[[1]][6])
                        results2.1a[counter,16]<-as.numeric(strsplit(result_raw2[i,7],split=":")[[1]][7])
                    }
                    if(substr(results2.1a[counter,3],2,2)=="G"){
                        results2.1a[counter,6]<-as.numeric(strsplit(result_raw2[i,8],split=":")[[1]][2])
                        results2.1a[counter,10]<-as.numeric(strsplit(result_raw2[i,8],split=":")[[1]][4])
                        results2.1a[counter,12]<-as.numeric(strsplit(result_raw2[i,8],split=":")[[1]][6])
                        results2.1a[counter,16]<-as.numeric(strsplit(result_raw2[i,8],split=":")[[1]][7])
                    }
                    if(substr(results2.1a[counter,3],2,2)=="T"){
                        results2.1a[counter,6]<-as.numeric(strsplit(result_raw2[i,9],split=":")[[1]][2])
                        results2.1a[counter,10]<-as.numeric(strsplit(result_raw2[i,9],split=":")[[1]][4])
                        results2.1a[counter,12]<-as.numeric(strsplit(result_raw2[i,9],split=":")[[1]][6])
                        results2.1a[counter,16]<-as.numeric(strsplit(result_raw2[i,9],split=":")[[1]][7])
                    }
                    for(j in 10:length(result_raw2[i,])){
                        if(result_raw2[i,j]!=""&&
                           paste("-",substr(results2.1a[counter,3],2,nchar(results2.1a[counter,3])),sep="")==strsplit(result_raw2[i,j],split=":")[[1]][1]){
                            results2.1a[counter,7]<-as.numeric(strsplit(result_raw2[i,j],split=":")[[1]][2])
                            results2.1a[counter,11]<-NA
                            results2.1a[counter,13]<-as.numeric(strsplit(result_raw2[i,j],split=":")[[1]][6])
                            results2.1a[counter,17]<-as.numeric(strsplit(result_raw2[i,j],split=":")[[1]][7])
                        }
                    }
                    results2.1a[counter,9]<-results2.1a[counter,7]/results2.1a[counter,8]
                    for(o in 6:length(result_raw2[i,])){
                        if(result_raw2[i,o]!=""&&
                           vcountPattern("+",strsplit(result_raw2[i,o],split=":")[[1]][1])==0){
                            results2.1a[counter,14]<-sum(results2.1a[counter,14],as.numeric(strsplit(result_raw2[i,o],split=":")[[1]][6]),na.rm=T)
                            results2.1a[counter,18]<-sum(results2.1a[counter,18],as.numeric(strsplit(result_raw2[i,o],split=":")[[1]][7]),na.rm=T)   
                        }
                    }
                    results2.1a[counter,15]<-results2.1a[counter,13]/results2.1a[counter,14]
                    results2.1a[counter,19]<-results2.1a[counter,17]/results2.1a[counter,18] 
                    counter<-counter+1
                }
            }
            if(length(results1)>0){
              results1<-rbind(results1,results2.1a)  
            }
            if(length(results1)==0){
              results1<-results2.1a
            }
        }

        input2.1b<-input2[nchar(input2[,3])==1&input2[,3]==substr(input2[,4],1,1),]
        if(length(input2.1b[,1])>0){
            message("Ins:",length(input2.1b[,1]))
            results2.1b<-data.frame(input2.1b[,c(1:5)],Nr_Ref=NA,Nr_Alt=NA,DP=NA,VAF=NA,BQ_REF=NA,BQ_ALT=NA,
                                    Nr_Ref_fwd=NA,Nr_Alt_fwd=NA,DP_fwd=NA,VAF_fwd=NA,
                                    Nr_Ref_rev=NA,Nr_Alt_rev=NA,DP_rev=NA,VAF_rev=NA)
            write.table(cbind(input2.1b[,1],input2.1b[,2],input2.1b[,2]),paste(dir1,"temp",start_pos,"_",end_pos,".txt",sep=""),
                        quote=F,row.names = F,sep="\t")
            system(paste("tail -n +2 ",dir1,"temp",start_pos,"_",end_pos,".txt > ",dir1,"temp",start_pos,"_",end_pos,"_2.txt",sep=""),intern = F)
            result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                     dir_bam,"/",samples[n,1],".bam -l ",dir1,"temp",start_pos,"_",end_pos,"_2.txt -w 0",sep=""),intern=T)
            result_raw2<-str_split_fixed(result_raw,pattern = "\t",n=Inf)
            counter<-1
            message("Evaluation")
            for(i in 1:length(result_raw2[,1])){
                if((i%%1000)==0){
                    message("i=",i," out of ",length(result_raw2[,1]))
                }
                if(counter<=length(input2.1b[,1])&&
                   (result_raw2[i,1]!=input2.1b[counter,1]||result_raw2[i,2]!=input2.1b[counter,2])){
                    counter<-counter+1
                }
                if(counter<=length(input2.1b[,1])&&
                   (result_raw2[i,1]==input2.1b[counter,1]&&result_raw2[i,2]==input2.1b[counter,2])){
                    results2.1b[counter,8]<-as.numeric(result_raw2[i,4])
                    if(results2.1b[counter,3]=="A"){
                        results2.1b[counter,6]<-as.numeric(strsplit(result_raw2[i,6],split=":")[[1]][2])
                        results2.1b[counter,10]<-as.numeric(strsplit(result_raw2[i,6],split=":")[[1]][4])
                        results2.1b[counter,12]<-as.numeric(strsplit(result_raw2[i,6],split=":")[[1]][6])
                        results2.1b[counter,16]<-as.numeric(strsplit(result_raw2[i,6],split=":")[[1]][7])
                    }
                    if(results2.1b[counter,3]=="C"){
                        results2.1b[counter,6]<-as.numeric(strsplit(result_raw2[i,7],split=":")[[1]][2])
                        results2.1b[counter,10]<-as.numeric(strsplit(result_raw2[i,7],split=":")[[1]][4])
                        results2.1b[counter,12]<-as.numeric(strsplit(result_raw2[i,7],split=":")[[1]][6])
                        results2.1b[counter,16]<-as.numeric(strsplit(result_raw2[i,7],split=":")[[1]][7])
                    }
                    if(results2.1b[counter,3]=="G"){
                        results2.1b[counter,6]<-as.numeric(strsplit(result_raw2[i,8],split=":")[[1]][2])
                        results2.1b[counter,10]<-as.numeric(strsplit(result_raw2[i,8],split=":")[[1]][4])
                        results2.1b[counter,12]<-as.numeric(strsplit(result_raw2[i,8],split=":")[[1]][6])
                        results2.1b[counter,16]<-as.numeric(strsplit(result_raw2[i,8],split=":")[[1]][7])
                    }
                    if(results2.1b[counter,3]=="T"){
                        results2.1b[counter,6]<-as.numeric(strsplit(result_raw2[i,9],split=":")[[1]][2])
                        results2.1b[counter,10]<-as.numeric(strsplit(result_raw2[i,9],split=":")[[1]][4])
                        results2.1b[counter,12]<-as.numeric(strsplit(result_raw2[i,9],split=":")[[1]][6])
                        results2.1b[counter,16]<-as.numeric(strsplit(result_raw2[i,9],split=":")[[1]][7])
                    }
                    for(j in 10:length(result_raw2[i,])){
                        if(result_raw2[i,j]!=""&&
                            paste("+",substr(results2.1b[counter,4],2,nchar(results2.1b[counter,4])),sep="")==strsplit(result_raw2[i,j],split=":")[[1]][1]){
                            results2.1b[counter,7]<-as.numeric(strsplit(result_raw2[i,j],split=":")[[1]][2])
                            results2.1b[counter,11]<-NA
                            results2.1b[counter,13]<-as.numeric(strsplit(result_raw2[i,j],split=":")[[1]][6])
                            results2.1b[counter,17]<-as.numeric(strsplit(result_raw2[i,j],split=":")[[1]][7])
                        }
                    }
                    results2.1b[counter,9]<-results2.1b[counter,7]/results2.1b[counter,8]
                    for(o in 6:length(result_raw2[i,])){
                        if(result_raw2[i,o]!=""&&
                           vcountPattern("+",strsplit(result_raw2[i,o],split=":")[[1]][1])==0){
                            results2.1b[counter,14]<-sum(results2.1b[counter,14],as.numeric(strsplit(result_raw2[i,o],split=":")[[1]][6]),na.rm=T)
                            results2.1b[counter,18]<-sum(results2.1b[counter,18],as.numeric(strsplit(result_raw2[i,o],split=":")[[1]][7]),na.rm=T)   
                        }
                    }
                    results2.1b[counter,15]<-results2.1b[counter,13]/results2.1b[counter,14]
                    results2.1b[counter,19]<-results2.1b[counter,17]/results2.1b[counter,18]
                }
                counter<-counter+1
            }
            if(length(results1)>0){
              results1<-rbind(results1,results2.1b)  
            }
            if(length(results1)==0){
              results1<-results2.1b
            }
        }
        
        #Rest2
        message("Indels")
        input2.2<-input2[(nchar(input2[,3])!=1&nchar(input2[,4])!=1)|(nchar(input2[,3])==1&input2[,3]!=substr(input2[,4],1,1)),]
        if(length(input2.2[,1])>0){
            results2.2<-data.frame(input2.2[,c(1:5)],Nr_Ref=NA,Nr_Alt=NA,DP=NA,VAF=NA,BQ_REF=NA,BQ_ALT=NA,
                                   Nr_Ref_fwd=NA,Nr_Alt_fwd=NA,DP_fwd=NA,VAF_fwd=NA,
                                   Nr_Ref_rev=NA,Nr_Alt_rev=NA,DP_rev=NA,VAF_rev=NA)
            
            for(i in 1:length(input2.2[,1])){
                message("Mutation ",i," of ",length(input2.2[,1]))
                #INDEL
                if(nchar(input2.2[i,3])==nchar(input2.2[i,4])){
                    message(input2.2[i,])
                }
                if(nchar(input2.2[i,3])!=nchar(input2.2[i,4])){
                    #2 ins am anfang und ende
                    if(vcountPattern(input2.2[i,3],input2.2[i,4])>0){
                        result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                                 dir_bam,"/",input2.2[i,5],".bam ",input2.2[i,1],":",(input2.2[i,2]-1),"-",(input2.2[i,2]+nchar(input2.2[i,3])-1)," -w 0",sep=""),intern=T)
                        result_raw2<-strsplit(result_raw,split = "\t")
                        if(length(result_raw2)>0&&length(result_raw2)==(nchar(input2.2[i,3])+1)){
                            dp<-rep(0,nchar(input2.2[i,3]))
                            ref<-rep(0,nchar(input2.2[i,3]))
                            alt<-rep(0,nchar(input2.2[i,3])+2)
                            bq_ref<-rep(0,nchar(input2.2[i,3]))
                            bq_alt<-NA
                            
                            for(j in 1:nchar(input2.2[i,3])){
                                dp[j]<-as.numeric(result_raw2[[j]][4])
                                if(substr(input2.2[i,3],j,j)=="A"){
                                    ref[j]<-as.numeric(strsplit(result_raw2[[j+1]][6],split=":")[[1]][2])
                                    bq_ref[j]<-as.numeric(strsplit(result_raw2[[j+1]][6],split=":")[[1]][4])
                                }
                                if(substr(input2.2[i,3],j,j)=="C"){
                                    ref[j]<-as.numeric(strsplit(result_raw2[[j+1]][7],split=":")[[1]][2])
                                    bq_ref[j]<-as.numeric(strsplit(result_raw2[[j+1]][7],split=":")[[1]][4])
                                }
                                if(substr(input2.2[i,3],j,j)=="G"){
                                    ref[j]<-as.numeric(strsplit(result_raw2[[j+1]][8],split=":")[[1]][2])
                                    bq_ref[j]<-as.numeric(strsplit(result_raw2[[j+1]][8],split=":")[[1]][4])
                                }
                                if(substr(input2.2[i,3],j,j)=="T"){
                                    ref[j]<-as.numeric(strsplit(result_raw2[[j+1]][9],split=":")[[1]][2])
                                    bq_ref[j]<-as.numeric(strsplit(result_raw2[[j+1]][9],split=":")[[1]][4])
                                }
                            }
                            position<-matchPattern(input2.2[i,3],input2.2[i,4])
                            indel1<-substr(input2.2[i,4],1,start(position[1])-1)
                            indel2<-substr(input2.2[i,4],end(position[1])+1,nchar(input2.2[i,4]))
                            for(j in 10:length(result_raw2[[1]])){
                                if(paste("+",indel1,sep="")==strsplit(result_raw2[[1]][j],split=":")[[1]][1]){
                                    alt[1]<-as.numeric(strsplit(result_raw2[[1]][j],split=":")[[1]][2])
                                }
                            }
                            #sonderfall homopolymer
                            if(alt[1]==0&&result_raw2[[1]][3]==indel1){
                                move<-2
                                result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                                         dir_bam,"/",input2.2[i,5],".bam ",input2.2[i,1],":",(input2.2[i,2]-move),"-",(input2.2[i,2]+nchar(input2.2[i,3])-1)," -w 0",sep=""),intern=T)
                                result_raw_hp<-strsplit(result_raw,split = "\t")
                                while(alt[1]==0&&result_raw_hp[[1]][3]==indel1){
                                    for(j in 10:length(result_raw_hp[[1]])){
                                        if(paste("+",indel1,sep="")==strsplit(result_raw_hp[[1]][j],split=":")[[1]][1]){
                                            alt[1]<-as.numeric(strsplit(result_raw_hp[[1]][j],split=":")[[1]][2])
                                        }
                                    }
                                    move<-move+1
                                    result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                                             dir_bam,"/",input2.2[i,5],".bam ",input2.2[i,1],":",(input2.2[i,2]-move),"-",(input2.2[i,2]+nchar(input2.2[i,3])-1)," -w 0",sep=""),intern=T)
                                    result_raw_hp<-strsplit(result_raw,split = "\t")
                                }
                                for(j in 10:length(result_raw_hp[[1]])){
                                    if(paste("+",indel1,sep="")==strsplit(result_raw_hp[[1]][j],split=":")[[1]][1]){
                                        alt[1]<-as.numeric(strsplit(result_raw_hp[[1]][j],split=":")[[1]][2])
                                    }
                                }
                            }
                            
                            for(j in 1:length(ref)){
                                alt[j+1]<-ref[j]
                            }
                            for(j in 10:length(result_raw2[[length(result_raw2)]])){
                                if(paste("+",indel2,sep="")==strsplit(result_raw2[[length(result_raw2)]][j],split=":")[[1]][1]){
                                    alt[length(alt)]<-as.numeric(strsplit(result_raw2[[length(result_raw2)]][j],split=":")[[1]][2])
                                }
                            }
                            results2.2[i,6]<-min(ref)
                            results2.2[i,7]<-min(alt)
                            results2.2[i,8]<-min(dp)
                            results2.2[i,9]<-results2.2[i,7]/results2.2[i,8]
                            results2.2[i,10]<-mean(bq_ref)
                            results2.2[i,11]<-bq_alt
                        }
                    }
                    if(vcountPattern(input2.2[i,3],input2.2[i,4])==0){
                        result_raw<-system(paste("bam-readcount"," -f ",reference," ",
                                                 dir_bam,"/",input2.2[i,5],".bam ",input2.2[i,1],":",input2.2[i,2],"-",(input2.2[i,2]+nchar(input2.2[i,3])-1)," -w 0",sep=""),intern=T)
                        result_raw2<-strsplit(result_raw,split = "\t")
                        if(length(result_raw2)>0&&length(result_raw2)==nchar(input2.2[i,3])){
                            dp<-rep(0,nchar(input2.2[i,3]))
                            ref<-rep(0,nchar(input2.2[i,3]))
                            alt<-rep(0,nchar(input2.2[i,3]))
                            bq_ref<-rep(0,nchar(input2.2[i,3]))
                            bq_alt<-"ComplexIndel"
                            
                            for(j in 1:nchar(input2.2[i,3])){
                                dp[j]<-as.numeric(result_raw2[[j]][4])
                                if(substr(input2.2[i,3],j,j)=="A"){
                                    ref[j]<-as.numeric(strsplit(result_raw2[[j]][6],split=":")[[1]][2])
                                    bq_ref[j]<-as.numeric(strsplit(result_raw2[[j]][6],split=":")[[1]][4])
                                }
                                if(substr(input2.2[i,3],j,j)=="C"){
                                    ref[j]<-as.numeric(strsplit(result_raw2[[j]][7],split=":")[[1]][2])
                                    bq_ref[j]<-as.numeric(strsplit(result_raw2[[j]][7],split=":")[[1]][4])
                                }
                                if(substr(input2.2[i,3],j,j)=="G"){
                                    ref[j]<-as.numeric(strsplit(result_raw2[[j]][8],split=":")[[1]][2])
                                    bq_ref[j]<-as.numeric(strsplit(result_raw2[[j]][8],split=":")[[1]][4])
                                }
                                if(substr(input2.2[i,3],j,j)=="T"){
                                    ref[j]<-as.numeric(strsplit(result_raw2[[j]][9],split=":")[[1]][2])
                                    bq_ref[j]<-as.numeric(strsplit(result_raw2[[j]][9],split=":")[[1]][4])
                                }
                            }
                            alt<-dp-ref
                            results2.2[i,6]<-min(ref)
                            results2.2[i,7]<-max(alt)
                            results2.2[i,8]<-min(dp)
                            results2.2[i,9]<-results2.2[i,7]/results2.2[i,8]
                            results2.2[i,10]<-mean(bq_ref)
                            results2.2[i,11]<-bq_alt
                        }
                    }
                }
            }
            if(length(results1)>0){
              results1<-rbind(results1,results2.2)  
            }
            if(length(results1)==0){
              results1<-results2.2
            }
        }
    }
    write.table(results1,paste(dir1,samples[n,1],".frequency.vcf",sep=""),row.names=F,sep="\t",quote=F)
}





