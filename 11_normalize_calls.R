library(BSgenome)
library(stringr)

indel_converter<-function(input){
    input[,4]<-as.character(input[,4])
    input[,5]<-as.character(input[,5])
    minus<-grep(pattern="-",input[,5],fixed = T)
    plus<-grep(pattern="+",input[,5],fixed = T)
    relevant<-c(minus,plus)
    message("Indel_converter: ",length(relevant)," Positionen")
    for(i in relevant){
        if((which(i==relevant)%%100)==0){
            message(which(i==relevant)," von ",length(relevant))
        }
        if(substr(input[i,5],1,1)=="-"){
            ref_new<-paste(input[i,4],substr(input[i,5],2,nchar(input[i,5])),sep="")
            alt_new<-input[i,4]
            input[i,4]<-ref_new
            input[i,5]<-alt_new
        }
        if(substr(input[i,5],1,1)=="+"){
            alt_new<-paste(input[i,4],substr(input[i,5],2,nchar(input[i,5])),sep="")
            input[i,5]<-alt_new
        }
    }
    return(input)
}

mnv_converter<-function(input){
  input[,4]<-as.character(input[,4])
  input[,5]<-as.character(input[,5])
  ncharRef<-nchar(input[,4])
  ncharAlt<-nchar(input[,5])
  relevant<-seq(1,length(input[,1]))
  relevant<-relevant[ncharRef>1&ncharAlt>1&ncharRef==ncharAlt]
  message("MNV_converter: ",length(relevant)," Positionen")
  if(length(relevant)>0){
      temp1<-input[relevant,]
      temp2<-input[setdiff(seq(1,length(input[,1])),relevant),]
      tempRef<-str_split_fixed(string = temp1[,4],pattern="",n = Inf)
      tempAlt<-str_split_fixed(string = temp1[,5],pattern="",n = Inf)
      temp1.2<-cbind(temp1[,1:3],V4=tempRef[,1],Alt=tempAlt[,1])
      temp1.2[,4]<-as.character(temp1.2[,4])
      temp1.2[,5]<-as.character(temp1.2[,5])
      temp1[,2]<-as.numeric(as.character(temp1[,2]))
      if(length(tempRef[1,])>1){
          for(i in 2:length(tempRef[1,])){
              add<-cbind(V1=as.character(temp1[,1]),V2=(temp1[,2]+i-1),V3=".",
                         V4=tempRef[,i],Alt=tempAlt[,i])
              add[,4]<-as.character(add[,4])
              add[,5]<-as.character(add[,5])
              add<-add[add[,5]!="",]
              temp1.2<-rbind(temp1.2,add)
          }
      }
      temp1.3<-temp1.2[temp1.2[,4]!=temp1.2[,5],]
      input_new<-rbind(temp2,temp1.3)
      return(input_new)
  }
  return(input)
}

check_alternative_bases<-function(input){
  input[,4]<-as.character(input[,4])
  input[,5]<-as.character(input[,5])
  relevant<-seq(1,length(input[,1]))
  relevant<-relevant[grep(pattern=",",input[,5],fixed=T)]
  message("check bases: ",length(relevant)," Positionen")
  temp<-str_split_fixed(string = input[,5],pattern=",",n = Inf)
  temp2<-cbind(input[,c(1:4)],Alt=temp[,1])
  temp2[,5]<-as.character(temp2[,5])
  if(length(temp[1,])>1){
    for(i in 2:length(temp[1,])){
      add<-cbind(input[,c(1:4)],Alt=temp[,i])
      add[,5]<-as.character(add[,5])
      add<-add[add[,5]!="",]
      temp2<-rbind(temp2,add)
    }
  }
  return(temp2)
}

string_diff_finder<-function(input){
    #start from the back
    input[,4]<-as.character(input[,4])
    input[,5]<-as.character(input[,5])
    ncharRef<-nchar(input[,4])
    ncharAlt<-nchar(input[,5])
    relevant<-seq(1,length(input[,1]))
    relevant<-relevant[ncharRef>1&ncharAlt>1]
    message("String_diff_finder: ",length(relevant)," Positionen")
    for(i in relevant){
        if((which(i==relevant)%%100)==0){
            message(which(i==relevant)," von ",length(relevant))
        }
        ref<-strsplit(substr(input[i,4],2,(ncharRef[i])),split="")[[1]]
        alt<-strsplit(substr(input[i,5],2,(ncharAlt[i])),split="")[[1]]
        flag_diff<-F
        j_ref<-length(ref)
        j_alt<-length(alt)
        while(flag_diff==F&&j_ref>0&&j_alt>0){
            if(ref[j_ref]!=alt[j_alt]){
                flag_diff<-T
            }
            if(flag_diff==F){
                j_ref<-j_ref-1
                j_alt<-j_alt-1
            }
        }
        input[i,4]<-substr(input[i,4],1,(j_ref+1))
        input[i,5]<-substr(input[i,5],1,(j_alt+1))
    }
    #now check the front
    ncharRef<-nchar(input[,4])
    ncharAlt<-nchar(input[,5])
    relevant<-seq(1,length(input[,1]))
    relevant<-relevant[ncharRef>1&ncharAlt>1]
    input[,2]<-as.numeric(as.character(input[,2]))
    for(i in relevant){
        if((which(i==relevant)%%1000)==0){
            message(which(i==relevant)," von ",length(relevant))
        }
        flag_diff<-F
        while(flag_diff==F&&nchar(input[i,4])>1&&nchar(input[i,5])>1){
            ref<-substr(input[i,4],1,1)[[1]]
            alt<-substr(input[i,5],1,1)[[1]]
            if(ref[1]==alt[1]){
                input[i,2]<-input[i,2]+1
                input[i,4]<-substr(input[i,4],2,nchar(input[i,4]))
                input[i,5]<-substr(input[i,5],2,nchar(input[i,5]))
            }
            if(ref[1]!=alt[1]){
                flag_diff<-T
            }
        }
    }
    return(input)
}


args = commandArgs()
args = args[length(args)]
args = strsplit(args, ",")[[1]]
gatk<-args[1]
platypus<-args[2]
varscan<-args[3]
lofreq<-args[4]
freebayes<-args[5]
snver<-args[6]
samtools<-args[7]
vardict<-args[8]

dir1<-paste(args[9],"/",sep="")
samples<-read.table(paste(dir1,"SampleNames.txt",sep=""),header=F)
snps=data.frame(x=1)

if(gatk==T){
    message("GATK")
    dir2<-paste(dir1,"/gatk/intermRes/calling_comp/",sep="") 
    for(i in 1:length(samples[,1])){
        message("Sample ",i," von ",length(samples[,1]))
        results<-data.frame(stringsAsFactors=F)
        data<-lapply(snps,function(x){
            tryCatch(read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F),error=function(e) NULL)
        })
        if(!is.null(data[[1]])){
            snps<-read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F,
                            colClasses = c(NULL,NULL,NULL,"character","character",NULL,NULL,NULL,NULL,NULL))
            snps<-read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F)
            temp2<-check_alternative_bases(snps)
            temp<-temp2
            temp2<-string_diff_finder(temp)
            results<-temp2 
        }
        write.table(results,paste(dir2,"/",samples[i,1],".gatk_raw2.vcf",sep=""),row.names=F,quote=F,sep="\t")
    }
}

if(platypus==T){
    message("Platypus")
    dir2<-paste(dir1,"/Platypus/calling_comp/",sep="")
    for(i in 1:length(samples[,1])){
        message("Sample ",i," von ",length(samples[,1]))
        results<-data.frame(stringsAsFactors=F)
                data<-lapply(snps,function(x){
                   tryCatch(read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F),error=function(e) NULL)
              })
             if(!is.null(data[[1]])){
        snps<-read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F,
                        colClasses = c(NULL,NULL,NULL,"character","character",NULL,NULL,NULL,NULL,NULL))
        snps<-read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F)
        temp2<-check_alternative_bases(snps)
        temp<-temp2
        temp2<-string_diff_finder(temp)
        results<-temp2 
        }
        write.table(results,paste(dir2,"/",samples[i,1],".Platypus_raw2.vcf",sep=""),row.names=F,quote=F,sep="\t")
    }
}

if(varscan==T){
    message("VarScan")
    dir2<-paste(dir1,"/varscan/calling_comp/",sep="")
    for(i in 1:length(samples[,1])){
        message("Sample ",i," von ",length(samples[,1]))
                results<-data.frame(stringsAsFactors=F)
               data<-lapply(snps,function(x){
                  tryCatch(read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F),error=function(e) NULL)
             })
            if(!is.null(data[[1]])){
        snps<-read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F,
                        colClasses = c(NULL,NULL,NULL,"character","character",NULL,NULL,NULL,NULL,NULL))
        snps<-read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F)
        temp2<-indel_converter(snps)
        temp3<-check_alternative_bases(temp2)
        temp<-temp3
        temp2<-string_diff_finder(temp)
        results<-temp2 
        }
        write.table(results,paste(dir2,"/",samples[i,1],".varscan_raw2.vcf",sep=""),row.names=F,quote=F,sep="\t")
    }
}

if(lofreq==T){
    message("LoFreq")
    dir2<-paste(dir1,"/lofreq/calling_comp/",sep="")
    for(i in 1:length(samples[,1])){
        message("Sample ",i," von ",length(samples[,1]))
        results<-data.frame(stringsAsFactors=F)
                data<-lapply(snps,function(x){
                    tryCatch(read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F),error=function(e) NULL)
                })
               if(!is.null(data[[1]])){
        snps<-read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F,
                        colClasses = c(NULL,NULL,NULL,"character","character",NULL,NULL,NULL,NULL,NULL))
        snps<-read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F)
        temp2<-check_alternative_bases(snps)
        temp<-temp2
        temp2<-string_diff_finder(temp)
        results<-temp2 
        }
        write.table(results,paste(dir2,"/",samples[i,1],".lofreq_raw2.vcf",sep=""),row.names=F,quote=F,sep="\t")
    }
}

if(snver==T){
    message("SNVer")
    dir2<-paste(dir1,"/SNVer/calling_comp/",sep="")
    for(i in 1:length(samples[,1])){
        message("Sample ",i," von ",length(samples[,1]))
        
        results<-data.frame(stringsAsFactors=F)
        data<-lapply(snps,function(x){
            tryCatch(read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F),error=function(e) NULL)
        })
        if(!is.null(data[[1]])){
        snps<-read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F,
                        colClasses = c(NULL,NULL,NULL,"character","character",NULL,NULL,NULL,NULL,NULL))
        snps<-read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F)
        temp2<-check_alternative_bases(snps)
        temp<-temp2
        temp2<-string_diff_finder(temp)
        results<-temp2 
        }
        write.table(results,paste(dir2,"/",samples[i,1],".snver_raw2.vcf",sep=""),row.names=F,quote=F,sep="\t")
    }
}

if(samtools==T){
    message("SAMtools")
    dir2<-paste(dir1,"/samtools/calling_comp/",sep="")
    for(i in 1:length(samples[,1])){
        message("Sample ",i," von ",length(samples[,1]))
        results<-data.frame(stringsAsFactors=F)
        data<-lapply(snps,function(x){
            tryCatch(read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F),error=function(e) NULL)
        })
        if(!is.null(data[[1]])){
        snps<-read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F,
                        colClasses = c(NULL,NULL,NULL,"character","character",NULL,NULL,NULL,NULL,NULL))
        snps<-read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F)
        temp2<-check_alternative_bases(snps)
        temp<-temp2
        temp2<-string_diff_finder(temp)
        results<-temp2 
        }
        write.table(results,paste(dir2,"/",samples[i,1],".samtools_raw2.vcf",sep=""),row.names=F,quote=F,sep="\t")
    }
}

if(vardict==T){
    message("VarDict")
    dir2<-paste(dir1,"/vardict/calling_comp/",sep="")
    for(i in 1:length(samples[,1])){
        message("Sample ",i," von ",length(samples[,1]))
        results<-data.frame(stringsAsFactors=F)
        data<-lapply(snps,function(x){
            tryCatch(read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F),error=function(e) NULL)
        })
        if(!is.null(data[[1]])){
        snps<-read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F,
                        colClasses = c(NULL,NULL,NULL,"character","character",NULL,NULL,NULL,NULL,NULL))
        snps<-read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F)
        temp2<-check_alternative_bases(snps)
        temp<-temp2
        temp2<-mnv_converter(temp)
        temp<-temp2
        temp2<-string_diff_finder(temp)
        results<-temp2 
        }
        write.table(results,paste(dir2,"/",samples[i,1],".vardict_raw2.vcf",sep=""),row.names=F,quote=F,sep="\t")
    }
}

if(freebayes==T){
    message("FreeBayes")
    dir2<-paste(dir1,"/freebayes/calling_comp/",sep="")
    for(i in 1:length(samples[,1])){
        message("Sample ",i," von ",length(samples[,1]))
                results<-data.frame(stringsAsFactors=F)
               data<-lapply(snps,function(x){
                   tryCatch(read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F),error=function(e) NULL)
               })
              if(!is.null(data[[1]])){
        snps<-read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F,
                        colClasses = c(NULL,NULL,NULL,"character","character",NULL,NULL,NULL,NULL,NULL))
        snps<-read.table(paste(dir2,samples[i,1],".target.vcf",sep=""),header=F)
        temp2<-check_alternative_bases(snps)
        temp<-temp2
        temp2<-mnv_converter(temp)
        temp<-temp2
        temp2<-string_diff_finder(temp)
        results<-temp2 
        }
        write.table(results,paste(dir2,"/",samples[i,1],".freebayes_raw2.vcf",sep=""),row.names=F,quote=F,sep="\t")
    }
}