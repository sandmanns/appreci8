args = commandArgs()
args = args[length(args)]
args = strsplit(args, ",")[[1]]
dir<-paste(args[1],"/documents/test/",sep="")

samples<-read.table(paste(dir,"../../SampleNames.txt",sep=""),header=F,stringsAsFactors=F)
helper2=data.frame(x=1)

for(i in 1:length(samples[,1])){
  helper2=data.frame(x=1)
    message("Sample: ",samples[i,1])
    input<-read.table(paste(dir,samples[i,1],"_temp.uniq.vcf",sep=""),header=F,stringsAsFactors = F)
    output<-data.frame(sample=rep(samples[i,1],length(input[,1])),chr=input[,1],pos=input[,2],
                                  ref=input[,3],alt=input[,4],ann=input[,5],GATK=NA,
                                  Platypus=NA,VarScan=NA,LoFreq=NA,FreeBayes=NA,SNVer=NA,SAMtools=NA,
                                  VarDict=NA)
    helper1<-read.table(paste(dir,samples[i,1],"_temp.uniq2.vcf",sep=""),header=F,stringsAsFactors = F,sep="Z")
    helper1<-helper1[,1]
    
    message("GATK")
    helper2<-read.table(paste(dir,samples[i,1],"_gatk.vcf",sep=""),header=F,stringsAsFactors = F,sep="Z")
    helper2<-helper2[,1]
    output[helper1%in%helper2,7]<-1
    
    message("Platypus")
    helper2<-read.table(paste(dir,samples[i,1],"_platypus.vcf",sep=""),header=F,stringsAsFactors = F,sep="Z")
    helper2<-helper2[,1]
    output[helper1%in%helper2,8]<-1
    
    message("VarScan")
    helper2<-read.table(paste(dir,samples[i,1],"_varscan.vcf",sep=""),header=F,stringsAsFactors = F,sep="Z")
    helper2<-helper2[,1]
    output[helper1%in%helper2,9]<-1

    message("LoFreq")
    helper2<-read.table(paste(dir,samples[i,1],"_lofreq.vcf",sep=""),header=F,stringsAsFactors = F,sep="Z")
    helper2<-helper2[,1]
    output[helper1%in%helper2,10]<-1
    
    message("FreeBayes")
    helper2<-read.table(paste(dir,samples[i,1],"_freebayes.vcf",sep=""),header=F,stringsAsFactors = F,sep="Z")
    helper2<-helper2[,1]
    output[helper1%in%helper2,11]<-1
    
    message("SNVer")
    helper2<-read.table(paste(dir,samples[i,1],"_snver.vcf",sep=""),header=F,stringsAsFactors = F,sep="Z")
    helper2<-helper2[,1]
    output[helper1%in%helper2,12]<-1
    
    message("SAMtools")
    helper2<-read.table(paste(dir,samples[i,1],"_samtools.vcf",sep=""),header=F,stringsAsFactors = F,sep="Z")
    helper2<-helper2[,1]
    output[helper1%in%helper2,13]<-1

    message("VarDict")
    helper2<-read.table(paste(dir,samples[i,1],"_vardict.vcf",sep=""),header=F,stringsAsFactors = F,sep="Z")
    helper2<-helper2[,1]
    output[helper1%in%helper2,14]<-1

    write.table(output,paste(dir,"../",samples[i,1],".results_raw.txt",sep=""),row.names = F,quote=F,sep="\t")
}


