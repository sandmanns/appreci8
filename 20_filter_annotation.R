library(stringr)

filterAnnotation<-function(input){
    temp<-str_split_fixed(string = input[,8],pattern=",",n = Inf)
    temp[grep("5_prime_UTR_variant",temp)]<-""
    temp[grep("3_prime_UTR_variant",temp)]<-""
    temp[grep("downstream_gene_variant",temp)]<-""
    temp[grep("upstream_gene_variant",temp)]<-""
    temp[grep("=intron_variant",temp)]<-""
    temp[grep("^intron_variant",temp)]<-""
    temp[grep("intergenic_variant",temp)]<-""
    temp[grep("intragenic_variant",temp)]<-""
    temp[grep("intergenic_region",temp)]<-""
    #temp[grep("synonymous_variant",temp)]<-""
    temp[grep("non_coding_exon_variant",temp)]<-""
    temp[grep("protein_protein_contact",temp)]<-""
    temp[grep("sequence_feature",temp)]<-""
    temp[grep("^[0-999]",temp)]<-""
    temp[grep("TF_binding_site_variant",temp)]<-""
    
    test<-seq(1,length(input[,1]))[rowSums(temp!="")>0]
    temp<-temp[test,]
    if(is.vector(temp)==T){
        temp<-t(data.frame(V1=temp))
    }
    input<-input[test,]
    temp<-apply(temp,1,paste,collapse=",")
    temp<-gsub(pattern="EFF=",replacement="",temp)
    while(length(grep(",,",temp))>0){
        temp<-gsub(pattern=",,",replacement=",",temp)
    }
    temp<-gsub(pattern=",$",replacement="",temp)
    temp<-gsub(pattern="^,",replacement="",temp)
    results<-cbind(input[,1:7],temp,input[,9])
    return(results)
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
input=data.frame(x=1)

if(gatk==T){
    message("GATK")
    dir2<-paste(dir1,"/gatk/intermRes/calling_comp/",sep="") 
    for(i in 1:length(samples[,1])){
        message("Sample ",i," von ",length(samples[,1]))
      filtered<-data.frame(stringsAsFactors=F)
      data<-lapply(input,function(x){
        tryCatch(read.table(paste(dir2,samples[i,1],".gatk_ann.vcf",sep=""),header=F),error=function(e) NULL)
      })
      if(!is.null(data[[1]])){
        input<-read.table(paste(dir2,samples[i,1],".gatk_ann.vcf",sep=""),header=F)
	input<-input[input[,8]!=".",]
        filtered<-filterAnnotation(input)
      }
        write.table(filtered,paste(dir2,"/",samples[i,1],".gatk_filtered.vcf",sep=""),row.names=F,quote=F,sep="\t")
    }
}

if(platypus==T){
    message("Platypus")
    dir2<-paste(dir1,"/Platypus/calling_comp/",sep="")
    for(i in 1:length(samples[,1])){
        message("Sample ",i," von ",length(samples[,1]))
      filtered<-data.frame(stringsAsFactors=F)
      data<-lapply(input,function(x){
        tryCatch(read.table(paste(dir2,samples[i,1],".Platypus_ann.vcf",sep=""),header=F),error=function(e) NULL)
      })
      if(!is.null(data[[1]])){
        input<-read.table(paste(dir2,samples[i,1],".Platypus_ann.vcf",sep=""),header=F)
	input<-input[input[,8]!=".",]
        filtered<-filterAnnotation(input)
      }
        write.table(filtered,paste(dir2,"/",samples[i,1],".Platypus_filtered.vcf",sep=""),row.names=F,quote=F,sep="\t")
    }
}

if(varscan==T){
    message("VarScan")
    dir2<-paste(dir1,"/varscan/calling_comp/",sep="")
    for(i in 1:length(samples[,1])){
        message("Sample ",i," von ",length(samples[,1]))
      filtered<-data.frame(stringsAsFactors=F)
      data<-lapply(input,function(x){
        tryCatch(read.table(paste(dir2,samples[i,1],".varscan_ann.vcf",sep=""),header=F),error=function(e) NULL)
      })
      if(!is.null(data[[1]])){
        input<-read.table(paste(dir2,samples[i,1],".varscan_ann.vcf",sep=""),header=F)
	input<-input[input[,8]!=".",]
        filtered<-filterAnnotation(input)
      }
        write.table(filtered,paste(dir2,"/",samples[i,1],".varscan_filtered.vcf",sep=""),row.names=F,quote=F,sep="\t")
    }
}

if(lofreq==T){
    message("LoFreq")
    dir2<-paste(dir1,"/lofreq/calling_comp/",sep="")
    for(i in 1:length(samples[,1])){
        message("Sample ",i," von ",length(samples[,1]))
      filtered<-data.frame(stringsAsFactors=F)
      data<-lapply(input,function(x){
        tryCatch(read.table(paste(dir2,samples[i,1],".lofreq_ann.vcf",sep=""),header=F),error=function(e) NULL)
      })
      if(!is.null(data[[1]])){
        input<-read.table(paste(dir2,samples[i,1],".lofreq_ann.vcf",sep=""),header=F)
	input<-input[input[,8]!=".",]
        filtered<-filterAnnotation(input)
      }
        write.table(filtered,paste(dir2,"/",samples[i,1],".lofreq_filtered.vcf",sep=""),row.names=F,quote=F,sep="\t")
    }
}

if(snver==T){
    message("SNVer")
    dir2<-paste(dir1,"/SNVer/calling_comp/",sep="")
    for(i in 1:length(samples[,1])){
        message("Sample ",i," von ",length(samples[,1]))
      filtered<-data.frame(stringsAsFactors=F)
      data<-lapply(input,function(x){
        tryCatch(read.table(paste(dir2,samples[i,1],".snver_ann.vcf",sep=""),header=F),error=function(e) NULL)
      })
      if(!is.null(data[[1]])){
        input<-read.table(paste(dir2,samples[i,1],".snver_ann.vcf",sep=""),header=F)
	input<-input[input[,8]!=".",]
        filtered<-filterAnnotation(input)
      }
        write.table(filtered,paste(dir2,"/",samples[i,1],".snver_filtered.vcf",sep=""),row.names=F,quote=F,sep="\t")
    }
}

if(samtools==T){
    message("SAMtools")
    dir2<-paste(dir1,"/samtools/calling_comp/",sep="")
    for(i in 1:length(samples[,1])){
        message("Sample ",i," von ",length(samples[,1]))
      filtered<-data.frame(stringsAsFactors=F)
      data<-lapply(input,function(x){
        tryCatch(read.table(paste(dir2,samples[i,1],".samtools_ann.vcf",sep=""),header=F),error=function(e) NULL)
      })
      if(!is.null(data[[1]])){
        input<-read.table(paste(dir2,samples[i,1],".samtools_ann.vcf",sep=""),header=F)
	input<-input[input[,8]!=".",]
        filtered<-filterAnnotation(input)
      }
        write.table(filtered,paste(dir2,"/",samples[i,1],".samtools_filtered.vcf",sep=""),row.names=F,quote=F,sep="\t")
    }
}

if(vardict==T){
    message("VarDict")
    dir2<-paste(dir1,"/vardict/calling_comp/",sep="")
    for(i in 1:length(samples[,1])){
        message("Sample ",i," von ",length(samples[,1]))
      filtered<-data.frame(stringsAsFactors=F)
      data<-lapply(input,function(x){
        tryCatch(read.table(paste(dir2,samples[i,1],".vardict_ann.vcf",sep=""),header=F),error=function(e) NULL)
      })
      if(!is.null(data[[1]])){
        input<-read.table(paste(dir2,samples[i,1],".vardict_ann.vcf",sep=""),header=F)
	input<-input[input[,8]!=".",]
        filtered<-filterAnnotation(input)
      }
        write.table(filtered,paste(dir2,"/",samples[i,1],".vardict_filtered.vcf",sep=""),row.names=F,quote=F,sep="\t")
    }
}

if(freebayes==T){
    message("FreeBayes")
    dir2<-paste(dir1,"/freebayes/calling_comp/",sep="")
    for(i in 1:length(samples[,1])){
        message("Sample ",i," von ",length(samples[,1]))
      filtered<-data.frame(stringsAsFactors=F)
      data<-lapply(input,function(x){
        tryCatch(read.table(paste(dir2,samples[i,1],".freebayes_ann.vcf",sep=""),header=F),error=function(e) NULL)
      })
      if(!is.null(data[[1]])){
        input<-read.table(paste(dir2,samples[i,1],".freebayes_ann.vcf",sep=""),header=F)
	input<-input[input[,8]!=".",]
        filtered<-filterAnnotation(input)
      }
        write.table(filtered,paste(dir2,"/",samples[i,1],".freebayes_filtered.vcf",sep=""),row.names=F,quote=F,sep="\t")
    }
}
