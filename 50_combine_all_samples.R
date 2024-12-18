library(stringr)
args = commandArgs()
args = args[length(args)]
args = strsplit(args, ",")[[1]]
dir2<-paste(args[1],"/",sep="")
dir1<-paste(dir2,"/documents/",sep="")
samples<-read.table(paste(dir2,"SampleNames.txt",sep=""),header=F,stringsAsFactors = F)

nr_alt<-as.numeric(args[2])
dp<-as.numeric(args[3])
vaf<-as.numeric(args[4])
low_bq<-as.numeric(args[5])
bq_diff<-as.numeric(args[6])

for(n in 1:length(samples[,1])){
    results<-read.table(paste(dir1,samples[n,1],".results_raw.txt",sep=""),header=T,stringsAsFactors =F,quote = "",sep="\t")
    frequency<-read.table(paste(dir1,samples[n,1],".frequency.vcf",sep=""),header=T,stringsAsFactors =F)
    results<-results[order(as.character(results[,2]),as.numeric(results[,3]),as.character(results[,4]),as.character(results[,5])),]
    frequency<-frequency[order(as.character(frequency[,1]),as.numeric(frequency[,2]),as.character(frequency[,3]),as.character(frequency[,4])),]
    exclude<-data.frame(VAF=rep(0,length(results[,1])),BQ=rep(0,length(results[,1])))
    
    message("Sample: ",samples[n,1]," Calls: ",length(results[,1]))
    #Filter for frequency
    temp<-cbind(frequency[,7]<nr_alt,frequency[,8]<dp,frequency[,9]<vaf)
    exclude[,1]<-rowSums(temp,na.rm = T)
    #Filter for low base quality
    temp<-cbind(as.numeric(frequency[,11]<low_bq),as.numeric(frequency[,10])-as.numeric(frequency[,11])>=bq_diff)
    exclude[,2]<-rowSums(temp,na.rm = T)
    
    results_filtered<-results[rowSums(exclude)==0,]
    frequency_filtered<-frequency[rowSums(exclude)==0,]
    
    if(length(results_filtered[,1])>0){
        output<-data.frame(Sample=results_filtered[,1],Chr=results_filtered[,2],Pos=results_filtered[,3],
                           Ref=results_filtered[,4],Alt=results_filtered[,5],
                           Gene=NA,ENST=NA,Type=NA,Mutation=NA,
                           Nr_Ref=frequency_filtered[,6],Nr_Alt=frequency_filtered[,7],
                           DP=frequency_filtered[,8],VAF=frequency_filtered[,9],
                           GATK=results_filtered[,7],Platypus=results_filtered[,8],
                           VarScan=results_filtered[,9],LoFreq=results_filtered[,10],
                           FreeBayes=results_filtered[,11],SNVer=results_filtered[,12],
                           SAMtools=results_filtered[,13],VarDict=results_filtered[,14],
                           Called=rowSums(results_filtered[,7:14],na.rm = T),
                           ESP6500=NA,G1000=NA,Cosmic=NA,Cosmic_location=NA,
                           Cosmic_occurence=NA,dbSNP_SNPs=NA,dbSNP_SNVs=NA,
                           PM_flag=NA,ExAC=NA,ClinVar_clinical=NA,
                           ClinVar_noimpact=NA,Provean_Prediction=NA,Provean_Score=NA,
                           BQ_Ref=frequency_filtered[,10],BQ_Alt=frequency_filtered[,11],
                           Nr_Ref_fwd=frequency_filtered[,12],Nr_Alt_fwd=frequency_filtered[,13],
                           DP_fwd=frequency_filtered[,14],VAF_fwd=frequency_filtered[,15],
                           Nr_Ref_rev=frequency_filtered[,16],Nr_Alt_rev=frequency_filtered[,17],
                           DP_rev=frequency_filtered[,18],VAF_rev=frequency_filtered[,19])
        
        ann<-str_split_fixed(results_filtered[,6],pattern=",",n=Inf)
        
        output[,6]<-str_split_fixed(ann[,1],pattern="\\|",n=Inf)[,6]
        output[,7]<-str_split_fixed(ann[,1],pattern="\\|",n=Inf)[,9]
        output[,8]<-str_split_fixed(ann[,1],pattern="\\(",n=Inf)[,1]
        output[,9]<-paste(str_split_fixed(ann[,1],pattern="\\|",n=Inf)[,3],
                          str_split_fixed(ann[,1],pattern="\\|",n=Inf)[,4],sep="|")
        
        if(length(ann[1,])>1){
            for(i in 2:length(ann[1,])){
                output[,6]<-paste(output[,6],str_split_fixed(ann[,i],pattern="\\|",n=Inf)[,6],sep=",")
                output[,7]<-paste(output[,7],str_split_fixed(ann[,i],pattern="\\|",n=Inf)[,9],sep=",")
                output[,8]<-paste(output[,8],str_split_fixed(ann[,i],pattern="\\(",n=Inf)[,1],sep=",")
                output[,9]<-paste(output[,9],paste(str_split_fixed(ann[,i],pattern="\\|",n=Inf)[,3],
                                                   str_split_fixed(ann[,i],pattern="\\|",n=Inf)[,4],sep="|"),sep=",")
            }
        }
        while(length(grep(",,",output[,6]))>0){
            output[,6]<-gsub(",,",",",output[,6])
            output[,7]<-gsub(",,",",",output[,7])
            output[,8]<-gsub(",,",",",output[,8])
        }
        output[,9]<-gsub(",\\|","",output[,9])
        output[,6]<-gsub(",$","",output[,6])
        output[,7]<-gsub(",$","",output[,7])
        output[,8]<-gsub(",$","",output[,8])
    }
    if(length(results_filtered[,1])==0){
        output<-data.frame(Sample=NA,Chr=NA,Pos=NA,Ref=NA,Alt=NA,
                           Gene=NA,ENST=NA,Type=NA,Mutation=NA,
                           Nr_Ref=NA,Nr_Alt=NA,DP=NA,VAF=NA,GATK=NA,Platypus=NA,
                           VarScan=NA,LoFreq=NA,FreeBayes=NA,SNVer=NA,
                           SAMtools=NA,VarDict=NA,Called=NA,
                           ESP6500=NA,G1000=NA,Cosmic=NA,Cosmic_location=NA,
                           Cosmic_occurence=NA,dbSNP_SNPs=NA,dbSNP_SNVs=NA,
                           PM_flag=NA,ExAC=NA,ClinVar_clinical=NA,
                           ClinVar_noimpact=NA,Provean_Prediction=NA,Provean_Score=NA,
                           BQ_Ref=NA,BQ_Alt=NA,Nr_Ref_fwd=NA,Nr_Alt_fwd=NA,
                           DP_fwd=NA,VAF_fwd=NA,Nr_Ref_rev=NA,Nr_Alt_rev=NA,
                           DP_rev=NA,VAF_rev=NA)
    }
    
    write.table(output,paste(dir1,samples[n,1],".results_filtered_V1_pre.txt",sep=""),row.names=F,sep="\t",quote=F)
}





