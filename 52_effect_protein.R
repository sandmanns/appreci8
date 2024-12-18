library(Biostrings)
library(stringr)
library(DelayedArray)
args = commandArgs()
args = args[length(args)]
args = strsplit(args, ",")[[1]]
dir2<-paste(args[1],"/",sep="")

dir1<-paste(dir2,"/documents/",sep="")

#prepare variants
input_temp<-read.table(paste(dir1,"results_filtered_V1.txt",sep=""),header=T,stringsAsFactors =F,sep="\t",quote = "")
input<-cbind(input_temp[,7],input_temp[,9])
output<-data.frame(Variant=NA)
for(i in 1:length(input[,1])){
    output[i,1]<-strsplit(strsplit(input[i,2],split="p\\.")[[1]][2],split="/")[[1]][1]
}

output2<-gsub("Arg","R",output[,1])
output2<-gsub("His","H",output2)
output2<-gsub("Lys","K",output2)
output2<-gsub("Asp","D",output2)
output2<-gsub("Glu","E",output2)
output2<-gsub("Ser","S",output2)
output2<-gsub("Thr","T",output2)
output2<-gsub("Asn","N",output2)
output2<-gsub("Gln","Q",output2)
output2<-gsub("Cys","C",output2)
output2<-gsub("Sec","U",output2)
output2<-gsub("Gly","G",output2)
output2<-gsub("Pro","P",output2)
output2<-gsub("Ala","A",output2)
output2<-gsub("Val","V",output2)
output2<-gsub("Ile","I",output2)
output2<-gsub("Leu","L",output2)
output2<-gsub("Met","M",output2)
output2<-gsub("Phe","F",output2)
output2<-gsub("Trp","W",output2)
output2<-gsub("Tyr","Y",output2)
output2<-gsub("*","M",output2,fixed = T)

enst<-c(rep(NA,length(input[,1])))
for(i in 1:length(input[,1])){
    enst[i]<-input[i,1]
}
output3<-data.frame(Variant=output2,ENST=enst,stringsAsFactors = F)

write.table(output3,paste(dir1,"results.var",sep=""),row.names=F,sep="\t",quote=F)

#Create ENSTs:
enst_temp<-input_temp[,c(6,7)]
enst_temp2<-enst_temp[nchar(enst_temp[,2])>10,]
enst_temp3<-unique(enst_temp2)
enst_temp4<-data.frame()
for(i in 1:length(enst_temp3[,1])){
    pre<-strsplit(enst_temp3[i,2],split=",")[[1]]
    enst_temp4<-rbind(enst_temp4,as.data.frame(pre))
}
enst<-unique(enst_temp4)
enst<-as.data.frame(enst[enst[,1]!="",])
names(enst)<-c("TranscriptIDs")
write.table(enst,paste(dir2,"/snpEff_ann/ENSTs.txt",sep=""),row.names=F,sep="\t",quote=F)

###prepare transcripts
peptide_info<-args[2]

transcripts<-read.table(paste(dir2,"/snpEff_ann/ENSTs.txt",sep=""),header=T,stringsAsFactors =F)
for(i in 1:length(transcripts[,1])){
    message("Transcript: ",i,"=i; ",transcripts[i,1])
    line_temp<-system(paste("grep -n \"",transcripts[i,1],"\" ",peptide_info,sep=""),intern = T)
    if(length(line_temp)>0){
        line<-strsplit(line_temp,split=":")[[1]][1]
        test<-system(paste("tail -n +",line," ",peptide_info," > ",dir1,"/temp.fa",sep=""),intern=T)
        line2_temp<-system(paste("grep -n -m 2 \"ENSP\" ",dir1,"/temp.fa",sep=""),intern=T)
        line2<-strsplit(line2_temp[2],split=":")[[1]][1]
        system(paste("head -n ",as.character(as.numeric(line2)-1)," ",dir1,"/temp.fa > ",dir2,"/snpEff_ann/",transcripts[i,1],".fasta",sep=""),intern=F)
    }
}

#provean!
#results_final<-read.table(paste(dir1,"results_provean.txt",sep=""),header=T,sep="\t",stringsAsFactors = F)
results_final<-data.frame(output3,provean_score=NA,prediction=NA)
for(i in 1:length(enst[,1])){
    message("ENST ",i," out of ",length(enst[,1]))
    ofInterest<-unique(results_final[grep(enst[i,1],results_final[,2]),1])
    ofInterest<-ofInterest[!is.na(ofInterest)]
    if(length(ofInterest)>0){
        write.table(ofInterest,paste(dir1,"provean_temp.var",sep=""),row.names=F,sep="\t",quote=F)
        system(paste("tail -n +2 ",dir1,"provean_temp.var > ",dir1,"provean_temp2.var",sep=""),intern=F)
        temp<-system(paste("/mnt/home2/sandmans/provean/provean-1.1.5/scripts/provean.sh -q ",dir2,"/snpEff_ann/",enst[i,1],".fasta -v ",dir1,"provean_temp2.var --num_threads 20",sep=""),intern=T)
        if(temp[5]!="No variations entered"&&vcountPattern("File open error",temp[5])==0){
            for(j in 12:length(temp)){
                found<-strsplit(temp[j],split="\\t")[[1]][1]==results_final[,1]
                results_final[!is.na(found)&found==T,3]<-paste(results_final[!is.na(found)&found==T,3],strsplit(temp[j],split="\\t")[[1]][2],sep=",")
                
            }
        }   
    }
}
write.table(results_final,paste(dir1,"results_provean_temp.txt",sep=""),row.names=F,sep="\t",quote=F)

temp<-as.data.frame(str_split_fixed(results_final[,3],pattern=",",n=Inf))
temp2<-temp
for(i in 1:length(temp[1,])){
    temp2[,i]<-as.numeric(as.character(temp[,i]))
}
temp3<-as.matrix.data.frame(temp2)
results_final[,3]<-rowMins(temp3,na.rm=T)
results_final[,3]<-gsub("Inf","NA",results_final[,3])
results_final[,4]<-rowMins(temp3,na.rm=T)<=(-2.5)
results_final[,4]<-gsub("FALSE","Tolerated",results_final[,4])
results_final[,4]<-gsub("TRUE","Deleterious",results_final[,4])
results_final[results_final[,3]=="NA",4]<-"NA"


write.table(results_final,paste(dir1,"results_provean.txt",sep=""),row.names=F,sep="\t",quote=F)

