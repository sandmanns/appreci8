library("VariantAnnotation")
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library("Rsamtools")
library(Biostrings)

args = commandArgs()
args = args[length(args)]
args = strsplit(args, ",")[[1]]
dir2<-paste(args[1],"/documents/",sep="")
samples<-read.table(paste(dir2,"../SampleNames.txt",sep=""),header=F,stringsAsFactors = F)

snps_pre<-data.frame()
for(i in 1:length(samples[,1])){
    temp<-read.table(paste(dir2,samples[i,1],".results_filtered_V1_pre.txt",sep=""),header=T,sep="\t",quote = "",stringsAsFactors = F)
    if(!is.na(temp[1,1])){
        snps_pre<-rbind(snps_pre,temp)   
    }
}
snps_pre<-snps_pre[order(as.character(snps_pre[,2]),as.numeric(snps_pre[,3]),snps_pre[,4],snps_pre[,5],snps_pre[,1]),]
snps<-cbind(snps_pre[,c(2,3)],ID=rep(".",length(snps_pre[,1])),snps_pre[,c(4,5,1,6,9,7)])
snps[,4]<-gsub("TRUE","T",snps[,4])
dir3<-as.character(args[2])


checkMutation<-function(ref1,alt1,ref2,alt2){
    if((ref1==ref2&&alt1==alt2)||
       (as.character(complement(DNAString(ref1)))==ref2&&
        as.character(complement(DNAString(alt1)))==alt2)){
        return(TRUE)
    }
    if(vcountPattern(",",alt1)>0&&vcountPattern(",",alt2)==0){
        return(FALSE)
    }
    if(vcountPattern(",",alt1)==0&&vcountPattern(",",alt2)>0){
        variants<-strsplit(alt2,",")[[1]]
        for(i in 1:length(variants)){
            if(variants[i]==alt1){
                return(TRUE)
            }
        }
        return(FALSE)
    }
    if(vcountPattern(",",alt1)>0&&vcountPattern(",",alt2)>0){
        variants1<-strsplit(alt1,",")[[1]]
        variants2<-strsplit(alt2,",")[[1]]
        found<-0
        for(i in 1:length(variants1)){
            for(j in 1:length(variants2)){
                if(variants1[i]==variants2[j]){
                    found<-found+1
                }
            }
            if(found1==length(variants1)){
                return(TRUE)
            }
            if(found1<length(variants1)){
                return(FALSE)
            }
        }
    }
    return(FALSE)
}


db_esp<-as.character(args[10])
ESP_1<-TabixFile(file=paste(dir3,db_esp,".chr1.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr1.snps_indels.vcf.gz.tbi",sep=""))
ESP_2<-TabixFile(file=paste(dir3,db_esp,".chr2.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr2.snps_indels.vcf.gz.tbi",sep=""))
ESP_3<-TabixFile(file=paste(dir3,db_esp,".chr3.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr3.snps_indels.vcf.gz.tbi",sep=""))
ESP_4<-TabixFile(file=paste(dir3,db_esp,".chr4.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr4.snps_indels.vcf.gz.tbi",sep=""))
ESP_5<-TabixFile(file=paste(dir3,db_esp,".chr5.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr5.snps_indels.vcf.gz.tbi",sep=""))
ESP_6<-TabixFile(file=paste(dir3,db_esp,".chr6.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr6.snps_indels.vcf.gz.tbi",sep=""))
ESP_7<-TabixFile(file=paste(dir3,db_esp,".chr7.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr7.snps_indels.vcf.gz.tbi",sep=""))
ESP_8<-TabixFile(file=paste(dir3,db_esp,".chr8.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr8.snps_indels.vcf.gz.tbi",sep=""))
ESP_9<-TabixFile(file=paste(dir3,db_esp,".chr9.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr9.snps_indels.vcf.gz.tbi",sep=""))
ESP_10<-TabixFile(file=paste(dir3,db_esp,".chr10.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr10.snps_indels.vcf.gz.tbi",sep=""))
ESP_11<-TabixFile(file=paste(dir3,db_esp,".chr11.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr11.snps_indels.vcf.gz.tbi",sep=""))
ESP_12<-TabixFile(file=paste(dir3,db_esp,".chr12.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr12.snps_indels.vcf.gz.tbi",sep=""))
ESP_13<-TabixFile(file=paste(dir3,db_esp,".chr13.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr13.snps_indels.vcf.gz.tbi",sep=""))
ESP_14<-TabixFile(file=paste(dir3,db_esp,".chr14.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr14.snps_indels.vcf.gz.tbi",sep=""))
ESP_15<-TabixFile(file=paste(dir3,db_esp,".chr15.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr15.snps_indels.vcf.gz.tbi",sep=""))
ESP_16<-TabixFile(file=paste(dir3,db_esp,".chr16.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr16.snps_indels.vcf.gz.tbi",sep=""))
ESP_17<-TabixFile(file=paste(dir3,db_esp,".chr17.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr17.snps_indels.vcf.gz.tbi",sep=""))
ESP_18<-TabixFile(file=paste(dir3,db_esp,".chr18.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr18.snps_indels.vcf.gz.tbi",sep=""))
ESP_19<-TabixFile(file=paste(dir3,db_esp,".chr19.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr19.snps_indels.vcf.gz.tbi",sep=""))
ESP_20<-TabixFile(file=paste(dir3,db_esp,".chr20.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr20.snps_indels.vcf.gz.tbi",sep=""))
ESP_21<-TabixFile(file=paste(dir3,db_esp,".chr21.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr21.snps_indels.vcf.gz.tbi",sep=""))
ESP_22<-TabixFile(file=paste(dir3,db_esp,".chr22.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chr22.snps_indels.vcf.gz.tbi",sep=""))
ESP_X<-TabixFile(file=paste(dir3,db_esp,".chrX.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chrX.snps_indels.vcf.gz.tbi",sep=""))
ESP_Y<-TabixFile(file=paste(dir3,db_esp,".chrY.snps_indels.vcf.gz",sep=""),index=paste(dir3,db_esp,".chrY.snps_indels.vcf.gz.tbi",sep=""))

db_g1000<-as.character(args[12])
db_g1000x<-as.character(args[13])
db_g1000y<-as.character(args[14])
G1000_1<-TabixFile(file=paste(dir3,"ALL.","chr1.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr1.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_2<-TabixFile(file=paste(dir3,"ALL.","chr2.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr2.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_3<-TabixFile(file=paste(dir3,"ALL.","chr3.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr3.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_4<-TabixFile(file=paste(dir3,"ALL.","chr4.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr4.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_5<-TabixFile(file=paste(dir3,"ALL.","chr5.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr5.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_6<-TabixFile(file=paste(dir3,"ALL.","chr6.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr6.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_7<-TabixFile(file=paste(dir3,"ALL.","chr7.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr7.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_8<-TabixFile(file=paste(dir3,"ALL.","chr8.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr8.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_9<-TabixFile(file=paste(dir3,"ALL.","chr9.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr9.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_10<-TabixFile(file=paste(dir3,"ALL.","chr10.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr10.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_11<-TabixFile(file=paste(dir3,"ALL.","chr11.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr11.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_12<-TabixFile(file=paste(dir3,"ALL.","chr12.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr12.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_13<-TabixFile(file=paste(dir3,"ALL.","chr13.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr13.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_14<-TabixFile(file=paste(dir3,"ALL.","chr14.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr14.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_15<-TabixFile(file=paste(dir3,"ALL.","chr15.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr15.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_16<-TabixFile(file=paste(dir3,"ALL.","chr16.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr16.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_17<-TabixFile(file=paste(dir3,"ALL.","chr17.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr17.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_18<-TabixFile(file=paste(dir3,"ALL.","chr18.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr18.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_19<-TabixFile(file=paste(dir3,"ALL.","chr19.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr19.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_20<-TabixFile(file=paste(dir3,"ALL.","chr20.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr20.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_21<-TabixFile(file=paste(dir3,"ALL.","chr21.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr21.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_22<-TabixFile(file=paste(dir3,"ALL.","chr22.",db_g1000,".vcf.gz",sep=""),index=paste(dir3,"ALL.","chr22.",db_g1000,".vcf.gz.tbi",sep=""))
G1000_X<-TabixFile(file=paste(dir3,db_g1000x,".vcf.gz",sep=""),index=paste(dir3,db_g1000x,".vcf.gz.tbi",sep=""))
G1000_Y<-TabixFile(file=paste(dir3,db_g1000y,".vcf.gz",sep=""),index=paste(dir3,db_g1000y,".vcf.gz.tbi",sep=""))

CosmicCoding<-TabixFile(file=paste(dir3,as.character(args[6]),".vcf.gz",sep=""),index=paste(dir3,as.character(args[6]),".vcf.gz.tbi",sep=""))
CosmicNoncoding<-TabixFile(file=paste(dir3,as.character(args[7]),".vcf.gz",sep=""),index=paste(dir3,as.character(args[7]),".vcf.gz.tbi",sep=""))
Cosmic<-TabixFile(file=paste(dir3,as.character(args[5]),".bed.gz",sep=""),index=paste(dir3,as.character(args[5]),".bed.gz.tbi",sep=""))

dbSNP_poly<-TabixFile(file=paste(dir3,as.character(args[9]),".vcf.gz",sep=""),index=paste(dir3,as.character(args[9]),".vcf.gz.tbi",sep=""))
dbSNP<-TabixFile(file=paste(dir3,as.character(args[8]),".vcf.gz",sep=""),index=paste(dir3,as.character(args[8]),".vcf.gz.tbi",sep=""))
EXAC<-TabixFile(file=paste(dir3,as.character(args[11]),".vcf.gz",sep=""),index=paste(dir3,as.character(args[11]),".vcf.gz.tbi",sep=""))
clinvar_rel<-TabixFile(file=paste(dir3,as.character(args[3]),".vcf.gz",sep=""),index=paste(dir3,as.character(args[3]),".vcf.gz.tbi",sep=""))
clinvar_nons<-TabixFile(file=paste(dir3,as.character(args[4]),".vcf.gz",sep=""),index=paste(dir3,as.character(args[4]),".vcf.gz.tbi",sep=""))


#######
#SNPs:#
#######
firstVariant<-T
for(i in 1:length(snps[,1])){
    if(i!=1){
        if(as.character(snps[i,1])==as.character(snps[i-1,1])&&
           as.numeric(snps[i,2])==as.numeric(snps[i-1,2])&&
           as.character(snps[i,4])==as.character(snps[i-1,4])&&
           as.character(snps[i,5])==as.character(snps[i-1,5])){
            firstVariant<-F
            snps[i,10]<-snps[i-1,10]
            snps[i,11]<-snps[i-1,11]
            snps[i,12]<-snps[i-1,12]
            snps[i,13]<-snps[i-1,13]
            snps[i,14]<-snps[i-1,14]
            snps[i,15]<-snps[i-1,15]
            snps[i,16]<-snps[i-1,16]
            snps[i,17]<-snps[i-1,17]
            snps[i,18]<-snps[i-1,18]
            snps[i,19]<-snps[i-1,19]
            snps[i,20]<-snps[i-1,20]
        }
    }
    if(firstVariant==T){
        chr<-as.character(snps[i,1])
        pos<-as.numeric(snps[i,2])
        param<-GRanges(seqnames=chr,ranges=IRanges(pos,pos))
        ref<-as.character(snps[i,4])
        alt<-as.character(snps[i,5])
        if(vcountPattern(",",alt)>0||vcountPattern(",",ref)>0){
            snps[i,10]<-"Manual inspection"
        }
        if(vcountPattern(",",alt)==0&&vcountPattern(",",ref)==0){
            if(chr=="1"){
                res_esp<-scanTabix(ESP_1,param=param)
                res_1000g<-scanTabix(G1000_1,param=param)
            }
            if(chr=="2"){
                res_esp<-scanTabix(ESP_2,param=param)
                res_1000g<-scanTabix(G1000_2,param=param)
            }
            if(chr=="3"){
                res_esp<-scanTabix(ESP_3,param=param)
                res_1000g<-scanTabix(G1000_3,param=param)
            }
            if(chr=="4"){
                res_esp<-scanTabix(ESP_4,param=param)
                res_1000g<-scanTabix(G1000_4,param=param)
            }
            if(chr=="5"){
                res_esp<-scanTabix(ESP_5,param=param)
                res_1000g<-scanTabix(G1000_5,param=param)
            }
            if(chr=="6"){
                res_esp<-scanTabix(ESP_6,param=param)
                res_1000g<-scanTabix(G1000_6,param=param)
            }
            if(chr=="7"){
                res_esp<-scanTabix(ESP_7,param=param)
                res_1000g<-scanTabix(G1000_7,param=param)
            }
            if(chr=="8"){
                res_esp<-scanTabix(ESP_8,param=param)
                res_1000g<-scanTabix(G1000_8,param=param)
            }
            if(chr=="9"){
                res_esp<-scanTabix(ESP_9,param=param)
                res_1000g<-scanTabix(G1000_9,param=param)
            }
            if(chr=="10"){
                res_esp<-scanTabix(ESP_10,param=param)
                res_1000g<-scanTabix(G1000_10,param=param)
            }
            if(chr=="11"){
                res_esp<-scanTabix(ESP_11,param=param)
                res_1000g<-scanTabix(G1000_11,param=param)
            }
            if(chr=="12"){
                res_esp<-scanTabix(ESP_12,param=param)
                res_1000g<-scanTabix(G1000_12,param=param)
            }
            if(chr=="13"){
                res_esp<-scanTabix(ESP_13,param=param)
                res_1000g<-scanTabix(G1000_13,param=param)
            }
            if(chr=="14"){
                res_esp<-scanTabix(ESP_14,param=param)
                res_1000g<-scanTabix(G1000_14,param=param)
            }
            if(chr=="15"){
                res_esp<-scanTabix(ESP_15,param=param)
                res_1000g<-scanTabix(G1000_15,param=param)
            }
            if(chr=="16"){
                res_esp<-scanTabix(ESP_16,param=param)
                res_1000g<-scanTabix(G1000_16,param=param)
            }
            if(chr=="17"){
                res_esp<-scanTabix(ESP_17,param=param)
                res_1000g<-scanTabix(G1000_17,param=param)
            }
            if(chr=="18"){
                res_esp<-scanTabix(ESP_18,param=param)
                res_1000g<-scanTabix(G1000_18,param=param)
            }
            if(chr=="19"){
                res_esp<-scanTabix(ESP_19,param=param)
                res_1000g<-scanTabix(G1000_19,param=param)
            }
            if(chr=="20"){
                res_esp<-scanTabix(ESP_20,param=param)
                res_1000g<-scanTabix(G1000_20,param=param)
            }
            if(chr=="21"){
                res_esp<-scanTabix(ESP_21,param=param)
                res_1000g<-scanTabix(G1000_21,param=param)
            }
            if(chr=="22"){
                res_esp<-scanTabix(ESP_22,param=param)
                res_1000g<-scanTabix(G1000_22,param=param)
            }
            if(chr=="X"){
                res_esp<-scanTabix(ESP_X,param=param)
                res_1000g<-scanTabix(G1000_X,param=param)
            }
            if(chr=="Y"){
                res_esp<-scanTabix(ESP_Y,param=param)
                res_1000g<-scanTabix(G1000_Y,param=param)
            }
            
            res_cosCC<-scanTabix(CosmicCoding,param=param)
            res_cosNC<-scanTabix(CosmicNoncoding,param=param)
            res_snp<-scanTabix(dbSNP_poly,param=param)
            res_snv<-scanTabix(dbSNP,param=param)
            res_exac<-scanTabix(EXAC,param=param)
            if(as.character(snps[i,1])!="X"&&as.character(snps[i,1])!="Y"){
                res_cvR<-scanTabix(clinvar_rel,param=param)
                res_cvN<-scanTabix(clinvar_nons,param=param)
            }
            message("Mutation ",i," out of ",length(snps[,1]),":")
            
            if(length(res_esp[[1]])==0||strsplit(res_esp[[1]],split="\\t")[[1]][2]!=pos||
               !checkMutation(ref,alt,strsplit(res_esp[[1]],split="\\t")[[1]][4],strsplit(res_esp[[1]],split="\\t")[[1]][5])){
                snps[i,10]<-"NO"
            }
            if(length(res_esp[[1]])>0&&strsplit(res_esp[[1]],split="\\t")[[1]][2]==pos&&
               checkMutation(ref,alt,strsplit(res_esp[[1]],split="\\t")[[1]][4],strsplit(res_esp[[1]],split="\\t")[[1]][5])){
                message("\t",length(res_esp[[1]])," matches in ESP6500")
                snps[i,10]<-strsplit(strsplit(strsplit(strsplit(res_esp[[1]],split="\\t")[[1]][8],split=";")[[1]][5],split="=")[[1]][2],split=",")[[1]][1]
            }
            if(length(res_1000g[[1]])==0||strsplit(res_1000g[[1]],split="\\t")[[1]][2]!=pos||
               !checkMutation(ref,alt,strsplit(res_1000g[[1]],split="\\t")[[1]][4],strsplit(res_1000g[[1]],split="\\t")[[1]][5])){
                snps[i,11]<-"NO"
            }
            if(length(res_1000g[[1]])>0&&strsplit(res_1000g[[1]],split="\\t")[[1]][2]==pos&&
               checkMutation(ref,alt,strsplit(res_1000g[[1]],split="\\t")[[1]][4],strsplit(res_1000g[[1]],split="\\t")[[1]][5])){
                message("\t",length(res_1000g[[1]])," matches in 1000G")
                snps[i,11]<-strsplit(strsplit(strsplit(res_1000g[[1]],split="\\t")[[1]][8],split=";")[[1]][9],split="=")[[1]][2]
            }
            flag<-T
            if(length(res_cosCC[[1]])==0&&length(res_cosNC[[1]])==0){
                snps[i,12]<-"NO"
                snps[i,13]<-"NO"
                snps[i,14]<-"O"
                flag<-F
            }
            snps[i,12]<-"NO"
            if(nchar(ref)==nchar(alt)){
                if(flag){
                    temp<-c()
                    if(length(res_cosCC[[1]])>0){
                        message("\t",length(res_cosCC[[1]])," matches in Cosmic (coding)")
                        for(j in 1:length(res_cosCC[[1]])){
                            if(!is.na(strsplit(strsplit(strsplit(res_cosCC[[1]],split="\\t")[[j]][8],split="=c\\.")[[1]][2],split=";")[[1]][1])&&
                               !is.na(strsplit(strsplit(as.character(snps[i,8]),split="\\|")[[1]][2],split="c\\.")[[1]][2])&&
                               strsplit(res_cosCC[[1]],split="\\t")[[j]][2]==pos&&
                               checkMutation(ref,alt,strsplit(res_cosCC[[1]],split="\\t")[[j]][4],strsplit(res_cosCC[[1]],split="\\t")[[j]][5])&&
                               strsplit(strsplit(strsplit(as.character(snps[i,8]),split="\\|")[[1]][2],split="c\\.")[[1]][2],split=",")[[1]][1]==strsplit(strsplit(strsplit(res_cosCC[[1]],split="\\t")[[j]][8],split="=c\\.")[[1]][2],split=";")[[1]][1]){
                                if(length(temp)>0){
                                    temp<-paste(temp,strsplit(res_cosCC[[1]],split="\\t")[[j]][3],sep=", ")
                                }
                                if(length(temp)==0){
                                    temp<-strsplit(res_cosCC[[1]],split="\\t")[[j]][3]
                                }
                                snps[i,13]<-"NA"
                                snps[i,14]<-0
                            }
                        }
                    }
                    if(length(res_cosNC[[1]])>0){
                        message("\t",length(res_cosNC[[1]])," matches in Cosmic (noncoding)")
                        for(j in 1:length(res_cosNC[[1]])){
                            if(!is.na(strsplit(strsplit(strsplit(res_cosNC[[1]],split="\\t")[[j]][8],split="=c\\.")[[1]][2],split=";")[[1]][1])&&
                               !is.na(strsplit(strsplit(as.character(snps[i,8]),split="\\|")[[1]][2],split="c\\.")[[1]][2])&&
                               strsplit(res_cosNC[[1]],split="\\t")[[j]][2]==pos&&
                               checkMutation(ref,alt,strsplit(res_cosNC[[1]],split="\\t")[[j]][4],strsplit(res_cosNC[[1]],split="\\t")[[j]][5])&&
                               strsplit(strsplit(strsplit(as.character(snps[i,8]),split="\\|")[[1]][2],split="c\\.")[[1]][2],split=",")[[1]][1]==strsplit(strsplit(strsplit(res_cosNC[[1]],split="\\t")[[j]][8],split="=c\\.")[[1]][2],split=";")[[1]][1]){
                                if(length(temp)>0){
                                    temp<-paste(temp,strsplit(res_cosNC[[1]],split="\\t")[[j]][3],sep=", ")
                                }
                                if(length(temp)==0){
                                    temp<-strsplit(res_cosNC[[1]],split="\\t")[[j]][3]
                                }
                                snps[i,13]<-"NA"
                                snps[i,14]<-0
                            }
                        }
                    }
                    if(length(temp)==0){
                        snps[i,13]<-"NA"
                        snps[i,14]<-"NA"
                    }
                    if(length(temp)==1){
                        snps[i,12]<-temp
                    }
                    if(length(temp)>1){
                        snps[i,12]<-temp[1]
                        for(j in 1:length(temp)){
                            if(vcountPattern(temp[j],snps[i,12])==0){
                                snps[i,12]<-paste(snps[i,12],temp[j],sep=", ")
                            }
                        }
                    }
                    ###advanced information
                    if(length(temp)>0){
                        output<-c()
                        for(j in 1:(vcountPattern(",",snps[i,12])+1)){
                            id<-strsplit(snps[i,12],split=", ")[[1]][j]
                            output<-c(output,system(paste("zcat ",dir3,"CosmicCompleteExport.bed.gz"," | grep -P \"",id,"\\t\"",sep=""),intern=T))
                            output<-c(output,system(paste("zcat ",dir3,"CosmicCompleteExport.fail.bed.gz"," | grep -P \"",id,"\\t\"",sep=""),intern=T))
                        }
                        if(length(output)==0){
                            snps[i,13]<-"Something is odd..."
                        }
                        if(length(output)==1){
                            snps[i,13]<-strsplit(output[[1]],split="\\t")[[1]][8]
                            snps[i,14]<-1
                        }
                        if(length(output)>1){
                            snps[i,13]<-strsplit(output[[1]],split="\\t")[[1]][8]
                            for(j in 1:length(output)){
                                if(vcountPattern(strsplit(output[[j]],split="\\t")[[1]][8],snps[i,13])==0){
                                    snps[i,13]<-paste(snps[i,13],strsplit(output[[j]],split="\\t")[[1]][8],sep=", ")
                                }
                            }
                            for(j in 1:(vcountPattern(",",snps[i,13])+1)){
                                if(j==1) {
                                    snps[i,14]<-sum(vcountPattern(strsplit(snps[i,13],", ")[[1]][j],output))
                                }
                                if(j>1){
                                    snps[i,14]<-paste(snps[i,14],sum(vcountPattern(strsplit(snps[i,13],", ")[[1]][j],output)),sep=", ")
                                }
                            }
                        }
                    }
                    
                }
            }
            if(nchar(ref)!=nchar(alt)){
                if(flag){
                    temp<-c()
                    no_aa<-F
                    if(length(res_cosCC[[1]])>0){
                        message("\t",length(res_cosCC[[1]])," matches in Cosmic (coding)")
                        for(j in 1:length(res_cosCC[[1]])){
                            if(!is.na(strsplit(strsplit(strsplit(res_cosCC[[1]],split="\\t")[[j]][8],split="=c\\.")[[1]][2],split=";")[[1]][1])&&
                               !is.na(strsplit(strsplit(as.character(snps[i,8]),split="\\|")[[1]][2],split="c\\.")[[1]][2])&&
                               strsplit(res_cosCC[[1]],split="\\t")[[j]][2]==pos&&
                               checkMutation(ref,alt,strsplit(res_cosCC[[1]],split="\\t")[[j]][4],strsplit(res_cosCC[[1]],split="\\t")[[j]][5])&&
                               strsplit(strsplit(as.character(snps[i,8]),split="\\|")[[1]][2],split="c\\.")[[1]][2]==strsplit(strsplit(strsplit(res_cosCC[[1]],split="\\t")[[j]][8],split="=c\\.")[[1]][2],split=";")[[1]][1]){
                                if(length(temp)>0){
                                    temp<-paste(temp,strsplit(res_cosCC[[1]],split="\\t")[[j]][3],sep=", ")
                                }
                                if(length(temp)==0){
                                    temp<-strsplit(res_cosCC[[1]],split="\\t")[[j]][3]
                                }
                                snps[i,13]<-"NA"
                                snps[i,14]<-0
                            }
                        }
                    }
                    if(length(res_cosNC[[1]])>0){
                        message("\t",length(res_cosNC[[1]])," matches in Cosmic (noncoding)")
                        for(j in 1:length(res_cosNC[[1]])){
                            if(!is.na(strsplit(strsplit(strsplit(res_cosNC[[1]],split="\\t")[[j]][8],split="=c\\.")[[1]][2],split=";")[[1]][1])&&
                               !is.na(strsplit(strsplit(as.character(snps[i,8]),split="\\|")[[1]][2],split="c\\.")[[1]][2])&&
                               strsplit(res_cosNC[[1]],split="\\t")[[j]][2]==pos&&
                               checkMutation(ref,alt,strsplit(res_cosNC[[1]],split="\\t")[[j]][4],strsplit(res_cosNC[[1]],split="\\t")[[j]][5])&&
                               strsplit(strsplit(as.character(snps[i,8]),split="\\|")[[1]][2],split="c\\.")[[1]][2]==strsplit(strsplit(strsplit(res_cosNC[[1]],split="\\t")[[j]][8],split="=c\\.")[[1]][2],split=";")[[1]][1]){
                                if(length(temp)>0){
                                    temp<-paste(temp,strsplit(res_cosNC[[1]],split="\\t")[[j]][3],sep=", ")
                                }
                                if(length(temp)==0){
                                    temp<-strsplit(res_cosNC[[1]],split="\\t")[[j]][3]
                                }
                                snps[i,13]<-"NA"
                                snps[i,14]<-0
                            }
                        }
                    }
                    if(length(temp)==0){
                        snps[i,13]<-"NA"
                        snps[i,14]<-"NA"
                    }
                    if(length(temp)==1){
                        snps[i,12]<-temp
                    }
                    if(length(temp)>1){
                        snps[i,12]<-temp[1]
                        for(j in 1:length(temp)){
                            if(vcountPattern(temp[j],snps[i,12])==0){
                                snps[i,12]<-paste(snps[i,12],temp[j],sep=", ")
                            }
                        }
                    }
                    ###advanced information
                    if(length(temp)>0){
                        output<-c()
                        for(j in 1:(vcountPattern(",",snps[i,12])+1)){
                            id<-strsplit(snps[i,12],split=", ")[[1]][j]
                            output<-c(output,system(paste("zcat ",dir3,"CosmicCompleteExport.bed.gz"," | grep -P \"",id,"\\t\"",sep=""),intern=T))
                            output<-c(output,system(paste("zcat ",dir3,"CosmicCompleteExport.fail.bed.gz"," | grep -P \"",id,"\\t\"",sep=""),intern=T))
                        }
                        if(length(output)==0){
                            snps[i,13]<-"Something is odd..."
                        }
                        if(length(output)==1){
                            snps[i,13]<-strsplit(output[[1]],split="\\t")[[1]][8]
                            snps[i,14]<-1
                        }
                        if(length(output)>1){
                            snps[i,13]<-strsplit(output[[1]],split="\\t")[[1]][8]
                            for(j in 1:length(output)){
                                if(vcountPattern(strsplit(output[[j]],split="\\t")[[1]][8],snps[i,13])==0){
                                    snps[i,13]<-paste(snps[i,13],strsplit(output[[j]],split="\\t")[[1]][8],sep=", ")
                                }
                            }
                            for(j in 1:(vcountPattern(",",snps[i,13])+1)){
                                if(j==1) {
                                    snps[i,14]<-sum(vcountPattern(strsplit(snps[i,13],", ")[[1]][j],output))
                                }
                                if(j>1){
                                    snps[i,14]<-paste(snps[i,14],sum(vcountPattern(strsplit(snps[i,13],", ")[[1]][j],output)),sep=", ")
                                }
                            }
                        }
                    }
                    if(no_aa){
                        snps[i,12]<-paste("???",snps[i,12],sep=", ")
                    }
                    
                }
            }
            
            snps[i,15]<-"NO"
            if(length(res_snp[[1]])>0){
              for(n in 1:length(res_snp[[1]])){
                if(strsplit(res_snp[[1]][n],split="\\t")[[1]][2]==pos&&
                   checkMutation(ref,alt,strsplit(res_snp[[1]][n],split="\\t")[[1]][4],strsplit(res_snp[[1]][n],split="\\t")[[1]][5])){
                  message("\t",length(res_snp[[1]])," matches in dbSNP (SNPs)")
                  snps[i,15]<-strsplit(res_snp[[1]][n],split="\\t")[[1]][3]
                }
              }
            }
            snps[i,16]<-"NO"
            snps[i,17]<-"NO"
            if(length(res_snv[[1]])>0){
              for(n in 1:length(res_snv[[1]])){
                if(strsplit(res_snv[[1]][n],split="\\t")[[1]][2]==pos&&
                   checkMutation(ref,alt,strsplit(res_snv[[1]][n],split="\\t")[[1]][4],strsplit(res_snv[[1]][n],split="\\t")[[1]][5])){
                  message("\t",length(res_snv[[1]])," matches in dbSNP (SNVs)")
                  snps[i,16]<-strsplit(res_snv[[1]][n],split="\\t")[[1]][3]
                  snps[i,17]<-"NO"
                  if(sum(strsplit(strsplit(res_snv[[1]][n],split="\\t")[[1]][8],split=";")[[1]]=="PM")>0){
                    snps[i,17]<-"YES"
                  }
                }
              }
            }
            if(length(res_exac[[1]])==0||strsplit(res_exac[[1]],split="\\t")[[1]][2]!=pos||
               !checkMutation(ref,alt,strsplit(res_exac[[1]],split="\\t")[[1]][4],strsplit(res_exac[[1]],split="\\t")[[1]][5])){
                snps[i,18]<-"NO"
            }
            if(length(res_exac[[1]])>0&&strsplit(res_exac[[1]],split="\\t")[[1]][2]==pos&&
               checkMutation(ref,alt,strsplit(res_exac[[1]],split="\\t")[[1]][4],strsplit(res_exac[[1]],split="\\t")[[1]][5])){
                message("\t",length(res_exac[[1]])," matches in ExAC")
                if(vcountPattern(pattern="AF",strsplit(strsplit(res_exac[[1]],split="\\t")[[1]][8],split=";")[[1]][12])>0){
                    if(vcountPattern(pattern=",",strsplit(strsplit(res_exac[[1]],split="\\t")[[1]][8],split=";")[[1]][12])==0){
                        snps[i,18]<-as.numeric(strsplit(strsplit(strsplit(res_exac[[1]],split="\\t")[[1]][8],split=";")[[1]][12],split="=")[[1]][2])
                    }
                    if(vcountPattern(pattern=",",strsplit(strsplit(res_exac[[1]],split="\\t")[[1]][8],split=";")[[1]][12])>0){
                        snps[i,18]<-max(as.numeric(strsplit(strsplit(strsplit(strsplit(res_exac[[1]],split="\\t")[[1]][8],split=";")[[1]][12],split="=")[[1]][2],split=",")[[1]]))
                    }
                }
                if(vcountPattern(pattern="AF",strsplit(strsplit(res_exac[[1]],split="\\t")[[1]][8],split=";")[[1]][13])>0){
                    if(vcountPattern(pattern=",",strsplit(strsplit(res_exac[[1]],split="\\t")[[1]][8],split=";")[[1]][13])==0){
                        snps[i,18]<-as.numeric(strsplit(strsplit(strsplit(res_exac[[1]],split="\\t")[[1]][8],split=";")[[1]][13],split="=")[[1]][2])
                    }
                    if(vcountPattern(pattern=",",strsplit(strsplit(res_exac[[1]],split="\\t")[[1]][8],split=";")[[1]][13])>0){
                        snps[i,18]<-max(as.numeric(strsplit(strsplit(strsplit(strsplit(res_exac[[1]],split="\\t")[[1]][8],split=";")[[1]][13],split="=")[[1]][2],split=",")[[1]]))
                    }
                }
            }
            if(as.character(snps[i,1])=="X"||as.character(snps[i,1])=="Y"||
               length(res_cvR[[1]])==0||strsplit(res_cvR[[1]],split="\\t")[[1]][2]!=pos||
               !checkMutation(ref,alt,strsplit(res_cvR[[1]],split="\\t")[[1]][4],strsplit(res_cvR[[1]],split="\\t")[[1]][5])){
                snps[i,19]<-"NO"
            }
            if(as.character(snps[i,1])!="X"&&as.character(snps[i,1])!="Y"&&
               length(res_cvR[[1]])>0&&strsplit(res_cvR[[1]],split="\\t")[[1]][2]==pos&&
               checkMutation(ref,alt,strsplit(res_cvR[[1]],split="\\t")[[1]][4],strsplit(res_cvR[[1]],split="\\t")[[1]][5])){
                message("\t",length(res_exac[[1]])," matches in ClinVar (relevant)")
                snps[i,19]<-"YES"
            }
            if(as.character(snps[i,1])=="X"||as.character(snps[i,1])=="Y"||
               length(res_cvN[[1]])==0||strsplit(res_cvN[[1]],split="\\t")[[1]][2]!=pos||
               !checkMutation(ref,alt,strsplit(res_cvN[[1]],split="\\t")[[1]][4],strsplit(res_cvN[[1]],split="\\t")[[1]][5])){
                snps[i,20]<-"NO"
            }
            if(as.character(snps[i,1])!="X"&&as.character(snps[i,1])!="Y"&&
               length(res_cvN[[1]])>0&&strsplit(res_cvN[[1]],split="\\t")[[1]][2]==pos&&
               checkMutation(ref,alt,strsplit(res_cvN[[1]],split="\\t")[[1]][4],strsplit(res_cvN[[1]],split="\\t")[[1]][5])){
                message("\t",length(res_exac[[1]])," matches in ClinVar (no impact)")
                snps[i,20]<-"YES"
            }
            
        }
    }
    firstVariant<-T
}

names(snps)<-c("chr","pos","ID","ref","alt","sample","gene","mutation","ENST",
               "ESP6500","1000G","Cosmic","Cosmic_location","Cosmic_occurence",
               "dbSNP (SNPs)","dbSNP (SNVs)","Clinically_flagged_dbSNP","ExAC",
               "ClinVar (clinical)","ClinVar (no impact)")

write.table(snps,paste(dir2,"Databases.txt",sep=""),row.names=F,sep="\t",quote=F) 

snps_pre[,23:33]<-snps[,10:20]
write.table(snps_pre,paste(dir2,"results_filtered_V1.txt",sep=""),row.names=F,sep="\t",quote=F)


