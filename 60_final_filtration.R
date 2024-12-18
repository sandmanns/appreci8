library(seqinr)
library(Biostrings)
library(stringr)
args = commandArgs()
args = args[length(args)]
args = strsplit(args, ",")[[1]]
dir2<-paste(args[1],"/",sep="")

nr_alt<-as.numeric(args[2])
dp<-as.numeric(args[3])
vaf<-as.numeric(args[4])
low_bq<-as.numeric(args[5])
bq_diff<-as.numeric(args[6])
nrsamples<-as.numeric(args[7])
limit_provean<-as.numeric(args[8])
limit_provean2<-as.numeric(args[9])
primerpositions<-args[10]

dir1<-paste(dir2,"/documents/",sep="")
#prepare variants
input_temp<-read.table(paste(dir1,"results_filtered_V1.txt",sep=""),header=T,stringsAsFactors =F,sep="\t",quote = "")
provean<-read.table(paste(dir1,"results_provean.txt",sep=""),header=T,stringsAsFactors =F,sep="\t")
input<-cbind(input_temp[,1:33],provean[,3:4],input_temp[,36:45])
if(primerpositions){
    primer<-read.table(paste(dir1,"Primer_positions.bed",sep=""),header=F,stringsAsFactors =F,sep="\t")   
}
input2<-input[order(input[,2],input[,3],input[,4],input[,5],input[,1]),]
input<-input2

original<-read.table(paste(dir1,"raw_calls_counted.txt",sep=""),header=F,stringsAsFactors =F,sep="\t",quote = "")
while(length(grep("^ ",original[,1]))>0){
    original[,1]<-gsub("^ ","",original[,1])
}
temp<-str_split_fixed(original[,1]," ",n=Inf)
original[,1]<-temp[,2]
original<-cbind(original,Counts=temp[,1])
original[,5]<-as.numeric(as.character(original[,5]))

nrsamples_high<-ceiling(length(unique(input[,1]))/2)
if(nrsamples_high==1){
    nrsamples_high<-2
}
exclude<-rep(0,length(input[,1]))
artifact_because<-data.frame(freq=rep(NA,length(input[,1])),bq=rep(NA,length(input[,1])),alignment=rep(NA,length(input[,1])),
                             nr_samples=rep(NA,length(input[,1])),nr_samples_similar=rep(NA,length(input[,1])),
                             nr_databases=rep(NA,length(input[,1])),polymorphism_db=rep(NA,length(input[,1])),
                             mutation_db=rep(NA,length(input[,1])),cosmic_nr=rep(NA,length(input[,1])),
                             Poly_freq=rep(NA,length(input[,1])))
message("Calls: ",length(input[,1]))

temp<-cbind(as.numeric(input[,10]<low_bq)>(dp-nr_alt),input[,36]<=low_bq)
artifact_because[,2]<-rowSums(temp,na.rm = T)
artifact_because[,2]<-gsub("2","low bq",artifact_because[,2])
artifact_because[,2]<-gsub("1",NA,artifact_because[,2])
artifact_because[,2]<-gsub("0",NA,artifact_because[,2])

input[,4]<-gsub("TRUE","T",input[,4])

#nr of samples
message("Detect nr of samples with the same mutation")
original_1<-original[original[,1]=="1",]
original_2<-original[original[,1]=="2",]
original_3<-original[original[,1]=="3",]
original_4<-original[original[,1]=="4",]
original_5<-original[original[,1]=="5",]
original_6<-original[original[,1]=="6",]
original_7<-original[original[,1]=="7",]
original_8<-original[original[,1]=="8",]
original_9<-original[original[,1]=="9",]
original_10<-original[original[,1]=="10",]
original_11<-original[original[,1]=="11",]
original_12<-original[original[,1]=="12",]
original_13<-original[original[,1]=="13",]
original_14<-original[original[,1]=="14",]
original_15<-original[original[,1]=="15",]
original_16<-original[original[,1]=="16",]
original_17<-original[original[,1]=="17",]
original_18<-original[original[,1]=="18",]
original_19<-original[original[,1]=="19",]
original_20<-original[original[,1]=="20",]
original_21<-original[original[,1]=="21",]
original_22<-original[original[,1]=="22",]
original_x<-original[original[,1]=="X",]
original_y<-original[original[,1]=="Y",]

for(i in 1:length(input[,1])){
    if(i%%1000==0){
        message("i=",i," von ",length(input[,1]))
    }
    if(input[i,2]=="1"){
        temp<-original_1[grep(paste("^",input[i,3],"$",sep=""),original_1[,2]),]
    }
    if(input[i,2]=="2"){
        temp<-original_2[grep(paste("^",input[i,3],"$",sep=""),original_2[,2]),]
    }
    if(input[i,2]=="3"){
        temp<-original_3[grep(paste("^",input[i,3],"$",sep=""),original_3[,2]),]
    }
    if(input[i,2]=="4"){
        temp<-original_4[grep(paste("^",input[i,3],"$",sep=""),original_4[,2]),]
    }
    if(input[i,2]=="5"){
        temp<-original_5[grep(paste("^",input[i,3],"$",sep=""),original_5[,2]),]
    }
    if(input[i,2]=="6"){
        temp<-original_6[grep(paste("^",input[i,3],"$",sep=""),original_6[,2]),]
    }
    if(input[i,2]=="7"){
        temp<-original_7[grep(paste("^",input[i,3],"$",sep=""),original_7[,2]),]
    }
    if(input[i,2]=="8"){
        temp<-original_8[grep(paste("^",input[i,3],"$",sep=""),original_8[,2]),]
    }
    if(input[i,2]=="9"){
        temp<-original_9[grep(paste("^",input[i,3],"$",sep=""),original_9[,2]),]
    }
    if(input[i,2]=="10"){
        temp<-original_10[grep(paste("^",input[i,3],"$",sep=""),original_10[,2]),]
    }
    if(input[i,2]=="11"){
        temp<-original_11[grep(paste("^",input[i,3],"$",sep=""),original_11[,2]),]
    }
    if(input[i,2]=="12"){
        temp<-original_12[grep(paste("^",input[i,3],"$",sep=""),original_12[,2]),]
    }
    if(input[i,2]=="13"){
        temp<-original_13[grep(paste("^",input[i,3],"$",sep=""),original_13[,2]),]
    }
    if(input[i,2]=="14"){
        temp<-original_14[grep(paste("^",input[i,3],"$",sep=""),original_14[,2]),]
    }
    if(input[i,2]=="15"){
        temp<-original_15[grep(paste("^",input[i,3],"$",sep=""),original_15[,2]),]
    }
    if(input[i,2]=="16"){
        temp<-original_16[grep(paste("^",input[i,3],"$",sep=""),original_16[,2]),]
    }
    if(input[i,2]=="17"){
        temp<-original_17[grep(paste("^",input[i,3],"$",sep=""),original_17[,2]),]
    }
    if(input[i,2]=="18"){
        temp<-original_18[grep(paste("^",input[i,3],"$",sep=""),original_18[,2]),]
    }
    if(input[i,2]=="19"){
        temp<-original_19[grep(paste("^",input[i,3],"$",sep=""),original_19[,2]),]
    }
    if(input[i,2]=="20"){
        temp<-original_20[grep(paste("^",input[i,3],"$",sep=""),original_20[,2]),]
    }
    if(input[i,2]=="21"){
        temp<-original_21[grep(paste("^",input[i,3],"$",sep=""),original_21[,2]),]
    }
    if(input[i,2]=="22"){
        temp<-original_22[grep(paste("^",input[i,3],"$",sep=""),original_22[,2]),]
    }
    if(input[i,2]=="X"){
        temp<-original_x[grep(paste("^",input[i,3],"$",sep=""),original_x[,2]),]
    }
    if(input[i,2]=="Y"){
        temp<-original_y[grep(paste("^",input[i,3],"$",sep=""),original_y[,2]),]
    }
    
    artifact_because[i,4]<-temp[intersect(grep(paste("^",input[i,4],"$",sep=""),temp[,3]),
                                          grep(paste("^",input[i,5],"$",sep=""),temp[,4])),5]
    artifact_because[i,5]<-sum(temp[,5])
    
}


#databases
message("Detect nr of databases in which a mutation is present")
for(i in 1:length(input[,1])){
    #ESP
    if(!is.na(input[i,23])&&input[i,23]!="NO"){
        artifact_because[i,6]<-sum(artifact_because[i,6],1,na.rm=T)
        if(as.numeric(input[i,23])>0.03){
            artifact_because[i,7]<-sum(artifact_because[i,7],1,na.rm=T)   
        }
    }
    #1000G
    if(!is.na(input[i,24])&&input[i,24]!="NO"&&input[i,24]!="0"){
        artifact_because[i,6]<-sum(artifact_because[i,6],1,na.rm=T)
        if(vcountPattern(pattern=",",input[i,24])==0){
          if(as.numeric(input[i,24])>0.001){
            artifact_because[i,7]<-sum(artifact_because[i,7],1,na.rm=T)   
          }
        }
        if(vcountPattern(pattern=",",input[i,24])>0){
          if(max(as.numeric(strsplit(split=",",input[i,24])[[1]]))>0.001){
            artifact_because[i,7]<-sum(artifact_because[i,7],1,na.rm=T)   
          }
        }
    }
    #Cosmic
    if(!is.na(input[i,25])&&input[i,25]!="NO"){
        artifact_because[i,9]<-sum(as.numeric(strsplit(input[i,27],split=", ")[[1]]))
        artifact_because[i,6]<-sum(artifact_because[i,6],1,na.rm=T)
        if(artifact_because[i,9]>20){
            artifact_because[i,8]<-sum(artifact_because[i,8],1,na.rm=T)
        }
        if(artifact_because[i,9]<=20){
            artifact_because[i,7]<-sum(artifact_because[i,7],1,na.rm=T)
        }
    }
    #dbSNP
    if(!is.na(input[i,28])&&input[i,28]!="NO"){
        artifact_because[i,6]<-sum(artifact_because[i,6],1,na.rm=T)
        artifact_because[i,7]<-sum(artifact_because[i,7],1,na.rm=T)
    }
    if(!is.na(input[i,29])&&input[i,29]!="NO"){
        artifact_because[i,6]<-sum(artifact_because[i,6],1,na.rm=T)
        if(input[i,30]=="YES"){
            artifact_because[i,8]<-sum(artifact_because[i,8],1,na.rm=T)
        }
        if(input[i,30]=="NO"){
            #if(!is.na(input[i,28])&&input[i,28]!="NO"){
            #   artifact_because[i,7]<-sum(artifact_because[i,7],1,na.rm=T)
            #}
            if(!is.na(input[i,28])&&input[i,28]=="NO"){
                artifact_because[i,8]<-sum(artifact_because[i,8],1,na.rm=T)
            }
        }
    }
    #ExAC
    if(!is.na(input[i,31])&&input[i,31]!="NO"){
        artifact_because[i,6]<-sum(artifact_because[i,6],1,na.rm=T)
        if(as.numeric(input[i,31])>0.0005){
            artifact_because[i,7]<-sum(artifact_because[i,7],1,na.rm=T)   
        }
    }
    #Clinvar clinical
    if(!is.na(input[i,32])&&input[i,32]!="NO"){
        artifact_because[i,6]<-sum(artifact_because[i,6],1,na.rm=T)
        artifact_because[i,8]<-sum(artifact_because[i,8],1,na.rm=T)
    }
    #Clinvar common
    if(!is.na(input[i,33])&&input[i,33]!="NO"){
        artifact_because[i,6]<-sum(artifact_because[i,6],1,na.rm=T)
        artifact_because[i,7]<-sum(artifact_because[i,7],1,na.rm=T)
    }
}

#tolerated and freq
message("Check for freq when tolerated")
for(i in 1:length(input[,1])){
    if(!is.na(input[i,13])&&((input[i,13]>=0.35&&input[i,13]<=0.65)||(input[i,13]>=0.85))){
        artifact_because[i,10]<-1
    }
    if(!is.na(input[i,13])&&(input[i,13]<0.35||(input[i,13]>0.65&&input[i,13]<0.85))){
        artifact_because[i,10]<-0
    }
}

#large number of samples and high VAF
message("Check for high freq when large number of samples")
i<-1
while(i<=length(input[,1])){
    subset<-input[i:(i+artifact_because[i,4]-1),]
    if(length(subset[,1])>nrsamples&&sum(!is.na(subset[,13]))==length(subset[,1])&&
       sum(subset[,13]>0.85)>=floor(0.9*length(subset[,1]))){
        artifact_because[i:(i+artifact_because[i,4]-1),10]<-2
    }
    i<-i+artifact_because[i,4]
}

#test for strand bias
strandbias<-rep(NA,length(input[,1]))
for(i in 1:length(input[,1])){
    if(!is.na(input[i,38])&&!is.na(input[i,39])&&!is.na(input[i,42])&&!is.na(input[i,43])){
        test<-fisher.test(x=matrix(c(input[i,38],input[i,39],input[i,42],input[i,43]),ncol=2))
        strandbias[i]<-test$p.value
    }
    if(primerpositions==T){
        if(nchar(input[i,4])==1){
            chr<-as.character(input[i,2])==as.character(primer[,1])
            start<-input[i,3]>primer[,2]
            end<-input[i,3]<=primer[,3]
            if(sum(rowSums(cbind(chr,start,end))==3)>0){
                strandbias[i]<-2
            }      
        }
        if(nchar(input[i,4])>1){
            flag<-F
            chr<-as.character(input[i,2])==as.character(primer[,1])
            start<-matrix(rep(F,length(primer[,1])*nchar(input[i,4])),ncol=nchar(input[i,4]))
            end<-matrix(rep(F,length(primer[,1])*nchar(input[i,4])),ncol=nchar(input[i,4]))
            for(j in 1:nchar(input[i,4])){
                start[,j]<-(input[i,3]+j-1)>primer[,2]
                end[,j]<-(input[i,3]+j-1)<=primer[,3]
                if(flag==F&&sum(rowSums(cbind(chr,start[,j],end[,j]))==3)>0){
                    strandbias[i]<-2
                    flag=T
                }
            }
        }
    }
}

#check for hotspots
hotspots<-read.table(paste(dir2,"/snpEff_ann/Hotspots.txt",sep=""),header=T,stringsAsFactors = F)
hotspot<-rep(NA,length(input[,1]))
for(i in 1:length(input[,1])){
    if(input[i,6]!=""){
        gene<-strsplit(input[i,6],split=",")[[1]][1]
        if(vcountPattern("p.",input[i,9])>0){
            flag<-F
            ref<-strsplit(strsplit(input[i,9],split="p\\.")[[1]][2],split="/")[[1]][1]
            ref1<-a(substr(ref,1,3))
            ref2<-a(substr(ref,(nchar(ref)-2),nchar(ref)))
            if(is.na(ref2)){
                ref_new<-paste(ref1,substr(ref,4,nchar(ref)),sep="")
                if(vcountPattern(pattern="*",ref,fixed = T)>0){
                    ref_new2<-paste(ref1,substr(ref,4,nchar(ref)-1),sep="")
                }
                if(vcountPattern(pattern="*",ref,fixed = T)==0){
                    ref_new2<-paste(ref1,substr(ref,4,nchar(ref)-2),sep="")
                }
            }
            if(!is.na(ref2)){
                ref_new<-paste(ref1,substr(ref,4,nchar(ref)-3),ref2,sep="")
                ref_new2<-paste(ref1,substr(ref,4,nchar(ref)-3),sep="")
            }
            if(length(grep(ref_new,hotspots[,2],fixed=T))!=0&&
               length(intersect(grep(ref_new,hotspots[,2],fixed=T),grep(gene,hotspots[,1],fixed=T)))!=0){
                line_1<-grep(gene,hotspots[,1],fixed=T)
                line_2<-grep(ref_new,hotspots[,2],fixed=T)
                line<-intersect(line_1,line_2)
                if((!is.na(hotspots[line,3])&&!is.na(input[i,13])&&
                    input[i,13]>=hotspots[line,3])||is.na(hotspots[line,3])){
                    hotspot[i]<-1   
                    flag<-T
                }
            }
            if(flag==F&&length(grep(ref_new2,hotspots[,2],fixed=T))!=0&&
               length(intersect(grep(ref_new2,hotspots[,2],fixed=T),grep(gene,hotspots[,1],fixed=T)))!=0){
                line_1<-grep(gene,hotspots[,1],fixed=T)
                line_2<-grep(ref_new2,hotspots[,2],fixed=T)
                line<-intersect(line_1,line_2)
                if((!is.na(hotspots[line,3])&&!is.na(input[i,13])&&
                    input[i,13]>=hotspots[line,3])||is.na(hotspots[line,3])){
                    hotspot[i]<-1   
                }
            }
        }
    }
}


output<-cbind(artifact_because,input,strandbias)

##complex filtration
output4<-cbind(output,Category=NA)
for(i in 1:length(output4[,1])){
    artifact_score<-0
    if(output4[i,4]>nrsamples){
        artifact_score<-artifact_score+2
    }
    if(output4[i,4]>nrsamples_high&&is.na(hotspot[i])){
        artifact_score<-artifact_score+2
    }
    if((nchar(output4[i,14])>1||nchar(output4[i,15])>1)&&
       output4[i,5]>output4[i,4]){
        artifact_score<-artifact_score+1
    }
    if((nchar(output4[i,14])>1||nchar(output4[i,15])>1)&&
       !is.na(output4[i,23])&&output4[i,23]<0.05){
        artifact_score<-artifact_score+1
    }
    if(!is.na(output4[i,10])&&output4[i,10]==2){
        artifact_score<-artifact_score+2
    }
    if(primerpositions==F||(!is.na(output4[i,56])&&output4[i,56]!=2)){
        if(!is.na(output4[i,56])&&output4[i,56]<0.001){
            artifact_score<-artifact_score+1
        }
        if(!is.na(output4[i,49])&&!is.na(output4[i,53])&&
           (output4[i,49]>=(nr_alt/2)&&output4[i,53]>=(nr_alt/2))&&
           output4[i,56]<0.001){
            artifact_score<-artifact_score-1
        }
        if(!is.na(output4[i,48])&&!is.na(output4[i,49])&&
           (output4[i,49]<=2&&output4[i,48]>=(dp-nr_alt)/2)&&output4[i,56]>=0.001){
            artifact_score<-artifact_score+1
        }
        if(!is.na(output4[i,48])&&!is.na(output4[i,49])&&output4[i,56]<0.001&&
           (output4[i,49]<=2&&output4[i,48]<(dp-nr_alt)/2)){
            artifact_score<-artifact_score-1
        }
        if(!is.na(output4[i,52])&&!is.na(output4[i,53])&&
           (output4[i,53]<=2&&output4[i,52]>=(dp-nr_alt)/2)&&output4[i,56]>=0.001){
            artifact_score<-artifact_score+1
        } 
        if(!is.na(output4[i,52])&&!is.na(output4[i,53])&&output4[i,56]<0.001&&
           (output4[i,53]<=2&&output4[i,52]<(dp-nr_alt)/2)){
            artifact_score<-artifact_score-1
        } 
    }
    if(primerpositions==T&&!is.na(output4[i,56])&&output4[i,56]==2){
        artifact_score<-artifact_score-1
    }
    if(!is.na(output4[i,23])&&output4[i,23]<0.02){
        artifact_score<-artifact_score+2
    }
    if(is.na(output4[i,6])&&!is.na(output4[i,23])&&output4[i,23]<0.10){
        artifact_score<-artifact_score+1
    }
    if(is.na(output4[i,6])&&output4[i,4]>nrsamples_high){
        artifact_score<-artifact_score+1
    }
    if((!is.na(output4[i,44])&&output4[i,44]<limit_provean)){
        artifact_score<-artifact_score-1
    }
    if(!is.na(output4[i,44])&&output4[i,44]>limit_provean2&&
       !is.na(output4[i,10])&&output4[i,10]==0){
        artifact_score<-artifact_score+1
    }
    if(output4[i,32]>=4){
        artifact_score<-artifact_score-1
    }
    if(output4[i,32]>=5){
        artifact_score<-artifact_score-1
    }
    if(output4[i,32]>=6){
        artifact_score<-artifact_score-1
    }
    if(output4[i,32]<=1){
        artifact_score<-artifact_score+1
    }
    if(!is.na(output4[i,2])&&output4[i,2]=="low bq"){
        artifact_score<-artifact_score+4
    }
    if(!is.na(output4[i,40])&&output4[i,40]=="YES"&&is.na(hotspot[i])){
        artifact_score<-artifact_score-1
    }
    if(!is.na(hotspot[i])){
        artifact_score<-artifact_score-3
    }
    if(!is.na(output4[i,27])&&!is.na(output4[i,28])&&!is.na(output4[i,31])&&output4[i,27]==1&&output4[i,28]==1&&output4[i,31]==1){
        artifact_score<-artifact_score-3
    }
    if(artifact_score>-1){
        output4[i,57]<-paste("Artifact (",artifact_score,")",sep="")
    }
    if(artifact_score<=-1&&is.na(hotspot[i])){
        output4[i,57]<-paste("Probably True (",artifact_score,")",sep="")
    }
    if(artifact_score<=-1&&!is.na(hotspot[i])){
        output4[i,57]<-paste("Hotspot (",artifact_score,")",sep="")
    }
    
    poly_score<-0
    cosmic_flag<-F
    if(!is.na(output4[i,4])&&output4[i,4]>nrsamples){
        poly_score<-poly_score+1
    }
    if(!is.na(output4[i,4])&&output4[i,4]==1){
        poly_score<-poly_score-1
    }
    if(!is.na(output4[i,7])&&output4[i,7]>=2){
        poly_score<-poly_score+1
    }
    if(!is.na(output4[i,7])&&output4[i,7]>=4){
        poly_score<-poly_score+1
    }
    if(!is.na(output4[i,8])&&output4[i,8]>=2){
        poly_score<-poly_score-1
    }
    if(is.na(output4[i,7])){
        poly_score<-poly_score-1
    }
    if(output4[i,32]>=6){
        poly_score<-poly_score+1
    }
    if(!is.na(output4[i,18])&&vcountPattern(pattern="inframe",output4[i,18])>0&&vcountPattern(pattern="stop_gained",output4[i,18],fixed = T)==0){
        poly_score<-poly_score+1
    }
    if(!is.na(output4[i,10])&&output4[i,10]==1){
        poly_score<-poly_score+1
    }
    if(!is.na(output4[i,45])&&output4[i,45]=="Tolerated"&&output4[i,44]>=limit_provean2){
        poly_score<-poly_score+1
    }
    if(!is.na(output4[i,44])&&output4[i,44]<=limit_provean||!is.na(output4[i,18])&&output4[i,18]=="stop_gained"){
        poly_score<-poly_score-1
    }
    if(!is.na(output4[i,40])&&output4[i,40]=="YES"){
        poly_score<-poly_score-2
    }
    if(!is.na(output4[i,9])&&output4[i,9]>100){
        cosmic_flag<-T
    }
    if(poly_score>=2&&cosmic_flag==T&&is.na(hotspot[i])){
        output4[i,57]<-paste(output4[i,57],"Likely Polymorphism",sep="")
    }
    if(is.na(hotspot[i])&&((poly_score>=2&&cosmic_flag==F)
                           ||poly_score>=3)){
        output4[i,57]<-"Polymorphism"
    }
    if((poly_score>=2&&cosmic_flag==T&&is.na(hotspot[i]))||
       (is.na(hotspot[i])&&((poly_score>=2&&cosmic_flag==F)||poly_score>=3))){
        if(!is.na(output4[i,23])&&output4[i,23]<=0.1){
            artifact_score<-artifact_score+5
            if(artifact_score>-1){
                output4[i,57]<-paste("Artifact (",artifact_score,")",sep="")
            }
            if(artifact_score<=-1&&is.na(hotspot[i])){
                output4[i,57]<-paste("Probably True (",artifact_score,")",sep="")
            }
            if(artifact_score<=-1&&!is.na(hotspot[i])){
                output4[i,57]<-paste("Hotspot (",artifact_score,")",sep="")
            }
        }
        if(!is.na(output4[i,23])&&output4[i,23]<=0.2){
            artifact_score<-artifact_score+2
            if(artifact_score>-1){
                output4[i,57]<-paste("Artifact (",artifact_score,")",sep="")
            }
            if(artifact_score<=-1&&is.na(hotspot[i])){
                output4[i,57]<-paste("Probably True (",artifact_score,")",sep="")
            }
            if(artifact_score<=-1&&!is.na(hotspot[i])){
                output4[i,57]<-paste("Hotspot (",artifact_score,")",sep="")
            }
        }
        if(!is.na(output4[i,18])&&vcountPattern(pattern="frameshift",output4[i,18])>0){
            artifact_score<-artifact_score+2
            if(artifact_score>-1){
                output4[i,57]<-paste("Artifact (",artifact_score,")",sep="")
            }
            if(artifact_score<=-1&&is.na(hotspot[i])){
                output4[i,57]<-paste("Probably True (",artifact_score,")",sep="")
            }
            if(artifact_score<=-1&&!is.na(hotspot[i])){
                output4[i,57]<-paste("Hotspot (",artifact_score,")",sep="")
            }
        }
    }
}


output4.1<-data.frame(output4)
write.table(output4.1[,c(11:56,4,5,57)],paste(dir1,"results_filtered_V2.txt",sep=""),row.names=F,sep="\t",quote=F)

output5<-output4.1[vcountPattern(output4.1[,57],pattern="Artifact")==0,]
output6<-output5[output5[,57]!="Polymorphism",]

output7<-output6[!(is.na(output6[,21]) & !is.na(output6[,20]) & !is.na(output6[,22]) & (output6[,22]-output6[,20])<nr_alt),]

write.table(output7[,c(11:56,4,5,57)],paste(dir1,"results_filtered_V3.txt",sep=""),row.names=F,sep="\t",quote=F)
