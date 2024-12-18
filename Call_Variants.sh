#########################################################################################
# Analysis information:													                #
#															                            #
# 1. dir		    -> Directory for the analysis								        #
# 2. dir_bam		-> Directory containing the bam and bai files						#
# 3. samplenames	-> File containing the sample names, one entry per line (txt)		#	
# 4. targetregion	-> File containing the target region (bed)							#
#########################################################################################

dir=$1
dir_bam=$2

scripts="$( dirname "$0" )"
source "$scripts/Config.sh"
$scripts/00_intro.sh
starttime=$(date +%s)
mkdir -p $dir/documents

samplenames=$dir/SampleNames.txt
targetregion=$dir/targetRegions/targetRegions.bed

#0) Perform alignment: 
#   we suggest use of BWA mem


#1) Perform actual variant calling:
if [ -f $samplenames ];
then
	samples=$(cat $samplenames)
else
	samples=$samplenames
fi


for sample in ${samples[*]} ; do
  $scripts/01_gatk.sh $dir $dir_bam $sample $targetregion
  if [ $? -ne 0 ]; then
    exit 1
  fi
  $scripts/02_platypus.sh $dir $dir_bam $sample
  if [ $? -ne 0 ]; then
    exit 2
  fi
  $scripts/03_varscan.sh $dir $dir_bam $sample
  if [ $? -ne 0 ]; then
    exit 3
  fi
  $scripts/04_lofreq.sh $dir $dir_bam $sample
  if [ $? -ne 0 ]; then
    exit 4
  fi
  $scripts/05_freebayes.sh $dir $dir_bam $sample
  if [ $? -ne 0 ]; then
    exit 5
  fi
  $scripts/06_snver.sh $dir $dir_bam $sample
  if [ $? -ne 0 ]; then
    exit 6
  fi
  $scripts/07_samtools.sh $dir $dir_bam $sample
  if [ $? -ne 0 ]; then
    exit 7
  fi
  $scripts/08_vardict.sh $dir $dir_bam $sample $targetregion
  if [ $? -ne 0 ]; then
    exit 8
  fi
done

#2) Postprocessing on variant caller level
for sample in ${samples[*]} ; do
   echo "Sample $sample"
   mkdir $dir/gatk/intermRes/calling_comp/
   bedtools intersect -a $dir/gatk/intermRes/$sample.rawMutations.vcf -b $dir/targetRegions/targetRegions.bed > $dir/gatk/intermRes/calling_comp/$sample.target.vcf

   mkdir $dir/Platypus/calling_comp/
   bedtools intersect -a $dir/Platypus/$sample.vcf -b $dir/targetRegions/targetRegions.bed > $dir/Platypus/calling_comp/$sample.target.vcf

   mkdir $dir/varscan/calling_comp/
   tail -n +2 $dir/varscan/snvs/$sample.txt | awk '{FS="\t" ; {print $1"\t"$2"\t"".""\t"$3"\t"$4"\t""1""\t"".""\t""""\t""sample"}}' > $dir/varscan/temp1.txt
   tail -n +2 $dir/varscan/indels/$sample.txt | awk '{FS="\t" ; {print $1"\t"$2"\t"".""\t"$3"\t"$4"\t""1""\t"".""\t""""\t""sample"}}' > $dir/varscan/temp2.txt
   cat $dir/vcf_header.txt $dir/varscan/temp1.txt $dir/varscan/temp2.txt > $dir/varscan/temp.vcf 
   bedtools intersect -a $dir/varscan/temp.vcf -b $dir/targetRegions/targetRegions.bed > $dir/varscan/calling_comp/$sample.target.vcf

   mkdir $dir/lofreq/calling_comp/
   bedtools intersect -a $dir/lofreq/indels/$sample.vcf -b $dir/targetRegions/targetRegions.bed > $dir/lofreq/calling_comp/$sample.target.vcf

   mkdir $dir/freebayes/calling_comp/
   bedtools intersect -a $dir/freebayes/$sample.vcf -b $dir/targetRegions/targetRegions.bed > $dir/freebayes/calling_comp/$sample.target.vcf

   mkdir $dir/samtools/calling_comp/
   bedtools intersect -a $dir/samtools/$sample.vcf -b $dir/targetRegions/targetRegions.bed > $dir/samtools/calling_comp/$sample.target.vcf

   mkdir $dir/SNVer/calling_comp/
   grep -v "#" $dir/SNVer/$sample.vcf.filter.vcf > $dir/SNVer/temp1.txt
   grep -v "#" $dir/SNVer/$sample.vcf.indel.filter.vcf > $dir/SNVer/temp2.txt
   cat $dir/vcf_header.txt $dir/SNVer/temp1.txt $dir/SNVer/temp2.txt > $dir/SNVer/temp.vcf 
   bedtools intersect -a $dir/SNVer/temp.vcf -b $dir/targetRegions/targetRegions.bed > $dir/SNVer/calling_comp/$sample.target.vcf

   mkdir $dir/vardict/calling_comp/
   tail -n +2 $dir/vardict/$sample.vcf | awk '{FS="\t" ; {print $3"\t"$4"\t"".""\t"$6"\t"$7"\t""1""\t"".""\t""""\t""sample"}}' > $dir/vardict/temp.txt
   cat $dir/vcf_header.txt $dir/vardict/temp.txt > $dir/vardict/temp.vcf 
   bedtools intersect -a $dir/vardict/temp.vcf -b $dir/targetRegions/targetRegions.bed > $dir/vardict/calling_comp/$sample.target.vcf

# output: *.target.vcf
done


gatk=TRUE
platypus=TRUE
varscan=TRUE
lofreq=TRUE
freebayes=TRUE
snver=TRUE
samtools=TRUE
vardict=TRUE

R --vanilla <./speedup/11_normalize_calls.R $gatk,$platypus,$varscan,$lofreq,$freebayes,$snver,$samtools,$vardict,$dir
#output: *.caller_raw2.vcf


#3) Annotate variants
for sample in ${samples[*]} ; do
   echo "Sample $sample"
   echo "GATK"
   cat $dir/gatk/intermRes/calling_comp/$sample.gatk_raw2.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""1""\t"".""\t""""\t""sample"}}' > $dir/gatk/intermRes/calling_comp/temp.vcf
   tail -n +2 $dir/gatk/intermRes/calling_comp/temp.vcf > $dir/gatk/intermRes/calling_comp/temp2.vcf
   cat $dir/vcf_header.txt $dir/gatk/intermRes/calling_comp/temp2.vcf > $dir/gatk/intermRes/calling_comp/$sample.gatk_raw.vcf
   java -jar $snpeff eff -formatEff GRCh37.75 -t -noStats $dir/gatk/intermRes/calling_comp/$sample.gatk_raw.vcf > $dir/gatk/intermRes/calling_comp/$sample.gatk_ann.vcf

   echo "Platypus"
   cat $dir/Platypus/calling_comp/$sample.Platypus_raw2.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""1""\t"".""\t""""\t""sample"}}' > $dir/Platypus/calling_comp/temp.vcf
   tail -n +2 $dir/Platypus/calling_comp/temp.vcf > $dir/Platypus/calling_comp/temp2.vcf
   cat $dir/vcf_header.txt $dir/Platypus/calling_comp/temp2.vcf > $dir/Platypus/calling_comp/$sample.Platypus_raw.vcf
   java -jar $snpeff eff -formatEff GRCh37.75 -t -noStats $dir/Platypus/calling_comp/$sample.Platypus_raw.vcf > $dir/Platypus/calling_comp/$sample.Platypus_ann.vcf

   echo "VarScan"
   cat $dir/varscan/calling_comp/$sample.varscan_raw2.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""1""\t"".""\t""""\t""sample"}}' > $dir/varscan/calling_comp/temp.vcf
   tail -n +2 $dir/varscan/calling_comp/temp.vcf > $dir/varscan/calling_comp/temp2.vcf
   cat $dir/vcf_header.txt $dir/varscan/calling_comp/temp2.vcf > $dir/varscan/calling_comp/$sample.varscan_raw.vcf
   java -jar $snpeff eff -formatEff GRCh37.75 -t -noStats $dir/varscan/calling_comp/$sample.varscan_raw.vcf > $dir/varscan/calling_comp/$sample.varscan_ann.vcf

   echo "LoFreq"
   cat $dir/lofreq/calling_comp/$sample.lofreq_raw2.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""1""\t"".""\t""""\t""sample"}}' > $dir/lofreq/calling_comp/temp.vcf
   tail -n +2 $dir/lofreq/calling_comp/temp.vcf > $dir/lofreq/calling_comp/temp2.vcf
   cat $dir/vcf_header.txt $dir/lofreq/calling_comp/temp2.vcf > $dir/lofreq/calling_comp/$sample.lofreq_raw.vcf
   java -jar $snpeff eff -formatEff GRCh37.75 -t -noStats $dir/lofreq/calling_comp/$sample.lofreq_raw.vcf > $dir/lofreq/calling_comp/$sample.lofreq_ann.vcf

   echo "FreeBayes"
   cat $dir/freebayes/calling_comp/$sample.freebayes_raw2.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""1""\t"".""\t""""\t""sample"}}' > $dir/freebayes/calling_comp/temp.vcf
   tail -n +2 $dir/freebayes/calling_comp/temp.vcf > $dir/freebayes/calling_comp/temp2.vcf
   cat $dir/vcf_header.txt $dir/freebayes/calling_comp/temp2.vcf > $dir/freebayes/calling_comp/$sample.freebayes_raw.vcf
   java -jar $snpeff eff -formatEff GRCh37.75 -t -noStats $dir/freebayes/calling_comp/$sample.freebayes_raw.vcf > $dir/freebayes/calling_comp/$sample.freebayes_ann.vcf

   echo "SNVer"
   cat $dir/SNVer/calling_comp/$sample.snver_raw2.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""1""\t"".""\t""""\t""sample"}}' > $dir/SNVer/calling_comp/temp.vcf
   tail -n +2 $dir/SNVer/calling_comp/temp.vcf > $dir/SNVer/calling_comp/temp2.vcf
   cat $dir/vcf_header.txt $dir/SNVer/calling_comp/temp2.vcf > $dir/SNVer/calling_comp/$sample.snver_raw.vcf
   java -jar $snpeff eff -formatEff GRCh37.75 -t -noStats $dir/SNVer/calling_comp/$sample.snver_raw.vcf > $dir/SNVer/calling_comp/$sample.snver_ann.vcf

   echo "SAMtools"
   cat $dir/samtools/calling_comp/$sample.samtools_raw2.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""1""\t"".""\t""""\t""sample"}}' > $dir/samtools/calling_comp/temp.vcf
   tail -n +2 $dir/samtools/calling_comp/temp.vcf > $dir/samtools/calling_comp/temp2.vcf
   cat $dir/vcf_header.txt $dir/samtools/calling_comp/temp2.vcf > $dir/samtools/calling_comp/$sample.samtools_raw.vcf
   java -jar $snpeff eff -formatEff GRCh37.75 -t -noStats $dir/samtools/calling_comp/$sample.samtools_raw.vcf > $dir/samtools/calling_comp/$sample.samtools_ann.vcf

   echo "Vardict"
   cat $dir/vardict/calling_comp/$sample.vardict_raw2.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""1""\t"".""\t""""\t""sample"}}' > $dir/vardict/calling_comp/temp.vcf
   tail -n +2 $dir/vardict/calling_comp/temp.vcf > $dir/vardict/calling_comp/temp2.vcf
   cat $dir/vcf_header.txt $dir/vardict/calling_comp/temp2.vcf > $dir/vardict/calling_comp/$sample.vardict_raw.vcf
   java -jar $snpeff eff -formatEff GRCh37.75 -t -noStats $dir/vardict/calling_comp/$sample.vardict_raw.vcf > $dir/vardict/calling_comp/$sample.vardict_ann.vcf

#   output: *_ann.vcf
done

#4) Filter annotated files
gatk=TRUE
platypus=TRUE
varscan=TRUE
lofreq=TRUE
freebayes=TRUE
snver=TRUE
samtools=TRUE
vardict=TRUE

R --vanilla <./speedup/20_filter_annotation.R $gatk,$platypus,$varscan,$lofreq,$freebayes,$snver,$samtools,$vardict,$dir

mkdir $dir/documents/test

#5) Generate combined output
for sample in ${samples[*]} ; do
echo "Sample $sample"

  tail -n +2 $dir/gatk/intermRes/calling_comp/$sample.gatk_filtered.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$4"\t"$5"\t"$8"\t""gatk"}}' > $dir/documents/test/$sample\_gatk1.vcf
  tail -n +2 $dir/Platypus/calling_comp/$sample.Platypus_filtered.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$4"\t"$5"\t"$8"\t""platypus"}}' > $dir/documents/test/$sample\_platypus1.vcf
  tail -n +2 $dir/varscan/calling_comp/$sample.varscan_filtered.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$4"\t"$5"\t"$8"\t""varscan"}}' > $dir/documents/test/$sample\_varscan1.vcf
  tail -n +2 $dir/lofreq/calling_comp/$sample.lofreq_filtered.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$4"\t"$5"\t"$8"\t""lofreq"}}' > $dir/documents/test/$sample\_lofreq1.vcf
  tail -n +2 $dir/freebayes/calling_comp/$sample.freebayes_filtered.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$4"\t"$5"\t"$8"\t""freebayes"}}' > $dir/documents/test/$sample\_freebayes1.vcf
  tail -n +2 $dir/SNVer/calling_comp/$sample.snver_filtered.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$4"\t"$5"\t"$8"\t""snver"}}' > $dir/documents/test/$sample\_snver1.vcf
  tail -n +2 $dir/samtools/calling_comp/$sample.samtools_filtered.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$4"\t"$5"\t"$8"\t""samtools"}}' > $dir/documents/test/$sample\_samtools1.vcf
  tail -n +2 $dir/vardict/calling_comp/$sample.vardict_filtered.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$4"\t"$5"\t"$8"\t""vardict"}}' > $dir/documents/test/$sample\_vardict1.vcf
  
  cat $dir/documents/test/$sample\_gatk1.vcf $dir/documents/test/$sample\_platypus1.vcf $dir/documents/test/$sample\_varscan1.vcf $dir/documents/test/$sample\_lofreq1.vcf $dir/documents/test/$sample\_freebayes1.vcf $dir/documents/test/$sample\_snver1.vcf $dir/documents/test/$sample\_samtools1.vcf $dir/documents/test/$sample\_vardict1.vcf > $dir/documents/test/$sample\_temp.vcf
  sort $dir/documents/test/$sample\_temp.vcf > $dir/documents/test/$sample\_temp.sorted.vcf
  cat $dir/documents/test/$sample\_temp.sorted.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"}}' > $dir/documents/test/$sample\_temp2.vcf
  uniq -c $dir/documents/test/$sample\_temp2.vcf > $dir/documents/test/$sample\_temp.counted.vcf
  uniq $dir/documents/test/$sample\_temp2.vcf > $dir/documents/test/$sample\_temp.uniq.vcf

  cat $dir/documents/test/$sample\_gatk1.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"}}' > $dir/documents/test/$sample\_gatk.vcf
  cat $dir/documents/test/$sample\_platypus1.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"}}' > $dir/documents/test/$sample\_platypus.vcf
  cat $dir/documents/test/$sample\_varscan1.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"}}' > $dir/documents/test/$sample\_varscan.vcf
  cat $dir/documents/test/$sample\_lofreq1.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"}}' > $dir/documents/test/$sample\_lofreq.vcf
  cat $dir/documents/test/$sample\_freebayes1.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"}}' > $dir/documents/test/$sample\_freebayes.vcf
  cat $dir/documents/test/$sample\_snver1.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"}}' > $dir/documents/test/$sample\_snver.vcf
  cat $dir/documents/test/$sample\_samtools1.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"}}' > $dir/documents/test/$sample\_samtools.vcf
  cat $dir/documents/test/$sample\_vardict1.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"}}' > $dir/documents/test/$sample\_vardict.vcf

  rm $dir/documents/test/$sample\_temp.vcf $dir/documents/test/$sample\_temp.less.vcf $dir/documents/test/$sample\_gatk1.vcf $dir/documents/test/$sample\_platypus1.vcf
  rm $dir/documents/test/$sample\_varscan1.vcf $dir/documents/test/$sample\_lofreq1.vcf $dir/documents/test/$sample\_freebayes1.vcf $dir/documents/test/$sample\_snver1.vcf
  rm $dir/documents/test/$sample\_samtools1.vcf $dir/documents/test/$sample\_vardict1.vcf

  cat $dir/documents/test/$sample\_temp.uniq.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4"\t"}}' > $dir/documents/test/$sample\_temp.uniq2.vcf
  rm $dir/documents/test/$sample\_temp2.vcf 
# output: *.results_raw.txt

done


R --vanilla <./speedup/30_combine_per_sample.R $dir
# output: *.results_raw.txt


#6) Determine Characteristics of reported calls
R --vanilla <./speedup/40_determine_frequencies.R $dir,$bam_readcount,$genome,$dir_bam
# output: *.frequency.txt


#7) Combine and filter
min_alt=20
min_dp=50
min_vaf=0.01
min_bq_alt=15
max_bq_diff=7
R --vanilla <./speedup/50_combine_all_samples.R $dir,$min_alt,$min_dp,$min_vaf,$min_bq_alt,$max_bq_diff
# output:: *.results_filtered_V1_pre.txt


#8) Databases
R --vanilla <./speedup/51_database_lookup.R $dir,$dir_databases,$db_ClinVar_cc,$db_ClinVar_nkmi,$db_Cosmic,$db_CosmicCoding,$db_CosmicNoncoding,$db_dbSNP,$db_dbSNP_poly,$db_ESP,$db_ExAC,$db_G1000,$db_G1000_X,$db_G1000_Y
# output Databases.txt and results_filtered_V1.txt


mkdir $dir/snpEff_ann
R --vanilla <./speedup/52_effect_protein.R $dir,$peptides
# output: results_provean.txt


#9) Final filtration
echo "" > $dir/documents/raw_calls.txt
for sample in ${samples[*]} ; do
   tail -n +2 $dir/documents/$sample.frequency.vcf | awk '{FS="\t" ; {print $1"\t"$2"\t"$3"\t"$4}}' > $dir/documents/temp.txt
   cat $dir/documents/raw_calls.txt $dir/documents/temp.txt > $dir/documents/raw_calls_temp.txt
   mv $dir/documents/raw_calls_temp.txt $dir/documents/raw_calls.txt
done

tail -n +2 $dir/documents/raw_calls.txt > $dir/documents/raw_calls_temp.txt
mv $dir/documents/raw_calls_temp.txt $dir/documents/raw_calls.txt
sort $dir/documents/raw_calls.txt > $dir/documents/raw_calls_sorted.txt
uniq -c $dir/documents/raw_calls_sorted.txt > $dir/documents/raw_calls_counted.txt

max_samples=3
provean_del=-3
provean_tol=-1.5
primer="FALSE"
R --vanilla <./speedup/60_final_filtration.R $dir,$min_alt,$min_dp,$min_vaf,$min_bq_alt,$min_bq_diff,$max_samples,$provean_del,$provean_tol,$primer
# output: results_filtered_V2.txt and results_filtered_V3.txt


if [ -e $dir/documents/results_filtered_V2.txt ] && [ -e $dir/documents/results_filtered_V3.txt ]; then
  echo "Variant Calling successful!"
  endtime=$(date +%s)
  echo "Elapsed time: " $(($endtime - $starttime)) " seconds"
fi


