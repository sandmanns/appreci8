scripts="$( dirname "$0" )"
source "$scripts/Config.sh"

dir=$1
dir_bam=$2
sample=$3
targetregion=$4

tempdir=$dir/gatk/bam
intermres=$dir/gatk/intermRes
mkdir -p $tempdir $intermres

echo "Processing sample ${sample} with GATK"
if [ ! -f ${tempdir}/${sample}_RecalData.csv ] ; then
  $java_path -jar $gatk -T BaseRecalibrator -I $dir_bam/${sample}.bam -R $genome --maximum_cycle_value 1500 --covariate ContextCovariate \
    --covariate CycleCovariate --covariate QualityScoreCovariate --covariate ReadGroupCovariate -knownSites $dbsnp1 -knownSites $knownindels1 \
    -knownSites $knownindels2 -nct 1 -o ${tempdir}/${sample}_RecalData.csv 
  if [ $? -ne 0 ]; then
    $java_path -jar $gatk -T BaseRecalibrator -I $dir_bam/${sample}.bam -R $genome --maximum_cycle_value 1500 --covariate ContextCovariate \
    --covariate CycleCovariate --covariate QualityScoreCovariate --covariate ReadGroupCovariate -knownSites $dbsnp1 -knownSites $knownindels1 \
    -knownSites $knownindels2 -nct 1 --defaultBaseQualities 0 -o ${tempdir}/${sample}_RecalData.csv 
    if [ $? -ne 0 ]; then
      echo "Error GATK: BaseRecalibrator"
      rm ${tempdir}/${sample}_RecalData.csv
      exit 1
    fi
  fi
fi

if [ ! -f ${tempdir}/${sample}.bam ] ; then
  $java_path -jar $gatk -T PrintReads -I $dir_bam/${sample}.bam --filter_bases_not_stored -R $genome -BQSR ${tempdir}/${sample}_RecalData.csv -o ${tempdir}/${sample}.bam 
  if [ $? -ne 0 ]; then
    $java_path -jar $gatk -T PrintReads -I $dir_bam/${sample}.bam --filter_bases_not_stored --defaultBaseQualities 0 -R $genome -BQSR ${tempdir}/${sample}_RecalData.csv -o ${tempdir}/${sample}.bam 
    if [ $? -ne 0 ]; then    
      echo "Error GATK: PrintReads"
      rm ${tempdir}/${sample}.bam
      exit 2
    fi
  fi
fi

if [ ! -f ${tempdir}/${sample}.bai ] ; then
  $java_path -jar $pic/BuildBamIndex.jar VALIDATION_STRINGENCY="LENIENT" INPUT=${tempdir}/${sample}.bam OUTPUT=${tempdir}/${sample}.bai 
  if [ $? -ne 0 ]; then
    echo "Error GATK: BuildBamIndex"
    rm ${tempdir}/${sample}.bai
    exit 3
  fi
fi

if [ ! -f ${intermres}/${sample}.rawMutations.vcf ] ; then
  downsample="1500"
  $java_path -jar $gatk -T HaplotypeCaller -R $genome -stand_call_conf 30.0 -stand_emit_conf 10.0 --dbsnp $dbsnp2 -L $targetregion --max_alternate_alleles 9 \
    --downsample_to_coverage $downsample -nct 1 -I ${tempdir}/${sample}.bam -o ${intermres}/${sample}.rawMutations.vcf 
  if [ $? -ne 0 ]; then
    echo "Error GATK: HaplotypeCaller"
    rm $OUTFILE
    exit 4
  fi
fi


