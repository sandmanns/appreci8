scripts="$( dirname "$0" )"
source "$scripts/Config.sh"

dir=$1
dir_bam=$2
sample=$3

mkdir -p $dir/SNVer

##define the minimum VAF as -b (default for appreci8 set to 0.01)
##define the minimum base quality as -bq (default for appreci8 set to 30)
##define the minimum Nr_Alt per strand as -a (default for appreci8 set to 2)

echo "Processing sample ${sample} with SNVer"
if [ ! -f $dir/SNVer/${sample}\.vcf.filter.vcf ] && [ ! -f $dir/SNVer/${sample}\.vcf.indel.filter.vcf ] ; then
  $java_path -jar $snver -i $dir_bam/${sample}\.bam -r $genome -b 0.01 -bq 30 -a 2 -o $dir/SNVer/${sample}\.vcf
  if [ $? -ne 0 ]; then
    echo "Error SNVer"
    rm $dir/SNVer/${sample}\.vcf.indel.raw.vcf $dir/SNVer/${sample}\.vcf.raw.vcf $dir/SNVer/${sample}\.vcf.filter.vcf $dir/SNVer/${sample}\.vcf.indel.filter.vcf
    exit 1
  fi
  rm $dir/SNVer/${sample}\.vcf.failed.log 
fi

