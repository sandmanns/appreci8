scripts="$( dirname "$0" )"
source "$scripts/Config.sh"

dir=$1
dir_bam=$2
sample=$3

mkdir -p $dir/samtools

echo "Processing sample ${sample} with SAMtools"
if [ ! -f $dir/samtools/${sample}\.bcf ] ; then
  samtools mpileup -q 1 -g -u -o $dir/samtools/${sample}\.bcf -f $genome $dir_bam/${sample}\.bam
  if [ $? -ne 0 ]; then
    echo "Error SAMtools: mpileup"
    rm $dir/samtools/${sample}\.bcf
    exit 1
  fi
fi

if [ ! -f $dir/samtools/${sample}\.vcf ] ; then
  bcftools call -vmO v -o $dir/samtools/${sample}\.vcf $dir/samtools/${sample}\.bcf
  if [ $? -ne 0 ]; then
    echo "Error SAMtools: call"
    rm $dir/samtools/${sample}\.vcf
    exit 2
  fi
fi


