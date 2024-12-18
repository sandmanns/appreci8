scripts="$( dirname "$0" )"
source "$scripts/Config.sh"

dir=$1
dir_bam=$2
sample=$3

mkdir -p $dir/Platypus

echo "Processing sample ${sample} with Platypus"
if [ ! -f $dir/Platypus/${sample}.vcf ] ; then
  python $platypus callVariants --bamFiles=$dir_bam/${sample}.bam --refFile $genome --output=$dir/Platypus/${sample}.vcf --filterDuplicates=0 --minFlank=0
  if [ $? -ne 0 ]; then
    echo "Error Platypus"
    rm $dir/Platypus/${sample}.vcf
    exit 1
  fi
fi


