scripts="$( dirname "$0" )"
source "$scripts/Config.sh"

dir=$1
dir_bam=$2
sample=$3

mkdir -p $dir/freebayes

##define the minimum VAF as -F (default for appreci8 set to 0.01)
##define the minimum Nr_Alt as -C (default for appreci8 set to 20)

echo "Processing sample ${sample} with FreeBayes"
if [ ! -f $dir/freebayes/${sample}\.vcf ] ; then
  $freebayes -F 0.01 -C 20 -f $genome $dir_bam/${sample}\.bam > $dir/freebayes/${sample}\.vcf
  if [ $? -ne 0 ]; then
    echo "Error FreeBayes"
    rm $dir/freebayes/${sample}\.vcf
    exit 1
  fi
fi

