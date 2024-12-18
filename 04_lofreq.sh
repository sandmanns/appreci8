scripts="$( dirname "$0" )"
source "$scripts/Config.sh"

dir=$1
dir_bam=$2
sample=$3

mkdir -p $dir/lofreq/indels

echo "Processing sample ${sample} with LoFreq"
if [ ! -f $dir/lofreq/indels/${sample}\.vcf ] ; then
  $lofreq call --call-indels -f $genome -o $dir/lofreq/indels/${sample}\.vcf $dir_bam/${sample}\.bam
  if [ $? -ne 0 ]; then
    echo "Error LoFreq"
    rm $dir/lofreq/indels/${sample}\.vcf
    exit 1
  fi
fi



