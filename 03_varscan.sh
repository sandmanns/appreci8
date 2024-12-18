scripts="$( dirname "$0" )"
source "$scripts/Config.sh"

dir=$1
dir_bam=$2
sample=$3

mkdir -p $dir/varscan/bcf $dir/varscan/snvs $dir/varscan/indels

echo "Processing sample ${sample} with VarScan"
if [ ! -f $dir/varscan/bcf/${sample}\.bcf ] ; then
  samtools mpileup -f $genome $dir_bam/${sample}\.bam > $dir/varscan/bcf/${sample}\.bcf
  if [ $? -ne 0 ]; then
    echo "Error VarScan: mpileup"
    rm $dir/varscan/bcf/${sample}\.bcf
    exit 1
  fi
fi

if [ ! -f $dir/varscan/snvs/${sample}\.txt ] ; then
  $java_path -jar $varscan mpileup2snp $dir/varscan/bcf/${sample}\.bcf > $dir/varscan/snvs/${sample}\.txt
  if [ $? -ne 0 ]; then
    echo "Error VarScan: mpileup2snp"
    rm $dir/varscan/snvs/${sample}\.txt
    exit 2
  fi
fi

if [ ! -f $dir/varscan/indels/${sample}\.txt ] ; then
  $java_path -jar $varscan mpileup2indel $dir/varscan/bcf/${sample}\.bcf > $dir/varscan/indels/${sample}\.txt
  if [ $? -ne 0 ]; then
    echo "Error VarScan: mpileup2indel"
    rm $dir/varscan/indels/${sample}\.txt
    exit 3
  fi
fi




