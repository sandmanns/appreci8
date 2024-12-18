scripts="$( dirname "$0" )"
source "$scripts/Config.sh"

dir=$1
dir_bam=$2
sample=$3
targetregion=$4

mkdir -p $dir/vardict

##define the minimum VAF as af_thr (default for appreci8 set to 0.01)
##define the minimum Nr_Alt as -r (default for appreci8 set to 20)
##define the minimum phred score for a base to be considered "good" as -q (default for appreci8 set to 30)
##define the maximum indel size as -I (default for appreci8 set to 10)

echo "Processing sample ${sample} with VarDict"
if [ ! -f $dir/vardict/${sample}\.vcf ] ; then
  af_thr="0.01"
  $vardict -C -G $genome -f $af_thr -q 30 -I 10 -r 20 -N ${sample} -b $dir_bam/${sample}.bam -h -c 1 -S 2 -E 3 -g 4 $targetregion > $dir/vardict/${sample}.vcf
  if [ $? -ne 0 ]; then
    echo "Error VarDict"
    rm $dir/vardict/${sample}\.vcf
    exit 1
  fi
fi

