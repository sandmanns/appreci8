#########################################################################################################################
# Required genome / database resources:											                                        #
#															                                                            #
# 1. db_ClinVar_cc	     -> ClinVar common and clinical ($db_ClinVar_cc\.vcf.gz and tbi will be imported)		        #
# 2. db_ClinVar_nkmi     -> ClinVar no known medical impact ($db_ClinVar_nkmi\.vcf.gz and tbi will be imported)   	    #
# 3. db_Cosmic		     -> Cosmic Complete Export ($db_Cosmic\.bed.gz and tbi will be imported)			            #
# 4. db_CosmicCoding	 -> Cosmic Coding Mutations ($db_CosmicCoding\.vcf.gz and tbi will be imported)		            #
# 5. db_CosmicNoncoding  -> Cosmic Non Coding Variants ($db_CosmicNoncoding\.vcf.gz and tbi will be imported)	        #
# 6. db_dbSNP		     -> dbSNP ($db_dbSNP\.vcf.gz and tbi will be imported)					                        #
# 7. db_dbSNP_poly	     -> dbSNP polymorphisms only ($db_dbSNP_poly\.vcf.gz and tbi will be imported)		            #
# 8. db_ESP		         -> ESP6500 ($db_ESP\.chr[1-22|X|Y].snps_indels.vcf.gz and tbi will be imported)	          	#
# 9. db_ExAC		     -> ExAC ($ExAC\.vcf.gz and tbi will be imported)						                        #
#10. db_G1000		     -> 1000 Genomes (ALL.chr[1-22].$db_G1000\.vcf.gz and tbi will be imported)			            #
#11. db_G1000_X		     -> 1000 Genomes chromosome X ($db_G1000_X\.vcf.gz and tbi will be imported			            #
#12. db_G1000_Y		     -> 1000 Genomes chromosome Y ($db_G1000_Y\.vcf.gz and tbi will be imported			            #
#13. dbsnp1		         -> dbSNP data (vcf)									                                        #
#14. dbsnp2		         -> dbSNP data regarding polymorphisms, i.e. excluding sites after version 129 (vcf)	        #
#15. dir_databases	     -> Directory containing information on the databases					                        #
#16. genome		         -> Human genome to be used (fasta)								                                #
#17. knownindels1	     -> Gold standard set of indels (vcf)							                                #
#18. knownindels2	     -> Gold standard set of indels (vcf)							                                #
#19. peptides		     -> AA sequence for all peptides (fa)							                                #
#########################################################################################################################

db_ClinVar_cc="clinvar_20190916"
db_ClinVar_nkmi="common_no_known_medical_impact_20160203"
db_Cosmic="CosmicCompleteExport"
db_CosmicCoding="CosmicCodingMuts"
db_CosmicNoncoding="CosmicNonCodingVariants"
db_dbSNP="dbsnp_b151_GRCh37_all"
db_dbSNP_poly="dbsnp_b151_GRCh37_common"
db_ESP="ESP6500SI-V2-SSA137.GRCh38-liftover"
db_ExAC="ExAC.r0.3.sites.vep"
db_G1000="phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes"
db_G1000_X="ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes"
db_G1000_Y="ALL.chrY.phase3_integrated_v2a.20130502.genotypes"
dbsnp1=/        #path to dbsnp_138.b37.vcf
dbsnp2=/        #path to dbsnp_138.b37.excluding_sites_after_129.vcf
dir_databases=/ #path to database directory 
genome=/        #path to Homo_sapiens.GRCh37.67.dna.chromosome.all.fasta
genome2=/       #path to Homo_sapiens.GRCh37.67.dna.chr.all.fasta
knownindels1=/  #path to Mills_and_1000G_gold_standard.indels.b37.vcf
knownindels2=/  #path to 1000G_phase1.indels.b37.vcf
peptides=/      #path to Homo_sapiens.GRCh37.67.pep.all.fa



#########################################################################################################################
# Required tools:													                                                    #
#########################################################################################################################

bam_readcount=/ #path to bam-readcount
bcftools=/      #path to bcftools
freebayes=/     #path to freebayes
gatk=/          #path to gatk-package-4.0.4.0-local.jar
java_path=/     #path to java7
java_path8=/    #path to java8
lofreq=/        #path to lofreq
picard=/        #path to picard-tools-1.118
platypus=/      #path to Platypus.py
python=/        #path to python2.7 (old version required for Platypus.py)
samtools=/      #path to samtools
snver=/         #path to SNVerIndividual.jar
snpeff=/        #path to snpEff.jar
vardict=/       #path to vardict
varscan=/       #path to VarScan.v2.3.9.jar



#########################################################################################################################
# Additional requirements:												                                                #
#															                                                            #
# 1. $dir/vcf_header.txt     -> containing a typical vcf header						                                    #
#                                                                                                                       #
# 2. $dir/snpEff_ann/Hotspots.txt     -> list of hotspot mutations containing gene name, mutation (on protein level)	#
#					                     and (if available) minimum, expected VAF					                    #
#                                        if no hotspots are known, provide an empty txt file                            #
# 3. optionally $dir/documents/Primer_positions.bed   -> bed file with known primer positions to consider when          #
#                                                        calculating strand bias                                        #
#########################################################################################################################
