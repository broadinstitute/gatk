(cd .. && exec ./gradlew assemble localJar)

mkdir -p gvcf bpres

for bam in $(ls bam-snippets/*bam)
do 
    sn=$(basename $bam .bam)
    ../gatk HaplotypeCaller \
        -I $bam \
        -O gvcf/$sn.gvcf \
        -R ref/PlasmoDB-54_Pfalciparum3D7_Genome.fasta \
        -L Pf3D7_01_v3:100390-100430 \
        -ERC GVCF \
        -bamout gvcf/$sn.bam
    ../gatk HaplotypeCaller \
        -I $bam \
        -O bpres/$sn.gvcf \
        -R ref/PlasmoDB-54_Pfalciparum3D7_Genome.fasta \
        -L Pf3D7_01_v3:100390-100430 \
        -ERC BP_RESOLUTION \
        -bamout bpres/$sn.bam
done

for tag in gvcf bpres
do
    ls $tag/*gvcf > $tag/input.list
    ../gatk CombineGVCFs \
        -V $tag/input.list \
        -O $tag/combined.gvcf \
        -R ref/PlasmoDB-54_Pfalciparum3D7_Genome.fasta
    ../gatk GenotypeGVCFs \
        -V $tag/combined.gvcf \
        -O $tag/genotyped.vcf \
        -R ref/PlasmoDB-54_Pfalciparum3D7_Genome.fasta
done

diff gvcf/genotyped.vcf bpres/genotyped.vcf | grep -v '#'
