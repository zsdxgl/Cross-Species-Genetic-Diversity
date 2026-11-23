1. Rads mapping
ref="reference.fa"
path=""
i="samples"
T=10
bwa mem -t 30 -M -R @RG\tID:${i}\tSM:${i} ${path}/${ref} ${path}/${i}_1.fastq.gz ${path}/${i}_2.fastq.gz |samtools view -@ ${T} -q 30  -b  |samtools sort --threads ${T} -o ${path}/bam/${i}.bam --reference ${path}/${ref} 
sambamba markdup --overflow-list-size 600000 --tmpdir=${path}/bam/ -r ${path}/bam/${i}.bam ${path}/bam/${i}.D.bam


2. Orthologous genes identifing
busco -i anak.reference.fa -m genome --augustus -o busco.out -f  -c 10  -l actinopterygii_odb10
##output: busco.out/run_actinopterygii_odb10/full_table.tsv
          busco.out/run_actinopterygii_odb10/augustus_output/gff  

3. 4D sites identifing
Reseqtools/iTools_Code/iTools Gfftools getCdsPep -Ref $ref -Gff ${i}.gff -OutPut output -4DSite
 

4. Shared single-copy genes identifing
./get_shared_single_copy_gene_4D_sites.pl
###example input:
anak.chr.length 
anak.reference.fa
anak.full_table.tsv（from busco）
anak.reseqtools_out/anak.4Dsite.gz
###example output:
anak.intervals.txt; anak.4Dsites.txt 

5. SNP calling and 4D site heterozygosity estimating
bcftools mpileup -Ou -q 30 -Q 20 -f anak.reference.fa -R anak.intervals.txt --rf 0x2 --ff 0x100 --ff DUP --ff 0x800 -b bam.list  --threads 10 -a "FORMAT/AD,FORMAT/DP" | bcftools call --threads 10 -V indels  -cO z -o anak.vcf.gz

perl calc_heterozygosity_4D_sites.pl anak.vcf.gz anak.4Dsites.txt >anak.4Dsite.het

6. whole-genome site heterozygosity estimating
perl calc_heterozygosity_fast.pl  --vcf /data1/guoliang/ZhuoLJ/vcf_wholeSNP/trachinotus_ovatus.vcf.gz --out /data1/guoliang/ZhuoLJ/vcf_wholeSNP/trachinotus_ovatus.hewholeSNP.txt   --total-sites 100000 --threads 20


