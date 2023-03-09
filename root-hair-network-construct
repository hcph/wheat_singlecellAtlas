####1 TF Footprint analysis
##prepare files
#extrat peaks
perl -p -i -e 's/,/\t/g' *csv
perl -p -i -e 's/\"//g' *.csv
perl -p -i -e 's/^M//g' *csv
awk '{print $2"\t"$3"\t"$4"\t"$11"\t"$14"\t"$15}' ground.meristem.csv >ground.meristem.bed
awk '{print $2"\t"$3"\t"$4"\t"$11"\t"$13"\t"$14}' epidermis.csv >epidermis.bed
#extract barcodes
write.csv(projHeme5@cellColData,"cellColdata_barcode_celltype.csv")
perl -p -i -e 's/^M//g' *_barcode.txt
#extract reads
#epidermis
samtools view -H possorted_bam.bam > SAM_header
samtools view possorted_bam.bam | LC_ALL=C grep -F -f epidermis_barcode.txt > epidermis_SAM_body
cat SAM_header epidermis_SAM_body > epidermis.sam
samtools view -b epidermis.sam > epidermis.bam
#groud_meristem
samtools view -H possorted_bam.bam > SAM_header
samtools view possorted_bam.bam | LC_ALL=C grep -F -f ground.meristem_barcode.txt > ground.meristem_SAM_body
cat SAM_header ground.meristem_SAM_body > ground.meristem.sam
samtools view -b ground.meristem.sam > ground.meristem.bam
samtools sort  epidermis.bam -o epidermis.sort.bam
samtools sort  ground.meristem.bam -o ground.meristem.sort.bam
samtools index epidermis.sort.bam
samtools index ground.meristem.sort.bam
#bam2bw
bamCoverage --bam epidermis.sort.bam -o epidermis.bw
bamCoverage --bam ground.meristem.sort.bam -o ground.meristem.bw
#TOBIAS detect TF footprint 
#Tn5 bias correct
#conda install mamba -n base -c conda-forge
#mamba create -n tobias python=3.7
#mamba install -c bioconda tobias
#conda activate tobias
cd /public/home/xinwang/ZLH
#cd /public/home/chaohe/ZLH/second_time
for i in epidermis_DEGs ground.meristem_DEG;
do
TOBIAS ATACorrect --bam /public/home/xinwang/snATAC/fq/ajbk_AK58_root/outs/"$i".sort.bam \
--genome /public/home/xinwang/snRNA/db/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta \
--peaks "$i".bed  \
--outdir ATACorrect_test --cores 8
done
#scorebigwig calcute TF foot print
for i in epidermis_DEGs ground.meristem_DEG;
do
TOBIAS FootprintScores --signal /public/home/xinwang/ZLH/ATACorrect_test/"$i".sort_corrected.bw \
--regions "$i".bed \
--output "$i".footprint.bw --cores 8
done
#detect TF footprint
#annotation
for i in epidermis_DEGs ground.meristem_DEG;
do
awk '{print $1"\t"$2"\t"$3"\t.\t.\t"$5"\t"$5}' "$i".bed >"$i"_annotation.bed
done
#BINDetect
for i in epidermis_DEGs ground.meristem_DEG;
do
TOBIAS BINDetect --motifs \
/public/home/chaohe/ATAC/final/JASPAR2020_CORE_plants_non-redundant_pfms_jaspar.txt \
--signals "$i".footprint.bw \
--genome /public/home/xinwang/snRNA/db/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta \
--peaks "$i"_annotation.bed \
--peak_header merged_peaks_annotated_header.txt \
--outdir "$i"_BINDetect_output \
--cond_names atac --cores 8
done
#CreateNetwork
mkdir "$i"_annotation
mv "$i"_BINDetect_output/*/beds/*bound.bed "$i"_annotation
cat "$i"_annotation/*.bed >"$i"_merged_bound.bed
#--motif-pvalue 1e-6  --bound-pvalue 0.001
awk '{print $4}' "$i"_merged_bound.bed >4R1.txt
perl -p -i -e 's/_MA*.*//g' 4R1.txt
perl -p -i -e 's/TraesCSchr/TraesCS/g' "$i"_merged_bound.bed
paste 4R1.txt "$i"_merged_bound.bed | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t.\t"$11"\t"$12"\t"$13}' >"$i"_tf_merged_boundR1.bed
#perl -p -i -e 's/_MA*.*//g' BINDetect_outputR1/*/beds/*_bound.bed
#TOBIAS CreateNetwork --TFBS BINDetect_outputR1/*/beds/*_bound.bed --origin motif2gene_mapping.txt 
cd /public/home/xinwang/ZLH
i="epidermis"
TOBIAS CreateNetwork --TFBS "$i"_tf_merged_boundR1.bed --origin /public/home/chaohe/ZLH/first_time/motif2gene_mapping.txt
i="ground.meristem"
mkdir "$i"_annotation
mv "$i"_BINDetect_output/*/beds/*bound.bed "$i"_annotation
cat "$i"_annotation/*.bed >"$i"_merged_bound.bed
awk '{print $4}' "$i"_merged_bound.bed >4R1.txt
perl -p -i -e 's/_MA*.*//g' 4R1.txt
perl -p -i -e 's/TraesCSchr/TraesCS/g' "$i"_merged_bound.bed
paste 4R1.txt "$i"_merged_bound.bed | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t.\t"$11"\t"$12"\t"$13}' >"$i"_tf_merged_boundR1.bed
mkdir groundMeristemNetwork
cd /public/home/xinwang/ZLH/groundMeristemNetwork
TOBIAS CreateNetwork --TFBS /public/home/xinwang/ZLH/"$i"_tf_merged_boundR1.bed --origin /public/home/chaohe/ZLH/first_time/motif2gene_mapping.txt
cat *path_edges.txt > epidermis_path_edges.txt
cat *paths.txt > epidermis_paths.txt

###SCENIC interfer co-expressed gene regulated network
#creat cis-target databases
#load software
cd /public/home/chaohe/create_cisTarget_databases
fasta_filename="left2right.fa"
motifs_dir="motifs_cb_format"
motifs_list_filename="motifs.lst"
db_prefix="3500-1500"
nbr_threads=20
cd /public/home/chaohe/create_cisTarget_databases
###Score all motifs at once and create rankings
/public/home/chaohe/create_cisTarget_databases/create_cistarget_motif_databases.py \
-f "${fasta_filename}" \
-M "${motifs_dir}" \
-m "${motifs_list_filename}" \
-o "${db_prefix}" \
-t "${nbr_threads}"


























