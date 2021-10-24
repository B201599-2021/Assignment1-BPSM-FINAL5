#!/bin/bash

#copying the files to the current working directory
cp /localdisk/data/BPSM/AY21/fastq/* .
cp /localdisk/data/BPSM/AY21/Tcongo_genome/* .
cp "/localdisk/home/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed" .

#saving the current working directory in a variable called cwd
cwd=$(pwd)

# running fastqc on all *.fq.gz files
echo "running fastqc"
fastqc *.fq.gz

#displaying the output for the fastqc results on the firefox web browser
echo "displaying fastqc results in a web browser"
for i in $(ls $cwd/*fastqc.html | rev | cut -c 6- | cut -c -22 | rev | uniq); do firefox $i.html; done;

#using bowtie2 to index the reference using the name reference_prefix
echo "using bowtie2 to index the reference using the name reference_prefix"
/localdisk/home/software/bowtie2-2.4.4/bowtie2-build TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz reference_prefix

#using a for loop to align all the fq.gz files with the refernce_prefix and storing the output in sam files
#using bowtie2 to align
echo "using bowtie2 to align the files and storing the output in sam files"
for i in $(ls $cwd/*.fq.gz | rev | cut -c 9- | rev | uniq); do bowtie2 -x reference_prefix -1 ${i}_1.fq.gz -2 ${i}_2.fq.gz -S ${i}.sam; done;

#converting the sam files to bam files using samtools
echo "converting the sam files to bam files using samtools"
for i in $( ls $cwd/*.sam | rev | cut -c 5- | rev | uniq); do samtools view -S -b ${i}.sam > ${i}.bam; done;

#sorting all the bam files using samtools with a for loop
echo "sorting the bam files using samtools"
for i in $(ls $cwd/*.bam | rev | cut -c 5- | rev | uniq); do samtools sort ${i}.bam -o ${i}.sorted.bam; done;

#indexing all the previously sorted bam files using samtools with a for loop
echo "indexing the bam files using samtools"
for i in $( ls $cwd/*.sorted.bam | rev | cut -c 12- | rev | uniq); do samtools index ${i}.sorted.bam; done;

#generating counts data of indexed cnt files using bedtools
echo "generating counts data of indexed cnt files using bedtools"
for i in $( ls $cwd/*.sorted.bam | rev | cut -c 5- | rev | uniq); do bedtools multicov -D -q 10 -bams ${i}.bam -bed TriTrypDB-46_TcongolenseIL3000_2019.bed > ${i}.bam.bedtools.cnt; done;

#clone1 group1 Clone1_time0h_uninduced- calculating the average for this group 
#using a for loop to find the file names for this group
for i in $(awk '{FS="\t"; if($5 == "Uninduced" && $2 == "Clone1" && $4 == 0){print $6;}}' 100k.fqfiles | rev | cut -c 9- | rev | uniq);
do echo $i; done
awk '{FS="\t"; OFS="\t"; {print $5,$NF;}}' 100k.C1-1-501.sorted.bam.bedtools.cnt >> c1g11.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C1-2-502.sorted.bam.bedtools.cnt >> c1g12.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C1-3-503.sorted.bam.bedtools.cnt >> c1g13.txt
paste c1g12.txt c1g13.txt > c1g123.txt
#pasting the gene descriptions along with the counts data in a single plain text tab-delimited file
paste c1g11.txt c1g123.txt > c1g1f.txt
#calculating the average per gene 
awk '{FS="\t"; OFS="\t"; {print $1,$2,$3,$4,($2 + $3 + $4)/3;}}' c1g1f.txt > averagec1g1.txt
echo "The average counts per gene for Clone 1 at Time 0h Uninduced is:"
cat averagec1g1.txt 

#clone1 group 2 Clone1_time24h_uninduced
#using a for loop to find the file names for this group
for i in $(awk '{FS="\t"; if($5 == "Uninduced" && $2 == "Clone1" && $4 == 24){print $6;}}' 100k.fqfiles | rev | cut -c 9- | rev | uniq);
do echo $i; done
awk '{FS="\t"; OFS="\t"; {print $5,$NF;}}' 100k.C1-4-516.sorted.bam.bedtools.cnt >> c1g21.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C1-5-517.sorted.bam.bedtools.cnt >> c1g22.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C1-6-518.sorted.bam.bedtools.cnt >> c1g33.txt
paste c1g22.txt c1g33.txt > c1g223.txt
#pasting the gene descriptions along with the counts data in a single plain text tab-delimited file
paste c1g21.txt c1g223.txt > c1g2f.txt
#calculating the average per gene 
awk '{FS="\t"; OFS="\t"; {print $1,$2,$3,$4,($2 + $3 + $4)/3;}}' c1g2f.txt > averagec1g2.txt
echo "The average counts per gene for Clone 1 at Time 24h Uninduced is:"
cat averagec1g2.txt 

#clone1 group 3 Clone1_time48h_uninduced
#using a for loop to find the file names for this group
for i in $(awk '{FS="\t"; if($5 == "Uninduced" && $2 == "Clone1" && $4 == 48){print $6;}}' 100k.fqfiles | rev | cut -c 9- | rev | uniq);
do echo $i; done
awk '{FS="\t"; OFS="\t"; {print $5,$NF;}}' 100k.C1-4-536.sorted.bam.bedtools.cnt >> c1g31.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C1-5-539.sorted.bam.bedtools.cnt >> c1g32.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C1-6-541.sorted.bam.bedtools.cnt >> c1g333.txt
paste c1g32.txt c1g333.txt > c1g323.txt
#pasting the gene descriptions along with the counts data in a single plain text tab-delimited file
paste c1g31.txt c1g323.txt > c1g3f.txt
#calculating the average per gene 
awk '{FS="\t"; OFS="\t"; {print $1,$2,$3,$4,($2 + $3 + $4)/3;}}' c1g3f.txt > averagec1g3.txt
echo "The average counts per gene for Clone 1 at Time 48h Uninduced is:"
cat averagec1g3.txt 

#clone1 group 4 Clone1_time24h_induced
#using a for loop to find the file names for this group
for i in $(awk '{FS="\t"; if($5 == "Induced" && $2 == "Clone1" && $4 == 24){print $6;}}' 100k.fqfiles | rev | cut -c 9- | rev | uniq);
do echo $i; done
awk '{FS="\t"; OFS="\t"; {print $5,$NF;}}' 100k.C1-1-513.sorted.bam.bedtools.cnt >> c1g41.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C1-2-514.sorted.bam.bedtools.cnt >> c1g42.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C1-3-515.sorted.bam.bedtools.cnt >> c1g43.txt
paste c1g42.txt c1g43.txt > c1g423.txt
#pasting the gene descriptions along with the counts data in a single plain text tab-delimited file
paste c1g41.txt c1g423.txt > c1g4f.txt
#calculating the average per gene 
awk '{FS="\t"; OFS="\t"; {print $1,$2,$3,$4,($2 + $3 + $4)/3;}}' c1g4f.txt > averagec1g4.txt
echo "The average counts per gene for Clone 1 at Time 24h induced is:"
cat averagec1g4.txt 

#clone1 group 5 Clone1_time48h_induced
#using a for loop to find the file names for this group
for i in $(awk '{FS="\t"; if($5 == "Induced" && $2 == "Clone1" && $4 == 48){print $6;}}' 100k.fqfiles | rev | cut -c 9- | rev | uniq);
do echo $i; done
awk '{FS="\t"; OFS="\t"; {print $5,$NF;}}' 100k.C1-1-533.sorted.bam.bedtools.cnt >> c1g51.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C1-2-534.sorted.bam.bedtools.cnt >> c1g52.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C1-3-535.sorted.bam.bedtools.cnt >> c1g53.txt
paste c1g52.txt c1g53.txt > c1g523.txt
#pasting the gene descriptions along with the counts data in a single plain text tab-delimited file
paste c1g51.txt c1g523.txt > c1g5f.txt
#Calculating the average per gene 
awk '{FS="\t"; OFS="\t"; {print $1,$2,$3,$4,($2 + $3 + $4)/3;}}' c1g5f.txt > averagec1g5.txt
echo "The average counts per gene for Clone 1 at Time 48h induced is:"
cat averagec1g5.txt


#clone2 group 1 Clone2_time0h_uninduced
#using a for loop to find the file names for this group
for i in $(awk '{FS="\t"; if($5 == "Uninduced" && $2 == "Clone2" && $4 == 0){print $6;}}' 100k.fqfiles | rev | cut -c 9- | rev | uniq);
do echo $i; done
awk '{FS="\t"; OFS="\t"; {print $5,$NF;}}' 100k.C2-1-504.sorted.bam.bedtools.cnt >> c2g11.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C2-2-505.sorted.bam.bedtools.cnt >> c2g12.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C2-3-506.sorted.bam.bedtools.cnt >> c2g13.txt
paste c2g12.txt c2g13.txt > c2g123.txt
#pasting the gene descriptions along with the counts data in a single plain text tab-delimited file
paste c2g11.txt c2g123.txt > c2g1f.txt
#calculating the average per gene 
awk '{FS="\t"; OFS="\t"; {print $1,$2,$3,$4,($2 + $3 + $4)/3;}}' c2g1f.txt > averagec2g1.txt
echo "The average counts per gene for Clone 2 at Time 0h Uninduced is:"
cat averagec2g1.txt

#clone2 group 2 Clone2_time24h_uninduced
#using a for loop to find the file names for this group
for i in $(awk '{FS="\t"; if($5 == "Uninduced" && $2 == "Clone2" && $4 == 24){print $6;}}' 100k.fqfiles | rev | cut -c 9- | rev | uniq);
do echo $i; done
awk '{FS="\t"; OFS="\t"; {print $5,$NF;}}' 100k.C2-4-522.sorted.bam.bedtools.cnt >> c2g21.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C2-5-523.sorted.bam.bedtools.cnt >> c2g22.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C2-6-524.sorted.bam.bedtools.cnt >> c2g23.txt
paste c2g22.txt c2g23.txt > c2g223.txt
#pasting the gene descriptions along with the counts data in a single plain text tab-delimited file
paste c2g21.txt c2g223.txt > c2g2f.txt
#calculating the average per gene 
awk '{FS="\t"; OFS="\t"; {print $1,$2,$3,$4,($2 + $3 + $4)/3;}}' c2g2f.txt > averagec2g2.txt
cat averagec2g2.txt
echo "The average counts per gene for Clone 2 at Time 24h Uninduced is:"
cat averagec2g2.txt

#clone2 group 3 Clone2_time48h_uninduced
#using a for loop to find the file names for this group
for i in $(awk '{FS="\t"; if($5 == "Uninduced" && $2 == "Clone2" && $4 == 48){print $6;}}' 100k.fqfiles | rev | cut -c 9- | rev | uniq);
do echo $i; done
awk '{FS="\t"; OFS="\t"; {print $5,$NF;}}' 100k.C2-4-545.sorted.bam.bedtools.cnt >> c2g31.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C2-5-546.sorted.bam.bedtools.cnt >> c2g32.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C2-6-548.sorted.bam.bedtools.cnt >> c2g33.txt
paste c2g32.txt c2g33.txt > c2g323.txt
#pasting the gene descriptions along with the counts data in a single plain text tab-delimited file
paste c2g31.txt c2g323.txt > c2g3f.txt
#calculating the average per gene 
awk '{FS="\t"; OFS="\t"; {print $1,$2,$3,$4,($2 + $3 + $4)/3;}}' c2g3f.txt > averagec2g3.txt
echo "The average counts per gene for Clone 2 at Time 48h Uninduced is:"
cat averagec2g3.txt 

#clone2 group 4 Clone2_time24h_induced
#using a for loop to find the file names for this group
for i in $(awk '{FS="\t"; if($5 == "Induced" && $2 == "Clone2" && $4 == 24){print $6;}}' 100k.fqfiles | rev | cut -c 9- | rev | uniq);
do echo $i; done
awk '{FS="\t"; OFS="\t"; {print $5,$NF;}}' 100k.C2-1-519.sorted.bam.bedtools.cnt >> c2g41.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C2-2-520.sorted.bam.bedtools.cnt >> c2g42.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C2-3-521.sorted.bam.bedtools.cnt >> c2g43.txt
paste c2g42.txt c2g43.txt > c2g423.txt
#pasting the gene descriptions along with the counts data in a single plain text tab-delimited file
paste c2g41.txt c2g423.txt > c2g4f.txt
#calculating the average per gene 
awk '{FS="\t"; OFS="\t"; {print $1,$2,$3,$4,($2 + $3 + $4)/3;}}' c2g4f.txt > averagec2g4.txt
echo "The average counts per gene for Clone 2 at Time 24h induced is:"
cat averagec2g4.txt 

#clone2 group 5 Clone2_time48h_induced
#using a for loop to find the file names for this group
for i in $(awk '{FS="\t"; if($5 == "Induced" && $2 == "Clone2" && $4 == 48){print $6;}}' 100k.fqfiles | rev | cut -c 9- | rev | uniq);
do echo $i; done
awk '{FS="\t"; OFS="\t"; {print $5,$NF;}}' 100k.C2-1-542.sorted.bam.bedtools.cnt >> c2g51.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C2-2-543.sorted.bam.bedtools.cnt >> c2g52.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.C2-3-544.sorted.bam.bedtools.cnt >> c2g53.txt
paste c2g52.txt c2g53.txt > c2g523.txt
#pasting the gene descriptions along with the counts data in a single plain text tab-delimited file
paste c2g51.txt c2g523.txt > c2g5f.txt
#calculating the average per gene 
awk '{FS="\t"; OFS="\t"; {print $1,$2,$3,$4,($2 + $3 + $4)/3;}}' c2g5f.txt > averagec2g5.txt
echo "The average counts per gene for Clone 2 at Time 48h induced is:"
cat averagec2g5.txt

#wildtype group 1 wildtype_time0h_uninduced
#using a for loop to find the file names for this group
for i in $(awk '{FS="\t"; if($5 == "Uninduced" && $2 == "WT" && $4 == 0){print $6;}}' 100k.fqfiles | rev | cut -c 9- | rev | uniq);
do echo $i; done
awk '{FS="\t"; OFS="\t"; {print $5,$NF;}}' 100k.WT-1-507.sorted.bam.bedtools.cnt >> c3g11.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.WT-2-511.sorted.bam.bedtools.cnt >> c3g12.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.WT-3-512.sorted.bam.bedtools.cnt >> c3g13.txt
paste c3g12.txt c3g13.txt > c3g123.txt
#pasting the gene descriptions along with the counts data in a single plain text tab-delimited file
paste c3g11.txt c3g123.txt > c3g1f.txt
#calculating the average per gene 
awk '{FS="\t"; OFS="\t"; {print $1,$2,$3,$4,($2 + $3 + $4)/3;}}' c3g1f.txt > averagec3g1.txt
echo "The average counts per gene for WildType at Time 0h Uninduced is:"
cat averagec3g1.txt

#wildtype group 2 wildtype_time24h_uninduced
#using a for loop to find the file names for this group
for i in $(awk '{FS="\t"; if($5 == "Uninduced" && $2 == "WT" && $4 == 24){print $6;}}' 100k.fqfiles | rev | cut -c 9- | rev | uniq);
do echo $i; done
awk '{FS="\t"; OFS="\t"; {print $5,$NF;}}' 100k.WT-4-530.sorted.bam.bedtools.cnt >> c3g21.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.WT-5-531.sorted.bam.bedtools.cnt >> c3g22.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.WT-6-532.sorted.bam.bedtools.cnt >> c3g23.txt
paste c3g22.txt c3g23.txt > c3g223.txt
#pasting the gene descriptions along with the counts data in a single plain text tab-delimited file
paste c3g21.txt c3g223.txt > c3g2f.txt
#calculating the average per gene 
awk '{FS="\t"; OFS="\t"; {print $1,$2,$3,$4,($2 + $3 + $4)/3;}}' c3g2f.txt > averagec3g2.txt
echo "The average counts per gene for WildType at Time 24h Uninduced is:"
cat averagec3g2.txt

#wildtype group 3 wildtype_time48h_uninduced
#using a for loop to find the file names for this group
for i in $(awk '{FS="\t"; if($5 == "Uninduced" && $2 == "WT" && $4 == 48){print $6;}}' 100k.fqfiles | rev | cut -c 9- | rev | uniq);
do echo $i; done
awk '{FS="\t"; OFS="\t"; {print $5,$NF;}}' 100k.WT-4-553.sorted.bam.bedtools.cnt >> c3g31.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.WT-5-554.sorted.bam.bedtools.cnt >> c3g32.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.WT-6-555.sorted.bam.bedtools.cnt >> c3g33.txt
paste c3g32.txt c3g33.txt > c3g323.txt
#pasting the gene descriptions along with the counts data in a single plain text tab-delimited file
paste c3g31.txt c3g323.txt > c3g3f.txt
#calculating the average per gene 
awk '{FS="\t"; OFS="\t"; {print $1,$2,$3,$4,($2 + $3 + $4)/3;}}' c3g3f.txt > averagec3g3.txt
echo "The average counts per gene for WildType at Time 48h Uninduced is:"
cat averagec3g3.txt 

#wildtype group 4 wildtype_time24h_induced
#using a for loop to find the file names for this group
for i in $(awk '{FS="\t"; if($5 == "Induced" && $2 == "WT" && $4 == 24){print $6;}}' 100k.fqfiles | rev | cut -c 9- | rev | uniq);
do echo $i; done
awk '{FS="\t"; OFS="\t"; {print $5,$NF;}}' 100k.WT-1-525.sorted.bam.bedtools.cnt >> c3g41.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.WT-2-526.sorted.bam.bedtools.cnt >> c3g42.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.WT-3-529.sorted.bam.bedtools.cnt >> c3g43.txt
paste c3g42.txt c3g43.txt > c3g423.txt
#pasting the gene descriptions along with the counts data in a single plain text tab-delimited file
paste c3g41.txt c3g423.txt > c3g4f.txt
#calculating the average per gene 
awk '{FS="\t"; OFS="\t"; {print $1,$2,$3,$4,($2 + $3 + $4)/3;}}' c3g4f.txt > averagec3g4.txt
echo "The average counts per gene for WildType at Time 24h induced is:"
cat averagec3g4.txt 

#wildtype group 5 wildtype_time48h_induced
#using a for loop to find the file names for this group
for i in $(awk '{FS="\t"; if($5 == "Induced" && $2 == "WT" && $4 == 48){print $6;}}' 100k.fqfiles | rev | cut -c 9- | rev | uniq);
do echo $i; done
awk '{FS="\t"; OFS="\t"; {print $5,$NF;}}' 100k.WT-1-550.sorted.bam.bedtools.cnt >> c3g51.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.WT-2-551.sorted.bam.bedtools.cnt >> c3g52.txt
awk '{FS="\t"; OFS="\t"; {print $NF;}}' 100k.WT-3-552.sorted.bam.bedtools.cnt >> c3g53.txt
paste c3g52.txt c3g53.txt > c3g523.txt
#pasting the gene descriptions along with the counts data in a single plain text tab-delimited file
paste c3g51.txt c3g523.txt > c3g5f.txt
#calculating the average per gene 
awk '{FS="\t"; OFS="\t"; {print $1,$2,$3,$4,($2 + $3 + $4)/3;}}' c3g5f.txt > averagec3g5.txt
echo "The average counts per gene for WildType at Time 48h induced is:"
cat averagec2g5.txt
