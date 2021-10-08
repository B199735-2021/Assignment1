# 1.Assign all files needed for analysis
# 1.1.Assign raw fastq file path
fq_path='/localdisk/data/BPSM/AY21/fastq'
# 1.2.Genome file And Bed file
fa_file='/localdisk/data/BPSM/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz'
bed_file='/localdisk/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed'
#dir=$(ls -l $fq_path | awk '{FS=" ";  {print $9;}}' | sort | grep "fq.gz")
# file_name=$(ls -l /localdisk/data/BPSM/AY21/fastq | awk '{FS=" ";  {print $9;}}' | sort | grep "fq.gz")
# count=$(ls -l /localdisk/data/BPSM/AY21/fastq | awk '{FS=" ";  {print $9;}}' | sort | grep "fq.gz" | wc -l)
# echo $(ls -l $fq_path | awk '{FS=" ";  {print $9;}}' | sort | grep "fq.gz" | awk '{print(substr($1,1,13))}' | sort | uniq) 
# ls -l $fq_path | awk '{FS=" ";  {print $9;}}' | sort | grep "fq.gz" | awk '{print(substr($1,1,13))}' | sort | uniq

# 2.Perpreration for mapping
mkdir genome
cp $fa_file ./genome
ls ./genome/ | gzip -d 
hisat2-build ./genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta ./genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome


# 3.
for file in $fq_path/*.fq.gz
do
  	echo ${file:0:45}_1.fq.gz
    # Perform quality control
    ## -t thread: 5
    ## output to fastqc folder  
    fastqc ${file:0:45}_1.fq.gz ${file:0:45}_2.fq.gz -t 5 -o ./fastqc/
done

# count=$(ls -l $fq_path | awk '{FS=" ";  {print $9;}}' | grep "fq.gz" | wc -l)
# #echo $dir
# for ((i=1;i<=count;i=i+2))
# do
#     echo $i
#     echo $(ls -l /localdisk/data/BPSM/AY21/fastq | awk '{FS=" ";  {print $9;}}' | sort | grep "fq.gz" | sed -n '{$i}p')
# done
# mkdir -p fastq
# fastqc /localdisk/data/BPSM/AY21/fastq/100k.C1-1-501_1.fq.gz /localdisk/data/BPSM/AY21/fastq/100k.C1-1-501_2.fq.gz -t 5 -o ./fastqc/
# fastqc --extract -o fastqc_outputs
# # wget -c http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
# # unzip Trimmomatic-0.39.zip
# # cd Trimmomatic-0.39

# # mkdir trim
# # java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -trimlog trim.log \
# # /localdisk/data/BPSM/AY21/fastq/100k.C1-1-501_1.fq.gz /localdisk/data/BPSM/AY21/fastq/100k.C1-1-501_2.fq.gz \
# # ./trim/100k.C1-1-501_1_paired.fq.gz ./trim/100k.C1-1-501_1_unpaired.fq.gz ./trim/100k.C1-1-501_2_paired.fq.gz ./trim/100k.C1-1-501_2_unpaired.fq.gz \
# # ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:true \
# # SLIDINGWINDOW:5:20 \
# # LEADING:3 \
# # TRAILING:3 \
# # MINLEN:36

# mkdir -p genome
# cp /localdisk/data/BPSM/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz ./genome
# gzip -d ./genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz
# hisat2-build ./genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta ./genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome

# #map
# mkdir -p map
# index=./genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome
# hisat2 -p 3 -x ./genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome -1 /localdisk/data/BPSM/AY21/fastq/100k.C1-1-501_1.fq.gz -2 /localdisk/data/BPSM/AY21/fastq/100k.C1-1-501_2.fq.gz \
# -S ./map/100k.C1-1-501.sam

# #sam to bam
# samtools sort -@ 1 -o ./map/100k.C1-1-501.bam ./map/100k.C1-1-501.sam
# # bam indexing
# samtools index ./map/100k.C1-1-501.bam ./map/100k.C1-1-501.bam.bai

# mkdir count
# bedtools multicov -D -q 10 -bams ./map/100k.C1-1-501.bam -bed /localdisk/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed > ./count/100k.C1-1-501.cnt
