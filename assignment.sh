# 1.Assign all files needed for analysis
# 1.1.Assign raw fastq file path
fq_path='/localdisk/data/BPSM/AY21/fastq'
#fq_path='/localdisk/home/s2172876/Assignment1/rawfastq'
# 1.2.Genome file And Bed file
fa_file='/localdisk/data/BPSM/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz'
bed_file='/localdisk/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed'
ls -l $fq_path | awk '{FS=" ";  {print $9;}}' | sort | grep "fq.gz" | awk '{print(substr($1,1,13))}' | sort | uniq >./sampleid.txt

# 2.Perpreration for mapping
mkdir -p genome
cp $fa_file ./genome/Genome.fasta.gz
gzip -d ./genome/Genome.fasta.gz
hisat2-build ./genome/Genome.fasta ./genome/Genome
index_path='./genome/Genome'

# 3. fastq to bam(.bai)
## Using a loop to analysis each pair of fastq data
## Using hisat2 to align reads, Using samtools to generate bam files and its index
mkdir -p map
mkdir -p fastqc
cat ${PWD}/sampleid.txt | while read id
do
    # 3.1 Perform quality control
    ## -t thread: 5
    ## output to fastqc folder
    fastqc ${fq_path}/${id}_1.fq.gz ${fq_path}/${id}_2.fq.gz -t 5 -o ${PWD}/fastqc/
    # 3.2 Perform reads alignment
    ## -p thread: 5
    ## sam file output to map folder
    ## log message redirects to map folder and save with filename mapeval(mapping-evaluation)
    hisat2 -p 5 -x $index_path -1 ${fq_path}/${id}_1.fq.gz -2 ${fq_path}/${id}_2.fq.gz -S ./map/${id}.sam 2> "./map/${id}.mapeval" 
    # 3.3 sam to bam
    samtools sort -@ 5 ./map/${id}.sam -o ./map/${id}.bam
    # 3.4 bam index
    samtools index ./map/${id}.bam ./map/${id}.bam.bai
done

#4. Generate Count data
mkdir -p count
bedtools multicov -bams ./map/*.bam -bed /localdisk/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed > count.txt
## match sample id to count.txt
### extract bam name in order
ls ./map/*.bam | cut -d "/" -f 3 | cut -d "." -f 2 | awk '{printf "%s\t",$1}' > colname
### insert colname for first five columns 
sed 's/^/id\tstart\tend\tlocation\tdescription\t/' colname > expression.txt
cat count.txt >> expression.txt


#5. Statics Mean Calculation
cp /localdisk/data/BPSM/AY21/fastq/100k.fqfiles .
## extract 2,4,5 colnames to identity group(samples in a group identical in 2_4_5 row)
cut 100k.fqfiles -f2,4,5 | sed 's/\t/_/g' | paste 100k.fqfiles - | cut - -f1,8 | sort >> group
## extract expression matrix
cut -f 6- expression.txt > count.txt
## expression matrix transpose in order to merge with grouping information
awk '{for(i=1;i<=NF;i++){a[FNR,i]=$i}}END{for(i=1;i<=NF;i++){for(j=1;j<=FNR;j++){printf a[j,i]"\t"}print ""}}' count.txt > Texpression.txt
## merging expression matrix and grouping files
IFS=$'\t'
join -j 1 group Texpression.txt > expression_group.txt
## remove processing files
rm group count.txt Texpression.txt
## sort files with grouping column to make sure samples in sample group list adjacent
sed 's/ /\t/g' expression_group.txt | sort -k2,2 > Texpression.txt
## transposed expression matrix transpose
awk '{for(i=1;i<=NF;i++){a[FNR,i]=$i}}END{for(i=1;i<=NF;i++){for(j=1;j<=FNR;j++){printf a[j,i]"\t"}print ""}}' Texpression.txt > expression_group.txt
## paste sample with gene descrpition
cut -f 1,2,3,4,5 expression.txt > output.txt
sed -i "1a\group\tgroup\tgroup\tgroup\tgroup" output.txt
paste output.txt expression_group.txt > test.txt

