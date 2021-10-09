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

# 3.
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
ls ./map/*.bam | cut -d "/" -f 3 | cut -d "." -f 1,2 | awk '{printf "%s\t",$1}' > colname
### insert colname for first five columns 
sed 's/^/id\tstart\tend\tlocation\tdescription\t/' colname > expression.txt
cat count.txt >> expression.txt


#5. Statics Mean Calculation
