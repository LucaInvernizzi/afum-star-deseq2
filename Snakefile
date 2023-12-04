SAMPLES = ["control_1", "control_2", "control_3", 
    "feStarved_1", "feStarved_2", "feStarved_3", 
    "h2o2Treated_1", "h2o2Treated_2", "h2o2Treated_3",
    "feStarved_h2o2Treated_1", "feStarved_h2o2Treated_2", "feStarved_h2o2Treated_3"]
# CHECK Params in STAR align step before running pipeline

rule all:
    input:
        expand("star_align/{samples}/Aligned.sortedByCoord.out.bam", samples = SAMPLES), 
        expand("data/samples/fastqc/{samples}_fastqc.html", samples = SAMPLES),
        expand("feature_counts/{samples}/{samples}_counts.txt", samples = SAMPLES),
        "results/counts_tab.txt"

rule bbduck_trim:
    input:
        "data/samples/fastq/{samples}.fastq.gz"
    output:
        "data/samples/trim/{samples}.fastq"
    shell: # bbduk.sh is like --help, but used advised command
        """
        bbduk.sh -Xmx1g in={input} out={output} \
        ref=/home/luca/mambaforge/envs/afum-star-deseq2/opt/bbmap-39.01-1/resources/adapters.fa \
        ktrim=r k=23 mink=11 hdist=1 tpe tbo \
        """

rule fastqc:
    input: # solvable by simply give list of input but try different strategy first
        "data/samples/trim/{samples}.fastq"
    output:
        html = "data/samples/fastqc/{samples}_fastqc.html",
        zip = "data/samples/fastqc/{samples}_fastqc.zip"
    shell:
        "fastqc {input} -o data/samples/fastqc"
        
rule star_index:
    input:
        fa = "data/afum_genome.fa",
        gff = "data/afum_ann.gff" # could be converted to gtf with agat or cufflinks
    output:
        directory("data/star_index")
    threads: 6
    shell:
        """
        mkdir -p {output} && \
        STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {output} \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gff} \
        --sjdbGTFtagExonParentTranscript Parent \
        --sjdbOverhang 100 \
        """ 

rule star_se:
    input:
        fq1="data/samples/trim/{sample}.fastq",
        idx="data/star_index",
    output:
        # see STAR manual for additional output files
        aln="star_align/{sample}/Aligned.sortedByCoord.out.bam",
        log="star_align/logs/{sample}/Log.out",
        log_final="star_align/logs/{sample}/Log.final.out",
        unmapped="star_align/{sample}/unmapped.fastq",
    log:
        "star_align/logs/se/{sample}.log",
    params:
        # optional parameters, maybe it can be split on more lines?
        extra="--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMattributes XS",
    threads: 6
    wrapper:
        "v3.0.0/bio/star/align"

rule featureCounts: # from the subread package
    input:
        aln="star_align/{sample}/Aligned.sortedByCoord.out.bam",
        gff = "data/afum_ann.gff"
    output:
        features = "feature_counts/{sample}/{sample}_counts.txt"
    shell:
        """
        featureCounts -t exon -g gene_id \
        -a {input.gff} \
        -o {output.features} \
        {input.aln}
        """

rule counts_tab:
    input:
        expand("feature_counts/{samples}/{samples}_counts.txt", samples = SAMPLES)
    output:
        "results/counts_tab.txt"
    script:
        "scripts/counts_tab.R"

# counts_tab.R is now ready to be used in a DESeq2 analysis