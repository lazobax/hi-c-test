rule all:
    input:
        "data/draft_genome_quast_report.html",
        "data/draft_genome_busco_summary.txt",
        "data/scaffold_quast_report.html",
        "data/scaffold_busco_summary.txt"
rule download_draft_genome:
    output:
        filename = "draft-genome-data.zip"
    conda:
        "envs/ncbi-datasets.yaml"
    shell:
        "datasets download genome accession GCA_017916325.1 --filename {output.filename} "

rule move_genome_data:
    input:
        file = "draft-genome-data.zip"
    output:
        "data/draft-genome.fna"
    shell:
        '''
        unzip {input.file} -d draft-genome-data
        mv draft-genome-data/ncbi_dataset/data/GCA_017916325.1/*.fna {output}
        rm -rf draft-genome-data
        rm {input.file}
        '''

rule download_hiC_reads:
    output:
        "data/SRR18440378_1.fastq",
        "data/SRR18440378_2.fastq"
    conda:
        "envs/sra-tools.yaml"
    shell:
        '''
        prefetch SRR18440378
        fasterq-dump --split-files --outdir data/ SRR18440378
        '''

rule align_reads:
    input:
        draft_genome = "data/draft-genome.fna",
        read1 = "data/SRR18440378_1.fastq",
        read2 = "data/SRR18440378_2.fastq"
    output:
        "data/alignment.sam"
    conda:
        "envs/alignment.yaml"
    shell:
        "bwa mem {input.draft_genome} {input.read1} {input.read2} > {output}"

rule sam_to_bam:
    input:
        "data/alignment.sam"
    output:
        "data/alignment.bam"
    conda:
        "envs/alignment.yaml"
    shell:
        "samtools view -S -b {input} > {output}"

rule sort_bam:
    input:
        "data/alignment.bam"
    output:
        "data/sorted_alignment.bam"
    conda:
        "envs/alignment.yaml"
    shell:
        "samtools sort {input} -o {output}"

rule bam_to_bed:
    input:
        "data/sorted_alignment.bam"
    output:
        "data/alignment.bed"
    conda:
        "envs/alignment.yaml"
    shell:
        "bedtools bamtobed -i {input} > {output}"

rule sort_bed:
    input:
        "data/alignment.bed"
    output:
        "data/sorted_alignment.bed"
    conda:
        "envs/alignment.yaml"
    shell:
        "sort -k 4 {input} > {output}"

rule generate_contig_lengths:
    input:
        contigs = "data/draft-genome.fna"  # Same as the genome file
    output:
        "data/draft-genome.fna.fai"  # Index file for contig lengths
    conda:
        "envs/alignment.yaml"
    shell:
        "samtools faidx {input.contigs}"

rule scaffold:
    input:
        genome = "data/draft-genome.fna",
        lengths = "data/draft-genome.fna.fai",
        alignment = "data/sorted_alignment.bed",

    conda:
        "envs/SALSA"
    output:
        "data/scaffolds_FINAL.fasta"

    shell:
        "python SALSA/run_pipeline.py -a {input.genome} -l {input.lengths} -b {input.alignment} -e GATC -o data/scaffolds "

rule download_busco_dataset:
    output:
        "data/busco_lineage_data/lineage_eukaryota_odb10"
    conda:
        "envs/qc.yaml"
    shell:
        "busco --download lineage_eukaryota_odb10 --out data/busco_lineage_data"

rule qc_draft:
    input:
        genome = "data/draft-genome.fna"
    output:
        quast_report = "data/draft_genome_quast_report.html",
        busco_summary = "data/draft_genome_busco_summary.txt"
    conda:
        "envs/qc.yaml"
    shell:
        """
        # Run BUSCO
        busco -i {input.genome} -o data/draft_genome_busco_summary -l data/busco_lineage_data/lineage_eukaryota_odb10 -m genome -c 8

        # Run QUAST
        quast.py {input.genome} -o data/draft_genome_quast_report --threads 8
        """

rule qc_scaffold:
    input:
        genome = "data/scaffolds_FINAL.fasta"
    output:
        quast_report = "data/scaffold_quast_report.html",
        busco_summary = "data/scaffold_busco_summary.txt"
    conda:
        "envs/qc.yaml"

    shell:
        """
        # Run BUSCO
        busco -i {input.genome} -o data/scaffold_busco_summary -l data/busco_lineage_data/lineage_eukaryota_odb10 -m genome -c 8

        # Run QUAST
        quast.py {input.genome} -o data/scaffold_quast_report --threads 8
        """
