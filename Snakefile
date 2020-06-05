SAMPLES, = glob_wildcards('input/{sample}_R1_001.fastq.gz')
PATTERNS=["R1_001", "R2_001"]
rule allout:
    input:
         expand("output/fastqc_out/{sample}_{pattern}_fastqc.html",sample =SAMPLES, pattern=PATTERNS),
         "output/multiqc_out/multiqc_report.html",
         expand("output/fastqc_trim_out/{sample}_{pattern}_fastqc.html", sample =SAMPLES, pattern=PATTERNS),
         "output/multiqc_trim_out/multiqc_report.html",
         expand('output/bbduk_out/{sample}_{pattern}_trimmed.fastq.gz', sample =SAMPLES, pattern=PATTERNS),
         directory('STAR_Output'),
         expand('output/STAR/{sample}_pass/Aligned.sortedByCoord.out.bam', sample =SAMPLES),
         expand('output/FeatureCounts/{sample}_feature_counts_s1.txt', sample =SAMPLES),
         expand('output/FeatureCounts/{sample}_feature_counts_s1.txt', sample =SAMPLES),
         'output/Collibri_counts.txt',
         'output/KAPA_counts.txt'
rule fastqc:
    input:
        "input/{sample}_{pattern}.fastq.gz"
    output:
        "output/fastqc_out/{sample}_{pattern}_fastqc.html"
    conda:
        "environment.yaml"
    shell:
        "fastqc {input} -o ./output/fastqc_out/"

rule multiqc:
    input:
        "output/fastqc_out/"
    output:
        "output/multiqc_out/multiqc_report.html"
    conda:
        "environment.yaml"
    shell:
        "multiqc output/fastqc_out -o output/multiqc_out"

rule bbduk:
    input:
        "input/{sample}_{pattern}.fastq.gz"
    output:
        "output/bbduk_out/{sample}_{pattern}_trimmed.fastq.gz"
    conda:
        "environment.yaml"
    shell:
        "bbduk.sh in={input} out={output} ref=input/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10"
rule fastqc_trimmed:
    input:
        "output/bbduk_out/{sample}_{pattern}.fastq.gz"
    output:
        "output/fastqc_trim_out/{sample}_{pattern}_fastqc.html"
    conda:
        "environment.yaml"
    shell:
        "fastqc {input} -o ./output/fastqc_trim_out/"

rule multiqc_trim:
    input:
        "output/fastqc_trim_out/"
    output:
        "output/multiqc_trim_out/multiqc_report.html"
    conda:
        "environment.yaml"
    shell:
        "multiqc output/fastqc_trim_out -o output/multiqc_trim_out"
rule index:
    input:
        fa = 'ref_files/chr19_20Mb.fa',
        gtf = 'ref_files/chr19_20Mb.gtf'
    output:
        'STAR_Output'
    threads: 4
    shell:
        'mkdir {output} && '
        'STAR --runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {output} '
        '--genomeFastaFiles {input.fa} '
        '--sjdbGTFfile {input.gtf} '
        '--sjdbOverhang 100'
rule pass:
        input:
             R1L1 = 'output/bbduk_out/{sample}_R1_001_trimmed.fastq.gz',
             R2L1 = 'output/bbduk_out/{sample}_R2_001_trimmed.fastq.gz',
             refdir = 'STAR_Output'
        params:
            outfirst = 'output/STAR',
            outdir = '{sample}_pass',
            id = '{sample}'
        output:
            'output/STAR/{sample}_pass/Aligned.sortedByCoord.out.bam'
        threads: 4
        shell:
            'cd {params.outfirst} && '
            'rm -rf {params.outdir} &&'
            'mkdir {params.outdir} && '
            'cd {params.outdir} && '
            'STAR --runThreadN {threads} '
            '--genomeDir ~/Transcriptomics_exercise/{input.refdir} '  ###~/Transcriptomics_exercise/STAR_Output'
            '--readFilesIn ~/Transcriptomics_exercise/{input.R1L1} ~/Transcriptomics_exercise/{input.R2L1} '
            '--readFilesCommand zcat '
            '--outSAMtype BAM SortedByCoordinate '
rule indexing_bam:
        input:
            'output/STAR/{sample}_pass/Aligned.sortedByCoord.out.bam',
        output:
            'output/STAR/{sample}_pass/Aligned.sortedByCoord.out.bam.bai'
        shell:
            'samtools index {input}'
rule featureCounts1:
        input:
            mapp_res='output/STAR/{sample}_pass/Aligned.sortedByCoord.out.bam',
            ref_gtf='ref_files/chr19_20Mb.gtf'
        output:
            'output/FeatureCounts/{sample}_feature_counts_s1.txt',
        shell:
            'featureCounts -p -t exon -g gene_id -a {input.ref_gtf} -o {output} {input.mapp_res} -s 1'
rule featureCounts2:
        input:
            mapp_res='output/STAR/{sample}_pass/Aligned.sortedByCoord.out.bam',
            ref_gtf='ref_files/chr19_20Mb.gtf'
        output:
            'output/FeatureCounts/{sample}_feature_counts_s2.txt',
        shell:
            'featureCounts -p -t exon -g gene_id -a {input.ref_gtf} -o {output} {input.mapp_res} -s 2'
rule feature_count_combined:
        input:
        #    mapp_res='output/STAR/{sample}_pass/Aligned.sortedByCoord.out.bam',
            ref_gtf='ref_files/chr19_20Mb.gtf'
        output:
            'output/Collibri_counts.txt',
            'output/KAPA_counts.txt'
        shell:
            'featureCounts -p -t exon -g gene_id -a ref_files/chr19_20Mb.gtf -o output/Collibri_counts.txt -M output/STAR/Collibri*/*.bam -s 1 &&'
            'featureCounts -p -t exon -g gene_id -a ref_files/chr19_20Mb.gtf -o output/KAPA_counts.txt -M output/STAR/KAPA*/*.bam -s 2'
