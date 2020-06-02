SAMPLES, = glob_wildcards('input/{sample}_R1_001.fastq.gz')
PATTERNS=["R1_001", "R2_001"]
rule allout:
    input:
    #     expand('input/{sample}_{pattern}.fastq.gz', sample=SAMPLES, pattern=PATTERNS),
         expand('output/bbduk_out/{sample}_{pattern}_trimmed.fastq.gz', sample =SAMPLES, pattern=PATTERNS),
#         expand('output/STAR/{sample}_pass/Aligned.sortedByCoord.out.bam', sample=SAMPLES)
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
    #    "input/{sample}_L001_R1_001.fastq.gz"
    #    "input/{sample}_L001_R2_001.fastq.gz"
        expand("input/{sample}_{pattern}.fastq.gz", sample= SAMPLES, pattern = PATTERNS)
    output:
        expand("output/bbduk_out/{sample}_{pattern}_trimmed.fastq.gz", sample=SAMPLES, pattern = PATTERNS)
        #"output/{sample}_L001_R1_001.fastq.gz"
        #"ouput/{sample}_L001_R2_001.fastq.gz"
    conda:
        "environment.yaml"
    shell:
        #"bbduk.sh ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10 " +
        #"in1={input[0]} in2={input[1]} out1={output[0]} out2={output[1]}"
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
        fa = 'ref_files/chr19_20Mb.fa', # reference FASTA file
        gtf = 'ref_files/chr19_20Mb.gtf' # GTF file
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
        #    R1L1 = expand('output/bbduk_out/{sample}_R1_001_trimmed.fastq.gz', sample=SAMPLES),
        #    R2L1 = expand('output/bbduk_out/{sample}_R2_001_trimmed.fastq.gz', sample=SAMPLES),
        #    R1L1 = expand('output/bbduk_out/{sample}_R1_001_trimmed.fastq.gz', sample=SAMPLES),
        #    R2L1 = expand('output/bbduk_out/{sample}_R2_001_trimmed.fastq.gz', sample=SAMPLES),
             R1L1 = 'output/bbduk_out/{sample}_R1_001_trimmed.fastq.gz',
             R2L1 = 'output/bbduk_out/{sample}_R2_001_trimmed.fastq.gz',
             refdir = 'STAR_Output'
        params:
            outfirst = 'output/STAR',
            outdir = '{sample}_pass',
            id = '{sample}'
        output:
            #'{sample}_pass/Aligned.sortedByCoord.out.bam'
            'output/STAR/{sample}_pass/Aligned.sortedByCoord.out.bam'
        threads: 4
        shell:
        #    'rm -rf {params.outdir} &&' # be careful with this. I don't know why, but Snakemake had problems without this cleaning.
            'cd {params.outfirst} && '
            'rm -rf {params.outdir} &&'
            'mkdir {params.outdir} && '
            'cd {params.outdir} && '
            'STAR --runThreadN {threads} '
            '--genomeDir ~/Transcriptomics_exercise/{input.refdir} '  ###~/Transcriptomics_exercise/STAR_Output'
            #'--genomeDir /STAR_Output'
        #    '--runMode alignReads'
            '--readFilesIn ~/Transcriptomics_exercise/{input.R1L1} ~/Transcriptomics_exercise/{input.R2L1} '
            '--readFilesCommand zcat '
            '--outSAMtype BAM SortedByCoordinate '
