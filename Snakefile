rule fastqc:
    input:
        "input/{sample}.fastq.gz"
    output:
        "output/fastqc_out/{sample}_fastqc.html"
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
        "input/{sample}.fastq.gz"
    output:
        "output/bbduk_out/{sample}_trimmed.fastq.gz"
    conda:
        "environment.yaml"
    shell:
        "bbduk.sh in={input} out={output} ref=input/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10"

rule fastqc_trimmed:
    input:
        "output/bbduk_out/{sample}.fastq.gz"
    output:
        "output/fastqc_trim_out/{sample}_fastqc.html"
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
        directory('STAR_Output')
    threads: 20
    shell:
        'mkdir {output} && '
        'STAR --runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {output} '
        '--genomeFastaFiles {input.fa} '
        '--sjdbGTFfile {input.gtf} '
        '--sjdbOverhang 100'
rule pass1:
        input:
            R1L1 = '~/Transcriptomics_exercise/output/bbduk_out/{sample1}_L001_R1_001_trimmed.fastq.gz', # may need adjustment if your fastq file name format is different
            R2L1 = '~/Transcriptomics_exercise/output/bbduk_out/{sample1}_L001_R2_001_trimmed.fastq.gz',
        ##    refdir = directory('~/Transcriptomics_exercise/STAR_Output')
        params:
            outdir = '~/Transcriptomics_exercise/output/pass1_out/{sample1}_pass1',
            ##rmbam = '~/Transcriptomics_exercise/output/pass1_out/{sample}_pass1/Aligned.out.bam'
        output:
            'output/pass1_out/{sample1}_pass1_out/Aligned.sortedByCoord.out.bam'
        threads: 20
        shell:
            'rm -rf {params.outdir} &&' # be careful with this. I don't know why, but Snakemake had problems without this cleaning.
            'mkdir {params.outdir} && ' # snakemake had problems finding output files with --outFileNamePrefix, so I used this approach instead
            'cd {params.outdir} && '
            'STAR --runThreadN {threads} '
            '--genomeDir ~/Transcriptomics_exercise/STAR_Output '
            '--readFilesIn {input.R1L1} {input.R2L1} '
            '--readFilesCommand zcat '
            '--outSAMtype BAM SortedByCoordinate '## && rm {params.rmbam} && cd ..'
