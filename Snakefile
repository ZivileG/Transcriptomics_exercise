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
