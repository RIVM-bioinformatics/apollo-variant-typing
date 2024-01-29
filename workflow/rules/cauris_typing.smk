rule cauris_annotate_vcf:
    input:
        vcf=lambda wildcards: SAMPLES[wildcards.sample]["vcf"],
        gff_ref=OUT + "/prepared_files/{sample}_ref.gff",
        fasta_ref=OUT + "/prepared_files/{sample}_ref.fasta",
    output:
        vcf=OUT + "/cauris_typing/annotated_vcf/{sample}.vcf",
    message:
        "Annotate VCF for {wildcards.sample}"
    container:
        "docker://staphb/bcftools:1.18"
    conda:
        "../envs/bcftools.yaml"
    threads: config["threads"]["bcftools"]
    resources:
        mem_gb=config["mem_gb"]["bcftools"],
    log:
        OUT + "/log/bcftools_csq/{sample}.log",
    shell:
        """
bcftools csq \
    --phase a \
    -f {input.fasta_ref} \
    -g {input.gff_ref} \
    {input.vcf} \
    > {output.vcf} \
    2>{log}
        """


rule cauris_annotated_vcf_to_table:
    input:
        vcf=OUT + "/cauris_typing/annotated_vcf/{sample}.vcf",
    output:
        tsv=OUT + "/cauris_typing/annotated_variants/{sample}.tsv",
    message:
        "Convert annotated variants to table for {wildcards.sample}"
    container:
        "docker://broadinstitute/gatk:4.3.0.0"
    conda:
        "../envs/gatk_picard.yaml"
    threads: config["threads"]["gatk"]
    resources:
        mem_gb=config["mem_gb"]["gatk"],
    log:
        OUT + "/log/cauris_annotated_vcf_to_table/{sample}.log",
    shell:
        """
gatk VariantsToTable \
-V {input.vcf} \
-F CHROM \
-F POS \
-F TYPE \
-F REF \
-F ALT \
-F DP \
-F AF \
-F BCSQ \
-O {output.tsv} 2>&1>{log}
        """


rule cauris_extract_aa_mutations:
    input:
        tsv=OUT + "/cauris_typing/annotated_variants/{sample}.tsv",
        aa_resistance_variants_csv=lambda wildcards: SAMPLES[wildcards.sample][
            "aa_resistance_variants_csv"
        ],
    output:
        tsv=OUT + "/cauris_typing/resistance_mutations/{sample}.tsv",
        full=OUT + "/cauris_typing/resistance_mutations/{sample}.full.tsv",
    message:
        "Extract AMR mutations for {wildcards.sample}"
    log:
        OUT + "/log/cauris_extract_amr_mutations/{sample}.log",
    shell:
        """
python workflow/scripts/extract_amr_mutations.py \
    --input {input.tsv} \
    --output {output.tsv} \
    --full-output {output.full} \
    --resistance_variants_csv {input.aa_resistance_variants_csv}
        """


rule cauris_bam_to_fastq:
    input:
        bam=lambda wildcards: SAMPLES[wildcards.sample]["bam"],
    output:
        r1=temp(OUT + "/fastq/{sample}.R1.fastq"),
        r2=temp(OUT + "/fastq/{sample}.R2.fastq"),
    message:
        "Convert {input.bam} to fastq for {wildcards.sample}"
    container:
        "docker://broadinstitute/picard:2.27.5"
    conda:
        "../envs/gatk_picard.yaml"
    threads: config["threads"]["picard"]
    resources:
        mem_gb=config["mem_gb"]["picard"],
    log:
        OUT + "/log/cauris_bam_to_fastq/{sample}.log",
    shell:
        """
java -jar /usr/picard/picard.jar SamToFastq \
    --INPUT {input.bam}\
    --FASTQ {output.r1} \
    --SECOND_END_FASTQ {output.r2} \
    2> {log}
        """


rule cauris_auriclass:
    input:
        r1=OUT + "/fastq/{sample}.R1.fastq",
    output:
        OUT + "/cauris_typing/auriclass/{sample}.tsv",
    message:
        "Run auriclass for {wildcards.sample}"
    container:
        "docker://quay.io/biocontainers/auriclass:0.5.3--pyhdfd78af_0"
    conda:
        "../envs/auriclass.yaml"
    threads: config["threads"]["auriclass"]
    resources:
        mem_gb=config["mem_gb"]["auriclass"],
    params:
        name="{sample}",
    log:
        OUT + "/log/cauris_auriclass/{sample}.log",
    shell:
        """
auriclass \
    -o {output} \
    -n {params.name} \
    {input.r1}
        """


rule combine_auriclas:
    input:
        expand(OUT + "/cauris_typing/auriclass/{sample}.tsv", sample=SAMPLES),
    output:
        OUT + "/cauris_typing/auriclass.tsv",
    message:
        "Combine auriclass results"
    log:
        OUT + "/log/combine_auriclas.log",
    shell:
        """
cat <(head -n 1 {input[0]}) \
    <(for file in {input}; do tail -n +2 $file; done) \
    > {output} \
    2> {log}
        """
