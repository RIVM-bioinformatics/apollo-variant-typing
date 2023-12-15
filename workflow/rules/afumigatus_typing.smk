rule afumigatus_annotate_vcf:
    input:
        vcf=lambda wildcards: SAMPLES[wildcards.sample]["vcf"],
        gff_ref=OUT + "/prepared_files/{sample}_ref.gff",
        fasta_ref=OUT + "/prepared_files/{sample}_ref.fasta",
    output:
        vcf=OUT + "/afumigatus_typing/annotated_vcf/{sample}.vcf",
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


rule afumigatus_annotated_vcf_to_table:
    input:
        vcf=OUT + "/afumigatus_typing/annotated_vcf/{sample}.vcf",
    output:
        tsv=OUT + "/afumigatus_typing/annotated_variants/{sample}.tsv",
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
        OUT + "/log/afumigatus_annotated_vcf_to_table/{sample}.log",
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
