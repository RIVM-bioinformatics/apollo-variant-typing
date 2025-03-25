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


rule afumigatus_compare_aa_mutations:
    input:
        tsv=OUT + "/afumigatus_typing/annotated_variants/{sample}.tsv",
        aa_resistance_variants_csv=lambda wildcards: SAMPLES[wildcards.sample][
            "aa_resistance_variants_csv"
        ],
    output:
        tsv=OUT + "/afumigatus_typing/resistance_mutations/aa/{sample}.aa.tsv",
        full=OUT + "/afumigatus_typing/resistance_mutations/aa/{sample}.aa.full.tsv",
    message:
        "Extract AMR mutations (amino acid based) for {wildcards.sample}"
    resources:
        mem_gb=config["mem_gb"]["compare"]
    log:
        OUT + "/log/afumigatus_compare_aa_mutations/{sample}.log",
    shell:
        """
python workflow/scripts/compare_aa_mutations.py \
    --input {input.tsv} \
    --output {output.tsv} \
    --full-output {output.full} \
    --resistance_variants_csv {input.aa_resistance_variants_csv}
        """


rule afumigatus_compare_nt_mutations:
    input:
        tsv=OUT + "/afumigatus_typing/annotated_variants/{sample}.tsv",
        nt_resistance_variants_csv=lambda wildcards: SAMPLES[wildcards.sample][
            "nt_resistance_variants_csv"
        ],
    output:
        tsv=OUT + "/afumigatus_typing/resistance_mutations/nt/{sample}.nt.tsv",
    message:
        "Extract AMR mutations (nucleotide based) for {wildcards.sample}"
    resources:
        mem_gb=config["mem_gb"]["compare"]
    log:
        OUT + "/log/afumigatus_compare_nt_mutations/{sample}.log",
    shell:
        """
python workflow/scripts/compare_nt_mutations.py \
    --input {input.tsv} \
    --output {output.tsv} \
    --resistance_variants_csv {input.nt_resistance_variants_csv}
        """


rule afumigatus_combine_aa_nt_mutations:
    input:
        aa=OUT + "/afumigatus_typing/resistance_mutations/aa/{sample}.aa.tsv",
        nt=OUT + "/afumigatus_typing/resistance_mutations/nt/{sample}.nt.tsv",
        aa_full=OUT + "/afumigatus_typing/resistance_mutations/aa/{sample}.aa.full.tsv",
    output:
        tsv=OUT + "/afumigatus_typing/resistance_mutations/{sample}.combined.tsv",
        full=OUT + "/afumigatus_typing/resistance_mutations/{sample}.combined.full.tsv",
    message:
        "Combine AMR mutations (amino acid and nucleotide based) for {wildcards.sample}"
    resources:
        mem_gb=config["mem_gb"]["compare"]
    log:
        OUT + "/log/afumigatus_combine_aa_nt_mutations/{sample}.log",
    shell:
        """
        python workflow/scripts/combine_aa_nt_reports.py -aa {input.aa} -nt {input.nt} -o {output.tsv}
        python workflow/scripts/combine_aa_nt_reports.py -aa {input.aa_full} -nt {input.nt} -o {output.full}
        """
