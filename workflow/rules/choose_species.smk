def choose_species(wildcards):
    if (SAMPLES[wildcards.sample]["genus"] == "aspergillus") & (
        SAMPLES[wildcards.sample]["species"] == "fumigatus"
    ):
        return [
            OUT + "/afumigatus/annotated_variants/{sample}.vcf",
            OUT + "/afumigatus/resistance_mutations/{sample}.tsv",
        ]
    elif (SAMPLES[wildcards.sample]["genus"] == "candida") & (
        SAMPLES[wildcards.sample]["species"] == "auris"
    ):
        return [
            OUT + "/cauris_typing/annotated_vcf/{sample}.vcf",
            OUT + "/cauris_typing/resistance_mutations/{sample}.tsv",
            OUT + "/cauris_typing/auriclass/{sample}.tsv",
        ]
    else:
        return OUT + "/typing_check/{sample}/no_typing_necessary.txt"


rule aggregate_species:
    input:
        choose_species,
    output:
        temp(OUT + "/typing_check/{sample}_done.txt"),
    message:
        "Checking correct typing ran for {wildcards.sample}"
    threads: 1
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        "touch {output}"


rule no_typing:
    input:
        vcf=lambda wildcards: SAMPLES[wildcards.sample]["vcf"],
    output:
        temp(OUT + "/typing_check/{sample}/no_typing_necessary.txt"),
    message:
        "Skipping typing step for {wildcards.sample}."
    threads: 1
    shell:
        """
        touch {output}
        """
