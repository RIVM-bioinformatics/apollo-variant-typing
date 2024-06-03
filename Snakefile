import yaml


sample_sheet = config["sample_sheet"]
with open(sample_sheet) as f:
    SAMPLES = yaml.safe_load(f)

for param in ["threads", "mem_gb"]:
    for k in config[param]:
        config[param][k] = int(config[param][k])

# print(SAMPLES)

OUT = config["output_dir"]


def check_if_species_present(samples_dict, genus, species):
    return any(
        [
            samples_dict[sample]["genus"] == genus
            and samples_dict[sample]["species"] == species
            for sample in samples_dict
        ]
    )


localrules:
    all,
    copy_sample_bam,
    copy_ref,
    copy_ref_gff,
    aggregate_species,
    no_typing,
    cauris_extract_aa_mutations,
    combine_auriclas,


include: "workflow/rules/choose_species.smk"
include: "workflow/rules/prepare_files.smk"
include: "workflow/rules/cauris_typing.smk"
include: "workflow/rules/afumigatus_typing.smk"


expected_output = []
expected_output.append(expand(OUT + "/typing_check/{sample}_done.txt", sample=SAMPLES))

if check_if_species_present(SAMPLES, "candida", "auris"):
    expected_output.append(OUT + "/cauris_typing/auriclass.tsv")


rule all:
    input:
        expected_output,
