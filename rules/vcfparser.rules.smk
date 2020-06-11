configfile: "config.yaml"

rule parser:
    input:
        "QTL_VCF_to_Table/QTL_Table.csv",

rule QTL_VCF_to_Table_Parser:
    input:
        expand("Homozygous_Filtered_VCF/Homozygous_{vcf_file_name}.vcf.gz", vcf_file_name=config["vcf_file"]),
    output:
        "QTL_VCF_to_Table/QTL_Table.csv",
    script:
        "../Scripts/QTL_Parser.R"
