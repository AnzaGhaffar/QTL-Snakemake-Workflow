configfile: "config.yaml"

rule filtering:
    input:
        expand("Homozygous_Filtered_VCF/Homozygous_{vcf_file_name}.vcf.gz", vcf_file_name=config["vcf_file"])

rule VCF_Homozygous_Filtering:
    input:
        expand("{vcf_file_name}.vcf.gz", vcf_file_name= config["vcf_file"]),
    output:
        expand("Homozygous_Filtered_VCF/Homozygous_{vcf_file_name}.vcf.gz", vcf_file_name=config["vcf_file"])
    params:
        filter_data= config["vcffilter"]["filters"]
    shell:
        "( bcftools view"
        "   {params.filter_data}"
        "   {input}"
        "   -O z"
        "   -o {output}"
        ")"
