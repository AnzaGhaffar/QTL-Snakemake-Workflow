configfile: "config.yaml"

rule computomics:
    input:
       expand("Allel_Frequency_Plots_Computomics/{vcf_file_name}.tsv", vcf_file_name=config["AllelPlots"]["vcffile"]),
       expand("Allel_Frequency_Plots_Computomics/{vcf_file_name}.pdf", vcf_file_name=config["AllelPlots"]["vcffile"]),
       expand("Allel_Frequency_Plots_Computomics/{vcf_file_name}_window.pdf", vcf_file_name=config["AllelPlots"]["vcffile"])
#       "Allel_Frequency_Plots_Computomics/freebayes_allele_freq_D2.pdf",
#       "Allel_Frequency_Plots_Computomics/freebayes_allele_freq_window_D2.pdf",




rule Allel_Frequency_Tsv_Generator:
    input:
        expand("{vcf_file_name}.{ext}",vcf_file_name= config["AllelPlots"]["vcffile"],ext=config["AllelPlots"]["ext"])
    output:
        expand("Allel_Frequency_Plots_Computomics/{vcf_file_name}.tsv", vcf_file_name=config["AllelPlots"]["vcffile"])
    params:
        samples=config["AllelPlots"]["Samples"],
        data=config["AllelPlots"]["data"]
    shell:
        """
        ( printf "#Samples:con-all,D2,D2_F2_tt,D2_F2_TT\nCHROM\\tPOS\\tREF\\tALT\\tRO\\tAO\\tGT\\tGQ\\tSampleRO\\tSampleAO\n" > {output} 
          bcftools view {input} {params.samples} | bcftools query {params.data} >> {output})"""


rule Allel_Plots:
    input:
        expand("Allel_Frequency_Plots_Computomics/{vcf_file_name}.tsv", vcf_file_name=config["AllelPlots"]["vcffile"])
    output:
        expand("Allel_Frequency_Plots_Computomics/{vcf_file_name}.pdf", vcf_file_name=config["AllelPlots"]["vcffile"]),
        expand("Allel_Frequency_Plots_Computomics/{vcf_file_name}_window.pdf", vcf_file_name=config["AllelPlots"]["vcffile"])
#        "Allel_Frequency_Plots_Computomics/freebayes_allele_freq_D2.pdf",
#        "Allel_Frequency_Plots_Computomics/freebayes_allele_freq_window_D2.pdf"
    script:
        "../Scripts/plot_allele_freqs.py"
