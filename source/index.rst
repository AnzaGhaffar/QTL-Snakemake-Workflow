.. _manual-main:

=========
QTL-Snakemake-Workflow
=========

.. image:: https://img.shields.io/conda/dn/bioconda/snakemake.svg?label=Bioconda
    :target: https://bioconda.github.io/recipes/snakemake/README.html

.. image:: https://img.shields.io/pypi/pyversions/snakemake.svg
    :target: https://www.python.org

.. image:: https://img.shields.io/pypi/v/snakemake.svg
    :target: https://pypi.python.org/pypi/snakemake

.. image:: https://img.shields.io/docker/cloud/build/snakemake/snakemake
       :target: https://hub.docker.com/r/snakemake/snakemake

.. image:: https://github.com/snakemake/snakemake/workflows/CI/badge.svg?branch=master
    :target: https://github.com/snakemake/snakemake/actions?query=branch%3Amaster+workflow%3ACI

.. image:: https://img.shields.io/badge/stack-overflow-orange.svg
    :target: https://stackoverflow.com/questions/tagged/snakemake

.. image:: https://img.shields.io/twitter/follow/johanneskoester.svg?style=social&label=Follow
    :target: https://twitter.com/search?l=&q=%23snakemake%20from%3Ajohanneskoester


.. .. raw:: html
          <span class="__dimensions_badge_embed__" data-doi="https://doi.org/10.1093/bioinformatics/bts480" data-legend="always" data-style="large_rectangle"></span><script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>

The QTL-Snakemake-Workflow consists of 3 main fuctions which helps to utilize
the QTLseqr tool for Next Generation sequencing bulk segregant analysis from 
a VCF file.

.. _main-getting-started:

---------------
Getting started
---------------

To run this package first we need to install the required libraries for the workflow which can be installed by running 

.. code-block:: python
    
    conda env create --name QTL_Test --file condaenv.yml

Here we can define any name we want instead of the QTL_Test. This command will make a conda environment to run this package.

Once we have created the environment we need to configure the **config.yaml** file with the required parameters so that it can take input 
for different rules. After this to run the workflow we can use the command

.. code-block:: python

    snakemake -pr

This command will run the workflow and generate the **output** in the **QTL_Plots** folder.


.. _manual-Work_Flow_Rules:

-------------
Work FLow Rules
-------------

The QTL workflow consists of three rules:
 1. VCF_Homozygous_Filtering
 2. QTL_VCF_to_Table_Parser
 3. QTL_Plotting

1. VCF_Homozygous_Filtering

This rule takes any vcf file and filters the parent for selecting only the homozygous SNPS in the vcf file.

.. code-block:: python

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

The **input** is a vcf file which needs to be defined in the config file under the **vcf_file** name and the **output**
of this rule will be saved into the **Homozygous_Filtered_VCF** folder with the same vcf file name just a **Homozygous**
title will be added in the start of the name so that we can distinguish between different files if the workflow is run 
for multiple vcf files. The filtering is done on the parents and the shell command takes input from the config file.

.. code-block:: python

     -i "GT[0]='0/0' && GT[1]='1/1' || GT[0]='1/1' && GT[1]='0/0'"

In the above script for filtering the parents which comes from the config file it is taking the parents and from our reference 
vcf file the parents are on **GT[0]** and **GT[1]**, but it should be adjusted according to the given vcf file.

2. QTL_VCF_to_Table_Parser

This rule is for converting a vcf file into a table format file which is required by the QTL tool if the VCF file is not generated
from **GATK**. This rule runs a R scripts which parses the vcf file.

.. code-block:: python

    input:
        expand("Homozygous_Filtered_VCF/Homozygous_{vcf_file_name}.vcf.gz", vcf_file_name=config["vcf_file"]),
    output:
        "QTL_VCF_to_Table/QTL_Table.csv",
    script:
        "Scripts/QTL_Parser.R"

It takes an **input** a vcf file from the **Homozygous_Filtered_VCF** folder with a specified name as **Homozygous_{vcf_file_name}.vcf.gz**
the **vcf_file_name** comes from the config file where the name of the vf file was given. The **output** is saved in **QTL_VCF_to_Table** folder
The **QTL_VCF_to_Table** script takes input from the config file where the names of the **High Bulk, Low Bulk** are defined along with that it also 
takes the parameter **Number_of_Chromosomes** from the config file.

3. QTL_Plotting

This rule runs the QTLseqr tool and generates the plots for **Gprime Analysis** and **QTLseq Analysis**. This rules runs an R script for generating the plots.

.. code-block:: python

    input:
       "QTL_VCF_to_Table/QTL_Table.csv"
    output:
       "QTL_Plots/DP_Filtering Data.pdf",
       "QTL_Plots/REF Frequency for Filtering Data.pdf",
       "QTL_Plots/SNP Index for Filtering Data.pdf",
       "QTL_Plots/GPrime Distribution with Hampel Outlier Filter.pdf",
       "QTL_Plots/GPrime Distribution with deltaSNP Outlier Filter.pdf",
       "QTL_Plots/SNP Density Plot.pdf",
       "QTL_Plots/Delta SNP Index Plot with Intervals.pdf",
       "QTL_Plots/GPrime Value Plot.pdf"       
    script:
        "Scripts/QTL_Plotting.R"

It takes **input** the csv file developed by the **QTL_VCF_to_Table_Parser** along with the parameters defined within the config file for filtering the SNPs
for better results and the **output** is saved into **QTL_Plots** folder.

 




