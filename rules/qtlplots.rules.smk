configfile: "config.yaml"

rule plotting:
    input:
       "QTL_Plots/DP_Filtering Data.pdf",
       "QTL_Plots/REF Frequency for Filtering Data.pdf",
       "QTL_Plots/SNP Index for Filtering Data.pdf",
       "QTL_Plots/GPrime Distribution with Hampel Outlier Filter.pdf",
       "QTL_Plots/GPrime Distribution with deltaSNP Outlier Filter.pdf",
       "QTL_Plots/SNP Density Plot.pdf",
       "QTL_Plots/Delta SNP Index Plot with Intervals.pdf",
       "QTL_Plots/GPrime Value Plot.pdf",

rule QTL_Plotting:
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
       "QTL_Plots/GPrime Value Plot.pdf",
    script:
       "../Scripts/QTL_Plotting.R"
