# Gene sets

This directory contains the gene sets we found to be locally similar. All files are tab-separated and of the form "Gene name" TAB "Set identifier".

|File                       |Description                                                                                                                                                                             |
|---------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|CoexFields.tsv             |Raw output of our algorithm, extracting clusters of significantly high local coexpression.                                                                                              |
|Continents.tsv             |A curated version of CoexFields.tsv, also including the low-coexpression compartment and naming the compartments using continent names.                                                 |
|CommunityCenters.tsv       |List of genes forming a local community of 5 genes around them (gene +/- 2 neighboring genes on the same chromosome) assigned to one of the 7 histone profile classes                   |
|MotifFields.tsv            |Clusters of similar transcription factor binding sites. Extracted by Jaccard distance.                                                                                                  |
|JaccardIndexMotifFields.tsv|A refined version of MotifFields utilizing the Jaccard index instead of distance. This results in greater statistical significance. This recent addition is not documented in the paper.|
|PromoterFields.tsv         |Clusters of locally similar promoter-site histone modification profiles.                                                                                                                |
|RoughFeatures.tsv          |Arbitrary compartmentalization of the yeast genome, designating "Hat" (nucleolus), "Internal" and "External" regions.                                                                   |
|TightCommunities.tsv       |Chromosomal segments containing a small number of genes each, having a low Shannon entropy in terms of histone modification profiles.                                                   |

