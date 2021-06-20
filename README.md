# ScerSeg - Segmentation of S. cerevisiae genome

<table border="0"><tr>

<td>
This repository serves as an accompanying resource to my master thesis, in which we segment the budding yeast genome into areas of local similarity.

[One page summary](https://github.com/mgeorgoulopoulos/ScerSeg/blob/main/OnePageSummary.docx) is a word file describing in short this work's methodology and findings.
</td>

<td width="40%">
<img src="https://raw.githubusercontent.com/mgeorgoulopoulos/ScerSeg/main/r/images/PromoterFields.png"/>
</td>

</tr></table>

Both the [thesis document](https://github.com/mgeorgoulopoulos/ScerSeg/blob/main/r/Final.html) and [full technical reference](https://github.com/mgeorgoulopoulos/ScerSeg/blob/main/r/chapters/index.html) are data-driven and automatically generated from the data. I highly recommend reading these, instead of the word document of my thesis, as these have interactive 3D plots, allowing you to explore the yeast genome in 3D. To view these files, follow the corresponding links and download the html files to your system. Then open locally.


# Potentially useful resources

Feel free to browse the [Results](https://github.com/mgeorgoulopoulos/ScerSeg/tree/main/Results) directory for comma- and tab-separated files containing my results. The [GeneSets](https://github.com/mgeorgoulopoulos/ScerSeg/tree/main/Results/GeneSets) subdirectory contains the resulting gene sets along with a brief description for each.

[ScerGenePositions.tsv](https://github.com/mgeorgoulopoulos/ScerSeg/blob/main/Results/ScerGenePositions.tsv) is a tab-separated file containing x,y,z coordinates of all S. cerevisiae genes, corresponding to (Duan et al. , 2010) three-dimensional model of the yeast genome. Added here with Dr. Duan's permission.


In [src](https://github.com/mgeorgoulopoulos/ScerSeg/tree/main/src), one can find a collection of C++ programs. Most of these implement the "sphere-test" algorithm outlined in the thesis document and can serve as a basis for adapting the algorithm to new data.  [CsvPatcher](https://github.com/mgeorgoulopoulos/ScerSeg/blob/main/src/CsvPatcher.cpp) is a simple program that can fill missing values in CSV files, based on the surrounding cells.

# Acknowledgments


Dr. Christoforos Nikolaou is the supervisor of this work. He proposed the sphere-sampling algorithm and provided direction and prioritization of the various tests we carried out, as well as a plethora of insights on statistics, genome organization and writing. His contribution can not be overstated.

Athanasia Stavropoulou offered the transcription factor motif signal and insight regarding the resulting motif clusters, coming from her own work of training a hidden Markov model on the same dataset. In addition, she helped during debugging of parsing the coexpression dataset and took an active role in the creative process by both participating in the yeast-related weekly meetings and by constantly exchanging ideas and advice.

Nikos Vakirlis has provided the gene species count and overarching taxon per gene dataset.

CG2 team has contributed by critically discussing every aspact of this work.
	
