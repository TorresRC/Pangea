# Pangea
Pangea: The Pan-Genome Analysis Pipeline

Pangea is a pipeline to calculate the entire gene collection of a set of prokaryotic organisms (Pan genome). It is based on an all-against-all algorithm to find shared and non-shared genes between samples. 

Pangea was designed to compare high diverse prokaryotic genomes using profile Hidden Markov Model (profile HMM) to make a position-specific scoring. Pangea builds a profile HMM for each gene that constitutes the pan genome. Then each profile is used to search the corresponding gene in each sample. During this search the profiles are iteratively populated with the corresponding homologous sequences. This allows each profile HMM to include all the variations that a gene presents between samples, increasing the robustness and therefore the accuracy of the search.

First, Pangea makes a single blast comparison (blastn for nucleotides or blastp for amino acids) between the features of the first two query samples or between the first query and a reference sample. Homologous sequences, according to an established e-value and a percentage of identity, are considered shared features, while features without any blast hit are considered non-shared. Sequences of shared features are extracted and aligned with muscle to build a profile HMM from the alignment using hmmbuil, while sequences of shared features are extracted and used to build profile HMM from a single sequence. Then, the profiles HMM are used to look for homologous features in the remaining samples according to the same percentage of identity and e-value used in the blast comparison. Only the best hit is considered homolog (shared feature) and its sequence is added to the corresponding profile HMM, while non-homologous features are re-analyzed to discard duplicates or fragmented sequences, which would generate less significant hits. Remaining non-homologous features are considered new elements of the pan genome and their sequence are used to create profiles HMM of a single sequence. All features are reported in a presence/absence table, which represents the composition of the entire pan genome, including the core genome, accessory genome and singletons.

## Installation
Pangea is a collection of perl scripts and it has been tested on Linux and MacOS. For its correct execution requires the previous installation of some third-party programs.

### Dependencies 
#### Third-party:
- [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [Muscle](https://www.drive5.com/muscle/)
- [HMMER](http://hmmer.org)
- [R](https://www.r-project.org)

#### Perl modules:
- [FindBin](https://perldoc.perl.org/FindBin.html)
- [List::MoreUtils](https://metacpan.org/pod/release/REHSACK/List-MoreUtils-0.416/lib/List/MoreUtils.pm)
- [Getopt::Long](https://perldoc.perl.org/Getopt/Long.html)

#### R packages:
- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
- [gplots](https://cran.r-project.org/web/packages/gplots/index.html)


## Input files
For now, Pangea works only with the outputs produced by [Prokka](https://github.com/tseemann/prokka/blob/master/README.md). All the annotated samples should be contained in the same directory, one folder per sample. It is highly recommended to use the --locustag option during the annotation process with prokka. For more information on the installation and execution of prokka check [here](https://github.com/tseemann/prokka/blob/master/README.md).

You have to create a text file that lists all the samples you want to analyze. The samples will be analyzed according to the order specified in the list.

Additionally you can include a reference sample. This can be declared at the beginning of your sample list or you can use the --trusted option to include it. See the Usage section.

## Outputs
Depending on the used options, the possible outputs include the following:

|Name|Type|Description|
|-|-|-|
|Blast|Dir|Contains all the blast databases used by the program|
|PanGenomeHmmDb|Dir|Contains a database of the pan genome used during the HMMER comparisons|
|ORFeomes|Dir|Contains the ORFeomes or Proteomes of the queries|
|ORFs|Dir|Contains all the analyzed features. One folder per feature|
|CoreSequences|Dir|Contains the core and the aligned-core fasta files of each query|
|\*UniqueBlastComparison.txt|File|Contains the blast report of the first comparison|




## Citation

Seemann T.  
*Prokka: rapid prokaryotic genome annotation*  
**Bioinformatics** 2014 Jul 15;30(14):2068-9.   
[PMID:24642063](http://www.ncbi.nlm.nih.gov/pubmed/24642063)





Website: https://torresrc.github.io/Pangea/
