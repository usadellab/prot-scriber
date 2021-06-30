# Prot Scriber

Assigns short human readable descriptions to query biological sequences using references. For this, it consumes sequence similarity search (Blast or Diamond) results in tabular format.

## Prot Scriber Experiments

This branch (`evaluation_R`) of the project contains R code that has been used to carry out experiments with the prot-scriber project. Particularly, different methods of assigning human readable descriptions to query sequences have been tried out and their performance assessed using standard prediction scoring methods. Based on the results of these experiments the actual implementation of prot-scriber in Rust has been deviced and designed.

## Installation

The `evaluation_R` branch of prot-scriber is implemented as a R package. A *nix operating system is required for a few of the implemented commands sent to the underlying OS via the `system` R function. To install this package have a look at the `DESCRIPTION` file where the dependencies in form of R packages are listed. Make sure, they are installed. Once this is the case you can install this R package with:

```sh
git clone https://github.com/usadellab/prot-scriber.git ./prot-scriber
cd ./prot-scriber
git checkout evaluation_R
R CMD INSTALL .
```

## Materials

Two sets of query amino acid sequences (proteomes) have been used:

* The _Phaseolus coccineus_ proteome
* The recently published "MetaEuk" [1] proteins

As reference of the molecular functions of the above query proteomes we generated MapMan4 [2] annotations with Mercator [3] and identified conserved protein domains [4] using the Pfam-A database of hidden markov models [5].

For each query protein the words found in these reference annotations were merged into a respective gold standard set of words in lowercase.

To generate the input for the various tested prot-scriber annotation methods sequence similarity searches [8,9] were carried out. The reference databases used were Uniprot Swissprot and trEMBL [6,7].

## Methods and results

### Sequence similarity searches

All sequence similarity searches were carried out using Diamond [9] with the following example command:

```sh
diamond blastp -p 10 --quiet \
  -d uniprot_sprot.fasta.dmnd \
  -q Pcoccineus.faa \
  -o Pcoccineus_vs_swisprot_blastp.txt \
  -f 6 qseqid sseqid qlen qstart qend slen sstart send bitscore stitle evalue;
```
Note that the above example carried out the Diamond / Blast searches for the _P. coccineus_ proteome using the Uniprot Swissprot reference database.

### Identification of reference conserved protein domains

Using hidden markov models and respective search tools conserved protein domains of the Pfam-A database were identified for all query proteins. Consider the following example command:

```sh
hmmscan --cpu 20 \
  --tblout Pcoccineus_vs_PfamA_hmmscan_out.tsv \
  Pfam-A.hmm Pcoccineus.faa;
```
Note that the above command searches the _P. coccineus_ query proteome for conserved domains in the Pfam-A database.

### Identification of MapMan4 reference annotations

All query proteins were used as input to the Mercator annotation tool in order to assign MapMan4 Bin annotations. We used the online web tool for this step [3].

### Assignment and performance evaluation of human readable descriptions (HRDs)

In the following the steps are described that were applied to generate human readable descriptions (HRDs) for the _P. coccineus_ query proteome using and testing different prot-scriber annotation methods. The best Blast (BB) approach was also tested. All resulting HRD annotations were evaluated using two performance scores, the F-Score (with beta parameter set to one) and the Matthew's Correlation Coefficient (MCC). 

The following scripts in the package's `exec` directory were carried out in the given order:

1. loadPcoccineusSeqSimSearchResults.R
2. loadPcoccineus_PfamA_Mercator4_results.R
3. annotateAndEvaluatePcoccineus.R
4. p_coccineus_queries_alignment_regions.R
5. p_coccineus_multi_region_HRDs.R
6. visualize_p_coccineus_results.R

Each script was invoked with the `Rscript` command and the respective arguments as stated in the scripts' manuals.

Intermediate results and scientific plots were generated. 

* `./data/p_coccineus_seq_sim_search.RData` contains the sequence similarity search results for the _P. coccineus_ proteome.
* `./data/p_coccineus_pfamA_Mercator4.RData` contains the Pfam-A annotations (HMMER3 results) and the Mercator (MapMan4) annotations for the _P. coccineus_ proteome.
* `./data/p_coccineus_mercator_pfamA_references.RData` contain the reference word sets obtained from the above two reference function predictions.
* `./data/p_coccineus_reference_words.RData` contain the processed reference word sets, from which uninformative words have been removed.
* `./data/p_coccineus_HRDs.RData` contain the different human readable descriptions assigned to the _P. coccineus_ query proteins by the competing methods. Evaluation scores are included.
* `./data/p_coccineus_alignment_regions.RData` contains information about those query proteins which had disjoint regions in their sequence for which the seq. sim. searches generated local alignments.
* `./data/p_coccineus_HRDs_alignment_regios.RData` contains HRDs generated by providing a HRD for each disjoint local alignment region separately and then concatonate those HRDs. The performance of these concatonated HRDs is also included.

The script `visualize_p_coccineus_results.R` created the following performance score distribution plots:

* p_coccineus_eval_f-scores_dists.pdf
* p_coccineus_eval_MCC_dists.pdf
* p_coccineus_eval_relative_f-scores_dists.pdf
* p_coccineus_eval_relative_MCC_dists.pdf
* p_coccineus_multi_region_HRDs_F_Score_diffs.pdf

and these two HTML documents that show all HRDs that were annotated to the different _P. coccineus_ query proteins given the different HRD annotation methods:

* p_coccineus_HRDs.html.bz2
* p_coccineus_HRDs_multi_region.html.bz2

## References

1. Levy Karin, E., Mirdita, M. & Söding, J. MetaEuk—sensitive, high-throughput gene discovery, and annotation for large-scale eukaryotic metagenomics. Microbiome 8, 48 (2020).
2. Schwacke, R. et al. MapMan4: a refined protein classification and annotation framework applicable to multi-omics data analysis. Molecular Plant (2019) doi:10.1016/j.molp.2019.01.003.
3. Mercator: a fast and simple web server for genome scale functional annotation of plant sequence data - LOHSE - 2013 - Plant, Cell & Environment - Wiley Online Library. http://onlinelibrary.wiley.com/doi/10.1111/pce.12231/full.
4. Eddy, S. R. Accelerated Profile HMM Searches. PLoS Comput. Biol. 7, e1002195 (2011).
5. Bateman, A. et al. The Pfam protein families database. Nucleic Acids Research 32, D138–D141 (2004).
6. Boeckmann, B. et al. The SWISS-PROT protein knowledgebase and its supplement TrEMBL in 2003. Nucleic Acids Res 31, 365–370 (2003).
7. Bairoch, A. & Apweiler, R. The SWISS-PROT protein sequence database and its supplement TrEMBL in 2000. Nucleic Acids Res. 28, 45–48 (2000).
8. McGinnis, S. & Madden, T. L. BLAST: at the core of a powerful and diverse set of sequence analysis tools. Nucleic Acids Res 32, W20–W25 (2004).
9. Buchfink, B., Xie, C. & Huson, D. H. Fast and sensitive protein alignment using DIAMOND. Nat Meth 12, 59–60 (2015).
