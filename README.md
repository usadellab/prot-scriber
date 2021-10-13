# Prot Scriber

Assigns short human readable descriptions to query biological sequences using references. For this, it consumes sequence similarity search (Blast or Diamond) results in tabular format.

## Implementation

This branch `faba-analysis` is based on the git branch `evaluation_R`.

## Annotation for Faba

This branch `faba-analysis` is just a temporary solution to enable annotation of the Faba query proteome with the current `R` implementation. At the time of writing and using this branch on the Faba query proteins, the `Rust` version of prot-scriber was not yet finished.

### Work protocol

The above "abuse" of the `R` code implies that the protocol used to assign short human readable descriptions to the Faba query proteins is something of a "hack" - apologies.

### Directory on the IBG-4 cluster

Find all data and scripts of this annotation run in directory `/mnt/data/asis/prot-scriber/20210618_Faba_annotations_proteins_v1.0`. Results and intermediates are stored in the sub-folder `./prot-scriber/data`.

### Step 1 - Sequence similarity searches

Diamond was used to search for significant homologs in the UniProt SwissProt and trEMBL public reference databases. This was done using the scripts:

* `Faba_vs_swissprot_oge_job.sh` and
* `Faba_vs_trembl_oge_job.sh`

### Step 2 - Load and parse the Blastp results

The respective `R` script, originally intended for the _P. coccineus_ annotation was (mis)-used. The script did not need any adjustment.
```sh
Rscript ./exec/loadPcoccineusSeqSimSearchResults.R -s ../Faba_vs_swissprot_blastp_out.txt -t ../Faba_vs_trembl_blastp_out.txt -n 40 -d ./data
```

### Step 3 - Assign short human readable descriptions to the Faba query proteome

The respective `R` script, originally intended for the _P. coccineus_ annotation was (mis)-used. Note, that the script was adjusted to work for Faba query proteins. 
```sh
Rscript ./exec/annotateAndEvaluatePcoccineus.R -d ./data/ -n 40
```

### Step 4 - Generate the output table

After step 3 has finished, there will be a result binary file `./data/p_coccineus_HRDs.RData`. Note, that even though the file says "p_coccineus" that does not mean the annotations are for the query proteome of _P. coccineus_, they are those generated for the Faba queryome.

In order to generate a simple table in which the query protein identifiers are mapped to the assigned short human readable description, the following `R` code was executed in an interactive `R` shell. Note, that the working directory _must_ be `/mnt/data/asis/prot-scriber/20210618_Faba_annotations_proteins_v1.0/prot-scriber`.

```r
load('./data/p_coccineus_HRDs.RData')
x <- "centered.inverse.inf.cntnt.mean"
write.table( pc.hrds[ pc.hrds$Method == x, c( 'Protein.ID', 'HRD' )],
  './data/faba_human_readable_descriptions.txt', sep="\t", row.names=FALSE,
  quote=F )
```

### Results

Find the result table in file `./data/faba_human_readable_descriptions.txt`.
