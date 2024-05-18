[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=prot_scriber)

# prot-scriber

Assigns short human readable descriptions (HRD) to query biological sequences using reference candidate descriptions. In this, `prot-scriber` consumes sequence similarity search (Blast or Diamond or similar) results in tabular format. A customized lexical analysis is carried out on the descriptions of these Blast Hits and a resulting HRD is assigned to the query sequences. 

`prot-scriber` can also apply the same methodology to produce HRDs for sets of biological sequences, i.e. gene families. 

## Quick Start

(This section is for the lazy :wink: [TL;DR](https://en.wikipedia.org/wiki/Wikipedia:Too_long;_didn%27t_read))

`prot-scriber` can be used to generate human readable descriptions (HRDs) for _either_ query biological _sequences_ (proteins) or for _gene-families_. We will walk you through both use cases below with ready to use example input files.

### Step 1 - Get `prot-scriber`

Depending on your operating system, download the ready to use executable from the table in section [\"Installation\"](#download-ready-to-use-executables).

### Step 2 - Sequence similarity searches with Blast or Diamond

Independent of your use-case, sequence or gene-family annotation, you need to run a sequence similarity search of your query biological sequences against reference databases. We recommend searching [UniProt](https://www.uniprot.org/downloads) Swissprot ([uniprot_sprot.fasta.gz](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz)) and trEMBL ([uniprot_trembl.fasta.gz](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz)). You will need to format (Blast `makeblastdb`, Diamond `diamond makedb`) these UniProt reference databases and search them with either Blast 
```sh
blastp -db uniprot_sprot.fasta -query my_prots.fasta -num_threads 10 -out my_prots_vs_sprot.txt -outfmt \"6 delim=<TAB> qacc sacc stitle\"
```
(Note that the above `<TAB>` actually needs to be a tab character. Typically you type that in with \"Ctrl+v\" followed by \"Tab\".)

or Diamond 
```sh
diamond blastp -p 10 --quiet -d uniprot_sprot.fasta.dmnd -q my_prots.fasta -o my_prots_vs_sprot.txt -f 6 qseqid sseqid stitle
```
(See [the manual](#manual) section \"2.2 Example Blast or Diamond commands\" for details). For a quick test run you can assume to have carried out the searches and use the example output tables below (all files are taken from this repository's [`misc`](https://github.com/usadellab/prot-scriber/tree/master/misc) directory):

**To generate HRDs for twelve example biological sequences (proteins) use:**
* [`Twelve_Proteins_vs_Swissprot_blastp.txt`](https://raw.githubusercontent.com/usadellab/prot-scriber/master/misc/Twelve_Proteins_vs_Swissprot_blastp.txt)
* [`Twelve_Proteins_vs_trembl_blastp.txt`](https://raw.githubusercontent.com/usadellab/prot-scriber/master/misc/Twelve_Proteins_vs_trembl_blastp.txt)

**To generate HRDs for two gene-families, comprising four and three proteins, respectively, use:**
* [`families.txt`](https://raw.githubusercontent.com/usadellab/prot-scriber/master/misc/families.txt)
* [`family_prots_vs_Swissprot.txt`](https://raw.githubusercontent.com/usadellab/prot-scriber/master/misc/family_prots_vs_Swissprot.txt)
* [`family_prots_vs_trEMBL.txt`](https://raw.githubusercontent.com/usadellab/prot-scriber/master/misc/family_prots_vs_trEMBL.txt)

Please read section \"2.3 Gene Family preparation and analysis\" of [the manual](#manual) for a recipy on how to cluster biological sequences into gene-families.

### Step 3 - Assign human readable descriptions (HRDs)

Note that on Windows the below would be \"`prot-scriber.exe`\" instead of just \"`prot-scriber`\".

**To annotate biological sequences, e.g. proteins, with HRDs, use:**

```sh
prot-scriber -s Twelve_Proteins_vs_Swissprot_blastp.txt -s Twelve_Proteins_vs_trembl_blastp.txt -o Twelve_Proteins_HRDs.txt
```
Find `prot-scriber`'s output in file `Twelve_Proteins_HRDs.txt`.

**To annotate gene-families with HRDs, use:**
```sh
prot-scriber -f families.txt -s family_prots_vs_Swissprot.txt -s family_prots_vs_trEMBL.txt -o families_HRDs.txt
```
Find `prot-scriber`'s output in file `families_HRDs.txt`.

## Installation

### Download ready to use executables

You can choose to download a pre-built binary, ready to be executed, from the table below, if you want the latest stable version. Other versions can be downloaded from the [Releases page](https://github.com/usadellab/prot-scriber/releases). Have a look at the below table to know which is the version you need for your operating system and platform. 

We strongly recommend to _rename_ the downloaded release file to a simple `prot-scriber` (or `prot-scriber.exe` on Windows).

Note that on Mac OS and Unix / Linux operating systems you need to make the downloaded and renamed `prot-scriber` file _executable_. In order to achieve this, open a terminal shell, navigate (`cd`) to the directory where you saved `prot-scriber`, and execute `chmod a+x ./prot-scriber`.


|Operating System|CPU-Architecture|Release-Name (click to download)|Comment|
|---|---|---|---|
|Windows 7 or higher|any|[windows_prot-scriber.exe](https://github.com/usadellab/prot-scriber/releases/download/latest-stable/x86_64-pc-windows-gnu_prot-scriber.exe)|to be used in a terminal (`cmd` or Power-Shell)|
|any GNU-Linux|any Intel x86, 64 bits|[x86_64-unknown-linux-gnu_prot-scriber](https://github.com/usadellab/prot-scriber/releases/download/latest-stable/x86_64-unknown-linux-gnu_prot-scriber)|requires libm.so.6 (compiled with glibc version 2.27) and libc.so.6 (compiled with glibc 2.18) installed as is the case e.g. in Ubuntu >= 22.04|
|any GNU-Linux|any aarch, 64 bits|[aarch64-unknown-linux-gnu_prot-scriber](https://github.com/usadellab/prot-scriber/releases/download/latest-stable/aarch64-unknown-linux-gnu_prot-scriber)|e.g. for Raspberry Pi; requires libm.so.6 (compiled with glibc version 2.27) and libc.so.6 (compiled with glibc 2.18) installed as is the case e.g. in Ubuntu >= 22.04|
|Apple / Mac OS|any Mac Computer with Mac OS 10|[x86_64-apple-darwin_prot-scriber](https://github.com/usadellab/prot-scriber/releases/download/latest-stable/x86_64-apple-darwin_prot-scriber)||

### Compilation from source code

#### Prerequisites 

`prot-scriber` is written in Rust. That makes it extremely performant and, once compiled for your operating system (OS), can be used on any machine with that particular OS environment. To compile it for your platform you first need to have Rust and cargo installed. Follow [the official instructions](https://www.rust-lang.org/tools/install).

#### Obtain the code

Download the [latest stable release of `prot-scriber` here](https://github.com/usadellab/prot-scriber/archive/refs/tags/latest-stable.zip).

Unzip it, e.g. by double clicking it or by using the command line:
```sh
unzip latest-stable.zip
```

#### Compile `prot-scriber`

Change into the directory of the downloaded `prot-scriber` code, e.g. in a Mac OS or Linux terminal `cd prot-scriber-latest-stable` after having unpacked the latest stable release (see above).

Now, compile `prot-scriber` with
```sh
cargo build --release
```

The above compilation command has generated an executable binary file in the current directory `./target/release/prot-scriber`. You can just go ahead and use `prot-scriber` now that you have compiled it successfully (see Usage below).

#### Global installation

If you are familiar with installing self compiled tools on a system wide level, this section will provide no news to you. It is convenient to make the compiled executable `prot-scriber` program available from anywhere on your system. To achieve this, you need to copy it to any place you typically have your programs installed, or add its directory to your, our all users' `$PATH` environment. In doing so, e.g. in case you are a system administrator, you make `prot-scriber` available for all users of your infrastructure. Make sure you and your group have executable access rights to the file. You can adjust these access right with `chmod ug+x ./target/release/prot-scriber`. You, and possibly other users of your system, are now ready to run `prot-scriber`.

### Install via bioconda
In case you are using [conda](https://docs.conda.io/en/latest/) to manage your pacakges, `prot-scriber` is available on [bioconda](https://anaconda.org/bioconda/prot-scriber). Download via

```
conda install -c bioconda prot-scriber
```
or create and activate a new conda environment via
```
conda create -n prot-scriber -c bioconda -c conda-forge prot-scriber
conda activate prot-scriber
```

## Usage

`prot-scriber` is a command line tool and _must_ be used in a terminal application. On Windows that will be `cmd` or PowerShell, on Mac OS X or any Linux / Unix system that will be a standard terminal shell.

### Manual

<details>
    <summary><b>Please read the manual of the latest stable version (<i>click to expand</i>).</b></summary>

```
prot-scriber version 0.1.4

PLEASE USE '--help' FOR MORE DETAILS!

prot-scriber assigns human readable descriptions (HRD) to query biological sequences or sets of them
(a.k.a gene-families).

USAGE:
    prot-scriber [OPTIONS] --output <output> --seq-sim-table <seq-sim-table>

OPTIONS:
    -a, --annotate-non-family-queries
            Use this option only in combination with --seq-families (-f), i.e. when prot-scriber is
            used to generate human readable descriptions for gene families. If in that context this
            flag is given, queries for which there are sequence similarity search (Blast) results
            but that are NOT member of a sequence family will receive an annotation (human readable
            description) in the output file, too. Default value of this setting is 'OFF' (false).

    -b, --blacklist-regexs <blacklist-regexs>
            A file with regular expressions (Rust syntax), one per line. Any match to any of these
            regular expressions causes sequence similarity search result descriptions ('stitle' in
            Blast terminology) to be discarded from the prot-scriber annotation process. If multiple
            --seq-sim-table (-s) args are provided make sure the --blacklist-regexs (-b) args appear
            in the correct order, e.g. the first -b arg will be used for the first -s arg, the
            second -b will be used for the second -s and so on. Set to 'default' to use the hard
            coded default. An example file can be downloaded here:
            https://raw.githubusercontent.com/usadellab/prot-scriber/master/misc/blacklist_stitle_regexs.txt
            - Note that this is an expert option.

    -c, --capture-replace-pairs <capture-replace-pairs>
            A file with pairs of lines. Within each pair the first line is a regular expressions
            (fancy-regex syntax) defining one or more capture groups. The second line of a pair is
            the string used to replace the match in the regular expression with. This means the
            second line contains the capture groups (fancy-regex syntax). These pairs are used to
            further filter the sequence similarity search result descriptions ('stitle' in Blast
            terminology). In contrast to the --filter-regex (-l) matches are not deleted, but
            replaced with the second line of the pair. Filtering is used to process descriptions
            ('stitle' in Blast terminology) and prepare the descriptions for the prot-scriber
            annotation process. If multiple --seq-sim-table (-s) args are provided make sure the
            --capture-replace-pairs (-c) args appear in the correct order, e.g. the first -c arg
            will be used for the first -s arg, the second -c will be used for the second -s and so
            on. Set to 'default' to use the hard coded default. An example file can be downloaded
            here:
            https://raw.githubusercontent.com/usadellab/prot-scriber/master/misc/capture_replace_pairs.txt
            - Note that this is an expert option.

    -d, --polish-capture-replace-pairs
            The last step of the process generating human readable descriptions (HRDs) for the
            queries (proteins or sequence families) is to 'polish' the selected HRDs. Polishing is
            done by iterative application of regular expressions (fancy-regex) and replace
            instructions (capture-replace-pairs). If you do not want to use the default polishing
            capture replace pairs specify a file in which pairs of lines are given. Of each pair the
            first line hold a regular expression (fancy-regex syntax) and the second the replacement
            instructions providing access to capture groups. Set to 'none' or provide an empty file,
            if you want to suppress polishing. If you want to have a template file for your custom
            polishing capture-replace-pairs please refer to
            https://raw.githubusercontent.com/usadellab/prot-scriber/master/misc/polish_capture_replace_pairs.txt
            - Note that this an expert option.

    -e, --header <header>
            Header of the --seq-sim-table (-s) arg. Separated by space (' ') the names of the
            columns in order of appearance in the respective table. Required and default columns are
            'qacc sacc stitle'. Note that this option only understands Blast terminology, i.e. even
            if you ran Diamond, please provide 'qacc' instead of 'qseqid' and 'sacc' instead of
            'sseqid'. Luckily 'stitle' is 'stitle' in Diamond, too. You can have additional columns
            that will be ignored, as long as the required columns appear in the correct order.
            Consider this example: 'qacc sacc evalue bitscore stitle'. If multiple --seq-sim-table
            (-s) args are provided make sure the --header (-e) args appear in the correct order,
            e.g. the first -e arg will be used for the first -s arg, the second -e will be used for
            the second -s and so on. Set to 'default' to use the hard coded default.

    -f, --seq-families <seq-families>
            A file in which families of biological sequences are stored, one family per line. Each
            line must have format 'fam-name TAB gene1,gene2,gene3'. Make sure no gene appears in
            more than one family.

    -g, --seq-family-gene-ids-separator <seq-family-gene-ids-separator>
            A regular expression (Rust syntax) used to split the list of gene-identifiers in the
            argument --seq-families (-f) gene families file. Default is '(\s*,\s*|\s+)'.

    -h, --help
            Print help information

    -i, --seq-family-id-genes-separator <seq-family-id-genes-separator>
            A string used as separator in the argument --seq-families (-f) gene families file. This
            string separates the gene-family-identifier (name) from the gene-identifier list that
            family comprises. Default is '<TAB>' ("\t").

    -l, --filter-regexs <filter-regexs>
            A file with regular expressions (Rust syntax), one per line. Any match to any of these
            regular expressions causes the matched sub-string to be deleted, i.e. filtered out.
            Filtering is used to process descriptions ('stitle' in Blast terminology) and prepare
            the descriptions for the prot-scriber annotation process. In case of UniProt sequence
            similarity search results (Blast result tables), this removes the Blast Hit identifier
            (`sacc`) from the description (`stitle`) and also removes the taxonomic information
            starting with e.g. 'OS=' at the end of the `stitle` strings. If multiple --seq-sim-table
            (-s) args are provided make sure the --filter-regexs (-l) args appear in the correct
            order, e.g. the first -l arg will be used for the first -s arg, the second -l will be
            used for the second -s and so on. Set to 'default' to use the hard coded default. An
            example file can be downloaded here:
            https://raw.githubusercontent.com/usadellab/prot-scriber/master/misc/filter_stitle_regexs.txt
            - Note that this is an expert option.

    -n, --n-threads <n-threads>
            The maximum number of parallel threads to use. Default is the number of logical cores.
            Required minimum is two (2). Note that at most one thread is used per input sequence
            similarity search result (Blast table) file. After parsing these annotation may use up
            to this number of threads to generate human readable descriptions.

    -o, --output <output>
            Filename in which the tabular output will be stored.

    -p, --field-separator <field-separator>
            Field-Separator of the --seq-sim-table (-s) arg. The default value is the '<TAB>'
            character. Consider this example: '-p @'. If multiple --seq-sim-table (-s) args are
            provided make sure the --field-separator (-p) args appear in the correct order, e.g. the
            first -p arg will be used for the first -s arg, the second -p will be used for the
            second -s and so on. You can provide '-p default' to use the hard coded default (TAB).

    -q, --center-inverse-word-information-content-at-quantile <center-inverse-word-information-content-at-quantile>
            The quantile (percentile) to be subtracted from calculated inverse word information
            content to center these values. Consequently, this must be a value between zero and one
            or literal 50, which is interpreted as mean instead of a quantile. Default is 5o,
            implying centering at the mean. Note that this is an expert option.

    -r, --description-split-regex <description-split-regex>
            A regular expression in Rust syntax to be used to split descriptions (`stitle` in Blast
            terminology) into words. Default is '([~_\-/|\;,':.\s]+)'. Note that this is an expert
            option.

    -s, --seq-sim-table <seq-sim-table>
            File in which to find sequence similarity search results in tabular format (SSST). Use
            e.g. Blast or Diamond to produce them. Required columns are: 'qacc sacc stitle' (Blast)
            or 'qseqid sseqid stitle' (Diamond). (See section '2. prot-scriber input preparation'
            for more details.) If the required columns, or more, appear in different order than
            shown here you must use the --header (-e) argument. If any of the input SSSTs uses a
            different field-separator than the '<TAB>' character, you must provide the
            --field-separator (-p) argument. You can provide multiple SSSTs, simply by repeating the
            -s argument, e.g. '-s queries_vs_swissprot_diamond_out.txt -s
            queries_vs_trembl_diamond_out.txt'. Providing multiple --seq-sim-table (-s) arguments
            might imply the order in which you give other arguments like --header (-e) and
            --field-separator (-p). See there for more details.

    -v, --verbose
            Print informative messages about the annotation process.

    -V, --version
            Print version information

    -w, --non-informative-words-regexs <non-informative-words-regexs>
            The path to a file in which regular expressions (regexs) are stored, one per line. These
            regexs are used to recognize non-informative words, which will only receive a minimun
            score in the prot-scriber process that generates human readable description. There is a
            default list hard-coded into prot-scriber. An example file can be downloaded here:
            https://raw.githubusercontent.com/usadellab/prot-scriber/master/misc/non_informative_words_regexs.txt
            - Note that this is an expert option.

    -x, --exclude-not-annotated-queries
            Exclude results from the output table that could not be annotated, i.e. 'unknown
            protein' or 'unknown sequence family', respectively.



MANUAL
======

1. Summary
----------
'prot-scriber' uses reference descriptions ('stitle' in Blast terminology) from sequence similarity
search results (Blast Hits) to assign short human readable descriptions (HRD) to query biological
sequences or sets of them (a.k.a gene, or sequence, families). In this, prot-scriber consumes
sequence similarity search (Blast, Diamond, or similar) results in tabular format. A customized
lexical analysis is carried out on the descriptions ('stitle' in Blast terminology) of these Blast
Hits and a resulting HRD is assigned to the query sequences or query families, respectively.

2. prot-scriber input preparation
---------------------------------
This sections explains how to run your favorite sequence similarity search tool, so that it produces
tabular results in the format prot-scriber needs them. You can run sequence similarity searches with
Blast [McGinnis, S. & Madden, T. L. BLAST: at the core of a powerful and diverse set of sequence
analysis tools. Nucleic Acids Res 32, W20–W25 (2004).] or Diamond [Buchfink, B., Xie, C. & Huson,
D. H. Fast and sensitive protein alignment using DIAMOND. Nat Meth 12, 59–60 (2015).]. Note that
there are other tools to carry out sequence similarity searches which can be used to generate the
input for prot-scriber. As long as you have a tabular text file with the three required columns
holding the query identifier, the subject ('Hit') identifier, and the subject ('Hit') description
('stitle' in Blast terminology) prot-scriber will accept this as input.
Depending on the type of your query sequences the search method and searched reference databases
vary. For amino acid queries search protein reference databases, for nucleotide query sequences
search nucleotide reference databases. If you have protein coding nucleotide query sequences you can
choose to either search protein reference databases using translated nucleotide queries with
'blastx' or 'diamond blastx' or search reference nucleotide databases with 'blastn' or 'diamond
blastn'. Note, that before carrying out any sequence similarity searches you need to format your
reference databases. This is achieved by either the 'makeblastdb' (Blast) or 'makedb' (Diamond)
commands, respectively. Please see the respective tool's (Blast or Diamond) manual for details on
how to format your reference sequence database.

2.1 A note on TAB characters
----------------------------
TAB is often used as a field separator, e.g. by default in Diamond sequence similarity search result
tables, or to separate gene-family identifiers from their respective gene-lists. Consequently,
prot-scriber has several arguments that could be a TAB, e.g. the --field-separator (-p) or the
--seq-family-id-genes-separator (-i) (please see below for more details on these arguments).
Unfortunately providing the TAB character as a command line argument can be tricky. It is even more
tricky to write it into a manual like this, because it appears as a blank whitespace and cannot
easily be distiunguished from other whitespace characters. We thus write '<TAB>' whenever we mean
the TAB character. To type it in the command line and provide it as an argument to prot-scriber you
can (i) either use $'\t' (e.g. -p $'\t') or (ii) hit Ctrl+v and subsequently hit the TAB key on your
keyboard (e.g. -p '	').

2.2 Which reference databases to search
---------------------------------------
For amino acid (protein) or protein coding nucleotide query sequences we recommend searching
UniProt's Swissprot and trEMBL. For nucleotide sequences UniRef100 and, or UniParc might be good
choices. Note that you can search _any_ database you deem to hold valuable reference sequences.
However, you might have to provide custom blacklist, filter, and capture-replace arguments for Blast
or Diamond output tables stemming from searches in these non UniProt databases (see section '3.
Technical manual' on the arguments --blacklist-regexs (-b), --filter-regexs (-l), and
--capture-replace-pairs (-c) for further details). If you want to search any NCBI reference
database, please see section 2.2.1 for more details.

2.2.1 NCBI reference databases
------------------------------
The National Center for Biotechnology Information (NCBI) has excellent reference databases to be
searched by Blast or Diamond, too. Note that NCBI and UniProt update each other's databases very
frequently. So, by searching UniProt only you should not loose information. Anyway, NCBI has e.g.
the popular non redundant ('NR') database. However, NCBI has a different description ('stitle' in
Blast terminology) format. To make sure prot-scriber parses sequence similarity search result (Blast
or Diamond) tables (SSSTs) correctly, you should use a tailored --filter-regexs (-l) argument. A
file containing such a list of regular expressions specifically tailored for parsing SSSTs produced
by searching NCBI reference databases, e.g. NR, is provided with prot-scriber. You can download it,
and edit it if neccessary, here:
https://raw.githubusercontent.com/usadellab/prot-scriber/master/misc/filter_stitle_regexs_NCBI_NR.txt

2.2.2 UniRef reference databases
------------------------------
The UniRef databases (UniProt Reference Clusters) provide clustered sets of sequences from the
UniProt Knowledgebase and selected UniParc records to obtain complete coverage of sequence space at
several resolutions (100%, 90% and 50% identity) while hiding redundant sequences. The UniRef100
database combines identical sequences and subfragments from any source organism into a single UniRef
entry (i.e. cluster). UniRef90 and UniRef50 are built by clustering UniRef100 sequences at the 90%
or 50% sequence identity levels. To make sure prot-scriber parses sequence similarity search result
(Blast or Diamond) tables (SSSTs) correctly, you should use a tailored --filter-regexs (-l)
argument. A file containing such a list of regular expressions specifically tailored for parsing
SSSTs produced by searching UniRef databases is provided with prot-scriber. You can download it, and
edit it if neccessary, here:
https://raw.githubusercontent.com/usadellab/prot-scriber/master/misc/filter_stitle_regexs_UniRef.txt

2.3 Example Blast or Diamond commands
-------------------------------------
Note that the following instructions on how to execute your sequence similarity searches with Blast
or Diamond only include the information - in terms of selected output table columns - absolutely
required by 'prot-scriber'. You are welcome, of course, to have more columns in your tabular output,
e.g. 'bitscore' or 'evalue' etc. Note that you need to search each of your reference databases with
a separate Blast or Diamond command, respectively.

2.3.1 Blast
-----------
Generate prot-scriber input with Blast as follows. The following example uses 'blastp', replace it,
if your query sequence type makes that necessary with 'blastn' or 'blastx'.

blastp -db <reference_database.fasta> -query <your_query_sequences.fasta> -num_threads
<how-many-do-you-want-to-use> -out <queries_vs_reference_db_name_blastout.txt> -outfmt "6
delim=<TAB> qacc sacc stitle"

It is important to note, that in the above 'outfmt' argument the 'delim' set to '<TAB>' means you
need to actually type in a TAB character. (We write '<TAB>' here, so you see something, not only
whitespace.) Typically you can type it by hitting Ctrl+Tab in the terminal.

2.3.2 Diamond
-------------
Generate prot-scriber input with Diamond as follows. The following example uses 'blastp', replace
it, if your query sequence type makes that necessary with 'blastn' or 'blastx'.

diamond blastp -p <how-many-threads-do-you-want-to-use> --quiet -d <reference-database.dmnd> -q
<your_query_sequences.fasta> -o <queries_vs_reference_db_name_diamondout.txt> -f 6 qseqid sseqid
stitle

Note that diamond by default uses the '<TAB>' character as a field-separator for its output tables.

2.4 Gene Family preparation and analysis
----------------------------------------
Assume you have the proteomes of eight crucifer plant species and want to cluster the respective
amino acid sequences into gene families. Note that the following example provides code to be
executed in a BASH Shell (also available on Windows). We provide a very basic procedure to perform
the clustering:

(i) "All versus all" Blast or Diamond

Assume all amino acid sequences of the eight example proteomes stored in a single file
'all_proteins.fasta'
Run:

diamond makedb --in all_proteins.fasta -d all_proteins.fasta

diamond blastp --quiet -p <how-many-threads-do-you-want-to-use?> -d all_proteins.fasta.dmnd -q
all_proteins.fasta -o all_proteins_vs_all.txt -f 6 qseqid sseqid pident

(ii) Run markov clustering

Note that 'mcl' is a command line tool implementing the original Markov Clustering algorithm [Stijn
van Dongen, A cluster algorithm for graphs. Technical Report INS-R0010, National Research Institute
for Mathematics and Computer Science in the Netherlands, Amsterdam, May 2000]. On most systems you
can install the 'mcl' binary using the respective package manager, e.g. 'sudo apt-get update && sudo
apt-get install -y mcl' (Debian / Ubuntu).

mcl all_proteins_vs_all.txt -o all_proteins_gene_clusters.txt --abc -I 2.0

(iii) Add gene family names to mcl output and filter out singleton clusters

Note that we use the GNU tools 'sed' and 'awk' to do some basic post-processing of the 'mcl' output.

sed -e 's/\t/,/g' all_proteins_gene_clusters.txt | awk -F "," 'BEGIN{i=1}{if (NF > 1){print
"Seq-Fam_" i "\t" $0; i=i+1}}' > all_proteins_gene_families.txt

Congratulations! You now have clustered your eight plant crucifer proteomes into gene families (file
'all_proteins_gene_families.txt').

(iv) Run prot-scriber

We assume that you ran either 'blastp' or 'diamond blastp' (see section 2.3 for details) to search
your selected reference databases with the 'all_proteins.fasta' queries. Here, we assume you have
searched UniProt's Swissprot and trEMBL databases.

prot-scriber -f all_proteins_gene_families.txt -s all_proteins_vs_Swissprot_blastout.txt -s
all_proteins_vs_trEMBL_blastout.txt -o all_proteins_gene_families_HRDs.txt
```

</details>

Note, that you can get the manual directly from `prot-scriber`. In the command prompt (`cmd` or PowerShell on Windows) or the Terminal (Mac OS, Linux, or Unix) use  
```sh
prot-scriber --help
```
to get it printed.

_Happy `prot-scribing`!_
    
## Speed and memory requirements
    
`prot-scriber` is **blazingly fast** and has **low memory requirements**. Consider the following two standard use cases, in which `prot-scriber` generated Human readable descriptions (HRDs) for (i) a single species and (ii) gene families.

### single species 

On a standard Laptop with 4 cores `prot-scriber` took approx. **7 seconds** and used a little under **50 MB RAM** to generate human readable descriptions for a complete plant proteome with Blast search Hits for 32,567 distinct query proteins (input: Blast result table from searches in UniProt Swissprot 66 MB, Blast result table from searches in UniProt trEMBL 144 MB)

### gene families 

On a standard Laptop with 4 cores `prot-scriber` took approx. **15 seconds** and used a little under **180 MB RAM** to generate human readable descriptions for 24,072 gene families with Blast search Hits for 71,610 distinct query proteins (input: Blast result table from searches in UniProt Swissprot 126 MB, Blast result table from searches in UniProt trEMBL 273 MB)
    
## Word-Cloud visualization of prot-scriber results
    
prot-scriber comes with a simple and small _R_ script to generate a word-cloud plot from any prot-scriber results. To use it, you must have R and the following packages installed (click [here](https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages) to learn how to install R packages):
* `RColorBrewer`                                                          
* `wordcloud`                                                             
* `wordcloud2`                                                            
* `htmlwidgets`                                                           
* `webshot`

You find the script `prot-scriber-word-cloud.R` in the `misc` directory or download it directly from [here](https://raw.githubusercontent.com/usadellab/prot-scriber/master/misc/prot-scriber-word-cloud.R).
    
In your Terminal (`cmd` or Power-Shell on Windows) you can invoke the script as follows:
    
```sh
Rscript prot-scriber-word-cloud.R input-prot-scriber-table.txt output-files-name
```

Note that the first argument is the output-table generated by prot-scriber and the second is a file-name, _without_ file extension (e.g. `.pdf`). Several output files will be created, two PDFs and one HTML.
    
_Happy word-clouding!_
    
## Development / Contribute

`prot-scriber` is open source. Please feel free and invited to contribute. 

### Preparation of releases (pre-compiled executables)

This repository is set up to use [GitHub Actions](https://github.com/features/actions) (see `.github/workflows/push.yml` for details). We use GitHub actions to trigger compilation of `prot-scriber` every time a Git Tag is pushed to this repository. So, if you, after writing new code and committing it, do the following in your local repo:

```sh
git tag -a 'version-Foo-Bar-Baz' -m "My fancy new version called Foo Bar Baz"
git push origin master --tags
```

GitHub will automatically compile `prot-scriber` and provide executable binary versions of prot-scriber for the platforms and operating systems mentioned above. The resulting binaries are then made available for download on the [releases page](https://github.com/usadellab/prot-scriber/releases).

In short, you do not need to worry about compiling your latest version to make it available for download and the different platforms and operating systems, GitHub Actions take care of this for you. 
