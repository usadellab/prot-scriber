# prot-scriber

Assigns short human readable descriptions (HRD) to query biological sequences using reference candidate descriptions. In this, `prot-scriber` consumes sequence similarity search (Blast or Diamond or similar) results in tabular format. A customized lexical analysis is carried out on the descriptions of these Blast Hits and a resulting HRD is assigned to the query sequences. 

`prot-scriber` can also apply the same methodology to produce HRDs for sets of biological sequences, i.e. gene families. 

## Installation

### Pre built binaries

You can choose to download a pre-built binary ready to be executed from the table below, if you want the latest stable version. Other versions can be downloaded from the [Releases page](https://github.com/usadellab/prot-scriber/releases). Have a look at the below table to know which is the version you need for your operating system and platform. 

We strongly recommend to _rename_ the downloaded release file to a simple `prot-scriber` (or `prot-scriber.exe` on Windows).

Note that on Mac OS and Unix / Linux operating systems you need to make the downloaded and renamed `prot-scriber` file _executable_. In order to achieve this, open a terminal shell, navigate (`cd`) to the directory where you saved `prot-scriber`, and execute `chmod a+x ./prot-scriber`.


|Operating System|CPU-Architecture|Release-Name (click to download)|Comment|
|---|---|---|---|
|Windows 7 or higher|any|[windows_prot-scriber.exe](https://github.com/usadellab/prot-scriber/releases/download/latest-stable/x86_64-pc-windows-gnu_prot-scriber.exe)|to be used in a terminal (`cmd` or Power-Shell)|
|any GNU-Linux|any Intel x86, 64 bits|[x86_64-unknown-linux-gnu_prot-scriber](https://github.com/usadellab/prot-scriber/releases/download/latest-stable/x86_64-unknown-linux-gnu_prot-scriber)||
|any GNU-Linux|any aarch, 64 bits|[aarch64-unknown-linux-gnu_prot-scriber](https://github.com/usadellab/prot-scriber/releases/download/latest-stable/aarch64-unknown-linux-gnu_prot-scriber)|e.g. for Raspberry Pi|
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

## Usage

`prot-scriber` is a command line tool and _must_ be used in a terminal application. On Windows that will be `cmd` or PowerShell, on Mac OS X or any Linux / Unix system that will be a standard terminal shell.

### Manual

<details>
    <summary><b>Please read the manual of the latest stable version (<i>click to expand</i>).</b></summary>

```
prot-scriber version 0.1.0

1. Summary
----------
'prot-scriber' uses reference descriptions ('stitle' in Blast terminology) from sequence similarity
search results (Blast Hits) to assign short human readable descriptions (HRD) to query biological
sequences or sets of them (a.k.a gene, or sequence, families). In this, prot-scriber consumes
sequence similarity search (Blast, Diamond, or similar) results in tabular format. A customized
lexical analysis is carried out on the descriptions ('stitle' in Blast terminology) of these Blast
Hits and a resulting HRD is assigned to the query sequences or query families, respectively.

2. 'prot-scriber' input preparation
-----------------------------------
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

2.1 Which reference databases to search
---------------------------------------
For amino acid (protein) or protein coding nucleotide query sequences we recommend searching
UniProt's Swissprot and trEMBL. For nucleotide sequences UniRef100 and, or UniParc might be good
choices. Note that you can search _any_ database you deem to hold valuable reference sequences.

2.2 Example Blast or Diamond commands
-------------------------------------
Note that the following instructions on how to execute your sequence similarity searches with Blast
or Diamond only include the information - in terms of selected output table columns - absolutely
required by 'prot-scriber'. You are welcome, of course, to have more columns in your tabular output,
e.g. 'bitscore' or 'evalue' etc. Note that you need to search each of your reference databases with
a separate Blast or Diamond command, respectively.

2.2.1 Blast
-----------
Generate prot-scriber input with Blast as follows. The following example uses 'blastp', replace it,
if your query sequence type makes that necessary with 'blastn' or 'blastx'.

blastp -db <reference_database.fasta> -query <your_query_sequences.fasta> -num_threads <how-many-do-
you-want-to-use> -out <queries_vs_reference_db_name_blastout.txt> -outfmt "6 delim=<TAB> qacc sacc
stitle"

It is important to note, that in the above 'outfmt' argument the 'delim' set to '<TAB>' means you
need to actually type in a TAB character. (We write '<TAB>' here, so you see something, not only
whitespace.) Typically you can type it by hitting Ctrl+Tab in the terminal.

2.2.2 Diamond
-------------
Generate prot-scriber input with Diamond as follows. The following example uses 'blastp', replace
it, if your query sequence type makes that necessary with 'blastn' or 'blastx'.

diamond blastp -p <how-many-threads-do-you-want-to-use> --quiet -d <reference-database.dmnd> -q
<your_query_sequences.fasta> -o <queries_vs_reference_db_name_diamondout.txt> -f 6 qseqid sseqid
stitle

Note that diamond by default uses the '<TAB>' character as a field-separator for its output tables.

2.3 Gene Family preparation and analysis
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

sed -e 's/\t/,/g' all_proteins_gene_clusters.txt | awk -F "," 'BEGIN{i=1}{if (NF > 1){print "Seq-
Fam_" i "\t" $0; i=i+1}}' > all_proteins_gene_families.txt

Congratulations! You now have clustered your eight plant crucifer proteomes into gene families (file
'all_proteins_gene_families.txt').

(iv) Run prot-scriber

We assume that you ran either 'blastp' or 'diamond blastp' (see section 2.2 for details) to search
your selected reference databases with the 'all_proteins.fasta' queries. Here, we assume you have
searched UniProt's Swissprot and trEMBL databases.

prot-scriber -f all_proteins_gene_families.txt -s all_proteins_vs_Swissprot_blastout.txt -s
all_proteins_vs_trEMBL_blastout.txt -o all_proteins_gene_families_HRDs.txt


3. Technical manual
-------------------

USAGE:
    prot-scriber [OPTIONS] --output <output> --seq-sim-table <seq-sim-table>

OPTIONS:
    -a, --annotate-non-family-queries
            Use this option only in combination with --seq-families (-f), i.e. when prot-scriber is
            used to generate human readable descriptions for gene families. If in that context this
            flag is given, queries for which there are sequence similarity search (Blast) results
            but that are NOT member of a sequence family will receive an annotation (human readable
            description) in the output file, too. Default value of this setting is 'OFF' (false).

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

    -s, --seq-sim-table <seq-sim-table>
            File in which to find sequence similarity search results in tabular format (SSST). Use
            e.g. Blast or Diamond to produce them. Required columns are: 'qacc sacc stitle' (Blast)
            or 'qseqid sseqid stitle' (Diamond). (See section '2. prot-scriber input preparation'
            for more details.) If the required columns, or more, appear in different order than
            shown here you must use the --header (-e) argument. If any of the input SSSTs uses a
            different field-separator than the '<TAB>' character, you must provide the --field-
            separator (-p) argument. You can provide multiple SSSTs, simply by repeating the -s
            argument, e.g. '-s queries_vs_swissprot_diamond_out.txt -s
            queries_vs_trembl_diamond_out.txt'. Providing multiple --seq-sim-table (-s) arguments
            might imply the order in which you give other arguments like --header (-e) and --field-
            separator (-p). See there for more details.

    -v, --verbose
            Print informative messages about the annotation process.

    -V, --version
            Print version information
```

</details>

Note, that you can get the manual directly from `prot-scriber`. In the command prompt (`cmd` or PowerShell on Windows) or the Terminal (Mac OS, Linux, or Unix) use  
```sh
prot-scriber --help
```
to get it printed.

_Happy `prot-scribing`!_

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
