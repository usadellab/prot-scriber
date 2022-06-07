use crate::annotation_process::AnnotationProcess;
use crate::seq_family::SeqFamily;
use regex::Regex;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Parses line by line of the argument file `path` in which sets of biological sequence
/// identifiers (a.k.a. gene families) are stored; one family per line. Each parsed family is
/// stored in the argument `annotation_process`.
///
/// # Arguments
///
/// * `path` - The valid path to the file holding the to be parsed gene families.
/// * `annotation_process` - The AnnotationProcess to be provided with the parsed gene families.
pub fn parse_seq_families_file(path: &str, annotation_process: &mut AnnotationProcess) {
    // Open stream to the gene families input file
    let file_path = path.to_string();
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    // read file line by line
    for (i, line) in reader.lines().enumerate() {
        let family_line = line.unwrap();

        // parse line. panic if malformatted, add to the annotation_process if OK
        match parse_seq_family(family_line, &annotation_process.seq_family_id_genes_separator, &annotation_process.seq_family_gene_ids_separator) {
            Ok((seq_fam_name, seq_fam_instance)) => {
                annotation_process.insert_seq_family(seq_fam_name, seq_fam_instance)
            }
            Err(e) => panic!("\n\n{:?} in file {:?} line <{:?}>. The expected format is \"<family-name>TABgene1,gene2,gene3,...\"\n\n", e, file_path, i),
        }
    }
}

/// Parses a single line read from a respective "gene family input file" (see
/// `parse_seq_families_file`). The argument `family: String` is split into a family identifier
/// (name) and the set of sequence identifiers the family comprises, is made of. Returns a
/// `Result<(String, SeqFamily), Box<dyn Error>>` holding either the parsed gene family or an
/// error.
///
/// # Arguments
///
/// * `family` - The single line (`String`) holding the gene family information
/// * `fam_id_from_gene_id_list_separator` - The character that separates a gene-family's
/// identifier from the list of gene-identifiers the family comprises.
/// * `gene_ids_separator_regex` - The string representation of a regular expression to be used to
/// split the list of gene-identifiers.
fn parse_seq_family(
    family: String,
    fam_id_from_gene_id_list_separator: &String,
    gene_ids_separator_regex: &String,
) -> Result<(String, SeqFamily), Box<dyn Error>> {
    // Split the line by argument `fam_id_from_gene_id_list_separator`. There should be more than 1
    // element (>=2), panic if not:
    let family_cols: Vec<&str> = family
        .trim()
        .split(fam_id_from_gene_id_list_separator)
        .map(|x| x.trim())
        .filter(|x| !x.is_empty())
        .collect();
    if family_cols.len() < 2 {
        return Err("Malformatted gene family".into());
    }

    let seq_fam_name: String = String::from(family_cols[0]); // the family name
    let mut seq_fam_instance = SeqFamily::new(); // the genes contained in that family

    // split the gene column using the default separator.
    // In case genes were separated by a <TAB> character
    // we need to re-join the remaining elements.
    let split_gene_id_list_regex = Regex::new(gene_ids_separator_regex).unwrap();
    let gene_cols: Vec<String> = split_gene_id_list_regex
        .split(&family_cols[1..].join("\t"))
        .map(|x| x.trim().to_string())
        .filter(|x| !x.is_empty())
        .collect();

    // set family genes
    seq_fam_instance.query_ids = gene_cols;
    // return OK
    Ok((seq_fam_name, seq_fam_instance))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::default::{SPLIT_GENE_FAMILY_GENES_REGEX, SPLIT_GENE_FAMILY_ID_FROM_GENE_SET};
    use std::path::Path;

    #[test]
    fn parses_seq_families_file() {
        let mut ap = AnnotationProcess::new();
        let p = Path::new("misc")
            .join("test_gene_families.txt")
            .to_str()
            .unwrap()
            .to_string();
        parse_seq_families_file(&p, &mut ap);
        assert_eq!(ap.seq_families.len(), 6)
    }

    #[test]
    fn parses_correct_lines_ok() {
        let line_1: String = "OG0023617	VFABAHed036352,VFABAHed036353".to_string();
        let line_2: String =
            "OG0023617	VFABAHed036490 ,	 VFABAHed036491  VFABAHed036353  ".to_string();
        let line_3: String = "OG0023617	VFABAHed036352	VFABAHed036353, VFABAHed036354".to_string();

        match parse_seq_family(
            line_1,
            &(*SPLIT_GENE_FAMILY_ID_FROM_GENE_SET).to_string(),
            &(*SPLIT_GENE_FAMILY_GENES_REGEX).to_string(),
        ) {
            Ok((seq_fam_name, seq_fam_instance)) => {
                assert_eq!(seq_fam_name, "OG0023617");
                assert_eq!(seq_fam_instance.query_ids.len(), 2);
            }
            Err(e) => println!("{}", e),
        }

        match parse_seq_family(
            line_2,
            &(*SPLIT_GENE_FAMILY_ID_FROM_GENE_SET).to_string(),
            &(*SPLIT_GENE_FAMILY_GENES_REGEX).to_string(),
        ) {
            Ok((seq_fam_name, seq_fam_instance)) => {
                assert_eq!(seq_fam_name, "OG0023617");
                assert_eq!(seq_fam_instance.query_ids.len(), 3);
            }
            Err(e) => println!("{}", e),
        }

        match parse_seq_family(
            line_3,
            &(*SPLIT_GENE_FAMILY_ID_FROM_GENE_SET).to_string(),
            &(*SPLIT_GENE_FAMILY_GENES_REGEX).to_string(),
        ) {
            Ok((seq_fam_name, seq_fam_instance)) => {
                assert_eq!(seq_fam_name, "OG0023617");
                assert_eq!(seq_fam_instance.query_ids.len(), 3);
            }
            Err(e) => println!("{}", e),
        }
    }

    #[test]
    fn parse_faulty_line_malformatted() {
        let line = "OG0023619|VFABAHed036490,VFABAHed036491".to_string();
        assert!(parse_seq_family(
            line,
            &(*SPLIT_GENE_FAMILY_ID_FROM_GENE_SET).to_string(),
            &(*SPLIT_GENE_FAMILY_GENES_REGEX).to_string(),
        )
        .is_err())
    }
}
