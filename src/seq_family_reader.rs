use super::default::SPLIT_GENE_FAMILY_GENES_REGEX;
use crate::annotation_process::AnnotationProcess;
use crate::seq_family::SeqFamily;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

pub fn parse_seq_families_file(path: &str, annotation_process: &mut AnnotationProcess) {
    // Open stream to the gene families input file
    let file_path = path.to_string();
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    // read file line by line
    for (i, line) in reader.lines().enumerate() {
        let family_line = line.unwrap();

        // parse line. panic if malformatted, add to the annotation_process if OK
        match parse_seq_family(family_line) {
            Ok((seq_fam_name, seq_fam_instance)) => {
                annotation_process.insert_seq_family(seq_fam_name, seq_fam_instance)
            }
            Err(e) => panic!("{:?} in file {:?} line <{:?}>. The expected format is \"<family-name>TABgene1,gene2,gene3,...\"", e, file_path, i),
        }
    }
}

fn parse_seq_family(family: String) -> Result<(String, SeqFamily), Box<dyn Error>> {
    // split the line by "\t". There should be more than 1 element (>=2), panic if not
    let family_cols: Vec<&str> = family.trim().split("\t").collect();
    if family_cols.len() < 2 {
        return Err("Malformatted gene family".into());
    }

    let seq_fam_name: String = String::from(family_cols[0]); // the family name
    let mut seq_fam_instance = SeqFamily::new(); // the genes contained in that family

    // split the gene column using the default separator.
    let gene_cols: Vec<String> = SPLIT_GENE_FAMILY_GENES_REGEX
        .split(family_cols[1])
        .map(String::from)
        .collect();
    // set family genes
    seq_fam_instance.query_ids = gene_cols;
    // return OK
    Ok((seq_fam_name, seq_fam_instance))
}

#[cfg(test)]
mod tests {
    use super::*;
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
        let line_2: String = "OG0023617	VFABAHed036352 VFABAHed036353".to_string();
        let line_3: String = "OG0023617	VFABAHed036352   VFABAHed036353, VFABAHed036354".to_string();

        match parse_seq_family(line_1) {
            Ok((seq_fam_name, seq_fam_instance)) => {
                assert_eq!(seq_fam_name, "OG0023617");
                assert_eq!(seq_fam_instance.query_ids.len(), 2);
            }
            Err(e) => println!("{}", e),
        }

        match parse_seq_family(line_2) {
            Ok((seq_fam_name, seq_fam_instance)) => {
                assert_eq!(seq_fam_name, "OG0023617");
                assert_eq!(seq_fam_instance.query_ids.len(), 2);
            }
            Err(e) => println!("{}", e),
        }

        match parse_seq_family(line_3) {
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
        assert!(parse_seq_family(line).is_err())
    }
}
