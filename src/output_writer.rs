use std::collections::HashMap;
use std::fs::write;

/// Parse human_readable_descriptions into string and save to a file.
///
/// # Arguments
///
/// * `file_path: String` - The file path for saving output.
/// * `human_readable_descriptions: HashMap<String, String>` - The generated human readable descriptions.
pub fn write_output_table(
    file_path: String,
    human_readable_descriptions: HashMap<String, String>,
) -> std::io::Result<()> {
    if human_readable_descriptions.keys().len() > 0 {
        let mut output = String::from("Annotee-Identifier\tHuman-Readable-Description");
        // Sorted by annotee, not in HashMap order. RandomState seeds every map afresh, so
        // iterating the map directly emitted the same annotations in a different line order on
        // every run: two output files with identical content had different checksums, and diffing
        // one run against another showed the whole table as changed.
        let mut annotees: Vec<&String> = human_readable_descriptions.keys().collect();
        annotees.sort();
        // stream write line after line
        for annotee_name in annotees {
            let annotation = &human_readable_descriptions[annotee_name];
            output.push_str(&(format!("\n{}\t{}", annotee_name, annotation)));
        }
        // add trailing newline for the last annotation
        output.push_str("\n");
        return write(file_path, output);
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::output_writer::write_output_table;
    use std::collections::HashMap;
    #[test]
    fn writer_test() {
        let mut human_readable_descriptions: HashMap<String, String> = HashMap::new();
        human_readable_descriptions.insert(
            "Seq-Family-1".to_string(),
            "alien devouring protein".to_string(),
        );
        human_readable_descriptions.insert(
            "Protein-123".to_string(),
            "human devouring protein".to_string(),
        );
        assert_eq!(
            write_output_table(
                "./target/result.txt".to_string(),
                human_readable_descriptions
            )
            .is_ok(),
            true
        );
    }
}
