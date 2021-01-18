use super::default::FILTER_REGEXS;
use super::models::Hit;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;

pub fn parse_table(path: &str) -> std::io::Result<()> {
    println!("path is {}", path);
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line_content = line.unwrap();
        println!("record is {:?}", line_content);
    }

    Ok(())
}

//pub fn parse_hit(line: &str) -> Hit {
//    let record: Vec<&str> = line.split('\t').collect();
//}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn use_filter_regexs() {
        println!("FILTER_REGEXS is: {:?}", *FILTER_REGEXS);
        assert_eq!(true, true);
    }

    #[test]
    fn parses_table() {
        let pt_res = parse_table(
            Path::new("misc")
                .join("Soltu.DM.09G022410.3_vs_UniProt_Sprot_blastp.txt")
                .to_str()
                .unwrap(),
        )
        .unwrap();
        assert_eq!(true, true);
    }
}
