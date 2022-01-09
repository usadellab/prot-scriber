#[cfg(test)]
mod tests {
    use std::env;
    use std::fs::File;
    use std::io::prelude::*;
    use std::io::BufReader;
    use std::thread;

    #[test]
    fn read_lines_test() {
        let path = env::var("PS_F1").unwrap();
        let file = File::open(path).unwrap();
        let mut reader = BufReader::new(file);
        let mut buf = String::new();
        loop {
            match reader.read_line(&mut buf) {
                Err(_) | Ok(0) => break,
                Ok(_) => {
                    print!(".");
                    buf.clear();
                }
            }
        }
    }

    #[test]
    fn lines_test() {
        let path = env::var("PS_F1").unwrap();
        let file = File::open(path).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let _record_line = line.unwrap();
            print!(",");
        }
    }
}
