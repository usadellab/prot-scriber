#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::env;
    use std::fs::File;
    use std::io::prelude::*;
    use std::io::BufReader;
    use std::sync::mpsc::{channel, Receiver, Sender};
    use std::{thread, time};

    #[derive(Debug, Clone, Default)]
    struct Foo {
        pub id: String,
        pub dict: HashMap<String, Vec<String>>,
    }

    impl Foo {
        fn new(id: &String, dict_cells: Vec<String>) -> Foo {
            let mut foo = Foo {
                id: id.clone(),
                dict: HashMap::new(),
            };
            foo.dict.insert(dict_cells[1].clone(), dict_cells);
            foo
        }
    }

    fn row_to_cells(row: &str) -> Vec<String> {
        row.split("\t")
            .into_iter()
            .map(|x| (*x).to_string())
            .collect()
    }

    fn parse_file(path: &String, tx: Sender<Vec<String>>) {
        let file = File::open(path).unwrap();
        let mut reader = BufReader::new(file);
        let mut buf = String::new();
        loop {
            match reader.read_line(&mut buf) {
                Err(_) | Ok(0) => break,
                Ok(_) => {
                    let foo_cells = row_to_cells(&buf);
                    tx.send(foo_cells);
                    buf.clear();
                }
            }
        }
    }

    #[test]
    fn test_channel() {
        let paths = env::var("TEST_FILES")
            .unwrap()
            .split(",")
            .into_iter()
            .map(|x| (*x).to_string())
            .collect::<Vec<String>>();
        let mut memory_db: HashMap<String, Foo> = HashMap::new();
        let (tx, rx): (Sender<Vec<String>>, Receiver<Vec<String>>) = channel();

        for p in paths {
            let tx_i = tx.clone();

            thread::spawn(move || {
                parse_file(&p, tx_i);
                println!("Finished parsing file '{}'.", p);
            });
        }
        drop(tx);

        for foo_cells in rx {
            let foo_id = foo_cells[0].clone();
            if memory_db.contains_key(&foo_id) {
                let foo = memory_db.get_mut(&foo_id).unwrap();
                foo.dict.insert(foo_cells[1].clone(), foo_cells);
            } else {
                let foo = Foo::new(&foo_id, foo_cells);
                memory_db.insert(foo_id, foo);
                thread::sleep(time::Duration::from_millis(100));
            }
            print!(".");
        }
        println!("memory_db now contains {} entries", memory_db.len());
    }
}
