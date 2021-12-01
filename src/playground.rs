#[cfg(test)]
mod tests {
    use std::sync::mpsc;
    use std::sync::{Arc, Mutex};
    use std::thread;

    #[test]
    fn parse_input_with_events() {
        let n_blast_tables = 2;

        let (tx, rx) = mpsc::channel();
        let v: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));

        for i in 0..n_blast_tables {
            let tx_i = tx.clone();
            let v_i = v.clone();

            thread::spawn(move || {
                let vals = vec![
                    format!("({}) hello", i),
                    env!("CARGO_MANIFEST_DIR").to_string(),
                    String::from("from"),
                    String::from("the"),
                    String::from("thread"),
                ];

                for val in vals {
                    let mut vec_i = v_i.lock().unwrap();
                    vec_i.push(val);
                }
                tx_i.send(format!("thread {}'s message", i)).unwrap();
            });
        }
        // Receiver will wait eternally, as long as tx has not moved and finished sending. So
        // manually drop it.
        drop(tx);

        for received in rx {
            println!("Got: {}", received);
        }

        println!("Data vector 'v':\n{:?}", *v.lock().unwrap());

        assert!(true);
    }
}
