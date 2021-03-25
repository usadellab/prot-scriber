#[macro_use]
extern crate lazy_static;

use clap::{App, Arg};

mod cluster;
mod default;
mod hit;
mod models;
mod query;
mod seq_sim_clustering;
mod seq_sim_table_reader;

fn main() {
    let matches = App::new("Prot-Scribr")
        .version("0.1.0")
        .about("Assigns short human readable descriptions to query biological sequences using references. For this, it consumes sequence similarity search (Blast or Diamond) results in tabular format. See below on how to run your favorite sequence similarity search tool.")
        .arg(
            Arg::new("output")
            .required(true)
            .takes_value(true)
            .short('o')
            .long("output")
            .about("Filename in which the tabular output will be stored."),
        )
        .arg(
            Arg::new("seq-sim-table")
            .required(true)
            .short('s')
            .takes_value(true)
            .long("seq-sim-table")
            .multiple(true)
            .about("File in which to find sequence similarity search results in tabular format. Use e.g. Blast or Diamond to produce them. Required columns are: 'qacc sacc qlen qstart qend slen sstart send bitscore stitle'. If given in different order than shown here use the --header arg."),
        )
        .arg(
            Arg::new("header columns of sequence similarity search result table")
            .short('h')
            .takes_value(true)
            .long("header")
            .multiple(true)
            .about("Header of the --seq-sim-table (-s) arg. Separated by space (' ') the names of the columns as they appear in the respective table. Required and default columns are 'qacc sacc qlen qstart qend slen sstart send bitscore stitle'. You can have additional columns that will be ignored. If multiple --seq-sim-table (-s) args are provided make sure the --header (-h) args appear in the correct order, e.g. the first -h arg will be used for the first -s arg, the second -h will be used for the second -s and so on."),
        )
        .arg(
            Arg::new("debug")
            .short('d')
            .long("debug")
            .about("Turn debugging information on."),
        )
        .get_matches();

    let seq_sim_tables: Vec<_> = matches.values_of("seq-sim-table").unwrap().collect();
    println!("Input seq-sim-tables are: {:?}", seq_sim_tables);

    if let Some(o) = matches.value_of("output") {
        println!("Value for output: {}", o);
    }

    println!("Welcome!");
}
