extern crate ross;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;

fn main(){
    
    //ross::parse_args();

    //let fastq_iter = ross::ReadFastqCarefully::new("/dev/stdin");
    let my_file = File::open("/dev/stdin").expect("Could not open file");
    let mut my_buffer=BufReader::new(my_file);
    let mut fastq_reader=ross::Reader::new(my_buffer);
    while let Some(seq_obj) = fastq_reader.read() {
      println!("{}\n{}\n+\n{}",seq_obj.id.trim(),seq_obj.seq,seq_obj.qual);
    }
    println!("Made it");
}

