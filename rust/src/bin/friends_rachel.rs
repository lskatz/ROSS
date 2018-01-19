extern crate ross;
extern crate statistical;

use std::fs::File;
use std::io::BufReader;

use ross::io::fastq;

use statistical::mean;

fn main(){
    
    //ross::parse_args();

    let my_file = File::open("/dev/stdin").expect("Could not open file");
    let my_buffer=BufReader::new(my_file);
    let mut fastq_reader=fastq::Reader::new(my_buffer);
    let mut read_length = Vec::new();
    let mut avg_qual = Vec::new();
    let mut num_bases :f32 = 0.0;
    let mut min_read_length :f32 = std::f32::MAX;
    let mut max_read_length :f32 = 0.0;
    let mut num_reads :usize = 0;
    while let Some(seq_obj) = fastq_reader.read_carefully() {
      num_reads+=1;

      // length of read
      let len = seq_obj.seq.len() as f32;
      read_length.push(len);
      num_bases += len;

      if len > max_read_length {
        max_read_length = len;
      }
      if len < min_read_length {
        min_read_length = len;
      }

      // quality
      let qual :Vec<f32> = seq_obj.qual.chars()
                                .map(|x| {x as u8 as f32-33.0})
                                .collect();
      avg_qual.push(mean(&qual));
    }
    
    println!("{}",vec![
        "read_length", "num_bases",
        "min_read_length", "max_read_length",
        "avg_quality", "num_reads", "PE?",
        "coverage",
    ].join("\t"));
    println!("{}",vec![
        mean(&read_length).to_string(),
        num_bases.to_string(),
        min_read_length.to_string(), max_read_length.to_string(),
        mean(&avg_qual).to_string(), num_reads.to_string(), 
        "0".to_string(),
    ].join("\t"));
}

