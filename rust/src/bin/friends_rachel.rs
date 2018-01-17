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
    while let Some(seq_obj) = fastq_reader.read_carefully() {
      // length of read
      let len = seq_obj.seq.len() as f32;
      read_length.push(len);

      // quality
      let qual :Vec<f32> = seq_obj.qual.chars()
                                .map(|x| {x as u8 as f32-33.0})
                                .collect();
      let avgQual = mean(&qual);
      println!("{:?}",avgQual);
      break;
    }
    
    println!("Mean read length {}", mean(&read_length));
}

