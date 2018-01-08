extern crate ross;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;

fn main(){
    
    ross::parse_args();

    //let stdin = std::io::stdin();
    let stdin = File::open("/dev/stdin").expect("Could not open /dev/stdin");
    let reader = BufReader::new(&stdin);
    let fastq_iter = ross::ReadFastqCarefully::new(&reader);

    for seq in fastq_iter {
      //println!("{}\n{}\n{}",seq.id.trim(),seq.seq,seq.qual.trim());
      println!("{}\n",seq.id.trim());
    }

    /*
    let reads :Vec<ross::Seq> = ross::read_fastq_carefully();

    println!(">{}<",reads[0].id);
    */
}

