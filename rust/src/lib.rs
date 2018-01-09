extern crate regex;
extern crate getopts;

use std::io::Read;
use std::io::BufRead;
use std::io::BufReader;
use std::env;
use regex::Regex;
use getopts::Options;

pub struct Seq {
  pub id:    String,
  pub seq:   String,
  pub qual: String,
}

pub fn parse_args () -> Vec<String> {

    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();
    let mut opts = Options::new();
    opts.optflag("h", "help", "Print this help menu.");
    let matches = match opts.parse(&args[1..]) {
      Ok(m) => { m }
      Err(error) => { panic!(error.to_string()) }
    };
    if matches.opt_present("h") { 
      print_usage(&program,opts);
      std::process::exit(1);
    }
    
    args
}

pub fn print_usage(program: &str, opts: Options) {
  let brief = format!("Usage: {} FILE [options]", program);
  print!("{}", opts.usage(&brief));
}

// non-iterator version of reading a fastq
pub fn read_fastq_carefully_old() -> Vec<Seq> {
    let mut seq_vector :Vec<Seq>= Vec::new();
    
    let re = Regex::new(r"(\s+)").expect("malformed regex");
    let fh_in = std::io::stdin();

    'entry: loop{
        let mut seq = Seq {
          id:    String::new(),
          seq:   String::new(),
          qual:  String::new(),
        };

        // Read the ID of the entry
        match fh_in.read_line(&mut seq.id) {
            Ok(n) => {
                // if we're expecting an ID line, but
                // there are zero bytes read, then we are
                // at the end of the file. Break.
                if n < 1 {
                    break 'entry;
                }
            }
            Err(error) => {
                panic!("ERROR: {}",error);
            }
        }

        // Read the DNA line of the entry and count
        // how long it is.
        'dna: loop{
            let mut buf = String::new();
            match fh_in.read_line(&mut buf) {
                Ok(n) => {
                    if n < 1 {
                        panic!("ERROR: incomplete entry (no seq line), seqid {}\nbuf {}", seq.id.trim(),buf);
                    }
                    // if we hit the qual line, then it is a single
                    // character, +
                    else if &buf[0..1] == "+" {
                        break 'dna;
                    }
                    else {
                        seq.seq.push_str(&buf);
                    }
                }
                Err(error) => {
                    panic!("ERROR while reading seq for ID {}: {}",seq.id.trim(),error);
                }
            }
        }
        // remove all whitespace
        seq.seq = re.replace_all(&seq.seq,"").into_owned();
        let read_length :usize=seq.seq.len(); 

        // Let's get the qual string next
        // https://stackoverflow.com/a/30679861
        let mut qual_length_counter=0;
        for byte in std::io::stdin().bytes(){
          let qual_char = byte
                          .expect("Reached the end of the file too early while reading qual bytes")
                          as char;
          if qual_char == '\n' {
            if qual_length_counter == read_length {
              break;
            }
            continue;
          }
          seq.qual.push(qual_char);
          qual_length_counter+=1;
        }
          
        //println!("{}\n{}\n+\n{}",id.trim(),seq.trim(),qual);
        seq_vector.push(seq);
    }

    return seq_vector;
}

pub struct ReadFastqCarefully{
  pub reader: std::io::BufReader,
  pub whitespace_regex: regex::Regex,
}

impl ReadFastqCarefully {
  pub fn new(reader: std::io::BufReader) -> ReadFastqCarefully {
    ReadFastqCarefully{
      reader: &reader,
      whitespace_regex: Regex::new(r"(\s+)").expect("malformed regex")
    }
  }
}

impl Iterator for ReadFastqCarefully {
  type Item = (Seq);

  fn next (&mut self) -> Option<Self::Item> {
    let mut seq = Seq {
      id:    String::new(),
      seq:   String::new(),
      qual:  String::new(),
    };

    // Read the ID of the entry
    match self.reader.read_line(&mut seq.id) {
        Ok(n) => {
            // if we're expecting an ID line, but
            // there are zero bytes read, then we are
            // at the end of the file. Break.
            if n < 1 {
                return None;
            }
        }
        Err(error) => {
            panic!("ERROR: {}",error);
        }
      
    }
    // Read the DNA line of the entry and count
    // how long it is.
    'dna: loop{
        let mut buf = String::new();
        match self.reader.read_line(&mut buf) {
            Ok(n) => {
                if n < 1 {
                    panic!("ERROR: incomplete entry (no seq line), seqid {}\nbuf {}", seq.id.trim(),buf);
                }
                // if we hit the qual line, then it is a single
                // character, +
                else if &buf[0..1] == "+" {
                    break 'dna;
                }
                else {
                    seq.seq.push_str(&buf);
                }
            }
            Err(error) => {
                panic!("ERROR while reading seq for ID {}: {}",seq.id.trim(),error);
            }
        }
    }
    // remove all whitespace
    seq.seq = self.whitespace_regex.replace_all(&seq.seq,"").into_owned();
    let read_length :usize=seq.seq.len(); 

    // Let's get the qual string next
    // https://stackoverflow.com/a/30679861
    let mut qual_length_counter=0;
    for byte in self.reader.bytes(){
      let qual_char = byte
                      .expect("Reached the end of the file too early while reading qual bytes")
                      as char;
      if qual_char == '\n' {
        continue;
      }
      seq.qual.push(qual_char);
      if qual_length_counter >= read_length {
        break;
      }
      qual_length_counter+=1;
    }

    Some(seq)
  }
}

