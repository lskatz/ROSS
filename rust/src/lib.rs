extern crate regex;

use std::io;
use std::io::prelude::*;

use regex::Regex;

/// A sequence struct that contains the ID, sequence, and quality cigar line
pub struct Seq {
  pub id:    String,
  pub seq:   String,
  pub qual: String,
}

/// A FastQ reader
pub struct Reader<R: io::Read> {
  reader:           io::BufReader<R>,
}

impl<R: io::Read> Reader<R>{
  pub fn new(reader: R) -> Self {
    Reader {
      reader: io::BufReader::new(reader),
    }
  }

  pub fn read(&mut self) -> Option<Seq> {
    let whitespace_regex = Regex::new(r"(\s+)").expect("malformed regex");

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
    seq.seq = whitespace_regex.replace_all(&seq.seq,"").into_owned();
    let read_length :usize=seq.seq.len(); 
    
    // build onto the qual line until it has the right
    // number of bytes.
    'qual: loop{
        let mut buf = String::new();
        match self.reader.read_line(&mut buf) {
            Ok(n) => {
                if n < 1 {
                    panic!("ERROR: incomplete entry (no qual line), seqid {}\nbuf {}", seq.id.trim(),buf);
                }
                else {
                  seq.qual.push_str(&buf);
                }
                seq.qual = whitespace_regex.replace_all(&seq.qual,"").into_owned();
                if seq.qual.len() >= read_length {
                  break 'qual;
                }
            }
            Err(error) => {
                panic!("ERROR while reading qual for ID {}: {}",buf.trim(),error);
            }
        }
    }

    Some(seq)
  }
}

