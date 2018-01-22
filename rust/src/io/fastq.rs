use std::io;
use std::io::prelude::*;
use io::seq::Seq;
use io::seq::Cleanable;

use regex::Regex;

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

  /// An alias for self.read_quickly().
  pub fn read(&mut self) -> Option<Seq> {
    return self.read_quickly();
  }

  /// Read a fastq entry but assume that there are only
  /// four lines per entry (id, seq, plus, qual).
  pub fn read_quickly(&mut self) -> Option<Seq> {
    let mut seq = Seq {
      id:    String::new(),
      seq:   String::new(),
      qual:  String::new(),
    };

    // Read the ID of the entry
    self.reader.read_line(&mut seq.id).expect("ERROR: could not read ID line");

    self.reader.read_line(&mut seq.seq).expect("ERROR: could not read sequence line");
    
    // burn the plus sign
    let mut _plus = String::new();
    self.reader.read_line(&mut _plus).expect("ERROR: plus sign line not found");

    self.reader.read_line(&mut seq.qual).expect("ERROR: could not read qual line");

    Some(seq)
  }

  /// Read a fastq entry in the most correct way possible,
  /// allowing for whitespace in seq and qual lines.
  pub fn read_carefully(&mut self) -> Option<Seq> {
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

