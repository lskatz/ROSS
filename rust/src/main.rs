extern crate regex;
extern crate getopts;
use std::io::Read;
use regex::Regex;

use std::env;
use getopts::Options;

fn main(){
    
    parse_args();

    reformat_stdin();

}

fn parse_args () -> Vec<String> {

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

fn print_usage(program: &str, opts: Options) {
  let brief = format!("Usage: {} FILE [options]", program);
  print!("{}", opts.usage(&brief));
}

fn reformat_stdin() {
    
    let re = Regex::new(r"(\s+)").expect("malformed regex");
    let stdin = std::io::stdin();

    'entry: loop{
        let mut id  :String = String::new();
        let mut seq :String = String::new();
        let mut qual:String = String::new();

        // Read the ID of the entry
        match stdin.read_line(&mut id) {
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
            match stdin.read_line(&mut buf) {
                Ok(n) => {
                    if n < 1 {
                        panic!("ERROR: incomplete entry (no seq line), seqid {}\nbuf {}", id.trim(),buf);
                    }
                    // if we hit the qual line, then it is a single
                    // character, +
                    else if &buf[0..1] == "+" {
                        break 'dna;
                    }
                    else {
                        seq.push_str(&buf);
                    }
                }
                Err(error) => {
                    panic!("ERROR while reading seq for ID {}: {}",id.trim(),error);
                }
            }
        }
        // remove all whitespace
        seq = re.replace_all(&seq,"").into_owned();
        let read_length :usize=seq.len(); 

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
          qual.push(qual_char);
          qual_length_counter+=1;
        }
          
        println!("{}\n{}\n+\n{}",id.trim(),seq.trim(),qual);
    }
}

