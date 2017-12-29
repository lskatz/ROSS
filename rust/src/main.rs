extern crate regex;
use std::io::Read;
use regex::Regex;

fn main() {
    let stdin = std::io::stdin(); //io::stdin();
    let whitespace_regex = Regex::new(r"\s").unwrap();
    //let qual_offset=33;

    'entry: loop{
        let mut buf :String="".to_string();
        let mut id  :String="".to_string();

        // Read the ID of the entry
        match stdin.read_line(&mut buf) {
            Ok(n) => {
                if n < 1 {
                    break;
                }
                //println!("{} bytes read",n);
                id = buf.clone();
                println!("{}", id.trim());
            }
            Err(error) => {
                println!("ERROR: {}",error);
            }
        }

        // Read the DNA line of the entry and count
        // how long it is.
        let mut dna :String="".to_string();
        'dna: loop{
            buf="".to_string();
            match stdin.read_line(&mut buf) {
                Ok(n) => {
                    //println!("==> {}",buf.chars().nth(0).unwrap());
                    if n < 1 {
                        println!("ERROR: incomplete entry, seqid {}\n{}", id.trim(),buf);
                        break 'entry;
                    }
                    // if we hit the qual line, then it is a single
                    // character, +.
                    else if buf.chars().nth(0).unwrap() == '+' {
                        println!("{}",dna.trim());
                        println!("+");
                        break 'dna;
                    }
                    else {
                        dna = format!("{}{}",dna.trim(),buf.trim());
                    }
                }
                Err(error) => {
                    println!("ERROR: {}",error);
                }
            }
        }
        let read_length :usize=dna.len();
        //println!("{}", read_length);

        // Let's get the qual string next
        let mut qual :String = "".to_string();
        // https://stackoverflow.com/a/30679861
        for _i in 0..read_length {
            let byte :Option<i32> = std::io::stdin()
                .bytes()
                .next()
                .and_then(|result| result.ok())
                .map(|byte| byte as i32);
            let qual_char = byte.unwrap().to_string();
            let qual_int  = qual_char.parse::<i16>().unwrap();// - qual_offset;
            let qual_u8 :u8 = qual_int as u8;
            let mut qual_encoding :char = qual_u8 as char;

            // Whitespace unacceptable here, so go to the next byte
            // if it is found.
            if whitespace_regex.is_match(&qual_encoding.to_string()) { 
                let byte :Option<i32> = std::io::stdin()
                    .bytes()
                    .next()
                    .and_then(|result| result.ok())
                    .map(|byte| byte as i32);
                let qual_char = byte.unwrap().to_string();
                let qual_int  = qual_char.parse::<i16>().unwrap();// - qual_offset;
                let qual_u8 :u8 = qual_int as u8;
                qual_encoding = qual_u8 as char;
            }

            qual = format!("{}{}",qual,qual_encoding);
        }
        // the next char in the file should be \n, so it
        // needs to be burned.
        let byte :Option<i32> = std::io::stdin()
            .bytes()
            .next()
            .and_then(|result| result.ok())
            .map(|byte| byte as i32);
        let qual_char = byte.unwrap().to_string();
        let qual_int  = qual_char.parse::<i16>().unwrap();// - qual_offset;
        let qual_u8 :u8 = qual_int as u8;
        let qual_encoding = qual_u8 as char;
        if &qual_encoding.to_string() != "\n" {
            println!("ERROR: newline not found where expected at the end of qual line");
        }

        println!("{}",qual);

    }

}

