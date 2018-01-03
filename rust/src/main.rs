extern crate regex;
use std::io::Read;
use regex::Regex;

fn main() {
    let stdin = std::io::stdin();
    let whitespace_regex = Regex::new(r"\s").unwrap();

    'entry: loop{
        let mut id   = String::new();
        let mut seq  = String::new();
        let mut qual = String::new();

        let mut buf  = String::new();

        // Read the ID of the entry
        match stdin.read_line(&mut buf) {
            Ok(n) => {
                if n < 1 {
                    break;
                }
                id = String::from(buf);
                id = String::from(id.trim());
                //println!("{}", id.trim());
            }
            Err(error) => {
                println!("ERROR: {}",error);
            }
        }

        // reset the buffer
        let mut buf = String::new();

        // Read the DNA line of the entry and count
        // how long it is.
        'dna: loop{
            buf="".to_string();
            match stdin.read_line(&mut buf) {
                Ok(n) => {
                    if n < 1 {
                        println!("ERROR: incomplete entry, seqid {}\n{}", id.trim(),buf);
                        break 'entry;
                    }
                    // if we hit the qual line, then it is a single
                    // character, +
                    else if buf.chars().nth(0).unwrap() == '+' {
                        seq = String::from(seq.trim());
                        break 'dna;
                    }
                    else {
                        //seq = seq + &buf;
                        seq.push_str(&buf);
                    }
                }
                Err(error) => {
                    println!("ERROR: {}",error);
                }
            }
        }
        let read_length :usize=seq.len();
        //println!("{}", read_length);
        //println!("=>{}<=", seq);

        // Let's get the qual string next
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

            qual.push(qual_encoding);
            //qual = format!("{}{}",qual,qual_encoding);
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

        println!("{}\n{}\n+\n{}",id,seq,qual);
        break;

    }
}

