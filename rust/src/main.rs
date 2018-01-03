use std::io::Read;

fn main() {
    let stdin = std::io::stdin();

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
            }
            Err(error) => {
                println!("ERROR: {}",error);
            }
        }

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
                        break 'dna;
                    }
                    else {
                        seq.push_str(&buf);
                    }
                }
                Err(error) => {
                    println!("ERROR: {}",error);
                }
            }
        }
        let read_length :usize=seq.len() - 1;

        // Let's get the qual string next
        // https://stackoverflow.com/a/30679861
        let mut qual_length_counter=0;
        for byte in std::io::stdin().bytes(){
          let qual_char = byte.unwrap() as char;
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

