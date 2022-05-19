# Random Operations on Sequences Suite - ROSS

Perform random operations on reads files.  Right now ROSS only supports fastq and fastq.gz format.

## Deprecation

🛑 ROSS has been deprecated for some time now due to the newer repo ROSS.
These executables are about 10-20x faster than ROSS and install very easily using `cargo`.
https://github.com/lskatz/fasten

## General usage

All scripts accept the parameters

* `--help`
* `--numcpus` Not all scripts will take advantage of numcpus.
* `--tempdir`

## Ross script descriptions

|script               |Description|    |
|---------------------|-----------|----|
|`friends_monica.pl`  | Trims and cleans a fastq file| ![Monica](/images/monica.jpg) |
|`friends_carol.pl`   | Convert any fastq file to a standard four-line-per-entry format| ![Carol](/images/carol.jpg) | 
|`friends_rachel.pl`  | Prints basic read metrics| ![Rachel](/images/rachel.jpg) |
|`friends_ung.pl`     | Determines paired-endedness| ![UNG](/images/UNG.png) |
|`friends_ross.pl`    | Makes sure a fastq file is in a standard format and is unbroken | ![Ross](/images/ross.png) | 
|`friends_phoebe.pl`  | Randomizes reads| ![Phoebe](/images/phoebe.png) |
|`friends_chandler.pl`| Pure perl kmer counting. No outside dependencies.| ![Chandler](/images/chander.png) |
|`friends_marcel.pl`  | Rescores reads based on kmer abundance | ![Marcel](/images/marcel.png) | 
|`friends_ursula.pl`  | Removes duplicate reads and/or downsamples reads| ![Ursula](/images/ursula.png) | 
|`friends_joey.pl`    | Shuffles or deshuffles paired end reads| ![Joey](/images/joey.png) |
|`friends_barry.pl`   | Joins overlapping paired ends together | ![Barry](/images/barry.png) |
|`friends_charles.pl` | Introduces Ns into a read set to help simulate increasing ambiguities | ![Charles](/images/charles.png) |
|`friends_roger.pl`   | Finds invalid characters in a fastq file and fixes them | ![Roger](/images/roger.png) |

## see also

ROSS.rs, a Rust implementation of ROSS https://github.com/lskatz/ROSS.rs
