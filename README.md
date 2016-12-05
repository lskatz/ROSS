#Random Operations on Sequences Suite - ROSS

Perform random operations on reads files.  Right now ROSS only supports fastq and fastq.gz format.

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

