#Random Operations on Sequences Suite - ROSS

Perform random operations on reads files.  Right now ROSS only supports fastq and fastq.gz format.

## General usage

All scripts accept the parameters

* `--help`
* `--numcpus` Not all scripts will take advantage of numcpus.
* `--tempdir`

##Converting scripts from CGP and other sources to ROSS

|script               |Description|    |
|---------------------|-----------|----|
|`friends_monica.pl`  | Trims and cleans a fastq file| ![Monica](/images/monica.jpg) |
|`friends_carol.pl`   | Convert any fastq file to a standard four-line-per-entry format| ![Carol](/images/carol.jpg) | 
|`friends_rachel.pl`  | Prints basic read metrics| ![Rachel](/images/rachel.jpg) |
|`friends_ung.pl`     | Determines paired-endedness| ![UNG](/images/UNG.png) |
|`friends_ross.pl`    | Makes sure a fastq file is in a standard format and is unbroken |
|`friends_phoebe.pl`  | Randomizes reads|
|`friends_chandler.pl`| Pure perl kmer counting. No outside dependencies.|
|`friends_marcel.pl`  | Rescores reads based on kmer abundance |
|`friends_ursula.pl`  | Removes duplicate reads and/or downsamples reads|
|`friends_joey.pl`    | Shuffles or deshuffles paired end reads|
|`friends_barry.pl`   | Joins overlapping paired ends together |

