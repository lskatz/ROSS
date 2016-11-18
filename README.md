#Random Operations on Sequences Suite - ROSS

Perform random operations on reads files.  Right now ROSS only supports fastq and fastq.gz format.

##Converting scripts from CGP and other sources to ROSS

|Initially imported|Formally called|Description|
|--------------------|-------|-------|
|`friends_monica.pl`  | `run_assembly_trimClean.pl`                  | Trims and cleans a fastq file|
|`friends_carol.pl`   | `run_assembly_convertMultiFastqToStandard.pl`| Convert any fastq file to a standard four-line-per-entry format|
|`friends_rachel.pl`  | `run_assembly_readMetrics.pl`                | Prints basic read metrics|
|`friends_ung.pl`     | `run_assembly_isFastqPE.pl`                  | Determines paired-endedness|
|`friends_ross.pl`| `ross_validateFastq.pl`                      | Makes sure a fastq file is in a standard format and is unbroken |
|`friends_phoebe.pl`  | N/A                                          | Randomizes reads|
|`friends_chandler.pl`| N/A                                          | Pure perl kmer counting. No outside dependencies.|
|`friends_marcel.pl`  | N/A                                          | Rescores reads based on kmer abundance |
|`friends_ursula.pl`  | `run_assembly_removeDuplicateReads.pl`       | Removes duplicate reads and/or downsamples reads|
|`friends_joey.pl`    | `run_assembly_shuffleReads.pl`               | Shuffles or deshuffles paired end reads|

