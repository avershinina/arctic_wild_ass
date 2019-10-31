# README
This repo complements the paper by Vershinina et al 201X "The case of an arctic wild ass highlights the utility of ancient DNA for validating problematic identifications in museum collections", considered for review at Molecular Ecology Resources.

## Reference-guided aseembly of ancient mitochondrial genome
This script is a pipeline that utilizes a set of bioinformatic tools to process anicent DNA sequencing data with the aim of assembling a complete mitochondrial genome using a reference. This script works only for paired end Illumina data.

## To run the pipeline
Inside of the script, a user would need to include the path to reference genome, path to software, flags and settings. All these are hard-coded inside of the script. To run the pipeline, edit all corresponding variables manually using any script-friendly text editor, such as sublime-text or gedit. Make the script executable: chmod +x mito_assembly.sh

To run the pipeline, input three variables: path to forward and reverse read files and arbitrary sample name (used to name output files). 
Run it as follows:
```bash
./mito_assembly.sh path/to/sample.R1.fastq.gz path/to/sample.R2.fastq.gz samplename
```

## Steps of the pipeline

1) Trim adapters and merge reads. 
**Note.** Since we are processing ancient DNA, we assume that endogenious sequences are mainly short (30-70bp). This, however, may not be true for well preserved DNA. User should check fragment length distribution before deciding if most reads in a read pair are overlapping and whether they should be merged or not. Not all of the sequences are overlapping enough to merge them. Thus, the current version of the script utilizes both merged (short) and unmerged (long) reads. User may be interested in using only merged reads (cases of extremely degraded samples). If that is the case, the script should be manually modified according to user's needs. 

2) Remove low complexity reads (such as stretches of AAAAs, ATATATA, etc).

3) Concatenate files: collect both merged and unmerged reads together.

4) Convert fastq into fasta, applying a phred quality filter.

5) Remove duplicated reads (with the aim to reduce % of PCR duplicates).

6) Run MIA: mapping-iterative-assembler using a reference mitochondrial genome.

## Dependencies
The following programs should be installed before running the pipeline:
* SEQPREP https://github.com/jeizenga/SeqPrep2
* MIA https://github.com/mpieva/mapping-iterative-assembler/tree/5a7fb5afad735da7b8297381648049985c599874
* PRINSEQ http://prinseq.sourceforge.net/
* BBMAP https://github.com/BioInfoTools/BBMap
* FASTX_TOOLKIT http://hannonlab.cshl.edu/fastx_toolkit/index.html
