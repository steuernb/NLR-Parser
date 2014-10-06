# NLR-Parser README

NLR-Parser is a tool to rapidly annotate the NLR complement from sequenced plant genomes.

The NLR-Parser refines the output of MAST and reliably annotates disease resistance genes encoding for nucleotide-binding leucine-rich repeat (NLR) proteins.



## Prerequisites
### MEME suite
The MEME suite is available at [http://meme.nbcr.net/meme/](http://meme.nbcr.net/meme/)

Don't worry about setting up the Apache webserver. You just need MAST, so the quick install is sufficient. 

### JRE 1.6
Make sure you have the Java Runtime Environments 1.6 or higher. Download from [http://java.com](http://java.com)

### NLR motif definitions
Download the meme.xml that contains the definitions from [here](). 
The motifs were published by [Jupe et al. (2012)](http://www.biomedcentral.com/1471-2164/13/75). The downloaded meme.xml is an input argument for MAST.

### 6Frame translator
If you intend to screen nucleotide sequences for NLRs, it might make sense to translate your sequence in all 6 reading frames. To ensure the full functionality of the NLR-Parser, please make sure the 6 aa-sequences only differ by a suffix and end with:

* _frame+0
* _frame+1
* _frame+2
* _frame-0
* _frame-1
* _frame-2 

For this you can use the TranslateSequence.jar, which is part of this software.

## Installation

Just download this file. Run it from the command line.

`java -jar NLR-Parser.jar -i <mast.xml> -o <output.mast.txt> [-s <splitpattern>] [-p <pvalue>] [-b <blastfile>] [-gh] [-a <sequence>]` 

## Input parameters
 
parameter | argument | description
---       |   ---    | ---
**-i**    | *STR*    | The location of the xml output of MAST
**-o**    | *STR*    | Location and name of the outputfile that will be generated by the NLR-Parser. Note that an existing file will be overwritten
**-s**    | *STR*    | The splitpattern to combine 6-frame-translated nucleotide sequences to one output. **default: "_frame"**
**-p**    | *float*  | P-value threshold. Motifs with a p-value above will be ignored by the NLR-Parser. **default: 1E-5**
**-b**    | *STR*    | Location of an optional blast file to detemine start and end of an NLR. The format has to be NCBI blast with outfmt 5 (xml).
**-a**    | *STR*    | Location of an optional amino acid sequence file. This file should be the same as the one subjected to MAST. Providing this file allows extraction of the NB-ARC domain of the NLR, e.g. for phylogenetic studies. File has to be fasta format.
**-g**    |          | Output gff format instead of a tsv.
**-h**    |          | Print help


### -s splitpattern
In case a nucleotide sequence has to be annotated, it should be translated into its 6 reading frames. The NLR-Parser can assume the sequence names for the 6 amino acid sequences are of a type <common-prefix><splitpattern><framespecific-suffix>. In that case it will report the combined result in one line with <common-prefix> in the first column. It is highly unlikely that a sequence will have motifs in one forward strand and in the reverse strand at the same time. This makes sense if you annotate genomic sequence and introns cause a "frameshift".

This is of course a pit-fall if your sequence of interest contains two NLRs on different strands. In those cases, please use the workaround -s $$, assuming that none of your identifiers contains a "$$".

### -b blastFile
By default, the start- and end-columns of the output point to the beginning of the first motif and the end of the last motif. It is extremely difficult to actually predict the correct start of an NLR gene. You can blast the sequences you put into MAST against a set of well curated NLRs to determine potential start and end. Use this parameter to add the blast to the NLR-Parser report. The format of the file has to be XML and it is assumed that the sequences that were also subjected to MAST are the query sequences in the blast file. For ncbi-blast, use the parameter '-outfmt 5' to generate XML.

### -g 
Generate a gff file rather than a tsv table with the NLR-Parser results. This option is under development. Feel free to try and send us comments.

### -a aminoacidfile.fasta
One column of the NLR-Parser output is the aminoacid sequence of the NB-ARC domain. This is usually the most conserved part of the NLR and can be used for phylogenetic studies. If you do not provide the complete amino acid sequence of the genes, this column is empty.

### -p pvalue
This is the threshold of the p-values of the individual motifs. Motifs with a p-value above this threshold are ignored by the NLR-Parser. The default is **1E-5**. 



## Tips
* MAST has an e-value threshold. Sequences with an evalue above that are not displayed. This evalue is dependent on the number of input sequences. If you run MAST on a really large file, add the parameter `-ev 10000000` to your call. 
* If you want to annotate large files like genomes, it makes sense to chop them in overlapping fragments. 



## Citation
* For using the NLR motifs, please cite [Jupe et al. (2012)](http://www.biomedcentral.com/1471-2164/13/75)
* A publication for the NLR-parser is in preparation. Until then, please use [http://github.com/steuernb/NLR-Parser](http://github.com/steuernb/NLR-Parser) or contact [us](mailto:burkhard.steuernagel@jic.ac.uk) and get us on board.