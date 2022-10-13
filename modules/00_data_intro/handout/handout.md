## Key Learning Outcomes

After completing this practical the trainee should be able to:

-  Recognize various genomics formats to represent DNA/RNA sequence data 

-  Read in sequence data in FASTA format using Python or BioPython modules 
    decide on filters and cutoffs for cleaning up data ready for
    downstream analysis


***
## Resources You’ll be Using

### Tools Used

BioPython SeqIO package:  
https://biopython.org/wiki/SeqIO
https://biopython.org/docs/1.75/api/Bio.SeqIO.QualityIO.html

Google Colaboratory:  
https://colab.research.google.com/
***
##Useful Links

FASTA format: 
https://www.ncbi.nlm.nih.gov/genbank/fastaformat

FASTQ format: 
https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#fastq-files


***
## Author Information

*Author(s):*    
Sonika Tyagi, Sarthak Chauhan, Navya Tyagi, Tyrone Chen

***

DNA sequencing of a biological specimen is performed to find the order of letters (or bases) A, C, G, and T in it to determine the code of DNA. Different organisms have DNA of different length, for example, human DNA is made of 3 billion pairs of individual bases. A highthroughput DNA sequencing experiment can generate miliions of copies of DNA sequences and these are provided as long text files. There are standard formats for presenting DNA sequence data. Two of the most common primary DNA data format are the FASTA (stands for FAST alignment) and FASTQ format.

## FASTA format 

A FASTA format of biological sequence has mainly two fields: a unique sequence identifier line that starts with a `>` sign and the actuall sequence begins from the seconf line.

```
>My_seq_ID
TAATGGCTCT
GGAAGCTCT
TGGCTCTAGA
```

## FASTQ format

Raw data from a highthroughput DNA sequencer is commonly stored in a FASTQ format. This format is similar to FASTA except it contains 4 different fiels and additionally have DNA quality information in it. Here, the sequence identifier line starts with an `@` sign, followed by the actual sequence, a line starting `+` can optionally ahve metadata informaiton, and finally a quality string of the ssame length as the DNA sequence.

```
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC
+SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9IC
```

The fouurth line int eh above format is a quality string where each letter represent a PHRED scale DNA base quality score reperesented on a scale of `0-40+`. 
The table below provides an indication to how quality scores related to base call accuracies.

|Phred Quality Score|	Probability of incorrect base call| Base call accuracy|
|:----------------|:-----------------------|:------|
|10|	1 in 10	|90%|
|20|	1 in 100|	99%|
|30|	1 in 1000|	99.9%|
|40|	1 in 10,000|	99.99%|
|50|	1 in 100,000|	99.999%|
|60|	1 in 1,000,000|	99.9999%|


**Table 1:** Phred quality scores are logarithmically linked to error probabilities

***
## Reading in the sequence file in FASTA format 

A file containing one or more DNA or RNA squences in FASTA format can be read using one of the following two methods. 

### 1) Using a custom python script

Open your Python notebook and type the following code in teh code cell and run.

    def parse_fasta(fh):
       fs=[]
       fss=[]
       seq=''
       for ln in fh:
           if ln[0]=='A' or ln[0]=='T' or ln[0]=='G' or ln[0]=='C':
              seq=seq+ln[:-2]
           else:
              fs.append(seq)
              seq=''
       for element in fs:
         if element!='':
           fss.append(element)
         else:
           fss=fss
       return fss
 
    fh=open('/Sequence.fasta')
    input=parse_fasta(fh)
    print(input)


### 2) Using Biopython utilities

``` 
 from Bio import SeqIO
 for sequence in SeqIO.parse("example.fasta", "fasta"):
     print(sequence.id)
     print(sequence.seq)
     print (len(sequence)
```

A similar function is available to read the FASTQ format as well:

```
from Bio import SeqIO
for record in SeqIO.parse("example.fastq", "fastq"):
    print("%s %s" % (record.id, record.seq))
```
Output will look like this:
```
SeqID1 CCCTTCTTGTCTTCAGCGTTTCTCC
SeqID2 TTGGCAGGCCAAGGCCGATGGATCA
SeqID3 GTTGCTTCTGGCGTGGGTGGGGGGG
```
