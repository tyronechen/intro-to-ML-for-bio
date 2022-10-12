## Key Learning Outcomes

After completing this practical the trainee should be able to:

-  Recognize various genomics formats to represent DNA/RNA sequence data 

-  Read in sequence data in FASTA format using Python or BioPython modules 
    decide on filters and cutoffs for cleaning up data ready for
    downstream analysis


***
## Resources You’ll be Using

### Tools Used

Python Package:  
http://

BioPython SeqIO package:  
https://biopython.org/wiki/SeqIO

Google Colaboratory:  
http://

***
##Useful Links

FASTA format: 
https://www.ncbi.nlm.nih.gov/genbank/fastaformat

FASTQ format: 
https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#fastq-files

FASTQ Encoding: 
(http://en.wikipedia.org/wiki/FASTQ_format#Encoding)

***
## Author Information

*Primary Author(s):*    
Sonika Tyagi: sonika.tyagi@monash.edu

*Contributor(s):*    
Sarthak Chauhan, Navya Tyagi, Tyrone Chen

***
## FASTA format 


## FASTQ format

Quality Score Bins | Mapped quality scores
:----------|:-------------
N (no call) | N (no call)
2-9 | 6  
10-19 | 15  
20-24 | 22
25-29 | 27
30-34 | 33
35-39 | 37
\>=40 | 40


**Table 1:** Novaseq Q-score bins mapping

***
## Reading in the sequence file in FASTA format 

A file containing one or more DNA or RNA squences in FASTA format can be read using one of the following two methods. 

1. Using a custom python script

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


2.  Using Biopython utilities

``` 
 from Bio import SeqIO
 for record in SeqIO.parse("example.fasta", "fasta"):
     print(record.id)
```
