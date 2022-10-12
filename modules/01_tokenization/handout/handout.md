
# Preprocessing Biological Sequence Data for Natural Language Processing


## Learning objectives:

The overall objective is to understand how biological sequence data is transformed into a machine-readable format for input into a machine learning pipeline.

- View the flow of data through a conventional preproceessing pipeline
- Learn common data transformations
- Familiarity with common representations
- Explore how large sequences are separated into smaller components
- Awareness of the ecosystem of available libraries


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


***
## Author Information

*Primary Author(s):*    
Tyrone Chen, Navya Tyagi, Naima Vahab

*Contributor(s):*    
Sonika Tyagi, Sarthak Chauhan 

***

## Workflow

Assuming you have obtained well-defined sequence data and their corresponding metadata:
1. Split each input sequence into "k-mers" or "tokens".
2. Denoise the data if needed.
3. Convert these tokens into a numeric representation.
4. Transform the numeric representation as needed, depending on your use case.


*NOTE*: in the context of this tutorial, we will be using the term `token` and `k-mer` interchangeably.

Let us take a single biological sequence as an example. Assume that this corresponds to a gene of interest:

```
input_sequence = "TAATGGCTCT"
```

### 1) Split each input sequence into "tokens".

It is impractical in most cases for the machine to obtain a representation of a whole sequence, for example a whole document or chromosome. On the other hand, the same is true for obtaining individual letters of text or genomic DNA. While in the case of some languages such as English, it is possible to exploit the property of space-separation, this is not applicable to genomic DNA. Therefore, in most cases an arbitrary-length split is performed.

```
def get_kmers(sequence: str, length: int):
  """
  Take a dna sequence as input and split the sequence into k-mers / tokens
  """
  return [sequence[x:x+length].upper() for x in range(len(sequence) - length + 1)]

# in this example, we will split on length 3 k-mers
length = 3
kmers = get_kmers(input_sequence, length)
print(kmers)
```

### 2) Denoise the data if needed

In some cases, we want to remove certain tokens, commonly tokens with low information content. In English, this is straightforward and involves filtering out "stopwords", for example `the`, `and` or `is`. In Biology, this is not always done, and will require review on a case-by-case basis.

```
def filter_kmers(tokens: str, stopwords: list):
  """
  Take an input dna sequence and list of stopwords to remove from the data.
  """
  return [x for x in tokens if x not in stopwords]

# in this example, let us pretend this list of k-mers have low information content
stopwords = ["TAA", "AAT"]
filtered_kmers = filter_kmers(kmers, stopwords)
print(filtered_kmers)
```

#### Step 1: Numeric encoding

- Ordinal encoding
- One-hot encoding

The first type of text encodings are ordinal encodings. With this, individual `tokens` are assigned an integer value as shown below:

```
from sklearn.preprocessing import OrdinalEncoder

def encode_ordinal(input_sequence: str, length: int):
  """
  Take a list of k-mers and perform ordinal encoding
  """
  tokens = [[x] for x in get_kmers(input_sequence, length)]
  encoder = OrdinalEncoder()
  return encoder.fit_transform(tokens)

ordinal = encode_ordinal(input_sequence, length)
print(ordinal)
```

> *NOTE*: As an alternative to ordinal encoding, it is also possible to encode tokens as one-hot encodings. This method is not practical to use in most cases due to its sparsity, but we demonstrate an example for completeness. Compare this to the above:

```
[[ 0 0 0 0 0 1 0 0]
 [ 1 0 0 0 0 0 0 0]
 [ 0 1 0 0 0 0 0 0]
 [ 0 0 0 0 0 0 0 1]
 [ 0 0 0 0 1 0 0 0]
 [ 0 0 0 1 0 0 0 0]
 [ 0 0 1 0 0 0 0 0]
 [ 0 0 0 0 0 0 1 0]]
```


##### Frequency tables

**Count vectorisation**

First, each token is counted and these values are assigned to a dictionary. On its own, these counts can be vectorised directly. However, it is often practical to weight tokens by their frequency as shown in the next step:

**TF-IDF (Term Frequency-Inverse Document Frequency)**

TF-IDF is analogous to a library-size correction in conventional short-read RNA-Seq, allowing different size sequence collections to be compared on similar scales. It is computed on the previous counts prior to vectorisation.

TF-IDF can be calculated as follows:

```
tf-idf(t, d) = tf(t, d) * log(N/(df + 1))
```

Where:
- t = token
- d = sequence
- N = total number of sequences in the library

Within `scikit-learn`, the counting process is abstracted and not shown directly. We demonstrate this intermediate step for reference:


```
# input
'TAA AAT ATG TGG GGC GCT CTC TCT'

# intermediate (unweighted matrix)
{'TAA': 1, 'AAT': 1, 'ATG': 1, 'TGG': 1, 'GGC': 1, 'GCT': 1, 'CTC': 1, 'TCT': 1}

# resulting matrix
```
```
from sklearn.feature_extraction.text import TfidfVectorizer

def tfidf(input_sequence: str, length: int):
  """
  Take a dna sequence and k-mer size as input, 
  output a matrix of TF-IDF features.
  """
  data = [" ".join(get_kmers(input_sequence, length))]
  tfidf = TfidfVectorizer()
  tfidf = tfidf.fit_transform(data)
  return tfidf.toarray()

tfidf_matrix = tfidf(input_sequence,length) 
print(tfidf_matrix)

```


Although the frequency-based approaches above are often effective on their own, one key piece of information is lost. Individual tokens are treated as independent occurrences, which is not realistic in both natural language and biology. To account for the semantic relationship between tokens, a different approach is required.

##### Embedding vector

A latent space is generated from all tokens in the corpus. This n-dimensional `embedding` accounts for every combination of tokens in the data. Taking a single sequence (i.e. group of tokens) and projecting this back onto the embedding returns a unique vector for each sequence. Examining these further shows that tokens which are closely related to each other are naturally recapitulated by the model. For example, in a real life English dataset, you may expect the following word pairs to be closely related to each other:


#### IMAGE HERE

```
import gensim
from gensim.models import Word2Vec

def create_word2vec(input_sequence: str, length: int):
  """
  Input a list of tokens to generate an embedding of word vectors
  """
  data = get_kmers(input_sequence, length)
  return Word2Vec([data], min_count=1)

def map_token(input_sequence: str, model: gensim.models.base_any2vec):
  """
  Input a token and embedding and map the token onto the embedding
  """
  return model[input_sequence]

embedding = create_word2vec(input_sequence, length)
vector = map_token('TAA', embedding)
print(vector)
```

While we demonstrate `gensim` as one example, many embedding algorithms exist. A comprehensive overview of each method is out of the scope of this tutorial, but we list a few popular methods for your reference:

- Word2Vec
- dna2vec (dna equivalent of the above)
- GloVe

# END
