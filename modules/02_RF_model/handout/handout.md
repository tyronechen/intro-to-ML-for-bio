# Building a Random Forest model (RF) to classify a set of sequences

The objective of this exercise is to build a classiciation model that is training on two groups of sequences with known labels i.e. "promoter" and "non-promoter". The model is then trained on a held out model to make prediction - whether the given sequence is a promoter or not.
## Learning objectives:

At the end of this practical you should be able to:

- Prepare input data forms for a ML model (RF)
- Investigate parameters used in the model
- Prepare training and test performance metrices and visualization


***
## Resources You’ll be Using

### Data

We have promoter seuqneces of Yeast as one group of data labeled as `1` and a similar synthetic negative data was generated and labelled as `0`.
The data is provided as a file called [`promoter.csv`](../data/promoter.csv) 

### Tools Used

Python Package:  
pandas, numpym sklearn, matplotlib, seaborn

Google Colaboratory:  
https://colab.research.google.com/

***
##Useful Links
Random Classifier
https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html

Confusion matrix:
https://scikit-learn.org/stable/modules/generated/sklearn.metrics.confusion_matrix.html

A read on ROC curves
https://aiineverything.blogspot.com/2021/08/misconception--around-roc-auc.html

***
## Author Information

*Primary Author(s):*    
Sonika Tyagi, Navya Tyagi

*Contributor(s):*    
Tyrone Chen, Sarthak Chauhan 

***

## Workflow

Assuming you have obtained well-defined sequence data and their corresponding metadata:

1.	Read in the data and split each input sequence into "k-mers" or "tokens".
2. 	Denoise the data if needed.[optional]
3. 	Convert these tokens into a numeric representation.
4. 	Split the data into `train`, `test`, and `validation` sets.
5. 	Train, optimize, and validate the model.



*NOTE*: in the context of this tutorial, we will be using the term `token` and `k-mer` interchangeably.



First lets load the required python packages:

```
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score, cross_validate
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report, plot_confusion_matrix
from pprint import pprint
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.feature_extraction.text import CountVectorizer

```

### 1) Read in the data and split each input sequence into "tokens".

We have prepared a data file (`.csv` format) contain sequences from `promoter` (label `1`) and `non-promoter` (label `0`) groups.

```
dna = pd.read_csv('/content/promoter.csv')
dna.head()
```

Lets check the size of the data:

```
color = [sns.xkcd_rgb['medium blue'], sns.xkcd_rgb['pale red']]
sns.countplot('labels',data = dna, palette = color)
plt.gca().set_ylabel('Samples (number of sequences)')

```

Here, is a resuable function to split sequences into k-mers:
```

def get_kmers(sequence: str, length: int):
  """
  Take a dna sequence as input and split the sequence into k-mers / tokens
  """
  return [sequence[x:x+length].upper() for x in range(len(sequence) - length + 1)]

```

Use the above function with the sequence data now stored in a dataframe named `dna`:
```
# in this example, we will split on length 5 k-mers
kmer_length = 
seq=dna.sequence[0]
dna['kmers'] = dna.apply(lambda x: get_kmers(x['sequence'], kmer_length), axis=1)
print(dna.head())
```
We have previously determined the length of `kmer=5` to be used here.

We will next join the individual `kmers` to long string of sentence equivalent. Use the following function:

```
def to_sentence(words):
  r"""This function formats the k-mers such that it is useful for language models where algorithm expects sequences of words"""
  initial=''
  for element in range(len(words)):
    words[element]=' '.join(words[element])
  return words

```
First, extract the `kmers` column as a list and then perform the joining opration.

```
dna_text = list(dna['kmers'])
#get kmers/words joined together with space between them
joined_kmers=to_sentence(dna_text)
print(len(joined_kmers))
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
stopwords = [AAAAA", "TTTTT"]
filtered_kmers = filter_kmers(kmers, stopwords)
print(filtered_kmers)
```

### 3) Convert data into numeric presentation

Here, we are using a frequuency transformation based approach.

The n-gram size of 4 is previously determined by testing
```
cv = CountVectorizer(ngram_range=(,), max_features=500) 
X = cv.fit_transform(joined_kmers)

#get feature names
feature=cv.get_feature_names_out()
# X and Y should have same length 

#separate labels
Y=np.array(dna['labels'])
print("Total data iterms:",X.shape)
print("Total data labels", Y.shape)
```

Take a peek into the count vector we just created:

```
cv.get_feature_names_out()[10:50]
```

You can visualise the token frequencies:

Funtion to visualize frequencies:

```
from yellowbrick.text import FreqDistVisualizer
def token_freq_plot(feature):
  visualizer = FreqDistVisualizer(features=feature, orient='v')
  visualizer.fit(X)
  visualizer.show()
```
```
token_freq_plot(feature)
```

### 4) Split the data into train, test, and validation sets

```
def split_dataset(X, Y, train_ratio, test_ratio, validation_ratio):
    x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size=1 - train_ratio)
    x_val, x_test, y_val, y_test = train_test_split(x_test, y_test, test_size=test_ratio/(test_ratio + validation_ratio))
    return x_train, y_train, x_test, y_test, x_val, y_val

```

```
# split the dataset into train, test and validation set using sklearn
train_ratio = 0.70
validation_ratio = 0.15
test_ratio = 0.15
# train is now 70% of the entire data set
# test is now 15% of the initial data set
# validation is now 15% of the initial data set
x_train, y_train, x_test, y_test, x_val, y_val=split_dataset(X, Y, train_ratio, test_ratio, validation_ratio) 
print("Training data:",x_train.shape)
print("Training data labels:",y_train.shape)
print("Test data:",x_test.shape)
print("Test data labels:",y_test.shape)
print("Validation data:",x_val.shape)
print("Validation data labels:",y_val.shape)
```


Functions to run training and perform assessment

```
def train_model(model, param, x_train, y_train, x_test):
     clf=model(n_estimators=param)
     # fit the training data into the model
     clf.fit(x_train, y_train)
     y_pred=clf.predict(x_test)
     y_probas=clf.predict_proba(x_test)
     return clf, y_pred, y_probas

def model_metrics(model, x_test, y_test, y_pred):
  # accuracy prediction
  accuracy = accuracy_score(y_test, y_pred)
  print("Accuracy: %.2f%%" % (accuracy * 100.0))
  # classification report
  print("Classification report:\n")
  print(classification_report(y_test, y_pred))
  # confusion matrix
  conf=confusion_matrix(y_test, y_pred)
  print("Confusion matrix:\n", conf)


def feature_imp(model,feature_names, n_top_features):
  feats=np.array(feature_names)
  importances = model.feature_importances_
  indices = np.argsort(importances)[::-1]
  plt.figure(figsize=(8,10))
  plt.barh(feats[indices][:n_top_features ], importances[indices][:n_top_features ])
  plt.xlabel("RF feature Importance ")
```


#### Random Forest classication

We have two steps here: 1) we train the model 2) we test it on unseen data
```
print('RF BASE MODEL')
## training the model
model=RandomForestClassifier
param=100
rf_base, y_pred, y_probas=train_model(model, param, x_train, y_train, x_test)
# model metrics
model_metrics(rf_base, x_test, y_test, y_pred)


##Feature imp plot
print("Feature Importance Plot:\n")
feature_imp(rf_base, feature, 10)

```



Lets plot a ROC-AUC plot:

```
from sklearn.metrics import RocCurveDisplay
ax = plt.gca()
rfc_disp = RocCurveDisplay.from_estimator(rf_base, x_test, y_test, ax=ax, alpha=0.8)
#rfc_disp.plot(ax=ax, alpha=0.8)
plt.show()
```
For this exercise we have run the model with default paramters but in real practice we will running and optimisation step where we will find the best paramter settings for pur case study.

We will discuss that topic on another day!

#END
