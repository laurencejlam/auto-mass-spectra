# Predicting fragments from non-domain-specific mass spectra. 


<img src="https://github.com/user-attachments/assets/856285e0-e0fe-4036-a81d-ae3e0c8b1ac6" width=50%  height=50%>

This is a full scientific workflow to create precise models to extract information from experimental, non-domain-specific mass spectra. The code presented can be adapted to work on alternative data types and other model architectures.

## Data source
The dataset was curated from MassBank NA (https://mona.fiehnlab.ucdavis.edu/). Specifically, experimental LC-MSMS data in JSON format was used. You can use your own datasets as well, but note that the following scripts will need to be adapted for formatting. Tokenised-file-example.json is a short example of the formatting used for the full dataset from MassBank NA.

## Annotate data 
The data needs to be annotated in order to be useful for supervised learning. In the case of mass spectra, this means that spectral peaks need to be annotated with the corresponding fragments. This is done with the following steps:
1. Run streamingjson.py to stream JSON file to a csv, extracting useful information. Streaming method is used to accomodate large datasets.
2. Run inchi_sorted.py to keep track of the rows of all unique parent molecules.
3. Run generate-fragments.py to create a list of all possible fragments (in SMILES) given the set of parent molecules. It uses generates fragments by iteratively breaking bonds in the molecule. Adjust max_cleaves to create more fragments (note that computation time increases exponentially).
4. (OPTIONAL) Run count-number-of-molecules.py to count number of unique parent molecules in dataset.
5. Run tokenize.py to tokenize the parent molecules and the fragments. All annotations must be tokenized to be used in modelling.
6. Run matching-peaks.py to calculate molecular weights and match the correct generated fragments to their corresponding peaks.

## Peak stitching and normalization.

<img src="https://github.com/user-attachments/assets/ddfaf334-fcab-491a-bc8f-635f7105ac12" width=45%  height=45%>

Run peak-stitching-and-normalization.py to stitch together multiple spectra that make up a single data point, as well as normalizing the data using a peak area method to analyse the relative impact of each peak on its spectrum.

## Run model
Run model.py to generate models using hyperparameter sets generated through monte carlo methods. This is a computationally efficient alternative to generating almost optimal hyperparameters. Hyperparameter.csv is an example of a small number of hyperparameter sets.

The model is a 1-D CNN model designed to train on peak intensity data in 1 m/z bins.

## Evaluation
This workflow should output evaluation metrics (and optionally the models themselves) for multiple binary-cross-entropy CNN model to predict fragments from mass spectra. Further selection must be done to select the optimal hyperparameter sets for each fragment.
