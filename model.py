import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers, models
from tensorflow.keras import metrics as keras_metrics
from tensorflow.keras import backend as K
#from tensorflow.keras import regularizers
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import json
import sklearn
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import gc



keras.backend.clear_session()

# Check if GPU is available and limit TensorFlow to allocate memory as needed
gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    try:
        # Set memory growth to avoid allocating all GPU memory at once
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
    except RuntimeError as e:
        print(e)

mpl.rcParams['figure.figsize'] = (12, 10)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# Load the JSON data
with open('normalised_data5-19.json', 'r') as file:
    data = json.load(file)

# Number of m/z values
totalmz = 1001

# Load the tokens from the CSV
tokens_df = pd.read_csv('updated_fragment_counts.csv')
# Assuming the tokens are in the third column and the CSV has headers
# Replace 'Token' with the actual column name if it has one

#This sets the hyperparameter set to start with.
start = 45

#This selects the hyperparameter rows to iterate through.
top_30_tokens = tokens_df.iloc[start:60, 2].tolist()

# Initialize an empty list to collect all results
all_results = []

# Path to save the results
results_csv_path = 'hyper/30to100-top100params-24dec.csv'

# Initialize the results DataFrame
if not os.path.exists(results_csv_path):
    results_df = pd.DataFrame(columns=[
        'token',
        'learning_rate', 'epochs', 'batch_size', 'dropout_rate', 'pos_weight', 
        'conv1_neurons', 'conv2_neurons', 'conv3_neurons', 'dense_neurons',
        'binary_crossentropy', 'mean_squared_error', 'true_positives', 'false_positives', 
        'true_negatives', 'false_negatives', 'binary_accuracy', 'precision', 'recall', 
        'auc', 'f1_score'
    ])
    results_df.to_csv(results_csv_path, index=False)
else:
    results_df = pd.read_csv(results_csv_path)

# Load hyperparameters from CSV
hyperparams_df = pd.read_csv('hyper/hyper70.csv')
#l2_factor = 0.0001  # Typical values: 0.001 or 0.0001
start = start-1
# Loop over each token
for token in top_30_tokens:
    start = start + 1
    print(f"Processing token: {token} {start}")
    
    # Initialize lists to store data
    intensity_values = {i: [] for i in range(1, totalmz)}  # Create intensity variables for m/z values from 1 to 1000
    class_labels = []
    
    # Process each entry in the JSON data
    for entry in data:
        spectrum = entry['Spectrum']
        rounded_mz_values = {int(round(float(peak['m/z value']))): float(peak['Rel. Intensity']) for peak in spectrum}
        
        # Set intensity values based on presence in the spectrum
        for i in range(1, totalmz):
            if i in rounded_mz_values:
                intensity_values[i].append(rounded_mz_values[i])
            else:
                intensity_values[i].append(0)  # Set intensity to 0 if no peak exists
        
        # Check if the current token exists in the spectrum
        token_exists = any(str(peak.get('Token')) == str(token) for peak in spectrum)
        class_labels.append(1 if token_exists else 0)  # Assign class label
    
    # Create DataFrame
    raw_df = pd.DataFrame(intensity_values)
    raw_df['Class'] = class_labels
    
    # Augment dataset to balance classes
    true_samples = raw_df[raw_df['Class'] == 1]
    other_samples = raw_df[raw_df['Class'] == 0]
    
    num_true_samples_needed = len(true_samples)
    if num_true_samples_needed == 0:
        print(f"No positive samples for token {token}. Skipping.")
        continue  # Skip tokens with no positive samples
    balanced_other_samples = other_samples.sample(n=num_true_samples_needed, random_state=42)
    balanced_df = pd.concat([true_samples, balanced_other_samples])
    
    # Shuffle the dataframe
    balanced_df = balanced_df.sample(frac=1).reset_index(drop=True)
    
    # Split into training, validation, and test sets
    train_df, temp_df = train_test_split(balanced_df, test_size=0.3, random_state=42)
    val_df, test_df = train_test_split(temp_df, test_size=0.5, random_state=42)
    
    # Prepare features and labels
    train_labels = np.array(train_df.pop('Class')).reshape(-1, 1)  # Reshape to (num_samples, 1)
    val_labels = np.array(val_df.pop('Class')).reshape(-1, 1)
    test_labels = np.array(test_df.pop('Class')).reshape(-1, 1)
    
    train_features = np.array(train_df)
    val_features = np.array(val_df)
    test_features = np.array(test_df)
    
    # Scaling
    scaler = StandardScaler()
    train_features = scaler.fit_transform(train_features)
    val_features = scaler.transform(val_features)
    test_features = scaler.transform(test_features)
    
    train_features = np.clip(train_features, -5, 5)
    val_features = np.clip(val_features, -5, 5)
    test_features = np.clip(test_features, -5, 5)
    
    # Model metrics
    METRICS = [
          keras.metrics.BinaryCrossentropy(name='cross_entropy'),
          keras.metrics.MeanSquaredError(name='Brier score'),
          keras.metrics.TruePositives(name='tp'),
          keras.metrics.FalsePositives(name='fp'),
          keras.metrics.TrueNegatives(name='tn'),
          keras.metrics.FalseNegatives(name='fn'), 
          keras.metrics.BinaryAccuracy(name='accuracy'),
          keras.metrics.Precision(name='precision', dtype='float32'),
          keras.metrics.Recall(name='recall', dtype='float32'),
          keras.metrics.AUC(name='auc'),
          keras.metrics.AUC(name='prc', curve='PR'),
    ]
    
    def weighted_binary_crossentropy(pos_weight):
        def loss(y_true, y_pred):
            bce = tf.keras.losses.binary_crossentropy(y_true, y_pred)
            weighted_bce = tf.reduce_mean(pos_weight * bce)
            return weighted_bce
        return loss
    
    def make_model(learning_rate, dropout_rate, pos_weight, 
                   conv1_neurons, conv2_neurons, conv3_neurons, dense_neurons, metrics=METRICS, output_bias=None):
        if output_bias is not None:
            output_bias = tf.keras.initializers.Constant(output_bias)
            

        model = tf.keras.Sequential([
            tf.keras.layers.Reshape((1000, 1), input_shape=(1000,)),
            tf.keras.layers.Conv1D(conv1_neurons, 3, activation='relu'),
            tf.keras.layers.MaxPooling1D(2),
            tf.keras.layers.Conv1D(conv2_neurons, 3, activation='relu'), #,                            kernel_regularizer=regularizers.l2(l2_factor)\
            tf.keras.layers.MaxPooling1D(2),
            tf.keras.layers.Conv1D(conv3_neurons, 3, activation='relu'),
            tf.keras.layers.MaxPooling1D(2),
            tf.keras.layers.Flatten(),
            tf.keras.layers.Dense(dense_neurons, activation='relu'),
            tf.keras.layers.Dropout(dropout_rate),
            tf.keras.layers.Dense(1, activation='sigmoid', bias_initializer=output_bias)
        ])
    
        model.compile(
            optimizer=tf.keras.optimizers.Adam(learning_rate=learning_rate),
            loss=tf.keras.losses.BinaryCrossentropy(), #loss=tf.binary_crossentropy(pos_weight=pos_weight),  #
            metrics=metrics)
        
        return model
    start_row_index = 0  # Replace with the desired starting index


    # Your code here
    # Loop over hyperparameters
    for index, row in hyperparams_df.iloc[start_row_index:].iterrows():
        learning_rate = float(row.iloc[0])
        epochs = int(row.iloc[1])
        batch_size = int(row.iloc[2])
        dropout_rate = float(row.iloc[3])
        pos_weight = int(row.iloc[4])
        conv1_neurons = int(row.iloc[5])
        conv2_neurons = int(row.iloc[6])
        conv3_neurons = int(row.iloc[7])
        dense_neurons = int(row.iloc[8])
        
        
        for run in range(1, 5):
            print(f"Run {run}, set {index +1}, token {token}, start {start}")
            # Build the model
            model = make_model(learning_rate, dropout_rate, pos_weight, 
                            conv1_neurons, conv2_neurons, conv3_neurons, dense_neurons)
            
            early_stopping = tf.keras.callbacks.EarlyStopping(
                monitor='val_prc', 
                verbose=0,
                patience=15,
                mode='max',
                restore_best_weights=True)
            
            try:
                history = model.fit(
                    train_features,
                    train_labels,
                    batch_size=batch_size,
                    epochs=epochs,
                    validation_data=(val_features, val_labels),
                    verbose = 0,
                    callbacks=[early_stopping])
            except tf.errors.ResourceExhaustedError as e:
                print(f"ResourceExhaustedError for token {token} with hyperparameters set {index+1}: {e}")
                # Clean up and continue
                K.clear_session()
                del model
                del history
                gc.collect()
                continue  # Skip to the next set of hyperparameters
            
            # Make predictions
            predictions = model.predict(test_features, batch_size=batch_size, verbose=0)
            predictions_binary = (predictions > 0.5).astype(int)
            
            # Evaluate metrics
            binary_crossentropy = keras_metrics.BinaryCrossentropy()
            mean_squared_error = keras_metrics.MeanSquaredError()
            true_positives = keras_metrics.TruePositives()
            false_positives = keras_metrics.FalsePositives()
            true_negatives = keras_metrics.TrueNegatives()
            false_negatives = keras_metrics.FalseNegatives()
            binary_accuracy = keras_metrics.BinaryAccuracy()
            precision = keras_metrics.Precision()
            recall = keras_metrics.Recall()
            auc = keras_metrics.AUC()
            
            binary_crossentropy.update_state(test_labels, predictions)
            mean_squared_error.update_state(test_labels, predictions)
            true_positives.update_state(test_labels, predictions_binary)
            false_positives.update_state(test_labels, predictions_binary)
            true_negatives.update_state(test_labels, predictions_binary)
            false_negatives.update_state(test_labels, predictions_binary)
            binary_accuracy.update_state(test_labels, predictions_binary)
            precision.update_state(test_labels, predictions_binary)
            recall.update_state(test_labels, predictions_binary)
            auc.update_state(test_labels, predictions)
            
            # Calculate F1 score from precision and recall
            precision_value = float(precision.result().numpy())
            recall_value = float(recall.result().numpy())
            f1_score = 2 * (precision_value * recall_value) / (precision_value + recall_value) if (precision_value + recall_value) > 0 else 0
            
            result = {
                'token': token,
                'learning_rate': learning_rate,
                'epochs': epochs,
                'batch_size': batch_size,
                'dropout_rate': dropout_rate,
                'pos_weight': pos_weight,
                'conv1_neurons': conv1_neurons,
                'conv2_neurons': conv2_neurons,
                'conv3_neurons': conv3_neurons,
                'dense_neurons': dense_neurons,
                'binary_crossentropy': binary_crossentropy.result().numpy(),
                'mean_squared_error': mean_squared_error.result().numpy(),
                'true_positives': true_positives.result().numpy(),
                'false_positives': false_positives.result().numpy(),
                'true_negatives': true_negatives.result().numpy(),
                'false_negatives': false_negatives.result().numpy(),
                'binary_accuracy': binary_accuracy.result().numpy(),
                'precision': precision.result().numpy(),
                'recall': recall.result().numpy(),
                'auc': auc.result().numpy(),
                'f1_score': f1_score
            }
            
            all_results.append(result)
            
            # Update the CSV file after each iteration
            results_df = pd.concat([results_df, pd.DataFrame([result])], ignore_index=True)
            results_df.to_csv(results_csv_path, index=False)
            
            #print(f"Finished training for token: {token} with hyperparameters set {index+1}")
            #print("Test results:", result)
        
            # Clear memory
            K.clear_session()
            del model
            del history
            del predictions
            del predictions_binary
            del binary_crossentropy
            del mean_squared_error
            del true_positives
            del false_positives
            del true_negatives
            del false_negatives
            del binary_accuracy
            del precision
            del recall
            del auc
            gc.collect()
        
    # Clear data variables for the current token
    del raw_df
    del balanced_df
    del train_df
    del val_df
    del test_df
    del train_features
    del val_features
    del test_features
    del train_labels
    del val_labels
    del test_labels
    del scaler
    gc.collect()

# Optionally, save all results at the end
# results_df = pd.DataFrame(all_results)
# results_df.to_csv(results_csv_path, index=False)
keras.backend.clear_session()
