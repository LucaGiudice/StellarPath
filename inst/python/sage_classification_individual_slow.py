# verify that we're using the correct version of StellarGraph for this notebook
import stellargraph as sg
import warnings
warnings.filterwarnings("ignore")
import os
os.environ['TF_CPP_MIN_LOG_LEVEL']='3'
import tensorflow as tf

try:
    sg.utils.validate_notebook_version("1.2.1")
except AttributeError:
    raise ValueError(
        f"This notebook requires StellarGraph version 1.2.1, but a different version {sg.__version__} is installed.  Please see <https://github.com/stellargraph/stellargraph/issues/1172>."
    ) from None
    
import networkx as nx
import pandas as pd
import stellargraph as sg
from stellargraph.mapper import GraphSAGENodeGenerator
from stellargraph.layer import GraphSAGE
from tensorflow.keras import layers, optimizers, losses, Model
from tensorflow.keras.callbacks import EarlyStopping
from sklearn import preprocessing
from sklearn.metrics import matthews_corrcoef

def node_classification(allNames,nodes_all,allEdges):
    #Cut the metainformation
    nodes = nodes_all[['IDs','Groups']]
    #Generate the Series with all the nodes indipendetly by the set in which belong
    nodes_serie = pd.Series(nodes['Groups'].values, name='Groups', index=nodes['IDs'])
    nodes_hot=nodes.set_index("IDs")
    #Hot encoding with all the nodes indipendetly by the set in which belong
    nodes_hot = pd.get_dummies(
        nodes_hot, columns=["Groups"]
    )
    #print("hot encoding of all nodes:\n",nodes_hot.head())

    #Generate the Series related to the node's sets
    train_nodes = nodes_all.loc[nodes_all['train_indxs'] == 1]
    train_nodes = train_nodes[['IDs','Groups']]
    train_subjects = pd.Series(train_nodes['Groups'].values, name='Groups', index=train_nodes['IDs'])

    val_nodes = nodes_all.loc[nodes_all['val_indxs'] == 1]
    val_nodes = val_nodes[['IDs','Groups']]
    val_subjects = pd.Series(val_nodes['Groups'].values, name='Groups', index=val_nodes['IDs'])

    test_nodes = nodes_all.loc[nodes_all['test_indxs'] == 1]
    test_nodes = test_nodes[['IDs','Groups']]
    test_subjects = pd.Series(test_nodes['Groups'].values, name='Groups', index=test_nodes['IDs'])

    #print("train subjs:\n",train_subjects.value_counts().to_frame())
    #print("test subjs:\n",test_subjects.value_counts().to_frame())
    #print("validation subjs:\n",val_subjects.value_counts().to_frame())

    #Generate the encoding from labels to bin and the targets lists
    train_v=train_subjects.values.tolist()
    test_v=test_subjects.values.tolist()
    val_v=val_subjects.values.tolist()

    train_v = [(x,) for x in train_v]
    test_v = [(x,) for x in test_v]
    val_v = [(x,) for x in val_v]

    target_encoding = preprocessing.MultiLabelBinarizer()
    train_targets = target_encoding.fit_transform(train_v)
    val_targets = target_encoding.transform(val_v)
    test_targets = target_encoding.transform(test_v)

    #print("testing list:\n",test_v)
    #print("testing targs:\n",test_targets)
    
    #Set model params
    batch_size = 50
    num_samples = [10, 5]

    #Instantiate resulting vectors
    testing_res_mcc = []
    validation_res_mcc = []
    predictions_val_ll = []
    predictions_test_ll = []
    network_names = []
    
    for count in range(0,len(allNames)):
        network_name = allNames[count]
        #print("PSN: \n",network_name)
        edges = allEdges[count]
        
        #Create stellargraph object
        G = sg.StellarGraph({"samples": nodes_hot}, {"connections": edges})
        
        #Generate the model, add the node sets
        generator = GraphSAGENodeGenerator(G, batch_size, num_samples)
        train_gen = generator.flow(train_subjects.index, train_targets)
        val_gen = generator.flow(val_subjects.index, val_targets)
        
        #Create the model
        gs = GraphSAGE(
            layer_sizes=[32, 32], generator=generator, bias=True, dropout=0,
        )
            
        x_inp, x_out = gs.in_out_tensors()
        predictions = layers.Dense(units=train_targets.shape[1], activation="softmax")(x_out)
        model = Model(inputs=x_inp, outputs=predictions)
        model.compile(
            optimizer=optimizers.Adam(learning_rate=0.01),
            loss=losses.categorical_crossentropy,
            metrics=["acc"],
        )

        #Set the approaching for reaching the end of training
        es_callback = EarlyStopping(monitor="val_acc", patience=10, restore_best_weights=True)

        #Train
        model.fit(
            train_gen,
            epochs=200,
            validation_data=val_gen,
            verbose=0,
            shuffle=False,
            callbacks=[es_callback],
        )
       
        #Predict validation set
        val_predictions = model.predict(val_gen)
        
        #Convert binary prediction to class
        val_predictions_squeezed = val_predictions.squeeze()
        transformer = preprocessing.Binarizer(threshold=0.5).fit(val_predictions_squeezed)
        val_predictions_bin = transformer.transform(val_predictions_squeezed)
        val_predictions = target_encoding.inverse_transform(val_predictions_bin)
            
        #Save
        val_mcc = matthews_corrcoef(val_v, val_predictions)
        validation_res_mcc.append(val_mcc)
        predictions_val_ll.append(val_predictions)
    
        #Add the testing set
        test_gen = generator.flow(test_subjects.index, test_targets)
        
        #Predict the testing set
        test_predictions = model.predict(test_gen)
        
        #Convert binary prediction to class
        test_predictions_squeezed = test_predictions.squeeze()
        transformer = preprocessing.Binarizer(threshold=0.5).fit(test_predictions_squeezed)
        test_predictions_bin = transformer.transform(test_predictions_squeezed)
        test_predictions = target_encoding.inverse_transform(test_predictions_bin)
        
        #Save
        test_mcc = matthews_corrcoef(test_v, test_predictions)
        testing_res_mcc.append(test_mcc)
        predictions_test_ll.append(test_predictions)
        network_names.append(network_name)

        
    #Summary of performances and predictions
    df_metrics = pd.DataFrame({"PSNs": network_names, "validation_MCC": validation_res_mcc, "test_MCC": testing_res_mcc})
    df_val_predictions = pd.DataFrame(predictions_val_ll)
    df_test_predictions = pd.DataFrame(predictions_test_ll)

    #Write
    #df_metrics.to_excel(metrics_path)
    #df_val_predictions.to_excel(val_preds_path)
    #df_test_predictions.to_excel(test_preds_path)

    #print(df_metrics)
    return df_metrics, df_val_predictions, df_test_predictions
