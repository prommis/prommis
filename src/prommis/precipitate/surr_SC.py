# Import statements
import os
import numpy as np
import pandas as pd

# Import IDAES libraries
from idaes.core.surrogate.sampling.data_utils import split_training_validation
from idaes.core.surrogate.pysmo_surrogate import PysmoPolyTrainer, PysmoSurrogate
from idaes.core.surrogate.plotting.sm_plotter import (
    surrogate_scatter2D,
    surrogate_parity,
    surrogate_residual,
)
# Import training data
np.set_printoptions(precision=6, suppress=True)

Species = ["Al", "Ca", "Fe", "Sc", "Y", "La", "Ce", "Pr", "Nd", "Sm", "Gd", "Dy"]
# Species = ["Al"]
csv_data = pd.read_csv("all_species.csv")
csv_data.columns.values[0:3] =["species", "pH","sep"]

for j in Species:
    csv_data1 = csv_data[csv_data.species  == j] 
    print(csv_data1)

    data = csv_data1.sample(n=6)

    input_data = data.iloc[:, 1:2]
    output_data = data.iloc[:, 2:3]

    # # Define labels, and split training and validation data
    input_labels = list(input_data.columns)
    output_labels =  list(output_data.columns) 

    n_data = data[input_labels[0]].size
    data_training, data_validation = split_training_validation(
        data, 0.8, seed=n_data
    )

    # Create PySMO trainer object
    trainer = PysmoPolyTrainer(
        input_labels=input_labels,
        output_labels=output_labels,
        training_dataframe=data_training,
    )

    var = output_labels
    trainer.config.extra_features=['pH*pH']
    # Set PySMO options
    trainer.config.maximum_polynomial_order = 3
    trainer.config.multinomials = True
    trainer.config.training_split = 0.8
    trainer.config.number_of_crossvalidations = 3

    # Train surrogate (calls PySMO through IDAES Python wrapper)
    poly_train = trainer.train_surrogate()

    # create callable surrogate object
    xmin, xmax = [0.1,8], [10,100]
    input_bounds = {input_labels[i]: (xmin[i], xmax[i]) for i in range(len(input_labels))}
    poly_surr = PysmoSurrogate(poly_train, input_labels, output_labels, input_bounds)
    # save model to JSON
    model = poly_surr.save_to_file('pysmo_poly_surrogate'+str(j)+'.json', overwrite=True)

    # # visualize with IDAES surrogate plotting tools
    # surrogate_scatter2D(poly_surr, data_training, filename="pysmo_poly_train_scatter2D.pdf")
    # surrogate_parity(poly_surr, data_training, filename="pysmo_poly_train_parity.pdf")
    # surrogate_residual(poly_surr, data_training, filename="pysmo_poly_train_residual.pdf")