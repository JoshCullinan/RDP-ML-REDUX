# RDPML: Machine Learning for Viral Recombination Detection

This repository contains the code and datasets used in my MSc Bioinformatics thesis at the University of Cape Town, focused on enhancing viral recombination detection through machine learning approaches.

## Project Overview

The project applies modern machine learning techniques to improve upon existing viral recombination detection methods. Using SANTA-SIM generated viral evolution data, we developed and evaluated multiple computational approaches against the current standard, RDP5. The study trained and tested several models including:

- Logistic regression
- Gradient boosting (LightGBM)
- Random forests
- Neural networks with both binary and position-selection architectures

The neural network employing position selection achieved the highest performance with a weighted AUC of 0.784, surpassing RDP5's baseline AUC of 0.739.

## Repository Contents

- Custom SANTA-SIM implementation for recombination tracking
- Dataset generation and preprocessing scripts
- Model training and evaluation code
- Trained model weights

## Requirements

- Python 3.12.3
- TensorFlow 2.18.0
- SciKit-Learn 1.5.2
- Full requirements listed in requirements.txt

## Dataset

The training dataset consists of 491 124 sequences for classification, derived from multiple viral evolution simulations under varying parameters. A separate test dataset based on SARS-CoV-2 was used for validation.

For access to these datasets please contact Joshua Cullinan at CLLJOS001@myuct.ac.za

## File Descriptions

> SantaSim Files
*santa.jar* -> modified SantaSim executable file
*Simulation.py* -> pipeline to parse the XML files and run santasim simultaneously with multiple variable inputs.

> Event Classifier 
*event_classifier.py* -> Written to process the raw simulated genomic sequence data to identify optimal minor and major parent sequences for recombinant sequences, used in the custom version of RDP5.
*event_classifier_pipeline.py* -> Pipeline for event_classifier

> RDP Pipeline
*RDP_pipeline.py* -> scans through a supplied directory and runs RDPCL.exe with the files generated from the custom version of SantaSim, outputs the raw training statistics. 

> ML Data Parser
*output_parser.py* -> Used to process all of the RDP5 statistics with the santa sim output files to create the datasets used for machine learning.

> Machine Learning
*tools.py* -> contains frequently used functions across all the Jupyter notebooks.
*RDPML.ipynb* -> Jupyter notebook for the initial models trained on the data.
*RDPML_BNN.ipynb* -> Jupyter notebook for the binary approach neural network.
*RDP_PSNN.ipynb* -> Jupyter notebook for the posistion selection neural network.

## License

This project is licensed under the MIT License - see the LICENSE file for details.