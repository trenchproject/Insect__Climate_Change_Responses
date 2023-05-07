# Insect_Responses_Climate_Change
This repository contains the data and scripts for the manuscript "Temperature sensitivity of fitness components across life cycles drives insect responses to climate change"

Folder Structure:
* Biological data: contains field census data from Benin and China and laboratory temperature response data
* Climate data: contains "Climate station data.xlsx" and climate data files for each site in the recent and future climates
* Documentation: contains ReadMe files for downloading Python and climate data as well as running all scripts in "Scripts" folder
* Model parameters: contains "Habitat temperature parameters.csv" and "Temperature response parameters.csv"
* Model predictions: contains model prediction files of climate change effects on fitness metrics/components and population dynamics
* Scripts: contains all scripts for running models and analyses, which are briefly described below:
  * DDE population dynamics.py: Python script for simulating insect population dynamics
  * Download future climate data.py: Python script for accessing and downloading future climate data
  * Fitness metrics and components.R: R script for plotting fitness metrics and components for conceptual Figure 1
  * Habitat temperatures.R: R script for fitting habitat temperature parameters (saved in "Model parameters" folder)
  * Population analyses.R: R script for quantifying climate change effects on population dynamics
  * Read climate data.R: R script for reading downloaded climate data and producing the climate data files in "Climate data" folder
  * Statistical analyses.R: R script for analyzing and plotting fitness metrics/components and population dynamics for Figures 3-5
  * Temperature responses.R: R script for fitting temperature response parameters (saved in "Model parameters" folder)
  * Time series.R: R code for plotting predicted population dynamics for case studies in Figure 2
  * TPC and model analyses.R: R script for quantifying climate change effects on fitness metrics/components
* Time series data: contains density-dependent time-series data predicted by the population model
* Time series DI data: contains density-independent ("DI") time-series data predicted by population model

## Setting up working directories and paths
Scripts in this repository produce and/or read data files that are saved within the folders of the repository. It is therefore important to download the entire repository with all files and folders having the same names, locations, and file extensions as in this repository. As long as the working directory in R or Python matches the main folder of the downloaded repository on the user's computer (i.e., the folder containing "Insect_responses_Climate_Change.rproj"), all paths _should_ work without any changes to the scripts; however, it may be necessary to explicitly specify the paths on the user's computer as detailed in the script's ReadMe file. 

## Using the scripts:
Each script does a specific task as detailed in its ReadMe. It is not necessary to run all scripts to see the results reported in the manuscript. To produce the analyses and plots reported in Figures 3-5, for example, it is only necessary to run "Statistical analyses.R". To produce the conceputal conceptual Figure 1, run "Fitness metrics and components.R" and to produce the time-series plots in Figure 2, run "Time series.R".

To run all scripts used in the manuscript in the correct order, see the order of the ReadMe files in the "Documentation" folder. In general, this involves the following steps:
1. Install Python in order to access future climate data and run DDE population models (see "ReadME1 Install Python")
2. Download recent climate data (see "ReadMe2 Download recent climate data")
3. Run "Download future climate data.py" to download the recent and future climate data
4. Run "Read climate data.R" to produce climate files in "climate data" folder
5. Run "Habitat temperatures.R" to estimate habitat temperature parameters
6. Run "Temperature responses.R" to estimate temperature response parameters
7. Run "DDE population dynamics.py" to simulate the model for four cases involving two factors - recent climate (recent = True) and future climate (recent = false) as well as with competition (competition = True) and without competition (competition = False)
8. Run "Fitness metrics and components.R" to produce the conceptual Figure 1
9. Run "Time series.R" to produce Figure 2
10. Run "TPC and model analyses.R" to quantify climate change effects on fitness metrics/components
11. Run "Population dynamics analyses.R" to quantify climate change effects on population dynamics
12. Run "Statistical analyses.R" to analyze model predictions and produce Figures 3-5
