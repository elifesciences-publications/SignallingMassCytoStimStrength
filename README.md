This code is associated with the paper from Ma et al., Stimulation strength controls the rate of initiation but not the molecular organization of TCR-induced signalling. eLife, 2020. http://doi.org/10.7554/eLife.53948



# SignallingCyTOFStimStrength

This repo contains code for analysis of mass cytometry data looking at T cell signalling after stimulating with different ligands as described in Ma et al, 2020.

## Data download

Data can be downloaded from the FlowRepository (https://flowrepository.org/), accession numbers FR-FCM-Z2CX and FR-FCM-Z2CP. 

## Scripts

flow_data_bubble_plot.R creates a bubble plot from flow cytometry analyses depicting: a) percentage of cells positive for each marker, and b) median fluorescent intensity of the positive population. Inputs are annotated within the function.

MC_NormToDNA.Rmd and MC_NormToTotalProtein.Rmd present functions to normalize the signal of chosen channels in each cell to the DNA signal in that cell or, in the case of phospho-proteins, to the signal from each corresponding total protein in that cell.

Timecourse_peptides_analysis.Rmd runs through batch correction and differential abundance analyses of FR-FCM-Z2CP presented in Ma et al. Running this script assumes data from FR-FCM-Z2CP has been downloaded into a directory called "data" created within this respository directory.

Timecourse_peptides_plots.Rmd generates plots from analyses in Timecourse_peptides_analysis.Rmd.

Timecourse_peptides_write_results.R writes a table of statstics from analyses in Timecourse_peptides_analysis.Rmd.

Timecourse_trajectories_full.Rmd constructs activation trajectories and analyzes the order of activation events in FR-FCM-Z2CP presented in Ma et al. Running this script assumes data from FR-FCM-Z2CP has been downloaded into a directory called "data" created within this respository directory.

Timecourse_trajectories_012.Rmd constructs activation trajectories and analyzes the order of activation events in biological replicates containing data for 0, 1 and 2 hours of stimulation from FR-FCM-Z2CX presented in Ma et al. Running this script assumes data from FR-FCM-Z2CX has been downloaded into a directory called "data1" created within this respository directory.


