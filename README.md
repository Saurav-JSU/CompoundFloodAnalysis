# CompoundFloodAnalysis

This repository provides an R-based analysis of the dependence between storm surge and precipitation. The analysis includes:

- Data processing
- Threshold calculation
- Empirical return period estimation
- Dependence measurement
- Copula fitting
- Joint return period computation

## Directory Structure

- `Data/`: Contains the input data files.
- `FIGURES/`: Contains the generated figures.
- `CSV/`: Contains the output CSV files.
- `R/`: Contains the R scripts for the analysis.

## Sample Data

A sample data file `sample_data.csv` is provided in the `Data` directory to help you get started with the analysis. The file includes columns for date, precipitation, and sea level.

## Running the Analysis

To run the analysis, make sure you have the necessary R libraries installed:

```r
install.packages(c("dplyr", "fitdistrplus", "extRemes", "copula", "ggplot2", "gridExtra"))