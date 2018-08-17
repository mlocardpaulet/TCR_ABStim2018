# TCR_ABStim2018
Phosphoproteomics analysis of antibody-based stimulation of TCR in CD4+ T cells

## Folder organisation:

The folder "RAW" contains all the search results from MaxQuant:
- "DeepProteome": deep proteome generated from an in-gel fractionated sample (SDS-PAGE) to have a maximum coverage in terms of the cells protein IDs.
- "PeptideSamples": samples analysed before phospho-enrichment.
- "TiO2": 10% of the samples analysed after TiO2 and before phospho-tyrosine IP.
- "pYIP": samples analysed after TiO2 and phospho-tyrosine IP.

## Data analysis:

The data was prepared and normalised with the scripts in "DataPreparation.Rmd/html".
