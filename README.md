# TCR_ABStim2018
Phosphoproteomics analysis of antibody-based stimulation of TCR in CD4+ T cells

## Folder organisation:

`"RAW"` contains all the search results from MaxQuant:
- "DeepProteome": deep proteome generated from an in-gel fractionated sample (SDS-PAGE) to have a maximum coverage in terms of the cells protein IDs.
- "PeptideSamples": samples analysed before phospho-enrichment.
- "TiO2": 10% of the samples analysed after TiO2 and before phospho-tyrosine IP.
- "pYIP": samples analysed after TiO2 and phospho-tyrosine IP.
- "iRTNorm": peptide tables from the MaxQuant search of all the phospho-enriched runs to normalise the signal for instrument variation using spiked-in peptides.

`"PhosphoSitePlus_20180817"` contains the data downloaded from the PhosPhositePlus database (phosphosite.org) the 17/08/2018.

`"RFunctions"` contains the functions used for the analysis.

`"Figures"` contains the ouptput figures.

`"RData"` contains the data at different steps of the processing, and mapping tables or other information needed for the analysis:
- "DBMapping": data used for the mapping of our phosphorylation sites to the PhosphoSitePlus data base (see the file "PhosphoDBMapping.Rmd" to get the details).
- "01_ParsedTables.RData": Phospho(STY) tables from MaxQuant after parsing through the function "MultiPhosphoParsingMQ" to get the ID of the multi-phosphorylated peptides.
- "02_psitesInfo.RData": list of all the phosphosites of the analysis with their different IDs.
- "03_psiteInfoPhosphoSitePlus.RData": list of all the IDs of the phosphosites of the data set, with the phosphositePlus references (Mouse and all organisms if no mouse found)(see the file "PhosphoDBMapping.Rmd" to get the details). 
- "04_ParsedTablesPSP.RData": Quantification data with the PhosphoSitePlus ID.


## Data analysis:

The data was prepared and normalised with the scripts in "DataPreparation.Rmd/html".
