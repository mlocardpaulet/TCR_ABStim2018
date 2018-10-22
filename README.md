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
- "05_ParsedTablesNormPSP.RData": Data normalised for technical instrument variation using the intensities of the iRT. Also, I removed some runs that overperformed. 
- "06_Prottab.RData": Protein data (raw).
- "07_TableBeforeStat.RData": The normalised and QC data of the phosphoproteome before statistical analysis.
- "08_StatResultsPhospho.RData": The output of the statistical analysis of the phosphoproteome.
- "09_Correlations.RData": The output of the correlation analysis of the significantly regulated sites.
- "10_Cluster.RData": The clusters of phospho-kinetics (only include the significantly regulated sites).
- "11_TableClusters.RData": Phosphoproteomic data with the results of the statistical analysis and the clustering.
- "12_PhosphoTableWithProteinWaring.RData": Phosphoproteomic data with the results of the statistical analysis and the clustering. I added a final column that indicates if the protein FC correlates with the phosphorylation sites FC on this protein.
- "13_PhosphoscoutMapping.RData": Mapping of the phospho-regulated sites to the PhosphoScout database.

`"KEAEnrichments"` contains the results of phosphorylation site specific enrichment using [KEA](http://www.maayanlab.net/KEA2/#) the 4th September 2018 (no background). I saved the enrichment results with the following labels:
- "KS" for Literature Based Kinase-Substrate Library with Phosphosites, followed by the cluster number.
- "BT" for Biological Terms Associated with Phosphosites from Literature Mining, followed by the cluster number.
- "AllReg" is when I used all the regulated sites as input.
- "Background" is when I used all the sites identified in the data set.

`"PANTHER"` contains the results of GO term enrichments based on the gene names using [Panther](http://pantherdb.org/) the 4th September 2018. I used the detected sites as background.

`"Phomics"` contains the results of GO term enrichments based on the phosphorylation sites using [Phomics](http://phomics.jensenlab.org/phospho_enrichment_uniprot) the 4th September 2018. I used the detected sites as background.

`"tsne"` contains the scripts for tsne prodution and clustering (work of Romain Roncagalli).

## Data analysis:

The data was prepared and normalised with the scripts in "DataPreparation.Rmd/html" with some quality control. See "DataPreparation_Proteome.Rmd/html" for the quality control of the proteins data.

The scripts used for the statistical testing of the data (including kinetics normalisation and replacement of missing values) are in the folder `StatisticalAnalysis`. IT CONTAINS THE MANUAL CORRECTION OF SEVERAL NAMING ISSUES IN THE ORIGINAL TABLES. The same folder contains the script used for the statistical analysis of the proteome.

The clustering of the phosphorylation sites is performed in the document `PhoshositesClustering`. I chose to perform tight clustering on the mean values of the loops. I scaled the kinetics (row-wise).

I made a PCA plot of the regulated phosphorylation sites in the `PCA` document.

I performed a pairwise correlation analysis on the regulated phosphorylation sites. This is in the folder `CorrelationAnalysis`.