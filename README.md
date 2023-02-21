# DEGGs
## Differentially Expressed Gene-Gene pairs
The DEGGs package test for differential gene-gene correlations across different groups of samples in count data from high-throughput 
sequencing assays.  
Specific gene-gene interactions can be explored and gene-gene pair regression plots can be interactively shown.   

### Installation instructions 
To install from Github please use the following on your R console  
`devtools::install_github("elisabettasciacca/DEGGs", build_vignettes = TRUE)`

### Example  
Load package and sample data   
`library(DEGGs)  
data("BRCA_metadata")  
data("BRCA_normCounts")`  

Generate specific gene-gene networks for each subtype  
`subnetworks_object <- generate_subnetworks(normalised_counts = BRCA_normCounts,  
                                           metadata = BRCA_metadata,  
                                           subgroup_variable = "SUBTYPE",  
                                           subgroups = c("BRCA_Her2",  
                                                         "BRCA_LumA"),  
                                           entrezIDs = TRUE,  
                                           convert_to_gene_symbols = TRUE,  
                                           cores = 2)` 
  
Visualise  
`View_interactive_subnetwork(subnetworks_object)`  
  
Get a table listing all the significant gene-gene interactions found in each subtype  
`extract_sig_deggs(subnetworks_object)`  
   
Print differential regression fits for a single gene-gene interaction through the `print_regressions` function  
`print_regressions(gene_A = "NOTCH2", gene_B = "DTX4",
                  deggs_object = subnetworks_object,
                  legend_position = "bottomright")`
