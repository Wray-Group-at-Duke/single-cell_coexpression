# single-cell_coexpression
Python and R code for counting cells that co-express two given genes and generating plots showing changes over developmental time and between species. 

## Work flow
Original inputs are scRNAseq count tables for 6 time-points from 2 species (12 in total).  
  
For 2 genes, counts from separate accessions were aggregated using **agg_gene_models.R**:  
  - *alx1* because the gene model is fragmented  
  - *pmar1* because tandem paralogues cannot be distinguished and are functionally equivalent  
  - **agg_alx1/pmar1.txt** contains the accession numbers used for aggregating  
  
Count tables with aggregated gene models are processed by **CountCoexpr.py**  
  
**CountCoexpr.py** expects to be located in a folder containing:  
  - Count tables in Separate folders: **Lv_count_tables** and **He_count_tables**  
  - The folder **Comparison_sets**, which contains folders for different gene sets  
  - The folder for each gene set contains:
    - **G-S_<gene_set>.txt**, a look-up table for human-readable gene names
    - **Lv_<gene_set>_gene_list.txt**, the Lv accession numbers of genes in the set
    - **He_<gene_set>_gene_list.txt**, the He accession numbers of genes in the set
  	 
When **CountCoexpr.py** runs, it generates:	
  - **Lv_<gene_set>_gene_pairs.txt**, a list of all permutations of Lv gene pairs  
  - **Lv_tx_ordered.csv**, a table with normalized coexpression values for all Lv gene pairs  
  - **He_<gene_set>_gene_pairs.txt**, a list of all permutations of He gene pairs  
  - **He_tx_ordered.csv**, a table with normalized coexpression values for all He gene pairs  
   
 Compression tables from Lv and He are co-analyzed with **Coexpr_plots.R**  
   - Takes as input **Lv_tx_ordered.csv** and **He_tx_ordered.csv**  
   - Generates plots of coexpression over developmental time for specified pairs of genes

testing

