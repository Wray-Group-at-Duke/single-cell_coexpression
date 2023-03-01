###############################################################
# count percentage of cells with co-expression of 2 genes
#    for use with scRNA-seq data from Lv and He projects
#    currently computes this for 1 gene set in 1 species
#    used to generate Fig # in Massri et al. 2022
# inputs: gene list, gene-symbol table, scRNA-seq count tables
# output: csv file of multistage coexpression percentages
# prior to running: 
#    generate appropriate gene_list.txt and G-S.txt files 
#    set species and gene_set values in step 0
#
# status: development
# 2022-AUG-27 Greg Wray
###############################################################

#### step 0: set-up

# import packages
import pandas as pd
import itertools as it

# set variables !!! change values as needed for each run
species = 'He'
gene_set = 'SkC'
time_points = ('6','9','12','16','20','24')

# set paths for file operations
read_path = './'+species+'_count_tables/'
write_path = './Comparison_sets/'+gene_set+'/'

# splash and confirm species and gene set
splash = 'computing % coexpression for '
print(splash+gene_set+' in '+species+' using time points',time_points)

#### step 1: construct list of gene pairs from gene set

# update progress
print('generating list of gene pairs')

# open file with list of genes in  gene set
gene_list = pd.read_csv(write_path+species+'_'+gene_set+'_gene_list.txt')

# create list of all pairwise combinations (note: generates tuples)
gene_tups = list(it.combinations(gene_list, 2))
gene_combs = [list(j) for j in gene_tups]
gene_pairs = {}
i = 0
for pair in gene_combs:
    gene_pairs[i] = pair
    i += 1
       
# write gene pair list (record-keeping, not used by program)
outfile_name = write_path+species+'_'+gene_set+'_gene_pairs.txt'
df = pd.DataFrame(gene_pairs)
df = df.T
df.to_csv(outfile_name, index = False)

#### step 2: open count tables, compute overlap for specified gene pairs 
 
# update progress
print('extracting counts and computing overlap')

# loop through count tables for each time point
curr_time_point = 0
cell_counts = list()
for t in time_points:
    curr_count_file = species+t+'hpf.csv'
 
    # read in data from one time point
    df = pd.read_csv(read_path+curr_count_file)
    print('input file read')    

    # gene names column to be df index to facilitate single row extractions
    df = df.set_index('X')

    # record number of rows (cells) for later normalization
    num_UMIs, num_cells = df.shape
    cell_counts.append(num_cells)

    # print progress info
    print('    read file: ', read_path+curr_count_file)
    print('    cell count: ', cell_counts[curr_time_point])

    # loop through pairs of genes
    out_list = []
    for i in gene_pairs:
        g1, g2 = gene_pairs[i]    

        # extract row with counts for each gene 
        gene_1 = df.loc[g1]
        gene_2 = df.loc[g2]    

        # set any count value over 1 to 1 and convert row to Series
        gene_1 = pd.Series([1 if i > 0 else 0 for i in gene_1])
        gene_2 = pd.Series([1 if j > 0 else 0 for j in gene_2])
        
        # add to create coexpression Series
        coexpr = gene_1 + gene_2
        
        # convert to boolean True if coexpression 
        coexpr = [True if k == 2 else False for k in coexpr]
        
        # tally coexpressing cells from row and normalize
        coexpr_count = coexpr.count(True)
        coexpr_perc = (coexpr_count / num_cells) * 100
        next_row = [g1, g2, coexpr_perc]
        
        # append value for gene pair to output DataFrame
        out_list.append(next_row)    
    
    # convert list to DataFrame, add column names, and write to file
    out_df = pd.DataFrame(out_list, columns = ['gene_1', 'gene_2', 'coexpr'])
    time_point_file = write_path+'Time_points/'+species+t+'coex.txt'
    out_df.to_csv(time_point_file, index = False, header = None)
    print('    wrote file:', time_point_file)
    
    # increment to next time point
    curr_time_point += 1
    
#### step 3: merge tables into time-series, rename genes, and clean up

# update progress
print('merging time points and cleaning up')

# construct list of input files to merge
input_files = []
j = 0
for k in time_points:
    input_files.append(species+time_points[j]+'coex.txt')
    j += 1

# open individual time point files and merge into DataFrame
file_num = 0
for i in input_files:
    cols = ['g1', 'g2', time_points[file_num]]
    r_path = write_path+'Time_points/'
    if file_num == 0:
        tp_0 = pd.read_csv(r_path+i, header = None, names = cols)
        coexpr_df = tp_0.set_index(['g1', 'g2'])       
    else:
        tp_n = pd.read_csv(r_path+i, header = None, names = cols)
        tp_n = tp_n.set_index(['g1', 'g2'])
        coexpr_df = pd.concat([coexpr_df, tp_n], axis = 1)
    file_num += 1  
    
# write to file and strip headers
out_file = write_path+species+'_tc_raw.csv'
coexpr_df.to_csv(out_file, header = None)
print('    wrote file: '+write_path+species+'_tc_raw.csv')

# read clean file for processing
coexpr_d2 = pd.read_csv(out_file, header = None)

# rename columns
new_col_names = {0:'gene_1', 1:'gene_2'}
entry = 0
for m in time_points:
    new_col_names[entry + 2] = time_points[entry]
    entry += 1
coexpr_d2 = coexpr_d2.rename(columns = new_col_names)
    
# read in gene name translation table
tr_df = pd.read_csv(write_path+'G-S_'+gene_set+'.txt')

# join gene names from translation table and clean up
coexpr_d2 = pd.merge(coexpr_d2, tr_df, left_on = 'gene_1', right_on = 'species_id')
coexpr_d2 = pd.merge(coexpr_d2, tr_df, left_on = 'gene_2', right_on = 'species_id')
coexpr_d2 = coexpr_d2.drop(['gene_1','gene_2','species_id_x','species_id_y'], axis=1)
new_df = coexpr_d2.rename(columns={'symbol_x':'gene_1','symbol_y':'gene_2'})


# ordinate gene names alphabetically by column (e.g., alx1, wnt8)
g1 = new_df.loc[:,'gene_1']
gc1 = g1.copy()
g2 = new_df.loc[:,'gene_2']
gc2 = g2.copy()
i = 0
while i < gc1.size:
    if gc1[i] > gc2[i]:
        gc1[i] = new_df.loc[i,'gene_2']
        gc2[i] = new_df.loc[i,'gene_1']
    i += 1    

# attach reordered gene name columns, remove old gene name columns   
new_df = new_df.rename(columns={'gene_1':'g_1','gene_2':'g_2'})
coexpr_d3 = pd.concat([new_df, gc1, gc2], axis = 1)
coexpr_d3 = coexpr_d3.drop(['g_1', 'g_2'], axis = 1)

# reorder columns
col_names = ['gene_1', 'gene_2']
tp = 0
for k in time_points:
    col_names.append(time_points[tp])
    tp += 1
coexpr_d3 = coexpr_d3[col_names]

# write out modified table to csv file
coexpr_d3.to_csv(write_path+species+'_tc_ordered.csv')
print('    wrote file: '+write_path+species+'_tc_ordered.csv')

# report completion
print('done ')




