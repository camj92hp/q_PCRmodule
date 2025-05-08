
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import ttest_ind

import matplotlib.pyplot as plt
import os


#Multiple data frames simultaneaous cleaning

def load_file(file):

    if os.path.exists(file):
        df = pd.read_csv(file, index_col=0)
        print("df loaded successfully ✅")
        return df
    else:
        print(f"File not found: {file} ❌")
    
    

    




def data_cleaning(sample_name_list, *args):
    #sample_name_list = the list that describe the treatment you have in your cells. For instance if your treatment is the mutation, then you should call it Non-mutant as control
    #mutant as treatment 
    
    #concat the data frame
    df = pd.concat(args)

    if df.empty:
        return pd.DataFrame()

    #remove some useless columns, 
    # If ther are any columns that you will not use in your analysis, or that are just empty, you can remove them here,)
    # Just add the column name on the list below
    columns_to_drop = ['Starting Quantity (SQ)', 'Log Starting Quantity', 'SQ Mean', 'Well Note', 'Cq Std. Dev', 'Set Point', 'Cq Mean', 'SQ Std. Dev']
    
    columns_to_drop = [col for col in columns_to_drop if col in df.columns]

    df = df.drop(columns_to_drop, axis=1)


    #check that the timepoints are proper

    #proper_name = ['0 HR', '0.5 HR', '1 HR', '2 HRS', '4 HRS', '6 HRS']
    proper_name = sample_name_list

    ####################################################################################################################################################
    # OPTIONAL STEP
    #If you know what are the usual mistakes you make when naming your samples you can run the line below,
    # I created a dictionary that would keep the common mistakes I make and replace them with the actual sample list item.  

    df.loc[:, 'Sample'].replace({'CTL': '0 HR', '30MIN LPS': '0.5 HR', '2HR': '2 HRS', '4HR': '4 HRS', '6HR': '6 HRS', 'NT CTL WT': '0 HR', 
                         'NT CTL STING KI': '0 HR', '4  HRS': "4 HRS", '1 HR LPS': "1 HR",  '2HR LPS' :"2 HRS",  '6HRS' : '6 HRS',
                         '4HR LPS': "4 HRS", '1HR': "1 HR", 'O HR': '0 HR', '30 MIN': '0.5 HR',
                          '6HR LPS' :"6 HRS", '0H':'0 HR' , '0HRS': '0 HR'},  inplace=True) 
    
    
    ######################################################################################################################################################
    # Remove all NA

    df= df.dropna(subset=['Sample', 'Cq'])

    # Check if you have the right number of sample groups, so that  you have an idea if your data needs to be double check
    samples = df['Sample'].unique()
    counter = 0
    for i in samples:

        if i in proper_name:
            counter +=1
    
    if counter == len(proper_name):
        return df

    else:
        print("Samples need to be double checked manually")
        return df



#This is the same as above except that this is for a single df. 
def data_cleaning_single(df, drop, sample_name_list):
    
    #concat the data frame
    #df = pd.concat(args)

    if df.empty:
        return pd.DataFrame()

    #remove some useless columns
    df= df.drop(drop, axis = 1 )

    

    #check that the timepoints are proper

    #proper_name = ['0 HR', '0.5 HR', '1 HR', '2 HRS', '4 HRS', '6 HRS']
    proper_name = sample_name_list


    df.loc[:, 'Sample'].replace({'CTL': '0 HR', '30MIN LPS': '0.5 HR', '2HR': '2 HRS', '4HR': '4 HRS', '6HR': '6 HRS', 'NT CTL WT': '0 HR', 
                         'NT CTL STING KI': '0 HR', '4  HRS': "4 HRS", '1 HR LPS': "1 HR",  '2HR LPS' :"2 HRS",  '6HRS' : '6 HRS',
                         '4HR LPS': "4 HRS", '1HR': "1 HR", 'O HR': '0 HR', '30 MIN': '0.5 HR',
                           '6HR LPS' :"6 HRS", '0H':'0 HR', '0HRS': '0 HR'  }, inplace=True) 
    
    
    
    df= df.dropna(subset=['Sample', 'Cq'])

    samples = df['Sample'].unique()
    counter = 0
    for i in samples:

        if i in proper_name:
            counter +=1
    
    if counter == len(proper_name):
        return df

    else:
        print("Samples need to be double checked manually")
        return df






        #the first thing i need to do is calculate delta_ct. This delata city is only for the wt control at time 0
# so is the average of the Gene of interest (GOI) at the control time for the wt  - average of ref gene  for the WT with no treatment 
 

def delta_ct(df, gene, cell,  ref = 'GAPDH', sample = '0 HR'):
    
    ref_df = df[(df['Target'] == ref) & (df['Sample']== sample) & (df['Biological Set Name']== cell)]
    
    gene_df = df[(df['Target'] == gene) & (df['Sample']== sample) & (df['Biological Set Name']== cell)]
    
    ave_ref = ref_df['Cq'].mean()

    ave_gene = gene_df['Cq'].mean()

    del_ct = ave_gene- ave_ref

    
    return del_ct




    #now we know that we have to go gene by gene to check the delta delta ct. IF WE ARE INTERESTED IN TREATMENT EFFECT WE MUST DO THIS ONE
# HERE WE JUST CALCULATE WHAT IS THE EFFECT OF your treatment  ON EVERY CELL LINE

def delta_delta_ct_treatment_effect(df, gene, cell, ref_cell, ref = 'GAPDH',):

    ref_df = df[(df['Target'] == ref) & (df['Biological Set Name']== cell)] 
    gene_df = df[(df['Target'] == gene) & (df['Biological Set Name']== cell)]

    group_ref = ref_df.groupby(['Sample'])
    
    mean_ref_ct = group_ref['Cq'].mean().reset_index()


    merged_df = pd.merge(gene_df, mean_ref_ct, on=['Sample'], suffixes=('_gene_of_interest', '_ref_gene'))
    
    #we are making a data frame where we can store our answers
    res_df = pd.DataFrame()
    
    res_df[['Sample']] = merged_df[['Sample']]
    res_df['Biological Set Name'] = cell
    res_df [['Target', 'Cq_gene_of_interest', 'Cq_ref_gene', ]] = merged_df [['Target', 'Cq_gene_of_interest', 'Cq_ref_gene', ]]
    res_df['delta_ct'] = delta_ct(df, gene, ref_cell, ref)
    
    

    res_df["delta delta ct"] = (merged_df['Cq_gene_of_interest']- merged_df['Cq_ref_gene'])- res_df['delta_ct']
    res_df['log_change'] = 2**-(res_df['delta delta ct'])
    res_df['Analysis type'] = 'treatment_effect'

    return res_df




    #now we know that we have to go gene by gene to check the delta delta ct. IF WE ARE INTERESTED IN gene EFFECT WE MUST DO THIS ONE
# HERE WE JUST CALCULATE WHAT IS THE EFFECT OF the gene mutation ON THE 6 HRS ON EVERY CELL LINE

#ref_cell =  tells you what cell line and time point use as reference

def delta_delta_ct_gene_effect(df, gene, cell, ref_cell = 'ATCC', treatment='6 HRS', ref="GAPDH"):

    ref_df = df[(df['Sample'] == treatment) & (df['Biological Set Name']== ref_cell) & (df['Target'] == ref) ] 
    gene_df = df[(df['Target'] == gene) & (df['Biological Set Name']== cell) &(df['Sample'] == treatment)]
    # print("ref sample df")
    print(ref_df)
    # print("target cell sample df")
    # print(gene_df)
   

    group_ref = ref_df.groupby(['Sample'])
    # print("group df")
    # print(group_ref)
   
    
    mean_ref_ct = group_ref['Cq'].mean().reset_index() #take the average of the ct of my ref gene
    print("mean_ref_ct")
    print(mean_ref_ct)


    merged_df = pd.merge(gene_df, mean_ref_ct, on=['Sample'], suffixes=('_gene_of_interest', '_ref_gene'))

    #print(merged_df)
    
    #we are making a data frame where we can store our answers
    res_df = pd.DataFrame()
    
    res_df[['Sample']] = merged_df[['Sample']]
    res_df['Biological Set Name'] = cell
    res_df [['Target', 'Cq_gene_of_interest', 'Cq_ref_gene', ]] = merged_df [['Target', 'Cq_gene_of_interest', 'Cq_ref_gene', ]]
    res_df['delta_ct_ref_cells'] = delta_ct(df, gene,ref_cell, ref, treatment)
    res_df['delta_ct_target_cells'] = merged_df['Cq_gene_of_interest']- merged_df['Cq_ref_gene']
    #print(delta_ct(df, gene, ref_cell,  ref, treatment))
    

    res_df["delta delta ct"] =res_df['delta_ct_target_cells'] - res_df['delta_ct_ref_cells']
    res_df['log_change'] = 2**-(res_df['delta delta ct'])
    res_df['Analysis type'] = 'gene_effect'
    

    return res_df



def delta_delta_ct_ctl_baseline(df, gene, cell, ref_cell = 'ATCC', ref_treatment= '0 HR', treatment='6 HRS', ref="GAPDH"):


    #this is to do the delta-ct of gene of intereste so i need gene of interest and house keeping gene
    ref_df = df[(df['Sample'] == treatment) & (df['Biological Set Name']== cell) & (df['Target'] == ref) ] 
    gene_df = df[(df['Target'] == gene) & (df['Biological Set Name']== cell) &(df['Sample'] == treatment)]
    print("ref sample df")
    print(ref_df)
    print("target cell sample df")
    print(gene_df)
   

    group_ref = ref_df.groupby(['Sample'])
   
    
    mean_ref_ct = group_ref['Cq'].mean().reset_index() #take the average of the ct of my ref gene
    # print("mean_ref_ct")
    # print(mean_ref_ct)


    merged_df = pd.merge(gene_df, mean_ref_ct, on=['Sample'], suffixes=('_gene_of_interest', '_ref_gene'))

    #print(merged_df)
    
    #we are making a data frame where we can store our answers
    res_df = pd.DataFrame()
    
    res_df[['Sample']] = merged_df[['Sample']]
    res_df['Biological Set Name'] = cell
    res_df [['Target', 'Cq_gene_of_interest', 'Cq_ref_gene', ]] = merged_df [['Target', 'Cq_gene_of_interest', 'Cq_ref_gene', ]]
    res_df['delta_ct_ref_cells'] = delta_ct(df, gene, ref_cell, ref, ref_treatment) ## this is for the control gene I just need house keeping
    
    res_df['delta_ct_target_cells'] = merged_df['Cq_gene_of_interest']- merged_df['Cq_ref_gene']
    print(delta_ct(df, gene, ref_cell,  ref, treatment))
    

    res_df["delta delta ct"] =res_df['delta_ct_target_cells'] - res_df['delta_ct_ref_cells']
    res_df['log_change'] = 2**-(res_df['delta delta ct'])
    res_df['Analysis type'] = 'Baseline is 0hr ctl attc'
    

    return res_df


def all_deltas_treatment_effect(df, ref = 'GAPDH'):

    #now we can use the delta make to do all the genes at onece

    genes = df['Target'].unique()
    cell_type = df['Biological Set Name'].unique()
    
    df_res = pd.DataFrame()
    
    for g in genes:

        for c in cell_type:


            if g != ref:
                data = delta_delta_ct_treatment_effect(df, g, c, ref )
                
                df_res = pd.concat([df_res, data])
    
    return df_res




def all_deltas_gene_effect(df, ref_cells= 'ATTC', ref = 'GAPDH'):

    #now we can use the delta make to do all the genes at onece

    genes = df['Target'].unique()
    cell_type = df['Biological Set Name'].unique()
    sample_type = df['Sample'].unique()


    df_res = pd.DataFrame()

    for g in genes:

        for c in cell_type:

            for s in sample_type:


                if (g != ref) :
                    data = delta_delta_ct_gene_effect(df, g, c, ref_cells, s, ref)
                    
                    df_res = pd.concat([df_res, data])

    return df_res

def all_deltas_ctl_baseline(df, ref_cells='ATTC', ref_treatment='0 HR', ref='GAPDH', verbose=False):
    """
    Computes ΔΔCt for all genes, cell types, and sample conditions,
    using ΔCt of the gene at 0 HR as the reference.
    """
    genes = df['Target'].unique()
    cell_types = df['Biological Set Name'].unique()
    sample_types = df['Sample'].unique()

    df_res = pd.DataFrame()
    counter=0
    for g in genes:
        if g == ref:  # Skip housekeeping gene
            continue

        for c in cell_types:
            for s in sample_types:
                counter+=1
                if verbose:
                    print(f"Processing Gene: {g}, Cell: {c}, Sample: {s},")
                    print(counter)
                
                # Ensure that ΔCt at 0 HR is always used as reference
                data = delta_delta_ct_ctl_baseline(df, gene=g, cell=c, ref_cell=ref_cells, 
                                                   ref_treatment=ref_treatment, treatment=s, ref=ref)
                
                if data is not None and not data.empty:
                    df_res = pd.concat([df_res, data], ignore_index=True)

    return df_res




def graph_qpcr(df):

    # Get a list of unique genes
    genes = df['Target'].unique()

    # Create a figure and subplots
    fig, axes = plt.subplots(len(genes), 1, figsize=(6, 30))

    # Ensure `axes` is always iterable
    if len(genes) == 1:
        axes = [axes]  # Wrap in a list to avoid indexing errors

    # Iterate over each gene and plot the qPCR data
    for i, gene in enumerate(genes):
        gene_data = df[df['Target'] == gene]  # Filter data for the gene
        ax = axes[i]  # Safe indexing

        sns.barplot(data=gene_data, x='Sample', y='log_change', hue='Biological Set Name', ax=ax, errorbar='sd')
        ax.set_title(gene)
        ax.set_xlabel('Sample')
        ax.set_ylabel('Log Change')

    plt.tight_layout()
    plt.show()






def transposed_values(df):

    genes = df['Target'].unique()
    timepoints = df['Sample'].unique()
    bio = df['Biological Set Name'].unique()

    #first we are gonna create dictionary with the new data structure we want to store the group valiew we want to transpose
    dic_final = dict()
    for g in genes:
        for t in timepoints:
            for b in bio:

                values = df[(df['Target']==g) & (df['Sample'] ==t) & (df['Biological Set Name']==b)]['log_change'].values

                if 'Target' in dic_final:
                    
                    dic_final['Target'].extend([g])
                
                if 'Sample' in dic_final:
                    dic_final['Sample'].extend([t])
                
                if 'Biological Set Name' in dic_final:
                    dic_final['Biological Set Name'].extend([b])
                
                if 'log_change' in dic_final:
                    dic_final['log_change'].extend([values])
                
                if 'values_lenght' in dic_final:
                    dic_final['values_lenght'].extend([len(values)])
                
                else:
                    dic_final['Target'] = [g]
                    dic_final['Sample']= [t]
                    dic_final['Biological Set Name'] = [b]
                    dic_final['log_change']=[values]
                    dic_final['values_lenght'] = [len(values)]

    #Now that we have a dictionary we have to create the DataFrame that will contain the rows to transpose. 
    T_df = pd.DataFrame(dic_final)

    def split_array(row):
        print(row)
        row = list(row) if isinstance(row, np.ndarray) else row  # Convert NumPy array to list
        length = len(row)
        if length >= 3:
            return row[:3]
        elif length == 2:
            print("only 2 values her")
            return row + [0]
        else:
            print("only 3 values her")
            return row + [0, 0]
    
    T_df[['log_change_1', 'log_change_2', 'log_change_3']] = T_df['log_change'].apply(split_array).apply(pd.Series)


    return T_df




def saveToExcel(df, Filename):
    df.to_excel(Filename, sheet_name='Sheet1', index=False)



