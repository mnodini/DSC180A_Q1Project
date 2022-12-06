import pandas as pd
import io
import os
import gzip
from sklearn.linear_model import LinearRegression
from collections import defaultdict
import statsmodels.api as sm
import matplotlib.pyplot as plt
import numpy as np
def get_data(path):
    return pd.read_csv(path)

def get_gene_data(path):
    '''
    reads files into a variable
    '''
    exp = pd.read_csv(path, sep='\t')
    exp_22 = exp[exp['Chr'] == '22']
    return exp_22

def get_pop_data(path):
    '''
    reads population data into a variable
    '''
    return pd.read_csv(path, sep = ' ')
def get_vcf_data(path):
    '''
    reads vcf data
    '''
    def get_vcf_names(vcf_path):
        with gzip.open(vcf_path, "rt") as ifile:
              for line in ifile:
                if line.startswith("#CHROM"):
                      vcf_names = [x for x in line.split('\t')]
                      break
        ifile.close()
        return vcf_names
    names = get_vcf_names(path)
    return pd.read_csv(path, compression='gzip', chunksize=50000, comment='#',low_memory=False, delim_whitespace=True, header=None, names=names)

def barplot_eqtl_counts_test(barplot_csv):
    '''
    Generates barplot from data dictionary and then saves the output in a png file.
    '''
    df = pd.read_csv(barplot_csv, sep = '\t')


    pops = list(df['Unnamed: 0'].values)
    counts = list(df['eQTL Count'].values)

    #generate plot
    fig = plt.figure(figsize = (10, 5))
    # creating the bar plot
    plt.bar(pops, counts, color ='maroon',
            width = 0.4)

    plt.xlabel("Populations")
    plt.ylabel("No. of cis-eQTLs")
    plt.title("eQTL Count for Gene ENSG00000249263.2 across 3 Populations")
    plt.show()
    plt.savefig('output/barplot_eQTLs_counts_test.png')
def barplot_eqtl_counts(gbrCount, tsiCount, finCount):
    data = {'FIN' : finCount, 'GBR' : gbrCount, 'TSI' : tsiCount}
    fig = plt.figure(figsize = (10, 5))
    # creating the bar plot
    plt.bar(data.keys(), data.values(), color ='maroon',
            width = 0.4)

    plt.xlabel("Populations")
    plt.ylabel("No. of cis-eQTLs")
    plt.title("eQTL Count for Gene ENSG00000249263.2 across 3 Populations")
    plt.show()
    dataframe = pd.DataFrame.from_dict(data,orient='index')
    dataframe.to_csv('output/barplot_eqtl_counts', sep = '\t')
    
def box_whisker_plot_test(box_whisker_plot_csv):
    df = pd.read_csv(box_whisker_plot_csv)
    #creating plot
    fig,ax = plt.subplots()
    df.plot(x='populations',y='slopes',kind='scatter',yerr='std_errs',ylabel='eQTL effect size',figsize=(5,5),fontsize=10,ax=ax)
    plt.xlabel('Populations')
    plt.savefig('output/box-whisker_test.png')

def create_target(vcf):
    target = []

    for i, j in enumerate(vcf):
        target.append(j)
        print('REMOVE THIS BREAK LATER')
        break
    return target
def regression(exp_22, target):
    '''
    Runs entire pipeline to generate relevant plots and tables
    '''
    def snp_val_mod(x):
        pos_1 = int(str(x)[0])
        pos_2 = int(str(x)[2])
        return pos_1 + pos_2 
    loc = None
    for c in exp_22['Coord'].values:
            if c in target['POS'].values:
                loc = c
                break

    expression = exp_22.loc[exp_22['Coord'].isin([loc])]
    #function to narrow down relevant columns and convert datatypes 
    expression = expression.squeeze().iloc[4:].to_frame()
    allele = target.loc[target['POS'].isin([loc])]
    allele = allele.squeeze().iloc[9:]
    allele.index = allele.index.str.split().str[0]
    allele = allele.to_frame().reset_index()
    allele = allele.reset_index().rename(columns = {'index' : 'sample_id', allele.columns[1]: 'snp_value'})
    allele['snp_value'] = allele['snp_value'].apply(snp_val_mod)
    allele = allele.iloc[:,1:].set_index('sample_id')[['snp_value']]
    merged = pd.merge(expression, allele, left_index=True, right_index=True)
    merged = merged.rename(columns={merged.columns[0]: 'expression', merged.columns[1]: 'allele'})
    reg = LinearRegression().fit(merged['allele'].values.reshape(-1, 1), merged['expression'])
    print('Coefficient: ' + str(reg.coef_) + '\n Intercept: ' + str(reg.intercept_))
    print('P val: ' + str(reg.score(merged['allele'].values.reshape(-1, 1), merged['expression'])))
    
def merge(pop_data, exp_22):
    pop_data.rename(columns = {'sample' : 'sample_id'}, inplace = True)
    exp_22_transposed = exp_22.melt(list(exp_22.columns[:4]), var_name='sample_id', value_name='Value')
    pop_exp_merged = exp_22_transposed.merge(pop_data, on = 'sample_id')
    return pop_exp_merged
def match_gbr(gbr_pop,target): 
    gbr_dict = defaultdict(list)
    gbr_coords = list(gbr_pop['Coord'].unique())
    #gbr dict 
    for coord in gbr_coords:
        for j in target:
            for pos in j['POS'].values:
                if (int(pos) < (int(coord) - 50000)) and (int(pos) > (int(coord) + 50000)): 
                    continue
                else:
                    gbr_dict[coord].append(pos)
                    
    return gbr_dict
def regression_pop(population, pop, gbr_dict,target):
    #function for modifying the values
    def snp_val_mod(x):
        pos_1 = int(str(x)[0])
        pos_2 = int(str(x)[2])
        return pos_1 + pos_2
    gene_coord = 17140518
    snp_ids = []
    p_vals = []
    marker = []
    slopes = []
    positions = []
    count = 0
    for pos in gbr_dict[gene_coord]:
        #Unique issue
        if(pos in [16423023,16687501,17106293,17106294,17137315,17271797]):
            continue
        count += 1
        if count > 5000:
            break
        for chunks in target:
            allele = chunks.loc[chunks['POS'].isin([pos])]
            if len(allele) == 1:
                break
        if len(allele) == 0: 
            count -= 1
            continue
        positions.append(pos)
        #print(allele['ID'])
        snp_ids.append(allele['ID'].values[0])
        marker.append('chr22:' + str(allele['POS'].values[0]))        
        expression = pop.loc[pop['Coord'].isin([gene_coord])]

        allele = allele.squeeze().iloc[9:]
        
        if type(allele) != pd.core.frame.DataFrame: 
            allele = allele.to_frame()
        #step by step
        allele = allele.reset_index()
        allele = allele.rename(columns = {'index' : 'sample_id', allele.columns[1]: 'snp_value'})
        allele['snp_value'] = allele['snp_value'].apply(snp_val_mod)
        merged = pd.merge(expression, allele, on = 'sample_id')
        X = sm.add_constant(merged['Value'].values)
        est = sm.OLS(merged['snp_value'].values.astype(float),X.astype(float))
        est2 = est.fit()
        p_val = est2.pvalues[1]
        slope = est2.params[1]
        std_err = est2.bse[1]
        p_vals.append(p_val) #append p vals
        slopes.append(slope) # append slope
        
    locus_zoom_gbr = {'snp_ids' : snp_ids, 'p_vals' : p_vals, 'slope' : slopes, 'pos' : positions}
    locus_zoom_gbr = pd.DataFrame.from_dict(locus_zoom_gbr)
    
    if population == 'GBR': 
        locus_zoom_gbr.to_csv('output/locusTableGBR.txt', sep="\t")
    elif population == 'TSI':
        locus_zoom_gbr.to_csv('output/locusTableTSI.txt', sep="\t")
    else: 
        locus_zoom_gbr.to_csv('output/locusTableFIN.txt', sep="\t")
    significant = len(locus_zoom_gbr.loc[locus_zoom_gbr['p_vals'] <= .05/5000])
    print('# of Significant eQTLs: ' + str(significant))
    return significant
def generate_boxplots(target,pop_exp_merged):
    def snp_val_mod(x):
        pos_1 = int(str(x)[0])
        pos_2 = int(str(x)[2])
        return pos_1 + pos_2
    #Get SNP data
    allele = target.iloc[0]
    allele = allele.squeeze().iloc[9:]
    allele.index = allele.index.str.split().str[0]
    allele = allele.reset_index().rename(columns = {'index' : 'sample_id', 0: 'snp_value'})
    allele['snp_value'] = allele['snp_value'].apply(snp_val_mod)
    #create population dictionary 
    pop_dfs = dict()
    for pop in pop_exp_merged['population'].unique():
        df = pop_exp_merged.loc[pop_exp_merged['population'] == pop]
        if(df.shape[0] > 0):
            pop_dfs[pop] = pop_exp_merged.loc[pop_exp_merged['population'] == pop]
    LR_dfs = dict()
    p_vals = dict()
    std_errs = dict()
    slopes = dict()
    #then merge SNP data with population dataframes
    target_gene = 'ENSG00000249263.2'
    for key in pop_dfs.keys():
        df = pop_dfs[key][pop_dfs[key]['TargetID'] == target_gene]
        merged = pd.merge(df,allele,on='sample_id')
        LR_dfs[key] = merged
    #export plots
    #YRI
    fig, ax = plt.subplots()
    LR_dfs['YRI'].boxplot(column='Value',by='snp_value',grid=False,figsize=(5,5),ax=ax)
    ax.set_xlabel('SNP Class')
    ax.set_ylabel('Genetic Expression')
    ax.set_title('Distribution of Genetic Expression for population YRI')
    fig.suptitle('')
    plt.savefig('output/Distribution of Genetic Expression for population YRI')
    #FIN
    fig, ax = plt.subplots()
    LR_dfs['FIN'].boxplot(column='Value',by='snp_value',grid=False,figsize=(5,5),ax=ax)
    ax.set_xlabel('SNP Class')
    ax.set_ylabel('Genetic Expression')
    ax.set_title('Distribution of Genetic Expression for population FIN')
    fig.suptitle('')
    plt.savefig('output/Distribution of Genetic Expression for population FIN')
    #GBR
    fig, ax = plt.subplots()
    LR_dfs['GBR'].boxplot(column='Value',by='snp_value',grid=False,figsize=(10,10),ax=ax)
    ax.set_xlabel('SNP Class')
    ax.set_ylabel('Genetic Expression')
    ax.set_title('Distribution of Genetic Expression for population GBR')
    fig.suptitle('')
    plt.savefig('output/Distribution of Genetic Expression for population GBR')

def generate_slope_plots(target,pop_exp_merged):
    def snp_val_mod(x):
        pos_1 = int(str(x)[0])
        pos_2 = int(str(x)[2])
        return pos_1 + pos_2
    #first get SNP data for all populations
    allele = target.iloc[0]
    allele = allele.squeeze().iloc[9:]
    allele.index = allele.index.str.split().str[0]
    allele = allele.reset_index().rename(columns = {'index' : 'sample_id', 0: 'snp_value'})
    allele['snp_value'] = allele['snp_value'].apply(snp_val_mod)
    #create population dictionary 
    pop_dfs = dict()
    for pop in pop_exp_merged['population'].unique():
        df = pop_exp_merged.loc[pop_exp_merged['population'] == pop]
        if(df.shape[0] > 0):
            pop_dfs[pop] = pop_exp_merged.loc[pop_exp_merged['population'] == pop]
    p_vals = dict()
    std_errs = dict()
    slopes = dict()
    #then merge SNP data with population dataframes
    target_gene = 'ENSG00000249263.2'
    for key in pop_dfs.keys():
        df = pop_dfs[key][pop_dfs[key]['TargetID'] == target_gene]
        merged = pd.merge(df,allele,on='sample_id')
        X = sm.add_constant(merged['Value'].values)
        est = sm.OLS(merged['snp_value'].values.astype(float),X.astype(float))
        est2 = est.fit()
        p_val = est2.pvalues[1]
        std_err = est2.bse[1]
        slop = est2.params[1]
        if not np.isnan(p_val):
            p_vals[key] = p_val
            std_errs[key] = std_err
            slopes[key] = slop
    #Plotting
    plot_dict = dict()
    plot_dict['populations'] = list(p_vals.keys())
    plot_dict['slopes'] = list(slopes.values())
    plot_dict['std_errs'] = list(std_errs.values())
    fig, ax = plt.subplots()
    pd.DataFrame.from_dict(plot_dict).plot(x='populations',y='slopes',kind='scatter',yerr='std_errs',ylabel='eQTL effect size',figsize=(5,5),fontsize=10,ax=ax)
    ax.set_xlabel('Populations')
    ax.set_title('eQTL Effect Size Across Populations')
    plt.savefig('output/eQTL Effect Size Across Populations')

def test_generate_boxplots(filepath):
    LR_dfs = dict()
    LR_dfs['YRI'] = pd.read_csv(filepath + 'YRI_Boxplot.csv')
    LR_dfs['FIN'] = pd.read_csv(filepath + 'FIN_Boxplot.csv')
    LR_dfs['GBR'] = pd.read_csv(filepath + 'GBR_Boxplot.csv')
    #export plots
    #YRI
    fig, ax = plt.subplots()
    LR_dfs['YRI'].boxplot(column='Value',by='snp_value',grid=False,figsize=(5,5),ax=ax)
    ax.set_xlabel('SNP Class')
    ax.set_ylabel('Genetic Expression')
    ax.set_title('Distribution of Genetic Expression for population YRI')
    fig.suptitle('')
    plt.savefig('output/Distribution of Genetic Expression for population YRI')
    #FIN
    fig, ax = plt.subplots()
    LR_dfs['FIN'].boxplot(column='Value',by='snp_value',grid=False,figsize=(5,5),ax=ax)
    ax.set_xlabel('SNP Class')
    ax.set_ylabel('Genetic Expression')
    ax.set_title('Distribution of Genetic Expression for population FIN')
    fig.suptitle('')
    plt.savefig('output/Distribution of Genetic Expression for population FIN')
    #GBR
    fig, ax = plt.subplots()
    LR_dfs['GBR'].boxplot(column='Value',by='snp_value',grid=False,figsize=(10,10),ax=ax)
    ax.set_xlabel('SNP Class')
    ax.set_ylabel('Genetic Expression')
    ax.set_title('Distribution of Genetic Expression for population GBR')
    fig.suptitle('')
    plt.savefig('output/Distribution of Genetic Expression for population GBR')
    print('Exported Box Plots')

def test_generate_slope_plots(filepath):
    df = pd.read_csv(filepath)
    fix,ax = plt.subplots()
    df.plot(x='populations',y='slopes',kind='scatter',yerr='std_errs',ylabel='eQTL effect size',figsize=(5,5),fontsize=10,ax=ax)
    ax.set_xlabel('Populations')
    ax.set_title('eQTL Effect Size Across Populations')
    plt.savefig('output/eQTL Effect Size Across Populations')
    
    
    
    
    