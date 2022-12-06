#!/usr/bin/env python

import sys
import os
import io
import json
import gzip
import pandas as pd
from etl import get_data
from etl import get_gene_data
from etl import get_vcf_data
from etl import barplot_eqtl_counts_test
from etl import barplot_eqtl_counts
from etl import box_whisker_plot_test
from etl import create_target
from etl import regression
from etl import merge
from etl import match_gbr
from etl import regression_pop
from etl import generate_boxplots
from etl import get_pop_data
from etl import generate_slope_plots
from etl import test_generate_boxplots
from etl import test_generate_slope_plots

def main(targets):
    if 'test' in targets:
        filepath = 'data/test_data/'
        print('Plotting barplot')
        barplot_eqtl_counts_test(filepath + 'barplot_eqtl_counts.csv')
        print('Plotting boxplots')
        test_generate_boxplots(filepath)
        print('Plotting eQTL effect size plot')
        test_generate_slope_plots(filepath + 'slope_data.csv')
        
    elif 'raw' in targets:
        #step one: read in csvs
        filepath = 'data/raw_data/'
        print('Reading in CSVs')
        gene_exp = get_gene_data(filepath + 'GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz')
        pop_data = get_pop_data(filepath + 'ALL_1000G_phase1integrated_v3.sample')
        vcf = get_vcf_data(filepath + 'ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz')
        #populate target
        target = create_target(vcf)                
        #step two: regression on one snp-gene pair
        print('Regression on One SNP-Gene Pair')
        regression(gene_exp, target[0])
        print('Merging population and gene expression datasets')
        #step three: merge population and gene expression datasets
        
        pop_exp_merged = merge(pop_data, gene_exp)
        
        print('Creating population specific subsets')
        #step four: create a separate population specific dataset
        tsi_pop = pop_exp_merged.loc[pop_exp_merged['population'] == 'TSI']
        gbr_pop = pop_exp_merged.loc[pop_exp_merged['population'] == 'GBR']
        fin_pop = pop_exp_merged.loc[pop_exp_merged['population'] == 'FIN']
        print('Populating dictionary with SNP-Gene matches')
        #step five: populate dictionary with snp-gene matches:
        gbr_dict = match_gbr(gbr_pop,target)
        print('Running regression for one gene in 3 populations')
        #step six: run a regression for one gene in three populations and output tables 
        
        gbr_eqtl_counts = regression_pop('GBR', gbr_pop,gbr_dict,target)
        tsi_eqtl_counts = regression_pop('TSI', gbr_pop,gbr_dict,target)
        fin_eqtl_counts = regression_pop('FIN', gbr_pop,gbr_dict,target)       
        print('Plotting Barplot')
        #create barplot from metrics from above
        barplot_eqtl_counts(gbr_eqtl_counts, tsi_eqtl_counts, fin_eqtl_counts)
        
        print('Plotting boxplots')
        #create boxplots
        generate_boxplots(target[0], pop_exp_merged)
        print('Plotting slope plots')
        #create slope plots
        generate_slope_plots(target[0],pop_exp_merged)

    return

if __name__ == '__main__':
    targets = sys.argv[1:]
    main(targets)
