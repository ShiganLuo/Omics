#! /usr/bin/env python3


import subprocess
import pandas as pd
import numpy as np
import os, logging
from collections import defaultdict
from multiprocessing import Pool

log = logging.getLogger(__name__)

def filterCoverage (cov_table, min_cov):

	log.info("\n+-----------------------------------+ \
		\n| Final deconvolution and filtering |\
		\n+-----------------------------------+")
	if min_cov < 1:
		cov_table_new = cov_table.div(cov_table.sum(axis=0), axis=1)
		filtered_list = list(cov_table_new[(cov_table_new.values < min_cov).all(1) & (~cov_table_new.index.str.contains('mito')) & (~cov_table_new.index.str.contains('plastid')) & (~cov_table_new.index.str.contains('eColi'))].index)
		log.info("{} clusters filtered out according to minimum coverage threshold: {:.2%} of total tRNA coverage.".format(len(filtered_list), min_cov))
	else:
		if min_cov == 1:
			log.warning("--min-cov set to 1: treating as integer of absolute coverage, not a fraction of mapped reads!")
		filtered_list = list(cov_table[(cov_table.values < min_cov).all(1) & (~cov_table.index.str.contains('mito')) & (~cov_table.index.str.contains('plastid')) & (~cov_table.index.str.contains('eColi'))].index)
		log.info("{} clusters filtered out according to minimum coverage threshold: {} total read coverage per isodecoder.".format(len(filtered_list), min_cov))

	filter_warning = False
	if len(filtered_list) / len(cov_table.index) >= 0.7:
		log.warning("70% or more clusters filtered out by --min-cov: consider reducing or assessing data quality and depth!")
		filter_warning = True

	return(filtered_list, filter_warning)

def getBamList (sampleGroups):

	with open(sampleGroups,"r") as sampleGroups:
		baminfo = defaultdict(list)
		bamlist = list()
		for line in sampleGroups:
			line = line.strip()
			currbam = str(line.split("\t")[0])
			condition = line.split("\t")[1]
			librarySize = int(line.split("\t")[2])
			baminfo[currbam] = [condition,librarySize]
			bamlist.append(currbam)

	return(baminfo, bamlist)

def getCoverage(sampleGroups, out_dir, control_cond, filtered_cov, unsplitCluster_lookup):

	log.info("\n+-----------------------------------+\
		\n| Calculating coverage and plotting |\
		\n+-----------------------------------+")

	baminfo = getBamList(sampleGroups)[0]
	cov_mean = list()

	for bam, info in baminfo.items():

		coverage = pd.read_csv(os.path.dirname(bam) + "/" + bam.split("/")[-1] + "_coverage.txt", sep = "	") 
		coverage = coverage[~coverage.isodecoder.isin(filtered_cov)]
		coverage.loc[~coverage['isodecoder'].str.contains("chr"), 'isodecoder'] = coverage['isodecoder'].str.split("-").str[:-1].str.join("-")
		coverage.set_index("isodecoder", inplace=True)
		coverage['aa'] = coverage.index.format()
		coverage.loc[coverage.aa.str.contains('chr'), 'aa'] = coverage[coverage.aa.str.contains('chr')].aa.str.split("-").str[-4]
		coverage.loc[coverage.aa.str.contains('mito'), 'aa'] = "mito" + coverage[coverage.aa.str.contains('mito')].aa.str.split("-").str[-3]
		coverage.loc[coverage.aa.str.contains('plastid'), 'aa'] = "plastid" + coverage[coverage.aa.str.contains('plastid')].aa.str.split("-").str[-3]
		coverage.loc[coverage.aa.str.contains('nmt') & ~coverage.aa.str.contains('chr'), 'aa'] = "nmt" + coverage[coverage.aa.str.contains('nmt') & ~coverage.aa.str.contains('chr')].aa.str.split("-").str[-3]
		coverage.loc[~coverage.aa.str.contains('mito') & ~coverage.aa.str.contains('plastid') & ~coverage.aa.str.contains('nmt') & ~coverage.aa.str.contains('chr'), 'aa'] = coverage[~coverage.aa.str.contains('mito') & ~coverage.aa.str.contains('plastid') & ~coverage.aa.str.contains('nmt') & ~coverage.aa.str.contains('chr')].aa.str.split("-").str[-3]
		coverage = coverage[['pos','cov','aa','bam']]
		coverage['condition'] = info[0]
		coverage['cov'] = coverage['cov'].astype(float)
		coverage['cov_norm'] = coverage['cov'] / info[1]

		coverage['bin'] = coverage.groupby([coverage.index])['pos'].transform(lambda x: pd.qcut(x, 25, labels=range(4,104,4)))
		cov_mean.append(coverage)
		os.remove(os.path.dirname(bam) + "/" + bam.split("/")[-1] + "_coverage.txt")	 

	cov_mean = pd.concat(cov_mean, axis = 0)
	cov_mean_gene = cov_mean.copy()
	cov_mean_gene['Cluster'] = cov_mean_gene.index.format()
	cov_mean_gene = cov_mean_gene[['Cluster','pos','cov','aa','condition','cov_norm','bin','bam']]
	cov_mean_gene = cov_mean_gene.groupby(['Cluster', 'bin', 'condition', 'bam']).mean()
	cov_mean_gene = cov_mean_gene.dropna()
	cov_mean_gene.reset_index(inplace = True)
	cov_mean_gene = cov_mean_gene[['Cluster','bin','condition','bam','pos','cov','cov_norm']]
	cov_mean_gene['Cluster'] = cov_mean_gene['Cluster'].replace(unsplitCluster_lookup)
	cov_mean_gene.to_csv(out_dir + "coverage_bygene.txt", sep = "\t", index = False)

	cov_mean_aa = cov_mean.groupby(['aa', 'condition', 'bam', 'pos']).sum()
	cov_mean_aa = cov_mean_aa.reset_index()
	cov_mean_aa = cov_mean_aa.groupby('aa').apply(lambda group: group.loc[group['pos'] != group['pos'].max()])
	cov_mean_aa = cov_mean_aa.drop(columns='aa')
	cov_mean_aa = cov_mean_aa.reset_index()
	cov_mean_aa = cov_mean_aa.drop(columns='level_1')
	cov_mean_aa['bin'] = cov_mean_aa.groupby(['aa','condition','bam'])['pos'].transform(lambda x: pd.qcut(x, 25, labels=range(4,104,4)))
	cov_mean_aa = cov_mean_aa.groupby(['aa', 'bin', 'condition', 'bam']).mean()
	cov_mean_aa = cov_mean_aa.dropna()
	cov_mean_aa	= cov_mean_aa.reset_index()
	cov_mean_aa.to_csv(out_dir + "coverage_byaa.txt", sep = "\t", index = False)
	
	if control_cond:
		cov_mean_aa_controlcond = cov_mean_aa[cov_mean_aa.condition == control_cond]
	else:
		first_cond = cov_mean_aa['condition'].unique()[0]
		cov_mean_aa_controlcond = cov_mean_aa[cov_mean_aa.condition == first_cond]
	bam = pd.unique(cov_mean_aa_controlcond['bam'])[0] 
	cov_mean_aa_controlcond = cov_mean_aa_controlcond[cov_mean_aa_controlcond.bam == bam]
	cov_ratios = dict()
	for aa, data in cov_mean_aa_controlcond.groupby('aa'):
		try:
			ratio = float(data[data.bin == 8]['cov_norm']) / float(data[data.bin == 92]['cov_norm'])
		except ZeroDivisionError:
			ratio = float(data[data.bin == 92]['cov_norm'])
		except TypeError:
			ratio = float(data[data.bin.astype(float) == data.bin.astype(float).nlargest(2).iloc[1]]['cov_norm'])
		cov_ratios[aa] = ratio
	sorted_aa = sorted(cov_ratios, key = cov_ratios.get)
	sorted_aa = "_".join(str(e) for e in sorted_aa)

	return(sorted_aa)

def plotCoverage(out_dir, mito_trnas, sorted_aa):
	
	script_path = os.path.dirname(os.path.realpath(__file__))
	command = ["Rscript", script_path + "/coveragePlot.R", out_dir + "coverage_bygene.txt", out_dir + "coverage_byaa.txt", out_dir, sorted_aa, mito_trnas]
	try:
		process = subprocess.Popen(command, stdout = subprocess.PIPE)
		while True:
			line = process.stdout.readline()
			if not line:
				break
			line = line.decode("utf-8")
			log.info(line.rstrip())
		exitcode = process.wait()
	except Exception as e:
		logging.error("Error in {}".format(command), exc_info=e)
		raise