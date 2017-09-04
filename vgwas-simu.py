#!/usr/bin/env python
##
##

import os
import sys
import getopt
import re
import random
import math
from scipy.stats import norm
import numpy as np
from parser import *
from phenotyper import *
from converter import *


class CommandLine:
	def __init__(self):
		self.file = None
		self.outfile= None
		self.num_qtn = 1
		self.input = 'G'
		self.het = 0
		self.qual = 0
		self.prev = 0.01
		self.snpfile = None
		self.maf_range_causal = [0.05,0.45]
		self.num_genotypes = None
		self.dominant=0
		self.tot_var=1.0
		self.base_avg=0.0
		self.sample=0
		self.sample_size=[1,1]
		self.min = float("inf")
		self.max = float("inf")
		self.ld = 0
		self.ldfile = None

		
class CommandlineError(Exception): 
    def __init__(self): 
        pass		


def usage():
	usage='''
	
vgwas-simu is a tool for simulating mean and variance effects in phenotype value calculation V1.0

-i				type of input file ("G" for GENOME, "M" for ms) (default: G)
-h				a binary value: either homozygous (0) or heterozygous (1) (default: 0)
-q				quantitative phenotypes (0) or qualitative case/control phenotypes (1)
--file			name of input file (In the format of ms or GENOME)
--snpfile		file with tab delimited effect sizes of SNPs on total variance
--ldfile		file with tab delimited linkage disequilibrium in square format as outputted by PLINK
--outfile		prefix for the output files in PLINK format (default: name of input file)
--dominance		co-dominance (0) or complete dominant (1) model used (default: 0)
--maf_r			MAF range for causal markers if SNP not specified in snpfile (upper and lower bound, separated by a comma, no space) (default: 0.05,0.45)
--base_avg		baseline average for phenotype (default: 0.0)
--tot_var		total variance for phenotype (default: 1.0)
--min			minimum possible value for phenotype (default: -inf)
--max			maximum possible value for phenotype (default: +inf)
--prev			disease prevalence in case/control studies (default: 0.01)
--sample		sample case/control outputs (control then case size, separated by a comma e.g 500,500)
'''
	sys.exit(usage)


def parse_commandline(commandline):
	"""parse commandline and return instance of commandline object"""
	try:
		short_options="i:h:q:"
		long_options="file=:snpfile=:ldfile=:outfile=:dominance=:maf_r=:base_avg=:tot_var=:prev=:sample=:min=:max="
		long_options=long_options.split(':')
		opts, args = getopt.getopt(commandline, short_options,long_options)
	except getopt.GetoptError:
		sys.stderr.write('Wrong options\n')
		sys.exit(2)

	#print opts,args

	cl = CommandLine()
	for option, argument in opts:
		try:
			option = re.sub('\-+','',option)
			if option in ['file']:
				cl.file = argument
			elif option in ['i']:
				cl.input = argument
				if cl.input not in ['M','G']:
					raise ValueError,'wrong input option'
			elif option in ['snpfile']:
				cl.snpfile = argument
			elif option in ['outfile']:
				cl.outfile = argument
			elif option in ['ldfile']:
				cl.ld = 1
				cl.ldfile = argument
			elif option in ['maf_r']:
				cl.maf_range_causal = map(float,argument.split(','))
				if len(cl.maf_range_causal)!=2:
					raise ValueError, 'invalid range'
			elif option in ['h']:
				cl.het = int(argument)
				if cl.het not in [1,0]:
					raise ValueError,'wrong option'
			elif option in ['q']:
				cl.qual = int(argument)
				if cl.qual not in [1,0]:
					raise ValueError,'wrong option'
			elif option in ['prev']:
				cl.prev = float(argument)
				if(cl.prev<0 or cl.prev>1):
					raise ValueError,'wrong option'
			elif option in ['sample']:
				cl.sample = 1
				cl.sample_size = map(float,argument.split(','))
				if len(cl.sample_size)!=2:
					raise ValueError, 'invalid range'
			elif option in ['dominance']:
				cl.dominant = int(argument)
				if cl.dominant not in [1,0]:
					raise ValueError,'wrong option'
			elif option in ['base_avg']:
				cl.base_avg = float(argument)
			elif option in ['tot_var']:
				cl.tot_var = float(argument)
				if cl.tot_var<0:
					raise ValueError,'wrong option'
			elif option in ['min']:
				cl.min = float(argument)
			elif option in ['max']:
				cl.max = float(argument)
			# include more code here to parse other options
		except ValueError, e:
			err_str=' '.join(['Invalid option:','-'+option,argument])
			raise ValueError, err_str
	return cl



def binary_search(liste,eingabe):
	maxindex = len(liste) - 1
	suche_erfolgreich = False
	index = 0
	while not suche_erfolgreich and index <= maxindex and index >=0:
		mitte = index + ((maxindex - index)/2)
	
		if mitte==maxindex:
			if liste[mitte]>=eingabe:
				return mitte
		if liste[mitte]<=eingabe<liste[mitte+1]:
			return mitte+1
		if mitte==0 and liste[0]>=eingabe:
			return 0
		elif eingabe < liste[mitte]:
			maxindex = mitte - 1
		else:
			index = mitte + 1


def get_phen(val,threshold):
	if(val > threshold): 
		return 2 
	else: 
		return 1


def run(cl):

	print 'Running vgwas-simu...'


	if cl.sample and not cl.qual:
		raise ValueError, 'Cannot use sampling if not qualitative phenotype'

	if cl.dominant and not cl.het:
		raise ValueError, 'Cannot use dominant model if not heterozygous'
		
	if not cl.file:
		raise ValueError, 'No input file (-f) defined'

	if not cl.snpfile:
		raise ValueError, 'No SNP file (--snpfile) defined'

	if cl.min != float("inf") and cl.max != float("inf") :
		if cl.min > cl.max:
			raise ValueError, 'Minimum cannot be greater than maximum'
			

	#parse input
	if cl.input=='M':
		genotypes_all,positions_all=parse_ms(cl.file,cl.het)
	elif cl.input=='G':
		genotypes_all,positions_all=parse_genome(cl.file,cl.het)
	else:
		raise ValueError



	causal_snps = parse_snps(cl.snpfile)
	print(causal_snps)

	if cl.ld:
		ld_map = parse_ld(cl.ldfile, len(causal_snps))
	else: 
		ld_map = None

	print 'Number of Simulations: ',len(genotypes_all)
	print 'Number of Causal SNPs: ',len(causal_snps)
		
	tot_eff = 0.0
	for snp in causal_snps:
# 		print snp		
		tot_eff = tot_eff + float(snp[1]) + float(snp[2])

	if( tot_eff > 1.0):
		raise ValueError, 'Total effect sizes cannot be greater than 1.'
		
	#now per simulation
	for sim in xrange(len(genotypes_all)):
		geno=genotypes_all[sim]
		pos=positions_all[sim]


		#create phenotypes

		indices,phen,mafs, alphas, phis = assign_phenotype_quant(geno, causal_snps, maf_range=cl.maf_range_causal, het=cl.het, dominant=cl.dominant, tot_var=cl.tot_var, base_avg=cl.base_avg, min = cl.min, max = cl.max, ld=cl.ld, ld_map = ld_map) 		
		

		#create case/control values
		if(cl.qual):
			new_mean = cl.base_avg + np.dot(mafs,alphas)
			new_var = cl.tot_var + np.dot(mafs,phis)
			threshold = norm.ppf(1-cl.prev, loc=new_mean, scale=math.sqrt(new_var))
			qual_phen = [ get_phen(i,threshold) for i in phen]

			
			if(cl.sample):
				qual_phen0 = np.array(np.where(np.array(qual_phen)==1))[0]
				qual_phen1 = np.array(np.where(np.array(qual_phen)==2))[0]
				print(len(qual_phen0), len(qual_phen1))
				if(len(qual_phen0) < cl.sample_size[0] or len(qual_phen1) < cl.sample_size[1]):
					raise ValueError, 'Not enough subjects for sampling'
			
				sample0=list(np.random.choice(qual_phen0, size = cl.sample_size[0]))
				sample1=list(np.random.choice(qual_phen1, size = cl.sample_size[1]))
				sample = sample0+sample1



		#write output	
		if not cl.outfile:
			cl.outfile=cl.file

		fname='%s%d' %(cl.outfile,sim)
		
		if(cl.qual):
			if(cl.sample):
				geno_sample = [geno[i] for i in sample]
				phen_sample = [qual_phen[i] for i in sample]
				convert2plink(geno_sample,pos,phen_sample,fname,het=cl.het)
			else:
				convert2plink(geno,pos,qual_phen,fname,het=cl.het)			
		else:
			convert2plink(geno,pos,phen,fname,het=cl.het)
				
		#write causal_file
		f_causal=open('%s.causal' %fname,'w')
		for j in xrange(len(mafs)):
			snp = causal_snps[j]
			l='%d\t%f\t%s\t%s\n' %(indices[j],mafs[j],snp[1], snp[2])
			f_causal.write(l)

		f_causal.close()


def main():
	if len(sys.argv)==1:
		usage()
	cl = parse_commandline(sys.argv[1:])
	run(cl)


if __name__ == '__main__':
	main()
