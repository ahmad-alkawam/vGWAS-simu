import random
import math
import datetime
import numpy as np


def calc_phi_cov(phi,alleles, ld_matrix):
	nsnps = len(phi)
	
	cov_val = 0.0
	for i in range(1, nsnps):
		for j in range(i+1, nsnps):
			if(alleles[i] != 0 and alleles[j] != 0):
				cov_val =  cov_val+0.5*(float(alleles[i])*float(phi[i]) + float(alleles[j])*float(phi[j]))*float(ld_matrix[i][j])
	return cov_val
			
def calc_alpha_cov(alpha,alleles, ld_matrix):
	nsnps = len(alpha)
	
	cov_val = 0.0
	for i in range(1, nsnps):
		for j in range(i+1, nsnps):
			if(alleles[i] != 0 and alleles[j] != 0):
				cov_val =  cov_val+0.5*(float(alleles[i])*float(alpha[i]) + float(alleles[j])*float(alpha[j]))*float(ld_matrix[i][j])
	return cov_val




def assign_phenotype_quant(marker, causal_snps, maf_range, het, dominant, base_avg, tot_var, min, max, ld, ld_map):


	mafs=[]

# 	calculate minor allele frequencies of genotypes:
	for index in xrange(len(marker[0])):
		maf = 0.0
		for geno in marker:
			if het:
				if geno[index]=='1':
					maf+=0.5/len(marker)
				elif geno[index]=='2':
					maf+=1.0/len(marker)
# 					if dominant:
# 						d_p+=1.0/len(marker)
			else:
				if geno[index]=='1':
					maf+=1.0/len(marker)
		mafs.append(maf)


	causal_indices=[]
	causal_mafs=[]
	mean_effects=[]
	var_effects=[]



# 	set causal snp info	
	for snp in causal_snps:
		snp[0] = int(snp[0])
		if snp[0] > -1:
			if(snp[0] > len(marker[0])):
				raise ValueError, 'SNP index out of range'

			if(snp[0] in causal_indices):
				raise ValueError, 'Redundant SNP index'
			
			causal_indices.append(snp[0])
			causal_mafs.append(mafs[snp[0]])
			mean_effects.append(float(snp[1]))
			var_effects.append(float(snp[2]))	

		elif snp[0] == -1:
			visited = []
			maf = 0.0
			while maf<=maf_range[0] or maf>=maf_range[1]:
				if len(visited) == len(marker[0]):
					raise RuntimeError, 'No SNP in MAF range found'

				causal_index=random.randint(0,len(marker[0]))

				if causal_index in causal_indices:
					continue

				if causal_index in visited:
					continue

				visited.append(causal_index)

				maf = mafs[causal_index]

			causal_indices.append(causal_index)
			causal_mafs.append(mafs[causal_index])
			mean_effects.append(float(snp[1]))
			var_effects.append(float(snp[2]))
		else:
			raise RuntimeError, 'Wrong index for SNP'


# 	calculate alpha, phi, and sigma
			
	tot_var
	alpha=[]
	phi=[]
	dom_prob = list(map(lambda x:x**2, causal_mafs))
	for i in xrange(len(causal_indices)):
		if het:
			if dominant:
				alpha.append(math.sqrt( (mean_effects[i] * tot_var) / ((dom_prob[i])*(1-dom_prob[i]))))				
# 				phi.append(math.sqrt( (var_effects[i] * tot_var) / ((dom_prob[i])*(1-dom_prob[i]))))				
				phi.append((var_effects[i] * tot_var) / (dom_prob[i]))
			else:
				alpha.append(math.sqrt( (mean_effects[i] * tot_var) / (2*(causal_mafs[i])*(1-causal_mafs[i]))))
# 				phi.append(math.sqrt( (var_effects[i] * tot_var) / (2*(causal_mafs[i])*(1-causal_mafs[i]))))
				phi.append((var_effects[i] * tot_var) / (2*(causal_mafs[i])))
		else:
			alpha.append(math.sqrt( (mean_effects[i] * tot_var) / ((causal_mafs[i])*(1-causal_mafs[i]))))
# 			phi.append(math.sqrt( (var_effects[i] * tot_var) / ((causal_mafs[i])*(1-causal_mafs[i]))))
			phi.append( (var_effects[i] * tot_var) / (causal_mafs[i]))


# 	if dominant:
# 		base_var = math.sqrt( tot_var*(1-(sum(mean_effects)+sum(var_effects)))) - np.dot(dom_prob,phi)
# 	else:
# 		max_effects = 1 - (np.dot(causal_mafs,phi))**2 / tot_var
# 		if( max_effects < (sum(mean_effects)+sum(var_effects))):
# 			raise RuntimeError, 'For these settings, sum of effects should be less than '+str(max_effects)+' current var: '+str((sum(mean_effects)+sum(var_effects)))
# 		
# 
# 		base_var = math.sqrt( tot_var*(1-(sum(mean_effects)+sum(var_effects)))) - np.dot(causal_mafs,phi)

	base_var = tot_var*(1-(sum(mean_effects)+sum(var_effects)))

# 	calculate the phenotype values


	phenotypes=[]   
	for geno in marker:             
		alleles=[]		
		for i in causal_indices:
			if het:
				if dominant:
					if geno[i] == '2':
						alleles.append(2)
					else:
						alleles.append(0)
				else:
					alleles.append(geno[i])
			else:
				if geno[i] == '1':
					alleles.append(2)
				else:
					alleles.append(0)			

		alleles=map(float,alleles)
		
		tot_alpha = np.dot(alleles, alpha)
		tot_phi = np.dot(alleles, phi)

		if(ld):
			cov_alpha = calc_alpha_cov(alleles=alleles, alpha = alpha, ld_matrix = ld_map)
			cov_phi = calc_phi_cov(alleles=alleles, phi = phi, ld_matrix = ld_map)
			tot_alpha = tot_alpha - cov_alpha
			tot_phi = tot_phi - cov_phi
			

		counter = 0
		if min != float("inf") and max != float("inf"):
			phen = base_avg + tot_alpha + random.gauss(0, math.sqrt( base_var+tot_phi))
			while( phen > max or phen < min):
				phen = base_avg + tot_alpha + random.gauss(0, math.sqrt( base_var+tot_phi))			
				counter=counter+1
				if counter > 10000:
					break
		elif min != float("inf"):
			phen = base_avg + tot_alpha + random.gauss(0, math.sqrt( base_var+tot_phi))
			while( phen < min):
				phen = base_avg + tot_alpha + random.gauss(0, math.sqrt( base_var+tot_phi))					
				counter=counter+1
				if counter > 10000:
					break
		elif max != float("inf"):
			phen = base_avg + tot_alpha + random.gauss(0, math.sqrt( base_var+tot_phi))
			while( phen > max):
				phen = base_avg + tot_alpha + random.gauss(0, math.sqrt( base_var+tot_phi))			
				counter=counter+1
				if counter > 10000:
					break
		else:
			phen = base_avg + np.dot(alleles, alpha) + random.gauss(0, math.sqrt( base_var+tot_phi))
		
		if counter > 10000:
			raise ValueError, 'Cannot find phenotype value within range'

		phenotypes.append(phen)
	
	return causal_indices,phenotypes,causal_mafs, alpha, phi






# 		for geno in marker:
# 			if geno[causal_index]=='1':
# 				maf+=1.0/len(marker)
# 
# 			
# 		mafs.append(maf)
# 		causal_indices.append(causal_index)		
# 			
# 	i=0
# 	while i<len(mean_eff) and not causal_index_pre:
# 		if i==0:			
# 			causal_indices=[]
# 			mafs=[]		
# 			dom_prob=[]
# 			visited = []
# 		i+=1
# 		maf=0.0
# 
# 		count=0
# 
# 		while maf<=maf_range[0] or maf>=maf_range[1]:       
# 
# 			if len(visited) == len(marker[0]):
# 				raise RuntimeError, 'No SNP in MAF range found'
# 				
# # 			count+=1
# # 
# # 			if count>=2*len(marker[0]):
# # 				raise RuntimeError, 'No SNP in MAF range found'
# 			
# 			maf=0.0
# 			epi_freq=0.0
# 			d_p = 0.0
# 
# # 			causal_index=random.randint(int(0.1*len(marker[0])),int(0.9*len(marker[0])))
# 
# 			causal_index=random.randint(0,len(marker[0]))
# 
# 			if causal_index in causal_indices:
# 				continue
# 
# 			if causal_index in visited:
# 				continue
# 
# 			visited.append(causal_index)
# 
# 			
# 
# 		dom_prob.append(d_p)
# 		mafs.append(maf)
# 		causal_indices.append(causal_index)
# 		#print causal_indices
# 
# 	phenotypes=[]   
# 	
# 
# # 	if (math.sqrt(1-sum(map(abs,mean_eff))))<0:
# # 		raise ValueError,'Illegal Variance Proportions'
# 	
# 	tot_var
# 	alpha=[]
# 	phi=[]
# 	dom_prob = list(map(lambda x:x**2, mafs))
# 	for i in xrange(len(causal_indices)):
# 		if diploid:
# 			if dominant:
# 				alpha.append(math.sqrt( (mean_eff[i] * tot_var) / ((dom_prob[i])*(1-dom_prob[i]))))				
# 				phi.append(math.sqrt( (var_eff[i] * tot_var) / ((dom_prob[i])*(1-dom_prob[i]))))				
# 			else:
# 				alpha.append(math.sqrt( (mean_eff[i] * tot_var) / (2*(mafs[i])*(1-mafs[i]))))
# 				phi.append(math.sqrt( (var_eff[i] * tot_var) / (2*(mafs[i])*(1-mafs[i]))))
# 		else:
# 			alpha.append(math.sqrt( (mean_eff[i] * tot_var) / ((mafs[i])*(1-mafs[i]))))
# 			phi.append(math.sqrt( (var_eff[i] * tot_var) / ((mafs[i])*(1-mafs[i]))))
# 
# 	if dominant:
# 		base_var = math.sqrt( tot_var*(1-(sum(mean_eff)+sum(var_eff)))) - np.dot(dom_prob,phi)
# 	else:
# 		base_var = math.sqrt( tot_var*(1-(sum(mean_eff)+sum(var_eff)))) - np.dot(mafs,phi)
# 
# 	for geno in marker:             
# 		alleles=[]
# 		
# 		for i in causal_indices:
# 			if dominant:
# 				if geno[i] == '2':
# 					alleles.append(1)
# 				else:
# 					alleles.append(0)
# 			else:
# 				alleles.append(geno[i])
# 
# 		alleles=map(float,alleles)
# 		
# 		phen = base_avg + np.dot(alleles, alpha) + random.gauss(0, math.sqrt( base_var+np.dot(alleles, phi)))
# 		phenotypes.append(phen)
# 	
# 
# 	return causal_indices,phenotypes,mafs    



