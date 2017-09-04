import re

def _add_lists(x,y):
	#use numpy instead?
	target=[]
	assert len(x)==len(y)
	
	for i in xrange(len(x)):
		target.append(x[i]+y[i])
		
	return target
		
def make_diploids(genotypes):

	dipl_genotypes=[]
	for i in xrange(0,len(genotypes),2):
		if i+1>=len(genotypes):
			continue
		x=map(int,genotypes[i])
		y=map(int,genotypes[i+1])
		
		z=_add_lists(x,y)
		dipl_genotypes.append(map(str,z))
		
	
	genotypes=dipl_genotypes
	return genotypes

def parse_ms(fname,diploid=0):
	f=open(fname)
	s=f.read()
	f.close()

	positions_all=[]
	genotypes_all=[]

	simulations=s.split('//')[1:]

	for sim in simulations:
		simlines=sim.split('\n')
		genotypes,positions=parse_ms_simwise(simlines,diploid)

		positions_all.append(positions)
		genotypes_all.append(genotypes)

	return genotypes_all,positions_all



def parse_ms_simwise(lines,diploid=0):

	positions=[]
	genotypes=[]

	sim_count=0

	for l in lines:
		if len(l)==0:
			continue
		if l[:3]=='pos':
			s=re.sub('positions: ','',l.rstrip())
			positions=map(float,s.split(' '))

		if l[0] in ['0','1']:
			genotypes.append(list(l.rstrip()))
					
	if diploid:
		genotypes=make_diploids(genotypes)
		

	return genotypes,positions

def parse_genome(fname,diploid=0):
	f=open(fname)
	s=f.read()
	f.close()

	positions_all=[]
	genotypes_all=[]

	simulations=s.split('GENOME')[1:]

	for sim in simulations:
		simlines=sim.split('\n')
		genotypes,positions=parse_genome_simwise(simlines,diploid)
		positions_all.append(positions)
		genotypes_all.append(genotypes)

	return genotypes_all,positions_all
	

def parse_genome_simwise(lines,diploid=0):

	positions=[]
	genotypes=[]

	sim_count=0
	pos_next=0

	for i in xrange(len(lines)):
		l=lines[i]
		if len(l)==0:
			continue
		if l[:3]=='SNP':
			pos_next=1
		elif pos_next:
			s=l.rstrip()
			positions=map(float,s.split(' '))
			pos_next=0

		if l[:3]=='POP':
			s=re.sub('POP.*: ','',l.rstrip())
			genotypes.append(list(s))


	for i in xrange(len(positions)):
		positions[i]=float(positions[i])
	
		
	if diploid:
		genotypes=make_diploids(genotypes)


	return genotypes,positions



def parse_snps(fname):
	f=open(fname)
	snps_all=f.read().split('\n')
	f.close()
	causal_snps = []

	for s in snps_all:
		snp = s.split('\t')
		if( len(snp) != 3):
			continue
		if(int(snp[0]) != -1):
			snp[0] = int(snp[0])-1
		causal_snps.append(snp)
		
	return causal_snps
	
def	parse_ld(ldfile, nsnps):
	ld_matrix = []
	f=open(ldfile)
	ld_all=f.read().split('\n')
	for l in ld_all:	
		ld_snp = l.split('\t')
		if( len(ld_snp) == nsnps):
			ld_matrix.append(ld_snp)
		else:
			if( len(ld_snp) != 1 or str(ld_snp[0]) != ""):
				raise ValueError, 'Incorrect dimensions for ld matrix'
	if( len(ld_matrix) != nsnps):
		raise ValueError, 'Incorrect dimensions for ld matrix'	
	f.close()
	return ld_matrix
	
	
	