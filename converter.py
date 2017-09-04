def convert2plink(marker,positions,phenotype,fname,het):

	f=open('%s.ped'%fname,'w')

	for i in xrange(len(marker)):
		l='%s 1 0 0 0 %s ' %(i,phenotype[i])
		s=' '.join(map(str,marker[i]))+'\n'
		#print s
		if het:
			s=s.replace('2','A A')
			s=s.replace('1','A T')
			#print s
			s=s.replace('0','T T')
			s=s.replace('NA','0 0')
		else:
			s=s.replace('1','A A')
			#print s
			s=s.replace('0','T T')
			s=s.replace('NA','0 0')
		#print s
		f.write(l+s)

	f.close()

	f=open('%s.map'%fname,'w')

	for i in xrange(len(positions)):
		l='1 snp%s 0 %s\n' %(i,positions[i])
		f.write(l)

	f.close()

	f=open('%s.pheno'%fname,'w')

	for i in xrange(len(phenotype)):
		l='%s 1 %s\n' %(i,str(phenotype[i]))
		f.write(l)

	f.close()
	

	return

