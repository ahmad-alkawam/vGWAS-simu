# vgwas-simu

vGWAS-Simu is a tool for simulating phenotypes
from genotype data. vGWAS-Simu reads the output of the coalescent
genotype simulators and produces quantitative and qualitative phenotype
values. vGWAS-Simu provides the user with the ability to induce both:
mean and variance shifts resulting from on one or more causal loci.
vGWAS-Simu supports both homozygous and heterozygous diploid genotypes.
The phenotype could be generated using either a co-dominance or a
complete dominance model. The resulting phenotypes are outputted in
PLINK (Purcell et al, 2007) format to allow further analysis of the
data.

Setup
=====

The software was coded in Python 2.6 and uses some of Python’s standard
libraries. It was tested under both Linux and MacOS, but since Python is
portable, it should run on Windows as well. To use the tool:

1.  Unpack the vGWAS-Simu package.

2.  Add execution permission to vGWAS-Simu by typing: chmod +x
    vgwas-simu.py

3.  Make vGWAS-Simu available across your entire system through typing:
    export PATH=\$PATH:\$PWD. If you want to make this setting
    permanent, go to your home-directory and add the line export
    PATH=\$PATH:&lt;path to simwas&gt; to your .bashrc file.

For Windows:

1.  Make sure Python is installed on your computer.Python can be
    downloaded from http://www.python.org/.

2.  Unpack the vGWAS-Simu archive.

3.  Add execution permission to vGWAS-Simu by right-clicking on the
    file, select ’Properties’ and then uncheck the box ’Read-only’

4.  In order to make vGWAS-Simu available across your entire system, you
    have to add vGWAS-Simu’s directory to your system’s PATH.

Tool Description
================

Input
-----

vGWAS-Simu reads an input file defined using the –file option. This file
can be the output of ms (Hudson, 2002), msHOT (Hellenthal and Stephens,
2007), msms (Ewing and Hermisson, 2010), and GENOME (Liang et al. 2007).
The flag -i specifies the format of the input file by specifying ’M’ for
ms, msHOT, and msms, and ’G’ for GENOME. Under the default setting, the
individuals are treated as homozygous, by setting -h 1 they are treated
as heterozygous. Heterozygous individuals are created by combining two
simulated chromosomes to a joint genotype. A detailed list of inputs
along with their default values is given in Table S1.

Output
------

vGWAS-Simu writes the genotypes and phenotypes in file formats
compatible for PLINK (Purcell et al, 2007). The prefix of the output
files can be defined by the user (–outfile). The output consists of
three files: 1) A ’.map’ file containing each genotype’s chromosome
location, marker ID, genetic distance, and physical position. 2) A
’.ped’ file that describes the pedigree information of the genotypes and
their individuals. Currently, our simulator uses this file to assign
phenotype values to each genotype profile. However, we plan to
incorporate pedigree information in the subsequent versions of our
simulator. 3) The third file our simulator outputs is a ‘.causal’ file
which includes the position, index, minor allele frequency (MAF) and the
effect of the causal marker(s).

Phenotype Generation
--------------------

Our tool provides the user with the ability to introduce both: mean and
variance discrepancies based on one or more loci. The tool is also
capable of supporting both homozygous and heterozygous genotypes. For
heterozygotes, the phenotype could be generated using either the
co-dominance or the complete dominance model. 

Using the default settings, the genotypes are treated as heterozygotes
under co-dominance. By unsetting the heterozygous flag, (-h 0) the
genotypes are treated as homozygotes. When running the simulator in the
heterozygous case, the simulator reads the genotype input file and
considers every two consecutive genotype strings as the two chromosomes
of one individual. However, in the homozygous case and since the same
allele exists on both chromosomes, the simulator considers every
genotype string as a separate entry.

The user can control the value of several important parameters through
the command line inputs, such as the dominance model, the number of
loci, their minor allele frequencies, and their respective effect sizes.
For example, complete dominance could be simulated by setting the
dominance flag to one (- -dominance 1). A detailed list of inputs along
with their default values is given in Table S1. In order to account for
causal loci that are in LD, our simulator implements the mean and
variance shift adjustments described in the previous section. In such a
case, the simulator requires an $r^2$ matrix in a square format, as
outputted by PLINK. This matrix could be specified using the (- -ldfile)
flag.

The simulator is capable of generating qualitative phenotypes for
case/control studies using the liability model described above. Finally,
in the qualitative phenotype simulation, our tool provides the user with
the option of sampling the case and control population to produce
samples with a specified size. In typical scenarios, prevalence rates
are usually low, which results in a substantially larger control
population as compared to the case population. The sampling option gives
the user the ability to set the size of the output samples.

For more information about vGWAS-Simu, check the manual.pdf file.
