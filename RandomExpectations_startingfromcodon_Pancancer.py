import os,math,sys, numpy as np
import networkx as nx
np.set_printoptions(threshold=np.nan)
import matplotlib
import pylab
import matplotlib.pyplot as plt
from collections import defaultdict
import glob
import random

###This method is used to create a matrix of the expected liklihood of each amino acid changing to every other amino acid for a set of proteins(e.g. how often would we expect an Alanine to Valine to occur in all GPCRs)
##Given the frequency of nucleotide substitutions (e.g. C to T) in a particular cancer type
##Given the TCGA mutation file for that cancer type
###Given the Fasta files for the proteins of interest
####Output is a 20X20 (aa x aa) matrix of expected mutation rates

def CalcCodontable(cancertype,AAcodondict,mutorder):
	###This portion reads in an externale nucleotide transition matrix and uses the codon table to calculate the liklihood of each codon transitioning to every other codon
	###Only calculates odds a of a snv and not multiple polymorphisms in the same codon
	#Used the specific nucleotide transition frequency for each cancer type
	#####Reading in the Nucleotide transition frequencies caluclated at a seperate time for the cancer type of interest
	nuclist=['A','U','C','G','-']
	nucmat = []
	rowcounter=0
	nucmat=np.zeros((5,5))
	for line in open('../ImagesOlfactory/'+cancertype+'/NucTrans_GPCR_Frequency_Matrix.txt'):
		parts=line.strip('\n').split(' ')
		print(parts)
		colcounter=0
		for i in parts:
			nucmat[rowcounter][colcounter]=i
			colcounter+=1
		rowcounter+=1
	print('nucmat',nucmat)#This matrix is 5x5 and is in the order of nuclist
	return nucmat

def rosettaparse():
	genename_dict=defaultdict(list)
	Rosettastone=open('../Rosettastone_olf.txt')
	for line in Rosettastone:
		namepart=line.strip('\n').split('\t')
		#print(namepart)
		HGNC=namepart[2]
		GIid=namepart[3]
		nmid=namepart[0]
		genename_dict[GIid]=HGNC
	print('genename_dict',genename_dict)
	return genename_dict

def parsefasta():
	Fasta_dict={}
	file=open('../../NucleotideSequences/GPCR_sequences.txt')
	for line in file:
		#print(line)
		parts=line.strip('\n').split('\t')
		gene=parts[0]
		sequence=parts[1]
		Fasta_dict[gene]=sequence
	return Fasta_dict

def maphighcovtogenes(genename_dict,highcoveragearray):
	#This script creates a dictionary indexed by HGNC ID when available, GIid when not.
	#For each protein in the alignment I use it creates an array out '.' when that alignment position is a blank for this protein and a numeric value when that protein has an aa at that alingment position
	#This numeric value counts what protein position occurs at that alignment
	#This dictionary therefore translates the alignment position into the amino acid position for each protein.
	Protein_dict={}
	sequence_array=[]
	old_protein="cat"
	loop_count1=0
	FASTA=open("../Alignments/Promalalign_WithOlfactory.fa")
	for line in FASTA:
		#print(line)
		if line[0]=='>':
			nametemp=line.split('_')
			try:
				Protein_name=str(nametemp[1])
			except:
				''
			if loop_count1>0:
				try:
					Protein_dict[genename_dict[old_protein]]=sequence_array
				except:
					Protein_dict[old_protein]=sequence_array
			sequence_array=[]
			old_protein=str(Protein_name)
			pos_count=1
			loop_count1=loop_count1+1
			print(Protein_name)
		elif line[0]=='&':
			Protein_dict[old_protein]=sequence_array
		else:
			for position in line:
				if position != '\n':
					if position != '.':
						place=pos_count
						pos_count=pos_count+1
					else:
						place=position
					try:
						sequence_array.append(place)
					except:
						sequence_array=[place,]
	FASTA.close()
	#print('Protein_dict',Protein_dict)
	
	genposconvertdic=defaultdict(list)
	for position in highcoveragearray:
		for gene in Protein_dict:
			genposconvertdic[gene].append(Protein_dict[gene][position])
	print('genposconvertdic',genposconvertdic)
	for i in highcoveragearray:
		print(i,Protein_dict['CXCR2'][i])
	return(genposconvertdic,Protein_dict)
		
	
	




def splitintocodons(Fasta_dict,AAcodondict):
	##This part takes in the fasta dictionary and counts the number of times each codon occurs in that gene
	codoncountdict=defaultdict(list)
	codontotaldict={}
	for gene in Fasta_dict:
		if gene not in codoncountdict:
			codoncountdict[gene]=[]
			#print(codoncountdict[gene])
			positioncounter=0
			starpos=0
			endpos=3
			fasta=Fasta_dict[gene]
			#print(fasta)
			#print(codonlist)
			while endpos<len(fasta):
				codon=fasta[starpos:endpos]
				#print(starpos,endpos,len(fasta),codon)
				codoncountdict[gene].append(codon)
				try:
					codontotaldict[codon]+=1
				except:
					codontotaldict[codon]=1
				#print(fasta)
				#print(codoncountdict[gene])
				starpos+=3
				endpos+=3
	return codoncountdict,codontotaldict

def countcodonfreq(highcoveragearray,AAcodonlist,genposconvertdic,codoncountdict,badgenearray):
	#####Counting the codon frequency at each of the 265 positions across all GPCRs (the codon composition of each position based on the alignment)
	codonfreqmatrix=np.zeros((len(highcoveragearray),len(AAcodonlist)))
	codonrepresentationdict=defaultdict(list)
	countarray=range(0,len(highcoveragearray))
	codoncountdict_pos={}#This dictionary will count how often each of the 61 codons is occurs across all GPCRs at all positions
	print(countarray)
	poscounter=0
	for pos in countarray:
		codcount=0
		for codon in AAcodonlist:
			#print(pos,codon)
			#print(pos,poscounter,codcount)
			for gene in genposconvertdic:
				if gene in codoncountdict and gene not in badgenearray:
					#print(gene,pos,genposconvertdic[gene])
					if genposconvertdic[gene][pos]!='.':
						#print(codoncountdict[gene][genposconvertdic[gene][pos]-1])
						if codoncountdict[gene][genposconvertdic[gene][pos]-1]==codon:
							codonrepresentationdict[pos].append(codon)
							#print('These codons are a match',genposconvertdic[gene][pos],codon)
							codonfreqmatrix[poscounter][codcount]+=1
							try:
								codoncountdict_pos[codon]+=1
							except:
								codoncountdict_pos[codon]=1
			codcount+=1
		poscounter+=1
	for codon in codoncountdict_pos:
		print(codon,codoncountdict_pos[codon])
	print('codonrepresentationdict',codonrepresentationdict)

	return(codonrepresentationdict,codonfreqmatrix,codoncountdict_pos)






def Obtainmutations(cancertype):
	###This part loops over the mutation file from the tcga for the given cancer type and counts the number of unique mutations
	###in each gene for all genes of interest. These genes are listed in the array GPCR_protlist
	GPCR_protlist=['OR2M5','OR2A25','OR6C74','OR14I1','OR5H15','OR6C68','OR5K3','OR6C6','OR51F1','OR2T8','OR4C46','OR5H14','OR6C70','OR6C65','OR5H1','OR6C75','OR5B21','OR2AG2','OR6C76','OR5K4','OR2AT4','OR11H12','GPR25','F2RL2','GPR31','P2RY10','CCRL2','CXCR6','OR7A17','GPR171','RRH','TAAR5','FFAR1','FFAR3','GPR182','GPR37','FFAR2','MLNR','GPR39','GALR2','HCRTR1','HCRTR2','OR1F1','OR2T1','OR10H2','OR10H3','OR7C2','OR1I1','GALR3','GPR37L1','GPR32','LGR5','OR2B3','OR2J3','OR2J2','OR7C1','OR7A10','OR2F2','OR6B1','OR4F21','OR2A4','S1PR2','OR5F1','OR6A2','OR2C1','NTSR2','GPR75','OR2H2','S1PR4','OPN1SW','OPN1LW','OPN1MW2','MAS1','ADRB2','RHO','CHRM2','CHRM4','ADRB1','HTR1A','CHRM5','ADRA2A','OR56A5','CHRM1','ADRB3','DRD2','TSHR','ADRA2B','ADRA2C','CHRM3','TACR2','S1PR1','FPR1','CNR1','DRD1','TBXA2R','DRD4','DRD5','LHCGR','FSHR','EDNRB','HRH2','CXCR1','CXCR2','FPR3','FPR2','ADRA1D','EDNRA','TACR1','PTAFR','F2R','NPY1R','HTR1D','HTR1B','HTR2A','HTR2C','NMBR','HTR1E','ADORA2A','ADORA2B','TACR3','BDKRB2','AVPR2','ADORA1','GRPR','AGTR1','OXTR','SSTR1','SSTR2','HTR1F','OR1E1','OR10J1','GNRHR','NTSR1','SSTR4','CCKAR','CCKBR','MC4R','CCR1','BRS3','CCR7','GPR183','CXCR5','SSTR3','MC5R','ADORA3','HTR7','CNR2','TRHR','OR1D2','PTGER1','SSTR5','ADRA1A','HRH1','ADRA1B','OPRM1','PTGER4','MAS1L','APLNR','DRD3','AVPR1A','OPRD1','OPRK1','OPRL1','P2RY2','HTR2B','CCR2','MC3R','PTGFR','PTGER3','PTGER2','PTGIR','LPAR6','GPR3','GPR1','CCR10','GPR4','XCR1','GPR6','BDKRB1','GALR1','GPR12','RGR','OR3A1','OR1E2','OR3A3','OR1G1','OR3A2','HTR5A','P2RY1','AVPR1B','MTNR1A','NPBWR1','NPBWR2','HCAR3','NPY2R','CX3CR1','MTNR1B','CXCR3','PRLHR','GPR15','AGTR2','NPY4R','HTR6','P2RY4','CCR3','CCR4','CCR5','CCR6','CCR8','CCR9','F2RL1','OR1D5','OR2B6','OR4D2','OR10A3','OR12D2','GPR85','CXCR4','MC2R','MC1R','PTGDR','GPR17','GPR50','OR5I1','OR2F1','HTR4','GPR18','GPR176','P2RY6','P2RY14','OR1Q1','OR4D1','OR1C1','OR8B8','OR7A5','LTB4R','GPR68','GPR19','NPY5R','GPR162','C3AR1','OR2B11','OR10J3','FFAR4','OR2G6','GPR139','OR5M10','OR4C11','OR2T5','OR2T2','OR2A2','OR52W1','OR4A47','OR10K2','OR52E8','OR6B2','OR7E24','GPR153','ADORA3','OPN5','NPSR1','OR2W3','VN1R4','GPR142','GPR141','MRGPRG','MRGPRE','GPR149','P2RY8','OR4N5','GPR65','GPR135','OR4N4','OR8I2','OR5AS1','OR8H3','OR6V1','OR8H2','OR2L13','OR2C3','GPR161','GPR26','PROKR2','VN1R2','OR5T1','OR2T33','OR2T12','OR8G5','OR2L5','OR2M7','OR2M3','OR2AK2','OR2L3','OR13H1','OR11H1','OR7G3','OR2Z1','OR7D4','OR7G2','OR7G1','OR1M1','OR10H4','OR10H5','OR4F17','OR4S1','OR4M2','OR4F15','OR4F6','OR5AU1','OR11G2','OR4E2','OR10G2','OR10G3','OR4K17','OR11H6','OR11H4','OR4M1','OR4N2','OR4K2','OR4K5','OR4K1','OR4K14','OR10AD1','OR6C4','OR2AP1','OR10P1','OR10A7','OR9K2','OR4D9','OR9Q2','OR52B6','OR52R1','OR51D1','OR5AP2','OR10W1','OR5B17','OR4B1','OR4X2','OR8J3','OR5T2','OR5T3','OR8H1','OR8K1','OR8B12','OR8A1','OR8B3','OR2D3','OR56A1','OR52L1','OR56A4','OR52E4','OR52N2','OR52N4','OR56B1','OR4D11','OR4D10','OR10V1','OR5AN1','OR5A2','OR5A1','OR4D6','OR52H1','OR52E2','OR51L1','OR51A4','OR51A2','OR51S1','OR51T1','OR51G2','OR51G1','OR52B4','OR52K2','OR52K1','OR52M1','OR52I1','OR5D16','OR5L2','OR5D18','OR5L1','OR5D14','OR5D13','OR4A15','OR4P4','OR4C16','OR4C15','OR6M1','OR8D4','OR4D5','OR6T1','OR10S1','OR10G4','OR10G9','OR10G8','OR10G7','OR4C13','OR8J1','OR5M9','OR5M3','OR5M8','OR5M1','OR5AR1','OR9G4','OR6Q1','OR1S2','OR10Q1','OR9Q1','OR9I1','OR13A1','OR1L6','OR1K1','OR5C1','OR1L4','OR1B1','OR1L8','OR1N2','OR1N1','OR1J4','OR1J2','OR1J1','OR13F1','OR13C4','OR13C3','OR13C8','OR13C5','OR13C2','OR13C9','OR2K2','OR13J1','OR9A2','OR2A12','OR2A1','OR9A4','GPR150','OR2Y1','OR13D1','OR5H6','OR5H2','OR6B3','OR6K6','OR11L1','OR2T34','OR2T35','OR10T2','OR10K1','OR10R2','OR6Y1','OR6P1','OR10X1','OR10Z1','OR6K2','OR6K3','OR6N1','OR6N2','OR2L8','OR13G1','OR2G3','OR2G2','OR6F1','OR2T10','OR2T4','OR2T11','OR2T29','OR2T3','OR2T27','OR4Q3','OR11H2','OR8S1','OR8U1','OR2L2','OR5J2','OR10AG1','OR4F5','OR4C3','OR6S1','OR4K15','OR4K13','OR4L1','OR5B3','OR4X1','OR8K5','OR8K3','OR52N1','OR56A3','OR52N5','OR51Q1','OR52J3','OR51F2','OR51A7','OR52I2','OR5W2','OR4A16','OR4C6','OR4S2','OR10A6','OR56B4','OR6X1','OR4A5','OR5R1','OR9G1','OR5AK2','OR1S1','OR1L3','OR1L1','OR2AE1','OR2V1','OR5K1','OR5K2','OR10J5','OR14A16','OR14C36','OR2T6','OR51E1','PROKR1','HCAR2','OXER1','MRGPRD','GPR152','GPBAR1','RXFP4','GPR151','GPR148','GPR119','RXFP2','OR8D1','OR5P2','OR5P3','LPAR1','GHSR','KISS1R','TAAR8','MCHR2','MRGPRF','GPR146','OR10C1','MRGPRX4','MRGPRX3','MRGPRX2','MRGPRX1','QRFPR','GPR101','GPR82','OXGR1','GPR78','OR5B12','OR5B2','OR2M4','OR2M2','OR2V2','OR2A7','OR2A14','OR2A5','OR4C12','OR4F4','OR7D2','OR5M11','OR8B4','OR8B2','OR6C1','OR52B2','OR52E6','F2RL3','TAAR6','TAAR1','S1PR3','LPAR4','GPR20','GPR21','GPR22','MCHR1','CMKLR1','P2RY13','SUCNR1','LGR4','HCAR1','GPR174','GPR87','GPR63','GPR62','GPR61','OR2B2','OR2H1','OR11A1','OR8D2','GPR88','VN1R1','NMUR2','NPFFR1','LPAR5','OPN3','OR2AG1','OR10A5','OR10A2','OR10A4','OR2D2','S1PR5','P2RY12','OR51E2','OR52A5','OR51V1','OR51B5','OR51B6','OR51M1','OR51I1','OR51I2','OR52D1','HRH4','NMUR1','LPAR2','LGR6','RXFP1','GPR35','LTB4R2','OR2S2','GPR84','GPR173','GPR27','CYSLTR2','RXFP3','GPR83','OR6C3','OR6C2','OR5AC2','TAAR2','OR1A1','C5AR2','LPAR3','OR14J1','OR5V1','OR12D3','OPN4','GPR160','OR52A1','UTS2R','GPR132','GPR34','CYSLTR1','GPR52','GPR55','OR2W1','OR10H1','OR1A2','HRH3','OR51B4','OR51B2','NPFFR2','GPR45','PTGDR2']
	path='../DATA/Mutations/'+cancertype+'/'
	checkuniquedict=defaultdict(list)
	mutcountdict={}
	for filename in glob.glob(os.path.join(path,'*EA.maf')):
		print('filename',filename)
		maffile=open(filename)
		linecounter=0
		freqmatrix=np.zeros((5,5))
		gpcrfreqmatrix=np.zeros((5,5))
		for line in maffile:
			if linecounter!=0:
				parts=line.strip('\n').split('\t')
				PatientBC=parts[15]
				Mutationlist=parts[39]
				gene=parts[38]
				EAlist=parts[40]
				EAparts=EAlist.split(';')
				EAscore=EAparts[0]
				Mutparts=Mutationlist.split(';')
				mutation=Mutparts[0]
				uniquecheck=(PatientBC,gene,mutation)
				if gene in GPCR_protlist and uniquecheck not in checkuniquedict[gene]:
					#print(gene,mutation,EAscore)
					if gene in mutcountdict:
						mutcountdict[gene]+=1
					else:
						mutcountdict[gene]=1
			linecounter+=1
	return mutcountdict
	
	
def adjustforcodonfreq(codoncount,codonmatrix,genemut):	
	#Takes in the standard codon transition frequency and multiplies each row by the number of times that codon occurs in the protein divided by the total number of codons
	#print(codoncount)
	#print(codonmatrix)
	#print(sum(codoncount))
	codoncount=np.array(codoncount)/float(sum(codoncount))
	#print(codoncount)
	codonnormalizedmatrix=np.transpose(np.transpose(np.array(codonmatrix))*np.transpose(np.array(codoncount)))
	#print(codonnormalizedmatrix)
	codonnormalizedmatrix=codonnormalizedmatrix*float(genemut)
	#print('After multiplying by the number of mutations', genemut,codonnormalizedmatrix)
	
	
	return(codonnormalizedmatrix,codoncount)
	
def codonstoaa(Expectmatrix,codonlist,AAcodondict):	
	###This method takes a matrix of codon to codon transitions and converts it to an amino acid to amino acid transition matrix
	##Used the codon table in order to sum all transitions involving the same amino acid.
	aatransfinal=np.zeros((20,20))
	rowcounter=0
	for aa1 in AAcodondict:
		colcounter=0
		for aa2 in AAcodondict:##For every iterative combination of amino acids
			aarunningtotal=0
			for codon1 in AAcodondict[aa1]:##in those two amino acids, loop over all combinations of codon codon1 and codon2
				for codon2 in AAcodondict[aa2]:
					copos1=0
					copos2=0
					codoncounter=0
					for i in codonlist:##Loop over list to obtain which position codon1 and codon2 are at.
						if i[:3]==codon1:
							copos1=codoncounter
						if i[:3]==codon2:
							copos2=codoncounter
						codoncounter+=1
					print('codon1',codon1,copos1,codonlist[copos1],'codon2',codon2,copos2,codonlist[copos2],Expectmatrix[copos1][copos2])
					aarunningtotal=aarunningtotal+Expectmatrix[copos1][copos2]#Using the codon positions from above, acquire the transition frequency between codon1 and codon2 from matrix
			aatransfinal[rowcounter][colcounter]=aarunningtotal
			colcounter+=1
		rowcounter+=1
	print('aatransfinal',aatransfinal)
	return(aatransfinal)
	

def calculatemutationfrequency(cancerlist,GPCR_protlist,Olf_protlist,NoOlf_protlist):
	###Loops over each cancer type and parses 2 files to count the total number of mutations observed in each cancertype in GPCRs
	###also counts the number of non-synonymous muations
	###Need to make functionallity for Olf or No_olf option
	nummutdict={}
	observedcountdict={}
	for cancer in cancerlist:
		try:
			file=open('../ImagesOlfactory/%s/Allmutations_Ignoreexpression_Aug8%s.txt'%(cancer,cancer))
			filesilent=open('../ImagesOlfactory/%s/Allmutations_Ignoreexpression_Silentmutations_%s.txt'%(cancer,cancer))
			check1=True
		except:
			check1=False
		if check1==True:
			mutcounter=0
			for line in file:
				if line.strip('\n')!='':
					parts=line.strip('\n').split('\t')
					protein=parts[2]
					EAscore=parts[4]
					mutcounter+=1
			observedcountdict[cancer]=mutcounter
			for line in filesilent:
				if line.strip('\n')!='':
					parts=line.strip('\n').split('\t')
					protein=parts[2]
					mutcounter+=1
			nummutdict[cancer]=mutcounter
	return (nummutdict,observedcountdict)
		
def mutsimprotein(mutatpositiondict,poslist,codonrepresentationdict,genposconvertdic,badgenearray,codoncountdict,codonmutdict,cancertype,mutorder,AAcodondict):
	###providing an example simulation of how mutations may be distributed in different proteins
	###due to the low number of mutations, this is only a potential example, there are 1000s of potential distributions
	###Only designed to represent the rough mutation frequency of individual proteins over all, not specific variations between specific proteins
	countarray=range(0,len(poslist))
	outfile=open('../ImagesOlfactory/%s/Sample_Expectedmutations_%s.txt'%(cancertype,cancertype),'w')
	for pos in countarray:
		nummutations=0
		while nummutations<round(mutatpositiondict[poslist[pos]]):
			###making array of codons at this position based on the representation
			randomcodon=random.choice(codonrepresentationdict[pos])
			protein=random.choice(genposconvertdic.keys())
			if protein not in badgenearray:
				position=genposconvertdic[protein][pos]
				if position!='.' and protein in codoncountdict:
					proteincodon=codoncountdict[protein][position-1]
					if randomcodon==proteincodon:
						print(pos,randomcodon,protein,position,proteincodon,codonmutdict[randomcodon])
						if codonmutdict[randomcodon]!=[]:
							substitution=random.choice(codonmutdict[randomcodon])
							for aa in mutorder:
								for aacodon in AAcodondict[aa]:
									if aacodon==randomcodon:
										oldaa=aa
							if oldaa!=substitution:
								mutation=''.join([oldaa,str(position),substitution])
								print(protein,position,randomcodon, oldaa,substitution,mutation,nummutations)
								outfile.write('%s\t%s\t%s\n'%(poslist[pos],protein,mutation))
								nummutations+=1

def createoutputfile(poslist,mutationfreqmatrix,observedmutationarray,mutorder,cancertype):
	###designed to take the final products of all the data and make a .txt file displaying expected vs observed mutation rate
	outfile=open('../ImagesOlfactory/%s/ObserveredVsExpectedMutationrate_%s.txt'%(cancertype,cancertype),'w')
	counter=0
	outfile.write('Position')
	for aa in mutorder:
		outfile.write('\t%s'%(str(aa)))
	outfile.write('\tExpected\tObserved\n')
	for pos in poslist:
		outfile.write(str(pos))
		for value in mutationfreqmatrix[counter]:
			outfile.write('\t%s'%(str(value)))
		outfile.write('\t%s'%(str(sum(mutationfreqmatrix[counter]))))
		outfile.write('\t%s\n'%(str(observedmutationarray[counter])))
		counter+=1
	
	
def acquireobsmutations(cancertype,poslist):
	###creating an ordered array of the number of observed mutations at each position for each cancer type
	###Need to build functionallity for variations in EA and olfactory status
	mutdict=defaultdict(list)
	print('cancertype',cancertype)
	file=open('../ImagesOlfactory/%s/PositionSignificance_%s.txt'%(cancertype,cancertype))
	for line in file:
		parts=line.strip('\n').split('\t')
		position=parts[1]
		counts=parts[2]
		mutdict[position]=counts
	observedmutationarray=[]
	for pos in poslist:
		if str(pos) in mutdict.keys():
			print(pos,mutdict[str(pos)])
			observedmutationarray.append(int(mutdict[str(pos)]))
		else:
			observedmutationarray.append(0)
			print('pos not found',pos)
	print('observedmutationarray',observedmutationarray)
	return(observedmutationarray)
	
	
def Main():
	highcoveragearray=[915, 917, 918, 919, 920, 921, 922, 923, 924, 1358, 1359, 1360, 1361, 1362, 1363, 1364, 1365, 1366, 1367, 1368, 1373, 1374, 1375, 1376, 1377, 1378, 1382, 1383, 1384, 1385, 1386, 1387, 1388, 1389, 1390, 1402, 1403, 1404, 1405, 1406, 1407, 1410, 1423, 1424, 1425, 1426, 1427, 1428, 1460, 1461, 1462, 1463, 1464, 1465, 1466, 1467, 1468, 1474, 1475, 1476, 1482, 1483, 1484, 1485, 1486, 1487, 1488, 1499, 1500, 1501, 1502, 1505, 1506, 1507, 1517, 1518, 1519, 1520, 1521, 1522, 1527, 1528, 1529, 1536, 1537, 1538, 1540, 1541, 1542, 1543, 1544, 1545, 1591, 1592, 1593, 1594, 1595, 1596, 1597, 1598, 1599, 1600, 1601, 1602, 1603, 1604, 1605, 1606, 1607, 1608, 1609, 1610, 1611, 1612, 1622, 1623, 1624, 1625, 1644, 1645, 1661, 1662, 1663, 1664, 1665, 1666, 1668, 1669, 1670, 1671, 1672, 1710, 1711, 1712, 1714, 1715, 1716, 1717, 1718, 1719, 1731, 1732, 1733, 1734, 1735, 1782, 1784, 1785, 1786, 1787, 1909, 1910, 1911, 1912, 1913, 1967, 1976, 1978, 1979, 1986, 1987, 1988, 1989, 1990, 2001, 2002, 2007, 2008, 2009, 2011, 2012, 2013, 2014, 2079, 2080, 2081, 2082, 2083, 2084, 2085, 2086, 2087, 2124, 2125, 2126, 2127, 2128, 2145, 2146, 2147, 2181, 2701, 2702, 2703, 2704, 2705, 2728, 2729, 2730, 2731, 2732, 2733, 2734, 2735, 2736, 2737, 2738, 2748, 2749, 2750, 2751, 2752, 2753, 2754, 2755, 2762, 2763, 2764, 2765, 2766, 2808, 2818, 2821, 2824, 2826, 2834, 2836, 2837, 2838, 2839, 2848, 2849, 2850, 2851, 2852, 2853, 2854, 2870, 2871, 2872, 2873, 2874, 2875, 2876, 2877, 2878, 2879, 2880, 2881, 2882, 2883, 2884, 2885, 2886, 2950, 2951, 2952, 2953, 2954, 2955, 2956, 2957, 2958, 2959,]
	AAcodondict={'A':['GCU','GCC','GCA','GCG'],'C':['UGU','UGC'],'D':['GAU','GAC'],'E':['GAA','GAG'],'F':['UUU','UUC'],'G':['GGU','GGC','GGA','GGG'],'H':['CAU','CAC'],'I':['AUU','AUC','AUA'],'K':['AAA','AAG'],'L':['UUA','UUG','CUU','CUC','CUA','CUG'],'M':['AUG'],'N':['AAU','AAC'],'P':['CCU','CCC','CCA','CCG'],'Q':['CAA','CAG'],'R':['AGA','AGG','CGU','CGC','CGA','CGG'],'S':['UCU','UCC','UCA','UCG','AGU','AGC'],'T':['ACU','ACC','ACA','ACG'],'V':['GUU','GUC','GUA','GUG'],'W':['UGG'],'Y':['UAU','UAC']}
	AAcodonlist=['GCU','GCC','GCA','GCG','UGU','UGC','GAU','GAC','GAA','GAG','UUU','UUC','GGU','GGC','GGA','GGG','CAU','CAC','AUU','AUC','AUA','AAA','AAG','UUA','UUG','CUU','CUC','CUA','CUG','AUG','AAU','AAC','CCU','CCC','CCA','CCG','CAA','CAG','AGA','AGG','CGU','CGC','CGA','CGG','UCU','UCC','UCA','UCG','AGU','AGC','ACU','ACC','ACA','ACG','GUU','GUC','GUA','GUG','UGG','UAU','UAC']
	mutorder=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	GPCR_protlist=['OR2M5','OR2A25','OR6C74','OR14I1','OR5H15','OR6C68','OR5K3','OR6C6','OR51F1','OR2T8','OR4C46','OR5H14','OR6C70','OR6C65','OR5H1','OR6C75','OR5B21','OR2AG2','OR6C76','OR5K4','OR2AT4','OR11H12','GPR25','F2RL2','GPR31','P2RY10','CCRL2','CXCR6','OR7A17','GPR171','RRH','TAAR5','FFAR1','FFAR3','GPR182','GPR37','FFAR2','MLNR','GPR39','GALR2','HCRTR1','HCRTR2','OR1F1','OR2T1','OR10H2','OR10H3','OR7C2','OR1I1','GALR3','GPR37L1','GPR32','LGR5','OR2B3','OR2J3','OR2J2','OR7C1','OR7A10','OR2F2','OR6B1','OR4F21','OR2A4','S1PR2','OR5F1','OR6A2','OR2C1','NTSR2','GPR75','OR2H2','S1PR4','OPN1SW','OPN1LW','OPN1MW2','MAS1','ADRB2','RHO','CHRM2','CHRM4','ADRB1','HTR1A','CHRM5','ADRA2A','OR56A5','CHRM1','ADRB3','DRD2','TSHR','ADRA2B','ADRA2C','CHRM3','TACR2','S1PR1','FPR1','CNR1','DRD1','TBXA2R','DRD4','DRD5','LHCGR','FSHR','EDNRB','HRH2','CXCR1','CXCR2','FPR3','FPR2','ADRA1D','EDNRA','TACR1','PTAFR','F2R','NPY1R','HTR1D','HTR1B','HTR2A','HTR2C','NMBR','HTR1E','ADORA2A','ADORA2B','TACR3','BDKRB2','AVPR2','ADORA1','GRPR','AGTR1','OXTR','SSTR1','SSTR2','HTR1F','OR1E1','OR10J1','GNRHR','NTSR1','SSTR4','CCKAR','CCKBR','MC4R','CCR1','BRS3','CCR7','GPR183','CXCR5','SSTR3','MC5R','ADORA3','HTR7','CNR2','TRHR','OR1D2','PTGER1','SSTR5','ADRA1A','HRH1','ADRA1B','OPRM1','PTGER4','MAS1L','APLNR','DRD3','AVPR1A','OPRD1','OPRK1','OPRL1','P2RY2','HTR2B','CCR2','MC3R','PTGFR','PTGER3','PTGER2','PTGIR','LPAR6','GPR3','GPR1','CCR10','GPR4','XCR1','GPR6','BDKRB1','GALR1','GPR12','RGR','OR3A1','OR1E2','OR3A3','OR1G1','OR3A2','HTR5A','P2RY1','AVPR1B','MTNR1A','NPBWR1','NPBWR2','HCAR3','NPY2R','CX3CR1','MTNR1B','CXCR3','PRLHR','GPR15','AGTR2','NPY4R','HTR6','P2RY4','CCR3','CCR4','CCR5','CCR6','CCR8','CCR9','F2RL1','OR1D5','OR2B6','OR4D2','OR10A3','OR12D2','GPR85','CXCR4','MC2R','MC1R','PTGDR','GPR17','GPR50','OR5I1','OR2F1','HTR4','GPR18','GPR176','P2RY6','P2RY14','OR1Q1','OR4D1','OR1C1','OR8B8','OR7A5','LTB4R','GPR68','GPR19','NPY5R','GPR162','C3AR1','OR2B11','OR10J3','FFAR4','OR2G6','GPR139','OR5M10','OR4C11','OR2T5','OR2T2','OR2A2','OR52W1','OR4A47','OR10K2','OR52E8','OR6B2','OR7E24','GPR153','ADORA3','OPN5','NPSR1','OR2W3','VN1R4','GPR142','GPR141','MRGPRG','MRGPRE','GPR149','P2RY8','OR4N5','GPR65','GPR135','OR4N4','OR8I2','OR5AS1','OR8H3','OR6V1','OR8H2','OR2L13','OR2C3','GPR161','GPR26','PROKR2','VN1R2','OR5T1','OR2T33','OR2T12','OR8G5','OR2L5','OR2M7','OR2M3','OR2AK2','OR2L3','OR13H1','OR11H1','OR7G3','OR2Z1','OR7D4','OR7G2','OR7G1','OR1M1','OR10H4','OR10H5','OR4F17','OR4S1','OR4M2','OR4F15','OR4F6','OR5AU1','OR11G2','OR4E2','OR10G2','OR10G3','OR4K17','OR11H6','OR11H4','OR4M1','OR4N2','OR4K2','OR4K5','OR4K1','OR4K14','OR10AD1','OR6C4','OR2AP1','OR10P1','OR10A7','OR9K2','OR4D9','OR9Q2','OR52B6','OR52R1','OR51D1','OR5AP2','OR10W1','OR5B17','OR4B1','OR4X2','OR8J3','OR5T2','OR5T3','OR8H1','OR8K1','OR8B12','OR8A1','OR8B3','OR2D3','OR56A1','OR52L1','OR56A4','OR52E4','OR52N2','OR52N4','OR56B1','OR4D11','OR4D10','OR10V1','OR5AN1','OR5A2','OR5A1','OR4D6','OR52H1','OR52E2','OR51L1','OR51A4','OR51A2','OR51S1','OR51T1','OR51G2','OR51G1','OR52B4','OR52K2','OR52K1','OR52M1','OR52I1','OR5D16','OR5L2','OR5D18','OR5L1','OR5D14','OR5D13','OR4A15','OR4P4','OR4C16','OR4C15','OR6M1','OR8D4','OR4D5','OR6T1','OR10S1','OR10G4','OR10G9','OR10G8','OR10G7','OR4C13','OR8J1','OR5M9','OR5M3','OR5M8','OR5M1','OR5AR1','OR9G4','OR6Q1','OR1S2','OR10Q1','OR9Q1','OR9I1','OR13A1','OR1L6','OR1K1','OR5C1','OR1L4','OR1B1','OR1L8','OR1N2','OR1N1','OR1J4','OR1J2','OR1J1','OR13F1','OR13C4','OR13C3','OR13C8','OR13C5','OR13C2','OR13C9','OR2K2','OR13J1','OR9A2','OR2A12','OR2A1','OR9A4','GPR150','OR2Y1','OR13D1','OR5H6','OR5H2','OR6B3','OR6K6','OR11L1','OR2T34','OR2T35','OR10T2','OR10K1','OR10R2','OR6Y1','OR6P1','OR10X1','OR10Z1','OR6K2','OR6K3','OR6N1','OR6N2','OR2L8','OR13G1','OR2G3','OR2G2','OR6F1','OR2T10','OR2T4','OR2T11','OR2T29','OR2T3','OR2T27','OR4Q3','OR11H2','OR8S1','OR8U1','OR2L2','OR5J2','OR10AG1','OR4F5','OR4C3','OR6S1','OR4K15','OR4K13','OR4L1','OR5B3','OR4X1','OR8K5','OR8K3','OR52N1','OR56A3','OR52N5','OR51Q1','OR52J3','OR51F2','OR51A7','OR52I2','OR5W2','OR4A16','OR4C6','OR4S2','OR10A6','OR56B4','OR6X1','OR4A5','OR5R1','OR9G1','OR5AK2','OR1S1','OR1L3','OR1L1','OR2AE1','OR2V1','OR5K1','OR5K2','OR10J5','OR14A16','OR14C36','OR2T6','OR51E1','PROKR1','HCAR2','OXER1','MRGPRD','GPR152','GPBAR1','RXFP4','GPR151','GPR148','GPR119','RXFP2','OR8D1','OR5P2','OR5P3','LPAR1','GHSR','KISS1R','TAAR8','MCHR2','MRGPRF','GPR146','OR10C1','MRGPRX4','MRGPRX3','MRGPRX2','MRGPRX1','QRFPR','GPR101','GPR82','OXGR1','GPR78','OR5B12','OR5B2','OR2M4','OR2M2','OR2V2','OR2A7','OR2A14','OR2A5','OR4C12','OR4F4','OR7D2','OR5M11','OR8B4','OR8B2','OR6C1','OR52B2','OR52E6','F2RL3','TAAR6','TAAR1','S1PR3','LPAR4','GPR20','GPR21','GPR22','MCHR1','CMKLR1','P2RY13','SUCNR1','LGR4','HCAR1','GPR174','GPR87','GPR63','GPR62','GPR61','OR2B2','OR2H1','OR11A1','OR8D2','GPR88','VN1R1','NMUR2','NPFFR1','LPAR5','OPN3','OR2AG1','OR10A5','OR10A2','OR10A4','OR2D2','S1PR5','P2RY12','OR51E2','OR52A5','OR51V1','OR51B5','OR51B6','OR51M1','OR51I1','OR51I2','OR52D1','HRH4','NMUR1','LPAR2','LGR6','RXFP1','GPR35','LTB4R2','OR2S2','GPR84','GPR173','GPR27','CYSLTR2','RXFP3','GPR83','OR6C3','OR6C2','OR5AC2','TAAR2','OR1A1','C5AR2','LPAR3','OR14J1','OR5V1','OR12D3','OPN4','GPR160','OR52A1','UTS2R','GPR132','GPR34','CYSLTR1','GPR52','GPR55','OR2W1','OR10H1','OR1A2','HRH3','OR51B4','OR51B2','NPFFR2','GPR45','PTGDR2']
	Olf_protlist=['OR10A2','OR10A3','OR10A4','OR10A5','OR10A6','OR10A7','OR10AD1','OR10AG1','OR10C1','OR10G2','OR10G3','OR10G4','OR10G7','OR10G8','OR10G9','OR10H1','OR10H2','OR10H3','OR10H4','OR10H5','OR10J1','OR10J3','OR10J5','OR10K1','OR10K2','OR10P1','OR10Q1','OR10R2','OR10S1','OR10T2','OR10V1','OR10W1','OR10X1','OR10Z1','OR11A1','OR11G2','OR11H1','OR11H12','OR11H2','OR11H4','OR11H6','OR11L1','OR12D2','OR12D3','OR13A1','OR13C2','OR13C3','OR13C4','OR13C5','OR13C8','OR13C9','OR13D1','OR13F1','OR13G1','OR13H1','OR13J1','OR14A16','OR14C36','OR14I1','OR14J1','OR1A1','OR1A2','OR1B1','OR1C1','OR1D2','OR1D5','OR1E1','OR1E2','OR1F1','OR1G1','OR1I1','OR1J1','OR1J2','OR1J4','OR1K1','OR1L1','OR1L3','OR1L4','OR1L6','OR1L8','OR1M1','OR1N1','OR1N2','OR1Q1','OR1S1','OR1S2','OR2A1','OR2A12','OR2A14','OR2A2','OR2A25','OR2A4','OR2A5','OR2A7','OR2AE1','OR2AG1','OR2AG2','OR2AK2','OR2AP1','OR2AT4','OR2B11','OR2B2','OR2B3','OR2B6','OR2C1','OR2C3','OR2D2','OR2D3','OR2F1','OR2F2','OR2G2','OR2G3','OR2G6','OR2H1','OR2H2','OR2J2','OR2J3','OR2K2','OR2L13','OR2L2','OR2L3','OR2L5','OR2L8','OR2M2','OR2M3','OR2M4','OR2M5','OR2M7','OR2S2','OR2T1','OR2T10','OR2T11','OR2T12','OR2T2','OR2T27','OR2T29','OR2T3','OR2T33','OR2T34','OR2T35','OR2T4','OR2T5','OR2T6','OR2T8','OR2V1','OR2V2','OR2W1','OR2W3','OR2Y1','OR2Z1','OR3A1','OR3A2','OR3A3','OR4A15','OR4A16','OR4A47','OR4A5','OR4B1','OR4C11','OR4C12','OR4C13','OR4C15','OR4C16','OR4C3','OR4C46','OR4C6','OR4D1','OR4D10','OR4D11','OR4D2','OR4D5','OR4D6','OR4D9','OR4E2','OR4F15','OR4F17','OR4F21','OR4F4','OR4F5','OR4F6','OR4K1','OR4K13','OR4K14','OR4K15','OR4K17','OR4K2','OR4K5','OR4L1','OR4M1','OR4M2','OR4N2','OR4N4','OR4N5','OR4P4','OR4Q3','OR4S1','OR4S2','OR4X1','OR4X2','OR51A2','OR51A4','OR51A7','OR51B2','OR51B4','OR51B5','OR51B6','OR51D1','OR51E1','OR51E2','OR51F1','OR51F2','OR51G1','OR51G2','OR51I1','OR51I2','OR51L1','OR51M1','OR51Q1','OR51S1','OR51T1','OR51V1','OR52A1','OR52A5','OR52B2','OR52B4','OR52B6','OR52D1','OR52E2','OR52E4','OR52E6','OR52E8','OR52H1','OR52I1','OR52I2','OR52J3','OR52K1','OR52K2','OR52L1','OR52M1','OR52N1','OR52N2','OR52N4','OR52N5','OR52R1','OR52W1','OR56A1','OR56A3','OR56A4','OR56A5','OR56B1','OR56B4','OR5A1','OR5A2','OR5AC2','OR5AK2','OR5AN1','OR5AP2','OR5AR1','OR5AS1','OR5AU1','OR5B12','OR5B17','OR5B2','OR5B21','OR5B3','OR5C1','OR5D13','OR5D14','OR5D16','OR5D18','OR5F1','OR5H1','OR5H14','OR5H15','OR5H2','OR5H6','OR5I1','OR5J2','OR5K1','OR5K2','OR5K3','OR5K4','OR5L1','OR5L2','OR5M1','OR5M10','OR5M11','OR5M3','OR5M8','OR5M9','OR5P2','OR5P3','OR5R1','OR5T1','OR5T2','OR5T3','OR5V1','OR5W2','OR6A2','OR6B1','OR6B2','OR6B3','OR6C1','OR6C2','OR6C3','OR6C4','OR6C6','OR6C65','OR6C68','OR6C70','OR6C74','OR6C75','OR6C76','OR6F1','OR6K2','OR6K3','OR6K6','OR6M1','OR6N1','OR6N2','OR6P1','OR6Q1','OR6S1','OR6T1','OR6V1','OR6X1','OR6Y1','OR7A10','OR7A17','OR7A5','OR7C1','OR7C2','OR7D2','OR7D4','OR7E24','OR7G1','OR7G2','OR7G3','OR8A1','OR8B12','OR8B2','OR8B3','OR8B4','OR8B8','OR8D1','OR8D2','OR8D4','OR8G5','OR8H1','OR8H2','OR8H3','OR8I2','OR8J1','OR8J3','OR8K1','OR8K3','OR8K5','OR8S1','OR8U1','OR9A2','OR9A4','OR9G1','OR9G4','OR9I1','OR9K2','OR9Q1','OR9Q2']
	NoOlf_protlist=['ADORA1', 'ADORA2A', 'ADORA2B', 'ADORA3', 'ADORA3', 'ADRA1A', 'ADRA1B', 'ADRA1D', 'ADRA2A', 'ADRA2B', 'ADRA2C', 'ADRB1', 'ADRB2', 'ADRB3', 'AGTR1', 'AGTR2', 'APLNR', 'AVPR1A', 'AVPR1B', 'AVPR2', 'BDKRB1', 'BDKRB2', 'BRS3', 'C3AR1', 'C5AR2', 'CCKAR', 'CCKBR', 'CCR1', 'CCR10', 'CCR2', 'CCR3', 'CCR4', 'CCR5', 'CCR6', 'CCR7', 'CCR8', 'CCR9', 'CCRL2', 'CHRM1', 'CHRM2', 'CHRM3', 'CHRM4', 'CHRM5', 'CMKLR1', 'CNR1', 'CNR2', 'CX3CR1', 'CXCR1', 'CXCR2', 'CXCR3', 'CXCR4', 'CXCR5', 'CXCR6', 'CYSLTR1', 'CYSLTR2', 'DRD1', 'DRD2', 'DRD3', 'DRD4', 'DRD5', 'EDNRA', 'EDNRB', 'F2R', 'F2RL1', 'F2RL2', 'F2RL3', 'FFAR1', 'FFAR2', 'FFAR3', 'FFAR4', 'FPR1', 'FPR2', 'FPR3', 'FSHR', 'GALR1', 'GALR2', 'GALR3', 'GHSR', 'GNRHR', 'GPBAR1', 'GPR1', 'GPR101', 'GPR119', 'GPR12', 'GPR132', 'GPR135', 'GPR139', 'GPR141', 'GPR142', 'GPR146', 'GPR148', 'GPR149', 'GPR15', 'GPR150', 'GPR151', 'GPR152', 'GPR153', 'GPR160', 'GPR161', 'GPR162', 'GPR17', 'GPR171', 'GPR173', 'GPR174', 'GPR176', 'GPR18', 'GPR182', 'GPR183', 'GPR19', 'GPR20', 'GPR21', 'GPR22', 'GPR25', 'GPR26', 'GPR27', 'GPR3', 'GPR31', 'GPR32', 'GPR34', 'GPR35', 'GPR37', 'GPR37L1', 'GPR39', 'GPR4', 'GPR45', 'GPR50', 'GPR52', 'GPR55', 'GPR6', 'GPR61', 'GPR62', 'GPR63', 'GPR65', 'GPR68', 'GPR75', 'GPR78', 'GPR82', 'GPR83', 'GPR84', 'GPR85', 'GPR87', 'GPR88', 'GRPR', 'HCAR1', 'HCAR2', 'HCAR3', 'HCRTR1', 'HCRTR2', 'HRH1', 'HRH2', 'HRH3', 'HRH4', 'HTR1A', 'HTR1B', 'HTR1D', 'HTR1E', 'HTR1F', 'HTR2A', 'HTR2B', 'HTR2C', 'HTR4', 'HTR5A', 'HTR6', 'HTR7', 'KISS1R', 'LGR4', 'LGR5', 'LGR6', 'LHCGR', 'LPAR1', 'LPAR2', 'LPAR3', 'LPAR4', 'LPAR5', 'LPAR6', 'LTB4R', 'LTB4R2', 'MAS1', 'MAS1L', 'MC1R', 'MC2R', 'MC3R', 'MC4R', 'MC5R', 'MCHR1', 'MCHR2', 'MLNR', 'MRGPRD', 'MRGPRE', 'MRGPRF', 'MRGPRG', 'MRGPRX1', 'MRGPRX2', 'MRGPRX3', 'MRGPRX4', 'MTNR1A', 'MTNR1B', 'NMBR', 'NMUR1', 'NMUR2', 'NPBWR1', 'NPBWR2', 'NPFFR1', 'NPFFR2', 'NPSR1', 'NPY1R', 'NPY2R', 'NPY4R', 'NPY5R', 'NTSR1', 'NTSR2', 'OPN1LW', 'OPN1MW2', 'OPN1SW', 'OPN3', 'OPN4', 'OPN5', 'OPRD1', 'OPRK1', 'OPRL1', 'OPRM1', 'OXER1', 'OXGR1', 'OXTR', 'P2RY1', 'P2RY10', 'P2RY12', 'P2RY13', 'P2RY14', 'P2RY2', 'P2RY4', 'P2RY6', 'P2RY8', 'PRLHR', 'PROKR1', 'PROKR2', 'PTAFR', 'PTGDR', 'PTGDR2', 'PTGER1', 'PTGER2', 'PTGER3', 'PTGER4', 'PTGFR', 'PTGIR', 'QRFPR', 'RGR', 'RHO', 'RRH', 'RXFP1', 'RXFP2', 'RXFP3', 'RXFP4', 'S1PR1', 'S1PR2', 'S1PR3', 'S1PR4', 'S1PR5', 'SSTR1', 'SSTR2', 'SSTR3', 'SSTR4', 'SSTR5', 'SUCNR1', 'TAAR1', 'TAAR2', 'TAAR5', 'TAAR6', 'TAAR8', 'TACR1', 'TACR2', 'TACR3', 'TBXA2R', 'TRHR', 'TSHR', 'UTS2R', 'VN1R1', 'VN1R2', 'VN1R4', 'XCR1']
	badgenearray=['PROKR1','OPN3']
	#GPCR_protlist=['CCRL2','CXCR6']
	nuclist=['A','U','C','G','-']
	nuclist=['A','U','C','G']
	cancerlist=['ACC','LUAD','SKCM','STAD','HNSC','COAD','UCEC','BLCA','BRCA','ESCA','GBM','KIRC','KIRP','LIHC','LUSC','PAAD','READ','SARC','STAD','THCA','UCEC','UCS']
	print(sorted(AAcodondict.items()))
	
	syncount=0
	nonsyncount=0
	#cancertype='{0}'.format(sys.argv[1])
	genename_dict=rosettaparse()
	Fasta_dict=parsefasta()
	print('Fasta_dict',Fasta_dict)
	codoncountdict,codontotaldict=splitintocodons(Fasta_dict,AAcodondict)
	genposconvertdic,Protein_dict=maphighcovtogenes(genename_dict,highcoveragearray)
	#print('codoncountdict',codoncountdict['ADRB2'][143],codoncountdict['ADRB2'][64],codoncountdict['ADRB2'][310])
	#print(genposconvertdic['HRH4'][103],genposconvertdic['HRH4'][25],genposconvertdic['HRH4'][242])
	#print('codoncountdict',codoncountdict['HRH4'][genposconvertdic['HRH4'][103]],codoncountdict['HRH4'][genposconvertdic['HRH4'][25]],codoncountdict['HRH4'][genposconvertdic['HRH4'][242]])
	
	for i in codoncountdict:
		try:
			''		
			print(i,codoncountdict[i][genposconvertdic[i][103]],codoncountdict[i][genposconvertdic[i][24]],codoncountdict[i][genposconvertdic[i][242]])
		except:
			print(i,'There was some error with this gene')
	
	codonrepresentationdict,codonfreqmatrix,codoncountdict_pos=countcodonfreq(highcoveragearray,AAcodonlist,genposconvertdic,codoncountdict,badgenearray)
	
	###Obtain the number of mutations in each cancer type
	nummutdict,observedcountdict=calculatemutationfrequency(cancerlist,GPCR_protlist,Olf_protlist,NoOlf_protlist)
	
	#Calculating the expected mutations frequency for each cancertype
	for cancertype in cancerlist:
		np.savetxt('../ImagesOlfactory/%s/codonfrequency.txt'%(cancertype),codonfreqmatrix)##Count of how frequently each codon occurs at each of the 265 positions of interest across all the GPCRs
		
		mutationdict=defaultdict(list)
		codonmutdict=defaultdict(list)
		
		nummut=nummutdict[cancertype]
		nucmat=CalcCodontable(cancertype,AAcodondict,mutorder)
		nucmat=nummut*nucmat
		print('nucmat',nucmat)
		nuccounter=0
		for i in nuclist:
			nuccounter1=0
			for j in nuclist:
				if i!=j:
					repcounter=0##This will count the number of repitions the method makes at a single nucleotide transition type
					numreps=round(nucmat[nuccounter][nuccounter1],0)
					print(i,j,'numreps',numreps)
					while repcounter<numreps:##Start method loop
						codon=random.choice(AAcodonlist)
						print(codon)
						codonpos=random.choice([0,1,2])
						if codon[codonpos]==i:
							newcodon=list(codon)
							newcodon[codonpos]=j
							newcodon=''.join(newcodon)
							for aa in AAcodondict:
								for aacodon in AAcodondict[aa]:
									if aacodon==newcodon:
										newaa=aa
									if aacodon==codon:
										oldaa=aa
							try:
								print(oldaa,codon,codonpos,codon[codonpos],i,newcodon,newaa)
							except:
								print(newcodon)
							codonmutdict[codon].append(newaa)
							repcounter+=1
				nuccounter1+=1
			nuccounter+=1
		print('codonmutdict',codonmutdict)#Structure is 61 codons for the keys and each key has a list of single letter codes for AA, number of times that codon was mutated to that AA including silent mutations
	
	
		#This section counts the number of non-silent mutations at each codon. 'total'=non silent len(codonmutdict[codon])-total = number of silent mutations
		scaledcodondict={}
		for codon in AAcodonlist:
			if codon in codonmutdict.keys():###some codons weren't mutated in each cancer type, so first check if those codons have any mutations
				for aa in AAcodondict:
					for aacodon in AAcodondict[aa]:
						if aacodon==codon:
							aminoacid=aa
				total=0
				for i in codonmutdict[codon]:
					if aminoacid!=i:
						nonsyncount+=1
						total+=1
					else:
						syncount+=1

				print(codon,'non-synonymous=',total,'all mutations=',len(codonmutdict[codon]))
				scaledcodondict[codon]=float(total)/float(codoncountdict_pos[codon]) ## number of non-synomous mutations / # times that codon occurs?  Gives the non-synomous mutation rate of each codon for a single codon
			else:#if not mutations were observed at that codon, setting the scaled frequency to 0
				scaledcodondict[codon]=0
		print('scaledcodondict',scaledcodondict)
		
		print('codonmutdict',codonmutdict)
		for i in codonmutdict:
			print(i,len(codonmutdict[i]),codoncountdict_pos[i],scaledcodondict[i])
		print('codontotaldict',codontotaldict)
		for aa in AAcodondict:#converting the codon mutation expectation to the mutation expectation for each amino acid (combine multiple codons that code for the same AA)
			total=0
			for codon in AAcodondict[aa]:
				total=total+scaledcodondict[codon]
			print(aa,total/float(len(AAcodondict[aa])))
	
	
		#How often each codon occurs at each positions * the expected mutation rate of that codon= The number of mutations at that position for that codon. Sum over all codons at that position for the expected mutation rate.
		print(sorted(AAcodondict.items()))
		mutationfreqmatrix=np.zeros((len(highcoveragearray),20))
		rowcounter=0
		for row in codonfreqmatrix:
			print(row)
			codoncounter=0
			for i in row:
				codon=AAcodonlist[codoncounter]
				aacounter=0
				#print(codon)
				for aa in mutorder:
					#print(aa)
					for aacodon in AAcodondict[aa]:
						if aacodon==codon:
							#print(aacodon,codon,aacounter)
							colpos=aacounter
					aacounter+=1
				#print('colpos',colpos)
				scale=scaledcodondict[codon]
				mutationfreqmatrix[rowcounter][colpos]=mutationfreqmatrix[rowcounter][colpos]+i*scale
				#print(i)
			
				codoncounter+=1	
			rowcounter+=1
	
		###Create dictionary of the total number of mutations at each position
		poslist=['25',	'26',	'27',	'28',	'29',	'30',	'31',	'32',	'33',	'34',	'35',	'36',	'37',	'38',	'39',	'40',	'41',	'42',	'43',	'44',	'45',	'46',	'47',	'48',	'49',	'50',	'51',	'52',	'53',	'54',	'55',	'56',	'57',	'58',	'59',	'60',	'61',	'62',	'63',	'64',	'65',	'66',	'67',	'68',	'69',	'70',	'71',	'72',	'73',	'74',	'75',	'76',	'77',	'78',	'79',	'80',	'81',	'82',	'83',	'84',	'86',	'87',	'88',	'89',	'90',	'91',	'92',	'93',	'94',	'95',	'96',	'98',	'99',	'100',	'101',	'102',	'103',	'104',	'105',	'106',	'107',	'108',	'109',	'110',	'111',	'112',	'113',	'114',	'115',	'116',	'117',	'118',	'119',	'120',	'121',	'122',	'123',	'124',	'125',	'126',	'127',	'128',	'129',	'130',	'131',	'132',	'133',	'134',	'135',	'136',	'137',	'138',	'139',	'140',	'141',	'142',	'143',	'144',	'147',	'148',	'149',	'150',	'151',	'152',	'153',	'154',	'155',	'156',	'157',	'158',	'159',	'160',	'161',	'162',	'163',	'164',	'165',	'166',	'167',	'168',	'169',	'170',	'171',	'172',	'173',	'174',	'175',	'176',	'177',	'178',	'182',	'183',	'184',	'185',	'186',	'194',	'197',	'198',	'199',	'200',	'201',	'202',	'203',	'204',	'205',	'206',	'207',	'208',	'209',	'210',	'211',	'212',	'213',	'214',	'215',	'216',	'217',	'218',	'219',	'220',	'221',	'222',	'223',	'224',	'225',	'226',	'227',	'228',	'229',	'230',	'231',	'263',	'264',	'265',	'266',	'267',	'271',	'272',	'273',	'274',	'275',	'276',	'277',	'278',	'279',	'280',	'281',	'282',	'283',	'284',	'285',	'286',	'287',	'288',	'289',	'290',	'291',	'292',	'293',	'294',	'300',	'301',	'302',	'304',	'305',	'306',	'307',	'308',	'309',	'310',	'311',	'312',	'313',	'314',	'315',	'316',	'317',	'318',	'319',	'320',	'321',	'322',	'323',	'324',	'325',	'326',	'327',	'328',	'329',	'330',	'331',	'332',	'333',	'334',	'335',	'336',	'337',	'338',	'339',	'340',	'341',	'342',	'343',	'344']
		mutatpositiondict={}
		poscounter=0
		for pos in poslist:
			mutatpositiondict[pos]=sum(mutationfreqmatrix[poscounter])
			poscounter+=1
		print(mutatpositiondict)

		
		###Normalize the mutation freqmatrix based on the number of predicted non-syn mutations compared to the number of observed
		mutationfreqmatrix=mutationfreqmatrix*float(observedcountdict[cancertype])/sum(sum(mutationfreqmatrix))
		np.savetxt('../ImagesOlfactory/%s/ExpectedMutFreq.txt'%(cancertype),mutationfreqmatrix)
		
		####obtaining the observedmutation rate at each position
		observedmutationarray=acquireobsmutations(cancertype,poslist)
		
		####Creating a nice output with the expected mutation rate of each position vs the observed mutation rate
		createoutputfile(poslist,mutationfreqmatrix,observedmutationarray,mutorder,cancertype)
		
		
		
		
		##Now that we have the number of mutations and type of mutations at each position, distribute these across all proteins in order to simulate an actual mutation pattern
		mutsimprotein(mutatpositiondict,poslist,codonrepresentationdict,genposconvertdic,badgenearray,codoncountdict,codonmutdict,cancertype,mutorder,AAcodondict)
	print(syncount,nonsyncount)
	
Main()