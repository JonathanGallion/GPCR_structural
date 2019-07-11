import matplotlib.pyplot as plt
import os,math,sys, numpy as np
from collections import defaultdict
import scipy
from scipy import stats
import random
from scipy.stats import mannwhitneyu
from scipy.stats import ks_2samp
from scipy.stats import ttest_ind





def readinEA():
	EAdict=defaultdict(list)
	file=open('../DATA/P07500_EAtable.txt')
	linecount=0
	for line in file:
		if linecount!=0:
		
			parts=line.strip('\n').split('\t')
			#print(parts[0])
			origaa=parts[0][0]
			position=parts[0][:]
			EAdict[position]=[parts[1],]
			for i in parts[2:]:
				EAdict[position].append(i)
			#print(EAdict[position])
		else:
			linecount=1
	return(EAdict)

def obtainhighcov():
	highcovlist=[]
	file=open('../NumberingConversiontoAlignment.txt')
	for line in file:
		#print(line)
		parts=line.strip('\n').split('\t')
		b2arpos=int(parts[1])
		try:
			highcov=parts[2]
		except:
			highcov='NA'
		if highcov=='HighCov':
			highcovlist.append(b2arpos)
	return(highcovlist)

def obtainETscores():
	counter=0
	file=open('../DATA/P07500.csv')
	etdict={}
	for line in file:
		if counter>0:
			#print(line)
			parts=line.strip('\n').split(',')
			print(parts)
			res=parts[1]
			coverage=parts[3]
			print(res,coverage)
			etdict[int(res)]=float(coverage)
		counter+=1

	return(etdict)

def getlistofsigpos(cancertype):
	listofsigpos=[]
	lowsigpos=[]
	file=open('../ImagesOlfactory/%s/PositionSignificance_%s.txt'%(cancertype,cancertype))
	for line in file:
		parts=line.strip('\n').split('\t')
		pos=parts[1]
		sig=float(parts[3])
		if pos!='.':
			pos=int(pos)
			if sig>0.99:
				listofsigpos.append(pos)
			if sig<0.01:
				lowsigpos.append(pos)
	return(listofsigpos,lowsigpos)


def acquiremutations(cancertype):
	eamutationdict=defaultdict(list)
	eamutlist=[]
	maffile=open('../ImagesOlfactory/%s/Allmutations_Ignoreexpression_Aug8%s.txt'%(cancertype,cancertype))
	for line in maffile:
		parts=line.strip('\n').split('\t')
		position=parts[1].strip(' ')
		expression=parts[7]
		EAscore=float(parts[4])
		gpcrcount=0
		if position!='.':
			if int(position) in eamutationdict:
				eamutationdict[int(position)].append(EAscore)
			else:
				eamutationdict[int(position)]=[EAscore,]
		eamutlist.append(EAscore)
	return(eamutationdict,eamutlist)
	



def etandeamain():
	aaconversiondict={'Ala':'A','Cys':'C','Asp':'D','Glu':'E','Phe':'F','Gly':'G','His':'H','Ile':'I','Lys':'K','Leu':'L','Met':'M','Asn':'N','Pro':'P','Gln':'Q','Arg':'R','Ser':'S','Thr':'T','Val':'V','Trp':'W','Tyr':'Y'}
	toaaarray=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	EAdict=readinEA()
	highcovlist=obtainhighcov()
	etdict=obtainETscores()
	
	cancerfile=open('../DATA/Mutations/filelist.txt')
	cancerlist=[]
	for line in cancerfile:
		cancer=line.strip('\n')
		if cancer!='':
			cancerlist.append(cancer)
	print(cancerlist)
	cancerlist=['COAD','HNSC','LUAD','LUSC','SKCM','STAD','UCEC']
	#cancerlist=['LUAD','SKCM']
	#print(EAdict)
	#print(etdict)
	#print('highcovlist',highcovlist)
	etallarray=[]
	for i in etdict:
		etallarray.append(float(etdict[i]))
		
	for cancertype in cancerlist:
		print('cancertype',cancertype)
		eamutationdict,eamutlist=acquiremutations(cancertype)
		listofsigpos,lowsigpos=getlistofsigpos(cancertype)
		#print(eamutationdict)
		ethyperarray=[]
		ethypoarray=[]
		earandarray=[]
		eahyperarray=[]
		eahypoarray=[]
		eaallarray=eamutlist
		for pos in EAdict:
			position=int(pos[1:])
			etvalue=etdict[int(position)]
			eaarray=EAdict[pos]
			if position in highcovlist:
			 	''
			for j in EAdict[pos]:
				if j!='-':
					earandarray.append(float(j))
		print('listofsigpos',listofsigpos)
		print('lowsigpos',lowsigpos)
		for pos in eamutationdict:#Getting EA the significantly mutated positions, hypo and hyper
			#print(pos)
			if pos in listofsigpos:
				''
				for ea in eamutationdict[pos]:
					if eahyperarray==[]:
						eahyperarray=[float(ea),]
					else:
						eahyperarray.append(float(ea))

			if pos in lowsigpos:
				for ea in eamutationdict[pos]:
					if eahyperarray==[]:
						eahypoarray=[float(ea),]
					else:
						eahypoarray.append(float(ea))
		for pos in listofsigpos:
			if pos in etdict:
				ethyperarray.append(float(etdict[pos]))
		for pos in lowsigpos:
			if pos in etdict:
				ethypoarray.append(float(etdict[pos]))
	

		etplist=[]
		ettuple=[ethyperarray,ethypoarray,etallarray]
		for i in ettuple:
			for j in ettuple:
				try:
					t,p=ttest_ind(i,j)
				except:
					p=1.00
				#t,p=mannwhitneyu(i,j)
				#t,p=ks_2samp(i,j)
				etplist.append(round(p,3))

		plot1=plt.boxplot(ettuple,whis=1,patch_artist=True,)
		plt.setp(plot1['boxes'],color='Black',linewidth=2,facecolor='lightgray')
		plt.setp(plot1['whiskers'],         # customise whisker appearence
				 color='black',       # whisker colour
				 linewidth=1)             # whisker thickness

		plt.setp(plot1['caps'],             # customize lines at the end of whiskers 
				 color='black',       # cap colour
				 linewidth=1)             # cap thickness

		plt.setp(plot1['fliers'],           # customize marks for extreme values
				 color='white',            # set mark colour
				 marker='o',                # maker shape
				 markersize=4)             # marker size

		plt.setp(plot1['medians'],          # customize median lines
				 color='Black',            # line colour
				 linewidth=2)             # line thickness
		plt.title('ET comparison of %s \n%s'%(cancertype,etplist[-3:]))
		plt.xticks([1,2,3],['Frequent Positions','Infrequent Positions','Random Control'])
		plt.savefig('../Paper/PotentialFigures/EAET/ET_%s.png'%(cancertype))
		plt.close()
		#plt.hist(ethyperarray,alpha=0.5)
		#plt.hist(ethypoarray,alpha=0.5)
		#plt.hist(etallarray,alpha=0.5)
		#plt.show()
		
		eaoutfile=open('../ImagesOlfactory/%s/EAforboxplot_%s.txt'%(cancertype,cancertype),'w')
		eatuple=[eahyperarray,eahypoarray,eaallarray,earandarray]
		for i in eatuple:
			for j in i:
				eaoutfile.write(str(j))
				eaoutfile.write('\t')
			eaoutfile.write('\n')
		eaplist=[]
		for i in eatuple:
			for j in eatuple:
				try:
					t,p=ttest_ind(i,j)
				except:
					p=1.00
				#t,p=mannwhitneyu(i,j)
				#t,p=ks_2samp(i,j)
				eaplist.append(round(p,5))
		plot1=plt.boxplot(eatuple,whis=1,patch_artist=True,)
		plt.setp(plot1['boxes'],color='Black',linewidth=2,facecolor='lightgray')
		plt.setp(plot1['whiskers'],         # customise whisker appearence
				 color='black',       # whisker colour
				 linewidth=1)             # whisker thickness

		plt.setp(plot1['caps'],             # customize lines at the end of whiskers 
				 color='black',       # cap colour
				 linewidth=1)             # cap thickness

		plt.setp(plot1['fliers'],           # customize marks for extreme values
				 color='white',            # set mark colour
				 marker='o',                # maker shape
				 markersize=4)             # marker size

		plt.setp(plot1['medians'],          # customize median lines
				 color='Black',            # line colour
				 linewidth=2)             # line thickness
		plt.title('EA comparison of %s\n%s'%(cancertype,eaplist[-4:]))
		plt.xticks([1,2,3,4],['Frequent Positions','Infrequent Positions','All Mutations','Random Control'])
		
		plt.ylim(0,100)
		plt.savefig('../Paper/PotentialFigures/EAET/EA_%s.png'%(cancertype))
		plt.close()
etandeamain()
