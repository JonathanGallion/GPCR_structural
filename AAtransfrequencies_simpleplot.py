import os,math,sys, numpy as np
import csv
import networkx as nx
np.set_printoptions(threshold=np.nan)
import matplotlib
import pylab
import xlrd # Import the package
import matplotlib.pyplot as plt
from collections import Counter
from collections import defaultdict
import glob
from pylab import *
import random
import operator


#####The idea behind this code is to show simplistically



def getlistofsigpos(cancertype):
	listofsigpos=[]
	lowsigpos=[]
	file=open('../ImagesOlfactory/%s/PositionSignificance_%s.txt'%(cancertype,cancertype))
	sigcounter=0
	sigdict={}
	for line in file:
		parts=line.strip('\n').split('\t')
		pos=parts[1]
		count=parts[2]
		sig=float(parts[3])
		if pos!='.':
			sigdict[int(pos)]=int(count)
		
			pos=int(pos)
			if sig>0.99 and sigcounter<50:
				listofsigpos.append(pos)
				
				sigcounter+=1
	sorted_x = sorted(sigdict.items(), key=operator.itemgetter(1))
	print(sorted_x)
	listofsigpos=[]
	for i in sorted_x[-10:]:
		listofsigpos.append(i[0])
	print(sorted(listofsigpos))
	return(sorted(listofsigpos))


def mutperpos(cancertype,listofsigpos):
	file=open('../ImagesOlfactory/%s/Allmutations_Ignoreexpression_Aug8%s.txt'%(cancertype,cancertype))
	GPCR_protlist=['OR2M5','OR2A25','OR6C74','OR14I1','OR5H15','OR6C68','OR5K3','OR6C6','OR51F1','OR2T8','OR4C46','OR5H14','OR6C70','OR6C65','OR5H1','OR6C75','OR5B21','OR2AG2','OR6C76','OR5K4','OR2AT4','OR11H12','GPR25','F2RL2','GPR31','P2RY10','CCRL2','CXCR6','OR7A17','GPR171','RRH','TAAR5','FFAR1','FFAR3','GPR182','GPR37','FFAR2','MLNR','GPR39','GALR2','HCRTR1','HCRTR2','OR1F1','OR2T1','OR10H2','OR10H3','OR7C2','OR1I1','GALR3','GPR37L1','GPR32','LGR5','OR2B3','OR2J3','OR2J2','OR7C1','OR7A10','OR2F2','OR6B1','OR4F21','OR2A4','S1PR2','OR5F1','OR6A2','OR2C1','NTSR2','GPR75','OR2H2','S1PR4','OPN1SW','OPN1LW','OPN1MW2','MAS1','ADRB2','RHO','CHRM2','CHRM4','ADRB1','HTR1A','CHRM5','ADRA2A','OR56A5','CHRM1','ADRB3','DRD2','TSHR','ADRA2B','ADRA2C','CHRM3','TACR2','S1PR1','FPR1','CNR1','DRD1','TBXA2R','DRD4','DRD5','LHCGR','FSHR','EDNRB','HRH2','CXCR1','CXCR2','FPR3','FPR2','ADRA1D','EDNRA','TACR1','PTAFR','F2R','NPY1R','HTR1D','HTR1B','HTR2A','HTR2C','NMBR','HTR1E','ADORA2A','ADORA2B','TACR3','BDKRB2','AVPR2','ADORA1','GRPR','AGTR1','OXTR','SSTR1','SSTR2','HTR1F','OR1E1','OR10J1','GNRHR','NTSR1','SSTR4','CCKAR','CCKBR','MC4R','CCR1','BRS3','CCR7','GPR183','CXCR5','SSTR3','MC5R','ADORA3','HTR7','CNR2','TRHR','OR1D2','PTGER1','SSTR5','ADRA1A','HRH1','ADRA1B','OPRM1','PTGER4','MAS1L','APLNR','DRD3','AVPR1A','OPRD1','OPRK1','OPRL1','P2RY2','HTR2B','CCR2','MC3R','PTGFR','PTGER3','PTGER2','PTGIR','LPAR6','GPR3','GPR1','CCR10','GPR4','XCR1','GPR6','BDKRB1','GALR1','GPR12','RGR','OR3A1','OR1E2','OR3A3','OR1G1','OR3A2','HTR5A','P2RY1','AVPR1B','MTNR1A','NPBWR1','NPBWR2','HCAR3','NPY2R','CX3CR1','MTNR1B','CXCR3','PRLHR','GPR15','AGTR2','NPY4R','HTR6','P2RY4','CCR3','CCR4','CCR5','CCR6','CCR8','CCR9','F2RL1','OR1D5','OR2B6','OR4D2','OR10A3','OR12D2','GPR85','CXCR4','MC2R','MC1R','PTGDR','GPR17','GPR50','OR5I1','OR2F1','HTR4','GPR18','GPR176','P2RY6','P2RY14','OR1Q1','OR4D1','OR1C1','OR8B8','OR7A5','LTB4R','GPR68','GPR19','NPY5R','GPR162','C3AR1','OR2B11','OR10J3','FFAR4','OR2G6','GPR139','OR5M10','OR4C11','OR2T5','OR2T2','OR2A2','OR52W1','OR4A47','OR10K2','OR52E8','OR6B2','OR7E24','GPR153','ADORA3','OPN5','NPSR1','OR2W3','VN1R4','GPR142','GPR141','MRGPRG','MRGPRE','GPR149','P2RY8','OR4N5','GPR65','GPR135','OR4N4','OR8I2','OR5AS1','OR8H3','OR6V1','OR8H2','OR2L13','OR2C3','GPR161','GPR26','PROKR2','VN1R2','OR5T1','OR2T33','OR2T12','OR8G5','OR2L5','OR2M7','OR2M3','OR2AK2','OR2L3','OR13H1','OR11H1','OR7G3','OR2Z1','OR7D4','OR7G2','OR7G1','OR1M1','OR10H4','OR10H5','OR4F17','OR4S1','OR4M2','OR4F15','OR4F6','OR5AU1','OR11G2','OR4E2','OR10G2','OR10G3','OR4K17','OR11H6','OR11H4','OR4M1','OR4N2','OR4K2','OR4K5','OR4K1','OR4K14','OR10AD1','OR6C4','OR2AP1','OR10P1','OR10A7','OR9K2','OR4D9','OR9Q2','OR52B6','OR52R1','OR51D1','OR5AP2','OR10W1','OR5B17','OR4B1','OR4X2','OR8J3','OR5T2','OR5T3','OR8H1','OR8K1','OR8B12','OR8A1','OR8B3','OR2D3','OR56A1','OR52L1','OR56A4','OR52E4','OR52N2','OR52N4','OR56B1','OR4D11','OR4D10','OR10V1','OR5AN1','OR5A2','OR5A1','OR4D6','OR52H1','OR52E2','OR51L1','OR51A4','OR51A2','OR51S1','OR51T1','OR51G2','OR51G1','OR52B4','OR52K2','OR52K1','OR52M1','OR52I1','OR5D16','OR5L2','OR5D18','OR5L1','OR5D14','OR5D13','OR4A15','OR4P4','OR4C16','OR4C15','OR6M1','OR8D4','OR4D5','OR6T1','OR10S1','OR10G4','OR10G9','OR10G8','OR10G7','OR4C13','OR8J1','OR5M9','OR5M3','OR5M8','OR5M1','OR5AR1','OR9G4','OR6Q1','OR1S2','OR10Q1','OR9Q1','OR9I1','OR13A1','OR1L6','OR1K1','OR5C1','OR1L4','OR1B1','OR1L8','OR1N2','OR1N1','OR1J4','OR1J2','OR1J1','OR13F1','OR13C4','OR13C3','OR13C8','OR13C5','OR13C2','OR13C9','OR2K2','OR13J1','OR9A2','OR2A12','OR2A1','OR9A4','GPR150','OR2Y1','OR13D1','OR5H6','OR5H2','OR6B3','OR6K6','OR11L1','OR2T34','OR2T35','OR10T2','OR10K1','OR10R2','OR6Y1','OR6P1','OR10X1','OR10Z1','OR6K2','OR6K3','OR6N1','OR6N2','OR2L8','OR13G1','OR2G3','OR2G2','OR6F1','OR2T10','OR2T4','OR2T11','OR2T29','OR2T3','OR2T27','OR4Q3','OR11H2','OR8S1','OR8U1','OR2L2','OR5J2','OR10AG1','OR4F5','OR4C3','OR6S1','OR4K15','OR4K13','OR4L1','OR5B3','OR4X1','OR8K5','OR8K3','OR52N1','OR56A3','OR52N5','OR51Q1','OR52J3','OR51F2','OR51A7','OR52I2','OR5W2','OR4A16','OR4C6','OR4S2','OR10A6','OR56B4','OR6X1','OR4A5','OR5R1','OR9G1','OR5AK2','OR1S1','OR1L3','OR1L1','OR2AE1','OR2V1','OR5K1','OR5K2','OR10J5','OR14A16','OR14C36','OR2T6','OR51E1','PROKR1','HCAR2','OXER1','MRGPRD','GPR152','GPBAR1','RXFP4','GPR151','GPR148','GPR119','RXFP2','OR8D1','OR5P2','OR5P3','LPAR1','GHSR','KISS1R','TAAR8','MCHR2','MRGPRF','GPR146','OR10C1','MRGPRX4','MRGPRX3','MRGPRX2','MRGPRX1','QRFPR','GPR101','GPR82','OXGR1','GPR78','OR5B12','OR5B2','OR2M4','OR2M2','OR2V2','OR2A7','OR2A14','OR2A5','OR4C12','OR4F4','OR7D2','OR5M11','OR8B4','OR8B2','OR6C1','OR52B2','OR52E6','F2RL3','TAAR6','TAAR1','S1PR3','LPAR4','GPR20','GPR21','GPR22','MCHR1','CMKLR1','P2RY13','SUCNR1','LGR4','HCAR1','GPR174','GPR87','GPR63','GPR62','GPR61','OR2B2','OR2H1','OR11A1','OR8D2','GPR88','VN1R1','NMUR2','NPFFR1','LPAR5','OPN3','OR2AG1','OR10A5','OR10A2','OR10A4','OR2D2','S1PR5','P2RY12','OR51E2','OR52A5','OR51V1','OR51B5','OR51B6','OR51M1','OR51I1','OR51I2','OR52D1','HRH4','NMUR1','LPAR2','LGR6','RXFP1','GPR35','LTB4R2','OR2S2','GPR84','GPR173','GPR27','CYSLTR2','RXFP3','GPR83','OR6C3','OR6C2','OR5AC2','TAAR2','OR1A1','C5AR2','LPAR3','OR14J1','OR5V1','OR12D3','OPN4','GPR160','OR52A1','UTS2R','GPR132','GPR34','CYSLTR1','GPR52','GPR55','OR2W1','OR10H1','OR1A2','HRH3','OR51B4','OR51B2','NPFFR2','GPR45','PTGDR2']
	Olf_protlist=['OR10A2','OR10A3','OR10A4','OR10A5','OR10A6','OR10A7','OR10AD1','OR10AG1','OR10C1','OR10G2','OR10G3','OR10G4','OR10G7','OR10G8','OR10G9','OR10H1','OR10H2','OR10H3','OR10H4','OR10H5','OR10J1','OR10J3','OR10J5','OR10K1','OR10K2','OR10P1','OR10Q1','OR10R2','OR10S1','OR10T2','OR10V1','OR10W1','OR10X1','OR10Z1','OR11A1','OR11G2','OR11H1','OR11H12','OR11H2','OR11H4','OR11H6','OR11L1','OR12D2','OR12D3','OR13A1','OR13C2','OR13C3','OR13C4','OR13C5','OR13C8','OR13C9','OR13D1','OR13F1','OR13G1','OR13H1','OR13J1','OR14A16','OR14C36','OR14I1','OR14J1','OR1A1','OR1A2','OR1B1','OR1C1','OR1D2','OR1D5','OR1E1','OR1E2','OR1F1','OR1G1','OR1I1','OR1J1','OR1J2','OR1J4','OR1K1','OR1L1','OR1L3','OR1L4','OR1L6','OR1L8','OR1M1','OR1N1','OR1N2','OR1Q1','OR1S1','OR1S2','OR2A1','OR2A12','OR2A14','OR2A2','OR2A25','OR2A4','OR2A5','OR2A7','OR2AE1','OR2AG1','OR2AG2','OR2AK2','OR2AP1','OR2AT4','OR2B11','OR2B2','OR2B3','OR2B6','OR2C1','OR2C3','OR2D2','OR2D3','OR2F1','OR2F2','OR2G2','OR2G3','OR2G6','OR2H1','OR2H2','OR2J2','OR2J3','OR2K2','OR2L13','OR2L2','OR2L3','OR2L5','OR2L8','OR2M2','OR2M3','OR2M4','OR2M5','OR2M7','OR2S2','OR2T1','OR2T10','OR2T11','OR2T12','OR2T2','OR2T27','OR2T29','OR2T3','OR2T33','OR2T34','OR2T35','OR2T4','OR2T5','OR2T6','OR2T8','OR2V1','OR2V2','OR2W1','OR2W3','OR2Y1','OR2Z1','OR3A1','OR3A2','OR3A3','OR4A15','OR4A16','OR4A47','OR4A5','OR4B1','OR4C11','OR4C12','OR4C13','OR4C15','OR4C16','OR4C3','OR4C46','OR4C6','OR4D1','OR4D10','OR4D11','OR4D2','OR4D5','OR4D6','OR4D9','OR4E2','OR4F15','OR4F17','OR4F21','OR4F4','OR4F5','OR4F6','OR4K1','OR4K13','OR4K14','OR4K15','OR4K17','OR4K2','OR4K5','OR4L1','OR4M1','OR4M2','OR4N2','OR4N4','OR4N5','OR4P4','OR4Q3','OR4S1','OR4S2','OR4X1','OR4X2','OR51A2','OR51A4','OR51A7','OR51B2','OR51B4','OR51B5','OR51B6','OR51D1','OR51E1','OR51E2','OR51F1','OR51F2','OR51G1','OR51G2','OR51I1','OR51I2','OR51L1','OR51M1','OR51Q1','OR51S1','OR51T1','OR51V1','OR52A1','OR52A5','OR52B2','OR52B4','OR52B6','OR52D1','OR52E2','OR52E4','OR52E6','OR52E8','OR52H1','OR52I1','OR52I2','OR52J3','OR52K1','OR52K2','OR52L1','OR52M1','OR52N1','OR52N2','OR52N4','OR52N5','OR52R1','OR52W1','OR56A1','OR56A3','OR56A4','OR56A5','OR56B1','OR56B4','OR5A1','OR5A2','OR5AC2','OR5AK2','OR5AN1','OR5AP2','OR5AR1','OR5AS1','OR5AU1','OR5B12','OR5B17','OR5B2','OR5B21','OR5B3','OR5C1','OR5D13','OR5D14','OR5D16','OR5D18','OR5F1','OR5H1','OR5H14','OR5H15','OR5H2','OR5H6','OR5I1','OR5J2','OR5K1','OR5K2','OR5K3','OR5K4','OR5L1','OR5L2','OR5M1','OR5M10','OR5M11','OR5M3','OR5M8','OR5M9','OR5P2','OR5P3','OR5R1','OR5T1','OR5T2','OR5T3','OR5V1','OR5W2','OR6A2','OR6B1','OR6B2','OR6B3','OR6C1','OR6C2','OR6C3','OR6C4','OR6C6','OR6C65','OR6C68','OR6C70','OR6C74','OR6C75','OR6C76','OR6F1','OR6K2','OR6K3','OR6K6','OR6M1','OR6N1','OR6N2','OR6P1','OR6Q1','OR6S1','OR6T1','OR6V1','OR6X1','OR6Y1','OR7A10','OR7A17','OR7A5','OR7C1','OR7C2','OR7D2','OR7D4','OR7E24','OR7G1','OR7G2','OR7G3','OR8A1','OR8B12','OR8B2','OR8B3','OR8B4','OR8B8','OR8D1','OR8D2','OR8D4','OR8G5','OR8H1','OR8H2','OR8H3','OR8I2','OR8J1','OR8J3','OR8K1','OR8K3','OR8K5','OR8S1','OR8U1','OR9A2','OR9A4','OR9G1','OR9G4','OR9I1','OR9K2','OR9Q1','OR9Q2']
	NoOlf_protlist=['ADORA1', 'ADORA2A', 'ADORA2B', 'ADORA3', 'ADORA3', 'ADRA1A', 'ADRA1B', 'ADRA1D', 'ADRA2A', 'ADRA2B', 'ADRA2C', 'ADRB1', 'ADRB2', 'ADRB3', 'AGTR1', 'AGTR2', 'APLNR', 'AVPR1A', 'AVPR1B', 'AVPR2', 'BDKRB1', 'BDKRB2', 'BRS3', 'C3AR1', 'C5AR2', 'CCKAR', 'CCKBR', 'CCR1', 'CCR10', 'CCR2', 'CCR3', 'CCR4', 'CCR5', 'CCR6', 'CCR7', 'CCR8', 'CCR9', 'CCRL2', 'CHRM1', 'CHRM2', 'CHRM3', 'CHRM4', 'CHRM5', 'CMKLR1', 'CNR1', 'CNR2', 'CX3CR1', 'CXCR1', 'CXCR2', 'CXCR3', 'CXCR4', 'CXCR5', 'CXCR6', 'CYSLTR1', 'CYSLTR2', 'DRD1', 'DRD2', 'DRD3', 'DRD4', 'DRD5', 'EDNRA', 'EDNRB', 'F2R', 'F2RL1', 'F2RL2', 'F2RL3', 'FFAR1', 'FFAR2', 'FFAR3', 'FFAR4', 'FPR1', 'FPR2', 'FPR3', 'FSHR', 'GALR1', 'GALR2', 'GALR3', 'GHSR', 'GNRHR', 'GPBAR1', 'GPR1', 'GPR101', 'GPR119', 'GPR12', 'GPR132', 'GPR135', 'GPR139', 'GPR141', 'GPR142', 'GPR146', 'GPR148', 'GPR149', 'GPR15', 'GPR150', 'GPR151', 'GPR152', 'GPR153', 'GPR160', 'GPR161', 'GPR162', 'GPR17', 'GPR171', 'GPR173', 'GPR174', 'GPR176', 'GPR18', 'GPR182', 'GPR183', 'GPR19', 'GPR20', 'GPR21', 'GPR22', 'GPR25', 'GPR26', 'GPR27', 'GPR3', 'GPR31', 'GPR32', 'GPR34', 'GPR35', 'GPR37', 'GPR37L1', 'GPR39', 'GPR4', 'GPR45', 'GPR50', 'GPR52', 'GPR55', 'GPR6', 'GPR61', 'GPR62', 'GPR63', 'GPR65', 'GPR68', 'GPR75', 'GPR78', 'GPR82', 'GPR83', 'GPR84', 'GPR85', 'GPR87', 'GPR88', 'GRPR', 'HCAR1', 'HCAR2', 'HCAR3', 'HCRTR1', 'HCRTR2', 'HRH1', 'HRH2', 'HRH3', 'HRH4', 'HTR1A', 'HTR1B', 'HTR1D', 'HTR1E', 'HTR1F', 'HTR2A', 'HTR2B', 'HTR2C', 'HTR4', 'HTR5A', 'HTR6', 'HTR7', 'KISS1R', 'LGR4', 'LGR5', 'LGR6', 'LHCGR', 'LPAR1', 'LPAR2', 'LPAR3', 'LPAR4', 'LPAR5', 'LPAR6', 'LTB4R', 'LTB4R2', 'MAS1', 'MAS1L', 'MC1R', 'MC2R', 'MC3R', 'MC4R', 'MC5R', 'MCHR1', 'MCHR2', 'MLNR', 'MRGPRD', 'MRGPRE', 'MRGPRF', 'MRGPRG', 'MRGPRX1', 'MRGPRX2', 'MRGPRX3', 'MRGPRX4', 'MTNR1A', 'MTNR1B', 'NMBR', 'NMUR1', 'NMUR2', 'NPBWR1', 'NPBWR2', 'NPFFR1', 'NPFFR2', 'NPSR1', 'NPY1R', 'NPY2R', 'NPY4R', 'NPY5R', 'NTSR1', 'NTSR2', 'OPN1LW', 'OPN1MW2', 'OPN1SW', 'OPN3', 'OPN4', 'OPN5', 'OPRD1', 'OPRK1', 'OPRL1', 'OPRM1', 'OXER1', 'OXGR1', 'OXTR', 'P2RY1', 'P2RY10', 'P2RY12', 'P2RY13', 'P2RY14', 'P2RY2', 'P2RY4', 'P2RY6', 'P2RY8', 'PRLHR', 'PROKR1', 'PROKR2', 'PTAFR', 'PTGDR', 'PTGDR2', 'PTGER1', 'PTGER2', 'PTGER3', 'PTGER4', 'PTGFR', 'PTGIR', 'QRFPR', 'RGR', 'RHO', 'RRH', 'RXFP1', 'RXFP2', 'RXFP3', 'RXFP4', 'S1PR1', 'S1PR2', 'S1PR3', 'S1PR4', 'S1PR5', 'SSTR1', 'SSTR2', 'SSTR3', 'SSTR4', 'SSTR5', 'SUCNR1', 'TAAR1', 'TAAR2', 'TAAR5', 'TAAR6', 'TAAR8', 'TACR1', 'TACR2', 'TACR3', 'TBXA2R', 'TRHR', 'TSHR', 'UTS2R', 'VN1R1', 'VN1R2', 'VN1R4', 'XCR1']
	positiondict={}
	for line in file:
		parts=line.strip('\n').split('\t')
		protein=parts[2]
		position=parts[1]
		EAscore=float(parts[4])
		mutto=parts[3][-1]
		expression=float(parts[7])
		gpcrcount=0
		if position!='.' and int(position) in listofsigpos:
			position=int(position)
			if position in positiondict:
				positiondict[position].append(mutto)
			else:
				positiondict[position]=[mutto,]
	print('positiondict',positiondict)
	return(positiondict)


def piechart():
	mutorder=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	cancerlist=['COAD','HNSC','LUAD','LUSC','SKCM','STAD','UCEC']
	cancerlist=['SKCM','LUAD']
	#cancerlist=['SKCM','LUAD']
	#cancerlist=['Thousandgenomes',]
	for cancertype in cancerlist:
		print('cancertype',cancertype)
		listofsigpos=getlistofsigpos(cancertype)
		#listofsigpos=[211,334,88,67,138,25,329,319,218,102,50,323,131]
		#listofsigpos=['915','1378','1423','1518','1603','1610','2083','2871','2875','2881']
		#listofsigpos=[915,1378,1423,1518,1603,1610,2083,2871,2875,2881]
		#[211,334,88,67,138,25,329,319,218,102,50,323,131]
		print('listofsigpos',listofsigpos)
		positiondict=mutperpos(cancertype,listofsigpos)
		print('positiondict',positiondict)
		transdict=defaultdict(list)
		labels=[]
		for aa in mutorder:
			countarray=np.zeros(len(positiondict))
			tempdict={}
			counter=0

			for i in sorted(positiondict):
				labels.append(i)
				for j in positiondict[i]:
					if j==aa:
						countarray[counter]+=1
				counter+=1
			#print('countarray',countarray)
			transdict[aa]=countarray
		print('transdict',transdict)
		X=range(len(listofsigpos))
		colorlist=['black','rosybrown','firebrick','sandybrown','lightgoldenrodyellow','green','skyblue','darkblue','blueviolet','crimson','aqua','palegreen','pink','hotpink','gold','saddlebrown','slategray','darkkhaki','bisque','orchid']
		colorcount=0
		temptuple=[]
		for i in transdict:
			temptuple.append(transdict[i])
			SKCMtuple=np.matrix(temptuple)
		enddone=np.array(np.transpose(SKCMtuple))
		labels=sorted(positiondict.keys())
		fractiontuple=[]
		for i in enddone:
			sumvalue=sum(i)
			temparray=[]
			for j in i:
				temparray.append(j/float(sumvalue))
			fractiontuple.append(sorted(temparray,reverse=True))	
		print('enddone',enddone)
		print('fractiontuple',fractiontuple)
		meanarray=[]
		stdarray=[]
		counter=0
		print('np.transpose(fractiontuple)',fractiontuple[counter],'sum',np.mean(fractiontuple[counter]))
		while counter<len(mutorder):
			'''temparray=np.zeros(len(listofsigpos))
			for i in np.transpose(fractiontuple)[0:counter]:
				temparray=temparray+i
				#print(i)
			print('counter',counter)
			print(temparray)
			meanvalue=np.mean(temparray)
			stdvalue=np.std(temparray)'''
			meanvalue=np.mean(np.transpose(fractiontuple)[counter])
			stdvalue=np.std(np.transpose(fractiontuple)[counter])
			
			
			meanarray.append(meanvalue)
			stdarray.append(stdvalue)
			print(counter,meanvalue,stdvalue,(np.transpose(fractiontuple)[counter]))
			counter+=1
		print(cancertype,'meanarray',meanarray)
		ind=np.arange(len(meanarray))
		plt.bar(ind,meanarray,yerr=stdarray,color="Black",alpha=0.6)
		plt.ylim(0,1)
		plt.ylabel('Fraction of all Mutations')
		plt.xlabel('Number of types of aa substitutions')
		plt.show()

piechart()