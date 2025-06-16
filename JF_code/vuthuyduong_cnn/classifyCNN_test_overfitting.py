#!/usr/bin/env python
# FILE: classifyCNN.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019

# cited from https://github.com/vuthuyduong/fungiclassifiers/blob/master/models/CNN/classifyCNN.py


#from sklearn.model_selection import StratifiedKFold
#from sklearn.preprocessing import StandardScaler
import sys
if sys.version_info[0] >= 3:
	unicode = str
import os, argparse
#from keras.models import Sequential
#from keras.layers import Dense
#from keras.layers import Convolution1D
#from keras.datasets import mnist
#from keras.layers import Dropout, Activation, Flatten
#from keras.layers import Convolution1D, MaxPooling1D
#from keras.utils import np_utils
#from keras import backend as K
from keras.models import load_model
import numpy as np
import json
from Bio import SeqIO
#import random
import multiprocessing
parser=argparse.ArgumentParser(prog='classifyCNN.py',  
							   usage="%(prog)s [options] -i fastafile -c classifiername -mp minproba -mc mincoverage -j variationjsonfilename",
							   description='''Script that classifies the sequences of the fasta files using CNN model. The classified sequences with a probability less than minproba will be verified using BLAST. The mincoverage is given for BLAST comparision. The json file is the the variation of the sequences within each group of the training dataset. This file is created during the training of the model, and used optionally for the verification of the classification. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-o','--out', required=False, help='The folder name containing the model and associated files.') #optional
parser.add_argument('-c','--classifier', required=False, help='the folder containing the classifier.')

args=parser.parse_args()
testfastafilename= args.input
modelname=args.classifier
reportfilename=args.out
out_path = "/output/folder/"
nproc=multiprocessing.cpu_count()

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def LoadConfig(modelname):
	if modelname[len(modelname)-1]=="/":
		modelname=modelname[:-1]
	basename=modelname
	if "/" in modelname:
		basename=modelname[modelname.rindex("/")+1:]
	configfilename=modelname + "/" + basename + ".config"
	classifiername=""
	jsonfilename=""
	classificationfilename=""
	classificationpos=0
	data_max=0
	k=6
	configfilename = "config_file.config"
	configfile=open(configfilename)
	for line in configfile:
		texts=line.split(": ")
		if texts[0]=="Classifier name":
			classifiername=texts[1].rstrip()
		if texts[0]=="K-mer number":
			k=int(texts[1].rstrip())
		if texts[0]=="Classification filename":
			classificationfilename=texts[1].rstrip()
		if texts[0]=="Data max":
			data_max=float(texts[1].rstrip())
		if texts[0]=="Column number to be classified":
			classificationpos=int(texts[1].rstrip())
		if texts[0]=="Classes filename":
			jsonfilename=texts[1].rstrip()
	return jsonfilename,classifiername,classificationfilename,classificationpos,k,data_max

def loadClassification(classificationfilename,classificationpos):
	#load classification
	allseqids=[]
	classificationfile= open(classificationfilename)
	classifications=[]
	for line in classificationfile:
		texts=line.split("\t")
		if line.startswith("#"):
			continue 
		seqid=texts[0].replace(">","").rstrip()
		classname=""
		if classificationpos < len(texts):
			classname=texts[classificationpos].rstrip()
		if classname !="":
			allseqids.append(seqid)
			classifications.append(classname)
	classificationfile.close()
	return classifications,allseqids

def loadData(matrixfilename,data_max,classifications,allseqids):
	#load vectors from the matrix
	testvectors= list(open(matrixfilename, "r"))
	testvectors=testvectors[1:]
	X=[]
	testseqIDs=[]
	testtaxa=[]
	for vector in testvectors:
		elements=vector.split(",")
		seqid=elements[0]
		testseqIDs.append(seqid)
		taxonname=""
		if seqid in allseqids:
			taxonname= classifications[allseqids.index(seqid)]
		testtaxa.append(taxonname)
		X.append(elements[1:])
	X=np.array(X,dtype=float)
	if data_max >0:
		X = X/data_max
	return testseqIDs,X,testtaxa

def SavePrediction(classdict,testseqIDs,testtaxa,pred_labels,probas,outputname):
	output=open(outputname,"w")
	output.write("#SequenceID\tGiven label\tPrediction\tFull classification\tProbability\n")
	i=0
	f2 = open("species_list.txt", "r")
	lines = f2.readlines()
	f2.close()
	line = lines[0]
	keys = line.split("\t")[:-1]
# 	keys=list(classdict.keys())
	for seqid in testseqIDs:
		proba =probas[i][pred_labels[i]]
		giventaxonname=testtaxa[i]
		predictedname =keys[pred_labels[i]]
		classification=classdict[predictedname]['classification']
# 		output.write((seqid + "\t" + giventaxonname + "\t"  + "pridicted" + "\t"+ classification + "\t" + str(proba) + "\n").encode('ascii', 'ignore')	)
		output.write(seqid + "\t" + giventaxonname + "\t"  + predictedname + "\t"+ classification + "\t" + str(proba) + "\n")
		i=i+1
	output.close()
	
	
	
	
	
def GetTaxonomicClassification(classificationpos,header,texts):
	classification=""
	p_s=len(texts)
	p_g=len(texts)
	p_f=len(texts)
	p_o=len(texts)
	p_c=len(texts)
	p_p=len(texts)
	p_k=len(texts)
	i=0
	for text in header.split("\t"):
		text=text.rstrip()
		if text.lower()=="species":
			p_s=i
		elif text.lower()=="genus":
			p_g=i	
		elif text.lower()=="family":
			p_f=i	
		elif text.lower()=="order":
			p_o=i	
		elif text.lower()=="class":
			p_c=i	
		elif text.lower()=="phylum":
			p_p=i	
		elif text.lower()=="kingdom":
			p_k=i	
		i=i+1 
	species="s__"
	genus="g__"
	family="f__"
	order="o__" 
	bioclass="c__"
	phylum="p__"
	kingdom="k__"
	if p_s< len(texts):
		species="s__" + texts[p_s].rstrip()
	if p_g< len(texts):
		genus="g__" + texts[p_g].rstrip()	
	if p_f< len(texts):
		family="f__" + texts[p_f].rstrip()	
	if p_o< len(texts):
		order="o__" + texts[p_o].rstrip()
	if p_c< len(texts):
		bioclass="c__" + texts[p_c].rstrip()	
	if p_p< len(texts):
		phylum="p__" + texts[p_p].rstrip()
	if p_k< len(texts):
		kingdom="k__" + texts[p_k].rstrip()	
	if classificationpos==p_s:
		classification=kingdom +";"+phylum +";"+bioclass +";"+order+";"+family + ";"+ genus+";"+species
	elif classificationpos==p_g:
		classification=kingdom +";"+phylum +";"+bioclass +";"+order+";"+family + ";"+ genus
	elif classificationpos==p_f:
		classification=kingdom +";"+phylum + ";"+bioclass +";"+order+";"+family 
	elif classificationpos==p_o:
		classification=kingdom +";"+phylum + ";"+bioclass + ";"+order
	elif classificationpos==p_c:
		classification=kingdom +";"+phylum + ";"+bioclass 
	elif classificationpos==p_p:
		classification=kingdom +";"+phylum
	elif classificationpos==p_k:
		classification=kingdom
	else:
		classification=texts[classificationpos]
	return classification
	
	
def load_class_dict(fastafilename,classificationfilename,classificationpos):
	#load seqrecords
	seqids=[]
	seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	#load classification
	classificationfile = open(classificationfilename, "r")
	classnames=[]
	level=""
	newseqrecords=[]
	classdict={}
	header=next(classificationfile)
	texts=header.split("\t")
	if classificationpos < len(texts):
		level=texts[classificationpos].rstrip()
	for line in classificationfile:
		texts=line.split("\t")
		seqid=texts[0].replace(">","").rstrip()
		if not seqid in seqrecords.keys():
			continue
		classname=""
		if classificationpos < len(texts):
			classname=texts[classificationpos].rstrip()
			#classname=unicode(classname,errors='ignore')
		if classname !="":
			newseqrecords.append(seqrecords[seqid])
			seqids.append(seqid)
			classnames.append(classname)
			classification=GetTaxonomicClassification(classificationpos,header,texts)
			if classname in classdict.keys():
				if len(classification) > len(classdict[classname]['classification']):
# 				if len(classification) > classdict[classname]['classification']:
					classdict[classname]['classification']=classification
				classdict[classname]['seqids'].append(seqid)
			else:	
				classdict.setdefault(classname,{})
				classdict[classname]['classification']=classification
				classdict[classname]['seqids']=[seqid]
	classificationfile.close()
	basename=GetBase(fastafilename)
	if "/" in basename:
		basename=basename[basename.rindex("/")+1:]
	
	return classdict
	
	
	
if __name__ == "__main__":
	path=sys.argv[0]
	path=path[:-(len(path)-path.rindex("/")-1)]
	#load config of the model
	classesfilename,classifiername,classificationfilename,classificationpos,k,data_max=LoadConfig(modelname)
	print(classesfilename)
# 	#load ref class dict
# 	classdict={}
# 	with open(classesfilename) as classesfile:
# 		classdict = json.load(classesfile)
# 	#load classification
	with open("species.json", 'r') as fp:
		classdict = json.load(fp)

	classifications,allseqids=loadClassification(classificationfilename,classificationpos)	
	#represent sequences of the test dataset as k-mer vector
	testfilename=GetBase(testfastafilename)
	matrixfilename=testfilename + "." + str(k) + ".matrix"
	command="python " + path + "fasta2matrix.py " +  str(k) + " " + testfastafilename + " " + matrixfilename
	print(command)
	# os.system(command)
	testseqIDs,testinputs,testtaxa=loadData(matrixfilename,data_max,classifications,allseqids)
	#load model
# 	new_model = keras.models.load_model(classesfilename)
	testinputs = testinputs.reshape(testinputs.shape + (1,))
	for i in range(20):
		model = load_model("/models/folder/" + "model_" + str(i) + ".keras")
		#predict labels for test dataset
		pred_labels = model.predict_classes(testinputs,verbose = 0)
		probas = model.predict_proba(testinputs)
		#save prediction
		basename=GetBase(classifiername)
		if "/" in basename:
			basename=basename[basename.rindex("/")+1:]
		if reportfilename==None or reportfilename=="":	
			reportfilename=GetBase(testfastafilename) +  "." + basename + "_" + str(i) + ".classified"
		reportfilename=GetBase(testfastafilename) +  "." + basename + "_" + str(i) + ".classified"
		SavePrediction(classdict,testseqIDs,testtaxa,pred_labels,probas,reportfilename)
# 		SavePrediction(None,testseqIDs,testtaxa,pred_labels,probas,reportfilename) # my code
		print("The result is saved in the file: " + reportfilename)
