import numpy 
import sklearn
import random
import math
import pylab
import bisect
import operator
from sklearn.metrics import confusion_matrix,roc_auc_score
from sklearn import linear_model,cross_validation,preprocessing
from sklearn.feature_selection import chi2
from sklearn import preprocessing
from sklearn import linear_model
from multiprocessing import Pool, Value
from ctypes import c_int
import time
import timeit
import gc

import sys
import os
#######################################
# Based on OV data Successful results #
#######################################

#last result
#Parallel time =  3615.33759713 (vostro i3)
#(27578, 2)
#[50, 0.71818181818181825]

#with open('/home/debajyoti/Data/gnexp_data.csv','r') as source:
#    data = [ (random.random(), line) for line in source ]
#data.sort()
#with open('/home/debajyoti/Data/gnexp_data.csv_shuffle','w') as target:
#    for _, line in data:
#        target.write( line )
numpy.seterr(all='ignore') 
filename = os.path.abspath('../output/adjacency_2.csv')
ncombi = 2000
print(str(sys.argv))

scale_limit = 100
print "LOADING DATA..."
readdata = numpy.loadtxt(os.path.abspath('../output/train.csv'),delimiter=',',skiprows=1)
#train_data = numpy.transpose(readdata)
#readdata = numpy.loadtxt('test.csv',delimiter=',',skiprows=1)
#test_data = train_data
readdata = numpy.transpose(readdata)
train_data, test_data, y_train, y_test = cross_validation.train_test_split(readdata, readdata[:,0], test_size=0.5, random_state=1)
#numpy.save(os.path.abspath('../output/train_subset.npy'), train_data)
#numpy.save(os.path.abspath('../output/valid_subset.npy'), test_data)
#sys.exit(100)
#data=numpy.loadtxt('/home/debajyoti/Data/gnexp_data.csv_shuffle',delimiter=',')
random.seed(1)
#random.shuffle(data)
#numpy.savetxt('test.out', data, delimiter=',')
#m=data.shape[0]
n=train_data.shape[1]
print "No of features:", n - 1
#c=(80*m)//100
print "No of training samples:", train_data.shape[0]
print "No of combinations per feature:", scale_limit

#y_train=train_data[:,0]
#y_test=test_data[:,0]

# std_scale = preprocessing.MinMaxScaler().fit(train_data)
# train_data = std_scale.transform(train_data)
# test_data = std_scale.transform(test_data)

result = numpy.array([[0,0.]])

final_res = None
counter = Value(c_int)
comb_counter = Value(c_int)

processors = int(sys.argv[1])
begin = 1#int(sys.argv[2])
end = n#int(sys.argv[3])

def adjacency_clust (filename):
 listoflist=[]
 fileclust = open(filename, 'r')
 for line in fileclust.readlines():
  row = line.split(',')
  maxlim = len(row)
  cols = [int(x) for x in row[0:maxlim-1]]
  listoflist.append(cols[0:maxlim])
  fileclust.close()
 return listoflist

def clust_size(clusters):
 length = []
 for k in range(1, len(clusters)):
  length.append(len(clusters[k]))
 return length


def roulette_select(population, fitnesses, num):
 #print population
 #print fitnesses
 #print num
 total_fitness = float(sum(fitnesses))
 rel_fitness = [f/total_fitness for f in fitnesses]
 #print rel_fitness
 # Generate probability intervals for each individual
 probs = [sum(rel_fitness[:i+1]) for i in range(len(rel_fitness))]
 #print probs
 # Draw new population
 new_population = []
 for n in xrange(num):
  r = random.random()
  for (i, individual) in enumerate(population):
   if r <= probs[i]:
    new_population.append(individual)
    break
 return new_population

def combinations(clusters,nos):
 indexset = []
 req = 30#random.randint(len(clusters),2*len(clusters))
 for i in range(ncombi):
  #indexrow = random.sample(range(1,n), random.randint(25,50))	#sample values range, random sample size
  ###################
  #read cluster list
  #read each line to choose in proportion
  
  #print req
  assorted = []
  indexrow = []
  ########
  items =  roulette_select(range(0,len(clusters)),nos, req)
  #print items
  ########
  for k in range(0, req):
   row = clusters[items[k]]
   #row = clusters[k]
   maxlim = len(row)
   #nos = int(math.ceil(float(maxlim)/float(n)*req))
   assorted.append(random.sample(row, 1))
  assorted.sort() #sort list to search exixtance of a gene in each combination
  map(indexrow.extend, assorted)
  indexset.append(indexrow)
 #sys.stdout.write("Training progress: %d%%   \r" % (comb_counter.value/(float)(ncombi)*100) )
 #sys.stdout.flush()
 return indexset
 
def final_logit (indices):
 global result
 #print indices
 x_train=train_data[:,indices]
 x_test=test_data[:,indices]

 # find logistic function
 logistic = linear_model.LogisticRegression(C=1.0)
 probas_=logistic.fit(x_train,y_train)
 C=logistic.coef_
 dlabel=logistic.predict(x_test)
 ptest=logistic.predict_proba(x_test)

 # find confusion matrix
 #mat=confusion_matrix(y_test, dlabel)
 #print 'the confusion matrix is: \n ', mat
 
 #fpr,tpr,threshold=sklearn.metrics.roc_curve(y_test,probas_)
 #auc=sklearn.metrics.auc(fpr, tpr)
 #print ptest
 rauc=roc_auc_score(y_test,dlabel)
 newrow = [rauc,dlabel]
 #print ' func ', newrow
 #result = numpy.vstack((result,newrow)) 
 return rauc

def mylogit (indices):
 global result
 x_train=train_data[:,indices]
 x_test=test_data[:,indices]

 # find logistic function
 logistic = linear_model.LogisticRegression(C=1.0)
 probas_=logistic.fit(x_train,y_train)
 C=logistic.coef_
 dlabel=logistic.predict(x_test)
 ptest=logistic.predict_proba(x_test)

 # find confusion matrix
 #mat=confusion_matrix(y_test, dlabel)
 #print 'the confusion matrix is: \n ', mat
 
 #fpr,tpr,threshold=sklearn.metrics.roc_curve(y_test,probas_)
 #auc=sklearn.metrics.auc(fpr, tpr)
 rauc=roc_auc_score(y_test,ptest[:,1])
 newrow = [len(indices),rauc]
 #print ' func ', newrow
 #result = numpy.vstack((result,newrow)) 
 return newrow

def algorithm(gene):
 counter.value+=1
 
 #print counter.value
 res_inc = []
 res_exc = []
 s=list(range(ncombi))
 random.seed(1)
 random.shuffle(s)
 p_val = 1
 for i in range(scale_limit):
  indexrow = list(indexset[s[i]])
  high = len(indexrow)
  #print gene, i, s[i]
  r = bisect.bisect_left(indexrow,gene,lo=0,hi=high)
  if r< high and indexrow[r] == gene:
   res_inc.append(final_res[s[i],:])
   indexrow.remove(gene)
   auc_t = mylogit(indexrow)
   res_exc.append(auc_t)
   #print "found", "auc=", final_res[s[i],1]
  else:
   res_exc.append(final_res[s[i],:])
   indexrow.append(gene)
   auc_t = mylogit(indexrow)
   res_inc.append(auc_t)
  	#print "not found"
 ## END APPROXIMATE ##
 #####################
 final_inc = numpy.array(res_inc)
 final_exc = numpy.array(res_exc)
 #print final_exc.shape
 #print final_inc[:,1], final_exc[:,1]
 #print len(final_inc[:,1])
 res = list(numpy.array(final_inc[:,1]) - numpy.array(final_exc[:,1]))
 #print res.count(0.0)
 med_inc = numpy.median(final_inc[:,1])
 med_exc = numpy.median(final_exc[:,1])
 p_val = 2.0
 #if (scale_limit - res.count(0.0)) >  20 :
  #print scale_limit - res.count(0.0)
  #if(med_inc>med_exc):
  #t,p_val=stats.ranksums(final_inc[:,1], final_exc[:,1])
  #t,p_val=stats.wilcoxon(final_inc[:,1], final_exc[:,1], correction=True)
 #print 't\tp_val'
 #print t, p_val
 wil_row_p = [gene,p_val,sum(res)]
 wil_row_p.extend(final_inc[:,1])
 wil_row_p.extend(final_exc[:,1])
 #print numpy.median(final_inc[:,1])-numpy.median(final_exc[:,1]),' ',p_val
 #sys.stdout.write("Training progress: %d%%   \r" % (comb_counter.value/(float)(ncombi)*100) )
 sys.stdout.write("Evaluation progress: %d%%	\r" % (counter.value/(float)(end-begin)*100) )
 sys.stdout.flush()
 return wil_row_p


if __name__ == '__main__':
 print '===Parallel==='
 
 #worker pool
 p = processors #processes

 #################
 ## PRE-TRAIN ####
 #################
 print "Reading clusters..."
 clusts = adjacency_clust(filename)
 #nos = meta_genes(clusts)
 nos = clust_size(clusts)
 #remove cluster with max length
 print "TRAINING ", ncombi," combintaions..."
 indexset = combinations(clusts,nos)
 numpy.savetxt(os.path.abspath('../output/indexset.comb'), indexset, delimiter=',')
 pool = Pool(p)
 result = pool.map(mylogit,indexset)
 pool.close()
 pool.join()
 final_res = numpy.array(result)
 numpy.savetxt(os.path.abspath('../output/auc.comb'), final_res[:,1], delimiter=',')
 ##########################################
 
 counter = Value('i', 0)
 tic = timeit.default_timer()
 # for each selected gene 
 p = processors
 geneset=numpy.array(range(begin,end))
 #print geneset
 print "Evaluating features..."
 numpy.savetxt(os.path.abspath('../output/gnsel_id'), geneset, delimiter=',')
 pool = Pool(processes=p)
 wil = pool.map(func=algorithm,iterable=geneset)
 pool.close()
 pool.join()
 #wil = algorithm(218)
 wil_p = numpy.array(wil)
 numpy.savetxt(os.path.abspath('../output/new_wil_p'), wil_p, delimiter=',')
 count = 0

 #wil_p=numpy.array(wil_row)
 toc = timeit.default_timer()
 print 'Parallel time = ', toc - tic
 filepath = os.path.abspath('../output/wil_p_')+str(begin)+'_'+str(end)
 numpy.savetxt(filepath, wil_p, delimiter=',')
 copyfile(filepath,os.path.abspath('../output/bigfile'))
 print wil_p.shape

