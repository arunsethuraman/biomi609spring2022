#This program reads a FASTQ file, computes the Q value distribution
#Thereon summarizes it as a boxplot
import array
import itertools
import matplotlib
import matplotlib.pyplot as plt
import statistics
from statistics import mean
from statistics import median
import numpy as np


numreads=0 #counter variable for total number of reads
readlength=0 #counter variable for length of each read
#Note I'm assuming here that all reads are of equal length

#Counter for number of total reads in the file, calculating length of reads
with open("assignment1.fastq","r") as f:
    for line in f.read().split("\n")[3::4]:
        numreads=numreads+1
        readlength=len(line)

print(readlength)
print(numreads)

qvalarray=[]

x=0
with open("assignment1.fastq","r") as f:
    for line in f.read().split("\n")[3::4]:
        z=[]
        for base in range(len(line)-1):
            qval=ord(line[base])-33
            z.append(float(qval))
        qvalarray.append(z)
        x=x+1

print(qvalarray)

q=np.array(qvalarray)
meanqvals=[]
medianqvals=[]
#Compute means across each column; i.e. across each nucleotide location
for x in range(readlength-1):
	meanqvals.append(mean(q[:,x]))
	medianqvals.append(median(q[:,x]))
#plt.boxplot(qvalarray)
#plt.show()
#plt.set_title('Q-value distribution')
fig, ax = plt.subplots(figsize=(5,5))

ax.set_title('Q-value distribution')
ax.plot(meanqvals,'C1',label='Mean')
ax.plot(medianqvals,'C2',label='Median')
ax.set_xlabel('Nucleotide position')
ax.set_ylabel('Q score')
ax.legend();

#plt.show()
plt.savefig('qval.png')                  
