import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import csv
import sys


infile = sys.argv[1]
outfile = sys.argv[2]

matrix = []

subfams = []
positions = []

count = 0
with open(infile,'r') as source:
    for line in source:
    	count+=1
    	
    	if count == 1:
    		positions = line.split('\t')
    		positions = positions[1:-2]
    	else:
    		fields = line.split('\t')
    		subfams.append(fields[0])
    		fields = fields[1:-2]
    		
    		for i in range(len(fields)):
    			fields[i] = float(fields[i])
    		matrix.append(fields)


for j in range(len(positions)):
	max = 0
	for i in range(len(subfams)):
		if(matrix[i][j] > max):
			max = matrix[i][j]
				
	for i in range(len(subfams)):
		if max:
			matrix[i][j] = matrix [i][j] / max
		else:
			#columns started out empty, no skip state either
			matrix[i][j] = matrix[i][j]
			
matrix = np.array(matrix)
		
fig, ax = plt.subplots(figsize = (200, 50))

im = ax.imshow(matrix, interpolation='nearest', aspect='auto')

ax.set_yticks(np.arange(len(subfams)))
ax.set_yticklabels(subfams)

plt.savefig(outfile)

