import csv
import numpy as np

distance = [-25.0,-20.0,-15.0,-10.0,-5.0,0.0,5.0,10.0,15.0,20.0,25.0]
value = 900.0
with open('out_adjusted_LSF_0001.csv','w', newline='') as file:
    writer = csv.writer(file, delimiter=',')
    writer.writerow(['distance','flux_aggregate'])
    for i in range(len(distance)):
        writer.writerow([distance[i],value])
