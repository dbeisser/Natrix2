#import pandas as pd
#import logging
#import dinopy

#import sys
#this script is to change silva taxonomy to  mothur input taxonomy format
#with open (sys.argv[1]) as f:
with open (snakemake.input[0]) as f: #read file
        tax=f.readlines()[1::]  #remove header

list=[]
for line in (tax): #split line by tab and save in list
        line=line.split("\t")
        list.append(line)


#file=open((str(sys.argv[2])), "w") #new file
file=open(snakemake.output[0], "w") #new file

for i in list:
        file.write(i[0] + "|" + i[15][:-1] + "|" + "size=" + i[1] + "\t" + i[15][:-1] + ";" + "\n") 
# reformat to mothur format and print data in mothur format

file.close
