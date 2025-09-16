import sys
import re
from optparse import OptionParser
import subprocess
import os
import time
import gzip
import math
import numpy as np
import random

class generate_bwamapjobs:

    def __init__(self, path='/data/lixingyu/condition_analysis/dcm_hcm/vcm/'):
        self.path = path
        self.name = []
        total_files=os.listdir(path)
        for i in total_files:
            if i.endswith('rawcounts.csv'):
                self.name.append(i[:-14])
        self.input1 = [open(path + CELLTYPE + '/' + i + '_cluster_info.'+CELLTYPE+'.txt', 'r') for i in self.name]
        self.input2 = [open(path +'/'+ i + '_rawcounts' + '.csv', 'r') for i in self.name]
        self.output = [open(path + CELLTYPE + '/' + i + '.pseudobulk.'+CELLTYPE+'.csv', 'w') for i in self.name]


    def generate_bwamapjobs(self):
        for file_index in range(len(self.name)):
            alltypes = {}
            myconvert = {}
            allsamples = []
            for line in self.input1[file_index]:
                line = line.rstrip()
                parts = line.split('\t')
                if parts[1] == CELLTYPE:
                    alltypes[parts[1]] = 0
                    myconvert[parts[0]] = parts[1]

                    allsamples.append(parts[0])

            print(f'{CELLTYPE} number cur_sample: ',len(allsamples))
            keepnames = []
            for c in alltypes:
                keepnames.append(c)
            keepnames.sort()
            print_names = []
            for m in keepnames:
                print_names.append(str(m) + '_' + self.name[file_index])


            if len(allsamples)<400:
                print('\t'.join(print_names), file=self.output[file_index])
                i = 0
                for line in self.input2[file_index]:
                    line = line.rstrip()
                    parts = line.split(',')

                    if i == 0:
                        allindex = []

                        for mysample in allsamples:
                            for j in range(1,len(parts)):
                                myname = parts[j]
                                if myname == mysample:
                                    allindex.append(j)


                    else:
                        mysum = 0
                        for k in allindex:
                            mysum += int(float(parts[k]))
                        print(parts[0] + '\t' + str(mysum),file=self.output[file_index])
                    i += 1
            else:
                allsamples_new = [0] * 1

                for i in range(0, 1):
                    allsamples_new[i] = []


                    allsamples_new[i]=random.sample(allsamples,k=400)

                keepnames = []
                for c in alltypes:
                    keepnames.append(c)

                keepnames.sort()

                print_names = []

                for m in keepnames:
                    for i in range(0, 1):
                        print_names.append(str(m) + '_' + self.name[file_index])
                print('\t'.join(print_names), file=self.output[file_index])
                # print >> self.output, '\t'.join(print_names)

                i = 0
                for line in self.input2[file_index]:
                    line = line.rstrip()
                    parts = line.split(',')

                    if i == 0:
                        allindex = [0] * 1
                        # print parts[1],"parts[1]"

                        for m in range(0, 1):
                            allindex[m] = []

                            for mysample in allsamples_new[m]:

                                for j in range(1,len(parts)):
                                    myname = parts[j]
                                    if myname == mysample:
                                        allindex[m].append(j)


                    else:

                        print_vals = []
                        for m in range(0, 1):
                            mysum = 0
                            for k in allindex[m]:
                                mysum += int(float(parts[k]))

                            print_vals.append(str(mysum))
                        print(parts[0]+ '\t' + '\t'.join(print_vals), file=self.output[file_index])
                        # print >> self.output, parts[0].split('"')[1] + '\t' + '\t'.join(print_vals)

                    i += 1


if __name__ == '__main__':
    CELLTYPE = 'vcm'
    model = generate_bwamapjobs()
    model.generate_bwamapjobs()
