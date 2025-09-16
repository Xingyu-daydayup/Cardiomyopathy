import sys
import re
from optparse import OptionParser
import subprocess
import os
import time


# /lustre2/home/tangqianzi/software/mirnylab-hiclib-8c146707b6f6/examples/pipeline_2017/
# 5% and 0.5

class generate:

    def __init__(self, path='/data/lixingyu/condition_analysis/dcm_hcm/vcm/'):
        self.path = path
        self.input = open(self.path+ CELLTYPE + '/combined_pseudobulk.'+ CELLTYPE +'.sample.'+ GENETYPE + '.xls', 'r')
        self.input2 = open(self.path+'mito.txt', 'r')
        self.output = open(self.path + CELLTYPE + '/combined_pseudobulk.'+  CELLTYPE +'.sample.filtered.'+ GENETYPE + '.xls', 'w')

    def generate(self):

        mitogenes = []
        for line in self.input2:
            line = line.rstrip()
            mitogenes.append(line)

        i = 0
        count = 0
        for line in self.input:
            line = line.rstrip()
            parts = line.split('\t')
            if i == 0:
                print(line, file=self.output)

            else:
                count_5per = 0
                allvals = []
                for j in range(1, len(parts)):
                    if float(parts[j]) > 0:
                        count_5per += 1
                    allvals.append(float(parts[j]))

                mymean = sum(allvals) / float(len(allvals))

                if float(count_5per) / float(len(parts) - 1) > 0.05 and mymean > 0.1 and (parts[0] not in mitogenes):
                    print(line,file= self.output)
                    # print >> self.output, line
                count += 1

            i += 1

        print(count)

        self.input.close()
        self.input2.close()
        self.output.close()



if __name__ == '__main__':
    CELLTYPE = 'vcm'
    GENETYPE = 'PG'
    model = generate()
    model.generate()

