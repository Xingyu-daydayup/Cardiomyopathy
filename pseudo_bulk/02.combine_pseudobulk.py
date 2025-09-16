import os
import re
import sys


def main():
    path = '/data/lixingyu/condition_analysis/dcm_hcm/vcm/'

    # input = open(path + 'gene2type.txt', 'r')
    output11 = open(path + CELLTYPE + '/combined_pseudobulk.' + CELLTYPE + '.sample.PG.xls', 'w')
    # output12 = open(path + CELLTYPE + '/combined_pseudobulk.info.' + CELLTYPE + '.sample.PG.xls', 'w')

    # output21 = open(path + CELLTYPE + '/combined_pseudobulk.' + CELLTYPE + '.sample.lncRNA.xls', 'w')
    # output22 = open(path + CELLTYPE + '/combined_pseudobulk.info.' + CELLTYPE + '.sample.lncRNA.xls', 'w')

    gene2type = {}
    samples = []
    total_files=os.listdir(path)
    for i in total_files:
        if i.endswith('.csv'):
            samples.append(i[:-14])

    # sampleinfo = {'control': 'control', 'dehp': 'treat'}

    # for line in input:
    #     line = line.rstrip()
    #     parts = line.split()
    #     gene2type[parts[0]] = parts[1]

    alldata = {}
    uniongenes = {}
    print_names = []
    keep_names = []

    for mysample in samples:
        alldata[mysample] = {}
        handle = open(path + CELLTYPE + '/' + mysample + '.pseudobulk.' + CELLTYPE + '.csv', 'r')
        k = 0
        for line in handle:
            line = line.strip()
            parts = line.split()
            if k == 0:
                keepindex = []
                keepnames = []
                for j in range(0, len(parts)):
                    keepindex.append(j + 1)
                    keepnames.append(parts[j])
                    keep_names.append(parts[j])
                print_names.append('\t'.join(keepnames))
            else:
                each = []
                for j in keepindex:
                    each.append(parts[j])
                alldata[mysample][parts[0]] = '\t'.join(each)
                if parts[0] not in uniongenes:
                    uniongenes[parts[0]] = 0
                uniongenes[parts[0]] += 1
            k += 1
        handle.close()
    output11.write('\t'.join(print_names) + '\n')
    # output21.write('\t'.join(print_names) + '\n')

    for mygene in uniongenes:
        if uniongenes[mygene] != len(samples):
            continue
        each = []
        for mysample in samples:
            each.append(alldata[mysample][mygene])

        output11.write(mygene + '\t' + '\t'.join(each) + '\n')

            # if gene2type[mygene] == 'lncRNA':
            #     output21.write(mygene + '\t' + '\t'.join(each) + '\n')

    #     else:
    #         newgenename = mygene.split('.')[0]
    #
    #         if newgenename not in gene2type:
    #             newgenename = re.sub('-', '_', newgenename)
    #             if newgenename not in gene2type:
    #                 continue
    #             if gene2type[newgenename] == 'protein_coding':
    #                 output11.write(mygene + '\t' + '\t'.join(each) + '\n')
    #
    #             if gene2type[newgenename] == 'lncRNA':
    #                 output21.write(mygene + '\t' + '\t'.join(each) + '\n')
    #
    # output12.write('Cell' + '\t' + 'Cell_Cycle' + '\t' + 'Locus' + '\t' + 'Replicate' + '\n')
    # output22.write('Cell' + '\t' + 'Cell_Cycle' + '\t' + 'Locus' + '\t' + 'Replicate' + '\n')


    # count_info = {}
    # for c in keep_names:
    #     myloc = sampleinfo[c.split('_')[-1]]
    #
    #     if myloc not in count_info:
    #         count_info[myloc] = 1
    #     else:
    #         count_info[myloc] += 1
    #     myrep = str(count_info[myloc])
    #
    #
    #     output12.write(c + '\t' + c + '\t' + 'S' + '\t' + myloc + '\t' + myrep + '\n')
    #     output22.write(c + '\t' + c + '\t' + 'S' + '\t' + myloc + '\t' + myrep + '\n')

    # input.close()
    output11.close()
    # output12.close()

    # output21.close()
    # output22.close()


if __name__ == '__main__':
    CELLTYPE = 'vcm'
    main()

