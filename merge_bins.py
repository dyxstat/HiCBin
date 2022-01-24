#!/usr/bin/env python3

import os
import io
import numpy as np
import copy
import argparse
import warnings
import numpy as npi
import pandas as pd
from sklearn import metrics
from Bio import SeqIO


warnings.filterwarnings("ignore")

def get_no_hidden_folder_list(wd):
    folder_list = []
    for each_folder in os.listdir(wd):
        if not each_folder.startswith('.'):
            folder_list.append(each_folder)

    folder_list_sorte = sorted(folder_list)
    return folder_list_sorte

######Input TAXAassign file#############
def main(pwd_folder , contam_file):
    bin_path = pwd_folder + '/BIN'
    sub_path = pwd_folder + '/SUB_BIN'
    mk = 'mkdir' + ' ' + pwd_folder + '/FINAL_BIN'
    os.system(mk)
    output_dir = pwd_folder + '/FINAL_BIN'

    bin_list = get_no_hidden_folder_list(bin_path)
    sub_list = get_no_hidden_folder_list(sub_path)

    cm = pd.read_csv(contam_file , sep = ',' , header=None)
    cm = cm.values[:,0]
    cm = list(cm)
    skip = 0
    global_index = 1
    for index1 in range(len(bin_list)):
        if index1 < 10:
            bin_name = 'BIN'+ '000' + str(index1)
        elif index1 >= 10 and index1 < 100:
            bin_name = 'BIN'+ '00' + str(index1)
        elif index1 >= 100 and index1 < 1000:
            bin_name = 'BIN'+ '0' + str(index1)
        else:
            bin_name = 'BIN'+str(index1)

        if bin_name in cm:
            skip += 1
            continue
        else:
            bin = bin_path  + '/' +bin_name + '.fa'
            output = output_dir + '/BIN' + str(global_index) + '.fa'
            mv = 'cp' + ' ' + bin + ' ' + output
            global_index += 1
            os.system(mv)
    for index2 in range(len(sub_list)):
        if index2 < 10:
            bin_name = 'SUB'+ '000' + str(index2)
        elif index2 >= 10 and index2 < 100:
            bin_name = 'SUB'+ '00' + str(index2)
        elif index2 >= 100 and index2 < 1000:
            bin_name = 'SUB'+ '0' + str(index2)
        else:
            bin_name = 'SUB'+str(index2)

        bin = sub_path + '/' + bin_name + '.fa'
        output = output_dir + '/BIN' + str(global_index) + '.fa'
        mv = 'cp' + ' ' + bin + ' ' + output
        global_index += 1
        os.system(mv)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("CHECKM",help="Contaminated_Bins_file (csv file)")
    parser.add_argument("OUTDIR",help="Output directory of HiCBin")
    args=parser.parse_args()
    main(args.OUTDIR , args.CHECKM)

        




