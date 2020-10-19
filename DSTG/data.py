import os
import numpy as np
import random
import pandas as pd
import time as tm
from operator import itemgetter
from sklearn.model_selection import train_test_split
import pickle as pkl
import scipy.sparse
from metrics import *
from gutils import *
from graph import *


#' data preperation
def input_data(DataDir):
    Link_Graph(outputdir='Infor_Data')
    DataPath1 = '{}/Pseudo_ST1.csv'.format(DataDir)
    DataPath2 = '{}/Real_ST2.csv'.format(DataDir)
    LabelsPath1 = '{}/Pseudo_Label1.csv'.format(DataDir)
    LabelsPath2 = '{}/Real_Label2.csv'.format(DataDir)

    #' read the data
    data1 = pd.read_csv(DataPath1, index_col=0, sep=',')
    data2 = pd.read_csv(DataPath2, index_col=0, sep=',')
    lab_label1 = pd.read_csv(LabelsPath1, header=0, index_col=False, sep=',')
    lab_label2 = pd.read_csv(LabelsPath2, header=0, index_col=False, sep=',')

    lab_data1 = data1.reset_index(drop=True)  #.transpose()
    lab_data2 = data2.reset_index(drop=True)  #.transpose()

    random.seed(123)
    p_data = lab_data1
    p_label = lab_label1

    temD_train, temd_test, temL_train, teml_test = train_test_split(
        p_data, p_label, test_size=0.1, random_state=1)
    temd_train, temd_val, teml_train, teml_val = train_test_split(
        temD_train, temL_train, test_size=0.1, random_state=1)

    print((temd_train.index == teml_train.index).all())
    print((temd_test.index == teml_test.index).all())
    print((temd_val.index == teml_val.index).all())
    data_train = temd_train
    label_train = teml_train
    data_test = temd_test
    label_test = teml_test
    data_val = temd_val
    label_val = teml_val

    data_train1 = data_train
    data_test1 = data_test
    data_val1 = data_val
    label_train1 = label_train
    label_test1 = label_test
    label_val1 = label_val

    train2 = pd.concat([data_train1, lab_data2])
    lab_train2 = pd.concat([label_train1, lab_label2])

    #' save objects

    PIK = "{}/datasets.dat".format(DataDir)
    res = [
        data_train1, data_test1, data_val1, label_train1, label_test1,
        label_val1, lab_data2, lab_label2
    ]

    with open(PIK, "wb") as f:
        pkl.dump(res, f)

    print('load data succesfully....')
