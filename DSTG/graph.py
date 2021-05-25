from gutils import *
import glob
import pandas as pd
import os
import itertools


def Link_Graph(outputdir):
    path0 = os.path.join(os.getcwd(), outputdir)

    #' import processed data
    files1 = glob.glob(path0 + "/ST_count/*.csv")
    files1.sort()
    count_list = []
    for df in files1:
        print(df)
        count_list.append(pd.read_csv(df, index_col=0))

    files2 = glob.glob(path0 + "/ST_norm/*.csv")
    files2.sort()
    norm_list = []
    for df in files2:
        print(df)
        norm_list.append(pd.read_csv(df, index_col=0))

    files3 = glob.glob(path0 + "/ST_scale/*.csv")
    files3.sort()
    scale_list = []
    for df in files3:
        print(df)
        scale_list.append(pd.read_csv(df, index_col=0))

    files4 = glob.glob(path0 + "/ST_label/*.csv")
    files4.sort()
    label_list = []
    for df in files4:
        print(df)
        label_list.append(pd.read_csv(df, index_col=0))

    fpath = os.path.join(path0, 'Variable_features.csv')
    features = pd.read_csv(fpath, index_col=False).values.flatten()

    N = len(count_list)
    if (N == 1):
        combine = pd.Series([(0, 0)])
    else:
        combin = list(itertools.product(list(range(N)), list(range(N))))
        index = [i for i, x in enumerate([i[0] < i[1] for i in combin]) if x]
        combine = pd.Series(combin)[index]

    link1 = Link_graph(count_list=count_list,
                       norm_list=norm_list,
                       scale_list=scale_list,
                       features=features,
                       combine=combine)

     #' ---- input data for link grpah 2 -----
    files1 = glob.glob(path0 + "/ST_count/*.csv")
    files1.sort()
    tem_count = pd.read_csv(files1[1], index_col=0)
    tem_count.columns = tem_count.columns.str.replace("mixt_", "rept_")

    files2 = glob.glob(path0 + "/ST_norm/*.csv")
    files2.sort()
    tem_norm = pd.read_csv(files2[1], index_col=0)
    tem_norm.columns = tem_norm.columns.str.replace("mixt_", "rept_")

    files3 = glob.glob(path0 + "/ST_scale/*.csv")
    files3.sort()
    tem_scale = pd.read_csv(files3[1], index_col=0)
    tem_scale.columns = tem_scale.columns.str.replace("mixt_", "rept_")

    files4 = glob.glob(path0 + "/ST_label/*.csv")
    files4.sort()
    tem_label = pd.read_csv(files4[1], index_col=0)
    tem_label.columns = tem_label.columns.str.replace("mixt_", "rept_")

    count_list2 = [count_list[1],tem_count]
    norm_list2 = [norm_list[1], tem_norm]
    scale_list2 = [scale_list[1], tem_scale]

    link2 = Link_graph(count_list=count_list2,
                       norm_list=norm_list2,
                       scale_list=scale_list2,
                       features=features,
                       combine=combine,
                       k_filter=100)

    graph1 = link1[0].iloc[:, 0:2].reset_index()
    graph1 = graph1.iloc[:,1:3] 
    graph1.to_csv('./Datadir/Linked_graph1.csv')

    graph2 = link2[0].iloc[:, 0:2].reset_index()
    graph2 = graph2.iloc[:,1:3]
    graph2.to_csv('./Datadir/Linked_graph2.csv')

    label1 = label_list[0]
    label1.to_csv('./Datadir/Pseudo_Label1.csv', index=False)

    label2 = label_list[1]
    label2.to_csv('./Datadir/Real_Label2.csv', index=False)
