# -*- coding: utf-8 -*-
'''
@Time    : 2025/4/18 20:07
@Author  : WY
@File    : calECscore.py
'''

import numpy as np
import pandas as pd
from math import exp

def cal_EC_score_6(mol_esp,pro_esp):
    if mol_esp*pro_esp >= 0:
        ## 设定静电碰撞得分在[-1,0.5]之间
        EC_local = 1 - (abs(mol_esp)+abs(pro_esp))/max(abs(mol_esp),abs(pro_esp),7)
        if EC_local >= 0:
            EC_local = 0.5*EC_local
        else:
            EC_local = EC_local
    else:
        if mol_esp + pro_esp >= 0:
            EC_local = 1- (mol_esp + pro_esp)/max(mol_esp,pro_esp,7)
        else:
            EC_local = 1- (mol_esp + pro_esp)/min(mol_esp,pro_esp,-7)
    return EC_local

def cal_EC_score_7(mol_esp,pro_esp):
    if mol_esp*pro_esp >= 0:
        EC_local = 0.5*exp(-(abs(mol_esp)*abs(pro_esp)))
    else:
        EC_local = exp(-1/(abs(mol_esp)+abs(pro_esp)))
    return EC_local

def cal_EC_score_8(mol_esp,pro_esp):
    EC_local = cal_EC_score_7(mol_esp,pro_esp)/2
    return EC_local

def cal_formula(mol_esp, pro_esp):
    max_esp = max(abs(mol_esp), abs(pro_esp))
    min_esp = min(abs(mol_esp), abs(pro_esp))
    if max_esp < 0.1:
        score = cal_EC_score_8(mol_esp, pro_esp)
    else:
        if max_esp <= 1:
            score = cal_EC_score_7(mol_esp, pro_esp)
        else:
            score = cal_EC_score_6(mol_esp, pro_esp)
    return score

def pro_tra(dis_list):
    lis1 = []
    min_dis = min(dis_list)
    for j, distance in enumerate(dis_list):
        if distance == min_dis:
            lis1.append([j, float(min_dis)])  # 添加最短距离对应下的蛋白顶点索引值
    return lis1[-1]  # 返回目标蛋白顶点的最短距离索引值、最短距离

def calDistance(lig_item, pro_df):
    dis_df = np.sqrt((pro_df[:, 0] - lig_item[0]) ** 2 + (pro_df[:, 1] - lig_item[1]) ** 2 + (pro_df[:, 2] - lig_item[2]) ** 2)
    return dis_df.tolist()

def func(lig_df, pro_df, pro_coors, lig_coors):
    all_result = []
    for i, lig_coor in enumerate(lig_coors):
        distance_list = calDistance(lig_coor, pro_coors)
        lig_info = lig_df.loc[i,:]
        pro_dis_list = pro_tra(distance_list)
        pro_info = pro_df.loc[pro_dis_list[0],:]
        pro_residue = pro_info['residueName']+str(pro_info['residueNumber'])
        score = cal_formula(lig_info['pqr_esp'], pro_info['pqr_esp'])
        all_result.append([lig_info['own_atom'], pro_info['own_atom'],pro_residue,score])
    result_df=pd.DataFrame(np.array(all_result),columns=['lig_atom','pro_atom','pro_residue','ECscore'])
    return result_df

def getECPairs(lig_file,pro_file):
    lig_df = pd.read_csv(lig_file)
    pro_df = pd.read_csv(pro_file)
    pro_coors = np.array(pro_df[['X', 'Y', 'Z']], dtype=float)
    lig_coors = np.array(lig_df[['X', 'Y', 'Z']], dtype=float)

    pairECscore = func(lig_df, pro_df, pro_coors, lig_coors)
    lig_atoms = set(pairECscore['lig_atom'].values.tolist())
    res_list = []
    for lig_atom in lig_atoms:
        raw_ECscore = pairECscore[pairECscore['lig_atom'] == lig_atom]
        mid_ECscore = pd.DataFrame(raw_ECscore['pro_residue'].value_counts())
        final_ECscore = raw_ECscore[raw_ECscore['pro_residue'] == list(mid_ECscore.index)[0]].reset_index(drop=True)
        score = np.mean(list(map(float,final_ECscore['ECscore'].values.tolist())))
        res_list.append([final_ECscore.loc[0, 'lig_atom'], final_ECscore.loc[0, 'pro_residue'], score])
    res_df = pd.DataFrame(np.array(res_list), columns=['lig_atom', 'pro_residue', 'ECscore'])
    res_df.sort_values(by='ECscore', inplace=True, ascending=False)
    return res_df
