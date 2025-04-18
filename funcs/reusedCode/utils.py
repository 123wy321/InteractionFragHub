# -*- coding: utf-8 -*-
'''
@Time    : 2025/4/18 20:05
@Author  : WY
@File    : utils.py
'''
import numpy as np
import pandas as pd
from rdkit import Chem


def mol_with_atomLabel(mol):
    atom_dict = {}
    for atom in mol.GetAtoms():
        atomType = atom.GetSymbol().upper()
        if atomType not in atom_dict.keys():
            atom_dict[atomType] = 1
            atom.SetProp('atomLabel', atomType + str(atom_dict[atomType]))
        else:
            atom_dict[atomType] += 1
            atom.SetProp('atomLabel', atomType + str(atom_dict[atomType]))
    return mol


def transformDF(data, columnsList):
    if len(data) > 0:
        data_dict = {}
        for item in data:
            if item.lig_atom_name not in data_dict.keys():
                data_dict[item.lig_atom_name] = [item.pro_res_name]
            else:
                data_dict[item.lig_atom_name].append(item.pro_res_name)
        data_list = []
        for key, value in data_dict.items():
            data_list.append([key, listToStr(list(set(value)))])
        data_df = pd.DataFrame(np.array(data_list), columns=columnsList)
        return data_df


def listToStr(list):
    res_str = ''
    if len(list) > 0:
        for index, item in enumerate(list):
            res_str += item
            if index < len(list) - 1:
                res_str += ','
    return res_str

def delete_attachment_H(smi):
    m=Chem.MolFromSmiles(smi)
    if m is not None:
        mw=Chem.RWMol(m)
        while Chem.MolToSmiles(mw.GetMol()).count('*')>0:
            for atom in mw.GetMol().GetAtoms():
                if atom.GetSymbol()=='*':
                    mw.ReplaceAtom(atom.GetIdx(),Chem.Atom(1))
                    break
        mw=Chem.RemoveHs(mw)
        return Chem.MolToSmiles(mw)
    else:
        return '[Xe]'

def delete_attachment(smi):
    m = Chem.MolFromSmiles(smi)
    if m is not None:
        mw = Chem.RWMol(m)
        atta_idx = []
        for i in range(0,mw.GetNumAtoms()):
            if mw.GetAtomWithIdx(i).GetSymbol() == '*':
                atta_idx.append(i)
        for i in reversed(atta_idx):
            mw.RemoveAtom(i)
        return Chem.MolToSmiles(mw)
    else:
        return '[Xe]'

def Conn(data1,data2,attachment): #smiles列表
    connected = []
    for i in range(len(data1)):
        if i == 0:
            for s2 in data2:
                connected.extend(Mol_Conn(data1[i],s2,attachment))
        else:
            len_tmp = len(connected)
            for j in range(len(connected)):
                if Chem.MolFromSmiles(connected[j]) is not None:
                    connected.extend(Mol_Conn(connected[j],data1[i],attachment))
            del connected[:len_tmp]
    return connected

def Mol_Conn(s1,s2,attachment): #需连接的片段的smiles,取代位点的表示
    conn = [] #有多个取代位点
    m1 = Chem.MolFromSmiles(s1)
    m2 = Chem.MolFromSmiles(s2)
    m = Chem.CombineMols(m1, m2)
    mw = Chem.RWMol(m)
    for i in range(0,m1.GetNumAtoms()):
        if m1.GetAtomWithIdx(i).GetSymbol() == attachment:  #判断是否是取代位点
            neighbor1_idx = m1.GetAtomWithIdx(i).GetNeighbors()[0].GetIdx() #取代位点邻位原子在分子中的序号
            for j in range(0,m2.GetNumAtoms()):
                if m2.GetAtomWithIdx(j).GetSymbol() == attachment:
                    neighbor2_idx = m2.GetAtomWithIdx(j).GetNeighbors()[0].GetIdx()
                    mw.AddBond(neighbor1_idx,neighbor2_idx+m1.GetNumAtoms(),Chem.BondType.SINGLE)
                    mw.RemoveAtom(j+m1.GetNumAtoms())
                    mw.RemoveAtom(i)
                    conn.append(Chem.MolToSmiles(mw)) #一个位点取代结束，进行存储
                    mw = Chem.RWMol(m) #刷新mw,保证下一个位点取代序号正确
    return conn