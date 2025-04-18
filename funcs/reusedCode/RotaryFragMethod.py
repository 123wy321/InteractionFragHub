# -*- coding: utf-8 -*-
'''
@Time    : 2025/4/18 19:49
@Author  : WY
@File    : RotaryFragMethod.py
'''
import time
import random
random.seed(42)
from rdkit import Chem
from rdkit.Chem import Lipinski
from itertools import combinations

def RotaryFrag(mol):
    t1 = time.time()
    single_bond_list = []
    double_bond_list = []
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        for neighbor in atom.GetNeighbors():
            j = neighbor.GetIdx()
            bond = mol.GetBondBetweenAtoms(i, j)
            bondType = str(bond.GetBondType())
            atom_neighbor=[an.GetSymbol() for an in atom.GetNeighbors() if an.GetSymbol() != 'H']
            neig_neighbor=[nn.GetSymbol() for nn in neighbor.GetNeighbors() if nn.GetSymbol() != 'H']
            if len(atom_neighbor)>1 or len(neig_neighbor)>1:
                if (j,i) not in single_bond_list and bondType == 'SINGLE' and not bond.IsInRing() and atom.GetSymbol() != 'H' and neighbor.GetSymbol() != 'H':  # 增加不在环上的且末端连接不是H的单键
                    single_bond_list.append((i,j))
                if (j,i) not in double_bond_list and bondType == 'DOUBLE' and not bond.IsInRing() and atom.GetSymbol() == 'C' and neighbor.GetSymbol() == 'C':  # 增加碳碳双键
                    double_bond_list.append((i,j))
    #筛选酰胺键
    react_bond_list=[]
    pop_list=[]
    for index,sbond in enumerate(single_bond_list):
        atom1,atom2=mol.GetAtomWithIdx(sbond[0]),mol.GetAtomWithIdx(sbond[1])
        if atom1.GetSymbol()=='C' and atom2.GetSymbol()=='N':
            for neighbor1 in atom1.GetNeighbors():
                j1 = neighbor1.GetIdx()
                bond11 = str(mol.GetBondBetweenAtoms(sbond[0], j1).GetBondType())
                if neighbor1.GetSymbol()=='O' and bond11=='DOUBLE':
                    pop_list.append(index)
                    react_bond_list.append(sbond)
        elif atom2.GetSymbol()=='C' and atom1.GetSymbol()=='N':
            for neighbor2 in atom2.GetNeighbors():
                j2 = neighbor2.GetIdx()
                bond12 = str(mol.GetBondBetweenAtoms(sbond[1], j2).GetBondType())
                if neighbor2.GetSymbol()=='O' and bond12=='DOUBLE':
                    pop_list.append(index)
                    react_bond_list.append(sbond)
    mid_single_bond_list=[bond for index, bond in enumerate(single_bond_list) if index not in pop_list]
    new_react_bond_list=double_bond_list+react_bond_list
    new_bond_list=[]

    #判断断键数量上限是否超过17
    if len(mid_single_bond_list)+len(new_react_bond_list)>17 and len(new_react_bond_list)<=17:
        new_single_bond_list=random.sample(mid_single_bond_list,17-len(new_react_bond_list))
        new_bond_list = new_single_bond_list + new_react_bond_list
    elif len(mid_single_bond_list)+len(new_react_bond_list)>17 and len(new_react_bond_list)>17:
        new_bond_list=random.sample(mid_single_bond_list+new_react_bond_list,17)
    elif len(mid_single_bond_list)+len(new_react_bond_list)<=17:
        new_single_bond_list=mid_single_bond_list
        new_bond_list=new_single_bond_list+new_react_bond_list
    bond_index_list=[mol.GetBondBetweenAtoms(bond[0], bond[1]).GetIdx() for bond in new_bond_list]

    blocks=[]
    blocks_smi=[]
    if len(bond_index_list)>0:
        all_index_combs = splitBonds(bond_index_list)
        for index_combs in all_index_combs:
            for index_comb in index_combs:
                frags=Chem.FragmentOnBonds(mol,index_comb)
                frags_list=Chem.GetMolFrags(frags,asMols=True)
                for frag in frags_list:
                    mol_heavyAtomNum = Lipinski.HeavyAtomCount(frag)  #计算重原子数
                    frag_smi=Chem.MolToSmiles(frag)
                    if (Chem.RemoveHs(frag).GetNumAtoms()-frag_smi.count('*')>1) and mol_heavyAtomNum<=20 and (frag_smi not in blocks_smi):  # 去除只有一个重原子且重原子数大于20的片段
                        blocks_smi.append(frag_smi)
                        blocks.append(frag)
    t2=time.time()
    interval=t2-t1
    return blocks,len(bond_index_list),interval

def splitBonds(bonds_index_list):
    all_index_combinations=[]
    for i in range(1,len(bonds_index_list)+1):
        indexs=[]
        for j in combinations(bonds_index_list,i):
            indexs.append(list(j))
        all_index_combinations.append(indexs)
    return all_index_combinations

