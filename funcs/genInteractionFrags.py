# -*- coding: utf-8 -*-
'''
@Time    : 2025/4/18 20:03
@Author  : WY
@File    : genInteractionFrags.py
'''
import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from reuseCode.utils import mol_with_atomLabel
from reuseCode.calECscore import getECPairs
from reusedCode.RotaryFragMethod import RotaryFrag
from interaction_components.plinteraction import get_interactions
from reusedCode.utils import *

dummy = Chem.MolFromSmiles('*')
Hydrogen = Chem.MolFromSmarts('[H]')
def getRotaryFrags(tab,lig):
    lig_mol = Chem.MolFromMolFile(lig, removeHs=False)
    lig_mol = mol_with_atomLabel(lig_mol)
    fragsList,bondsNum,interval = RotaryFrag(lig_mol)
    print(f'[{tab}] finish split, bonds Num:[{bondsNum}],frag Num:[{len(fragsList)}],calculate Time:{interval}!')
    return fragsList

def genECFrags(workpath,lig,pro_atomTmesh,lig_atomTmesh):
    tab=lig.split('_lig')[0]
    fragsList=getRotaryFrags(tab,workpath+'/'+lig)  # RotaryFrag方法
    frags_dict={}
    frags_data_list=[]
    if len(fragsList) > 0:
        ECpairs=getECPairs(workpath+'/'+lig_atomTmesh,workpath+'/'+pro_atomTmesh)
        for frag in fragsList:
            flag=True
            ECscore=0
            pos_res_list, neg_res_list=[],[]
            for atom in frag.GetAtoms():
                if atom.GetSymbol()!='*':
                    atomType=atom.GetProp('atomLabel')
                    if len(ECpairs[ECpairs['lig_atom']==atomType]) > 0:
                        score=float(ECpairs[ECpairs['lig_atom']==atomType]['ECscore'].values.tolist()[0])
                        ECscore+=score
                        if score>=0:
                            pos_res_list.append(
                                ECpairs[ECpairs['lig_atom']==atomType]['pro_residue'].values.tolist()[0])
                        else:
                            neg_res_list.append(
                                ECpairs[ECpairs['lig_atom']==atomType]['pro_residue'].values.tolist()[0])
                    else:
                        flag=False
                        break
            if flag==True:  # 获取所有原子均有静电打分的片段
                rep_frag=AllChem.ReplaceSubstructs(frag, dummy, Hydrogen, replaceAll=True)[0]
                frags_dict[rep_frag]=[frag, ECscore, list(set(pos_res_list)), list(set(neg_res_list))]
                frags_data_list.append([tab, Chem.MolToSmiles(frag), Chem.MolToSmiles(rep_frag), str(ECscore),listToStr(list(set(pos_res_list))), listToStr(list(set(neg_res_list)))])
    print(f'{tab} split Frags: {len(frags_dict)}')
    if len(frags_data_list)>0:
        frags_df=pd.DataFrame(np.array(frags_data_list),columns=['PDBsource', 'Smiles', 'rmvDumSmiles', 'ECscore', 'posResidues','negResidues'])
        frag_path=os.path.join(tmpdir,f'{tab}_ECfrags.csv')
        frags_df.to_csv(frag_path,index=False)
        return frag_path
    else:
        error_path=os.path.join(tmpdir,f'ECError.txt')
        with open(error_path,'w') as f:
            f.write('无法生成静电互补片段集')

def genHbondFrags(workpath,lig,pro):
    tab=lig.split('_lig')[0]
    ligand=mol_with_atomLabel(Chem.MolFromMolFile(workpath+'/'+lig,removeHs=False))
    protein=Chem.MolFromPDBFile(workpath+'/'+pro,removeHs=False)
    result=get_interactions(protein, ligand)
    interactions=result.interactions
    hbonds_frags_list=[]

    hbonds_ldon=interactions.hbonds_ldon
    hbonds_pdon=interactions.hbonds_pdon
    all_hbonds=hbonds_ldon+hbonds_pdon

    hbonds_df=transformDF(all_hbonds, ['ligAtom', 'proRes'])

    if isinstance(hbonds_df, pd.DataFrame):
        fragsList, bondsNum, interval=RotaryFrag(ligand)  # RotaryFrag方法
        if len(fragsList) > 0:
            for frag in fragsList:
                hbonds_res,hydrophobic_res=[],[]
                rep_frag=AllChem.ReplaceSubstructs(frag, dummy, Hydrogen, replaceAll=True)[0]
                for atom in frag.GetAtoms():
                    if atom.GetSymbol()!='*':
                        atomType=atom.GetProp('atomLabel')
                        if isinstance(hbonds_df, pd.DataFrame) and len(hbonds_df[hbonds_df['ligAtom']==atomType])>0:
                            hbonds_res.append(hbonds_df[hbonds_df['ligAtom']==atomType]['proRes'].tolist()[0])
                if len(hbonds_res) > 0:
                    hbonds_frags_list.append([tab, Chem.MolToSmiles(frag), Chem.MolToSmiles(rep_frag), listToStr(list(set(hbonds_res)))])
    # 生成含氢键作用的片段的csv文件
    if len(hbonds_frags_list)>0:
        hbonds_frags_df=pd.DataFrame(np.array(hbonds_frags_list),columns=['PDBsource', 'Smiles', 'rmvDumSmiles', 'HbondResidues'])
        hbond_frag_path=os.path.join(tmpdir,f'{tab}_Hbondsfrags.csv')
        hbonds_frags_df.to_csv(hbond_frag_path, index=False)
        return hbond_frag_path
    else:
        error_path=os.path.join(tmpdir,f'HbondError.txt')
        with open(error_path,'w') as f:
            f.write('无法生成氢键片段集')

def genHydrophobicFrags(workpath,lig,pro):
    tab=lig.split('_lig')[0]
    ligand=mol_with_atomLabel(Chem.MolFromMolFile(workpath+'/'+lig, removeHs=False))
    protein=Chem.MolFromPDBFile(workpath+'/'+pro, removeHs=False)
    result=get_interactions(protein, ligand)
    interactions=result.interactions
    hydrophobic_frags_list=[]

    hydrophobic_contacts=interactions.hydrophobic_contacts

    hydrophobic_df=transformDF(hydrophobic_contacts,['ligAtom', 'proRes'])

    if isinstance(hydrophobic_df, pd.DataFrame):
        fragsList, bondsNum, interval = RotaryFrag(ligand)  # RotaryFrag方法
        print(f'[{tab}] finish split, bonds Num:[{bondsNum}],frag Num:[{len(fragsList)}],calculate Time:{interval}!')
        if len(fragsList) > 0:
            for frag in fragsList:
                hbonds_res, hydrophobic_res = [], []
                rep_frag = AllChem.ReplaceSubstructs(frag, dummy, Hydrogen, replaceAll=True)[0]
                for atom in frag.GetAtoms():
                    if atom.GetSymbol() != '*':
                        atomType = atom.GetProp('atomLabel')
                        if isinstance(hydrophobic_df, pd.DataFrame) and len(hydrophobic_df[hydrophobic_df['ligAtom'] == atomType]) > 0:
                            hydrophobic_res.append(hydrophobic_df[hydrophobic_df['ligAtom'] == atomType]['proRes'].tolist()[0])
                if len(hydrophobic_res) > 0:
                    hydrophobic_frags_list.append([tab, Chem.MolToSmiles(frag), Chem.MolToSmiles(rep_frag),listToStr(list(set(hydrophobic_res)))])
    # 生成含疏水作用的片段的csv文件
    if len(hydrophobic_frags_list) > 0:
        hydrophobic_frags_df = pd.DataFrame(np.array(hydrophobic_frags_list),columns=['PDBsource', 'Smiles', 'rmvDumSmiles', 'HydrophobicResidues'])
        hydrophobic_frag_path=os.path.join(tmpdir,f'{tab}_Hydrophobicfrags.csv')
        hydrophobic_frags_df.to_csv(hydrophobic_frag_path, index=False)
        return hydrophobic_frag_path
    else:
        error_path=os.path.join(tmpdir,f'HydrophobicError.txt')
        with open(error_path,'w') as f:
            f.write('无法生成疏水片段集')