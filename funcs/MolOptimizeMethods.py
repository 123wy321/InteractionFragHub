# -*- coding: utf-8 -*-
'''
@Time    : 2025/4/18 20:32
@Author  : WY
@File    : MolOptimizeMethods.py
'''
import os
from random import sample

import pandas as pd
from rdkit import Chem
from rdkit.Chem import QED, Draw

from reusedCode.utils import delete_attachment, Conn, Mol_Conn, delete_attachment_H
from reusedCode.highlightMolReplace import highlightMolsReplace
def fragmention(mol,hit):
    '''
    input -> 完整分子mol , 需替换子结构对应的索引列表hit
    output -> 根据hit打碎mol,返回除了hit以外的碎片
    '''
    mw = Chem.RWMol(mol)
    link = []
    for h in hit:
        for a in mol.GetAtomWithIdx(h).GetNeighbors():
            if a.GetIdx() not in hit:
                link.append([h,a.GetIdx()]) #link[i][1]为连接位点
    #在连接处打碎每一个分子
    for l in link:
        atom1 = l[0]
        atom2 = l[1]
        mw.RemoveBond(atom1,atom2)
        mw.AddAtom(Chem.Atom(0))
        mw.AddBond(atom1,mw.GetNumAtoms()-1,Chem.BondType.SINGLE)
        mw.AddAtom(Chem.Atom(0))
        mw.AddBond(atom2,mw.GetNumAtoms()-1,Chem.BondType.SINGLE)
    frags = Chem.MolToSmiles(mw).split('.')
    frags_without_att = [delete_attachment(s) for s in frags]
    h = Chem.MolFragmentToSmiles(mol, hit) #需替换的子结构,也即hit对应的smi
    h_idx = frags_without_att.index(h)
    n_atta = frags[h_idx].count('*')
    frags.pop(h_idx)
    return frags,n_atta

def replaceFrag(raw_smi,replace_part_smi,input_file):
    mol=Chem.MolFromSmiles(raw_smi)
    patt=Chem.MolFromSmiles(replace_part_smi)
    hit=list(mol.GetSubstructMatch(patt))
    frags_for_connection,n_atta=fragmention(mol,hit)
    choices=[]
    if input_file.name.endswith('csv'):
        raw_smiles=list(set(pd.read_csv(input_file.name)['Smiles'].tolist()))[:200]
        for smi in raw_smiles:
            if smi.count('*')>=n_atta:
                choices.append(smi)
    elif input_file.name.endswith('sdf'):
        suppls=Chem.SDMolSupplier(input_file.name)
        for m in suppls:
            smi=Chem.MolToSmiles(m)
            if smi.count('*')>=n_atta:
                choices.append(smi)
    connected=Conn(frags_for_connection, choices, '*')
    res=list(set([delete_attachment(s) for s in connected]))
    outputPath=os.path.join(tmpdir,'replaceFrag_Molecules.sdf')
    writer=Chem.SDWriter(outputPath)
    genMols=[]
    for i,r in enumerate(res):
        genMol=Chem.MolFromSmiles(r)
        if genMol is not None:
            genMol.SetProp('_Name',str(i))
            ringNum=genMol.GetRingInfo().NumRings()
            writer.write(genMol)
            if ringNum>=3:
                genMols.append(genMol)
    image=highlightMolsReplace(genMols,frags_for_connection)
    writer.close()
    #高亮前10个代表分子的替换结构
    return outputPath,image

def connectFrag(filepath,files):
    conns={}
    for i,f in enumerate(files.split(',')):
        if f.endswith('.sdf'):
            suppls=Chem.SDMolSupplier(filepath+'/'+f)
            conns[i]=[Chem.MolToSmiles(s) for s in suppls]
        elif f.endswith('csv'):
            smiles=list(set(pd.read_csv(filepath+'/'+f)['Smiles'].tolist()))
            conns[i]=sample(smiles,30)
    frags=list(conns.values())
    #两个文件组合
    genSmis=[]
    if len(frags)==2:
        for f1 in frags[0]:
            for f2 in frags[1]:
                if f1.count('*')>0 and f2.count('*')>0:
                    conn=Mol_Conn(f1, f2,'*')
                    genSmis.extend([delete_attachment_H(c) for c in conn])
    #多个文件组合
    else:
        tempFrags=[]
        for f3 in frags[0]:
            for f4 in frags[1]:
                if f3.count('*')>0 and f4.count('*')>0:
                    conn1=Mol_Conn(f3, f4, '*')
                    tempFrags.extend(conn1)
        for fs in frags[2:]:
            midSmis=[]
            for tf in tempFrags:
                for f5 in fs:
                    if tf.count('*') > 0 and f5.count('*') > 0:
                        conn2=Mol_Conn(tf,f5,'*')
                        midSmis.extend(conn2)
            tempFrags=[]
            tempFrags.extend(midSmis)
        genSmis.extend([delete_attachment_H(tfs) for tfs in tempFrags])
    new_gemSmis=set(genSmis)
    mol_path=os.path.join(tmpdir,'connectFrag_Molecules.sdf')
    writer=Chem.SDWriter(mol_path)
    index,j=1,0
    genMols=[]
    for gs in new_gemSmis:
        gm=Chem.MolFromSmiles(gs)
        if gm is not None:
            gm.SetProp('_Name',str(index))
            writer.write(gm)
            index+=1
            if j<10 and 0.4<=QED.qed(gm)<=0.8:
                genMols.append(gm)
                j+=1
    writer.close()
    image=Draw.MolsToGridImage(genMols, molsPerRow=5, subImgSize=(500, 500))
    return mol_path,image

