# -*- coding: utf-8 -*-
'''
@Time    : 2025/4/18 17:13
@Author  : WY
@File    : BreakMolecules.py
'''
import os
import random
import shutil
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools

from InterFragHub.funcs.reusedCode.RotaryFragMethod import RotaryFrag

def breakMols(file_obj):
    global img
    csvFile=os.getcwd()+'/output/breakMols/frags.csv'
    if file_obj.name.endswith('mol'):
        mol=Chem.MolFromMolFile(file_obj.name)
        smis=[]
        ligSmile=Chem.MolToSmiles(mol)
        if mol is not None:
            blocks, blocks_list, interval=RotaryFrag(mol)
            for b in blocks:
                smis.append(Chem.MolToSmiles(b))
        data=pd.DataFrame(np.array(smis),columns=['Smiles'])
        data['ligandSmile']=ligSmile
        data.to_csv(csvFile,index=False)
        mols = []
        for s in random.sample(smis,10):
            mols.append(Chem.MolFromSmiles(s))
        img=Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(300, 300), legends=[str(x+1) for x in range(10)])
    elif file_obj.name.endswith('sdf'):
        suppl=Chem.SDMolSupplier(file_obj.name)
        smis=[]
        for i, s in enumerate(suppl):
            if s is not None:
                ligSmile=Chem.MolToSmiles(s)
                blocks, blocks_list, interval = RotaryFrag(s)
                for b in blocks:
                    smis.append([Chem.MolToSmiles(b),ligSmile])
        data=pd.DataFrame(smis,columns=['Smiles','ligandSmile'])
        data.to_csv(csvFile,index=False)
        mols=[]
        for smi in random.sample(data['Smiles'].tolist(),10):
            mols.append(Chem.MolFromSmiles(smi))
        img=Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(300, 300), legends=['' for x in mols])
    return csvFile,img

