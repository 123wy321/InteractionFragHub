# -*- coding: utf-8 -*-
'''
@Time    : 2025/4/18 17:13
@Author  : WY
@File    : BreakMolecules.py
'''
import os
import shutil
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

tmpdir=os.getcwd()
def breakMols(file_obj):
    print('临时文件夹地址：{}'.format(tmpdir))
    print('上传文件的地址：{}'.format(file_obj.name))  #输出上传后的文件在gradio中保存的绝对地址

    #将文件复制到临时目录中
    shutil.copy(file_obj.name,tmpdir)
    if file_obj.name.endswith('mol'):
        mol=Chem.MolFromMolFile(file_obj.name)
        smiles=[]
        if mol is not None:
            blocks, blocks_list, interval=RotaryFrag(mol)
            for b in blocks:
              smiles.append(Chem.MolToSmiles(b))
        outputPath=os.path.join(tmpdir,'frag.csv')
        pd.DataFrame(np.array(smiles),columns=['Smiles']).to_csv(outputPath,index=False)

        data=pd.read_csv(outputPath)
        smis = data['Smiles'].tolist()[:10]
        mols = []
        for s in smis:
            mols.append(Chem.MolFromSmiles(s))
        img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(300, 300), legends=['' for x in mols])
        return outputPath,img

