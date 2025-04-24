# -*- coding: utf-8 -*-
'''
@Time    : 2025/4/19 16:21
@Author  : WY
@File    : findFrags.py
'''
import os
import shutil

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D


def findFrags(res,interaction):
    #查找片段集文件
    filename=f'{interaction}_{res}.csv'
    fragpath=fragsPath+f'{interaction}/'+filename
    shutil.copy(fragpath, tmpdir)
    outputPath = os.path.join(tmpdir, filename)
    options = rdMolDraw2D.MolDrawOptions()
    options.legendFontSize = 30  # 设置 legend 字体大小

    #显示代表分子
    data = pd.read_csv(fragpath)
    data.drop_duplicates(subset=['Smiles'],keep='first',inplace=True)
    data['HeavyAtoms']=data['Smiles'].map(calHeavyAtoms)
    data['RingNums']=data['Smiles'].map(calRingNums)
    data['IsP']=data['Smiles'].map(findPAtom)
    result=data[(data['HeavyAtoms']>10)&(data['IsP']==True)&(data['RingNums'].isin([1,2,3]))]
    df=result.sample(n=20).iloc[:,1:3]
    indexs=[i for i in range(1,21)]
    df.insert(loc=0,column='Index',value=indexs)
    smis = df['Smiles'].tolist()[:12]
    mols = []
    for s in smis:
        mols.append(Chem.MolFromSmiles(s))
    img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(500, 500), legends=[str(x) for x in range(1,13)],drawOptions=options)
    return outputPath,img,df

