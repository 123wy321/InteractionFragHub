# -*- coding: utf-8 -*-
'''
@Time    : 2025/4/18 20:25
@Author  : WY
@File    : searchFrags.py
'''
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import uuid
def selectNotNoneMol(smile):
    mol=Chem.MolFromSmiles(smile)
    if mol is not None:
        return 'True'
    else:
        return 'False'
    
def findFrags(res,interaction):
    #查找片段集文件
    filename=f'{interaction}_{res}.csv'
    fragpath=os.getcwd()+f'/data/{interaction}/'+filename
    options = rdMolDraw2D.MolDrawOptions()
    options.legendFontSize = 30  # 设置 legend 字体大小

    #显示代表分子
    data = pd.read_csv(fragpath)
    data.drop_duplicates(subset=['Smiles'],keep='first',inplace=True)
    data['SelectNotMol']=data['Smiles'].map(selectNotNoneMol)
    mid_data=data[data['SelectNotMol']=='True']
    mid_data1=mid_data.sample(n=12)
    unique_name=str(uuid.uuid4())
    sdfpath=os.getcwd()+f'/output/searchFrags/{unique_name}.sdf'
    pngpath=os.getcwd()+f'/output/searchFrags/{unique_name}.png'

    #生成sdf文件
    writer=Chem.SDWriter(sdfpath)
    mols=[]
    for index,row in mid_data1.iterrows():
        mol=Chem.MolFromSmiles(row['Smiles'])
        mol.SetProp('_Name',str(index+1))
        mol.SetProp('PDBsource',row['PDBsource'])
        mols.append(mol)
        writer.write(mol)
    writer.close()

    #生成png文件
    img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(500, 500), legends=[str(x) for x in range(1,13)],drawOptions=options)
    img.save(pngpath)

    #生成df文件
    mid_data1.insert(0,'Index',range(1,len(mid_data1)+1))
    df=mid_data1[['Index','PDBsource','Smiles']]
    df=mid_data1[['PDBsource','Smiles']]
    return sdfpath,img,df

