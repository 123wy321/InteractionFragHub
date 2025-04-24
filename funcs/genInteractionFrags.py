# -*- coding: utf-8 -*-
'''
@Time    : 2025/4/18 20:03
@Author  : WY
@File    : genInteractionFrags.py
'''
import os
import subprocess
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from biopandas.pdb import PandasPdb
from tqdm import tqdm

from funcs.interaction_components.plinteraction import get_interactions
from funcs.reusedCode.RotaryFragMethod import RotaryFrag
from funcs.reusedCode.calECscore import getECPairs
from funcs.reusedCode.utils import listToStr, mol_with_atomLabel, transformDF, updateResName
# from InterFragHub.funcs.reusedCode.utils import mol_with_atomLabel
# from InterFragHub.funcs.reusedCode.calECscore import getECPairs
# from InterFragHub.funcs.reusedCode.RotaryFragMethod import RotaryFrag
# from InterFragHub.funcs.interaction_components.plinteraction import get_interactions
# from InterFragHub.funcs.reusedCode.utils import *

dummy = Chem.MolFromSmiles('*')
Hydrogen = Chem.MolFromSmarts('[H]')

def getRotaryFrags(tab,lig):
    lig_mol = Chem.MolFromMolFile(lig, removeHs=False)
    lig_mol = mol_with_atomLabel(lig_mol)
    fragsList,bondsNum,interval = RotaryFrag(lig_mol)
    print(f'[{tab}] finish split, bonds Num:[{bondsNum}],frag Num:[{len(fragsList)}],calculate Time:{interval}!')
    return fragsList

def creatPqrFile(ligfile_path, pdbfile_path):
    script_dir = "ESP_DNN"
    command1 = (
        f"cd {script_dir} && "
        f'conda run -n esp-dnn-env python -m esp_dnn.predict -m ligand -i {ligfile_path}'
    )
    # 使用 subprocess 执行命令
    process = subprocess.Popen(''.join(command1), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True)
    for line in iter(process.stdout.readline, ''):
        print(line, end='')
    process.stdout.close()
    process.wait()
    command2 = (
        f"cd {script_dir} && "
        f'conda run -n esp-dnn-env python -m esp_dnn.predict -m protein -i {pdbfile_path}'
    )
    # 使用 subprocess 执行命令
    process = subprocess.Popen(''.join(command2), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True)
    for line in iter(process.stdout.readline, ''):
        print(line, end='')
    process.stdout.close()
    process.wait()
    return 'progress end'

def creatTomeshFile(pqr_filename, tmesh_filename):
    print(pqr_filename, tmesh_filename)
    script_dir = "/data/ywang/genECFrags/esp-surface-generator"
    command1 = [
        f"cd {script_dir} && "
        f'conda run -n InteractionFragHub node cli.js {pqr_filename} {tmesh_filename}'
    ]
    # 使用 subprocess 执行命令
    process = subprocess.Popen(''.join(command1), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True)
    for line in iter(process.stdout.readline, ''):
        print(line, end='')
    process.stdout.close()
    process.wait()
    return 'progress end'

def lig_pqr_file_process(path):
    for pqr_file in os.listdir(path):
        if pqr_file.endswith('.pqr'):
            try:
                f = open(path+'/'+pqr_file, 'r+')
                f_list = f.readlines()
                if f_list==[]:  # 排除无法计算电荷的分子
                    continue
                else:
                    for i in range(len(f_list)):
                        f_list[i] = f_list[i].split(' ')
                        while '' in f_list[i]:
                            f_list[i].remove('')  # 去除空字符串
                        f_list[i] = [j for j in f_list[i] if j.strip() != '_']  # 去除'_'
                        f_list[i][-1] = f_list[i][-1][:-1]  # 去除最后一个元素的\n
                        if len(f_list[i]) != 10:
                            if len(f_list[i][5:-2]) == 2:  # 将坐标规范化(两个坐标合并) ①
                                new_coord = []
                                for coord in f_list[i][5:-2]:
                                    if coord.count('.') > 1:
                                        before_str = coord[:coord.find('.')] + '.' + coord[coord.find('.') + 1:coord.find('.') + 4]
                                        after_str = coord[coord.find('.') + 4:coord.rfind('.')] + '.' + coord[coord.rfind('.') + 1:]
                                        new_coord.extend([before_str, after_str])
                                    elif coord.count('.') == 1:
                                        new_coord.append(coord)
                                f_list[i][5:-2] = new_coord
                            if len(f_list[i][5:-2]) == 1:  # 将坐标规范化(三个坐标合并) ②
                                new_coord1 = []
                                raw_coords = f_list[i][5:-2][0].split('.')
                                before_str1 = raw_coords[0]+'.'+raw_coords[1][:3]
                                mid_str1 = raw_coords[1][3:]+'.'+raw_coords[2][:3]
                                after_str1 = raw_coords[2][3:] + '.' + raw_coords[3]
                                new_coord1.extend([before_str1, mid_str1, after_str1])
                                f_list[i][5:-2] = new_coord1
                    f_list_df = pd.DataFrame(f_list, columns=['recordName', 'atomIndex', 'atomName', 'residueName','residueNumber', 'X', 'Y', 'Z', 'charge', 'radius'])
                    new_f_list_df = updateResName(f_list_df,4)
                    new_f_list_df.to_csv(path+'/'+pqr_file.split('.')[0] + '_pqr.csv', index=False)
                f.close()
            except Exception as e:
                print(pqr_file + f' has error:{e}')

def pro_pqr_file_process(path):
    for pqr_file in os.listdir(path):
        if pqr_file.endswith('.pqr'):
            try:
                f = open(path+'/'+pqr_file, 'r+')
                f_list = f.readlines()
                if f_list==[]:  # 排除无法计算电荷的分子
                    continue
                else:
                    for i in range(len(f_list)):
                        f_list[i] = f_list[i].split(' ')
                        while '' in f_list[i]:
                            f_list[i].remove('')  # 去除空字符串
                        f_list[i] = [j for j in f_list[i] if j.strip() != '_']  # 去除'_'
                        f_list[i][-1] = f_list[i][-1][:-1]  # 去除最后一个元素的\n
                        if len(f_list[i]) != 11:
                            if len(f_list[i][3]) != 3:  # 处理原子和残基种类相连情况 ①
                                mid1 = f_list[i][2]
                                f_list[i][2] = mid1[:-3]
                                f_list[i].insert(3, mid1[-3:])
                            if len(f_list[i][4]) > 1:  # 处理链名和残基序号相连情况 ②
                                mid = f_list[i][4]
                                f_list[i][4] = mid[0]
                                f_list[i].insert(5, mid[1:])
                            if len(f_list[i][6:-2])==2:  # 将坐标规范化(先按两个坐标合并情况处理) ③
                                new_coord = []
                                for coord in f_list[i][6:-2]:
                                    if coord.count('.') > 1:
                                        before_str = coord[:coord.find('.')] + '.' + coord[coord.find('.') + 1:coord.find('.') + 4]
                                        after_str = coord[coord.find('.') + 4:coord.rfind('.')] + '.' + coord[coord.rfind('.') + 1:]
                                        new_coord.extend([before_str, after_str])
                                    elif coord.count('.') == 1:
                                        new_coord.append(coord)
                                f_list[i][6:-2] = new_coord
                            if len(f_list[i][6:-2])==1:  # 将坐标规范化(三个坐标合并)  ④
                                new_coord1 = []
                                raw_coords = f_list[i][6:-2][0].split('.')
                                before_str1 = raw_coords[0]+'.'+raw_coords[1][:3]
                                mid_str1 = raw_coords[1][3:]+'.'+raw_coords[2][:3]
                                after_str1 = raw_coords[2][3:] + '.' + raw_coords[3]
                                new_coord1.extend([before_str1, mid_str1, after_str1])
                                f_list[i][6:-2] = new_coord1
                    f_list_df = pd.DataFrame(f_list,columns=['recordName', 'atomIndex', 'atomName', 'residueName', 'chainName','residueNumber', 'X', 'Y', 'Z', 'charge', 'radius'])
                    new_f_list_df=updateResName(f_list_df,5)
                    new_f_list_df.to_csv(path+'/'+pqr_file.split('.')[0]+'_pqr.csv', index=False)
                f.close()
            except Exception as e:
                print(pqr_file+f' has error:{e}')

def tmesh_file_process(path):
    for tmesh_file in os.listdir(path):
        if tmesh_file.endswith('.tmesh'):
            try:
                f = open(path+'/'+tmesh_file, 'r+')
                f_list = f.readlines()
                if f_list == []:  # 排除无法计算电荷的分子
                    continue
                else:
                    vertex_num = int(f_list[0])
                    for i in range(1, vertex_num + 1):
                        f_list[i] = f_list[i].split(' ')
                        while '' in f_list[i]:
                            f_list[i].remove('')  # 去除空字符串
                        f_list[i] = [j for j in f_list[i] if j.strip() != '_']  # 去除'_'
                        f_list[i][-1] = f_list[i][-1][:-1]  # 去除最后一个元素的\n
                    f_list_df = pd.DataFrame(f_list[1:vertex_num + 1],columns=['X', 'Y', 'Z', 'nX', 'nY', 'nZ', 's', 'pqr_esp', 'QChem_esp'])
                    f_list_df.to_csv(path+'/'+tmesh_file.split('.')[0] + '_tmesh.csv', index=False)
                f.close()
            except:
                print(tmesh_file+' has error')

def deleteTmeshFiles(path):
    for file in tqdm(os.listdir(path)):
        if file.endswith('.tmesh'):
            os.remove(path+'/'+file)

def lig_atomTmesh_file_process(path):
    pdb_path,pqr_path,temsh_path,tab='','','',''
    for file in os.listdir(path):
        if file.endswith('.mol.pdb'):
            pdb_path=path+'/'+file
            tab=file.split('.mol.pdb')[0]
        elif file.endswith('_pqr.csv'):
            pqr_path=path+'/'+file
        elif file.endswith('_tmesh.csv'):
            temsh_path=path+'/'+file
    try:
        pqr_data=generateNewPQR(pdb_path,pqr_path,'ligand')
        tmesh_data=pd.read_csv(temsh_path)
        new_tmesh_data=tmesh_data.copy(deep=False)
        for index,row in tmesh_data.iterrows():
            data=calWeight(row,pqr_data)
            new_tmesh_data.loc[index,'own_atom']=data['atomName'].values.tolist()[0]
        new_tmesh_data1 = new_tmesh_data.drop(['nX', 'nY', 'nZ', 's', 'QChem_esp'],axis=1)
        new_tmesh_data1.to_csv(path+'/'+tab+'_atomTmesh.csv',index=False)
        print(f'success generate {tab}_atomTmesh.csv!')
        return path+'/'+tab+'_atomTmesh.csv'
    except Exception as e:
        return f'{tab}_atomTmesh has error:{e}'

def pro_atomTmesh_file_process(path):
    pdb_path,pqr_path,temsh_path,tab='','','',''
    for file in os.listdir(path):
        if file.endswith('.pdb'):
            pdb_path=path+'/'+file
            tab=file.split('.pdb')[0]
        elif file.endswith('_pqr.csv'):
            pqr_path=path+'/'+file
        elif file.endswith('_tmesh.csv'):
            temsh_path=path+'/'+file
    try:
        pqr_data=generateNewPQR(pdb_path,pqr_path,'protein')
        tmesh_data=pd.read_csv(temsh_path)
        new_tmesh_data=tmesh_data.copy(deep=False)
        for index, row in tmesh_data.iterrows():
            data = calWeight(row, pqr_data)
            new_tmesh_data.loc[index, 'own_atom'] = data['atomName'].values.tolist()[0]
            new_tmesh_data.loc[index, 'residueName'] = data['residueName'].values.tolist()[0]
            new_tmesh_data.loc[index, 'residueNumber'] = str(data['residueNumber'].values.tolist()[0])
        new_tmesh_data1 = new_tmesh_data.drop(['nX', 'nY', 'nZ', 's', 'QChem_esp'],axis=1)
        new_tmesh_data1.to_csv(path+'/'+tab+'_atomTmesh.csv',index=False)
        print(f'success generate {tab}_atomTmesh.csv!')
        return path+'/'+tab+'_atomTmesh.csv'
    except Exception as e:
        return f'{tab}_atomTmesh has error:{e}'
    
def generateNewPQR(pdb_path,pqr_path,kind):
    pqr_data=pd.read_csv(pqr_path)
    pdb_data=PandasPdb().read_pdb(pdb_path)
    if kind=='ligand':
        raw_data=pdb_data.df['HETATM']
    elif kind=='protein':
        raw_data=pdb_data.df['ATOM']
    new_pqr_data=pqr_data.copy(deep=False)
    new_pqr_data.drop('atomIndex',axis=1,inplace=True)
    for index,row in pqr_data.iterrows():
        try:
            if len(raw_data[(raw_data['x_coord']==row['X'])&(raw_data['y_coord']==row['Y'])&(raw_data['z_coord']==row['Z'])])!=0:
                atomName=raw_data.loc[(raw_data['x_coord']==row['X'])&(raw_data['y_coord']==row['Y'])&(raw_data['z_coord']==row['Z']),'atom_name'].values.tolist()[0]
                new_pqr_data.loc[index,'atomName']=atomName
            else:
                new_pqr_data.drop(index=index,inplace=True)
        except Exception as e:
            print(e)
    new_pqr_data.insert(loc=1,column='atomIndex',value=[i for i in range(1,len(new_pqr_data)+1)])

    return new_pqr_data

def calWeight(tmesh,pqr_data):
    mid_pqr_data=pqr_data.copy(deep=False)
    mid_pqr_data['weight']=1-(((tmesh['X']-mid_pqr_data['X'])**2+(tmesh['Y']-mid_pqr_data['Y'])**2+(tmesh['Z']-mid_pqr_data['Z'])**2).apply(np.sqrt)/mid_pqr_data['radius'])
    mid_pqr_data['pqr_esp']=tmesh['pqr_esp']
    data=mid_pqr_data[mid_pqr_data['weight'] == mid_pqr_data['weight'].max()]
    return data
    
def genECFrags(lig,pro):
    ligpath,propath=lig.name,pro.name
    tab=ligpath.split('\\')[-1].split('_lig')[0]
    file_lig,file_pro=ligpath.split('\\')[-1],propath.split('\\')[-1]

    #生成pqr文件
    creatPqrFile(ligpath,propath)
    lig_pqrfile_path=os.getcwd()+f'/output/genInterFrags/{file_lig}.pdb.pqr'
    pro_pqrfile_path=os.getcwd()+f'/output/genInterFrags/{file_pro}.pdb.pqr'
    lig_tmeshfile_path=os.getcwd()+f'/output/genInterFrags/{file_lig[:-4]}.tmesh'
    pro_tmeshfile_path=os.getcwd()+f'/output/genInterFrags/{file_pro[:-4]}.tmesh'
    #生成tmesh文件
    creatTomeshFile(lig_pqrfile_path, lig_tmeshfile_path)
    creatTomeshFile(pro_pqrfile_path, pro_tmeshfile_path)
    #生成atomtmesh文件
    ligs=os.getcwd()+f'/output/genInterFrags/lig_file'
    pros=os.getcwd()+f'/output/genInterFrags/pro_file'
    savepath=os.getcwd()+f'/output/genInterFrags'
    lig_pqr_file_process(ligs)
    pro_pqr_file_process(pros)
    tmesh_file_process(ligs)
    tmesh_file_process(pros)
    lig_atomTmesh_file_process(ligs)
    pro_atomTmesh_file_process(pros)

    ligAtomTmesh,proAtomTmesh=f'{tab}_lig_addH_atomTmesh.csv',f'{tab}_poc_addH_atomTmesh.csv'
    fragsList=getRotaryFrags(tab,ligpath)  # RotaryFrag方法
    frags_dict={}
    frags_data_list=[]
    if len(fragsList) > 0:
        ECpairs=getECPairs(ligs+'/'+ligAtomTmesh,pros+'/'+proAtomTmesh)
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
        EC_frag_path=savepath+f'/{tab}_ECfrags.csv'
        frags_df.to_csv(EC_frag_path,index=False)
        return EC_frag_path
    else:
        return '无法生成静电互补片段集'

def genHbondFrags(lig,pro):
    ligpath,propath=lig.name,pro.name
    tab=ligpath.split('\\')[-1].split('_lig')[0]
    ligand=mol_with_atomLabel(Chem.MolFromMolFile(ligpath,removeHs=False))
    protein=Chem.MolFromPDBFile(propath,removeHs=False)
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
        hbond_frag_path=os.getcwd()+f'/output/genInterFrags/{tab}_Hbondsfrags.csv'
        hbonds_frags_df.to_csv(hbond_frag_path, index=False)
        return hbond_frag_path
    else:
        print('无法生成氢键片段集')

def genHydrophobicFrags(lig,pro):
    ligpath,propath=lig.name,pro.name
    tab=ligpath.split('\\')[-1].split('_lig')[0]
    ligand=mol_with_atomLabel(Chem.MolFromMolFile(ligpath, removeHs=False))
    protein=Chem.MolFromPDBFile(propath, removeHs=False)
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
        hydrophobic_frag_path=os.getcwd()+f'/output/genInterFrags/{tab}_Hydrophobicfrags.csv'
        hydrophobic_frags_df.to_csv(hydrophobic_frag_path, index=False)
        return hydrophobic_frag_path
    else:
        print('无法生成疏水片段集')