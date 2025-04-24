# -*- coding: utf-8 -*-
'''
@Time    : 2025/4/18 20:32
@Author  : WY
@File    : MolOptimizeMethods.py
'''
import os
from random import sample
import subprocess
import uuid
import pandas as pd
from rdkit import Chem
from rdkit.Chem import QED, Draw,PandasTools
from rdkit.Chem.Draw import rdMolDraw2D
from funcs.reusedCode.highlightMolReplace import highlightNewMolsReplace
from funcs.reusedCode.utils import delete_attachment, Conn, Mol_Conn, delete_attachment_H, selectNotNoneMol
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
    unique_name=str(uuid.uuid4())
    sdfPath=os.getcwd()+f'/output/molOptimize/{unique_name}.sdf'
    writer=Chem.SDWriter(sdfPath)
    for i,r in enumerate(res):
        genMol=Chem.MolFromSmiles(r)
        if genMol is not None:
            genMol.SetProp('_Name',str(i))
            writer.write(genMol)
    writer.close()
    df=PandasTools.LoadSDF(sdfPath)
    image=highlightNewMolsReplace(df.sample(10)['ROMol'],frags_for_connection)
    #高亮前10个代表分子的替换结构
    return sdfPath,image

def connectFrag(files):
    filespath=[f.name for f in files]
    conns={}
    for i,f in enumerate(filespath):
        suppls=Chem.SDMolSupplier(f)
        conns[i]=[Chem.MolToSmiles(s) for s in suppls]
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
    unique_name=str(uuid.uuid4())
    sdfPath=os.getcwd()+f'/output/molOptimize/{unique_name}.sdf'
    writer=Chem.SDWriter(sdfPath)
    for gs in new_gemSmis:
        gm=Chem.MolFromSmiles(gs)
        if gm is not None:
            writer.write(gm)
    writer.close()
    df=PandasTools.LoadSDF(sdfPath)
    image=Draw.MolsToGridImage(df.sample(10)['ROMol'], molsPerRow=5, subImgSize=(500, 500))
    return sdfPath,image

def linkGenImages(smiPath):
    options = rdMolDraw2D.MolDrawOptions()
    options.legendFontSize = 30  # 设置 legend 字体大小
    smiFile=open(smiPath,'r')
    smis=[]
    warheads=''
    for line in smiFile:
        if line.startswith('Warheads:'):
            warheads=line.strip('\n')[11:-2].split('|')
        if not line.startswith('Molecules:') and not line.startswith('Warheads:'):
            smi=line.strip('\n')
            if len(smis)<100:
                smis.append(smi)
    data=pd.DataFrame(smis,columns=['Smiles'])
    data['SelectNotMol']=data['Smiles'].map(selectNotNoneMol)
    data1=data[data['SelectNotMol']=='True'].sample(10)
    mols=[Chem.MolFromSmiles(s) for s in data1['Smiles'].tolist()]
    img=highlightNewMolsReplace(mols,warheads)
    return img

def Linkinvent(warheads):
    input_uuid=str(uuid.uuid4())
    output_uuid=str(uuid.uuid4())

    #生成warheads的smi文件
    input_smiPath=os.getcwd()+f'/output/molOptimize/{input_uuid}.smi'
    input_smifile=open(input_smiPath, 'w')
    input_smifile.write(f"1\t1_1|1_2\t{warheads}\n")
    input_smifile.close()

    python_path="/data/Software/anaconda3/envs/reinvent.v3.2/bin/python"
    script_dir="/data/ywang/Reinvent"
    command = (
        f"cd {script_dir} && "
        "conda run -n reinvent.v3.2 --no-capture-output &&"
        f"{python_path} genMols_input_local.py --base_config configs/config.json --run_config configs/LinkInvent_genmols_config.json --warheads_file {input_smiPath} --save_folder {os.getcwd()+f'/output/molOptimize'}--output_uuid {output_uuid}"
    )
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,text=True)
    for line in iter(process.stdout.readline, ''):
        if 'report' in line:
            progress = progress + 10
            # socketio.emit("progress_update", {"progress": progress})  # 通过 WebSocket 发送更新
            print(progress)
        print(line, end='')
    process.stdout.close()
    process.wait()
    # input_smiPath=os.getcwd()+f'/output/molOptimize/{output_uuid}.smi'
    input_smiPath=os.getcwd()+f'/output/molOptimize/1020e835-78aa-4717-be4c-9649eefb43bd.smi'
    image=linkGenImages(input_smiPath)
    return input_smiPath,image
