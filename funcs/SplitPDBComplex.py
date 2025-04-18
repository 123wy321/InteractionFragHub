# -*- coding: utf-8 -*-
'''
@Time    : 2025/4/18 20:00
@Author  : WY
@File    : SplitPDBComplex.py
'''
import os
import pandas as pd
from Bio.PDB import PDBParser
from biopandas.pdb import PandasPdb
from rdkit import Chem
from rdkit.Chem import AllChem

tmpdir=os.getcwd()
def splitPDB(pdb,ligandName):
    #先测试一个配体
    pdbName=pdb.name.split('/')[-1][:-4]
    if ',' not in ligandName:
        ligandNameList=[ligandName]
    else:
        ligandNameList=ligandName.split(',')
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdbName, pdb.name)
    lig_dict={}
    chainInfo=[]
    for model in structure:
        for chain in model:
            chainInfo.append(chain.id)
            for residue in chain:
                if residue.id[0][2:] in ligandNameList:
                    if residue.id[0][2:] not in lig_dict.keys():
                        lig_dict[residue.id[0][2:]]=[chain.id]
                    elif residue.id[0][2:] in lig_dict.keys():
                        lig_dict[residue.id[0][2:]].append(chain.id)
    flag=True
    for key,value in lig_dict.items():
        if len(value)!=len(chainInfo):
            flag=False
            break
    # 获取蛋白并生成文件
    ppdb1 = PandasPdb()
    ppdb1.read_pdb(pdb.name)
    protein_biop = ppdb1.df['ATOM']
    if flag == True:
        right_pro = protein_biop[(protein_biop['chain_id'] == chainInfo[0])]   #配体在每条链上都出现，取一条蛋白链;
    else:
        right_pro = protein_biop   #存在配体在某条链上未出现，取所有的蛋白链
    ppdb1.df['ATOM'] = right_pro
    ppdb1.df['HETATM'] = pd.DataFrame()
    ppdb1.df['ANISOU'] = pd.DataFrame()
    ppdb1.df['OTHERS'] = ppdb1.df['OTHERS'][(ppdb1.df['OTHERS']['record_name'] == 'CONECT') | (ppdb1.df['OTHERS']['record_name'] == 'END')]
    pro_pdb_path=os.path.join(tmpdir,f'{pdbName}_pro.pdb')
    ppdb1.to_pdb(pro_pdb_path)
    print(f'[protein]success generate {pdbName}_pro.pdb')
    lig_mols,poc_pdbs=[],[]
    #获取配体并生成文件
    for ligName,chainId in lig_dict.items():
        ppdb = PandasPdb()
        ppdb.read_pdb(pdb.name)
        ligand_biop = ppdb.df['HETATM']
        right_lig = ligand_biop[(ligand_biop['residue_name']==ligName)&(ligand_biop['chain_id']==chainId[0])]
        ppdb.df['HETATM'] = right_lig
        ppdb.df['ATOM'] = pd.DataFrame()
        ppdb.df['ANISOU'] = pd.DataFrame()
        ppdb.df['OTHERS'] = ppdb.df['OTHERS'][(ppdb.df['OTHERS']['record_name'] == 'CONECT') | (ppdb.df['OTHERS']['record_name'] == 'END')]
        lig_pdb_path=os.path.join(tmpdir,f'{pdbName}_{ligName}_{chainId[0]}_lig.pdb')
        ppdb.to_pdb(lig_pdb_path)
        lig_mols.append(lig_pdb_path)
        mol=AllChem.MolFromPDBFile(lig_pdb_path)
        if mol!=None:
            lig_mol_path=os.path.join(tmpdir,f'{pdbName}_{ligName}_{chainId[0]}_lig.mol')
            Chem.MolToMolFile(mol,lig_mol_path)
            lig_mols.append(lig_mol_path)
            print(f'[ligand]success generate {pdbName}_{ligName}_{chainId[0]}_lig.mol&pdb')
        else:
            print(f'{pdbName}_{ligName}_{chainId[0]}_lig.pdb can not convert into mol file')
    return pro_pdb_path,lig_mols

