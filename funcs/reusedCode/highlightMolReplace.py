# -*- coding: utf-8 -*-
'''
@Time    : 2025/4/18 20:29
@Author  : WY
@File    : highlightMolReplace.py
'''
import io
from rdkit import Chem
from PIL import Image as PILImage
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from funcs.reusedCode.utils import delete_attachment_H

def highlightOneMolReplace(raw_smi,replace_part_smi):
    mol=Chem.MolFromSmiles(raw_smi)
    patt=Chem.MolFromSmiles(replace_part_smi)
    hit_atoms=mol.GetSubstructMatches(patt)[0]
    hit_bonds=[]
    for bond in patt.GetBonds():
        aid1=hit_atoms[bond.GetBeginAtomIdx()]
        aid2=hit_atoms[bond.GetEndAtomIdx()]
        hit_bonds.append(mol.GetBondBetweenAtoms(aid1,aid2).GetIdx())
    drawer=rdMolDraw2D.MolDraw2DCairo(300, 300)
    colors=(0.96, 0.96, 0)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, highlightAtoms=list(hit_atoms), highlightBonds=hit_bonds,highlightAtomColors={i:colors for i in hit_atoms},highlightBondColors={i:colors for i in hit_bonds})
    drawer.FinishDrawing()
    img_bytes=drawer.GetDrawingText()
    image = PILImage.open(io.BytesIO(img_bytes))
    return image

def highlightNewMolsReplace(genMols,frags_for_connection):
    options = rdMolDraw2D.MolDrawOptions()
    options.legendFontSize = 30  # 设置 legend 字体大小
    highlightAtoms=[]
    new_genMols=[]
    for gm in genMols:
        num=0
        deleteIndexs=[]
        atomIndexs=[atom.GetIdx() for atom in gm.GetAtoms()]
        for i, ffc in enumerate(frags_for_connection):
            ffc_mol=Chem.MolFromSmiles(delete_attachment_H(ffc))
            ffc_indexs=gm.GetSubstructMatches(ffc_mol)
            if len(ffc_indexs)==1:
                deleteIndexs.extend(list(ffc_indexs[0]))
                num+=1
        if num==len(frags_for_connection):
            hit_atoms=list(set(atomIndexs)-set(deleteIndexs))
            highlightAtoms.append(hit_atoms)
            new_genMols.append((gm))
    colors=(0, 1, 0)
    image=Draw.MolsToGridImage(new_genMols,molsPerRow=5,subImgSize=(500, 500),highlightAtomLists=highlightAtoms,highlightAtomColors=[{i:colors for i in ha} for ha in highlightAtoms],legends=[str(x) for x in range(1,len(genMols)+1)],drawOptions=options)
    return image

