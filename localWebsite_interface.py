# -*- coding: utf-8 -*-
'''
@Time    : 2025/4/18 17:09
@Author  : WY
@File    : localWebsite_interface.py
'''
import os
import gradio as gr
from funcs.BreakMolecules import breakMols
from funcs.SplitPDBComplex import splitPDB
from funcs.genInteractionFrags import genECFrags,genHbondFrags,genHydrophobicFrags
from funcs.reusedCode.highlightMolReplace import highlightOneMolReplace
from funcs.searchFrags import findFrags
from funcs.MolOptimizeMethods import Linkinvent, replaceFrag,connectFrag

res_types=['GLY','ALA','SER','THR','CYS','VAL','LEU','ILE','MET','PRO','PHE','TYR','TRP','ASP','GLU','ASN','GLN','HIS','LYS','ARG']
interaction_types=['ECFrags','HydrophobicFrags','HbondFrags']

if __name__ == '__main__':
    with gr.Blocks(css='./static/css/style.css') as demo:
        gr.Image(value='./static/picture/logo.png',width=300,show_download_button=False,container=False)
        #分子打碎
        with gr.Tab(label="Break Molecule", elem_id="body1"):
            gr.Markdown('Use RotaryFrag to break Molecules.',elem_id='markdown1')
            with gr.Row():
                with gr.Column(scale=1):
                    inputs = gr.components.File(label='Molecule File(.sdf,.mol)',elem_id='component1')
                with gr.Column(scale=1):
                    outputs = gr.components.File(label='Fragment File(.csv)')
            with gr.Row():
                breakdown = gr.Button(value='Break Molecule',elem_id="btn1")
            with gr.Row():
                example1 = gr.Image(label='Display Fragments')
            breakdown.click(fn=breakMols, inputs=inputs, outputs=[outputs, example1])
        #分离蛋白复合物
        with gr.Tab(label="Split PDB Complex", elem_id="body2"):
            gr.Markdown('Please enter PDB File, then click Split button to generate PDB file of protein and MOL file of ligand.',elem_classes='markdown')
            with gr.Row():
                with gr.Column(scale=1):
                    inputs = gr.components.File(label='Upload PDB file of protein complex')
                with gr.Column(scale=1):
                    ligand_name = gr.Textbox(label='Enter ligand names',placeholder='For example: L7H')
            with gr.Row():
                breakdown = gr.Button(value='Split PDB Complex',elem_id="btn2")
            with gr.Row():
                with gr.Column(scale=1):
                    pro_outputs = gr.components.File(label='Protein file(.pdb)')
                with gr.Column(scale=1):
                    lig_outputs = gr.components.File(label='Ligand file(.mol&.pdb)')
            breakdown.click(fn=splitPDB, inputs=[inputs, ligand_name], outputs=[pro_outputs, lig_outputs])
        #生成相互作用片段集
        with gr.Tab(label="Generate Interaction Fragments", elem_id="body3"):
            gr.Markdown('Click to generate different interaction frags.')
            with gr.Tab(label="Generate EC Fragments", elem_classes="body3_mini_1"):
                with gr.Row():
                    with gr.Column(scale=1):
                        lig1=gr.components.File(label='Ligand File(.mol)')
                    with gr.Column(scale=1):
                        pro1=gr.components.File(label='Protein File(.pdb)')
                with gr.Row():
                    genFrags1=gr.Button(value='Generate EC Fragments',elem_id="btn3")
                with gr.Row():
                    frags1=gr.Files(label='Generated ECFrags Files(.csv)')   
                genFrags1.click(fn=genECFrags,inputs=[lig1,pro1],outputs=frags1)     
            with gr.Tab(label="Generate Hbond Fragments", elem_classes="body3_mini_2"):
                with gr.Row():
                    with gr.Column(scale=1):
                        lig2=gr.components.File(label='Ligand File(.mol)')
                    with gr.Column(scale=1):
                        pro2=gr.components.File(label='Protein File(.pdb)')
                with gr.Row():
                    genFrags2=gr.Button(value='Generate Hbond Fragments',elem_id="btn4")
                with gr.Row():
                    frags2=gr.Files(label='Generated HbondFrags Files(.csv)')
                genFrags2.click(fn=genHbondFrags, inputs=[lig2,pro2], outputs=frags2)
            with gr.Tab(label="Generate Hydrophobic Fragments", elem_classes="body3_mini_3"):
                with gr.Row():
                    with gr.Column(scale=1):
                        lig3=gr.components.File(label='Ligand File(.mol)')
                    with gr.Column(scale=1):    
                        pro3=gr.components.File(label='Protein File(.pdb)')
                with gr.Row():
                    genFrags3=gr.Button(value='Generate Hydrophobic Fragments',elem_id="btn5")
                with gr.Row():
                    frags3=gr.Files(label='Generated HydrohobicFrags Files(.csv)')
                genFrags3.click(fn=genHydrophobicFrags, inputs=[lig3,pro3], outputs=frags3)
        #查找片段
        with gr.Tab(label="Search Fragments", elem_id="body4"):
            gr.Markdown('Please select residue and interaction types,then click Search button to get the needed fragment file.')
            with gr.Row():
                with gr.Column(scale=1):
                    res = gr.Dropdown(res_types, label="Residue types")
                with gr.Column(scale=1):
                    interaction = gr.Dropdown(interaction_types, label='Interaction types')
            with gr.Row():
                search = gr.Button(value='SEARCH',elem_id="btn6")
            with gr.Row():
                outputs = gr.components.File(label='Interaction fragment file')
            with gr.Row():
                with gr.Column(scale=1):
                    example2 = gr.Image(label='Display Fragments')
                with gr.Column(scale=1):
                    gr.Markdown('Fragments Information')
                    mol_df = gr.DataFrame()
            search.click(fn=findFrags, inputs=[res, interaction], outputs=[outputs, example2,mol_df])
        #基于片段的分子优化
        with gr.Tab(label="Optimize molecules based on Fragments", elem_id="body5"):
            gr.Markdown('Select different methods to generate molecules.')
            with gr.Tab(label="Replace Fragments", elem_classes="body5_mini_1"):
                with gr.Row():
                    with gr.Column(scale=1):
                        raw_smi = gr.Textbox(label="Primitive Molecule Smile", elem_classes="input-box")
                        replace_part_smi = gr.Textbox(label="Substitute Smiles", elem_classes="input-box")
                        highlight = gr.Button(value='Highlight Substitute Part',elem_id="btn7")
                    with gr.Column(scale=1):
                        with gr.Group():
                            with gr.Row():
                                highlightMol = gr.Image(label='Relationship of primitive molecule and substitute part',elem_id="output_image",width=300,height=300)
                                input_file = gr.components.File(label='Upload fragments sdf file')
                with gr.Row():
                    with gr.Column(scale=1):
                        replace_outputs = gr.components.File(label='Replace molecules sdf file')
                with gr.Row():
                    replace = gr.Button(value='Replace Fragments',elem_id="btn8")
                with gr.Row():
                    example3 = gr.Image(label='Display replaced molecules')
                highlight.click(fn=highlightOneMolReplace, inputs=[raw_smi, replace_part_smi], outputs=highlightMol)
                replace.click(fn=replaceFrag, inputs=[raw_smi, replace_part_smi, input_file],outputs=[replace_outputs, example3])
            with gr.Tab(label="Link Fragments", elem_classes="body5_mini_2"):
                with gr.Row():
                    with gr.Column(scale=1):
                        input_file1=gr.File(label='Upload 2 or 3 fragments sdf file',file_count='multiple')
                    with gr.Column(scale=1):
                        output_file1 = gr.components.File(label='Link molecules sdf file')
                with gr.Row():
                    connect = gr.Button(value='Link Fragments',elem_id="btn9")
                with gr.Row():
                    example4 = gr.Image(label='Display linked molecules')
                connect.click(fn=connectFrag, inputs=[input_file1], outputs=[output_file1, example4])
            with gr.Tab(label="LinkInvent", elem_classes="body5_mini_3"):
                with gr.Row():
                    with gr.Column(scale=1):
                        warheads = gr.Textbox(label='Warheads File',placeholder='For example: *N[C@H](C)c1cc(N)cc(C(F)(F)F)c1|*C1CCOC1')
                        linkinvent = gr.Button(value='Linker Design',elem_id="btn10")
                    with gr.Column(scale=1):
                        output_file2 = gr.components.File(label='Generate molecules sdf file')
                with gr.Row():
                    example5=gr.Image(label='Display generated molecules')
                linkinvent.click(fn=Linkinvent,inputs=[warheads],outputs=[output_file2, example5])
        demo.launch(share=True)
