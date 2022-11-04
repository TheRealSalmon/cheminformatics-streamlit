import streamlit as st

import py3Dmol
from stmol import showmol
from rdkit import Chem
from rdkit.Chem.AllChem import EmbedMolecule, MMFFOptimizeMolecule

from rdkit_utils import smi_to_svg, smi_to_canon_smi


def st_display_mol_with_smi(smiles, img_size=(100, 100), labels=[]):
    columns = [col for col in st.columns(len(smiles))]
    if len(labels) == 0: labels = [''] * len(smiles)
    for smi, label, col in zip(smiles, labels, columns):
        col.image(smi_to_svg(smi, img_size=img_size))
        if label != '':
            col.markdown(f'&nbsp;&nbsp;**{label}**, {smi}', unsafe_allow_html=True)
        else:
            col.markdown(f'&nbsp;&nbsp;{smi}', unsafe_allow_html=True)

def st_match_smi_to_mol(smiles, img_size=(100, 100), labels=[], idx_offset=0):
    n_smi = len(smiles)
    columns = [col for col in st.columns(n_smi)]
    keys = [f'smi{i+idx_offset}' for i in range(n_smi)]
    if len(labels) == 0: labels = [''] * len(smiles)
    for smi, label, col, key in zip(smiles, labels, columns, keys):
        col.text_input('SMILES:', key=key)
        col.image(smi_to_svg(smi, img_size=img_size))
        if label != '':
            col.markdown(f'&nbsp;&nbsp;**{label}**', unsafe_allow_html=True)
        try: 
            if smi_to_canon_smi(st.session_state[key]) == smi_to_canon_smi(smi):
                col.success('Correct!')
        except: pass

def st_display_3dmol(mol):
    mol = Chem.Mol(mol)
    mol = Chem.AddHs(mol)
    EmbedMolecule(mol)
    MMFFOptimizeMolecule(mol)
    view = py3Dmol.view(width=500, height=400)
    view.addModel(Chem.MolToMolBlock(mol), 'mol')
    view.setBackgroundColor('black')
    view.setStyle({'stick': {'colorscheme': 'pinkCarbon', 'radius': 0.25}})
    view.setViewStyle(
        {"style": "outline", "color": 'hotpink', "width": 0.04}
    )
    view.zoomTo()
    showmol(view)
