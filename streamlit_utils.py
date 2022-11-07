import streamlit as st

import py3Dmol
from stmol import showmol
from rdkit import Chem
from rdkit.Chem.AllChem import EmbedMolecule, MMFFOptimizeMolecule

from rdkit_utils import smi_to_svg, smi_to_canon_smi


def st_display_mol_with_smi(smiles, substruct='', use_smiles=False, display_substruct=True, img_size=(150, 100), labels=[]):
    if substruct != '' and display_substruct:
        if use_smiles:
            st.markdown(f'##### substruct SMILES: `{substruct.replace(")", "&rpar;")}`')
        else:
            st.markdown(f'##### substruct SMARTS: `{substruct.replace(")", "&rpar;")}`')
    columns = [col for col in st.columns(len(smiles))]
    if len(labels) == 0: labels = [''] * len(smiles)
    for smi, label, col in zip(smiles, labels, columns):
        col.image(smi_to_svg(smi, img_size=img_size, substruct=substruct, use_smiles=use_smiles))
        if label != '':
            col.markdown(f'&nbsp;&nbsp;**{label}**, {smi.replace(")", "&rpar;")}', unsafe_allow_html=True)
        else:
            col.markdown(f'&nbsp;&nbsp;{smi.replace(")", "&rpar;")}', unsafe_allow_html=True)

def st_match_smi_to_mol(smiles, img_size=(150, 100), labels=[], idx_offset=0):
    n_smi = len(smiles)
    columns = [col for col in st.columns(n_smi)]
    key_idxs = [i+idx_offset for i in range(n_smi)]
    if len(labels) == 0: labels = [''] * len(smiles)
    for smi, label, col, key_idx in zip(smiles, labels, columns, key_idxs):
        smi_key = f'smi_{key_idx}'
        button_key = f'button_match_smi_to_mol{key_idx}'
        col.text_input('SMILES:', key=smi_key)
        col.image(smi_to_svg(smi, img_size=img_size))
        if label != '':
            col.markdown(f'&nbsp;&nbsp;**{label}**', unsafe_allow_html=True)
        try: 
            if smi_to_canon_smi(st.session_state[smi_key]) == smi_to_canon_smi(smi):
                col.success('Correct!')
        except: pass
        if col.button('Answer', key=button_key): col.info(smi)

def st_check_smarts(smiles_in, smiles_out, answer, prompt, img_size=(150, 100), idx_offset=0):
    st.write('')
    st.markdown(f"""
        ##### {prompt}
    """)
    smt_key = f'smt_{idx_offset}'
    button_key = f'button_check_smarts_{idx_offset}'
    smt = st.text_input('SMARTS:', key=smt_key)
    substruct = Chem.MolFromSmarts(smt)
    st.markdown('***Match these molecules:***')
    st_display_mol_with_smi(smiles_in, substruct=smt, display_substruct=False, img_size=img_size)
    st.markdown("""***Don't match these molecules:***""")
    st_display_mol_with_smi(smiles_out, substruct=smt, display_substruct=False, img_size=img_size)
    if (
        sum([Chem.MolFromSmiles(smi).HasSubstructMatch(substruct) for smi in smiles_in]) == len(smiles_in) and
        sum([Chem.MolFromSmiles(smi).HasSubstructMatch(substruct) for smi in smiles_out]) == 0
    ):
        st.success('Correct!')
    if st.button('Answer', key=button_key): st.info(answer)

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
