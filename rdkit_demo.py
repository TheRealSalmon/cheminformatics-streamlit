import streamlit as st
import numpy as n
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image

from rdkit_utils import mol_to_svg


st.text_input('SMILES:', value='C', key='smi')
# mol = Chem.MolFromSmiles('C1CC1C(=O)NC2=CC=CC(=C2)[C@H](C)C3=NC=NC(=C3)NC4=CC=CC(=C4(C))C(F)(F)F')
mol = Chem.MolFromSmiles(st.session_state.smi)

im = mol_to_svg(mol)

st.image(im)

if st.checkbox('Show dataframe'):
    chart_data = pd.DataFrame(
       np.random.randn(20, 3),
       columns=['a', 'b', 'c'])

    chart_data
