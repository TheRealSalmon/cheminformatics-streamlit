import streamlit as st
from PIL import Image
from rdkit import Chem
from rdkit.Chem.AllChem import EmbedMolecule, MMFFOptimizeMolecule
from stmol import showmol
# from streamlit_ace import st_ace
import py3Dmol

from streamlit_utils import st_display_3dmol
from rdkit_utils import mol_to_svg

st.set_page_config(
    page_title="home",
    # page_icon=img,
    layout="wide",
)

st.markdown("""
    # Interactive Cheminformatics

    If you're here, you're probably already interested in cheminformatics.
    But it can be rather hard to learn. Tutorials are spread out all over
    the place and documentation can be incredibly hard to read.

    With this webapp I will help you teach yourself how to use
    computers to represent molecules, reactions, and even your chemical
    intuition.

    Here is an example of what I mean by "interactive". Give it a try! Edit the
    SMILES or try writing your own.
""")

LEFT, RIGHT = st.columns(2)

LEFT.text_input(
    'SMILES:',
    value='CC(C)(C#CC1=NC(=C(C=C1)C2=C3C(=C(C=C2)Cl)C(=NN3CC(F)(F)F)NS(=O)(=O)C)C(CC4=CC(=CC(=C4)F)F)NC(=O)CN5C6=C(C7CC7C6(F)F)C(=N5)C(F)(F)F)S(=O)(=O)C',
    key='smi'
)
try:
    mol = Chem.MolFromSmiles(st.session_state.smi)
    LEFT.image(mol_to_svg(mol))
    with RIGHT: st_display_3dmol(mol)
except:
    RIGHT.markdown('invalid molecule')

# code = st_ace(
#     language='python',
#     theme='twilight',
#     font_size=12,
# )
 