# from streamlit_ace import st_ace
import streamlit as st
from rdkit import Chem

from rdkit_utils import mol_to_svg
from streamlit_utils import st_display_3dmol

st.set_page_config(
    page_title="home",
    # page_icon=img,
    layout="wide",
)

st.markdown(
    """
    # Interactive Cheminformatics

    **What is cheminformatics?** To me, cheminformatics is the magic of
    expressing a chemical intuition in a computer-readable format.

    So much work has been done in cheminformatics, including:
    * how can we represent molecules?
    * how can we quantify the "similarity" of molecules?
    * how do we estimate molecular properties?

    However, despite the exciting progress in cheminformatics, there is a
    relative shortage of high-quality educational resources. Tutorials are
    spread out all over the place, documentation can be opaque, and there
    certainly isn't enough practice materials for beginners.

    ---

    With this webapp I will introduce you to the field of cheminformatics, one
    concept at a time. This curriculum, however, is entirely interactive and you
    are meant and learn and "play" your way through these topics!

    To get started, open the left side-bar and start working through a topic
    that interests you! Right now we have:
    * strings: learn how molecules, substructures, and reactions can be encoded
    in a sequence of letters, numbers, and symbols

    I plan to regularly release fresh content until I feel that I've adequately
    covered the major topics in modern cheminformatics.

    ---

    Here is an example of what I mean by "interactive". Give it a try! Edit the
    SMILES or try writing your own.
"""
)

LEFT, RIGHT = st.columns(2)

LEFT.text_input(
    "SMILES:",
    value="CC(C)(C#CC1=NC(=C(C=C1)C2=C3C(=C(C=C2)Cl)C(=NN3CC(F)(F)F)NS(=O)(=O)C"
    ")C(CC4=CC(=CC(=C4)F)F)NC(=O)CN5C6=C(C7CC7C6(F)F)C(=N5)C(F)(F)F)S(=O)(=O)C",
    key="smi",
)
try:
    mol = Chem.MolFromSmiles(st.session_state.smi)
    LEFT.image(mol_to_svg(mol))
    with RIGHT:
        st_display_3dmol(mol)
except ValueError:
    RIGHT.markdown("invalid molecule")

# code = st_ace(
#     language='python',
#     theme='twilight',
#     font_size=12,
# )
