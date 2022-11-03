import streamlit as st

import pandas as pd

# from rdkit_utils import smi_to_svg

def page():
    st.markdown("""
        # SMILES

        SMILES are a Simplified Molecular Input Line Entry System. In other
        words, it's a clever way to represent 2D molecules as a 1D
        string of letters and numbers.

        We will start with an overview of the syntax, and then move on to
        examples of greater and greater complexity.

        ## Atoms

        As with molecules, the most important building block of a SMILES is the
        atom. Atoms are defined by their atomic symbol.
    """)

    smiles = ['C', 'O', '[H]', '[Zn]', 'c1ccccc1']
    columns = [col for col in st.columns(5)]
    for smi, col in zip(smiles, columns):
        col.image(smi_to_svg(smi, img_size=(100, 100)))
        col.markdown(f'&nbsp;&nbsp;{smi}', unsafe_allow_html=True)

    st.markdown("""
        In general, just the atomic symbol is sufficient. However, for hydrogen
        and elements with two-letter atomic symbols, brackets `[H]` are
        required.

        Atomic symbols usually have their first letter capitalized. However, the
        first letter may be lowercase if it is aromatic like this benzene.

        SMILES also allows us to specify a few atom-level properties such as:
        * isotope
        * &#8203;# of hydrogens
        * formal charge
    """)

    smiles = ['C', '[13C]', '[13CH4]', '[CH3+]', '[Fe+2]',]
    columns = [col for col in st.columns(5)]
    for smi, col in zip(smiles, columns):
        col.image(smi_to_svg(smi, img_size=(100, 100)))
        col.markdown(f'&nbsp;&nbsp;{smi}', unsafe_allow_html=True)

    st.markdown("""
        As we see, specifying atom-level properties requires the use of
        brackets.
        * The isotope is specified by the numbers that PRECEDE the atomic symbol, ie. `[13C]`.
        * The # of hydrogens is specified by adding `H` to the end of the atomic symbol and putting the # of hydrogens after the `H`, ie. `[CH3]`.
        * The formal charge is specified by adding a `+` or `-` after the atom symbol/# of hydrogens. For multiple charges, an additional number can be put after the charge, ie. `[Fe+2]`.

        ### Your turn!
        """)

    smiles = ['N', '[H]', '[He]', '[OH-]', '[18OH2]']
    columns = [col for col in st.columns(5)]
    keys = [f'smi{i}' for i in range(5)]
    for smi, col, key in zip(smiles, columns, keys):
        col.text_input
        col.write(smi) 
