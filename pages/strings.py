import streamlit as st
from PIL import Image

from streamlit_utils import st_display_mol_with_smi, st_match_smi_to_mol


def intro_page():
    st.markdown("""
        # Molecules as Strings

        The first place to start in cheminformatics is with a good attitude and
        a SMILES ... or a SMARTS or a SMIRKS? Maybe you're feeling a little
        InChI? And how do you make an InChI into an InChi key?

        There is a plethora of confusing and overlapping string representation
        that all play their own role in cheminformatics. But no worries, we'll
        take them on one-by-one.

        Navigate to the sidebar and select a *sub-topic* to continue.
    """)

def smiles_page():
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
    st_display_mol_with_smi(smiles)

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
    st_display_mol_with_smi(smiles)

    st.markdown("""
        As we see, specifying atom-level properties requires the use of
        brackets.
        * The isotope is specified by the numbers that PRECEDE the atomic symbol, ie. `[13C]`.
        * The # of hydrogens is specified by adding `H` to the end of the atomic symbol and putting the # of hydrogens after the `H`, ie. `[CH3]`.
        * The formal charge is specified by adding a `+` or `-` after the atom symbol/# of hydrogens. For multiple charges, an additional number can be put after the charge, ie. `[Fe+2]`.

        ### Your turn!

        Write a SMILES string that matches the molecule depicted below.
        """)

    smiles = ['N', '[H]', '[He]', '[OH-]', '[18OH2]']
    st_match_smi_to_mol(smiles)

    st.markdown("""
        ## Bonds

        After learning to represent atoms, next we come to bonds. Bonds are far
        simpler than atoms as they come only in 4 categories:
        * single
        * double
        * triple
        * aromatic
    """)

    smiles = ['C-C', 'C=C', 'C#C', 'C:C']
    st_display_mol_with_smi(smiles)

    st.markdown("""
        And there we go! We know enough now to take a few steps forward...

        ## Linear molecules

        The easiest molecules to represent in SMILES are linear molecules. We
        will cover branching in the next section but let's make sure we've
        covered the basics. Here are some common solvents represented as SMILES:
    """)

    smiles = ['CCO', 'CCCCCC']
    labels = ['ethanol', 'n-hexane']
    st_display_mol_with_smi(smiles, img_size=(200,100), labels=labels)

    st.markdown("""
        ### Your turn!
    """)

    smiles = ['CCOCC', 'CC#N', 'COCCOC']
    labels = ['diethyl ether', 'acetonitrile', 'glyme']
    st_match_smi_to_mol(smiles, labels=labels, idx_offset=5)

subpages = [
    ('0 - intro', intro_page),
    ('1 - SMILES', smiles_page),
]

st.sidebar.markdown('## sub-topics')
subpage = st.sidebar.selectbox(
    'select page',
    subpages,
    format_func=lambda x: x[0]
)
st.sidebar.markdown("---")
subpage[1]()
