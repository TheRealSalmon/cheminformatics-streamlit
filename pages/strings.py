import streamlit as st
from PIL import Image

from streamlit_utils import st_display_mol_with_smi, st_match_smi_to_mol


def intro_page():
    st.markdown("""
        # Molecules as Strings

        The first place to start in cheminformatics is with a good attitude and
        a SMILES ... or a SMARTS or a SMIRKS? Maybe you're feeling a little
        InChI? And how do you make an InChI into an InChi key?

        There is a plethora of confusing and overlapping string representations
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

        Just like molecules, the most important building block of a SMILES is
        the atom. Atoms are defined by their atomic symbol.
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
        And there we go! We have the basic language needed to describe simple
        molecules.

        ## Linear molecules

        The easiest molecules to represent in SMILES are linear molecules. Here
        are some common solvents represented as SMILES:
    """)

    smiles = ['CCO', 'CCCCCC']
    labels = ['ethanol', 'n-hexane']
    st_display_mol_with_smi(smiles, img_size=(200,100), labels=labels)

    st.markdown("""
        SMILES don't respect "order" so `CCO` or `OCC` both represent ethanol.
        As you see, you don't need to place a `-` between atoms to place a
        single bond, a single bond will be placed by default. So `C-C-O` and
        `O-C-C` are the same.

        ### Your turn!
    """)

    smiles = ['CCOCC', 'COCCOC', 'CC#N']
    labels = ['diethyl ether', 'glyme', 'acetonitrile']
    st_match_smi_to_mol(smiles, labels=labels, idx_offset=5)

    st.markdown("""
        ## Branched molecules

        As we know, molecules are not always linear and so we need a way to
        represent branches in molecules. This is done by wrapping the branching
        atoms in parentheses.
    """)

    smiles = ['ClC(Cl)Cl', 'CN(C)C=O']
    labels = ['chloroform', 'DMF']
    st_display_mol_with_smi(smiles, img_size=(200,100), labels=labels)

    st.markdown("""
        While we represent chloroform as `ClC(Cl)Cl`, it would be just as valid
        to put `C` first as `C(Cl)(Cl)Cl`. Furthermore, you can put multiple
        atoms within parentheses, ie. DMF can be both `CN(C)C=O` or `CN(C=O)C`.

        ### Your turn!
    """)

    smiles = ['CS(=O)C', 'CCOC(=O)C']
    labels = ['DMSO', 'ethyl acetate']
    st_match_smi_to_mol(smiles, labels=labels, idx_offset=8)

    st.markdown("""
        ## Put a ring on it

        Lastly, we need to be able to represent rings. The way SMILES represents
        rings is to use numbers *outside* of brackets to denote connection
        points.
    """)

    smiles = ['CCC', 'C1CC1', '[CH-]1CC1', 'C1CC1C', 'C1C(C)C1']
    st_display_mol_with_smi(smiles, img_size=(100,100))

    st.markdown("""
        As you see, we "closed" the propane by putting `1` after the first and
        last C of the chain. As illustrated by the deprotonated cyclopropane,
        the `1` goes outside of the bracket. Lastly, rings can be combined with
        branching, either by adding atoms to the ring's closing atoms or by
        using parentheses.

        A special topic within rings is the proper representation of aromatic
        rings. Aromatic rings are unique in that there are two valid ways to
        represent them.
    """)

    smiles = ['C1=CC=CC=C1', 'c1ccccc1', 'N1=CNC=C1', 'n1c[nH]cc1']
    labels = ['kekulized', 'not kekulized'] * 2
    st_display_mol_with_smi(smiles, img_size=(150,100), labels=labels)

    st.markdown("""
        You will likely find that using the *not kekulized* form is more
        convenient as you don't need to remember where the double bonds in the
        ring are placed. However, if you have an implicit hydrogen on one of the
        ring atoms, you *will* need to specify it by replacing `n` with `[nH]`.
        Otherwise, you will get an error telling you that RDKit is unable to
        kekulize the molecule.

        Representing a molecule with two rings is done by using both `1` and
        `2`. Note that `C12` is interpreted as "C one" and "C two" rather than
        "C twelve". However this is pretty complicated and hopefully you rarely
        need to try and write bi/tricyclic SMILES by hand.
    """)
    smiles = ['C1CCCCN1c2ncccc2', 'c1cccc(cc[nH]2)c12', 'C12CC1C2']
    st_display_mol_with_smi(smiles, img_size=(150,100))

    st.markdown("""
        ### Your turn!
    """)

    smiles = ['C1CC=CCC1', 'N1C(=O)CC1', 'c1c(C)[nH]c(OC)c1']
    st_match_smi_to_mol(smiles, idx_offset=10)

    st.markdown("""
        ## Conclusion

        Congratulations on reaching the end! Hopefully this was a useful primer
        on SMILES strings. SMILES may not impress you quite yet, but the more
        you use them, the more you will appreciate the simplicity and power of
        representing complex molecules as "simple" 1D sequences of letters,
        numbers, and symbols. Without such a representation, cheminformatics
        would be far more complicated and difficult.
    """)


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
