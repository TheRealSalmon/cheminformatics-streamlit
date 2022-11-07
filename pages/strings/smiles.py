import streamlit as st

from streamlit_utils import (
    st_display_mol_with_smi,
    st_match_smi_to_mol,
)


def smiles_page():
    IDX_OFFSET = 0

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
    st_display_mol_with_smi(smiles, display_smiles=True)

    st.markdown("""
        In general, just the atomic symbol is sufficient. However, for hydrogen
        and elements with two-letter atomic symbols, brackets `[H]` are
        required.

        Atomic symbols usually have their first letter capitalized. However, the
        first letter may be lowercase if it is aromatic like this benzene.
    """)

    st.info("""
        Hydrogen is unique as the only one-letter atom that cannot exist outside
        of a bracket. This is because SMILES treats hydrogen as a second class
        citizen and it usually exists as an accessory to a heavy atom.
    """)

    st.markdown("""
        SMILES also allows us to specify a few atom-level properties such as:
        * isotope
        * &#8203;# of hydrogens
        * formal charge
    """)

    smiles = ['C', '[13C]', '[13CH4]', '[CH3+]', '[Fe+2]',]
    st_display_mol_with_smi(smiles, display_smiles=True)

    st.markdown("""
        As we see, specifying atom-level properties requires the use of
        brackets.
        * The isotope (superscript numbers) is specified by the numbers that PRECEDE the atomic symbol, ie. `[13C]`.
        * The # of hydrogens is specified by adding `H#` to the end of the atomic symbol and where # is the number of hydrogens, ie. `[CH3]`.
        * The formal charge is specified by adding a `+` or `-` after the atomic symbol/# of hydrogens. For multiple charges, an additional number can be put after the charge, ie. `[Fe+2]`.
    """)

    st.info("""
        When an atom is specified without brackets, Hs are added implicitly. So
        `C` gets translated to `[CH4]` automatically. However, once we add
        brackets, we take complete control of the SMILES string and so we need
        to manually add the Hs ourselves.
    """)

    st.markdown("""
        ### Your turn!

        Write a SMILES string that matches the molecule depicted below.
        """)

    st.info("""
        Don't forget that SMILES are CaPs-SeNsItIvE. Uppercase means aliphatic/
        non-aromatic while lowercase means aromatic.
    """)

    smiles = ['N', '[H]', '[He]', '[OH-]', '[18OH2]']
    st_match_smi_to_mol(smiles, img_size=(100,100), idx_offset=IDX_OFFSET)
    IDX_OFFSET += len(smiles)

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
    st_display_mol_with_smi(smiles, display_smiles=True)

    st.markdown("""
        ### Your turn!
    """)

    smiles = ['C=O', 'C#N', 'CC', 'C=C=C']
    st_match_smi_to_mol(smiles, img_size=(125,100), idx_offset=IDX_OFFSET)
    IDX_OFFSET += len(smiles)

    st.markdown("""
        And there we go! We have the basic language needed to describe simple
        molecules.
    """)

    st.info("""
        You may notice that `CC` and `C-C` are the same thing. This is because
        SMILES will automatically guess that the missing bond between `C` and
        `C` is a single bond. This is the same for `cc` and `c:c`. SMILES will
        guess that the bond between the two aromatic atoms is an aromatic bond.
    """)

    st.markdown("""
        ## Linear molecules

        The easiest molecules to represent in SMILES are linear molecules as
        they can be represented as a simple sequence of atoms and bonds. Here
        are some common solvents represented as SMILES:
    """)

    smiles = ['CCO', 'CCCCCC']
    labels = ['ethanol', 'n-hexane']
    st_display_mol_with_smi(smiles, display_smiles=True, img_size=(200,100), labels=labels)

    st.info("""SMILES don't respect "order" so `CCO` or `OCC` are equivalent.""")

    st.markdown("""
        ### Your turn!
    """)

    smiles = ['CCOCC', 'COCCOC', 'CC#N']
    labels = ['diethyl ether', 'glyme', 'acetonitrile']
    st_match_smi_to_mol(smiles, labels=labels, idx_offset=IDX_OFFSET)
    IDX_OFFSET += len(smiles)

    st.markdown("""
        ## Branched molecules

        As we know, molecules are not always linear and so we need a way to
        represent branches in molecules. This is done by wrapping the branching
        atom(s) in parentheses.
    """)

    smiles = ['ClC(Cl)Cl', 'CN(C)C=O']
    labels = ['chloroform', 'DMF']
    st_display_mol_with_smi(smiles, display_smiles=True, img_size=(200,100), labels=labels)

    st.markdown("""
        While we represent chloroform as `ClC(Cl)Cl`, it would be just as valid
        to put `C` first as `C(Cl)(Cl)Cl`. Furthermore, you can put multiple
        atoms within parentheses, ie. DMF can be both `CN(C)C=O` or `CN(C=O)C`.

        ### Your turn!
    """)

    smiles = ['CS(=O)C', 'CCOC(=O)C']
    labels = ['DMSO', 'ethyl acetate']
    st_match_smi_to_mol(smiles, labels=labels, idx_offset=IDX_OFFSET)
    IDX_OFFSET += len(smiles)

    st.markdown("""
        ## Put a ring on it

        Lastly, we need to be able to represent rings. The way SMILES represents
        rings is to use numbers *outside* of brackets to denote connection
        points.
    """)

    smiles = ['CCC', 'C1CC1', '[CH-]1CC1', 'C1CC1C', 'C1C(C)C1']
    st_display_mol_with_smi(smiles, display_smiles=True, img_size=(100,100))

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
    st_display_mol_with_smi(smiles, display_smiles=True, img_size=(150,100), labels=labels)

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
    st_display_mol_with_smi(smiles, display_smiles=True, img_size=(150,100))

    st.markdown("""
        ### Your turn!
    """)

    smiles = ['C1CC=CCC1', 'N1C(=O)CC1', 'c1c(C)[nH]c(OC)c1']
    st_match_smi_to_mol(smiles, idx_offset=IDX_OFFSET)
    IDX_OFFSET += len(smiles)

    st.markdown("""
        ## Chirality

        The last topic to be covered with SMILES is chirality. If the chiral
        center has a hydrogen, you must specify it as well with `@H` or `@@H`.
        If the chiral center has no hydrogens, you specify chirality with just a
        `@` or `@@`. 
    """)

    smiles = ['CC(F)CC', 'C[C@H](F)CC', 'C[C@@H](F)CC', 'C[C@](F)(Cl)CC', 'C[C@@](Cl)(F)CC']
    st_display_mol_with_smi(smiles, display_smiles=True, img_size=(150,100))

    st.markdown("""
        ### Your turn!

        This one is extra tough, don't worry too much about getting it right but
        feel free to give it a try.
    """)

    smiles = ['OC(=O)[C@H](C)N', 'OC(=O)[C@@H](C)N', 'OC(=O)[C@H](CC(C)C)N']
    st_match_smi_to_mol(smiles, idx_offset=IDX_OFFSET)
    IDX_OFFSET += len(smiles)

    st.markdown("""
        ## Conclusion

        Congratulations on reaching the end! Hopefully this was a useful primer
        on SMILES strings. SMILES may not impress you quite yet, but the more
        you use them, the more you will appreciate the simplicity and power of
        representing complex molecules as "simple" 1D sequences of letters,
        numbers, and symbols. Without such a representation, cheminformatics
        would be far more complicated and difficult.
    """)