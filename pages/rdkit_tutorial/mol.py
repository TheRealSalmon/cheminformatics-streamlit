import streamlit as st

from streamlit_utils import st_display_mol_with_smi


def mol_page():
    st.markdown(
        """
        # RDKit Mol

        In RDKit, the Mol is a c++ graph object wrapped in some Python. First,
        we need to import the package.
        ```python
        from rdkit import Chem
        ```
        Next, you can create an RDKit Mol. This is often done in one of a few
        ways.
        ```python
        mol = Chem.MolFromSmiles(smi_str)
        mol = Chem.MolFromMolFile(mdl_file_path, removeHs=False)
        ```
        You can also read a list of Mols from an SD File.
        ```python
        with Chem.SDMolSupplier(sd_file_path, removeHs=False) as sd_fh:
            mols = [mol for mol in sd_fh if mol is not None]
        ```
        Congrats! You now have an RDKit Mol. If you are in a Jupyter Notebook,
        you will find it useful to import the IPythonConsole module so that Mol
        objects display as images.
        ```python
        from rdkit.Chem.Draw import IPythonConsole
        ```
        """
    )

    st_display_mol_with_smi(["n1c(C)cccc1"])

    st.markdown(
        """
        ## Common operations.
        1. Add and remove hydrogens.
        ```python
        mol_no_hs = Chem.RemoveHs(mol)
        mol_with_hs = Chem.AddHs(mol)
        ```
        2. Substructure searches.
        ```python
        smt = Chem.MolFromSmarts()'[CH3]')
        mol = Chem.MolFromSmiles('CCC')

        assert mol.HasSubstructMatch(smt)
        assert mol.GetSubstructMatches(smt) == ((0,), (2,),)
        ```
        3. Calculate descriptors.
        ```python
        from rdkit.Chem.Descriptors import MolWt, MolLogP
        from rdkit.Chem.rdMolDescriptors import (
            CalcTPSA,
            CalcNumRotatableBonds,
            CalcNumLipinskiHBA,
        )

        n_ha = mol.GetNumHeavyAtoms()
        mw = MolWt(mol)
        clogp = MolLogP(mol)
        tpsa = CalcTPSA(mol)
        n_rot_bonds = CalcNumRotatableBonds(mol)
        n_hba = CalcNumLipinskiHBA(mol)
        ```
        4. Add conformer(s). These operations are done in-place.
        ```python
        from rdkit.Chem.AllChem import EmbedMolecule, EmbedMultipleConfs

        # conf = mol.GetConformer()  # this will fail
        EmbedMolecule(
            mol,
            useSmallRingTorsions=True,
            useMacrocycleTorsions=True,
            ETversion=2,
        )
        conf = mol.GetConformer()  # this will now succeed

        EmbedMultipleConfs(
            mol,
            numConfs=100,
            pruneRmsThresh=0.5,
            numThreads=1,
            useSmallRingTorsions=True,
            useMacrocycleTorsions=True,
            ETversion=2,
        )
        ```

        ## datamol
        RDKit has a steep learning curve and it can take a while until you feel
        comfortable using it. That's where datamol comes in. `datamol` is what
        happens when you adapt RDKit to fit modern (early 2020s) Python
        paradigms. I would still highly recommend learning RDKit, as datamol is
        just a wrapper for RDKit, but learning datamol first maybe helpful.

        ## Conclusions
        This is just the basics of RDKit Mols. You'll learn more by playing
        around with them rather than reading more. Good luck!
        """
    )
