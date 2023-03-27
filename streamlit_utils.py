import py3Dmol
import streamlit as st
from rdkit import Chem
from rdkit.Chem.AllChem import EmbedMolecule, MMFFOptimizeMolecule
from stmol import showmol

from rdkit_utils import (
    run_reaction,
    smi_to_canon_smi,
    smi_to_svg,
    smirks_to_svg,
)


def st_display_mol_with_smi(
    smiles,
    display_smiles=False,
    substruct="",
    use_smiles=False,
    display_substruct=True,
    img_size=(150, 100),
    labels=[],
):
    if substruct != "" and display_substruct:
        if use_smiles:
            st.markdown(f"substruct SMILES: `{substruct}`")
        else:
            st.markdown(f"substruct SMARTS: `{substruct}`")
    n_cols = 4 if len(smiles) < 4 else len(smiles)
    columns = [col for col in st.columns(n_cols)]
    if len(labels) == 0:
        labels = [""] * len(smiles)
    for smi, label, col in zip(smiles, labels, columns):
        col.image(
            smi_to_svg(
                smi,
                img_size=img_size,
                substruct=substruct,
                use_smiles=use_smiles,
            )
        )
        if label != "":
            col.markdown(f"&nbsp;&nbsp;**{label}**")
        if display_smiles:
            col.markdown(f"&nbsp;&nbsp;`{smi}`")
        if substruct != "":
            mol = Chem.MolFromSmiles(smi)
            if use_smiles:
                substruct_mol = Chem.MolFromSmiles(substruct)
            else:
                substruct_mol = Chem.MolFromSmarts(substruct)
            if mol.HasSubstructMatch(substruct_mol):
                col.markdown("*Match!*")
            else:
                col.markdown("*No match!*")


def st_match_smi_to_mol(smiles, img_size=(150, 100), labels=[], idx_offset=0):
    n_smi = len(smiles)
    n_cols = 4 if len(smiles) < 4 else len(smiles)
    columns = [col for col in st.columns(n_cols)]
    key_idxs = [i + idx_offset for i in range(n_smi)]
    if len(labels) == 0:
        labels = [""] * len(smiles)
    for smi, label, col, key_idx in zip(smiles, labels, columns, key_idxs):
        smi_key = f"smi_{key_idx}"
        button_key = f"button_match_smi_to_mol{key_idx}"
        col.text_input("SMILES:", key=smi_key)
        col.image(smi_to_svg(smi, img_size=img_size))
        if label != "":
            col.markdown(f"&nbsp;&nbsp;**{label}**")
        try:
            if smi_to_canon_smi(st.session_state[smi_key]) == smi_to_canon_smi(
                smi
            ):
                col.success("Correct!")
        except ValueError:
            pass
        if col.button("Answer", key=button_key):
            col.info(f"`{smi}`")


def st_check_smarts(
    smiles_in, smiles_out, answer, prompt, img_size=(150, 100), idx_offset=0
):
    st.write("")
    st.markdown(
        f"""
        ##### {prompt}
    """
    )
    smt_key = f"smt_{idx_offset}"
    button_key = f"button_check_smarts_{idx_offset}"
    smt = st.text_input("SMARTS:", key=smt_key)
    substruct = Chem.MolFromSmarts(smt)
    st.markdown("***Match these molecules:***")
    st_display_mol_with_smi(
        smiles_in, substruct=smt, display_substruct=False, img_size=img_size
    )
    st.markdown("""***Don't match these molecules:***""")
    st_display_mol_with_smi(
        smiles_out, substruct=smt, display_substruct=False, img_size=img_size
    )
    if (
        sum(
            [
                Chem.MolFromSmiles(smi).HasSubstructMatch(substruct)
                for smi in smiles_in
            ]
        )
        == len(smiles_in)
        and sum(
            [
                Chem.MolFromSmiles(smi).HasSubstructMatch(substruct)
                for smi in smiles_out
            ]
        )
        == 0
    ):
        st.success("Correct!")
    if st.button("Answer", key=button_key):
        st.info(f"`{answer}`")


def st_display_reaction(
    reactants=[],
    reagents=[],
    products=[],
    reactants_are_smarts=False,
    reagents_are_smarts=False,
    products_are_smarts=False,
    img_size=(400, 100),
):
    def smi_to_smt(smi):
        return Chem.MolToSmarts(Chem.MolFromSmiles(smi))

    if not reactants_are_smarts:
        reactants = [smi_to_smt(smi) for smi in reactants]
    if not reagents_are_smarts:
        reagents = [smi_to_smt(smi) for smi in reagents]
    if not products_are_smarts:
        products = [smi_to_smt(smi) for smi in products]
    smirks = f'{".".join(reactants)}>{".".join(reagents)}>{".".join(products)}'
    st.image(smirks_to_svg(smirks, img_size=img_size))


def st_run_reaction(smk, reactants, display_smk=True, img_size=(400, 100)):
    if display_smk:
        st.markdown(f"SMIRKS: `{smk}`")
    products = run_reaction(smk, reactants)
    if [Chem.MolToSmiles(product) for product in products] == [""]:
        st_display_reaction(reactants=reactants, products=["*"])
        st.error("no reaction")
    else:
        st_display_reaction(
            reactants=reactants,
            products=[Chem.MolToSmarts(product) for product in products],
            products_are_smarts=True,
        )


def st_check_smirks(
    reactantss, productss, answer, prompt, img_size=(400, 100), idx_offset=0
):
    st.write("")
    st.markdown(
        f"""
        ##### {prompt}
    """
    )
    n_schemes = len(reactantss)
    smk_key = f"smk_{idx_offset}"
    button_key = f"button_check_smirks_{idx_offset}"
    smk = st.text_input("SMIRKS:", key=smk_key)
    for i, reactants_products in enumerate(zip(reactantss, productss)):
        reactants, products = reactants_products
        if products == [""]:
            st_display_reaction(
                reactants=reactants, products="*", products_are_smarts=True
            )
            st.error("no reaction")
        else:
            st_display_reaction(reactants=reactants, products=products)
        if smk != "":
            st.markdown("***your reaction:***")
            st_run_reaction(smk, reactants, display_smk=False)
            products = [smi_to_canon_smi(product) for product in products]
            output_products = run_reaction(smk, reactants)
            if sorted(products) == sorted(
                [
                    Chem.MolToSmiles(output_product)
                    for output_product in output_products
                ]
            ):
                st.success("Correct!")
        if i < n_schemes - 1:
            st.markdown("---")
    if st.button("Answer:", key=button_key):
        st.info(answer)


def st_display_3dmol(mol):
    mol = Chem.Mol(mol)
    mol = Chem.AddHs(mol)
    EmbedMolecule(mol)
    MMFFOptimizeMolecule(mol)
    view = py3Dmol.view(width=500, height=400)
    view.addModel(Chem.MolToMolBlock(mol), "mol")
    view.setBackgroundColor("black")
    view.setStyle({"stick": {"colorscheme": "pinkCarbon", "radius": 0.25}})
    view.setViewStyle({"style": "outline", "color": "hotpink", "width": 0.04})
    view.zoomTo()
    showmol(view)
