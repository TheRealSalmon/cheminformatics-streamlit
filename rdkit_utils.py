from rdkit import Chem
from rdkit.Chem.AllChem import ReactionFromSmarts
from rdkit.Chem.Draw.rdMolDraw2D import MolDrawOptions, MolDraw2DSVG, SetDarkMode


DRAW_OPTIONS = MolDrawOptions()
SetDarkMode(DRAW_OPTIONS)
DRAW_OPTIONS.addStereoAnnotation = True
DRAW_OPTIONS.atomLabelDeuteriumTritium = True
DRAW_OPTIONS.simplifiedStereoGroupLabel = True
DRAW_OPTIONS.fixedBondLength = 22
DRAW_OPTIONS.fixedFontSize = 14
DRAW_OPTIONS.annotationFontScale = 0.6
DRAW_OPTIONS.setBackgroundColour((0.05, 0.05, 0.05))


def mol_to_svg(mol, img_size=(500, 300), show_atom_idx=False, highlight_atoms=None):
    svg_drawer = MolDraw2DSVG(*img_size)
    svg_drawer.SetDrawOptions(DRAW_OPTIONS)
    if show_atom_idx:
        svg_drawer.drawOptions().addAtomIndices = True
    if highlight_atoms is None:
        svg_drawer.DrawMolecule(mol)
    else:
        color = (0.35, 0.3, 0.3)
        atom_colors = {atom_idx: color for atom_idx in range(mol.GetNumAtoms())}
        bond_colors = {bond_idx: color for bond_idx in range(mol.GetNumBonds())}
        highlight_bonds = [bond.GetIdx() for bond in mol.GetBonds() if bond.GetBeginAtomIdx() in highlight_atoms and bond.GetEndAtomIdx() in highlight_atoms]
        svg_drawer.DrawMolecule(
            mol,
            highlightAtoms=highlight_atoms,
            highlightBonds=highlight_bonds,
            highlightAtomColors=atom_colors,
            highlightBondColors=bond_colors
        )
    svg_drawer.FinishDrawing()
    return svg_drawer.GetDrawingText()

def smi_to_svg(smi, img_size=(500, 300), show_atom_idx=False, substruct='', use_smiles=False):
    mol = Chem.MolFromSmiles(smi)
    if substruct == '':
        return mol_to_svg(mol, img_size)
    if use_smiles:
        substruct_matches = mol.GetSubstructMatches(Chem.MolFromSmiles(substruct))
    else:
        substruct_matches = mol.GetSubstructMatches(Chem.MolFromSmarts(substruct))
    substruct_idxs = set()
    for substruct_match in substruct_matches:
        substruct_idxs.update(substruct_match)
    return mol_to_svg(mol, img_size, show_atom_idx, substruct_idxs)

def smirks_to_svg(
    smirks, img_size=(500, 300)
):
    svg_drawer = MolDraw2DSVG(*img_size)
    svg_drawer.SetDrawOptions(DRAW_OPTIONS)
    reaction = ReactionFromSmarts(smirks)
    svg_drawer.DrawReaction(reaction)
    svg_drawer.FinishDrawing()
    return svg_drawer.GetDrawingText()

def smi_to_canon_smi(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))

def run_reaction(smk, reactants):
    reaction = ReactionFromSmarts(smk)
    reactants = [Chem.MolFromSmiles(reactant) for reactant in reactants]
    outcomes = reaction.RunReactants(reactants)
    if len(outcomes) == 0:
        return [Chem.MolFromSmiles('')]
    elif len(outcomes) == 1:
        return outcomes[0]
    num_products = len(outcomes[0])
    unique_products = set()
    for outcome in outcomes:
        unique_products.update([Chem.MolToSmiles(product) for product in outcome])
    if len(unique_products) == num_products:
        return [Chem.MolFromSmiles(product) for product in unique_products]
    elif len(unique_products) > num_products:
        raise ValueError('obtained more than one unique product')
    elif len(unique_products) < num_products:
        raise ValueError('you should never see this error')
