from rdkit import Chem
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG, SetDarkMode


def mol_to_svg(mol, img_size=(500, 300), show_atom_idx=False, highlight_atoms=None):
    svg_drawer = MolDraw2DSVG(*img_size)
    SetDarkMode(svg_drawer.drawOptions())
    if show_atom_idx:
        svg_drawer.drawOptions().addAtomIndices = True
    svg_drawer.drawOptions().addStereoAnnotation = True
    svg_drawer.drawOptions().atomLabelDeuteriumTritium = True
    svg_drawer.drawOptions().simplifiedStereoGroupLabel = True
    svg_drawer.drawOptions().fixedBondLength = 22
    svg_drawer.drawOptions().fixedFontSize = 14
    svg_drawer.drawOptions().annotationFontScale = 0.6
    svg_drawer.drawOptions().setBackgroundColour((0.05, 0.05, 0.05))
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

def smi_to_canon_smi(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))
