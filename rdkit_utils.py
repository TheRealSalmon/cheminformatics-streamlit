from rdkit import Chem
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG, SetDarkMode


def mol_to_svg(mol, img_size=(500, 300), show_atom_idx=False):
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
    svg_drawer.DrawMolecule(mol)
    svg_drawer.FinishDrawing()
    return svg_drawer.GetDrawingText()

def smi_to_svg(smi, img_size=(500, 300), show_atom_idx=False):
    return mol_to_svg(Chem.MolFromSmiles(smi), img_size, show_atom_idx)

def smi_to_canon_smi(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))
