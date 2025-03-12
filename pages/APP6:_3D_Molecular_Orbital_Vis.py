import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
import json

def generate_3d_coordinates(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    return Chem.MolToMolBlock(mol)

def create_viewer(mol_block):
    viewer_html = f"""
    <div id="molView" style="height: 400px; width: 600px;"></div>
    <script src="https://3Dmol.org/build/3Dmol.js"></script>
    <script>
    let viewer = new $3Dmol.createViewer(document.getElementById('molView'), {{
        backgroundColor: 'white'
    }});
    viewer.addModel(`{mol_block}`, 'sdf');
    viewer.setStyle({{}}, {{stick: {{}}}});
    viewer.zoomTo();
    viewer.render();
    </script>
    """
    return viewer_html

def main():
    st.title('3D Molecule Viewer')
    smiles_input = st.text_input("Enter a SMILES string:", "CCO")  # Example for ethanol

    if st.button("Visualize Molecule"):
        mol_block = generate_3d_coordinates(smiles_input)
        viewer_html = create_viewer(mol_block)
        st.components.v1.html(viewer_html, height=420)

if __name__ == "__main__":
    main()
