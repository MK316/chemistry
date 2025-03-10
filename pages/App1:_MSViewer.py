import streamlit as st
import py3Dmol
from rdkit import Chem
from rdkit.Chem import Draw

# Get SMILES input from user
smiles = st.text_input("Enter a SMILES string:")

if smiles:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        st.image(Draw.MolToImage(mol), caption="2D Structure")
        
        # 3D Model using py3Dmol
        xyz = Chem.MolToMolBlock(mol)
        viewer = py3Dmol.view(width=400, height=400)
        viewer.addModel(xyz, "sdf")
        viewer.setStyle({"stick": {}})
        viewer.zoomTo()
        viewer.show()
