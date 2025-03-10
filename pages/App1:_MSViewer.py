import streamlit as st
import py3Dmol
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
from io import BytesIO

st.title("Molecule Viewer")

# Get SMILES input from user
smiles = st.text_input("Enter a SMILES string:")

if smiles:
    mol = Chem.MolFromSmiles(smiles)
    
    if mol:
        # Generate 2D structure
        img = Draw.MolToImage(mol)
        buf = BytesIO()
        img.save(buf, format="PNG")
        st.image(buf.getvalue(), caption="2D Structure")

        # 3D Model using py3Dmol
        xyz = Chem.MolToM
