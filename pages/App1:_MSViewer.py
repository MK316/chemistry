import streamlit as st
import py3Dmol
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

# Get SMILES input from user
smiles = st.text_input("Enter a SMILES string:")

if smiles:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        # Generate 2D coordinates for visualization
        AllChem.Compute2DCoords(mol)
        st.image(Draw.MolToImage(mol), caption="2D Structure")
        
        # Generate 3D coordinates for visualization
        AllChem.EmbedMolecule(mol)
        xyz = Chem.MolToMolBlock(mol)

        # Display 3D model using py3Dmol
        viewer = py3Dmol.view(width=400, height=400)
        viewer.addModel(xyz, "sdf")
        viewer.setStyle({"stick": {}})
        viewer.zoomTo()
        st.components.v1.html(viewer.show(), height=400, width=400)
    else:
        st.error("Invalid SMILES string. Please try again.")
