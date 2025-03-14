import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import py3Dmol
import pandas as pd
import tempfile

# Create tabs
tabs = st.tabs(["🔎 App overview", "🐾 Practice with APP"])

with tabs[0]:
    st.markdown("""
    #### Purpose:
    + To generate and visualize the 2D and 3D structures of chemical compounds using SMILES notation.  
   
    #### Key Features:
   - 🧪 **SMILES Input** – Allow users to input a SMILES string to generate molecular structures.
   - 🖼️ **2D Visualization** – Display a 2D structural diagram of the molecule.
   - ✅ **Structure Validation** – Provide feedback if the input SMILES is invalid.
    """)

with tabs[1]:
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
            mol_block = Chem.MolToMolBlock(mol)

            # Create 3D viewer using py3Dmol
            viewer = py3Dmol.view(width=400, height=400)
            viewer.addModel(mol_block, "sdf")
            viewer.setStyle({"stick": {}})
            viewer.zoomTo()

            # Use write_html to generate complete HTML
            with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmpfile:
                viewer.write_html(tmpfile.name)
                html_file = tmpfile.name

            # Read HTML file and display it in Streamlit
            with open(html_file, "r") as file:
                html_content = file.read()

            st.components.v1.html(html_content, height=400, width=400)

        else:
            st.error("Invalid SMILES string. Please try again.")

    # Test Cases Table
    test_cases = {
        "SMILES String": [
            "C", 
            "CCO", 
            "CC(=O)OC1=CC=CC=C1C(=O)O", 
            "C1=CC=CC=C1", 
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "CC1=C2C(=CC=C1)C3C(OC2C4C3(OC(=O)C4=C)C5=C6C(=C7C(=C5)OCO7)C(=O)O6)OC(=O)CCO"
        ],
        "Molecule": [
            "Methane", 
            "Ethanol", 
            "Aspirin", 
            "Benzene", 
            "Ibuprofen",
            "Taxol"
        ]
    }
    
    # Convert to pandas DataFrame and display
    df = pd.DataFrame(test_cases)
    st.markdown("### 🧪 Test Cases")
    st.table(df)
