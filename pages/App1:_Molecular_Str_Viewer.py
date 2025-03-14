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
    # Function to reset the input field and clear images
    def reset_input():
        st.session_state.smiles_input = ""

    # Column layout for reset button and text input
    col1, col2 = st.columns([1, 4])

    with col1:
        st.button("Reset", on_click=reset_input)

    with col2:
        # Using a default key for text_input to avoid conflicts
        smiles = st.text_input("Enter a SMILES string:", value=st.session_state.get('smiles_input', ''), key="smiles_input")

    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            AllChem.Compute2DCoords(mol)
            col1, col2 = st.columns([1, 1])
            with col1:
                st.image(Draw.MolToImage(mol), caption="2D Structure")
            with col2:
                # Placeholder for molecule name resolution
                st.write("Molecule details would appear here.")
        else:
            st.error("Invalid SMILES string. Please try again.")

    # Display test cases as clickable text for copying
    test_cases = {
        "SMILES String": [
            "C", "CCO", "CC(=O)OC1=CC=CC=C1C(=O)O", "C1=CC=CC=C1", "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "CC1=C2C(=CC=C1)C3C(OC2C4C3(OC(=O)C4=C)C5=C6C(=C7C(=C5)OCO7)C(=O)O6)OC(=O)CCO"
        ],
        "Molecule": ["Methane", "Ethanol", "Aspirin", "Benzene", "Ibuprofen", "Taxol"]
    }

    df = pd.DataFrame(test_cases)
    st.dataframe(df)  # Displays the dataframe with headers
