import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import py3Dmol
import pandas as pd
import tempfile

# Initialize session state variables if not already set
if 'reset' not in st.session_state:
    st.session_state['reset'] = False
if 'smiles' not in st.session_state:
    st.session_state['smiles'] = ""

# Function to clear input and reset images
def reset_form():
    st.session_state.smiles = ""
    st.session_state.reset = True

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
    col1, col2 = st.columns([1, 4])
    with col1:
        st.button("Reset", on_click=reset_form)
    with col2:
        st.session_state.smiles = st.text_input("Enter a SMILES string:", key="smiles", value=st.session_state.smiles)

    smiles = st.session_state.smiles
    if smiles and not st.session_state.reset:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Generate 2D coordinates for visualization
            AllChem.Compute2DCoords(mol)

            # Create two columns to display image and name side-by-side
            col1, col2 = st.columns([1, 1])

            with col1:
                st.image(Draw.MolToImage(mol), caption="2D Structure")

            with col2:
                molecule_name = None
                for i, row in enumerate([
                    "C", "CCO", "CC(=O)OC1=CC=CC=C1C(=O)O", 
                    "C1=CC=CC=C1", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", 
                    "CC1=C2C(=CC=C1)C3C(OC2C4C3(OC(=O)C4=C)C5=C6C(=C7C(=C5)OCO7)C(=O)O6)OC(=O)CCO"
                ]):
                    if smiles == row:
                        molecule_name = [
                            "Methane", "Ethanol", "Aspirin", 
                            "Benzene", "Ibuprofen", "Taxol"
                        ][i]
                        break
                if molecule_name:
                    st.markdown(f"<h2 style='color:#4CAF50; font-size: 28px;'>{molecule_name}</h2>", unsafe_allow_html=True)
                else:
                    st.markdown(f"<h2 style='color:#4CAF50; font-size: 28px;'>Unknown Compound</h2>", unsafe_allow_html=True)

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
    elif st.session_state.reset:
        st.session_state.reset = False

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

    st.markdown("### 🧪 Test Cases")
    df = pd.DataFrame(test_cases)
    st.dataframe(df)
