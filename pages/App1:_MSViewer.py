import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import tempfile
import pandas as pd

tabs = st.tabs(["🍐 App overview", "🍏 Practice with APP"])

with tabs[0]:
st.caption("Overview")

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
    
            # Write to temporary file for 3D visualization
            with tempfile.NamedTemporaryFile(delete=False, suffix=".sdf") as tmpfile:
                tmpfile.write(mol_block.encode())
                sdf_file = tmpfile.name
            
            # HTML block to render with 3Dmol.js
            html_block = f"""
            <div style="width: 400px; height: 400px;">
                <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
                <script>
                    let viewer = $3Dmol.createViewer($("#3dmoldiv"), {{
                        width: 400,
                        height: 400,
                        backgroundColor: "white"
                    }});
                    viewer.addModel(`{mol_block}`, "sdf");
                    viewer.setStyle({{"stick": {{}}}});
                    viewer.zoomTo();
                    viewer.render();
                </script>
            </div>
            """
    
            # Display HTML using Streamlit components
            st.components.v1.html(html_block, height=400, width=400)
    
        else:
            st.error("Invalid SMILES string. Please try again.")
    
    # Test Cases Table
    test_cases = {
        "SMILES String": [
            "C",
            "CCO",
            "CC(=O)OC1=CC=CC=C1C(=O)O",
            "C1=CC=CC=C1",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
        ],
        "Molecule": [
            "Methane",
            "Ethanol",
            "Aspirin",
            "Benzene",
            "Ibuprofen"
        ]
    }
    
    # Convert to pandas DataFrame and display
    df = pd.DataFrame(test_cases)
    
    st.markdown("### 🧪 Test Cases")
    st.table(df)
