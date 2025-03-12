import streamlit as st
import py3Dmol

def show_molecule(mol_id):
    xyzview = py3Dmol.view(query='cid:{}'.format(mol_id))
    xyzview.setStyle({'stick': {}})
    xyzview.zoomTo()
    return xyzview.show()

def main():
    st.title('3D Molecular Orbital Viewer')
    st.write("This app demonstrates the 3D structure of molecules. Enter a PubChem compound ID to view its molecular structure.")
    
    # User input for the molecule's PubChem CID
    mol_id = st.text_input("Enter PubChem CID:", "1234")  # Example CID for butane

    # Button to show the molecule
    if st.button('Show Molecule'):
        st.write("Rendering 3D model for CID: {}".format(mol_id))
        show_molecule(mol_id)

if __name__ == "__main__":
    main()
