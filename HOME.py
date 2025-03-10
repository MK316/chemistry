import streamlit as st

# Set the page title and icon
st.set_page_config(page_title="Chemistry Learning Apps", page_icon="🧪")

# Title and introduction
st.title("🔬 Chemistry Learning Apps")
st.write(
    """
    Welcome to the Chemistry Learning Apps!  
    This platform is designed to help future chemistry teachers understand how to integrate technology into their teaching.  
    Below are the available tools:
    """
)

# Overview of available tools
st.markdown("### 🚀 Available Tools:")
st.markdown("1. **Molecular Structure Viewer** – Generate and visualize 3D structures of chemical compounds using SMILES.")
st.markdown("2. **Interactive Periodic Table** – Explore the periodic table and element properties interactively.")
st.markdown("3. **Chemical Equation Balancer** – Automatically balance chemical reactions using algebraic methods.")

# Footer
st.caption("This platform was created to demonstrate how Streamlit can be used to build educational apps for future chemistry teachers.")

# Add an image or logo (Optional)
st.image("https://github.com/MK316/chemistry/raw/main/images/chemistry01.png", 
         caption="Enhance Chemistry Learning with Technology", 
         width=400)  # Adjust width value as needed
