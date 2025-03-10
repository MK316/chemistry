import streamlit as st

# Set the page title and icon
st.set_page_config(page_title="Chemistry Learning Apps", page_icon="🧪")

# Title and introduction
st.markdown("### 🔬 Chemistry Learning App Gallery")
st.caption(
    """
    Welcome to the Chemistry Learning Apps! This platform is designed to help future chemistry teachers understand how to integrate technology into their teaching. This platform was created to demonstrate how Streamlit can be used to build educational apps for future chemistry teachers. (Since Mar. 10, 2025)  
    """)

# ✅ Add images using columns
col1, col2 = st.columns([1, 1])  # Divide into two equal columns

with col1:
    st.image("https://github.com/MK316/chemistry/raw/main/images/chemistry01.png", 
             caption="Enhance Chemistry Learning with Teacher-designed Applications", 
             width=300)  # Adjust width value as needed

with col2:
    st.image("https://github.com/MK316/chemistry/raw/main/images/chemi25.png", 
             caption="GNU Chemistry Edu Workshop 2025", 
             width=200)  # Adjust width value as needed


st.markdown("---")
# Overview of available tools
st.markdown("#### 🚀 Available Tools:")
st.markdown("1. **Molecular Structure Viewer** – Generate and visualize 2D structures of chemical compounds using SMILES.")
st.markdown("2. **Interactive Periodic Table** – Explore the periodic table and element properties interactively.")
st.markdown("3. **States of Matter Simulation** – To visualize how molecules behave in different states of matter.")
st.markdown("4. **pH and Acid-Base Indicator** – To understand acids, bases, and neutral solutions based on pH levels.  ")
st.markdown("5. **Chemical Equation Balancer** – Automatically balance chemical reactions using algebraic methods.")
# Footer

st.caption("©️Github: MK316")


