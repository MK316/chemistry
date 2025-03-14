import streamlit as st

# Set the page title and icon
st.set_page_config(page_title="Chemistry Learning Apps", page_icon="🧪")

# Title and introduction
st.markdown("### 🔬 Chemistry Learning App Gallery")
st.caption(
    """
    Welcame to the Chemistry Learning Apps! This platform is designed to help future chemistry teachers understand how to integrate technology into their teaching. This platform was created to demonstrate how Streamlit can be used to build educational apps for future chemistry teachers. (Since Mar. 10, 2025)  
    """)

# ✅ Add images using columns
col1, col2 = st.columns([1, 1])  # Divide into two equal columns

with col1:
    st.image("https://github.com/MK316/chemistry/raw/main/images/chemistry01.png", 
             caption="Enhance Chemistry Learning with Teacher-designed Applications", 
             width=300)  # Adjust width value as needed

with col2:
    st.image("https://github.com/MK316/chemistry/raw/main/images/chemi25.png", 
             caption="GNU ChemEdu Workshop 2025", 
             width=150)  # Adjust width value as needed


# Footer

st.caption("©️Github: MK316")


