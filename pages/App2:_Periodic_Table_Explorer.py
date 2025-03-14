import streamlit as st
import pandas as pd
import plotly.express as px

tabs = st.tabs(["🔎 App overview", "🐾 Practice with APP"])

with tabs[0]:
    st.markdown("""
    #### Purpose:
    + Help students understand the periodic table structure and trends.
    + Allow interactive exploration of element properties.
    #### Key Features:
    + Display an interactive periodic table.
    + Allow students to click on an element to see detailed properties (atomic number, mass, electron configuration, etc.).
    + Highlight trends (e.g., atomic size, electronegativity) using color gradients.


    """)

with tabs[1]:
    st.markdown("#### 🌀 Periodic Table")
    # Load periodic table data
    df = pd.read_csv('https://raw.githubusercontent.com/Bowserinator/Periodic-Table-JSON/master/PeriodicTableCSV.csv')
    
    # Create a scatter plot for the periodic table
    fig = px.scatter(df, x="xpos", y="ypos", text="symbol", color="category",
                     hover_name="name", size_max=40, width=900, height=400)
    
    fig.update_traces(marker=dict(size=30), textposition='middle center')
    
    st.plotly_chart(fig)
    
    # Display details when an element is selected
    selected_element = st.selectbox("Choose an element", df['name'])
    if selected_element:
        element_data = df[df['name'] == selected_element].iloc[0]
        st.write(element_data)
    
    st.info("Data source from https://github.com/Bowserinator/")
