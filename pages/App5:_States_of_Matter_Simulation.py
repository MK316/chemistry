import streamlit as st
import matplotlib.pyplot as plt
import numpy as np


tabs = st.tabs(["Overview", "Practice with APP"])

with tabs[0]:
    st.markdown("""

    #### Purpose:
    + Help students understand the three main states of matter: solid, liquid, and gas.
    + Demonstrate how molecules behave in different states through animation.
    #### Key Features:
    + Allow user to select a state of matter.
    + Display dynamic particle animation (movement, spacing) to reflect molecular behavior.
    + Explain key properties of each state (e.g., density, shape, volume).


    """)

with tabs[1]:

    # Define function to plot molecular behavior
    def plot_particles(state):
        plt.figure(figsize=(5, 5))
        
        if state == "Solid":
            x = np.random.uniform(0.1, 0.9, 50)
            y = np.random.uniform(0.1, 0.9, 50)
            plt.scatter(x, y, s=200, color="blue")  # Densely packed
            plt.title("Solid - Fixed shape and volume")
            
        elif state == "Liquid":
            x = np.random.uniform(0.1, 0.9, 30)
            y = np.random.uniform(0.1, 0.9, 30)
            plt.scatter(x, y, s=150, color="green")  # Less dense
            plt.title("Liquid - No fixed shape, fixed volume")
            
        elif state == "Gas":
            x = np.random.uniform(0.1, 0.9, 20)
            y = np.random.uniform(0.1, 0.9, 20)
            plt.scatter(x, y, s=100, color="red")  # Freely moving
            plt.title("Gas - No fixed shape or volume")
            
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.axis("off")
        st.pyplot(plt)
    
    st.title("🧪 States of Matter Simulation")
    
    # State selection
    state = st.radio("Select a state of matter:", ["Solid", "Liquid", "Gas"])

# Plot based on state
plot_particles(state)

# Explanation
if state == "Solid":
    st.markdown("**Solid:** Molecules are closely packed together, and the shape and volume are fixed.")
elif state == "Liquid":
    st.markdown("**Liquid:** Molecules are less tightly packed, the volume is fixed, but the shape can change.")
else:
    st.markdown("**Gas:** Molecules move freely, and both the shape and volume can change.")
