import streamlit as st
import plotly.express as px
import pandas as pd
import numpy as np
import time

# Create tabs
tabs = st.tabs(["🔎 App overview", "🐾 Practice with APP"])

# Overview Tab
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

# Practice Tab
with tabs[1]:
    st.title("🧪 States of Matter Simulation")

    # State selection
    state = st.radio("Select a state of matter:", ["Solid", "Liquid", "Gas"])

    # Create empty container for the simulation
    plot_area = st.empty()

    # Function to create particle simulation
    def simulate_particles(num_particles, size, speed):
        # Initialize particle positions
        x = np.random.uniform(0, 1, num_particles)
        y = np.random.uniform(0, 1, num_particles)

        for _ in range(100):  # Number of frames for animation
            # Update positions randomly to simulate motion
            x += np.random.uniform(-speed, speed, num_particles)
            y += np.random.uniform(-speed, speed, num_particles)
            
            # Keep particles within bounds
            x = np.clip(x, 0, 1)
            y = np.clip(y, 0, 1)

            # Create DataFrame for Plotly
            df = pd.DataFrame({'x': x, 'y': y})
            
            # Create scatter plot with Plotly
            fig = px.scatter(df, x="x", y="y", size_max=size, size=[size]*num_particles)
            fig.update_traces(marker=dict(color="blue" if state == "Solid" else 
                                          "green" if state == "Liquid" else "red"))
            fig.update_layout(width=400, height=400, 
                              xaxis=dict(visible=False), 
                              yaxis=dict(visible=False))

            # Display plot in real time
            plot_area.plotly_chart(fig, use_container_width=True)

            # Short delay to control animation speed
            time.sleep(0.1)

    # Simulate based on state
    if state == "Solid":
        simulate_particles(num_particles=50, size=15, speed=0.005)
        st.markdown("**Solid:** Molecules are closely packed together, and the shape and volume are fixed.")
    elif state == "Liquid":
        simulate_particles(num_particles=30, size=12, speed=0.02)
        st.markdown("**Liquid:** Molecules are less tightly packed, the volume is fixed, but the shape can change.")
    else:
        simulate_particles(num_particles=20, size=10, speed=0.05)
        st.markdown("**Gas:** Molecules move freely, and both the shape and volume can change.")

