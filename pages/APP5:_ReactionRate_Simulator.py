import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

def calculate_reaction_rate(k, initial_concentration, time):
    """ Simple first-order reaction rate calculation.
    Calculates the concentration of a reactant over time using the first-order kinetics formula:
    [A] = [A]₀ * exp(-kt)
    where [A]₀ is the initial concentration, k is the rate constant, and t is the time.
    """
    return initial_concentration * np.exp(-k * time)

def main():
    st.title('💫 Interactive Reaction Rate Simulator')
    st.write(
        """
        ## Overview
        This application simulates the reaction rate of a first-order chemical reaction. 
        Use the sidebar to input the initial concentration of the reactant, the rate constant, 
        and the time range over which you want to observe the reaction.
        
        The output is a plot that shows how the concentration of the reactant decreases over time.
        
        ### How to Use
        1. Enter the **Initial Concentration** of the reactant in moles per liter (mol/L).
        2. Input the **Rate Constant** of the reaction in per second (1/s).
        3. Adjust the **Time Range** slider to set the maximum time for the reaction observation.
        4. Click the **Calculate and Plot** button to generate the reaction rate plot.
        """
    )

    # Sidebar for inputs
    st.sidebar.header("Enter Reaction Parameters")
    initial_concentration = st.sidebar.number_input('Initial Concentration (mol/L)', min_value=0.0, value=1.0, step=0.1)
    rate_constant = st.sidebar.number_input('Rate Constant (1/s)', min_value=0.0, value=0.1, step=0.01)
    max_time = st.sidebar.slider('Time Range (s)', 0, 100, 50)

    # Main area for output display
    if st.sidebar.button('Calculate and Plot'):
        st.header("Reaction Rate Plot")
        time_values = np.linspace(0, max_time, 500)
        concentration_values = calculate_reaction_rate(rate_constant, initial_concentration, time_values)

        fig, ax = plt.subplots()
        ax.plot(time_values, concentration_values)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Concentration (mol/L)')
        ax.set_title('Concentration vs. Time')
        st.pyplot(fig)

if __name__ == "__main__":
    main()
