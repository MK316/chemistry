import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

def calculate_reaction_rate(k, initial_concentration, time):
    """ Simple first-order reaction rate calculation """
    return initial_concentration * np.exp(-k * time)

def main():
    st.title('Interactive Reaction Rate Simulator')

    with st.tabs(["Input", "Output"]):
        with st.tab("Input"):
            # User inputs for the reaction
            st.header("Enter Reaction Parameters")
            initial_concentration = st.number_input('Initial Concentration (mol/L)', min_value=0.0, value=1.0, step=0.1)
            rate_constant = st.number_input('Rate Constant (1/s)', min_value=0.0, value=0.1, step=0.01)
            max_time = st.slider('Time Range (s)', 0, 100, 50)

        with st.tab("Output"):
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
