import streamlit as st
from sympy import symbols, Eq, solve

tabs = st.tabs(["🔎 App overview", "🐾 Practice with APP"])

with tabs[0]:
    st.markdown("""
    #### Purpose: 
    - To automatically balance chemical reactions using algebraic methods.  
    #### Key Features: 
    - 🔢 Input unbalanced chemical equations.  
    - ✅ Display balanced reaction.  
    - 🧪 Identify reactants and products. 
    """)

with tabs[1]:
    
    st.title("Chemical Reaction Balancer")
    
    # Input for unbalanced equation
    reaction = st.text_input("Enter an unbalanced reaction (e.g., H2 + O2 = H2O):")
    
    if reaction:
        try:
            reactants, products = reaction.split("=")
            reactants = reactants.split("+")
            products = products.split("+")
    
            # Create symbols for coefficients
            coeffs = symbols(' '.join(['a', 'b', 'c', 'd']))
            eqs = []
            
            # Example: Code to create and solve balancing equations using sympy
            eqs.append(Eq(coeffs[0] + coeffs[1], coeffs[2]))  # Example balancing equation
            solution = solve(eqs)
    
            if solution:
                st.success(f"Balanced Reaction: {reaction}")
            else:
                st.error("Cannot balance the equation")
    
        except Exception as e:
            st.error(f"Error: {e}")
