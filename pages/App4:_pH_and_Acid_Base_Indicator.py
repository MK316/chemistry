import streamlit as st

tabs = st.tabs(["Overview", "Practice with APP"])
with tabs[0]:

    st.markdown("""
    
    #### Purpose:  
    + Teach students the concept of pH and how acids and bases work.
    + Show how color changes reflect pH levels using common indicators (e.g., litmus paper).
    #### Key Features:
    + Allow the user to input a pH value.
    + Display whether the solution is acidic, neutral, or basic.
    + Simulate color changes based on the pH level using a dynamic color display.
    """)
with tabs[1]:
  
    # Function to determine pH level and color
    def get_ph_info(ph):
        if ph < 0 or ph > 14:
            return "Invalid pH value", "white"
        elif ph < 7:
            return "Acidic", "#ff4b4b"  # Red for acidic
        elif ph == 7:
            return "Neutral", "#90ee90"  # Green for neutral
        else:
            return "Basic", "#6495ed"  # Blue for basic
    
    st.title("🌡️ pH and Acid-Base Indicator")
    
    # pH value input
    ph = st.slider("Select a pH value:", min_value=0.0, max_value=14.0, value=7.0)
    
    # Get pH level and color
    ph_status, color = get_ph_info(ph)
    
    # Display result
    st.markdown(f"### The solution is **{ph_status}**.")
    st.markdown(
        f'<div style="width:100%; height:50px; background-color:{color};"></div>',
        unsafe_allow_html=True
    )
    
    # Add a note about pH ranges
    st.info("🔎 **pH Guide:**\n- Acidic: pH < 7\n- Neutral: pH = 7\n- Basic: pH > 7")
