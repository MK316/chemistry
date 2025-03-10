import streamlit as st

# Define a dictionary of fixed questions and answers
qa_pairs = {
    "What time is the class?": "The class starts at 9:00 AM on Monday and Wednesday.",
    "When is the assignment due?": "The next assignment is due on Friday at 11:59 PM.",
    "Is there a make-up class?": "Yes, the make-up class is scheduled for next Tuesday at 10:00 AM.",
    "What’s the exam schedule?": "The midterm exam is on March 25th, and the final exam is on May 15th.",
    "When are office hours?": "Office hours are from 2:00 PM to 4:00 PM on Thursdays."
}

# Streamlit app layout
st.title("🎓 Class Info Chatbot")

# Dropdown for selecting a question
selected_question = st.selectbox(
    "Select a question about the class:",
    options=["Select a question..."] + list(qa_pairs.keys())
)

# Display the answer if a valid question is selected
if selected_question != "Select a question...":
    response = qa_pairs.get(selected_question)
    if response:
        st.success(response)

# Optional: Button to reset the dropdown
if st.button("Reset"):
    st.rerun()

# Add a helpful footer or note
st.markdown("---")
st.caption("💡 Tip: If you don't see your question listed, try asking during class or office hours!")
