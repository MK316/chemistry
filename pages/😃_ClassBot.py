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

# User input
question = st.text_input("Ask a question about the class:")

# Check for matching answers
if question:
    response = qa_pairs.get(question)
    
    if response:
        st.success(response)
    else:
        st.warning("Sorry, I can only answer questions about class time, assignments, make-up classes, and exams.")

# Display possible questions to guide the user
st.markdown("### 📌 Example Questions:")
for q in qa_pairs.keys():
    st.markdown(f"- {q}")

