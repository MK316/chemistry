import streamlit as st
from transformers import pipeline

# Load the Hugging Face pipeline (no API key required)
generator = pipeline("text-generation", model="gpt2")

# Set page title and icon
st.set_page_config(page_title="Free AI Chatbot", page_icon="🤖")

# Title and subtitle
st.title("🤖 AI Chatbot")
st.write("Ask me anything!")

# Initialize session state to store conversation history
if "messages" not in st.session_state:
    st.session_state.messages = []

# Display previous messages
for message in st.session_state.messages:
    with st.chat_message(message["role"]):
        st.markdown(message["content"])

# User input
if user_input := st.chat_input("Type your message..."):
    # Add user message to session state
    st.session_state.messages.append({"role": "user", "content": user_input})
    with st.chat_message("user"):
        st.markdown(user_input)

    # Generate AI response using GPT-2
    with st.spinner("Thinking..."):
        response = generator(user_input, max_length=100, num_return_sequences=1)[0]["generated_text"]

    # Add AI response to session state
    ai_message = response[len(user_input):].strip()
    st.session_state.messages.append({"role": "assistant", "content": ai_message})
    with st.chat_message("assistant"):
        st.markdown(ai_message)

# Button to clear the chat history
if st.button("Clear Chat"):
    st.session_state.messages = []
    st.rerun()
