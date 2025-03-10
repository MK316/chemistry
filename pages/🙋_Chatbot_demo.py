import streamlit as st
from transformers import pipeline

# Set page title and icon
st.set_page_config(page_title="Free AI Chatbot", page_icon="🤖")

# Title and subtitle
st.title("🤖 AI Chatbot")
st.write("Ask me nothing! (This is just a demo. I hope you don't expect much.)")

# Model selection with radio buttons
model_option = st.radio(
    "Select an AI model:",
    options=[
        "GPT-2 (Hugging Face)", 
    ],
    index=0  # Default to the first option
)

# Initialize the generator based on the selected model
if model_option == "GPT-2 (Hugging Face)":
    st.write("**GPT-2 model selected. Generating text...**")
    generator = pipeline("text-generation", model="gpt2")
elif model_option == "Talk to me (GPT-2)":
    st.write("**Talk to me (GPT-2) model selected. Generating text...**")
    generator = pipeline("text-generation", model="https://chatgpt.com/g/g-m8RLcYhOz-talktome")
else:
    st.error("Selected model is not available. Please choose a valid model.")
    generator = None  # This sets generator to None if no valid model is selected.

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

    # Check if generator is initialized before generating response
    if generator is not None:
        with st.spinner("Thinking..."):
            try:
                response = generator(user_input, max_length=100, num_return_sequences=1)[0]["generated_text"]
                ai_message = response[len(user_input):].strip()
            except Exception as e:
                ai_message = f"An error occurred: {str(e)}"
    else:
        ai_message = "No model selected or model initialization failed."

    # Add AI response to session state
    st.session_state.messages.append({"role": "assistant", "content": ai_message})
    with st.chat_message("assistant"):
        st.markdown(ai_message)

# Button to clear the chat history
if st.button("Clear Chat"):
    st.session_state.messages = []
    st.rerun()
