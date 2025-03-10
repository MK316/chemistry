import streamlit as st
from transformers import pipeline

# Set page title and icon
st.set_page_config(page_title="Free AI Chatbot", page_icon="🤖")

# Title and subtitle
st.title("🤖 AI Chatbot")
st.write("Ask me anything!")

# Model selection with radio buttons
model_option = st.radio(
    "Select an AI model:",
    options=[
        "GPT-2 (Hugging Face)", 
        "Llama 2 (Meta, 7B)",
        "GPT-Neo (EleutherAI, 1.3B)",
        "GPT-J (EleutherAI, 6B)",
        "Bloom (BigScience, 1.7B)"
    ],
    index=0  # Default to the first option
)

# Load the selected model
if model_option == "GPT-2 (Hugging Face)":
    st.write("**GPT-2 model selected. Generating text...**")
    generator = pipeline("text-generation", model="gpt2")

elif model_option == "Llama 2 (Meta, 7B)":
    st.write("**Llama 2 model selected. Generating text...**")
    generator = pipeline("text-generation", model="meta-llama/Llama-2-7b-hf")

elif model_option == "GPT-Neo (EleutherAI, 1.3B)":
    st.write("**GPT-Neo model selected. Generating text...**")
    generator = pipeline("text-generation", model="EleutherAI/gpt-neo-1.3B")

elif model_option == "GPT-J (EleutherAI, 6B)":
    st.write("**GPT-J model selected. Generating text...**")
    generator = pipeline("text-generation", model="EleutherAI/gpt-j-6B")

elif model_option == "Bloom (BigScience, 1.7B)":
    st.write("**Bloom model selected. Generating text...**")
    generator = pipeline("text-generation", model="bigscience/bloom-1b7")

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

    # Generate AI response based on selected model
    with st.spinner("Thinking..."):
        try:
            response = generator(user_input, max_length=100, num_return_sequences=1)[0]["generated_text"]
            ai_message = response[len(user_input):].strip()
        except Exception as e:
            ai_message = f"An error occurred: {str(e)}"

    # Add AI response to session state
    st.session_state.messages.append({"role": "assistant", "content": ai_message})
    with st.chat_message("assistant"):
        st.markdown(ai_message)

# Button to clear the chat history
if st.button("Clear Chat"):
    st.session_state.messages = []
    st.rerun()
