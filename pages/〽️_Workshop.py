import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import io
import math
import random
import os
from PIL import Image
import soundfile as sf

tab1, tab2, tab3, tab4, tab5 = st.tabs(["📖 Lecture slides", "Web links", "🌀App2: Complex waves", "💦 Quiz", "🌀 Videos"])

# CSS to adjust the alignment of the dropdown to match the buttons
st.markdown("""
    <style>
    .stSelectbox div[data-baseweb="select"] {
        margin-top: -30px;  /* Adjust this value to align with the buttons */
    }
    </style>
    """, unsafe_allow_html=True)

# Set up the path to the slides folder
slides_path = "pages/slides/"  # Ensure this is correct relative to your app's location
slide_files = sorted([f for f in os.listdir(slides_path) if f.endswith(".png")])
num_slides = len(slide_files)

# Initialize session state variables if they do not exist

# Initialize session state variables if they do not exist
if "slide_index" not in st.session_state:
    st.session_state.slide_index = 0  # Start with the first slide

# Check if there are slides in the folder
if num_slides == 0:
    st.error("No slides found in the specified folder.")
    st.stop()  # Stop the app if there are no slides


# Function to load and display the image based on the current index with resizing
def display_image():
    slide_path = os.path.join(slides_path, slide_files[st.session_state.slide_index])
    image = Image.open(slide_path)
    
    # Set your desired width for resizing
    desired_width = 1200  # Adjust this value as needed
    aspect_ratio = image.height / image.width
    new_height = int(desired_width * aspect_ratio)
    resized_image = image.resize((desired_width, new_height))

    st.image(resized_image, caption=f"Slide {st.session_state.slide_index + 1} of {num_slides}")


#############################
with tab1:

    # Arrange 'Start', 'Previous', 'Next', and 'Slide Selector' in a single row
    col1, col2, col3, col4 = st.columns([1, 1, 1, 5])
    
    with col1:
        if st.button("⛳", key="start", help="Reset to the first slide"):
            st.session_state.slide_index = 0

    with col2:
        if st.button("◀️", key="previous", help="Go back to the previous slide"):
            if st.session_state.slide_index > 0:
                st.session_state.slide_index -= 1
            else:
                st.warning("This is the first slide.")

    with col3:
        if st.button("▶️", key="next", help="Go to the next slide"):
            if st.session_state.slide_index < num_slides - 1:
                st.session_state.slide_index += 1
            else:
                st.warning("Final slide")

    with col4:
        # Display slide selector dropdown
        selected_slide = st.selectbox("",
                                      options=[f"Slide {i + 1}" for i in range(num_slides)],
                                      index=st.session_state.slide_index)

        # Update slide index if dropdown selection changes
        selected_slide_index = int(selected_slide.split()[-1]) - 1
        if selected_slide_index != st.session_state.slide_index:
            st.session_state.slide_index = selected_slide_index

    # Display the image
    display_image()

    st.markdown("---")
    st.caption("Introduction: Slides 1~8")
    st.caption("Chapter 01: Slides 9~52")

########################
with tab2:
    st.markdown("### 💻 Useful Platforms for Coding and Application Design")

    st.write("""
    Here are some key platforms useful for education and coding-based application designs:
    """)

    platforms = {
        "Python": {
            "link": "https://www.python.org/",
            "description": "Python is a versatile programming language widely used in education and app development due to its readability and vast library support. (Python은 읽기 쉽고 광범위한 라이브러리를 지원하여 교육 및 앱 개발에 널리 사용되는 다용도 프로그래밍 언어입니다.)"
        },
        "Colab": {
            "link": "https://colab.research.google.com/",
            "description": "Google Colab allows users to write and execute Python code in a browser, supporting collaboration and GPU access for deep learning tasks. (Google Colab은 브라우저에서 Python 코드를 작성하고 실행할 수 있도록 하며, 딥러닝 작업을 위한 협업 및 GPU 접근을 지원합니다.)"
        },
        "GitHub": {
            "link": "https://github.com/",
            "description": "GitHub is a code hosting platform for version control and collaboration, allowing teams to work together on projects. (GitHub은 버전 관리 및 협업을 위한 코드 호스팅 플랫폼으로, 팀이 프로젝트에서 협력할 수 있도록 합니다.)"
        },
        "Streamlit": {
            "link": "https://streamlit.io/",
            "description": "Streamlit is an open-source Python framework for building and sharing data-based web applications quickly. (Streamlit은 데이터 기반 웹 애플리케이션을 빠르게 구축하고 공유할 수 있는 오픈소스 Python 프레임워크입니다.)"
        },
        "Gradio": {
            "link": "https://gradio.app/",
            "description": "Gradio allows developers to create customizable user interfaces for machine learning models and applications. (Gradio는 개발자가 머신러닝 모델 및 애플리케이션을 위한 사용자 인터페이스를 쉽게 만들 수 있도록 합니다.)"
        },
        "Huggingface": {
            "link": "https://huggingface.co/",
            "description": "Hugging Face is a leader in natural language processing (NLP) and machine learning model hosting and sharing. (Hugging Face는 자연어 처리(NLP) 및 머신러닝 모델 호스팅 및 공유 분야의 선두 주자입니다.)"
        },
        "ChatGPT": {
            "link": "https://chat.openai.com/",
            "description": "ChatGPT is a conversational AI model from OpenAI, widely used for language-based learning and automation. (ChatGPT는 OpenAI의 대화형 AI 모델로, 언어 기반 학습 및 자동화에 널리 사용됩니다.)"
        },
        "DeepSeek": {
            "link": "https://www.deepseek.com/",
            "description": "DeepSeek is a platform for deep learning research and model deployment, useful for AI-based applications. (DeepSeek는 딥러닝 연구 및 모델 배포를 위한 플랫폼으로, AI 기반 애플리케이션에 유용합니다.)"
        }
    }

    for platform, info in platforms.items():
        st.markdown(f"🌱 **[{platform}]({info['link']})** — {info['description']}")

    

########################
with tab3:
    st.markdown("### 🎶 Create a Complex Wave")

    # User selects number of sine wave components
    num_components = st.number_input("Number of sine waves to combine", min_value=1, max_value=5, value=2)

    frequencies = []
    amplitudes = []

    # User inputs for multiple frequencies and amplitudes
    for i in range(num_components):
        col1, col2 = st.columns(2)
        with col1:
            freq = st.number_input(f"Frequency {i+1} (Hz)", min_value=1, max_value=5000, value=440 if i == 0 else 880, step=1)
            frequencies.append(freq)
        with col2:
            amp = st.slider(f"Amplitude {i+1}", min_value=0.1, max_value=1.0, value=0.5, step=0.05)
            amplitudes.append(amp)

    duration = 1.0  # Fixed duration of 1 second
    sampling_rate = 44100  # Standard audio sampling rate

    # Generate time values
    t = np.linspace(0, duration, int(sampling_rate * duration), endpoint=False)

    # Initialize complex wave
    complex_wave = np.zeros_like(t)

    # Define a color palette for different waves
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]  # Blue, Orange, Green, Red, Purple

    # Create a figure with subplots for individual waves and final wave
    fig, axes = plt.subplots(num_components + 1, 1, figsize=(6, 2 * (num_components + 1)))

    for i in range(num_components):
        # Generate each sine wave
        sine_wave = amplitudes[i] * np.sin(2 * np.pi * frequencies[i] * t)
        complex_wave += sine_wave  # Add to the complex wave
        
        # Plot individual sine wave with different colors
        axes[i].plot(t[:1000], sine_wave[:1000], color=colors[i % len(colors)], linewidth=1.5)
        axes[i].set_title(f"Sine Wave {i+1}: {frequencies[i]} Hz, Amplitude: {amplitudes[i]}", fontsize=10)
        axes[i].set_xlabel("Time (s)", fontsize=9)
        axes[i].set_ylabel("Amplitude", fontsize=9)
        axes[i].tick_params(axis='both', labelsize=8)
        axes[i].grid(True)

    # Normalize the complex wave to prevent clipping
    complex_wave /= np.max(np.abs(complex_wave))

    # Plot the final complex wave
    axes[-1].plot(t[:1000], complex_wave[:1000], color="black", linewidth=1.5)  # Final wave in black
    axes[-1].set_title("Final Complex Waveform (Sum of All Components)", fontsize=10)
    axes[-1].set_xlabel("Time (s)", fontsize=8)
    axes[-1].set_ylabel("Amplitude", fontsize=8)
    axes[-1].tick_params(axis='both', labelsize=8)
    axes[-1].grid(True)

    # Adjust spacing to prevent overlapping
    plt.subplots_adjust(hspace=0.8)

    # Display all plots
    st.pyplot(fig)

    # Save the wave as a temporary audio file
    audio_buffer = io.BytesIO()
    sf.write(audio_buffer, complex_wave, sampling_rate, format='WAV')
    audio_buffer.seek(0)

    # Provide a download button for the generated sound
    st.audio(audio_buffer, format='audio/wav')
    st.download_button(label="Download Complex Wave File", data=audio_buffer, file_name="complex_wave.wav", mime="audio/wav")

#################################################


# Sample chemistry questions defined outside of any session-state-dependent function
questions = [
    {"question": "What is the chemical formula for water?",
     "options": ["H2O", "H2O2", "CO2", "O2"],
     "answer": "H2O"},
    {"question": "What is the pH of pure water at 25°C?",
     "options": ["3", "5", "7", "8"],
     "answer": "7"},
    {"question": "Which element has the chemical symbol 'O'?",
     "options": ["Gold", "Oxygen", "Silver", "Iron"],
     "answer": "Oxygen"}
]

def quiz_app():
    st.title('💦 Chemistry Basics Quiz')
    
    # Dropdown to select the question
    question_titles = [q["question"] for q in questions]  # Extracting questions for the dropdown
    selected_question = st.selectbox("Select a question:", question_titles)

    # Find the selected question in the list
    question = next((q for q in questions if q["question"] == selected_question), None)
    
    if question:
        # Display the options as radio buttons
        user_answer = st.radio("Choose an answer:", question["options"], key="answer")

        # Button to check the answer
        if st.button("Check answer"):
            if user_answer == question["answer"]:
                st.success("Correct! 🎉")
            else:
                st.error(f"Wrong! The correct answer is {question['answer']}.")
with tab4:
    quiz_app()
    

#########################

with tab5:

    video_url = "https://www.youtube.com/embed/XLfQpv2ZRPU?si=5zKkYufSdvLbsCp3"

    st.markdown("#### 1. Video: Understanding physical aspect of sound")
    st.markdown(
        f"""
        <iframe width="560" height="315" src="{video_url}" 
        title="YouTube video player" frameborder="0" 
        allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" 
        referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
        """,
        unsafe_allow_html=True
    )

    

    video_url2 = "https://www.youtube.com/embed/rYrdiQckGhw?si=nod5AuttYchOzH5X"

    st.markdown("#### 2. Fun experiments with sound")
    st.markdown(
        f"""
        <iframe width="560" height="315" src="{video_url2}" 
        title="YouTube video player" frameborder="0" 
        allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" 
        referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
        """,
        unsafe_allow_html=True
    )
  

  


