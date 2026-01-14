
import os
from google import genai
from google.genai import types
from dotenv import load_dotenv

load_dotenv()
api_key = os.getenv("GOOGLE_API_KEY")

def main():
    # 1. Load project context
    load_dotenv()
    api_key = os.getenv("GOOGLE_API_KEY")
    
    # 2. Configure Client (Passing proxy explicitly if environment is stubborn)
    # The manual suggests nurion-dm nodes handle proxying best
    client = genai.Client(
        api_key=api_key,
        http_options=types.HttpOptions(
            client_args={'proxy': 'http://proxy.ksc.re.kr:8080'}
        )
    )

    print(f"🚀 Running AI Analysis in: {os.getcwd()}")
    
    try:
        response = client.models.generate_content(
            model="gemini-2.0-flash", 
            contents="Explain the VASP convergence criteria for a Cu(111) surface."
        )
        print(f"\nAI Response:\n{response.text}")
    except Exception as e:
        print(f"❌ Connection failed. Ensure you are on 'nurion-dm' node. Error: {e}")

if __name__ == "__main__":
    main()
