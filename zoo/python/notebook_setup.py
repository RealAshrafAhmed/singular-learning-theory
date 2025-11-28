import os
import subprocess
import sys
from pathlib import Path

# --- CONFIGURATION ---
# Replace these values with your actual GitHub repository details
REPO_URL = "https://github.com/RealAshrafAhmed/singular-learning-theory"
REPO_NAME = "singular-learning-theory"  # The name of your local repository folder
REQUIREMENTS_FILE = "requirements.txt"
# ---------------------


def is_colab_environment():
    """Detects if the notebook is running in Google Colab."""
    return "google.colab" in sys.modules


def setup_environment():
    """
    Sets up the environment by cloning the repo and installing dependencies
    if running in Colab and the repository hasn't been cloned yet.
    """
    # 1. Check if running in Google Colab
    if not is_colab_environment():
        print(
            "Running in local Jupyter environment. Assuming data and dependencies are available."
        )
        # Optionally, you can add local-specific setup here if needed.
        return False  # Indicate that Colab setup was skipped

    print("Running in Google Colab environment.")

    # 2. Define expected path for the repository
    repo_path = Path(f"/content/{REPO_NAME}")

    # 3. Check if the repository has already been cloned
    if repo_path.is_dir():
        print(f"Repository '{REPO_NAME}' already cloned. Skipping clone.")
        # Change directory into the repository path
        os.chdir(f"{repo_path}/code/python")
    else:
        print(f"Cloning repository from {REPO_URL}...")
        try:
            # Clone the repository
            subprocess.run(
                ["git", "clone", REPO_URL], check=True, capture_output=True, text=True
            )
            print("Cloning successful.")
            # Change directory into the newly cloned repository path
            os.chdir(f"{repo_path}/code/python")
        except subprocess.CalledProcessError as e:
            print(f"Error cloning repository: {e.stderr}")
            print("Please ensure the REPO_URL is correct and the repository is public.")
            return False

    # 4. Check for and install dependencies
    if Path(REQUIREMENTS_FILE).is_file():
        print(f"Installing dependencies from {REQUIREMENTS_FILE}...")
        try:
            # Use -q (quiet) for cleaner output in Colab
            subprocess.run(
                [sys.executable, "-m", "pip", "install", "-r", REQUIREMENTS_FILE, "-q"],
                check=True,
            )
            print("Dependencies installed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error installing dependencies: {e.stderr}")
            return False
    else:
        print(
            f"Warning: {REQUIREMENTS_FILE} not found. Skipping dependency installation."
        )

    # 5. Add the repository directory to Python's path
    # This allows you to import local modules within your cloned repo.
    if str(repo_path) not in sys.path:
        sys.path.insert(0, str(repo_path))
        print(f"Added '{repo_path}' to sys.path.")

    return True  # Indicate that Colab setup was performed


if __name__ == "__main__":
    setup_environment()
    # The main notebook logic should follow the setup.
    # In the actual notebook, you will use %run to execute this.
