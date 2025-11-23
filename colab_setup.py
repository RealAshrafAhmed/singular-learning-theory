"""
Colab Setup Utility for GitHub Repository Integration

This module provides utilities to set up Google Colab environments for working
with private GitHub repositories. It handles SSH authentication, repository cloning,
dependency installation, and Google Drive mounting.

Usage:
    # In your Colab notebook, first mount Drive and add this script to path:
    from google.colab import drive
    drive.mount('/content/drive')
    import sys
    sys.path.insert(0, '/content/drive/MyDrive')  # Adjust path as needed
    
    import colab_setup
    paths = colab_setup.setup_repo(
        github_user="YourUsername",
        repo_name="your-repo",
        branch="main"  # optional
    )
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import Dict, Optional, Tuple


def is_colab() -> bool:
    """
    Check if the code is running in Google Colab.
    
    Returns:
        bool: True if running in Colab, False otherwise
    """
    try:
        import google.colab
        return True
    except ImportError:
        return False


def setup_ssh(ssh_key_name: str = 'ash@colab') -> bool:
    """
    Set up SSH authentication for GitHub using Colab secrets.
    
    This function:
    - Retrieves the SSH private key from Colab userdata/secrets
    - Fixes Windows line ending issues
    - Sets proper permissions
    - Adds GitHub to known_hosts
    - Tests the SSH connection
    
    Args:
        ssh_key_name: Name of the SSH key in Colab secrets (default: 'ash@colab')
    
    Returns:
        bool: True if SSH setup was successful, False otherwise
    """
    try:
        from google.colab import userdata
        
        print("🔑 Setting up SSH authentication...")
        
        # Get SSH key from secrets
        try:
            ssh_key = userdata.get(ssh_key_name)
        except Exception as e:
            print(f"❌ Failed to retrieve SSH key '{ssh_key_name}' from secrets")
            print(f"   Make sure you've added your SSH private key to Colab secrets")
            print(f"   Error: {e}")
            return False
        
        # Clean the key - handle various line ending formats
        if '\\n' in ssh_key:
            ssh_key = ssh_key.replace('\\n', '\n')
        
        # Remove carriage returns (Windows line endings)
        ssh_key = ssh_key.replace('\r\n', '\n')
        ssh_key = ssh_key.replace('\r', '')
        
        # Ensure proper format
        ssh_key = ssh_key.strip()
        if not ssh_key.endswith('\n'):
            ssh_key = ssh_key + '\n'
        
        # Create SSH directory with proper permissions
        ssh_dir = Path.home() / '.ssh'
        ssh_dir.mkdir(mode=0o700, exist_ok=True)
        
        # Write SSH private key
        key_file = ssh_dir / 'id_rsa'
        key_file.write_text(ssh_key)
        key_file.chmod(0o600)
        
        # Add GitHub to known hosts
        known_hosts = ssh_dir / 'known_hosts'
        result = subprocess.run(
            ['ssh-keyscan', 'github.com'],
            capture_output=True,
            text=True,
            check=False
        )
        if result.returncode == 0:
            known_hosts.write_text(result.stdout)
        
        # Test SSH connection
        test_result = subprocess.run(
            ['ssh', '-T', 'git@github.com'],
            capture_output=True,
            text=True,
            check=False
        )
        
        if 'successfully authenticated' in test_result.stderr:
            print("✅ SSH authentication successful!")
            return True
        else:
            print("⚠️  SSH test response:")
            print(test_result.stderr)
            return 'successfully authenticated' in test_result.stderr.lower()
            
    except Exception as e:
        print(f"❌ SSH setup failed: {e}")
        return False


def clone_repo(
    github_user: str,
    repo_name: str,
    branch: Optional[str] = None,
    target_dir: Optional[str] = None,
    use_ssh: bool = True
) -> Optional[Path]:
    """
    Clone a GitHub repository.
    
    Args:
        github_user: GitHub username or organization
        repo_name: Repository name
        branch: Branch to checkout (optional, defaults to repo's default branch)
        target_dir: Directory to clone into (optional, defaults to /content/{repo_name})
        use_ssh: Use SSH for cloning (default: True). Set to False for public repos
    
    Returns:
        Path: Path to the cloned repository, or None if cloning failed
    """
    try:
        # Set target directory
        if target_dir is None:
            target_dir = f'/content/{repo_name}'
        
        target_path = Path(target_dir)
        
        # Remove existing directory if it exists
        if target_path.exists():
            print(f"📁 Removing existing directory: {target_path}")
            subprocess.run(['rm', '-rf', str(target_path)], check=True)
        
        # Construct clone URL
        if use_ssh:
            repo_url = f'git@github.com:{github_user}/{repo_name}.git'
        else:
            repo_url = f'https://github.com/{github_user}/{repo_name}.git'
        
        print(f"📥 Cloning {github_user}/{repo_name}...")
        
        # Clone repository
        clone_cmd = ['git', 'clone', repo_url, str(target_path)]
        if branch:
            clone_cmd.extend(['--branch', branch])
        
        result = subprocess.run(
            clone_cmd,
            capture_output=True,
            text=True,
            check=False
        )
        
        if result.returncode != 0:
            print(f"❌ Clone failed:")
            print(result.stderr)
            return None
        
        print(f"✅ Repository cloned to: {target_path}")
        return target_path
        
    except Exception as e:
        print(f"❌ Repository cloning failed: {e}")
        return None


def install_requirements(repo_path: Path, requirements_file: str = 'requirements.txt') -> bool:
    """
    Install Python requirements from a requirements file.
    
    Args:
        repo_path: Path to the repository
        requirements_file: Name of the requirements file (default: 'requirements.txt')
    
    Returns:
        bool: True if installation was successful, False otherwise
    """
    try:
        req_file = repo_path / requirements_file
        
        if not req_file.exists():
            print(f"ℹ️  No {requirements_file} found, skipping dependency installation")
            return True
        
        print(f"📦 Installing dependencies from {requirements_file}...")
        
        result = subprocess.run(
            ['pip', 'install', '-q', '-r', str(req_file)],
            capture_output=True,
            text=True,
            check=False
        )
        
        if result.returncode != 0:
            print(f"⚠️  Some dependencies may have failed to install:")
            print(result.stderr)
            return False
        
        print("✅ Dependencies installed successfully!")
        return True
        
    except Exception as e:
        print(f"❌ Dependency installation failed: {e}")
        return False


def mount_drive(mount_point: str = '/content/drive') -> bool:
    """
    Mount Google Drive.
    
    Args:
        mount_point: Where to mount the drive (default: '/content/drive')
    
    Returns:
        bool: True if mounting was successful, False otherwise
    """
    try:
        from google.colab import drive
        
        if Path(mount_point).exists() and any(Path(mount_point).iterdir()):
            print(f"✅ Drive already mounted at {mount_point}")
            return True
        
        print(f"💾 Mounting Google Drive at {mount_point}...")
        drive.mount(mount_point)
        print("✅ Drive mounted successfully!")
        return True
        
    except Exception as e:
        print(f"❌ Drive mounting failed: {e}")
        return False


def setup_output_dir(output_dir: Optional[str] = None) -> Path:
    """
    Create and return the output directory path.
    
    Args:
        output_dir: Custom output directory (optional, defaults to Drive/MyDrive/outputs)
    
    Returns:
        Path: Path to the output directory
    """
    if output_dir is None:
        output_dir = '/content/drive/MyDrive/outputs'
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    print(f"📂 Output directory: {output_path}")
    
    return output_path


def setup_repo(
    github_user: str,
    repo_name: str,
    branch: Optional[str] = None,
    ssh_key_name: str = 'ash@colab',
    requirements_file: str = 'requirements.txt',
    output_dir: Optional[str] = None,
    mount_drive_flag: bool = True,
    use_ssh: bool = True
) -> Dict[str, Path]:
    """
    Complete setup for a GitHub repository in Google Colab.
    
    This is the main function that orchestrates the entire setup process:
    1. Checks if running in Colab
    2. Sets up SSH authentication (if use_ssh=True)
    3. Clones the repository
    4. Installs requirements
    5. Adds repo to Python path
    6. Mounts Google Drive
    7. Creates output directory
    
    Args:
        github_user: GitHub username or organization
        repo_name: Repository name
        branch: Branch to checkout (optional)
        ssh_key_name: Name of SSH key in Colab secrets (default: 'ash@colab')
        requirements_file: Name of requirements file (default: 'requirements.txt')
        output_dir: Custom output directory (optional)
        mount_drive_flag: Whether to mount Google Drive (default: True)
        use_ssh: Use SSH for cloning (default: True). Set to False for public repos
    
    Returns:
        Dict[str, Path]: Dictionary containing:
            - 'repo': Path to the cloned repository
            - 'output': Path to the output directory
            - 'drive': Path to the mounted drive (if mounted)
    
    Example:
        >>> paths = setup_repo(
        ...     github_user="RealAshrafAhmed",
        ...     repo_name="singular-learning-theory",
        ...     branch="main"
        ... )
        >>> print(f"Working in: {paths['repo']}")
        >>> print(f"Saving to: {paths['output']}")
    """
    print("=" * 60)
    print("🚀 Starting Colab Setup")
    print("=" * 60)
    
    # Check if in Colab
    if not is_colab():
        print("⚠️  Not running in Google Colab - skipping Colab-specific setup")
        return {}
    
    paths = {}
    
    # Change to /content directory
    os.chdir('/content')
    
    # Setup SSH if using SSH authentication
    if use_ssh:
        if not setup_ssh(ssh_key_name):
            print("⚠️  SSH setup failed - attempting to continue anyway...")
    
    # Clone repository
    repo_path = clone_repo(github_user, repo_name, branch, use_ssh=use_ssh)
    if repo_path is None:
        raise RuntimeError("Failed to clone repository")
    
    paths['repo'] = repo_path
    
    # Change to repo directory
    os.chdir(repo_path)
    
    # Install requirements
    install_requirements(repo_path, requirements_file)
    
    # Add repo to Python path
    if str(repo_path) not in sys.path:
        sys.path.insert(0, str(repo_path))
        print(f"✅ Added {repo_path} to Python path")
    
    # Mount Drive
    if mount_drive_flag:
        if mount_drive():
            paths['drive'] = Path('/content/drive')
    
    # Setup output directory
    output_path = setup_output_dir(output_dir)
    paths['output'] = output_path
    
    print("\n" + "=" * 60)
    print("✅ Setup Complete!")
    print("=" * 60)
    print(f"📁 Repository: {paths['repo']}")
    print(f"💾 Output Dir: {paths['output']}")
    if 'drive' in paths:
        print(f"🗂️  Drive: {paths['drive']}")
    print("=" * 60)
    
    return paths


# Quick setup function for common use case
def quick_setup(github_user: str, repo_name: str, branch: Optional[str] = None) -> Dict[str, Path]:
    """
    Simplified setup with common defaults.
    
    Args:
        github_user: GitHub username or organization
        repo_name: Repository name
        branch: Branch to checkout (optional)
    
    Returns:
        Dict[str, Path]: Dictionary with 'repo', 'output', and optionally 'drive' paths
    """
    return setup_repo(github_user, repo_name, branch)
