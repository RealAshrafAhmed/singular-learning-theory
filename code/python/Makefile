# .PHONY: declares targets that do not correspond to actual files.
# This ensures make runs the commands even if a file with the same name exists.
.PHONY: all venv install clean activate

# Define variables for the virtual environment directory and executables
VENV_DIR = .venv
PYTHON_EXECUTABLE = $(VENV_DIR)/bin/python3
PIP_EXECUTABLE = $(VENV_DIR)/bin/pip

# Default target: runs 'install-pymc' which in turn depends on 'venv'
all: install

# Target to create the virtual environment
# Checks if the .venv directory already exists to avoid recreating it.
venv:
	@echo "Creating virtual environment in $(VENV_DIR)..."
	@if [ ! -d "$(VENV_DIR)" ]; then \
		python3 -m venv $(VENV_DIR); \
		echo "Virtual environment created."; \
	else \
		echo "Virtual environment already exists."; \
	fi

# Target to install PyMC within the virtual environment
# It depends on 'venv' to ensure the environment is created first.
# Commands are run in a subshell, so we directly use the venv's pip executable.
install: venv
	@echo "Upgrading pip and installing dependencies from requirements.txt..."
	# Ensure pip is up-to-date within the virtual environment
	@$(PIP_EXECUTABLE) install --upgrade pip
	# Install packages listed in requirements.txt
	@$(PIP_EXECUTABLE) install -r requirements.txt
	@echo "Dependencies installed successfully."

dev: venv
	@echo "Upgrading pip and installing dev dependencies..."
	# Ensure pip is up-to-date within the virtual environment
	@$(PIP_EXECUTABLE) install --upgrade pip
	# Install packages listed in requirements-dev.txt
	@$(PIP_EXECUTABLE) install -r requirements-dev.txt
	@echo "Development dependencies installed successfully."

# Target to provide instructions for manual activation
# Makefiles run commands in subshells, so 'source' commands won't persist
# in your main terminal session. This target tells you how to activate manually.
activate:
	@echo "To activate the virtual environment, run the following command in your terminal:"
	@echo "source $(VENV_DIR)/bin/activate"

# Target to clean up the virtual environment
clean:
	@echo "Removing virtual environment directory $(VENV_DIR)..."
	@rm -rf $(VENV_DIR)
	@echo "Virtual environment removed."

format:
	@$(PYTHON_EXECUTABLE) -m ruff format normal-mixtures-1d

lab: dev install
	@$(PYTHON_EXECUTABLE) -m ipykernel install --user --name=$(VENV_DIR) --display-name="$(VENV_DIR)"
	@$(PYTHON_EXECUTABLE) -m jupyter lab

mixnorm_ksigma_rlctx: install
	@$(PYTHON_EXECUTABLE) -m estimation.main mixnorm_ksigma rlctx 5000 .\
		--pymc-draws 10000\
		--pymc-chains 4\
		--pymc-cores 4\
		--pymc-tune 4000
#	    --pymc-progressbar\
