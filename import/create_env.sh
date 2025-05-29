#!/bin/bash

# Creates an 'omero' mircomamba environment with omero-py
# omero-cli-render and omero-cli-zarr ready for importing
# ome.zarr images

# Check if micromamba is already installed
if command -v micromamba &> /dev/null; then
    echo "micromamba is already installed"
else
    # Detect OS
    if [[ "$(uname)" == "Darwin" ]]; then
        OS="osx-64"
    elif [[ "$(uname)" == "Linux" ]]; then
        OS="linux-64"
    else
        echo "Unsupported operating system"
        exit 1
    fi

    # Install micromamba
    cd "$HOME"
    echo "Installing micromamba for $OS..."
    curl -Ls "https://micro.mamba.pm/api/micromamba/$OS/latest" | tar -xvj bin/micromamba

    # Set up micromamba
    mkdir -p "$HOME/.local/bin"
    mv bin/micromamba "$HOME/.local/bin/"

    # Initialize micromamba
    eval "$("$HOME/.local/bin/micromamba" shell hook -s bash)"

    # Add micromamba to PATH
    export PATH="$HOME/.local/bin:$PATH"

    # Configure micromamba
    micromamba shell init -s bash

    echo "Micromamba has been installed and configured."
    echo "Please restart your shell or run: source ~/.bashrc"
fi

# Create omero environment and install packages
echo "Creating omero environment..."
micromamba create -n omero -c conda-forge -c ome omero-py git -y

# Activate the environment
eval "$(micromamba shell hook -s bash)"
micromamba activate omero
echo "omero environment has been created and activated."

pip install -U omero-cli-render
echo "omero-cli-render has been installed."

# Clone and install omero-cli-zarr
cd "$HOME"
git clone https://github.com/dominikl/omero-cli-zarr.git
cd omero-cli-zarr
git checkout extinfo
pip install .
cd ..
rm -rf omero-cli-zarr
echo "omero-cli-zarr has been installed."

echo "Installation complete!"
echo "You can activate the environment in the future with: micromamba activate omero"
