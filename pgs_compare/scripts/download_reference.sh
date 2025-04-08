#!/bin/bash

# Define directory
DIR="${1:-data/reference}"

# Create directory
mkdir -p "$DIR"

# Download PGS compiled reference build
echo "Downloading PGS compiled reference build (HGDP+1kGP)..."
wget -O "$DIR/pgsc_HGDP+1kGP_v1.tar.zst" "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/resources/pgsc_HGDP+1kGP_v1.tar.zst"

if [ $? -eq 0 ]; then
    echo "Successfully downloaded HGDP+1kGP reference build"
else
    echo "Failed to download HGDP+1kGP reference build"
fi

# Download 1000 Genomes reference build
echo "Downloading 1000 Genomes reference build..."
wget -O "$DIR/pgsc_1000G_v1.tar.zst" "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/resources/pgsc_1000G_v1.tar.zst"

if [ $? -eq 0 ]; then
    echo "Successfully downloaded 1000 Genomes reference build"
else
    echo "Failed to download 1000 Genomes reference build"
fi

echo "Download process complete." 