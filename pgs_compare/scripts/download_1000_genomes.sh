#!/bin/bash

# Define directory
DIR="${1:-data/1000_genomes}"

# Create directories
mkdir -p "$DIR"

# Base URL for 1000 Genomes Phase 3 data
BASE_URL="ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

# Download metadata
echo "Downloading metadata..."
wget -P "$DIR" "$BASE_URL/integrated_call_samples_v3.20130502.ALL.panel"

# Process .panel to create sex info file
echo "Generating sex info file..."
echo -e "#IID\tSEX" > "$DIR/sex_info.txt"
awk 'NR>1 {if($4=="male") print $1"\t1"; else if($4=="female") print $1"\t2"; else print $1"\t0"}' \
    "$DIR/integrated_call_samples_v3.20130502.ALL.panel" >> "$DIR/sex_info.txt"

# Function to download and process a chromosome
process_chromosome() {
    local chr=$1
    local file=$2
    
    echo "Downloading chromosome ${chr} to $DIR..."
    wget -P "$DIR" "${BASE_URL}/${file}"
    
    echo "Converting to plink for chr_${chr}..."
    if [ "$chr" = "X" ]; then
        plink2 --vcf "$DIR/${file}" \
               --update-sex "$DIR/sex_info.txt" \
               --split-par b37 \
               --make-pgen \
               --out "$DIR/chr_${chr}"
    else
        plink2 --vcf "$DIR/${file}" \
               --update-sex "$DIR/sex_info.txt" \
               --make-pgen \
               --out "$DIR/chr_${chr}"
    fi
    
    # Clean up .vcf.gz file
    rm "$DIR/${file}"
}

# Process autosomes 1 to 22
for CHR in {1..22}; do
    FILE="ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    process_chromosome "${CHR}" "${FILE}"
done

# Process non-autosomal chromosomes
process_chromosome "MT" "ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz"
process_chromosome "X" "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"
process_chromosome "Y" "ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz"

echo "Download complete." 