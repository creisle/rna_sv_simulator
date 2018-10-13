#!/bin/bash
set -eu

# Download Flux and unpack it.
# Download ensembl assembly and annotation.

ENSEMBL_ROOT="ftp://ftp.ensembl.org/pub/release-69"
ENSEMBL_DNA="$ENSEMBL_ROOT"/"fasta/homo_sapiens/dna"
ENSEMBL_GTF="$ENSEMBL_ROOT"/"gtf/homo_sapiens/Homo_sapiens.GRCh37.69.gtf.gz"
FLUX_DOWNLOAD="http://artifactory.sammeth.net/artifactory/barna/barna/barna.simulator/1.2.1/flux-simulator-1.2.1.tgz"

FLUX_BUNDLE="flux-bundle"
GENOME_DIR="$FLUX_BUNDLE"/"genome"

# Create directories and whatnot
mkdir -p $GENOME_DIR

# Get flux and unpack
if [ ! -f $(basename $FLUX_DOWNLOAD) ]; then
    wget $FLUX_DOWNLOAD
    tar -xvf $(basename $FLUX_DOWNLOAD)    
    echo "Flux Simulator unpacked at $PWD"
else
    echo "Flux Simulator already downloaded."
fi

# Get chromosomes (.fasta)
pushd $GENOME_DIR/
for i in {1..22} X Y MT $(cat ../ensembl69_patches.txt | tr '\n' ' '); do
    CHROMOSOME="Homo_sapiens.GRCh37.69.dna.chromosome.""$i"".fa.gz"
    if [ ! -f $i".fa" ]; then
	wget "$ENSEMBL_DNA"/"$CHROMOSOME" \
	    && gunzip $CHROMOSOME \
	    && mv $(basename $CHROMOSOME | sed 's|\(.*\)\..*|\1|g') $i".fa" \
	    && echo "Chromosome $i downloaded." \
	    &       

    else
	echo "Chromosome $i already downloaded."	
    fi

done

## Get upplementary files/nonchromosomes.
NONCHROMOSOMAL_GZ="Homo_sapiens.GRCh37.69.dna.nonchromosomal.fa.gz"
NONCHROMOSOMAL_FILE=$(basename $NONCHROMOSOMAL_GZ | sed 's|\(.*\)\..*|\1|g')
wget "$ENSEMBL_DNA"/"$NONCHROMOSOMAL_GZ"
gunzip $NONCHROMOSOMAL_GZ
python ../../split_fasta_by_contig.py $NONCHROMOSOMAL_FILE

popd

# Get annotation (.gtf)
pushd $FLUX_BUNDLE
if [ ! -f $(basename $ENSEMBL_GTF | sed 's|\(.*\)\..*|\1|g' ) ]; then
    # Remove symbolic link if stale.
    rm -f annotation.gtf

    wget $ENSEMBL_GTF
    gunzip $(basename $ENSEMBL_GTF)
    ln -s $(basename $ENSEMBL_GTF | sed 's|\(.*\)\..*|\1|g') annotation.gtf
    echo "GTF file downloaded at $PWD"
else
    echo "GTF annotation already downloaded."    
fi
popd
