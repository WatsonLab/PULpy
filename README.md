# PULpy
Open prediction of Polysaccharide Utilisation Loci (PUL)

# Create conda env
```sh
conda env create -f envs/PULpy.yaml
source activate PULpy
```

# Get Pfam data
```sh
# Pfam

mkdir pfam_data && cd pfam_data
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz
gunzip Pfam-A.hmm.gz Pfam-A.hmm.dat.gz active_site.dat.gz
hmmpress Pfam-A.hmm
cd ..
```

# Get DBCAN data
```sh
mkdir dbcan_data && cd dbcan_data
wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/hmmscan-parser.sh
wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/dbCAN-fam-HMMs.txt
hmmpress dbCAN-fam-HMMs.txt
chmod 755 hmmscan-parser.sh
cd ..
```

# Edit config.json if you need to....

# Make scripts executable

```sh
chmod -R 755 scripts
```

# Input file

PULpy is designed to process genomes downloaded from NCBI, therefore expects them in a format:

```
"genomes/{id}_genomic.fna.gz")
```

That is, in a directory called "genomes", with a unique_id followed by "_genomic.fna.gz"

# Run it
```sh
snakemake --use-conda
```


