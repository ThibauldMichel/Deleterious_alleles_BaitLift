#!/bin/bash
#SBATCH --job-name="test GTF"
#SBATCH --partition=medium
#SBATCH --export=ALL

source /mnt/apps/users/tmichel/conda/etc/profile.d/conda.sh
conda activate perl_env

BMAS_GFF="/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.gff"
BMAS_FASTA="/home/tmichel/projects/rbge/tmichel/reference_genomes/Bmas.fa"
OUTPUT_DIR="${1:-.}"
ORG="Bmas"

WORK_DIR="$OUTPUT_DIR/sift4g_build_work"
mkdir -p "$WORK_DIR"
GTF_FILE="$WORK_DIR/${ORG}.gene.gtf"

echo ""
echo "    [3b-a] Converting GFF to GTF format..."

agat_convert_sp_gff2gtf.pl --gff "$BMAS_GFF" -o "$GTF_FILE"
if [ $? -ne 0 ] || [ ! -f "$GTF_FILE" ]; then
    echo "[!] AGAT GFF->GTF conversion failed." >&2
    exit 1
fi

# Fix biotype attributes — ensembl_gene_format_to_ucsc.pl requires
# gene_biotype "protein_coding" on every feature line
echo "    Fixing biotype attributes..."
sed -i 's/original_biotype "mrna"/gene_biotype "protein_coding"/g' "$GTF_FILE"
sed -i '/\texon\t/s/$/ gene_biotype "protein_coding";/' "$GTF_FILE"
sed -i '/\tCDS\t/s/$/ gene_biotype "protein_coding";/' "$GTF_FILE"


echo "    Verifying line formats:"
echo "    exon:"
grep "exon" "$GTF_FILE" | head -2 | sed 's/^/      /'
echo "    CDS:"
grep $'\tCDS\t' "$GTF_FILE" | head -2 | sed 's/^/      /'
echo "    transcript:"
grep "transcript" "$GTF_FILE" | grep -v "^#" | head -2 | sed 's/^/      /'

echo "    Feature counts:"
cut -f3 "$GTF_FILE" | sort | uniq -c | sed 's/^/      /'

echo "    Compressing FASTA and GTF..."
FASTA_GZ="$WORK_DIR/${ORG}.fa.gz"
GTF_GZ="$WORK_DIR/${ORG}.gene.gtf.gz"

gzip -c "$BMAS_FASTA" > "$FASTA_GZ"
gzip -c "$GTF_FILE"   > "$GTF_GZ"

if [ ! -f "$FASTA_GZ" ] || [ ! -f "$GTF_GZ" ]; then
    echo "[!] Failed to compress FASTA or GTF." >&2
    exit 1
fi

echo "    Done."
