# ==============================================================================
# Step 5b — Remove WARNING_REF_DOES_NOT_MATCH_GENOME variants
# ==============================================================================
echo ""
echo "[5b] Filtering WARNING_REF_DOES_NOT_MATCH_GENOME variants..."

SNPEFF_VCF_FILTERED="${OUTPUT_DIR}/all_annotated_snpeff.filtered.vcf"

bcftools filter \
    -e 'INFO/ANN ~ "WARNING_REF_DOES_NOT_MATCH_GENOME"' \
    -O v \
    -o "$SNPEFF_VCF_FILTERED" \
    "$SNPEFF_VCF"

if [ $? -ne 0 ]; then
    echo "[!] Filtering step failed." >&2
    exit 1
fi

# Quick sanity check
N_BEFORE=$(bcftools view -H "$SNPEFF_VCF" | wc -l)
N_AFTER=$(bcftools view -H "$SNPEFF_VCF_FILTERED" | wc -l)
echo "    Variants before filtering : $N_BEFORE"
echo "    Variants after filtering  : $N_AFTER"
echo "    Removed                   : $((N_BEFORE - N_AFTER))"
