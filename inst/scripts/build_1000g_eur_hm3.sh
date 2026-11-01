#!/usr/bin/env bash
# Build a 1000 Genomes phase 3 EUR (n = 503) PLINK reference restricted to the
# Pan-UKB HapMap3 SNP universe (UKBB.EUR.l2.ldscore.gz, 1,094,844 SNPs, GRCh37).
# All inputs come from public S3 buckets reachable from the cloud sandbox.
# Usage: build_1000g_eur_hm3.sh OUTDIR [CHROMS...]
set -euo pipefail
OUT=${1:?outdir}; shift || true
CHROMS=${@:-$(seq 1 22)}
mkdir -p "$OUT"; cd "$OUT"
KG=https://1000genomes.s3.amazonaws.com/release/20130502
[ -x plink ] || { curl -sf -o plink.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip && unzip -oq plink.zip plink; }
[ -s panel.txt ] || curl -sf -o panel.txt $KG/integrated_call_samples_v3.20130502.ALL.panel
awk '$3=="EUR"{print $1,$1}' panel.txt > eur.keep
[ -s UKBB.EUR.l2.ldscore.gz ] || curl -sf -o UKBB.EUR.l2.ldscore.gz https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.EUR.l2.ldscore.gz
zcat UKBB.EUR.l2.ldscore.gz | awk 'NR>1{print $1,$3,$3,$2}' > hm3.range
for c in $CHROMS; do
  [ -s eur_hm3_chr$c.bed ] && continue
  vcf=ALL.chr$c.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
  for i in 1 2 3 4; do curl -sf -o $vcf $KG/$vcf && break; sleep $((2**i)); done
  ./plink --vcf $vcf --keep eur.keep --extract range hm3.range --snps-only just-acgt \
    --biallelic-only strict --set-missing-var-ids @:#:\$1:\$2 --make-bed \
    --out eur_hm3_chr$c --memory 4000 > /dev/null
  rm -f $vcf
  echo "chr$c $(wc -l < eur_hm3_chr$c.bim) SNPs"
done
