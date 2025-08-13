#!/usr/bin/env bash

mkdir tcga_cna_gistic
cd tcga_cna_gistic

base="https://gdac.broadinstitute.org/runs/analyses__latest/data/"
folders=$(curl -s "$base" \
  | grep -oP '(?<=href=")[^"]+/' \
  | sed 's:/$::' \
  | grep '-' |
  | grep -v 'COADREAD' \
  | grep -v 'GBMLGG' \
  | grep -v 'KIPAN' \
  | grep -v 'STES')


for f in $folders; do
  url="${base}/${f}/20160128/gdac.broadinstitute.org_${f}.CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz"
  echo "Downloading ${f}..."
  wget -c "$url" || echo "Warning: failed to download $f"
done

for f in *.tar.gz; do
  tar -xvzf "$f" && rm "$f"
done

