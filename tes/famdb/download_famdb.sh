#!/usr/bin/env bash

BASE_URL="https://www.dfam.org/releases/current/families/FamDB"
OUTDIR="${1:-.}"

mkdir -p "$OUTDIR"
cd "$OUTDIR" || exit 1

FILES=(
  README.txt
  dfam39_full.0.h5.gz   dfam39_full.0.h5.gz.md5
  dfam39_full.1.h5.gz   dfam39_full.1.h5.gz.md5
  dfam39_full.2.h5.gz   dfam39_full.2.h5.gz.md5
  dfam39_full.3.h5.gz   dfam39_full.3.h5.gz.md5
  dfam39_full.4.h5.gz   dfam39_full.4.h5.gz.md5
  dfam39_full.5.h5.gz   dfam39_full.5.h5.gz.md5
  dfam39_full.6.h5.gz   dfam39_full.6.h5.gz.md5
  dfam39_full.7.h5.gz   dfam39_full.7.h5.gz.md5
  dfam39_full.8.h5.gz   dfam39_full.8.h5.gz.md5
  dfam39_full.9.h5.gz   dfam39_full.9.h5.gz.md5
  dfam39_full.10.h5.gz  dfam39_full.10.h5.gz.md5
  dfam39_full.11.h5.gz  dfam39_full.11.h5.gz.md5
  dfam39_full.12.h5.gz  dfam39_full.12.h5.gz.md5
  dfam39_full.13.h5.gz  dfam39_full.13.h5.gz.md5
  dfam39_full.14.h5.gz  dfam39_full.14.h5.gz.md5
  dfam39_full.15.h5.gz  dfam39_full.15.h5.gz.md5
  dfam39_full.16.h5.gz  dfam39_full.16.h5.gz.md5
)

for FILE in "${FILES[@]}"; do
  wget -c --progress=bar:force -O "$FILE" "$BASE_URL/$FILE"
  if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to download $FILE" >&2
    exit 1
  fi
done

echo ""
echo "=== Verifying MD5 checksums ==="
FAIL=0
for MD5FILE in *.gz.md5; do
  [[ -f "$MD5FILE" ]] || continue
  if md5sum -c "$MD5FILE"; then
    echo "OK: $MD5FILE"
  else
    echo "FAIL: $MD5FILE" >&2
    FAIL=1
  fi
done

if [[ $FAIL -eq 0 ]]; then
  echo "All checksums passed."
else
  echo "One or more checksums FAILED." >&2
  exit 1
fi
