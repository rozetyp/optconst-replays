#!/bin/sh
# Fetch the two certificates at the commits cited on the constant pages and check their SHA-256.
set -e
mkdir -p certs && cd certs
[ -d c1c ] || git clone -q https://github.com/DennisKoshta/c1c-autocorrelation-certificate.git c1c
[ -d c86 ] || git clone -q https://github.com/DennisKoshta/c86-sss-trace-certificate.git c86
(cd c1c && git checkout -q a4ee94ef29ff)
(cd c86 && git checkout -q 90c57242b2f43aa942257154be6731c29898b537)
shasum -a 256 c1c/K24576_0.4103333.json c86/certificate/c86_A19d7_robust_1.8084990.json
