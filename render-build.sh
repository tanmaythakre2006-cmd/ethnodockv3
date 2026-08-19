#!/usr/bin/env bash
# Exit on any error
set -o errexit

echo "============================================================"
echo "🌿 EthnoDock Pro • Flagship Production Deployment on Render"
echo "============================================================"

# Upgrade pip & install production wheels
python -m pip install --upgrade pip setuptools wheel
pip install -r requirements.txt

echo "============================================================"
echo "✅ EthnoDock Pro Flagship Environment Prepared Successfully!"
echo "============================================================"
