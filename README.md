---
title: EthnoDock Pro • In-Silico Pharmacognosy
emoji: 🌿
colorFrom: green
colorTo: indigo
sdk: streamlit
sdk_version: 1.30.0
app_file: app.py
pinned: false
license: mit
---

# 🌿 EthnoDock Pro • In-Silico Pharmacognosy & Molecular Docking

An advanced computational ethnopharmacology suite connecting classical **Traditional Chinese Medicine (TCM)** pharmacopoeias (*Shennong Bencaojing*, *Bencao Gangmu*, *Shanghan Lun*) with **Western macromolecular structural docking** (AutoDock Vina, RCSB PDB, UniProt, and RDKit ADMET/PAINS).

## 🚀 Key Modules:
- **Phase 1: Botanical Specimen & Chemical Characterization:** Dual-view macroscopic specimen photographs alongside microscopic 2D molecular structures.
- **Phase 2: Target Cavity Setup:** Automated pocket detection (`smart_cavity_finder`) and AutoDock grid initialization.
- **Phase 3: AutoDock Vina Simulation & 3D WebGL Studio:** 9-pose conformational matrix, thermodynamic inhibition constants ($K_i$), and high-contrast 3Dmol.js WebGL interaction viewer with dashed non-covalent cylinders and 3D residue badges.
- **Phase 4: Classical Paozhi (炮制) Transformation Audit:** Real-time chemical reaction modeling of thermal hydrolysis and detoxification for toxic raw herbs (*Aconitum*, *Pinellia*, *Rhubarb*, *Strychnos*).
- **Phase 5: Semi-Synthetic Bioisostere Lead Optimization:** Rational analog generation with medicinal chemistry rationale.
- **Phase 6: ADMET & Regulatory Dossier Export:** Lipinski Ro5 compliance, PAINS structural alert audit, and 1-click Nature/ACS-grade scientific research dossier export.
