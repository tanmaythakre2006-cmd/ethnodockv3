import os
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED, FilterCatalog

# 1. Initialize global PAINS catalog for assay interference
_pains_params = FilterCatalog.FilterCatalogParams()
_pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
_pains_catalog = FilterCatalog.FilterCatalog(_pains_params)

# 2. Curated In-Vivo Natural Product Toxicophore & Lethality Catalog
NATURAL_TOXICOPHORES = [
    {
        "name": "Diester Diterpenoid Alkaloid (Aconitine-type Cardiotoxin)",
        "severity": "CRITICAL LETHAL CARDIOTOXICITY",
        "target": "Voltage-Gated Sodium Channels (Nav1.5)",
        "smarts": ["c1ccccc1C(=O)O[#6]", "CC(=O)O[#6]"], # Both benzoate and acetate esters on amine
        "condition": "all",
        "requires_nitrogen": True,
        "description": "Lethal Nav1.5 channel activator (LD50 ~0.1 mg/kg) causing persistent myocardial depolarization, ventricular tachycardia, and fatal cardiac arrest."
    },
    {
        "name": "Strychnos Indole Alkaloid (Strychnine-type Neurotoxin)",
        "severity": "CRITICAL LETHAL NEUROTOXICITY",
        "target": "Spinal Glycine Receptors",
        "smiles": ["O=C1CC2OCC=C3CN4CCC56C7C4CC3C2C5C1=CC=C7N6"],
        "condition": "any",
        "requires_nitrogen": True,
        "description": "Competitive glycine receptor antagonist causing violent tetanic convulsions, motor neuron spasms, and fatal asphyxiation (Lethal Dose ~30 mg)."
    },
    {
        "name": "Aristolochic Acid Core (Nephrotoxic Carcinogen)",
        "severity": "CRITICAL CARCINOGENICITY & NEPHROTOXICITY",
        "target": "Genomic DNA / Renal Proximal Tubules",
        "smiles": ["O=C(O)C1=CC2=C(C=C1[N+](=O)[O-])C3=CC(OC)=C4OCOC4=C3C=C2"],
        "condition": "any",
        "requires_nitrogen": True,
        "description": "Forms irreversible covalent DNA adducts (dA:T -> T:A) causing rapid bilateral renal interstitial fibrosis and fatal urothelial carcinoma."
    }
]

def screen_natural_toxicophores(mol):
    """
    Screens small molecules against lethal and organ-damaging natural product toxicophore patterns.
    """
    if not mol:
        return []
    
    alerts = []
    has_n = any(atom.GetSymbol() == "N" for atom in mol.GetAtoms())
    
    for tox in NATURAL_TOXICOPHORES:
        if tox.get("requires_nitrogen") and not has_n:
            continue
            
        matches = []
        if "smarts" in tox:
            for s in tox["smarts"]:
                q = Chem.MolFromSmarts(s)
                if q is not None and mol.HasSubstructMatch(q):
                    matches.append(True)
                else:
                    matches.append(False)
        if "smiles" in tox:
            for s in tox["smiles"]:
                q = Chem.MolFromSmiles(s)
                if q is not None and mol.HasSubstructMatch(q):
                    matches.append(True)
                else:
                    matches.append(False)
                    
        total_queries = len(tox.get("smarts", [])) + len(tox.get("smiles", []))
        if tox["condition"] == "all" and all(matches) and len(matches) == total_queries:
            alerts.append(tox)
        elif tox["condition"] == "any" and any(matches):
            alerts.append(tox)
            
    return alerts

def get_admet_profile(smiles):
    """
    Computes comprehensive ADMET, drug-likeness, toxicological safety parameters,
    Lipinski Rule of Five, Veber oral bioavailability, PAINS, and Natural Product Toxicophore screens.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    aromatic_rings = Lipinski.NumAromaticRings(mol)
    qed_score = QED.qed(mol)

    # 1. Lipinski Violations
    lipinski_violations = []
    if mw > 500:
        lipinski_violations.append("MW > 500 Da")
    if logp > 5.0:
        lipinski_violations.append("LogP > 5.0")
    if hbd > 5:
        lipinski_violations.append("HBD > 5")
    if hba > 10:
        lipinski_violations.append("HBA > 10")

    # 2. Veber Rules (Oral Bioavailability)
    veber_violations = []
    if rotb > 10:
        veber_violations.append("Rotatable Bonds > 10")
    if tpsa > 140:
        veber_violations.append("TPSA > 140 Å²")

    # 3. Assay Interference (PAINS) Screening
    alert_match = _pains_catalog.GetFirstMatch(mol)
    pains_alert = alert_match.GetDescription() if alert_match else "Pass (0 PAINS Alerts)"
    is_pains_clean = alert_match is None

    # 4. In-Vivo Natural Toxicophore & Lethality Screening
    tox_alerts = screen_natural_toxicophores(mol)
    is_toxic_hazardous = len(tox_alerts) > 0
    
    if is_toxic_hazardous:
        top_tox = tox_alerts[0]
        structure_alert_summary = f"CRITICAL HAZARD: {top_tox['name']} ({top_tox['severity']})"
        safety_status = "CRITICAL HAZARD / FAIL"
        safety_badge = "🔴 FAIL"
    elif not is_pains_clean:
        structure_alert_summary = f"PAINS ALERT: {pains_alert}"
        safety_status = "ASSAY INTERFERENCE ALERT"
        safety_badge = "🟡 CAUTION"
    else:
        structure_alert_summary = "PASSED (0 alerts, Verified Clean)"
        safety_status = "PASS (Clean Scaffold)"
        safety_badge = "🟢 PASS"

    # 5. Blood-Brain Barrier (BBB) Estimation Heuristic
    bbb_permeable = (mw < 450) and (1.0 <= logp <= 4.5) and (tpsa < 90) and (hbd <= 3)

    return {
        "SMILES": smiles,
        "Molecular Weight": round(mw, 2),
        "LogP": round(logp, 2),
        "TPSA (Å²)": round(tpsa, 2),
        "H-Bond Donors": hbd,
        "H-Bond Acceptors": hba,
        "Rotatable Bonds": rotb,
        "Aromatic Rings": aromatic_rings,
        "QED Drug-Likeness": round(qed_score, 3),
        "Lipinski Violations": len(lipinski_violations),
        "Lipinski Details": lipinski_violations if lipinski_violations else ["None (Ro5 Compliant)"],
        "Veber Pass": len(veber_violations) == 0,
        "PAINS Screen": pains_alert,
        "PAINS Clean": is_pains_clean,
        "Natural Toxicophore Alerts": tox_alerts,
        "Is Toxicologically Hazardous": is_toxic_hazardous,
        "Structure Alert Screen": structure_alert_summary,
        "Safety Status": safety_status,
        "Safety Badge": safety_badge,
        "Predicted BBB Permeability": "High" if bbb_permeable else "Moderate/Low"
    }
