from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED, FilterCatalog

# Initialize global PAINS catalog to avoid re-allocating on every query
_pains_params = FilterCatalog.FilterCatalogParams()
_pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
_pains_catalog = FilterCatalog.FilterCatalog(_pains_params)

def get_admet_profile(smiles):
    """
    Computes comprehensive ADMET, drug-likeness, and toxicology parameters
    including Lipinski Rule of Five, Veber rules, QED, and PAINS structural alerts.
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

    # 3. PAINS Screening
    alert_match = _pains_catalog.GetFirstMatch(mol)
    pains_alert = alert_match.GetDescription() if alert_match else "Pass (0 PAINS Alerts)"
    is_pains_clean = alert_match is None

    # 4. Blood-Brain Barrier (BBB) Estimation Heuristic
    # Typically: MW < 400, LogP 2-4, TPSA < 90, HBD < 3
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
        "Predicted BBB Permeability": "High" if bbb_permeable else "Moderate/Low"
    }
