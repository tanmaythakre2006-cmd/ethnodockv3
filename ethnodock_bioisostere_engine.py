import re
from rdkit import Chem
from rdkit.Chem import AllChem

TRANSFORMATION_RULES = [
    {
        "name": "Fluorination (-OH -> -F)",
        "smarts": "[OX2H]",
        "replace_smiles": "F",
        "rationale": "Replaces vulnerable hydroxyl with fluorine to block Phase II glucuronidation while mimicking oxygen electronegativity."
    },
    {
        "name": "Chlorination (-OH -> -Cl)",
        "smarts": "[OX2H]",
        "replace_smiles": "Cl",
        "rationale": "Increases lipophilicity (LogP) and engages in halogen bonding with receptor aromatic residues."
    },
    {
        "name": "O-Methylation (-OH -> -OCH3)",
        "smarts": "[OX2H]",
        "replace_smiles": "OC",
        "rationale": "A common metabolic transformation that increases membrane permeability and removes polar H-bond donor penalties."
    },
    {
        "name": "Amination (-OH -> -NH2)",
        "smarts": "[OX2H]",
        "replace_smiles": "N",
        "rationale": "Converts hydroxyl to amine to create a stronger basic H-bond donor center for acidic binding pockets."
    },
    {
        "name": "Amidation (-COOH -> -CONH2)",
        "smarts": "C(=O)[OX2H1]",
        "replace_smiles": "C(=O)N",
        "rationale": "Neutralizes anionic carboxylate charge, reducing renal clearance and enhancing cellular target entry."
    },
    {
        "name": "Esterification (-COOH -> -COOCH3)",
        "smarts": "C(=O)[OX2H1]",
        "replace_smiles": "C(=O)OC",
        "rationale": "Classic prodrug strategy to enhance lipophilicity and gastrointestinal absorption."
    },
    {
        "name": "N-Methylation (-NH2 -> -NHCH3)",
        "smarts": "[NX3H2]",
        "replace_smiles": "NC",
        "rationale": "Enhances metabolic stability against monoamine oxidases and optimizes steric cavity fit."
    }
]

def generate_tcm_derivatives(smiles, max_variants=8):
    """
    Generates rationally designed bioisosteric derivatives and semi-synthetic analogs
    from a parent TCM phytochemical SMILES string.
    
    Returns a list of dicts with variant SMILES, modification name, and medicinal chemistry rationale.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return []

    parent_smiles = Chem.MolToSmiles(mol)
    results = []
    seen_smiles = {parent_smiles}

    for rule in TRANSFORMATION_RULES:
        try:
            pattern = Chem.MolFromSmarts(rule["smarts"])
            replacement = Chem.MolFromSmiles(rule["replace_smiles"])
            
            if pattern and replacement and mol.HasSubstructMatch(pattern):
                replacements = Chem.ReplaceSubstructs(mol, pattern, replacement, replaceAll=False)
                for r in replacements:
                    try:
                        Chem.SanitizeMol(r)
                        var_smiles = Chem.MolToSmiles(r)
                        if var_smiles not in seen_smiles:
                            seen_smiles.add(var_smiles)
                            results.append({
                                "variant_smiles": var_smiles,
                                "name": rule["name"],
                                "rationale": rule["rationale"],
                                "parent_smiles": parent_smiles
                            })
                            if len(results) >= max_variants:
                                return results
                    except Exception:
                        pass
        except Exception as e:
            print(f"Rule {rule['name']} processing note: {e}")

    return results
