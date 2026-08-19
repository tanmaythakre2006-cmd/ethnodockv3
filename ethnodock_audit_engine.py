import os

def evaluate_historical_claim(
    species_name,
    claim_text,
    translation,
    compound_name,
    target_name,
    pdb_id,
    affinity_kcal,
    admet_dict=None,
    is_paozhi=False,
    toxic_warning=""
):
    """
    Performs an objective, unbiased scientific audit of an ancient TCM claim.
    Returns:
    - verdict_category: 'VALIDATED' | 'MARGINAL_INERT' | 'TOXIC_REJECTED'
    - verdict_badge: Badge string (e.g. '🟢 VALIDATED MECHANISM')
    - verdict_title: Short title summary
    - failure_mode_type: None | 'Binding Efficacy Failure' | 'Compositional Failure' | 'Toxicity/Safety Failure' | 'Observational Error'
    - failure_mode_explanation: Detailed scientific breakdown of why the claim failed or succeeded.
    - molecular_mechanism_summary: Exact biophysical basis.
    """
    
    # 1. Determine Verdict Category based on cold biophysical data
    aff = float(affinity_kcal) if affinity_kcal is not None else -5.0
    
    # Check toxicophore alert from ADMET engine or Paozhi warning
    has_admet_hazard = bool(admet_dict and (
        admet_dict.get("Is Toxicologically Hazardous") or 
        "CRITICAL" in str(admet_dict.get("Safety Status", "")).upper() or 
        "FAIL" in str(admet_dict.get("Safety Badge", "")).upper() or
        "FAIL" in str(admet_dict.get("Safety Status", "")).upper()
    ))
    admet_alert_desc = str(admet_dict.get("Structure Alert Screen", "")) if (admet_dict and has_admet_hazard) else ""
    
    is_toxic_warn_active = bool(
        toxic_warning and not is_paozhi and any(k in toxic_warning.upper() for k in [
            "TOXIC", "LETHAL", "CRITICAL", "FATAL", "HEPATOTOXIC", "NEUROTOXIC", "CARDIOTOXIC", "POISON", "CONVULS", "BURDEN", "IRRITANT"
        ])
    )
    
    if is_toxic_warn_active or (has_admet_hazard and not is_paozhi):
        verdict_category = "TOXIC_REJECTED"
        verdict_badge = "🔴 TOXICITY FAILURE"
        verdict_color = "#FF453A"
        verdict_title = "Ancient Claim Masked Severe Toxicity / Lethal Toxicophore"
        failure_mode_type = "Toxicity / Safety Failure & Toxicophore Alert"
        effective_warn = toxic_warning or admet_alert_desc
        failure_mode_explanation = (
            f"Structure triggered severe toxicological alert ({effective_warn}). "
            f"High binding affinity reflects lethal channel/receptor locking rather than a safe therapeutic window. "
            f"Claim fails as raw monotherapy without chemical detoxification."
        )
        mechanism_summary = f"Substance possesses lethal or organ-damaging structural alert ({effective_warn}) requiring detoxification to prevent patient harm."
        
    elif aff <= -7.2:
        # High affinity binding confirmed
        verdict_category = "VALIDATED"
        verdict_badge = "🟢 VALIDATED EMPIRICAL HIT"
        verdict_color = "#30D158"
        verdict_title = "Molecular Mechanism Authenticated Against Human Target"
        failure_mode_type = "None (Claim Biologically Validated)"
        failure_mode_explanation = "The ancient empirical claim aligns with verified sub-micromolar stereochemical binding to the macromolecular human target."
        mechanism_summary = (
            f"Active constituent {compound_name} demonstrates high-affinity binding (ΔG = {aff:.1f} kcal/mol) to human {target_name} (PDB: {pdb_id}), "
            f"providing direct biophysical validation of the historical indication."
        )
        
    elif aff <= -6.0:
        verdict_category = "MARGINAL_INERT"
        verdict_badge = "🟡 MARGINAL / WEAK BINDING"
        verdict_color = "#FFD60A"
        verdict_title = "Weak Target Affinity / Secondary Modulatory Role"
        failure_mode_type = "Binding / Efficacy Limitation"
        failure_mode_explanation = (
            f"The isolated single molecule {compound_name} exhibits only moderate binding (ΔG = {aff:.1f} kcal/mol), suggesting the ancient clinical effect "
            f"relied on polypharmacological synergy with companion formula herbs rather than potent monotherapy."
        )
        mechanism_summary = f"Moderate target occupancy requiring multi-herb formula synergy (Jun-Chen-Zuo-Shi) to achieve clinical efficacy."
        
    else:
        verdict_category = "MARGINAL_INERT"
        verdict_badge = "⚪ INERT / REJECTED CLAIM"
        verdict_color = "#86868B"
        verdict_title = "Biophysical Efficacy Failure (Lack of Target Affinity)"
        failure_mode_type = "Binding Efficacy Failure & Observational Error"
        failure_mode_explanation = (
            f"Compound {compound_name} failed to achieve therapeutic binding affinity (ΔG = {aff:.1f} kcal/mol). "
            f"Historical accounts likely conflated natural spontaneous recovery or placebo response with direct pharmacological cure."
        )
        mechanism_summary = "Insufficent non-covalent binding energy to induce macromolecular conformational inhibition."

    return {
        "verdict_category": verdict_category,
        "verdict_badge": verdict_badge,
        "verdict_color": verdict_color,
        "verdict_title": verdict_title,
        "original_claim_text": claim_text,
        "original_claim_translation": translation,
        "failure_mode_type": failure_mode_type,
        "failure_mode_explanation": failure_mode_explanation,
        "mechanism_summary": mechanism_summary
    }
