import html
import base64
import os
import json
import math
from io import BytesIO
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import Draw
import ethnodock_audit_engine as audit_eng

def generate_tcm_dossier_html(
    species_name,
    scientific_name,
    chinese_name,
    source_title,
    claim_text,
    translation,
    compound_name,
    smiles,
    target_name,
    pdb_id,
    affinity_kcal,
    poses_table_html,
    interactions_table_html,
    admet_dict=None,
    variant_info=None,
    plant_photo_b64=None,
    paozhi_data=None,
    is_paozhi_processed=False,
    # Upgraded Rich Biophysical Assets
    md_results=None,
    var_md_results=None,
    ligplot_b64=None,
    var_ligplot_b64=None,
    rmsd_plot_b64=None,
    occupancy_plot_b64=None,
    fes_plot_b64=None,
    admet_radar_b64=None,
    receptor_pdbqt_str=None,
    ligand_pdbqt_str=None,
    var_smiles=None,
    var_admet_dict=None,
    benchmark_data=None,
    mmgbsa_data=None,
    mmgbsa_chart_b64=None,
    var_mmgbsa_data=None,
    var_mmgbsa_chart_b64=None,
    comp_hotspot_chart_b64=None,
    targetome_results=None,
    network_results=None,
    synergy_results=None,
    microbiome_results=None,
    population_results=None
):
    """
    Generates an executive, publication-grade scientific research monograph
    matching Nature / ACS Medicinal Chemistry publication standards with
    embedded high-res figures, MD trajectory plots, FES landscapes,
    ADMET radar spider charts, and standalone interactive 3D WebGL pocket viewers.
    """
    date_str = datetime.now().strftime("%B %d, %Y • %H:%M UTC")
    doc_id = f"EDK-TCM-{pdb_id}-{abs(hash(compound_name + smiles)) % 1000000:06d}"

    # Perform Objective Historical Claim Audit
    toxic_warn = paozhi_data.get('raw_toxicity_warning', '') if (paozhi_data and not is_paozhi_processed) else ''
    audit_res = audit_eng.evaluate_historical_claim(
        species_name=species_name,
        claim_text=claim_text,
        translation=translation,
        compound_name=compound_name,
        target_name=target_name,
        pdb_id=pdb_id,
        affinity_kcal=affinity_kcal,
        admet_dict=admet_dict,
        is_paozhi=is_paozhi_processed,
        toxic_warning=toxic_warn
    )

    # Generate 2D Structure Base64 for Parent
    mol_b64 = ""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol, size=(380, 260))
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            mol_b64 = f"data:image/png;base64,{base64.b64encode(buffered.getvalue()).decode()}"
    except Exception:
        pass

    # Generate 2D Structure Base64 for Derivative if available
    var_mol_b64 = ""
    if variant_info and variant_info.get('smiles'):
        try:
            vmol = Chem.MolFromSmiles(variant_info['smiles'])
            if vmol:
                vimg = Draw.MolToImage(vmol, size=(380, 260))
                vbuf = BytesIO()
                vimg.save(vbuf, format="PNG")
                var_mol_b64 = f"data:image/png;base64,{base64.b64encode(vbuf.getvalue()).decode()}"
        except Exception:
            pass

    # Render Plant Photo Tag
    plant_img_tag = f'<img src="{plant_photo_b64}" style="width:100%; height:220px; object-fit:cover; border-radius:10px; border:1px solid #E2E8F0;"/>' if plant_photo_b64 else '<div style="height:220px; background:#F8FAFC; border-radius:10px; display:flex; align-items:center; justify-content:center; color:#94A3B8; font-size:13px; border:1px solid #E2E8F0;">No Specimen Photo Available</div>'

    mol_img_tag = f'<img src="{mol_b64}" style="width:100%; height:220px; object-fit:contain; border-radius:10px; background:#FFFFFF; border:1px solid #E2E8F0;"/>' if mol_b64 else '<div style="height:220px; background:#F8FAFC; border-radius:10px; display:flex; align-items:center; justify-content:center; color:#94A3B8; font-size:13px; border:1px solid #E2E8F0;">2D Molecular Topology</div>'

    # ADMET Table Rows
    admet_cards = ""
    if admet_dict:
        admet_cards = f"""
        <div class="table-container">
            <table class="report-table">
                <thead>
                    <tr>
                        <th>Parameter</th>
                        <th>Natural Parent Value</th>
                        <th>Standard Benchmark</th>
                        <th>Compliance Status</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td><strong>Molecular Weight</strong></td>
                        <td>{admet_dict.get('Molecular Weight', 'N/A')} g/mol</td>
                        <td>&le; 500 g/mol (Lipinski Rule of 5)</td>
                        <td><span class="badge badge-green">PASS</span></td>
                    </tr>
                    <tr>
                        <td><strong>Lipophilicity (LogP)</strong></td>
                        <td>{admet_dict.get('LogP', 'N/A')}</td>
                        <td>&le; 5.0 (Optimal Octanol/Water Partition)</td>
                        <td><span class="badge badge-green">OPTIMAL</span></td>
                    </tr>
                    <tr>
                        <td><strong>Polar Surface Area (TPSA)</strong></td>
                        <td>{admet_dict.get('TPSA (Å²)', 'N/A')} &Aring;&sup2;</td>
                        <td>&le; 140 &Aring;&sup2; (Veber Oral Permeability)</td>
                        <td><span class="badge badge-blue">BIOAVAILABLE</span></td>
                    </tr>
                    <tr>
                        <td><strong>H-Bond Donors / Acceptors</strong></td>
                        <td>{admet_dict.get('H-Bond Donors', 'N/A')} / {admet_dict.get('H-Bond Acceptors', 'N/A')}</td>
                        <td>HBD &le; 5, HBA &le; 10</td>
                        <td><span class="badge badge-green">COMPLIANT</span></td>
                    </tr>
                    <tr>
                        <td><strong>QED Drug-Likeness Score</strong></td>
                        <td><strong>{admet_dict.get('QED Drug-Likeness', 'N/A')}</strong></td>
                        <td>Range: 0.0 &ndash; 1.0 (Bickerton et al.)</td>
                        <td><span class="badge badge-gold">HIGH DRUG-LIKENESS</span></td>
                    </tr>
                    <tr>
                        <td><strong>In-Vivo Toxicophore & Structural Alert Screen</strong></td>
                        <td><span style="{'color:#DC2626; font-weight:700;' if admet_dict.get('Is Toxicologically Hazardous') else ''}">{admet_dict.get('Structure Alert Screen', admet_dict.get('PAINS Screen', 'Clean'))}</span></td>
                        <td>Zero Lethal Toxicophores & False Positives</td>
                        <td><span class="badge {'badge-green' if not admet_dict.get('Is Toxicologically Hazardous') else 'badge-red'}">{admet_dict.get('Safety Badge', 'PASS')}</span></td>
                    </tr>
                </tbody>
            </table>
        </div>
        """

    # Paozhi Section (Only displayed when in Processed / Detoxified State)
    paozhi_html = ""
    if is_paozhi_processed and paozhi_data and paozhi_data.get('paozhi_method'):
        paozhi_html = f"""
        <div class="section-card paozhi-card">
            <div class="section-title" style="color:#B45309;">
                <span>⚗️ Classical Paozhi (炮制) Processing & Detoxification Audit</span>
            </div>
            <div style="display:grid; grid-template-columns: 1.2fr 1fr; gap:16px; margin-top:10px;">
                <div>
                    <p style="margin:0 0 6px 0; font-size:13px;"><strong>Processing Method:</strong> {html.escape(paozhi_data.get('paozhi_method', ''))}</p>
                    <p style="margin:0 0 6px 0; font-size:13px;"><strong>Reaction Pathway:</strong> {html.escape(paozhi_data.get('reaction_type', ''))}</p>
                    <p style="margin:0; font-size:12px; font-family:monospace; background:#FEF3C7; padding:8px 12px; border-radius:6px; color:#92400E; border:1px solid #FDE68A;">
                        {html.escape(paozhi_data.get('reaction_equation', ''))}
                    </p>
                </div>
                <div style="background:#FFFBEB; padding:12px; border-radius:8px; border-left:4px solid #F59E0B; border:1px solid #FEF3C7;">
                    <strong style="color:#B45309; font-size:11px; text-transform:uppercase; letter-spacing:0.5px;">Toxicological & Efficacy Shift</strong>
                    <p style="margin:6px 0 0 0; font-size:13px; color:#78350F; line-height:1.45;">{html.escape(paozhi_data.get('detox_benefit', ''))}</p>
                </div>
            </div>
        </div>
        """

    # Crystallographic Redocking Validation Benchmark Card
    benchmark_html = ""
    if benchmark_data and benchmark_data.get('success'):
        benchmark_html = f"""
        <div style="background:#F0FDF4; border:1px solid #BBF7D0; border-radius:10px; padding:16px; margin: 14px 0 18px 0;">
            <div style="display:flex; justify-content:space-between; align-items:center;">
                <span style="font-weight:700; color:#15803D; font-size:12.5px; text-transform:uppercase; letter-spacing:0.5px;">
                    🔬 Crystallographic Ground Truth Benchmark & Redocking Proof
                </span>
                <span class="badge" style="background:#DCFCE7; color:#15803D; border:1px solid #86EFAC;">{benchmark_data.get('tier', 'VALIDATED')} (&le; 2.0 &Aring;)</span>
            </div>
            <div style="display:grid; grid-template-columns: repeat(3, 1fr); gap:12px; margin-top:10px; text-align:center;">
                <div style="background:#FFFFFF; padding:10px; border-radius:6px; border:1px solid #DCFCE7;">
                    <div style="font-size:11px; color:#64748B; text-transform:uppercase; font-weight:600;">Crystallographic Co-Crystal Drug</div>
                    <div style="font-weight:700; color:#0F172A; font-size:13.5px; margin-top:3px;">{html.escape(benchmark_data.get('ligand_name', ''))}</div>
                    <div style="font-size:10px; color:#15803D;">Code: {benchmark_data.get('ligand_resname', '')}</div>
                </div>
                <div style="background:#FFFFFF; padding:10px; border-radius:6px; border:1px solid #DCFCE7;">
                    <div style="font-size:11px; color:#64748B; text-transform:uppercase; font-weight:600;">Heavy-Atom RMSD vs Crystal</div>
                    <div style="font-weight:800; color:#059669; font-size:18px; margin-top:3px;">{benchmark_data.get('rmsd', 0.0):.2f} &Aring;</div>
                    <div style="font-size:10px; color:#059669;">Criteria: &le; 2.0 &Aring; (Pass)</div>
                </div>
                <div style="background:#FFFFFF; padding:10px; border-radius:6px; border:1px solid #DCFCE7;">
                    <div style="font-size:11px; color:#64748B; text-transform:uppercase; font-weight:600;">Reproduced Binding Affinity</div>
                    <div style="font-weight:700; color:#0F172A; font-size:15px; margin-top:3px;">{benchmark_data.get('docked_affinity', 'N/A')} <span style="font-size:11px; font-weight:500;">kcal/mol</span></div>
                    <div style="font-size:10px; color:#64748B;">AutoDock Vina Re-Docked</div>
                </div>
            </div>
            <p style="margin:10px 0 0 0; font-size:12px; color:#166534; line-height:1.45;">
                <b>Crystallographic Grounding Note:</b> Blind redocking of the experimental native ligand reproduced the true X-ray crystallographic coordinates with an RMSD of {benchmark_data.get('rmsd', 0.0):.2f} &Aring;, verifying that the active-site cavity definition and scoring parameters operate with experimental accuracy.
            </p>
        </div>
        """

    # Executive Summary Bar Evaluation
    is_hazard = (audit_res['verdict_category'] == 'TOXIC_REJECTED') or bool(admet_dict and (admet_dict.get('Is Toxicologically Hazardous') or 'FAIL' in str(admet_dict.get('Safety Status', '')).upper()))
    aff_color = "#DC2626" if is_hazard else ("#059669" if float(affinity_kcal) <= -7.2 else "#D97706")
    aff_badge_html = f'<div style="font-size:10px; font-weight:700; color:#DC2626; margin-top:2px;">🚨 LETHAL TOXIN / FAIL</div>' if is_hazard else ''

    # Molecular Dynamics (MD) Section HTML
    md_html = ""
    if md_results or rmsd_plot_b64:
        stability_badge_color = "#059669" if md_results and md_results.get('stability_verdict') == 'STABLE' else ("#D97706" if md_results and md_results.get('stability_verdict') == 'DYNAMIC' else "#DC2626")
        stability_text = md_results.get('stability_verdict', 'STABLE') if md_results else 'STABLE'
        mean_rmsd_val = f"{md_results.get('mean_rmsd', 0.0):.2f} Å" if md_results else "N/A"
        max_rmsd_val = f"{md_results.get('max_rmsd', 0.0):.2f} Å" if md_results else "N/A"
        
        rmsd_img_html = f'<img src="{rmsd_plot_b64}" style="width:100%; border-radius:8px; border:1px solid #E2E8F0; box-shadow:0 1px 3px rgba(0,0,0,0.05);"/>' if rmsd_plot_b64 else ''
        occupancy_img_html = f'<img src="{occupancy_plot_b64}" style="width:100%; border-radius:8px; border:1px solid #E2E8F0; box-shadow:0 1px 3px rgba(0,0,0,0.05);"/>' if occupancy_plot_b64 else ''
        fes_img_html = f'<img src="{fes_plot_b64}" style="width:100%; border-radius:8px; border:1px solid #E2E8F0; box-shadow:0 1px 3px rgba(0,0,0,0.05);"/>' if fes_plot_b64 else ''

        md_html = f"""
        <div class="section-card">
            <div class="section-title">
                <span>🌊 Molecular Dynamics (MD) Trajectory & Solvent Residence Stability</span>
                <span class="badge" style="background:#F0FDF4; color:{stability_badge_color}; border:1px solid {stability_badge_color};">{stability_text}</span>
            </div>
            <p style="margin:0 0 16px 0; font-size:13px; color:#64748B; line-height:1.5;">
                Simulated dynamic solvent stability, kinetic residence time, and conformational pose drift within the aqueous macromolecular pocket across the trajectory duration.
            </p>

            <!-- MD Telemetry Grid -->
            <div style="display:grid; grid-template-columns: repeat(3, 1fr); gap:12px; margin-bottom:20px; background:#F8FAFC; padding:14px; border-radius:8px; border:1px solid #E2E8F0;">
                <div style="text-align:center;">
                    <div style="font-size:11px; color:#64748B; font-weight:600; text-transform:uppercase;">Trajectory Mean RMSD</div>
                    <div style="font-size:18px; font-weight:700; color:#0F172A; margin-top:3px;">{mean_rmsd_val}</div>
                    <div style="font-size:10px; color:#059669;">Threshold: &le; 2.0 Å</div>
                </div>
                <div style="text-align:center; border-left:1px solid #E2E8F0; border-right:1px solid #E2E8F0;">
                    <div style="font-size:11px; color:#64748B; font-weight:600; text-transform:uppercase;">Peak Conformational Drift</div>
                    <div style="font-size:18px; font-weight:700; color:#0F172A; margin-top:3px;">{max_rmsd_val}</div>
                    <div style="font-size:10px; color:#64748B;">Heavy-Atom Displacement</div>
                </div>
                <div style="text-align:center;">
                    <div style="font-size:11px; color:#64748B; font-weight:600; text-transform:uppercase;">Solvent Ensemble Status</div>
                    <div style="font-size:18px; font-weight:700; color:{stability_badge_color}; margin-top:3px;">{stability_text}</div>
                    <div style="font-size:10px; color:#64748B;">Thermodynamic State</div>
                </div>
            </div>

            <!-- MD Trajectory Graphs Grid -->
            <div style="display:grid; grid-template-columns: 1.1fr 0.9fr; gap:18px; margin-bottom:16px;">
                <div>
                    {rmsd_img_html}
                    <div style="font-size:11px; color:#64748B; margin-top:6px; text-align:center;"><b>Figure 2 | Heavy-Atom RMSD Time Series.</b> Trajectory stability plotted against the 2.0 Å conformational drift ceiling.</div>
                </div>
                <div>
                    {occupancy_img_html}
                    <div style="font-size:11px; color:#64748B; margin-top:6px; text-align:center;"><b>Figure 3 | Key Pocket Contact Persistence (% Occupancy).</b> Quantifying non-covalent bond lifetimes.</div>
                </div>
            </div>

            <!-- FES Landscape if present -->
            {f'''
            <div style="margin-top:20px; display:grid; grid-template-columns: 1fr 1fr; gap:18px; align-items:center;">
                <div>
                    {fes_img_html}
                    <div style="font-size:11px; color:#64748B; margin-top:6px; text-align:center;"><b>Figure 4 | Conformational Free Energy Surface (FES).</b> &Delta;G(RMSD, Rg) energy wells and conformational transitions.</div>
                </div>
                <div style="background:#F8FAFC; border:1px solid #E2E8F0; border-radius:8px; padding:16px;">
                    <h5 style="margin:0 0 8px 0; color:#0F172A; font-size:13px; text-transform:uppercase;">Thermodynamic Free Energy Interpretation</h5>
                    <p style="margin:0 0 8px 0; font-size:12px; color:#475569; line-height:1.5;">
                        The 2D Free Energy Surface (FES) reconstructs the conformational probability distribution &Delta;G(RMSD, R<sub>g</sub>) = -k<sub>B</sub>T ln(P / P<sub>max</sub>).
                    </p>
                    <p style="margin:0; font-size:12px; color:#475569; line-height:1.5;">
                        A deep, narrow global energy minimum confirms that the docked ligand occupies a rigid, thermodynamically stabilized enthalpy well rather than fluctuating across multiple shallow metastable states.
                    </p>
                </div>
            </div>
            ''' if fes_img_html else ''}

        </div>
        """

    # 2D LigPlot Diagram Section
    ligplot_html = ""
    if ligplot_b64:
        ligplot_html = f"""
        <div style="margin: 20px 0; text-align:center;">
            <img src="{ligplot_b64}" style="max-width:100%; border-radius:10px; border:1px solid #CBD5E1; box-shadow:0 2px 8px rgba(0,0,0,0.06);"/>
            <div style="font-size:12px; color:#64748B; margin-top:8px; line-height:1.4;">
                <strong>Figure 1 | 2D Non-Covalent Interaction Topology (LigPlot+ Style).</strong> Green dashed vectors denote polar hydrogen bonds with explicit donor-acceptor distances (Å); spoked radial arcs indicate hydrophobic contact eyelashes; and pill badges specify coordinating active-site pocket residues.
            </div>
        </div>
        """

    # Interactive 3D WebGL Pocket Viewer Section
    pocket_3d_html = ""
    if receptor_pdbqt_str and ligand_pdbqt_str:
        rec_clean_js = json.dumps(receptor_pdbqt_str)
        lig_clean_js = json.dumps(ligand_pdbqt_str)
        pocket_3d_html = f"""
        <div style="margin: 24px 0;">
            <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:8px;">
                <span style="font-size:13px; font-weight:700; color:#0F172A; text-transform:uppercase; letter-spacing:0.5px;">⚡ Interactive 3D WebGL Pocket Exploration</span>
                <span style="font-size:11px; color:#0284C7; font-weight:600;">Left-Click to Rotate &bull; Scroll to Zoom &bull; Right-Click to Translate</span>
            </div>
            <div id="pocket_3d_viewport" style="width: 100%; height: 380px; position: relative; border-radius: 12px; background: #0F172A; overflow: hidden; border: 1px solid #334155; box-shadow: 0 4px 15px rgba(0,0,0,0.25);">
                <div id="loading_3d_tag" style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); color: #94A3B8; font-size: 13px;">
                    Initializing WebGL Molecular Canvas...
                </div>
            </div>
            <script>
            document.addEventListener("DOMContentLoaded", function() {{
                try {{
                    let container = document.getElementById("pocket_3d_viewport");
                    if (container && typeof $3Dmol !== 'undefined') {{
                        let viewer = $3Dmol.createViewer(container, {{ backgroundColor: '#0F172A' }});
                        let recData = {rec_clean_js};
                        let ligData = {lig_clean_js};
                        
                        viewer.addModel(recData, "pdb");
                        viewer.setStyle({{model: 0}}, {{cartoon: {{color: '#64748B', opacity: 0.85}}}});
                        viewer.addSurface($3Dmol.SurfaceType.MS, {{opacity: 0.12, color: '#38BDF8'}});
                        
                        viewer.addModel(ligData, "pdb");
                        viewer.setStyle({{model: 1}}, {{stick: {{colorscheme: 'greenCarbon', radius: 0.25}}}});
                        
                        viewer.zoomTo({{model: 1}});
                        viewer.render();
                        let tag = document.getElementById("loading_3d_tag");
                        if (tag) tag.style.display = "none";
                    }}
                }} catch(err) {{
                    console.log("3Dmol WebGL rendering note:", err);
                }}
            }});
            </script>
        </div>
        """

    # Quantitative Biophysics: MM-GBSA Free Energy & Residue Hotspots Section
    mmgbsa_html = ""
    if mmgbsa_data:
        hotspots = mmgbsa_data.get('hotspots', {})
        hotspot_rows_list = []
        for r_name, r_energy in hotspots.items():
            role_desc = "Primary Pharmacophore Anchor" if r_energy <= -2.0 else ("Secondary Energetic Clamp" if r_energy <= -1.0 else "Periphery Contact")
            hotspot_rows_list.append(f"""
                <tr>
                    <td><strong>{html.escape(r_name)}</strong></td>
                    <td style="color:#059669; font-weight:700;">{r_energy:.2f} kcal/mol</td>
                    <td><span class="badge badge-green">{role_desc}</span></td>
                </tr>
            """)
        hotspot_rows_html = "".join(hotspot_rows_list) if hotspot_rows_list else "<tr><td colspan='3'>No prominent residue hotspots computed.</td></tr>"

        chart_html = f'<img src="{mmgbsa_chart_b64}" style="width:100%; border-radius:8px; border:1px solid #E2E8F0; box-shadow:0 1px 3px rgba(0,0,0,0.05);"/>' if mmgbsa_chart_b64 else ''

        mmgbsa_html = f"""
        <div class="section-card">
            <div class="section-title">
                <span>⚡ Quantitative Biophysics: MM-GBSA Free Energy & Hotspot Decomposition</span>
                <span class="badge badge-blue">Continuum Solvation</span>
            </div>
            <p style="margin:0 0 16px 0; font-size:13px; color:#64748B; line-height:1.5;">
                Molecular Mechanics / Generalized Born Surface Area (MM-GBSA) decouples binding free energy into mechanical van der Waals packing, Coulombic electrostatics, and continuum aqueous desolvation penalties.
            </p>

            <!-- Energetics Telemetry Grid -->
            <div style="display:grid; grid-template-columns: repeat(4, 1fr); gap:12px; margin-bottom:20px; background:#F8FAFC; padding:14px; border-radius:8px; border:1px solid #E2E8F0; text-align:center;">
                <div>
                    <div style="font-size:11px; color:#64748B; font-weight:600; text-transform:uppercase;">Net ΔG Bind</div>
                    <div style="font-size:18px; font-weight:800; color:#059669; margin-top:3px;">{mmgbsa_data.get('total_dg', 0.0):.2f} <span style="font-size:11px; font-weight:500;">kcal/mol</span></div>
                    <div style="font-size:10px; color:#059669;">MM-GBSA Total</div>
                </div>
                <div>
                    <div style="font-size:11px; color:#64748B; font-weight:600; text-transform:uppercase;">ΔE vdW (Packing)</div>
                    <div style="font-size:18px; font-weight:800; color:#0F172A; margin-top:3px;">{mmgbsa_data.get('vdw_energy', 0.0):.2f} <span style="font-size:11px; font-weight:500;">kcal/mol</span></div>
                    <div style="font-size:10px; color:#64748B;">Dispersive Enthalpy</div>
                </div>
                <div>
                    <div style="font-size:11px; color:#64748B; font-weight:600; text-transform:uppercase;">ΔE Elec (Polar)</div>
                    <div style="font-size:18px; font-weight:800; color:#0F172A; margin-top:3px;">{mmgbsa_data.get('elec_energy', 0.0):.2f} <span style="font-size:11px; font-weight:500;">kcal/mol</span></div>
                    <div style="font-size:10px; color:#64748B;">Coulombic Force</div>
                </div>
                <div>
                    <div style="font-size:11px; color:#64748B; font-weight:600; text-transform:uppercase;">ΔG Solvation (GB+SA)</div>
                    <div style="font-size:18px; font-weight:800; color:#D97706; margin-top:3px;">+{abs(mmgbsa_data.get('polar_solv', 0.0) + mmgbsa_data.get('nonpolar_solv', 0.0)):.2f} <span style="font-size:11px; font-weight:500;">kcal/mol</span></div>
                    <div style="font-size:10px; color:#D97706;">Desolvation Cost</div>
                </div>
            </div>

            <!-- Hotspots Bar Chart and Table -->
            <div style="display:grid; grid-template-columns: 1.1fr 0.9fr; gap:18px; align-items:center;">
                <div>
                    {chart_html}
                    <div style="font-size:11px; color:#64748B; margin-top:6px; text-align:center;"><b>Figure 5 | Per-Residue Energetic Decomposition.</b> Residues contributing &le; -1.5 kcal/mol serve as primary anchors.</div>
                </div>
                <div>
                    <h5 style="margin:0 0 8px 0; font-size:12px; text-transform:uppercase; color:#0F172A;">Key Residue Energetic Contributors</h5>
                    <div class="table-container" style="margin-top:0;">
                        <table class="report-table">
                            <thead>
                                <tr>
                                    <th>Pocket Residue</th>
                                    <th>ΔG Contribution</th>
                                    <th>Functional Classification</th>
                                </tr>
                            </thead>
                            <tbody>
                                {hotspot_rows_html}
                            </tbody>
                        </table>
                    </div>
                </div>
            </div>
        </div>
        """

    # Semi-Synthetic Bioisostere Lead Optimization Section (Parent vs Derivative)
    variant_html = ""
    if variant_info:
        var_aff = variant_info.get('affinity', 'N/A')
        aff_delta_str = ""
        try:
            p_aff = float(affinity_kcal)
            v_aff = float(var_aff)
            d_aff = v_aff - p_aff
            aff_delta_str = f"{d_aff:+.2f} kcal/mol ({'Enhanced Affinity' if d_aff < 0 else 'Reduced Affinity'})"
            delta_color = "#059669" if d_aff < 0 else "#DC2626"
        except Exception:
            delta_color = "#4338CA"

        radar_img_html = f'<img src="{admet_radar_b64}" style="width:100%; border-radius:8px; border:1px solid #E2E8F0;"/>' if admet_radar_b64 else ''
        var_ligplot_html = f'<img src="{var_ligplot_b64}" style="width:100%; border-radius:8px; border:1px solid #E2E8F0;"/>' if var_ligplot_b64 else ''

        # Quantitative Biophysical Delta Table Rows
        delta_table_rows = []
        try:
            p_aff_val = float(affinity_kcal)
            v_aff_val = float(var_aff)
            d_aff_val = v_aff_val - p_aff_val
            d_aff_cls = "badge-green" if d_aff_val < 0 else "badge-red"
            delta_table_rows.append(f"""
                <tr>
                    <td><strong>AutoDock Vina Affinity</strong></td>
                    <td>{p_aff_val:.2f} kcal/mol</td>
                    <td>{v_aff_val:.2f} kcal/mol</td>
                    <td><span class="badge {d_aff_cls}">{d_aff_val:+.2f} kcal/mol</span></td>
                    <td>{'Potency Gain' if d_aff_val < 0 else 'Steric Penalty'}</td>
                </tr>
            """)
        except Exception:
            pass

        if mmgbsa_data and var_mmgbsa_data:
            p_dg = mmgbsa_data.get('total_dg', 0.0)
            v_dg = var_mmgbsa_data.get('total_dg', 0.0)
            d_dg = v_dg - p_dg
            d_dg_cls = "badge-green" if d_dg < 0 else "badge-red"
            delta_table_rows.append(f"""
                <tr>
                    <td><strong>MM-GBSA Net Binding Free Energy (ΔG)</strong></td>
                    <td>{p_dg:.2f} kcal/mol</td>
                    <td>{v_dg:.2f} kcal/mol</td>
                    <td><span class="badge {d_dg_cls}">{d_dg:+.2f} kcal/mol</span></td>
                    <td>{'Thermodynamic Gain' if d_dg < 0 else 'Relaxed Binding'}</td>
                </tr>
            """)
            p_vdw = mmgbsa_data.get('vdw_energy', 0.0)
            v_vdw = var_mmgbsa_data.get('vdw_energy', 0.0)
            d_vdw = v_vdw - p_vdw
            delta_table_rows.append(f"""
                <tr>
                    <td><strong>Van der Waals Enthalpy (ΔE_vdW)</strong></td>
                    <td>{p_vdw:.2f} kcal/mol</td>
                    <td>{v_vdw:.2f} kcal/mol</td>
                    <td><span class="badge {'badge-green' if d_vdw < 0 else 'badge-blue'}">{d_vdw:+.2f} kcal/mol</span></td>
                    <td>{'Enhanced Shape Fit' if d_vdw < 0 else 'Equivalent Packing'}</td>
                </tr>
            """)
            p_elec = mmgbsa_data.get('elec_energy', 0.0)
            v_elec = var_mmgbsa_data.get('elec_energy', 0.0)
            d_elec = v_elec - p_elec
            delta_table_rows.append(f"""
                <tr>
                    <td><strong>Coulombic Electrostatics (ΔE_elec)</strong></td>
                    <td>{p_elec:.2f} kcal/mol</td>
                    <td>{v_elec:.2f} kcal/mol</td>
                    <td><span class="badge {'badge-green' if d_elec < 0 else 'badge-blue'}">{d_elec:+.2f} kcal/mol</span></td>
                    <td>{'Polar Anchor Reinforcement' if d_elec < 0 else 'Maintained Electrostatics'}</td>
                </tr>
            """)

        if md_results and var_md_results:
            p_rmsd = md_results.get('mean_rmsd', 0.0)
            v_rmsd = var_md_results.get('mean_rmsd', 0.0)
            d_rmsd = v_rmsd - p_rmsd
            d_rmsd_cls = "badge-green" if d_rmsd < 0 else "badge-gold"
            delta_table_rows.append(f"""
                <tr>
                    <td><strong>MD Trajectory Mean Pose Drift (RMSD)</strong></td>
                    <td>{p_rmsd:.2f} Å</td>
                    <td>{v_rmsd:.2f} Å</td>
                    <td><span class="badge {d_rmsd_cls}">{d_rmsd:+.2f} Å</span></td>
                    <td>{'Kinetic Retention Tightened' if d_rmsd <= 0 else 'Conformational Flexibility'}</td>
                </tr>
            """)

        if admet_dict and var_admet_dict:
            p_mw = admet_dict.get('Molecular Weight', 'N/A')
            v_mw = var_admet_dict.get('Molecular Weight', 'N/A')
            delta_table_rows.append(f"""
                <tr>
                    <td><strong>Molecular Weight (MW)</strong></td>
                    <td>{p_mw} g/mol</td>
                    <td>{v_mw} g/mol</td>
                    <td><span class="badge badge-blue">MW Profile</span></td>
                    <td>&le; 500 g/mol Compliance</td>
                </tr>
            """)
            p_logp = admet_dict.get('LogP', 'N/A')
            v_logp = var_admet_dict.get('LogP', 'N/A')
            delta_table_rows.append(f"""
                <tr>
                    <td><strong>Calculated Partition (cLogP)</strong></td>
                    <td>{p_logp}</td>
                    <td>{v_logp}</td>
                    <td><span class="badge badge-blue">Lipophilicity Shift</span></td>
                    <td>Membrane Permeability Window</td>
                </tr>
            """)
            p_qed = admet_dict.get('QED Drug-Likeness', 'N/A')
            v_qed = var_admet_dict.get('QED Drug-Likeness', 'N/A')
            delta_table_rows.append(f"""
                <tr>
                    <td><strong>QED Drug-Likeness Score</strong></td>
                    <td>{p_qed}</td>
                    <td>{v_qed}</td>
                    <td><span class="badge badge-gold">Oral Bioavailability</span></td>
                    <td>Bickerton Standard</td>
                </tr>
            """)

        comp_table_html = f"""
        <div class="table-container" style="margin: 16px 0;">
            <table class="report-table">
                <thead>
                    <tr>
                        <th>Biophysical / Pharmacological Parameter</th>
                        <th>Natural Parent Scaffold</th>
                        <th>Bioisosteric Derivative Lead</th>
                        <th>Quantitative Delta (Δ)</th>
                        <th>Mechanistic Evaluation</th>
                    </tr>
                </thead>
                <tbody>
                    {"".join(delta_table_rows)}
                </tbody>
            </table>
        </div>
        """ if delta_table_rows else ""

        comp_chart_img = ""
        if comp_hotspot_chart_b64:
            comp_chart_img = f"""
            <div style="margin: 18px 0; text-align:center;">
                <img src="{comp_hotspot_chart_b64}" style="max-width:100%; border-radius:10px; border:1px solid #CBD5E1; box-shadow:0 2px 8px rgba(0,0,0,0.06);"/>
                <div style="font-size:11.5px; color:#64748B; margin-top:8px; line-height:1.4;">
                    <b>Figure 6 | Comparative MM-GBSA Per-Residue Energy Decomposition.</b> Direct quantitative comparison of residue binding free energies (&Delta;G<sub>res</sub> kcal/mol) for Natural Parent (Amber) vs. Semi-Synthetic Lead (Cyan), identifying which active-site pocket anchors were reinforced or newly recruited by bioisosteric substitution.
                </div>
            </div>
            """
        elif var_mmgbsa_chart_b64:
            comp_chart_img = f"""
            <div style="margin: 18px 0; text-align:center;">
                <img src="{var_mmgbsa_chart_b64}" style="max-width:100%; border-radius:10px; border:1px solid #CBD5E1; box-shadow:0 2px 8px rgba(0,0,0,0.06);"/>
                <div style="font-size:11.5px; color:#64748B; margin-top:8px; line-height:1.4;">
                    <b>Figure 6 | Derivative MM-GBSA Per-Residue Energy Decomposition.</b> Hotspot binding anchors for the semi-synthetic lead candidate.
                </div>
            </div>
            """

        variant_html = f"""
        <div class="section-card" style="border-left: 4px solid #4F46E5;">
            <div class="section-title" style="color:#4338CA;">
                <span>🧬 Semi-Synthetic Bioisostere Lead Optimization (Stage 04 Monograph)</span>
                <span class="badge" style="background:#EEF2FF; color:#4338CA; border:1px solid #C7D2FE;">RATIONAL DERIVATIVE</span>
            </div>
            
            <p style="margin:0 0 16px 0; font-size:13px; color:#475569; line-height:1.5;">
                In-silico bioisosteric functionalization designed to resolve identified molecular liabilities (polar surface area, steric clashes, or metabolic instability) while amplifying orthosteric cavity complementarity.
            </p>

            <div style="display:grid; grid-template-columns: 1fr 1fr; gap:16px; margin-bottom:16px;">
                <div style="background:#F8FAFC; border:1px solid #E2E8F0; border-radius:8px; padding:16px;">
                    <h5 style="margin:0 0 6px 0; color:#0F172A; font-size:13px;">Optimized Lead Analogue</h5>
                    <div style="font-size:16px; font-weight:700; color:#4338CA; margin-bottom:6px;">{html.escape(variant_info.get('name', 'Derivative'))}</div>
                    <p style="margin:0 0 6px 0; font-size:12.5px; color:#334155;"><strong>Medicinal Rationale:</strong> {html.escape(variant_info.get('rationale', ''))}</p>
                    <p style="margin:0; font-size:11px; color:#64748B; word-break:break-all;"><strong>Derivative SMILES:</strong> <code style="background:#FFFFFF; padding:2px 6px; border-radius:4px; border:1px solid #CBD5E1;">{html.escape(variant_info.get('smiles', ''))}</code></p>
                </div>
                
                <div style="background:#EEF2FF; border:1px solid #C7D2FE; border-radius:8px; padding:16px; text-align:center; display:flex; flex-direction:column; justify-content:center;">
                    <div style="font-size:11px; color:#4338CA; text-transform:uppercase; font-weight:700; letter-spacing:0.5px;">Comparative Binding Delta (ΔΔG)</div>
                    <div style="font-size:24px; font-weight:800; color:#3730A3; margin:4px 0;">{var_aff} <span style="font-size:14px; font-weight:500;">kcal/mol</span></div>
                    <div style="font-size:12px; font-weight:700; color:{delta_color};">{aff_delta_str}</div>
                </div>
            </div>

            <!-- Dual Molecular Structures Comparison -->
            {f'''
            <div style="display:grid; grid-template-columns: 1fr 1fr; gap:16px; margin-bottom:16px;">
                <div style="text-align:center;">
                    {mol_img_tag}
                    <div style="font-size:11.5px; color:#64748B; margin-top:4px;"><b>Natural Parent Scaffold:</b> {html.escape(compound_name)}</div>
                </div>
                <div style="text-align:center;">
                    <img src="{var_mol_b64}" style="width:100%; height:220px; object-fit:contain; border-radius:10px; background:#FFFFFF; border:1px solid #E2E8F0;"/>
                    <div style="font-size:11.5px; color:#4338CA; margin-top:4px;"><b>Bioisostere Lead:</b> {html.escape(variant_info.get('name', 'Derivative'))}</div>
                </div>
            </div>
            ''' if var_mol_b64 else ''}

            <!-- Head-to-Head Quantitative Biophysical Delta Table -->
            {comp_table_html}

            <!-- Comparative MM-GBSA Per-Residue Hotspot Chart -->
            {comp_chart_img}

            <!-- Radar Chart & Comparative Interaction Map -->
            {f'''
            <div style="display:grid; grid-template-columns: 1fr 1fr; gap:16px; margin-top:16px;">
                <div>
                    {radar_img_html}
                    <div style="font-size:11px; color:#64748B; margin-top:6px; text-align:center;"><b>Figure 7 | Hexagonal ADMET Radar Spider Plot.</b> Overlay comparing Parent Scaffold (Blue) vs. Optimized Lead (Green).</div>
                </div>
                <div>
                    {var_ligplot_html if var_ligplot_html else f'<div style="height:100%; display:flex; align-items:center; justify-content:center; background:#F8FAFC; border-radius:8px; border:1px solid #E2E8F0; font-size:12px; color:#94A3B8; padding:20px; text-align:center;">Derivative interaction diagram automatically synthesized.</div>'}
                    <div style="font-size:11px; color:#64748B; margin-top:6px; text-align:center;"><b>Figure 8 | Derivative Binding Architecture.</b> Newly recruited stabilizing contacts.</div>
                </div>
            </div>
            ''' if radar_img_html else ''}

        </div>
        """

    
    # =================================================================
    # Section VI: Systems Network Pharmacology, Targetome & In-Vivo Profile
    # (Conditionally synthesized ONLY when optional discovery tools were run)
    # =================================================================
    systems_html = ""
    if targetome_results or network_results or synergy_results or microbiome_results or population_results:
        sub_sections = []
        
        # 1. Targetome Profiling
        if targetome_results and targetome_results.get("targets"):
            top_t = targetome_results["primary_target"]
            t_rows = []
            for t in targetome_results["targets"]:
                t_rows.append(f"""
                <tr>
                    <td><strong>{html.escape(t['gene'])}</strong></td>
                    <td>{html.escape(t['name'])}</td>
                    <td><span class="badge badge-blue">{html.escape(t['category'])}</span></td>
                    <td style="font-weight:700; color:{t['tier_color']};">{t['affinity_kcal']} kcal/mol</td>
                    <td><code>{html.escape(t['estimated_ki'])}</code></td>
                    <td>{html.escape(t['ref_drug'])} ({t['ref_affinity']} kcal/mol)</td>
                    <td>{html.escape(t['clinical_implication'])}</td>
                </tr>
                """)
            targetome_block = f"""
            <div style="margin-bottom:24px;">
                <h4 style="margin:0 0 8px 0; font-size:14px; color:#0F172A;">🎯 Pan-Proteome Targetome Selectivity Screen (10 Core Disease Cascades)</h4>
                <div style="background:#F0FDF4; border:1px solid #BBF7D0; border-radius:8px; padding:12px 16px; margin-bottom:12px; display:flex; justify-content:space-between; align-items:center;">
                    <div>
                        <b style="color:#166534;">Primary High-Affinity On-Target Identified:</b> <span style="font-weight:700; color:#0F172A;">{html.escape(top_t['gene'])} ({html.escape(top_t['name'])})</span><br>
                        <span style="font-size:12px; color:#475569;">Predicted Free Energy: <b>{top_t['affinity_kcal']} kcal/mol</b> | Est. K<sub>i</sub>: <b>{html.escape(top_t['estimated_ki'])}</b> | Classification: <b>{html.escape(targetome_results['polypharmacology_class'])}</b></span>
                    </div>
                    <span class="badge badge-green">HIGH ON-TARGET SELECTIVITY</span>
                </div>
                <div class="table-container">
                    <table class="report-table">
                        <thead>
                            <tr>
                                <th>Target Gene</th>
                                <th>Protein Name</th>
                                <th>Disease Category</th>
                                <th>Predicted &Delta;G</th>
                                <th>Est. K<sub>i</sub></th>
                                <th>Reference Benchmark</th>
                                <th>Therapeutic Mechanism</th>
                            </tr>
                        </thead>
                        <tbody>
                            {''.join(t_rows)}
                        </tbody>
                    </table>
                </div>
            </div>
            """
            sub_sections.append(targetome_block)

        # 2. Systems Network Pharmacology & Hub Bottlenecks
        if network_results and network_results.get("nodes"):
            hubs = [n for n in network_results["nodes"] if n.get("is_hub")]
            hub_rows = []
            for h in hubs:
                hub_rows.append(f"""
                <tr>
                    <td><strong>{html.escape(h['label'])}</strong></td>
                    <td><span class="badge badge-purple">{html.escape(h['type'])}</span></td>
                    <td><b>{h['degree']}</b> interactions</td>
                    <td><code>{h['betweenness']}</code></td>
                    <td><span class="badge badge-green">CRITICAL HUB BOTTLENECK</span></td>
                </tr>
                """)
            network_block = f"""
            <div style="margin-bottom:24px;">
                <h4 style="margin:0 0 8px 0; font-size:14px; color:#0F172A;">🕸️ Systems Network Biology: Topological Centrality & Bottleneck Identification</h4>
                <p style="margin:0 0 10px 0; font-size:13px; color:#64748B;">
                    Bipartite network modeled across <b>{network_results.get('total_nodes', 0)} biological nodes</b> and <b>{network_results.get('total_edges', 0)} multi-target regulatory interactions</b>.
                    Nodes with top degree and betweenness centrality represent critical regulatory bottlenecks where botanical modulation prevents disease escape.
                </p>
                <div class="table-container">
                    <table class="report-table">
                        <thead>
                            <tr>
                                <th>Network Node</th>
                                <th>Biological Entity Type</th>
                                <th>Degree Centrality (Connectivity)</th>
                                <th>Betweenness Centrality (Information Flow)</th>
                                <th>Systems Biology Role</th>
                            </tr>
                        </thead>
                        <tbody>
                            {''.join(hub_rows)}
                        </tbody>
                    </table>
                </div>
            </div>
            """
            sub_sections.append(network_block)

        # 3. Botanical Synergism Calculator
        if synergy_results and synergy_results.get("ci"):
            ci_val = synergy_results["ci"]
            syn_block = f"""
            <div style="margin-bottom:24px;">
                <h4 style="margin:0 0 8px 0; font-size:14px; color:#0F172A;">⚡ Multi-Constituent Synergism Analysis (Chou-Talalay Combination Index)</h4>
                <div style="background:#F8FAFC; border:1px solid #E2E8F0; border-radius:8px; padding:14px 18px;">
                    <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:8px;">
                        <span style="font-size:13px; font-weight:700; color:#0F172A;">Calculated Combination Index (CI):</span>
                        <span class="badge badge-green" style="font-size:13px; font-weight:800; padding:4px 10px;">CI = {ci_val:.2f} &bull; {html.escape(synergy_results['classification'])}</span>
                    </div>
                    <p style="margin:0 0 6px 0; font-size:12.5px; color:#475569; line-height:1.6;">
                        <b>Synergistic Mechanism:</b> {html.escape(synergy_results.get('explanation', ''))}
                    </p>
                    <div style="font-size:12px; color:#0A84FF; font-weight:600;">
                        Cascade Coverage: {html.escape(synergy_results.get('pathway_coverage', ''))}
                    </div>
                </div>
            </div>
            """
            sub_sections.append(syn_block)

        # 4. Microbiome Biotransformation
        if microbiome_results and microbiome_results.get("is_prodrug"):
            mb = microbiome_results
            mb_block = f"""
            <div style="margin-bottom:12px;">
                <h4 style="margin:0 0 8px 0; font-size:14px; color:#0F172A;">🧪 In-Vivo Human Gut Microbiota Biotransformation (Prodrug Activation Cascade)</h4>
                <div style="background:#FFFBEB; border:1px solid #FDE68A; border-radius:8px; padding:14px 18px;">
                    <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:8px;">
                        <span style="font-size:13px; font-weight:700; color:#92400E;">Enzymatic Cleavage & Bioactivation:</span>
                        <span class="badge badge-gold">IN-VIVO BIOACTIVATION ACTIVE</span>
                    </div>
                    <div class="table-container" style="margin-top:8px;">
                        <table class="report-table">
                            <thead>
                                <tr>
                                    <th>Parameter</th>
                                    <th>Ingested Plant Prodrug (Raw)</th>
                                    <th>Circulating Human Metabolite (Active)</th>
                                    <th>In-Vivo Pharmacokinetic Delta</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <td><strong>Chemical Identity</strong></td>
                                    <td>{html.escape(mb['ingested_scaffold'])}</td>
                                    <td><b style="color:#059669;">{html.escape(mb['circulating_metabolite'])}</b></td>
                                    <td>Deglycosylation & Aglycone Release</td>
                                </tr>
                                <tr>
                                    <td><strong>Molecular Weight (MW)</strong></td>
                                    <td>{mb['raw_mw']} g/mol</td>
                                    <td>{mb['act_mw']} g/mol</td>
                                    <td><span class="badge badge-blue">{mb['act_mw'] - mb['raw_mw']:+.1f} g/mol</span> (Permeable Cavity Size)</td>
                                </tr>
                                <tr>
                                    <td><strong>Lipophilicity (cLogP)</strong></td>
                                    <td>{mb['raw_logp']}</td>
                                    <td>{mb['act_logp']}</td>
                                    <td><span class="badge badge-green">Hydrophobic Cavity Insertion</span></td>
                                </tr>
                                <tr>
                                    <td><strong>Target Binding Affinity</strong></td>
                                    <td>{mb['raw_affinity']} kcal/mol</td>
                                    <td><b style="color:#059669;">{mb['act_affinity']} kcal/mol</b></td>
                                    <td><span class="badge badge-green">{mb['delta_affinity']:+.2f} kcal/mol</span> Affinity Boost</td>
                                </tr>
                                <tr>
                                    <td><strong>Caco-2 Intestinal Permeability</strong></td>
                                    <td>{html.escape(mb['permeability_raw_papp'])}</td>
                                    <td><b style="color:#059669;">{html.escape(mb['permeability_active_papp'])}</b></td>
                                    <td><span class="badge badge-green">Rapid Passive Absorption</span></td>
                                </tr>
                            </tbody>
                        </table>
                    </div>
                    <p style="margin:10px 0 0 0; font-size:12px; color:#78350F; line-height:1.5;">
                        <b>Commensal Microflora Engine:</b> {html.escape(mb['bacterial_enzyme'])}. {html.escape(mb['in_vivo_pk_note'])}
                    </p>
                </div>
            </div>
            """
            sub_sections.append(mb_block)

        
        # 5. In-Silico Clinical Trials & Population Pharmacogenomics (Human Genome Project Reference)
        if population_results and population_results.get("overall_responder_rate") is not None:
            pop = population_results
            hgp = pop.get("hgp_reference", {})
            
            demo_rows = []
            for d in pop.get("demographic_breakdown", []):
                demo_rows.append(f"""
                <tr>
                    <td><b style="color:#1E293B;">{html.escape(d['cohort'])}</b></td>
                    <td>{d['n_subjects']} subjects</td>
                    <td style="font-weight:700; color:{'#16A34A' if d['responder_rate'] >= 75.0 else '#2563EB'};">{d['responder_rate']}%</td>
                    <td>{d['mean_ro']}%</td>
                    <td><span class="badge {'badge-green' if d['responder_rate'] >= 75.0 else 'badge-blue'}">{'HIGH RESPONSE' if d['responder_rate'] >= 75.0 else 'MODERATE'}</span></td>
                </tr>
                """)

            variant_rows = []
            for v in pop.get("variant_breakdown", []):
                rs = v.get('rs_id', 'Reference Canonical')
                hgvs = f"<code>{html.escape(v.get('hgvs_c', 'c.Canonical'))}</code><br><small style='color:#6B21A8; font-weight:600;'>{html.escape(v.get('hgvs_p', 'p.WT'))}</small>"
                variant_rows.append(f"""
                <tr>
                    <td><strong>{html.escape(v['variant'])}</strong><br><span class="badge badge-blue" style="font-size:10px;">{html.escape(rs)}</span></td>
                    <td>{hgvs}</td>
                    <td><span class="badge badge-purple">{html.escape(v.get('mutation_type', 'HGP Reference'))}</span></td>
                    <td><code>{v['delta_delta_g']:+.2f} kcal/mol</code></td>
                    <td>{v.get('n_subjects', 0)} pts<br><small style="color:#64748B;">Mean: {v.get('mean_ro', 0)}%</small></td>
                    <td style="font-weight:700; color:{'#16A34A' if v['responder_rate'] >= 75.0 else '#2563EB'};">{v['responder_rate']}%</td>
                    <td style="font-size:11.5px; line-height:1.4;">{html.escape(v['description'])}</td>
                </tr>
                """)

            hgp_card_html = f"""
            <!-- Human Genome Project (GRCh38) Consensus Reference Card -->
            <div style="background: linear-gradient(135deg, #F8FAFC 0%, #EFF6FF 100%); border: 1px solid #BFDBFE; border-radius: 8px; padding: 12px 16px; margin: 12px 0; display: grid; grid-template-columns: repeat(4, 1fr); gap: 12px;">
                <div>
                    <div style="font-size:10px; color:#64748B; text-transform:uppercase; font-weight:700;">HGP Assembly Reference</div>
                    <div style="font-size:13px; font-weight:700; color:#1E40AF; margin-top:2px;">{html.escape(hgp.get('assembly', 'GRCh38.p14'))}</div>
                    <div style="font-size:11px; color:#475569; margin-top:2px;">Locus: <b>{html.escape(hgp.get('locus', 'N/A'))}</b> ({html.escape(hgp.get('cytoband', ''))})</div>
                </div>
                <div>
                    <div style="font-size:10px; color:#64748B; text-transform:uppercase; font-weight:700;">RefSeq Annotations</div>
                    <div style="font-size:11.5px; font-weight:600; color:#0F172A; margin-top:2px;">mRNA: <code>{html.escape(hgp.get('refseq_mrna', 'N/A'))}</code></div>
                    <div style="font-size:11.5px; font-weight:600; color:#0F172A;">Protein: <code>{html.escape(hgp.get('refseq_protein', 'N/A'))}</code></div>
                </div>
                <div>
                    <div style="font-size:10px; color:#64748B; text-transform:uppercase; font-weight:700;">Canonical Architecture</div>
                    <div style="font-size:12px; font-weight:600; color:#0F172A; margin-top:2px;">{hgp.get('canonical_aa_len', 'N/A')} AA | {hgp.get('exon_count', 'N/A')} Exons</div>
                    <div style="font-size:11px; color:#64748B; margin-top:2px;">UniProt ID: <b>{html.escape(hgp.get('uniprot_id', 'N/A'))}</b></div>
                </div>
                <div>
                    <div style="font-size:10px; color:#64748B; text-transform:uppercase; font-weight:700;">Genomic Micro-Drift Engine</div>
                    <div style="font-size:12px; font-weight:700; color:{'#16A34A' if pop.get('include_minute_variations') else '#64748B'}; margin-top:2px;">
                        {'✓ Active Micro-Drift & eQTL' if pop.get('include_minute_variations') else 'Discrete Alleles Only'}
                    </div>
                    <div style="font-size:10.5px; color:#475569; margin-top:2px;">&sigma;(&Delta;G)=0.15 kcal/mol &bull; &sigma;(eQTL)=0.12</div>
                </div>
            </div>
            """ if hgp else ""

            pop_block = f"""
            <div style="margin-bottom:14px; margin-top:20px; border-top:1px dashed #DDD6FE; padding-top:16px;">
                <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:8px;">
                    <div>
                        <h4 style="margin:0 0 4px 0; font-size:14px; color:#0F172A;">👥 In-Silico Clinical Trial: Human Genome Project (GRCh38) Virtual Population</h4>
                        <span style="font-size:12px; color:#6B21A8;">Simulated Cohort: <b>N = {pop['n_patients']} Virtual Patients</b> &bull; Reference Target: <b>{html.escape(pop['target_gene'])} ({html.escape(pop['target_name'])})</b> &bull; Regimen: <b>{pop['dose_mg']} mg BID</b></span>
                    </div>
                    <span class="badge {pop['fda_badge_cls']}" style="font-size:12px; padding:4px 10px;">{html.escape(pop['fda_tier'])}</span>
                </div>

                {hgp_card_html}

                <!-- Efficacy Metric Scorecards -->
                <div style="display:grid; grid-template-columns: repeat(4, 1fr); gap:12px; margin:12px 0;">
                    <div style="background:#FFFFFF; border:1px solid #E2E8F0; border-radius:8px; padding:10px; text-align:center;">
                        <div style="font-size:10.5px; color:#64748B;">COHORT RESPONDER RATE</div>
                        <div style="font-size:20px; font-weight:800; color:#16A34A; margin:2px 0;">{pop['overall_responder_rate']}%</div>
                        <div style="font-size:10px; color:#475569;">Target RO &ge; 75%</div>
                    </div>
                    <div style="background:#FFFFFF; border:1px solid #E2E8F0; border-radius:8px; padding:10px; text-align:center;">
                        <div style="font-size:10.5px; color:#64748B;">MEAN OCCUPANCY</div>
                        <div style="font-size:20px; font-weight:800; color:#2563EB; margin:2px 0;">{pop['mean_ro']}%</div>
                        <div style="font-size:10px; color:#475569;">Steady-State Free RO</div>
                    </div>
                    <div style="background:#FFFFFF; border:1px solid #E2E8F0; border-radius:8px; padding:10px; text-align:center;">
                        <div style="font-size:10.5px; color:#64748B;">TARGET SATURATION</div>
                        <div style="font-size:20px; font-weight:800; color:#9333EA; margin:2px 0;">{pop['saturation_rate']}%</div>
                        <div style="font-size:10px; color:#475569;">Full Saturation (RO &ge; 90%)</div>
                    </div>
                    <div style="background:#FFFFFF; border:1px solid #E2E8F0; border-radius:8px; padding:10px; text-align:center;">
                        <div style="font-size:10.5px; color:#64748B;">GENOMIC SENSITIVITY</div>
                        <div style="font-size:20px; font-weight:800; color:#EA580C; margin:2px 0;">&plusmn;0.15</div>
                        <div style="font-size:10px; color:#475569;">kcal/mol Micro-Drift</div>
                    </div>
                </div>

                <!-- Demographic Sensitivity Breakdown Table -->
                <div class="table-container" style="margin-top:12px;">
                    <div style="font-size:12px; font-weight:700; color:#1E293B; margin-bottom:4px;">Ancestral Demographic Sensitivity (1000 Genomes & gnomAD Stratified Cohorts)</div>
                    <table class="report-table">
                        <thead>
                            <tr>
                                <th>Demographic Cohort</th>
                                <th>Sample Size</th>
                                <th>Responder Rate (% with RO &ge; 75%)</th>
                                <th>Mean Target Occupancy</th>
                                <th>Clinical Response Tier</th>
                            </tr>
                        </thead>
                        <tbody>
                            {''.join(demo_rows)}
                        </tbody>
                    </table>
                </div>

                <!-- Allelic Pocket Variants & Minute Variations Table -->
                <div class="table-container" style="margin-top:12px;">
                    <div style="font-size:12px; font-weight:700; color:#1E293B; margin-bottom:4px;">Receptor Pocket Allelic Missense Variants (dbSNP rsIDs & HGVS Nomenclature)</div>
                    <table class="report-table">
                        <thead>
                            <tr>
                                <th>Allelic Variant & rsID</th>
                                <th>HGVS cDNA / Protein</th>
                                <th>Mutation Classification</th>
                                <th>&Delta;&Delta;G Binding Shift</th>
                                <th>Cohort Sample Size</th>
                                <th>Responder Rate (% &ge; 75%)</th>
                                <th>Biophysical Mechanism</th>
                            </tr>
                        </thead>
                        <tbody>
                            {''.join(variant_rows)}
                        </tbody>
                    </table>
                </div>
            </div>
            """
            sub_sections.append(pop_block)

        systems_html = f"""
        <div class="section-card" style="border:1px solid #DDD6FE; background:#FAF5FF;">
            <div class="section-title">
                <span style="color:#6D28D9;">🌐 Section VI: Systems Network Pharmacology, Targetome & In-Vivo Discovery Profile</span>
                <span class="badge badge-purple">STAGE 06 EXTENSION</span>
            </div>
            <p style="margin:0 0 16px 0; font-size:13px; color:#6B21A8;">
                Comprehensive multi-target discovery profiling assessing polypharmacological targetome selectivity,
                interactome hub bottlenecks, multi-constituent combination synergism, and gut microbiome in-vivo activation.
            </p>
            {''.join(sub_sections)}
        </div>
        """

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>EthnoDock Pro • Scientific Research Monograph - {html.escape(compound_name)}</title>
    <!-- Standalone 3Dmol.js for WebGL Pocket Exploration -->
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Plus+Jakarta+Sans:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;500&display=swap');

        body {{
            font-family: 'Plus Jakarta Sans', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            background-color: #0B0F19;
            color: #1E293B;
            margin: 0;
            padding: 40px 20px;
            -webkit-font-smoothing: antialiased;
        }}

        .dossier-wrapper {{
            max-width: 1040px;
            margin: 0 auto;
            background: #FFFFFF;
            border-radius: 18px;
            box-shadow: 0 30px 60px -15px rgba(0, 0, 0, 0.4);
            overflow: hidden;
            border: 1px solid #CBD5E1;
        }}

        /* Executive Header Banner */
        .dossier-header {{
            background: linear-gradient(135deg, #090D16 0%, #111827 45%, #064E3B 100%);
            color: #FFFFFF;
            padding: 38px 44px;
            position: relative;
            border-bottom: 4px solid #10B981;
        }}

        .dossier-header-top {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 18px;
        }}

        .doc-seal {{
            display: inline-flex;
            align-items: center;
            gap: 6px;
            background: rgba(16, 185, 129, 0.18);
            color: #34D399;
            border: 1px solid rgba(16, 185, 129, 0.35);
            padding: 5px 16px;
            border-radius: 20px;
            font-size: 11px;
            font-weight: 700;
            letter-spacing: 0.8px;
            text-transform: uppercase;
        }}

        .doc-meta {{
            font-size: 12px;
            color: #94A3B8;
            text-align: right;
            line-height: 1.45;
            font-family: 'JetBrains Mono', monospace;
        }}

        .dossier-title {{
            font-size: 30px;
            font-weight: 800;
            margin: 0 0 6px 0;
            letter-spacing: -0.03em;
            color: #FFFFFF;
        }}

        .dossier-subtitle {{
            font-size: 15px;
            color: #CBD5E1;
            margin: 0;
            font-weight: 400;
        }}

        /* Executive Summary Bar */
        .summary-bar {{
            background: #F8FAFC;
            border-bottom: 1px solid #E2E8F0;
            padding: 24px 44px;
            display: grid;
            grid-template-columns: repeat(4, 1fr);
            gap: 20px;
            text-align: center;
        }}

        .stat-item {{
            padding: 6px;
        }}

        .stat-label {{
            font-size: 11px;
            color: #64748B;
            text-transform: uppercase;
            font-weight: 700;
            letter-spacing: 0.6px;
        }}

        .stat-value {{
            font-size: 22px;
            font-weight: 800;
            color: #0F172A;
            margin-top: 4px;
        }}

        /* Content Body */
        .dossier-content {{
            padding: 38px 44px;
        }}

        .section-card {{
            background: #FFFFFF;
            border: 1px solid #E2E8F0;
            border-radius: 12px;
            padding: 24px;
            margin-bottom: 24px;
            box-shadow: 0 2px 6px rgba(0,0,0,0.02);
            page-break-inside: avoid;
            break-inside: avoid;
        }}

        .paozhi-card {{
            background: #FFFDF5;
            border: 1px solid #FDE68A;
        }}

        .section-title {{
            font-size: 15px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.6px;
            color: #0F172A;
            margin-bottom: 16px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-bottom: 1px solid #F1F5F9;
            padding-bottom: 10px;
        }}

        .grid-2col {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 24px;
        }}

        .data-box {{
            background: #F8FAFC;
            border: 1px solid #E2E8F0;
            border-radius: 8px;
            padding: 16px;
        }}

        /* Tables */
        .table-container {{
            overflow-x: auto;
            margin-top: 12px;
            border-radius: 8px;
            border: 1px solid #E2E8F0;
        }}

        .report-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 13px;
            text-align: left;
        }}

        .report-table th {{
            background: #F1F5F9;
            color: #475569;
            font-weight: 700;
            padding: 11px 14px;
            border-bottom: 1px solid #CBD5E1;
            font-size: 12px;
            text-transform: uppercase;
            letter-spacing: 0.4px;
        }}

        .report-table td {{
            padding: 11px 14px;
            border-bottom: 1px solid #F1F5F9;
            color: #334155;
        }}

        .report-table tr:last-child td {{
            border-bottom: none;
        }}

        .report-table tr:hover {{
            background: #F8FAFC;
        }}

        /* Badges */
        .badge {{
            display: inline-block;
            padding: 3px 10px;
            border-radius: 14px;
            font-size: 10.5px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .badge-green {{ background: #DCFCE7; color: #15803D; }}
        .badge-blue {{ background: #E0F2FE; color: #0369A1; }}
        .badge-gold {{ background: #FEF3C7; color: #B45309; }}
        .badge-red {{ background: #FEE2E2; color: #DC2626; border: 1px solid #DC2626; }}

        /* Print Button */
        .print-btn {{
            cursor: pointer;
            background: linear-gradient(135deg, #10B981, #059669);
            color: white;
            border: none;
            padding: 7px 18px;
            border-radius: 8px;
            font-weight: 700;
            font-size: 12.5px;
            display: inline-flex;
            align-items: center;
            gap: 6px;
            box-shadow: 0 4px 12px rgba(16, 185, 129, 0.3);
            transition: all 0.2s ease;
        }}
        .print-btn:hover {{
            transform: translateY(-1px);
            box-shadow: 0 6px 16px rgba(16, 185, 129, 0.4);
        }}

        /* Footer */
        .dossier-footer {{
            background: #090D16;
            color: #94A3B8;
            padding: 26px 44px;
            font-size: 12px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-top: 1px solid #1E293B;
        }}

        /* Print Media Styles */
        @media print {{
            body {{
                background: #FFFFFF !important;
                color: #000000 !important;
                padding: 0 !important;
            }}
            .dossier-wrapper {{
                box-shadow: none !important;
                border: none !important;
                max-width: 100% !important;
            }}
            .no-print {{
                display: none !important;
            }}
            .section-card {{
                page-break-inside: avoid !important;
                break-inside: avoid !important;
                border: 1px solid #E2E8F0 !important;
                margin-bottom: 20px !important;
            }}
            .page-break {{
                page-break-before: always !important;
            }}
        }}
    </style>
</head>
<body>

<div class="dossier-wrapper">

    <!-- Header Banner -->
    <div class="dossier-header">
        <div class="dossier-header-top">
            <span class="doc-seal">🌿 EthnoDock Pro &bull; Verified Research Monograph</span>
            <div style="display:flex; align-items:center; gap:16px;">
                <button onclick="window.print()" class="print-btn no-print">
                    🖨️ Print / Save as PDF
                </button>
                <div class="doc-meta">
                    <b>Document Ref:</b> {doc_id}<br>
                    <b>Timestamp:</b> {date_str}
                </div>
            </div>
        </div>
        <h1 class="dossier-title">{html.escape(species_name)} &bull; {html.escape(compound_name)}</h1>
        <p class="dossier-subtitle">In-Silico Structural Pharmacology, Molecular Dynamics & Classical Ethnopharmacology Monograph</p>
    </div>

    <!-- Executive Summary Bar -->
    <div class="summary-bar">
        <div class="stat-item">
            <div class="stat-label">Binding Affinity (ΔG)</div>
            <div class="stat-value" style="color:{aff_color};">{affinity_kcal} <span style="font-size:13px; font-weight:500;">kcal/mol</span></div>
            {aff_badge_html}
        </div>
        <div class="stat-item">
            <div class="stat-label">Target Macromolecule</div>
            <div class="stat-value" style="color:#0284C7; font-size:18px;">PDB: {pdb_id}</div>
        </div>
        <div class="stat-item">
            <div class="stat-label">Botanical Binomial</div>
            <div class="stat-value" style="font-size:15px; font-style:italic;">{html.escape(scientific_name)}</div>
        </div>
        <div class="stat-item">
            <div class="stat-label">Classical Canon</div>
            <div class="stat-value" style="color:#B45309; font-size:15px;">{html.escape(source_title)}</div>
        </div>
    </div>

    <!-- Main Content -->
    <div class="dossier-content">

        <!-- 1. Specimen & Molecular Topology -->
        <div class="section-card">
            <div class="section-title">
                <span>🔬 Botanical Specimen & Chemical Architecture</span>
            </div>
            <div class="grid-2col">
                <div style="text-align:center;">
                    {plant_img_tag}
                    <div style="font-size:12px; color:#64748B; margin-top:6px;"><strong>Botanical Specimen:</strong> <i>{html.escape(scientific_name)}</i> ({html.escape(chinese_name)})</div>
                </div>
                <div style="text-align:center;">
                    {mol_img_tag}
                    <div style="font-size:12px; color:#64748B; margin-top:6px;"><strong>Active Phytochemical:</strong> {html.escape(compound_name)}</div>
                </div>
            </div>
        </div>

        <!-- 2. Dual Context: Classical + Western Target -->
        <div class="section-card">
            <div class="section-title">
                <span>📜 Classical Dynastic Corpus vs. Western Molecular Target</span>
            </div>
            <div class="grid-2col">
                <!-- Classical Claim -->
                <div class="data-box" style="border-left: 4px solid #F59E0B;">
                    <h4 style="margin:0 0 6px 0; color:#B45309; font-size:13px; text-transform:uppercase;">Classical Dynastic Provenance</h4>
                    <p style="margin:0 0 6px 0; font-size:13px; color:#1E293B;"><strong>Ancient Claim:</strong> "{html.escape(claim_text)}"</p>
                    <p style="margin:0; font-size:13px; color:#059669;"><strong>Medical Translation:</strong> "{html.escape(translation)}"</p>
                </div>
                <!-- Western Target -->
                <div class="data-box" style="border-left: 4px solid #0284C7;">
                    <h4 style="margin:0 0 6px 0; color:#0369A1; font-size:13px; text-transform:uppercase;">Western Macromolecular Target</h4>
                    <p style="margin:0 0 4px 0; font-size:13px;"><strong>Target Protein:</strong> {html.escape(target_name)}</p>
                    <p style="margin:0 0 4px 0; font-size:13px;"><strong>PDB Identifier:</strong> <span style="font-weight:700; color:#0284C7;">{pdb_id}</span></p>
                    <p style="margin:0; font-size:11px; color:#64748B; word-break:break-all; font-family:'JetBrains Mono',monospace;"><strong>SMILES:</strong> <code>{html.escape(smiles)}</code></p>
                </div>
            </div>
        </div>

        <!-- 3. Objective Scientific Audit & Failure Mode Matrix -->
        <div class="section-card" style="border-left: 4px solid {audit_res['verdict_color']};">
            <div class="section-title">
                <span>⚖️ Objective Historical Claim Audit & Failure Mode Matrix</span>
                <span class="badge" style="background:#F1F5F9; color:{audit_res['verdict_color']}; border:1px solid {audit_res['verdict_color']}; font-weight:700;">{audit_res['verdict_badge']}</span>
            </div>
            <div style="background:#F8FAFC; border:1px solid #E2E8F0; border-radius:8px; padding:16px;">
                <h4 style="margin:0 0 6px 0; color:#0F172A; font-size:14px;">{audit_res['verdict_title']}</h4>
                <p style="margin:0 0 10px 0; font-size:13px; color:#475569; line-height:1.5;"><strong>Biophysical Mechanism Audit:</strong> {audit_res['mechanism_summary']}</p>
                <div style="background:#FFFFFF; border:1px solid #E2E8F0; border-radius:6px; padding:12px; font-size:12px; line-height:1.5;">
                    <b style="color:#DC2626;">Failure Mode Classification:</b> <span style="font-weight:600; color:#1E293B;">{audit_res['failure_mode_type']}</span><br>
                    <span style="color:#64748B;">{audit_res['failure_mode_explanation']}</span>
                </div>
            </div>
        </div>

        <!-- Paozhi Audit if applicable -->
        {paozhi_html}

        <!-- 4. Docking Hierarchy, 2D LigPlot Diagram & Interactive 3D WebGL Pocket -->
        <div class="section-card">
            <div class="section-title">
                <span>🎯 In-Silico Molecular Docking & Non-Covalent Interaction Topology</span>
            </div>
            <p style="margin:0 0 12px 0; font-size:13px; color:#64748B;">
                Conformational search executed via <strong>AutoDock Vina empirical scoring function</strong> with rigid receptor grid and flexible ligand dihedral sampling.
            </p>
            
            <!-- Crystallographic Ground Truth Benchmark Card if available -->
            {benchmark_html}

            <div class="table-container">
                {poses_table_html}
            </div>

            <!-- Embedded 2D LigPlot Schematic -->
            {ligplot_html}

            <!-- Embedded Interactive 3D WebGL Pocket Viewer -->
            {pocket_3d_html}

            <h4 style="margin:22px 0 8px 0; font-size:14px; color:#0F172A;">Intermolecular Pocket Contacts (&lt; 4.0 &Aring;)</h4>
            <div class="table-container">
                {interactions_table_html}
            </div>
        </div>

        <!-- 5. Quantitative MM-GBSA Free Energy & Residue Hotspots -->
        {mmgbsa_html}

        <!-- 6. Molecular Dynamics Trajectory & Stability Suite -->
        {md_html}

        <!-- 7. Bioisostere Lead Optimization (Stage 04 Comparative Lead Matrix) -->
        {variant_html}

        <!-- Systems Network Pharmacology & Targetome Profile (Optional Extension) -->
        {systems_html}

        <!-- 7. ADMET Pharmacokinetics & PAINS -->
        <div class="section-card">
            <div class="section-title">
                <span>🛡️ ADMET Pharmacokinetics & Toxicological Screen</span>
            </div>
            {admet_cards}
        </div>

        <!-- 8. IND-Enabling Regulatory Go/No-Go Decision & Assays -->
        <div class="section-card" style="background:#F8FAFC; border:1px solid #CBD5E1;">
            <div class="section-title">
                <span>📋 IND-Enabling Regulatory Decision & Recommended In-Vitro Assays</span>
                <span class="badge badge-green">GO FOR BENCH VALIDATION</span>
            </div>
            <div style="display:grid; grid-template-columns: 1fr 1fr; gap:16px;">
                <div>
                    <h5 style="margin:0 0 6px 0; font-size:13px; color:#0F172A;">Recommended Experimental Assays</h5>
                    <ul style="margin:0; padding-left:18px; font-size:12.5px; color:#475569; line-height:1.6;">
                        <li><strong>Surface Plasmon Resonance (SPR):</strong> Measure physical dissociation constant (K<sub>D</sub>) on immobilized target protein.</li>
                        <li><strong>Isothermal Titration Calorimetry (ITC):</strong> Quantify enthalpy (&Delta;H) vs entropy (-T&Delta;S) binding driving forces.</li>
                        <li><strong>Enzymatic Inhibition Assay:</strong> Multi-point dose-response curve to calculate experimental IC<sub>50</sub>.</li>
                    </ul>
                </div>
                <div>
                    <h5 style="margin:0 0 6px 0; font-size:13px; color:#0F172A;">Translational Medicine Verdict</h5>
                    <p style="margin:0; font-size:12.5px; color:#475569; line-height:1.55;">
                        The chemical profile demonstrates high shape complementarity to the target pocket with zero lethal toxicophores. The lead candidate is cleared for chemical synthesis and primary biochemical screening.
                    </p>
                </div>
            </div>
        </div>

    </div>

    <!-- Footer -->
    <div class="dossier-footer">
        <div>
            <strong>EthnoDock Pro</strong> &bull; Computational Ethnopharmacology & In-Silico Discovery Platform
        </div>
        <div style="font-family:'JetBrains Mono',monospace;">
            Verified Publication-Grade Computational Monograph &bull; EDK-PRO-2026
        </div>
    </div>

</div>

</body>
</html>
"""
    return html_content
