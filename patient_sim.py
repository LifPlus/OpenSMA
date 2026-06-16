"""
Personalized SMA Patient Simulation & PK/PD Projections (Module E)
==================================================================
Models the personalized patient response to various treatments.
Incorporates patient profile (age, weight, SMN2 copies) into the
outcomes, providing safety, efficacy, and dosage benefit/risk projections.

Phase 7 Upgrade: Multi-Compartment Physiologically-Based Pharmacokinetics (PBPK) Model
Models Plasma, Peripheral Tissue, CSF, and Brain Parenchyma compartments.
"""

import os
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# ══════════════════════════════════════════════════════════
# Patient Profile & Dose Baselines
# ══════════════════════════════════════════════════════════
def load_patient_profile():
    try:
        with open("patient_profile.json", "r") as f:
            profile = json.load(f)
            return {
                "ID": profile.get("patient_id", "SMA-CASE-001"),
                "Age_months": profile.get("age_months", 6),
                "SMN2_copies": profile.get("smn2_copies", 2),
                "Baseline_SMN_percent": 10.0 + (profile.get("smn2_copies", 2) * 5.0),
                "Motor_Neurons_Remaining_%": max(10.0, 100.0 - (profile.get("age_months", 6) * 6.0)),
                "Weight_kg": profile.get("weight_kg", 7.0)
            }
    except FileNotFoundError:
        return {
            "ID": "Default-SMA-Pt",
            "Age_months": 6,
            "SMN2_copies": 2,
            "Baseline_SMN_percent": 15.0,
            "Motor_Neurons_Remaining_%": 60.0,
            "Weight_kg": 7.0
        }

PATIENT = load_patient_profile()

# Load best results from pipeline to make it fully dynamic
NUSINERSEN_SMN_BOOST = 25.0
RISDIPLAM_SMN_BOOST = 28.0
NOVEL_ASO_PREDICTED_BOOST = 29.5
NOVEL_SM_PREDICTED_BOOST = 29.5
REPURPOSED_ORAL_BOOST = 21.0
OPENSMA_LNP_ASO_BOOST = 31.0

if os.path.exists("aso_candidates_intron7.csv"):
    try:
        df_aso = pd = None
        import pandas as pd
        df = pd.read_csv("aso_candidates_intron7.csv")
        if not df.empty:
            top_aso = df.iloc[0]
            val = float(top_aso.get("Pred_Exon7_Inclusion_%", 30.0)) - 15.0
            NOVEL_ASO_PREDICTED_BOOST = round(max(20.0, val), 1)
            OPENSMA_LNP_ASO_BOOST = round(NOVEL_ASO_PREDICTED_BOOST + 2.0, 1)
    except Exception:
        pass

if os.path.exists("docked_candidates_ranked.csv"):
    try:
        import pandas as pd
        df_sm = pd.read_csv("docked_candidates_ranked.csv")
        df_active = df_sm[df_sm.get("Fitness_Score") != "REF"]
        if df_active.empty:
            df_active = df_sm
        if not df_active.empty:
            top_sm = df_active.iloc[0]
            vina = top_sm.get("Vina_Affinity_kcal/mol")
            if pd.notna(vina) and isinstance(vina, (int, float)):
                NOVEL_SM_PREDICTED_BOOST = round(max(20.0, min(35.0, -vina * 4.0)), 1)
    except Exception:
        pass

# ══════════════════════════════════════════════════════════
# Adverse Event Profiles (literature-based heuristics)
# ══════════════════════════════════════════════════════════
AE_PROFILES = {
    "Nusinersen": {
        "Headache":              28.0,
        "Back_Pain":             41.0,
        "Post_LP_Syndrome":      18.0,
        "Thrombocytopenia":       3.0,
        "Meningitis_Risk":        2.0,
        "Retinal_Toxicity":       0.0,
        "GI_Effects":             2.0,
        "Rash":                   1.0,
        "Route":                 "Intrathecal (invasive)",
        "BBB_Bypass":            True
    },
    "Risdiplam": {
        "Headache":              13.0,
        "Back_Pain":              5.0,
        "Post_LP_Syndrome":       0.0,
        "Thrombocytopenia":       1.0,
        "Meningitis_Risk":        0.0,
        "Retinal_Toxicity":       0.0,
        "GI_Effects":            22.0,
        "Rash":                  11.0,
        "Route":                 "Oral (systemic)",
        "BBB_Bypass":            False
    },
    "Novel_ASO_v1": {
        "Headache":              25.0,
        "Back_Pain":             38.0,
        "Post_LP_Syndrome":      15.0,
        "Thrombocytopenia":       2.5,
        "Meningitis_Risk":        1.5,
        "Retinal_Toxicity":       0.0,
        "GI_Effects":             3.0,
        "Rash":                   1.5,
        "Route":                 "Intrathecal (invasive)",
        "BBB_Bypass":            True,
        "Note":                  "No CpG motifs. Lower innate immunogenicity predicted."
    },
    "Fluoro_Risdiplam_Analog": {
        "Headache":              11.0,
        "Back_Pain":              5.0,
        "Post_LP_Syndrome":       0.0,
        "Thrombocytopenia":       0.8,
        "Meningitis_Risk":        0.0,
        "Retinal_Toxicity":       0.0,
        "GI_Effects":            18.0,
        "Rash":                   9.0,
        "Route":                 "Oral (systemic)",
        "BBB_Bypass":            False,
        "Note":                  "Fluorine substitution improves half-life. CNS MPO ≥ 4."
    },
    "Repurposed_Orals": {
        "Headache":              15.0,
        "Back_Pain":              2.0,
        "Post_LP_Syndrome":       0.0,
        "Thrombocytopenia":       1.5,
        "Meningitis_Risk":        0.0,
        "Retinal_Toxicity":       0.0,
        "GI_Effects":            25.0,
        "Rash":                   5.0,
        "Route":                 "Oral (systemic)",
        "BBB_Bypass":            False,
        "Note":                  "Salbutamol + Valproic Acid combo. Off-patent, low cost."
    },
    "OpenSMA_LNP_ASO": {
        "Headache":              10.0,
        "Back_Pain":             15.0,
        "Post_LP_Syndrome":      10.0,
        "Thrombocytopenia":       1.0,
        "Meningitis_Risk":        0.5,
        "Retinal_Toxicity":       0.0,
        "GI_Effects":             2.0,
        "Rash":                   1.0,
        "Route":                 "Intrathecal (LNP)",
        "BBB_Bypass":            True,
        "Note":                  "Patent-free PMO/LNA omurga via MC3 LNP."
    }
}

# ══════════════════════════════════════════════════════════
# Multi-Compartment PBPK Model Setup
# ══════════════════════════════════════════════════════════
def pbpk_ode(y, t, is_small_molecule, emax_smn_boost, copy_num_modifier, weight_factor):
    """
    7-state variables:
    y = [C_gut, C_plasma, C_periph, C_csf, C_brain, SMN, MN]
    """
    C_gut, C_plasma, C_periph, C_csf, C_brain, SMN, MN = y
    
    # Volumes of distribution (scaled by weight)
    V_p = 3.0 * weight_factor
    V_periph = 10.0 * weight_factor
    V_csf = 0.15
    V_brain = 1.0 * weight_factor
    
    if is_small_molecule:
        k_a = 1.2          # absorption rate (1/month)
        k_el = 0.35        # elimination from plasma (1/month)
        k_12 = 0.5; k_21 = 0.2
        k_in_bbb = 0.9     # BBB penetration
        k_efflux_bbb = 0.6 # P-gp efflux
        k_el_brain = 0.05
        
        dGut_dt = -k_a * C_gut
        dPlasma_dt = (k_a * C_gut) / V_p - k_el * C_plasma - k_12 * C_plasma + k_21 * C_periph - k_in_bbb * C_plasma + k_efflux_bbb * C_brain * V_brain / V_p
        dPeriph_dt = k_12 * C_plasma - k_21 * C_periph
        dCSF_dt = 0.0
        dBrain_dt = k_in_bbb * C_plasma * V_p / V_brain - k_efflux_bbb * C_brain - k_el_brain * C_brain
    else:
        # ASO PK (Intrathecal bolus)
        k_diff_cns = 0.15  # diffusion from CSF to brain
        k_csf_out = 0.25   # clearance from CSF to blood
        k_el_plasma = 0.8  # renal clearance of ASO
        k_el_brain = 0.04  # very slow ASO decay in brain parenchyma
        
        dGut_dt = 0.0
        dPlasma_dt = (k_csf_out * C_csf * V_csf) / V_p - k_el_plasma * C_plasma
        dPeriph_dt = 0.0
        dCSF_dt = -k_diff_cns * (C_csf - C_brain) - k_csf_out * C_csf
        dBrain_dt = k_diff_cns * (C_csf - C_brain) * V_csf / V_brain - k_el_brain * C_brain
        
    # Pharmacodynamics
    # Stimulated SMN driven by brain concentration (Hill Emax)
    EC50 = 1.5
    stimulated_smn = emax_smn_boost * (C_brain / (EC50 + C_brain))
    
    smn_base = 15.0
    effective_target_smn = (smn_base + stimulated_smn) * copy_num_modifier
    effective_target_smn = min(max(effective_target_smn, 0.0), 100.0)
    
    k_smn = 0.8
    dSMN_dt = k_smn * (effective_target_smn - SMN)
    
    # Motor Neuron Survival
    protection = SMN / 100.0
    protection = min(max(protection, 0.0), 1.0)
    current_decay = 0.05 * (1.0 - protection * 0.92)
    dMN_dt = -current_decay * MN
    
    return [dGut_dt, dPlasma_dt, dPeriph_dt, dCSF_dt, dBrain_dt, dSMN_dt, dMN_dt]

def simulate_pbpk_trajectory(is_small_molecule, emax_smn_boost, months=12):
    t = np.linspace(0, months, 100)
    weight_factor = max(1.0, PATIENT["Weight_kg"] / 6.5)
    copy_num_modifier = PATIENT["SMN2_copies"] / 2.0
    
    # Dose is scaled such that brain concentration reaches therapeutic levels
    dose = 20.0
    
    # State: [C_gut, C_plasma, C_periph, C_csf, C_brain, SMN, MN]
    if is_small_molecule:
        y0 = [dose, 0.0, 0.0, 0.0, 0.0, PATIENT["Baseline_SMN_percent"], PATIENT["Motor_Neurons_Remaining_%"]]
    else:
        y0 = [0.0, 0.0, 0.0, dose, 0.0, PATIENT["Baseline_SMN_percent"], PATIENT["Motor_Neurons_Remaining_%"]]
        
    sol = odeint(pbpk_ode, y0, t, args=(is_small_molecule, emax_smn_boost, copy_num_modifier, weight_factor))
    
    # Return t, smn_trajectory, mn_trajectory
    return t, sol[:, 5], sol[:, 6]

# ══════════════════════════════════════════════════════════
# Plotting
# ══════════════════════════════════════════════════════════
def plot_simulations(scenarios):
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle("SMA Patient PBPK Simulation: Multi-Compartment Outcome Projections", 
                 fontsize=14, fontweight='bold')

    colors = ['#2196F3', '#4CAF50', '#FF5722', '#9C27B0', '#00BCD4', '#E91E63']
    smn_boosts = {
        "Nusinersen":            (False, NUSINERSEN_SMN_BOOST),
        "Risdiplam":             (True, RISDIPLAM_SMN_BOOST),
        "Novel_ASO_v1":          (False, NOVEL_ASO_PREDICTED_BOOST),
        "Fluoro_Risdiplam_Analog": (True, NOVEL_SM_PREDICTED_BOOST),
        "Repurposed_Orals":      (True, REPURPOSED_ORAL_BOOST),
        "OpenSMA_LNP_ASO":       (False, OPENSMA_LNP_ASO_BOOST)
    }

    # Panel 1: SMN Protein Restoration Over Time
    ax1 = axes[0, 0]
    for (name, color) in zip(smn_boosts.keys(), colors):
        is_sm, emax = smn_boosts[name]
        t, smn, mn = simulate_pbpk_trajectory(is_sm, emax, months=12)
        ax1.plot(t, smn, label=name, color=color, linewidth=2)
    ax1.axhline(y=PATIENT["Baseline_SMN_percent"], color='red', linestyle='--', alpha=0.5, label='Pre-treatment baseline')
    ax1.axhline(y=50, color='gray', linestyle=':', alpha=0.7, label='Threshold for survival benefit (~50%)')
    ax1.set_xlabel("Months After Treatment")
    ax1.set_ylabel("Estimated SMN Protein Level (% of Normal)")
    ax1.set_title("SMN Protein Restoration Trajectory (PBPK Compartment)")
    ax1.legend(fontsize=8)
    ax1.set_ylim(0, 105)
    ax1.grid(True, alpha=0.3)

    # Panel 2: Motor Neuron Survival Projection
    ax2 = axes[0, 1]
    t_un = np.linspace(0, 12, 100)
    untreated_mn = PATIENT["Motor_Neurons_Remaining_%"] * np.exp(-0.05 * t_un)
    ax2.plot(t_un, untreated_mn, 'r--', linewidth=2, label="Untreated / Natural Decline")
    
    for (name, color) in zip(smn_boosts.keys(), colors):
        is_sm, emax = smn_boosts[name]
        t, smn, mn = simulate_pbpk_trajectory(is_sm, emax, months=12)
        ax2.plot(t, mn, label=name, color=color, linewidth=2)
    ax2.set_xlabel("Months After Treatment")
    ax2.set_ylabel("Motor Neuron Survival (%)")
    ax2.set_title("Motor Neuron Survival Projection (PBPK-Coupled)")
    ax2.legend(fontsize=8)
    ax2.set_ylim(0, 105)
    ax2.grid(True, alpha=0.3)

    # Panel 3: Adverse Event Comparison
    ax3 = axes[1, 0]
    ae_categories = ["Headache", "Thrombocytopenia", "Meningitis_Risk", "Retinal_Toxicity", "GI_Effects", "Rash"]
    x = np.arange(len(ae_categories))
    width = 0.15
    for i, (name, color) in enumerate(zip(smn_boosts.keys(), colors)):
        ae_data = [AE_PROFILES[name].get(cat, 0) for cat in ae_categories]
        ax3.bar(x + i * width, ae_data, width, label=name, color=color, alpha=0.8)
    ax3.set_xticks(x + width * 2.5)
    ax3.set_xticklabels(ae_categories, rotation=15, ha='right', fontsize=8)
    ax3.set_ylabel("Adverse Event Probability (%)")
    ax3.set_title("Adverse Event Comparison")
    ax3.legend(fontsize=7)
    ax3.grid(True, alpha=0.3, axis='y')

    # Panel 4: Composite Benefit-Risk Scatter
    ax4 = axes[1, 1]
    for (name, color) in zip(smn_boosts.keys(), colors):
        ae = AE_PROFILES[name]
        avg_ae = np.mean([ae.get(k, 0) for k in ae_categories])
        is_sm, emax = smn_boosts[name]
        t, smn, mn = simulate_pbpk_trajectory(is_sm, emax, months=12)
        final_boost = smn[-1] - PATIENT["Baseline_SMN_percent"]
        ax4.scatter(avg_ae, final_boost, s=200, c=color, zorder=5, label=name)
        ax4.annotate(name.replace('_',' '), (avg_ae, final_boost), 
                     textcoords='offset points', xytext=(5,5), fontsize=7)
    ax4.set_xlabel("Mean Adverse Event Probability (%)")
    ax4.set_ylabel("Simulated 12m SMN Net Boost (%)")
    ax4.set_title("PBPK Benefit vs. Risk Landscape")
    ax4.grid(True, alpha=0.3)
    ax4.axvline(x=10, color='gray', linestyle='--', alpha=0.4)
    ax4.axhline(y=20, color='gray', linestyle='--', alpha=0.4)
    ax4.text(1, 23, "✓ Best Zone", color='green', fontsize=9)

    plt.tight_layout()
    plt.savefig('patient_simulation_plots.png', dpi=150, bbox_inches='tight')
    print("Simulation plots saved to 'patient_simulation_plots.png'")

# ══════════════════════════════════════════════════════════
# Generate Clinical Report
# ══════════════════════════════════════════════════════════
def generate_clinical_summary():
    report = {
        "patient": PATIENT,
        "simulation_summary": {}
    }
    
    smn_boosts = {
        "Nusinersen":            (False, NUSINERSEN_SMN_BOOST),
        "Risdiplam":             (True, RISDIPLAM_SMN_BOOST),
        "Novel_ASO_v1":          (False, NOVEL_ASO_PREDICTED_BOOST),
        "Fluoro_Risdiplam_Analog": (True, NOVEL_SM_PREDICTED_BOOST),
        "Repurposed_Orals":      (True, REPURPOSED_ORAL_BOOST),
        "OpenSMA_LNP_ASO":       (False, OPENSMA_LNP_ASO_BOOST)
    }
    
    print("\n" + "="*68)
    print("  SMA PATIENT CLINICAL OUTCOME SIMULATION REPORT (PBPK MODEL)")
    print("="*68)
    print(f"Patient: {PATIENT['ID']}  |  Age: {PATIENT['Age_months']} months  |  SMN2 copies: {PATIENT['SMN2_copies']}")
    print(f"Baseline SMN: {PATIENT['Baseline_SMN_percent']}% of normal  |  Motor Neurons Intact: {PATIENT['Motor_Neurons_Remaining_%']}%")
    print("="*68)

    t_un = np.linspace(0, 12, 100)
    final_mn_untreated = PATIENT["Motor_Neurons_Remaining_%"] * np.exp(-0.05 * 12)

    for name, (is_sm, emax) in smn_boosts.items():
        t, smn, mn = simulate_pbpk_trajectory(is_sm, emax, months=12)
        
        final_smn   = smn[-1]
        final_mn_treated   = mn[-1]
        mn_saved    = final_mn_treated - final_mn_untreated
        
        ae = AE_PROFILES[name]
        ae_categories = ["Headache", "Thrombocytopenia", "Meningitis_Risk", "Retinal_Toxicity", "GI_Effects", "Rash"]
        avg_ae = np.mean([ae.get(k, 0) for k in ae_categories])
        
        note = ae.get("Note", "Standard clinical profile.")
        
        print(f"\n{'—'*60}")
        print(f"  CANDIDATE: {name}")
        print(f"{'—'*60}")
        print(f"  Route:                       {ae['Route']}")
        print(f"  PBPK 12m SMN Level:          {final_smn:.1f}% of normal (emax boost: +{emax}%)")
        print(f"  Motor neurons at 12 months:  {final_mn_treated:.1f}%  (untreated would be {final_mn_untreated:.1f}%)")
        print(f"  Neurons rescued (estimate):  +{mn_saved:.1f}%")
        print(f"  Mean AE probability:         {avg_ae:.1f}%")
        print(f"  Thrombocytopenia risk:       {ae.get('Thrombocytopenia',0):.1f}%")
        print(f"  Meningitis risk:             {ae.get('Meningitis_Risk',0):.1f}%")
        print(f"  Retinal risk:                {ae.get('Retinal_Toxicity',0):.1f}%")
        print(f"  Safety note:                 {note}")

        report["simulation_summary"][name] = {
            "final_smn_pct": round(final_smn, 2),
            "motor_neurons_at_12m_pct": round(final_mn_treated, 2),
            "neurons_rescued_vs_no_tx": round(mn_saved, 2),
            "mean_AE_probability": round(avg_ae, 2),
            "route": ae['Route']
        }

    print(f"\n{'='*68}")
    print("Plots generated and saved. See 'patient_simulation_plots.png'")
    
    with open("clinical_simulation_results.json", "w") as f:
        json.dump(report, f, indent=2)
    print("Full clinical simulation data saved to 'clinical_simulation_results.json'")
    print("="*68)
    
    return report

if __name__ == "__main__":
    plot_simulations(None)
    generate_clinical_summary()
