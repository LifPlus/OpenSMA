"""
OpenSMA — Unified End-to-End Therapeutic Discovery Pipeline
=============================================================
Master pipeline that orchestrates ALL OpenSMA modules in the
correct dependency order and generates a unified final report.

Run order:
  Phase 1  — De Novo Molecule Generation      (molecule_optimizer.py)
  Phase 2  — Smart ASO Design                 (aso_designer.py)
  Phase 3  — CRISPR/Base Editing Design       (crispr_designer.py)
  Phase 4  — ADMET Screening                  (admet_screener.py)
  Phase 5  — ASO Toxicity Screen              (aso_toxicity.py)
  Phase 6  — AutoDock Vina 3D Docking         (docking_engine.py)
  Phase 7  — Real Tox21 ML Toxicity           (toxicity_model.py)
  Phase 8  — ODE Calibration & Off-Target     (ode_calibration_and_screening.py)
  Phase 9  — Final Composite Scoring          (final_scorer.py)
  Phase 10 — Patient Clinical Simulation      (patient_sim.py)
  Phase 11 — Full Cure Combination Sim        (full_cure_sim.py)

Output:
  - opensma_final_report.json  (machine-readable)
  - opensma_final_report.txt   (human-readable summary)
"""

import os
import sys
import json
import time
import subprocess

# ──────────────────────────────────────────────────────────────
# Color helpers for terminal output
# ──────────────────────────────────────────────────────────────
GREEN  = "\033[92m"
YELLOW = "\033[93m"
RED    = "\033[91m"
CYAN   = "\033[96m"
BOLD   = "\033[1m"
RESET  = "\033[0m"

def banner(title, color=CYAN):
    width = 70
    print(f"\n{color}{BOLD}{'═' * width}{RESET}")
    print(f"{color}{BOLD}  {title}{RESET}")
    print(f"{color}{BOLD}{'═' * width}{RESET}\n")

def step(num, name):
    print(f"{YELLOW}{BOLD}[Phase {num:02d}]{RESET} {name}...")

def ok(msg="Done"):
    print(f"  {GREEN}✓ {msg}{RESET}")

def fail(msg):
    print(f"  {RED}✗ {msg}{RESET}")

def run_module(script_name, phase_num, phase_name):
    """Run a Python module script as a subprocess."""
    step(phase_num, phase_name)
    start = time.time()
    result = subprocess.run(
        [sys.executable, script_name],
        capture_output=True, text=True
    )
    elapsed = time.time() - start
    
    if result.returncode == 0:
        ok(f"Completed in {elapsed:.1f}s")
        return True, result.stdout
    else:
        # Show only the last 3 lines of error
        err_lines = result.stderr.strip().splitlines()
        short_err = "\n".join(err_lines[-3:]) if err_lines else "Unknown error"
        fail(f"Failed ({elapsed:.1f}s): {short_err}")
        return False, result.stderr

# ──────────────────────────────────────────────────────────────
# Pipeline Execution
# ──────────────────────────────────────────────────────────────
def run_full_pipeline():
    banner("OpenSMA — Unified Therapeutic Discovery Pipeline")
    print(f"  Patient profile: {BOLD}patient_profile.json{RESET}")
    
    # Load patient info early for the summary
    try:
        with open("patient_profile.json") as f:
            patient = json.load(f)
        print(f"  Patient ID  : {patient.get('patient_id', 'Unknown')}")
        print(f"  SMN2 copies : {patient.get('smn2_copies', 2)}")
        print(f"  Age         : {patient.get('age_months', '?')} months")
        print(f"  SMA Type    : {patient.get('sma_type', '?')}")
    except FileNotFoundError:
        patient = {}
        print(f"  {YELLOW}patient_profile.json not found — using defaults{RESET}")
    
    print(f"\n  Starting {BOLD}12-phase pipeline{RESET}...\n")
    
    pipeline_steps = [
        ("molecule_optimizer.py",            1, "De Novo Small Molecule Generation (Genetic Algorithm)"),
        ("aso_designer.py",                  2, "Smart ASO Design (Generative + Splicing Estimator)"),
        ("crispr_designer.py",               3, "CRISPR/Base Editing gRNA Design"),
        ("admet_screener.py",                4, "ADMET Pharmacokinetic Screening"),
        ("aso_toxicity.py",                  5, "ASO Toxicity Pre-Screen"),
        ("docking_engine.py",                6, "AutoDock Vina 3D Molecular Docking"),
        ("toxicity_model.py",                7, "Tox21 ML Toxicity Prediction (GBM / 1999 molecules)"),
        ("ode_calibration_and_screening.py", 8, "ODE Calibration (Finkel 2017, Baranello 2021) + Genome Off-Target"),
        ("final_scorer.py",                  9, "Composite Final Scoring & Ranking"),
        ("patient_sim.py",                  10, "Personalized Patient Clinical Simulation"),
        ("full_cure_sim.py",                11, "Full Cure Protocol Simulation (5-Year ODE)"),
        ("pocket_dynamics.py",               12, "Metropolis Monte Carlo 3D Binding Stability Simulator"),
    ]
    
    results = {}
    all_passed = True
    
    for script, num, name in pipeline_steps:
        if not os.path.exists(script):
            print(f"  {RED}[Phase {num:02d}]{RESET} {name} — {RED}SKIPPED{RESET}: {script} not found")
            results[f"phase_{num:02d}"] = {"status": "SKIPPED", "name": name}
            continue
        
        success, output = run_module(script, num, name)
        results[f"phase_{num:02d}"] = {
            "status": "PASS" if success else "FAIL",
            "name": name,
            "script": script,
        }
        if not success:
            all_passed = False
    
    # ──────────────────────────────────────────────────────────────
    # Collect Key Results
    # ──────────────────────────────────────────────────────────────
    banner("Collecting Results from All Phases", CYAN)
    
    final_report = {
        "pipeline_status": "PASS" if all_passed else "PARTIAL",
        "patient": patient,
        "phase_results": results,
        "key_outputs": {}
    }
    
    # Collect docking results
    if os.path.exists("docked_candidates_ranked.csv"):
        import pandas as pd
        df = pd.read_csv("docked_candidates_ranked.csv")
        numeric = df[df["Vina_Affinity_kcal/mol"].apply(lambda x: str(x).lstrip('-').replace('.','').isdigit())]
        if len(numeric) > 0:
            best_aff = numeric["Vina_Affinity_kcal/mol"].min()
            ok(f"Best Vina Docking Affinity: {best_aff} kcal/mol")
            final_report["key_outputs"]["best_docking_affinity_kcal"] = float(best_aff)
    
    # Collect ASO results  
    if os.path.exists("aso_candidates_intron7.csv"):
        import pandas as pd
        df = pd.read_csv("aso_candidates_intron7.csv")
        best_aso = df.iloc[0]
        ok(f"Best ASO Smart Score: {best_aso.get('Smart_Score', '?')}")
        ok(f"Best ASO Predicted Exon7 Inclusion: {best_aso.get('Pred_Exon7_Inclusion_%', '?')}%")
        final_report["key_outputs"]["best_aso_score"] = float(best_aso.get("Smart_Score", 0))
        final_report["key_outputs"]["best_exon7_inclusion"] = float(best_aso.get("Pred_Exon7_Inclusion_%", 0))
    
    # Collect ODE calibration
    if os.path.exists("calibrated_ode_params.json"):
        with open("calibrated_ode_params.json") as f:
            calib = json.load(f)
        ok(f"Nusinersen ODE k_in = {calib.get('Nusinersen', {}).get('k_in', '?')} | RMSE = {calib.get('Nusinersen', {}).get('SMN_RMSE_vs_literature', '?')}")
        ok(f"Risdiplam  ODE k_in = {calib.get('Risdiplam', {}).get('k_in', '?')} | RMSE = {calib.get('Risdiplam', {}).get('SMN_RMSE_vs_literature', '?')}")
        final_report["key_outputs"]["calibrated_ode"] = calib
    
    # Collect genomic off-target
    if os.path.exists("genome_off_target_report.csv"):
        import pandas as pd
        df = pd.read_csv("genome_off_target_report.csv")
        safe = df[df["Risk"].str.contains("LOW")]
        flagged = df[df["Risk"].str.contains("HIGH")]
        ok(f"Off-Target screening: {len(safe)} safe, {len(flagged)} flagged as HIGH RISK")
        final_report["key_outputs"]["safe_sequences"] = len(safe)
        final_report["key_outputs"]["flagged_sequences"] = len(flagged)
    
    # Collect full cure functional outcomes
    if os.path.exists("full_cure_functional_outcomes.json"):
        with open("full_cure_functional_outcomes.json") as f:
            fco = json.load(f)
        walk_pct = fco.get("Full_Cure_Protocol", {}).get("month_60", {}).get("Walk_Independently_%", "?")
        survival = fco.get("Full_Cure_Protocol", {}).get("month_60", {}).get("Survival_at_24m_%", "?")
        ok(f"Full Cure: Walk Independently (5yr) = {walk_pct}%")
        ok(f"Full Cure: Survival Indicator     (5yr) = {survival}%")
        final_report["key_outputs"]["full_cure_walk_pct_5yr"] = walk_pct
        final_report["key_outputs"]["full_cure_survival_5yr"] = survival
    
    # Collect Tox21 predictions
    if os.path.exists("tox21_predictions.csv"):
        import pandas as pd
        df = pd.read_csv("tox21_predictions.csv")
        final_report["key_outputs"]["tox21_predictions"] = df.to_dict(orient="records")
        
    # Collect MMC stability results
    if os.path.exists("mmc_stability_results.csv"):
        import pandas as pd
        df = pd.read_csv("mmc_stability_results.csv")
        if not df.empty:
            best_mmc = df.iloc[0]
            ok(f"Best MMC Pocket Stability: {best_mmc.get('MMC_Stability', '?')} (RMSD: {best_mmc.get('MMC_Final_RMSD_A', '?')} Å)")
            final_report["key_outputs"]["best_mmc_stability_rmsd"] = float(best_mmc.get("MMC_Final_RMSD_A", 0.0))
            final_report["key_outputs"]["best_mmc_stability_class"] = str(best_mmc.get("MMC_Stability", ""))
    
    # ──────────────────────────────────────────────────────────────
    # Save Final Report
    # ──────────────────────────────────────────────────────────────
    with open("opensma_final_report.json", "w") as f:
        json.dump(final_report, f, indent=2)
    
    # Write human-readable summary
    with open("opensma_final_report.txt", "w", encoding="utf-8") as f:
        f.write("=" * 70 + "\n")
        f.write("  OpenSMA — Unified Discovery Pipeline Final Report\n")
        f.write("=" * 70 + "\n\n")
        
        f.write(f"Pipeline Status: {final_report['pipeline_status']}\n\n")
        
        f.write("PATIENT PROFILE\n")
        for k, v in patient.items():
            f.write(f"  {k}: {v}\n")
        
        f.write("\nPHASE RESULTS\n")
        for phase, info in results.items():
            status_sym = "✓" if info["status"] == "PASS" else ("⊘" if info["status"] == "SKIPPED" else "✗")
            f.write(f"  {status_sym} {phase}: {info['name']} [{info['status']}]\n")
        
        f.write("\nKEY SCIENTIFIC OUTPUTS\n")
        for k, v in final_report["key_outputs"].items():
            if not isinstance(v, (dict, list)):
                f.write(f"  {k}: {v}\n")
    
    # Final summary banner
    banner("Pipeline Complete", GREEN if all_passed else YELLOW)
    status_str = f"{GREEN}ALL PHASES PASSED{RESET}" if all_passed else f"{YELLOW}COMPLETED WITH SOME FAILURES{RESET}"
    print(f"  Status : {status_str}")
    print(f"  Reports: {BOLD}opensma_final_report.json{RESET}")
    print(f"           {BOLD}opensma_final_report.txt{RESET}")
    print(f"\n  Generated outputs in this directory:")
    
    output_files = [
        ("admet_screened_molecules.csv",   "De novo small molecules (ADMET-filtered)"),
        ("aso_candidates_intron7.csv",     "Smart ASO candidates"),
        ("crispr_grna_candidates.csv",     "CRISPR gRNA designs"),
        ("docked_candidates_ranked.csv",   "3D Vina docking affinity scores"),
        ("tox21_predictions.csv",          "Real Tox21 toxicity predictions"),
        ("calibrated_ode_params.json",     "Literature-calibrated PK/PD parameters"),
        ("genome_off_target_report.csv",   "Human genome off-target risk report"),
        ("final_ranked_candidates.csv",    "Composite-scored final ranking"),
        ("patient_simulation_plots.png",   "Patient PK/PD simulation plots"),
        ("full_cure_simulation.png",       "5-year Full Cure Protocol trajectories"),
        ("full_cure_functional_outcomes.json", "Motor milestone probabilities"),
        ("mmc_stability_results.csv",      "Metropolis Monte Carlo 3D binding stability results"),
    ]
    
    for fname, desc in output_files:
        exists = "✓" if os.path.exists(fname) else " "
        print(f"    [{exists}] {fname:<45} {desc}")

if __name__ == "__main__":
    run_full_pipeline()
