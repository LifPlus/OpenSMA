# OpenSMA: Pipeline Architecture and Discovery Logic

## 1. System Overview
OpenSMA is an AI-driven computational biology pipeline designed to identify therapeutic candidates for Spinal Muscular Atrophy. The architecture follows a multi-stage discovery flow from genomic analysis to clinical simulation.

## 2. Pipeline Stages

### Stage 1: Genomic Mapping (`fetch_smn_sequences.py`)
*   **Logic:** Retrieval of *SMN1* and *SMN2* consensus sequences from NCBI.
*   **Action:** Identifies the critical C→T difference in Exon 7 and maps the Intron 7 Splicing Silencer (ISS-N1).
*   **Output:** Baseline FASTA sequences for ASO and CRISPR targeting.

### Stage 2: ASO Sliding-Window Discovery (`aso_designer.py`)
*   **Logic:** Thermodynamic modeling of DNA-RNA hybridization.
*   **Heuristic:** Calculates Melting Temperature (Tm) using the Nearest Neighbor model.
*   **Selection:** Scans for sequences that displace native repressors (hnRNP A1) with higher affinity than current standards.

### Stage 3: Small Molecule Optimization (`molecule_optimizer.py` & `admet_screener.py`)
*   **Logic:** Physicochemical property refinement of pyridazine scaffolds.
*   **Metrics:** LogP for CNS penetration, PAINS filter for toxicity, and CNS MPO score (targeting > 4.0).
*   **Refinement:** Automatic generation of SMILES analogs with predictive ADMET screening.

### Stage 4: Toxicological Filtering (`aso_toxicity.py`)
*   **Logic:** Sequence-based safety heuristics.
*   **Filters:** Elimination of CpG motifs (TLR9 immune activation), G-quadruplexes (off-target risk), and high GC/Poly-run stretches (cytotoxicity).

### Stage 5: CRISPR/Base Editing Design (`crispr_designer.py`)
*   **Logic:** PAM recognition and gRNA efficiency scoring.
*   **Heuristic:** Doench 2016 Rule Set 2 scoring for on-target activity.
*   **Outcome:** Identification of "Base-Editing Suitable" windows to convert *SMN2* to *SMN1* without double-strand breaks.

### Stage 6: Clinical Simulation (`full_cure_sim.py`)
*   **Logic:** Differential equation modeling of SMN protein restoration vs. motor neuron decay.
*   **Outcome:** 5-year longitudinal projection of functional milestones (breathing, sitting, walking).

---

## 3. The "OpenSMA Synergy" Logic
The core discovery of this project is the **Temporal Synergy Protocol**:
1.  **Immediate Rescue (ASO + Small Molecule):** Prevents further neurodegeneration within weeks.
2.  **Episomal Boost (Gene Addition):** Provides massive SMN surge.
3.  **Permanent Consolidation (Base Editing):** Rewrites the patient's DNA to ensure SMN levels never drop back to disease levels.

---
## 4. How to Extend This Project
Future researchers can extend this pipeline by:
- Importing real VCF files to the `patient_analyzer.py` module.
- Running MD simulations on the ASO-RNA binding complexes.
- Validating the Fluoro-Risdiplam analog via Docking simulations (AutoDock Vina).

---
*Created by OpenSMA AI Pipeline | Engineering Documentation*
