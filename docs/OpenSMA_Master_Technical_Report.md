# OpenSMA: End-to-End Computational Drug Discovery & Systems Biology Simulation Platform Master Technical Report

---

## FOREWORD AND PROJECT SCOPE
This document serves as the master technical report detailing the architectural, algorithmic, mathematical, and biological foundations of **OpenSMA**, an in-silico drug discovery and personalized systems biology platform designed to screen therapies for Spinal Muscular Atrophy (SMA). Built on open-source libraries (RDKit, Biopython, AutoDock Vina, Scikit-Learn, SciPy), the platform accelerates pre-clinical screening and simulates longitudinal patient-specific response kinetics.

---

## 1. INTRODUCTION TO SPINAL MUSCULAR ATROPHY (SMA) BIOLOGY

### 1.1. Molecular Pathology and Genetics
Spinal Muscular Atrophi (SMA) is a recessive neuromuscular disorder characterized by the progressive degeneration of anterior horn motor neurons in the spinal cord, leading to proximal muscle wasting and respiratory failure. The primary cause of the disease is the homozygous deletion or loss-of-function mutation of the *SMN1* (Survival Motor Neuron 1) gene located on chromosome 5q13.2.

The human genome contains an evolutionary duplication of *SMN1* called *SMN2*. However, *SMN2* cannot fully compensate for *SMN1* loss due to a single C $\to$ T transition at position +6 in Exon 7.

```
Genomic Comparison:
SMN1 Exon 7: ... GCT ATC ATA ATT TTG TTC ATC [C] GA CAA AAT ... (SRSF1 binding ESE)
                                          ▼
SMN2 Exon 7: ... GCT ATC ATA ATT TTG TTC ATC [T] GA CAA AAT ... (SRSF1 inactive / hnRNP A1 active)
```

This single nucleotide transition:
1.  Disrupts an Exonic Splicing Enhancer (ESE) recognized by splicing activator **SRSF1** (SF2/ASF).
2.  Creates an Exonic Splicing Silencer (ESS) that recruits **hnRNP A1** to the nearby Intronic Splicing Silencer N1 (**ISS-N1**) region in Intron 7.

Consequently, ~90% of *SMN2* transcription skips Exon 7, producing a truncated, unstable protein isoform ($\Delta7$ SMN) that undergoes rapid proteasomal degradation. Only ~10% results in functional full-length (FL) SMN protein.

### 1.2. Therapeutic Modalities and Splicing modulators
Current FDA-approved therapies target *SMN2* splicing or replace the *SMN1* gene:
*   **Nusinersen (Spinraza):** An intrathecally delivered Antisense Oligonucleotide (ASO) blocking the ISS-N1 silencer, promoting Exon 7 inclusion.
*   **Risdiplam (Evrysdi):** A small molecule oral modulator stabilizing the transient U1 snRNP-RNA complex.
*   **Onasemnogene Abeparvovec (Zolgensma):** An AAV9-mediated gene replacement therapy.

---

## 2. OPENSMA PIPELINE ARCHITECTURE

The platform implements an integrated 12-phase pipeline spanning chemical space generation, thermodynamic analysis, structural docking, machine learning safety filters, systems pharmacology, and cohort simulations.

```
                                  [ Start Pipeline ]
                                          │
                  ┌───────────────────────┴───────────────────────┐
                  ▼                                               ▼
         [ Phase 1: Small Molecule ]                      [ Phase 2: ASO Design ]
         (Genetic Algorithm Scaffolds)                   (Sliding Window / Tm NN)
                  │                                               │
                  ▼                                               ▼
         [ Phase 4: ADMET Screening ]                    [ Phase 5: ASO Toxicity ]
         (Lipinski, CNS MPO, hERG, PAINS)                (CpG, G-Quadruplex, Poly-run)
                  │                                               │
                  └───────────────────────┬───────────────────────┘
                                          ▼
                              [ Phase 3: CRISPR / BE ]
                              (SpCas9, SaCas9, CBE gRNA)
                                          │
                                          ▼
                              [ Phase 6: PDB 8R62 Docking ]
                              (Meeko Conversion + Vina API)
                                          │
                                          ▼
                              [ Phase 7: Tox21 ML Screening ]
                              (Morgan FP + Descriptors + GBM)
                                          │
                                          ▼
                              [ Phase 8: ODE & Off-Target ]
                              (scipy.optimize + Bio.pairwise2)
                                          │
                                          ▼
                              [ Phase 9: Composite Scoring ]
                              (Efficacy/Safety/Delivery Matrix)
                                          │
                                          ▼
                              [ Phase 10: Patient Profile ]
                              (patient_profile.json dynamic PK)
                                          │
                                          ▼
                              [ Phase 11: Full Cure Sim ]
                              (5-Year Longitudinal Prognosis)
                                          │
                                   [ End Pipeline ]
```

---

## 3. ALGORITHMS AND MATHEMATICAL FORMULATIONS

### 3.1. Phase 1: Small Molecule Optimization (`molecule_optimizer.py`)
This module performs scaffold-based modifications on the core structure of Risdiplam using a Genetic Algorithm (GA).

#### Mutation Rules (Reaction SMARTS)
1.  **Fluorination:** `[c:1][H]>>[c:1]F` (Improves metabolic stability).
2.  **Methylation:** `[c:1][H]>>[c:1]C` (Modulates lipophilicity).
3.  **Methoxylation:** `[c:1][H]>>[c:1]OC` (Modulates polar surface area).
4.  **Chlorination:** `[c:1][H]>>[c:1]Cl` (Halogen binding potential).
5.  **Aromatic Nitrogen Swap:** `[c:1][H]>>[n:1]` (Alters ring basicity).

#### Composite Fitness Function
$$\text{Fitness} = (\text{QED} \times 100) + S_{\text{BBB\_bonus}} - S_{\text{penalties}} + S_{\text{3D\_stability}}$$

Where:
-   **QED:** Quantitative Estimate of Druglikeness (Bickerton et al., 2012).
-   **$S_{\text{BBB\_bonus}}$:** $+20$ for $2.0 \le \text{LogP} \le 5.0$, $+20$ for $\text{TPSA} < 90\text{ Å}^2$, and $+10$ for $\text{MW} < 500\text{ Da}$.
-   **$S_{\text{penalties}}$:** $-10$ for $\text{HBD} > 3$, $-10$ for $\text{HBA} > 7$, $-20$ for $\text{LogP} > 5.5$, and $-30$ for $\text{MW} > 600\text{ Da}$.
-   **$S_{\text{3D\_stability}}$:** RDKit MMFF94 force-field geometry optimization. Conformer strain energy $< 150 \text{ kcal/mol}$ adds $+15$, while $> 300 \text{ kcal/mol}$ subtracts $-10$. Failed 3D embedding subtracts $-15$.

### 3.2. Phase 2: Antisense Oligonucleotide (ASO) Design (`aso_designer.py`)
Scans Intron 7 targeting the ISS-N1 region (sequence: `5'-YTTTATAATGY-3'`) using a sliding window.

#### Hybridization Thermodynamics (Freier 1986)
$$\Delta G^\circ_{37\text{, total}} = \Delta G^\circ_{\text{init}} + \sum_{i=1}^{L-1} \Delta G^\circ_{\text{doublet}}(i, i+1) + N_{\text{term\_AU}} \cdot \Delta G^\circ_{\text{term\_AU}}$$

Doublets like **GG/CC** ($\Delta G^\circ_{37} = -3.1 \text{ kcal/mol}$) or **GC/CG** ($\Delta G^\circ_{37} = -3.4 \text{ kcal/mol}$) are summed with helix nucleation ($\Delta G^\circ_{\text{init}} = +3.4 \text{ kcal/mol}$) and terminal AU pair penalty ($\Delta G^\circ_{\text{term\_AU}} = +0.5 \text{ kcal/mol}$).

### 3.3. Phase 3: CRISPR/gRNA Design (`crispr_designer.py`)
Designs gRNAs targeting the *SMN2* genomic region to edit the +6 site (T $\to$ C) or disrupt splicing silencers. Supports **SpCas9** (PAM: `NGG`) and **SaCas9** (PAM: `NNGRRT`). Scores guides based on GC content ($40\% - 60\%$) and homopolymer run penalties.

### 3.4. Phase 4: ADMET Screening (`admet_screener.py`)
Filters GA-generated small molecules:
-   **Lipinski's Rule of 5:** $\text{MW} \le 500$, $\text{LogP} \le 5$, $\text{HBD} \le 5$, $\text{HBA} \le 10$.
-   **CNS MPO (Central Nervous System Multiparameter Optimization):** Score $\ge 4.0$ (scale 0-6).
-   **PAINS Filters:** Substructure matching to exclude pan-assay interference compounds.

### 3.5. Phase 5: ASO Toxicity Filter (`aso_toxicity.py`)
ASO candidate filters:
-   **Poly-run limit:** Max 4 consecutive identical bases (prevents off-target structural traps).
-   **CpG Island count:** Flagged if CpG density $> 0.15$ (minimizes toll-like receptor mediated immune stimulation).
-   **G-Quadruplex motifs:** Reject sequences forming $G_xN_yG_xN_zG_xN_wG_x$ (prevents non-antisense tertiary aggregates).

### 3.6. Phase 6: PDB 8R62 Molecular Docking (`docking_engine.py`)
Performs ligand preparation with Meeko and docking using AutoDock Vina against the SMN2 pre-mRNA pocket (PDB 8R62). Pocket dimensions: center $(22.1, 4.3, 11.2)$, size $(20, 20, 20)\text{ Å}$.

### 3.7. Phase 7: Machine Learning Toxicity Model (`toxicity_model.py`)
Trained on 1,999 Tox21 compounds for Mitochondrial Membrane Potential Disruption (SR-MMP). Utilizes a Gradient Boosting Classifier combining Morgan Fingerprints (ECFP4, 1024-bits) and 8 physical descriptors, achieving a cross-validated **ROC-AUC of 0.875**.

### 3.8. Phase 8: PBPK ODE Calibration & Genomic Off-Target (`ode_calibration_and_screening.py`)
Calibrates 5-compartment PBPK equations against clinical trial literature (Finkel 2017 & Baranello 2021) using `scipy.optimize.curve_fit`. Fits Nusinersen bulk clearance and Risdiplam absorption coefficients. Performs local sequence alignments (Smith-Waterman) to filter off-target matches ($>80\%$ homology for gRNA, $>85\%$ for ASOs).

### 3.9. Phase 9: Composite Scoring and Selection (`final_scorer.py`)
Selects top therapeutics using a normalized matrix balancing binding energy, thermodynamic affinity, safety indices, and PBPK brain exposure profile.

### 3.10. Phase 10: Patient Profile Generation (`patient_analyzer.py`)
Generates virtual patient specifications (`patient_profile.json`) incorporating age (infant vs adult), weight, and genetic cofactor parameters to scale volume compartments dynamically.

### 3.11. Phase 11: Full Cure 5-Year Sim (`full_cure_sim.py`)
Couples PBPK concentration profiles with Hill pharmacodynamics and motor neuron survival rates:

$$\frac{d[MN]}{dt} = -r_{\text{decay}} \cdot \left(1.0 - \text{protection}(t) \cdot \eta \right) \cdot [MN]$$

Integrates Markov transition matrices over 5 years (60 months) to estimate cohort outcomes (survival and independent walking probability).

---

## 4. EXPERIMENTAL VALIDATION WET-LAB PROTOCOLS

To bridge computational predictions with clinical reality, we have prepared a detailed wet-lab protocol (`OpenSMA_Lab_Protocol.md`):
-   **ASO-LNP Formulation:** Synthesis using microfluidics (DLin-MC3-DMA, DSPC, Cholesterol, DMG-PEG2000 at 50:10:38.5:1.5 molar ratio).
-   **GM03813 Transfection:** Lipofectamine/LNP delivery into patient-derived fibroblasts.
-   **Splicing Assay:** Semi-quantitative RT-PCR using Exon 6 and Exon 8 flanking primers, resolving bands at **263 bp** (full length) and **209 bp** ($\Delta7$).
-   **SMI-32 Immunofluorescence:** Quantification of primary motor neuron survival in co-culture systems.

---

## 5. RESEARCH CONCLUSIONS AND ROADMAP

The OpenSMA computational pipeline bridges structural molecular design with clinical systems biology. By offering defensive open-source disclosures of therapeutic candidates and detailed synthesis pathways, the platform bypasses conventional patent barriers and reduces drug discovery cost overheads. Future development aims to extend this multi-compartment simulation model to other splicing diseases, such as Duchenne Muscular Dystrophy (DMD).
