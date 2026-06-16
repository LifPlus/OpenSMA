# OpenSMA LinkedIn and Reddit Posting Templates

This file contains optimized copy for sharing the project on professional networks and scientific subreddits.

---

## 1. LinkedIn Post Template (Professional & DeSci Focused)

**Title:** Open-Source Splicing Modulation & Clinical Simulation: OpenSMA Platform is Live!

**Body:**
> How can we accelerate and democratize pre-clinical drug discovery for rare neuromuscular diseases like Spinal Muscular Atrophy (SMA), where commercial therapies cost millions of dollars?
> 
> I am thrilled to share that we have open-sourced **OpenSMA**, an end-to-end computational drug discovery pipeline and systems biology simulator for SMA.
> 
> **What is OpenSMA?**
> 🧬 **Smart ASO Design:** A sliding window scanner that computes hybridization free energy ($\Delta G^\circ_{37}$) of ASOs targeting the pre-mRNA ISS-N1 intron 7 region using Freier 1986 Nearest-Neighbor doublets.
> 🧪 **3D Molecular Docking & Pocket Dynamics:** AutoDock Vina molecular docking against the experimental PDB 8R62 structure, coupled with 1000-step Metropolis Monte Carlo dynamics simulating thermal stability at body temperature ($310\text{ K}$).
> 📊 **Physiologically-Based Pharmacokinetics (PBPK):** A 5-compartment continuous mass-transfer ODE system. PK parameters are fitted and calibrated against Phase 3 clinical trials (Finkel 2017 & Baranello 2021) with high correlation (RMSE of $0.33$ and $0.28$).
> 📈 **Full Cure Splicing PD & Monte Carlo Cohorts:** Long-term virtual trial projection (10,000 patients) for a multi-modal combination therapy combining transient ASO induction, oral modulators, and CRISPR Base Editing.
> 
> **Why Open Science and DeSci?**
> Proprietary drug discovery models and high pricing structures slow down rare disease research. OpenSMA releases all pipeline code, trained machine learning classifiers (Tox21 mitochondrial safety models), and most importantly, a detailed **Wet-Lab Experimental Validation Protocol** under the MIT license. Our goal is to bridge dry-lab predictions with direct, reproducible laboratory testing.
> 
> Check out the repository, run the simulations, and contribute:
> 🔗 [INSERT YOUR GITHUB REPO LINK HERE]
> 
> #DeSci #OpenScience #ComputationalBiology #Bioinformatics #SMATreatment #Biotech #RareDisease

---

## 2. Reddit r/bioinformatics & r/desci Post Template (Technical & Algoritmik)

**Title:** [Showcase] OpenSMA: An open-source 12-phase pipeline for SMA drug discovery - Generative GA, RNA Nearest-Neighbor Thermodynamics, Metropolis Monte Carlo, and 5-Compartment PBPK ODEs.

**Body:**
> Hi everyone,
> 
> I wanted to share **OpenSMA**, an end-to-end computational drug discovery and systems biology platform we built for Spinal Muscular Atrophy (SMA), now fully open-sourced under the MIT license.
> 
> The repository contains a complete 12-phase pipeline written in Python, integrating physics-based simulations with clinical trial-calibrated pharmacokinetics.
> 
> ### Key Technical Features:
> 1.  **De Novo Ligand Generation (GA):** Uses RDKit to perform scaffold-based mutations on known splicing modulators (like Risdiplam) using specific Reaction SMARTS (fluorination, methylation, nitrogen insertion). Optimization relies on a composite fitness function scoring QED, CNS MPO, and MMFF94 force-field strain energy.
> 2.  **ASO Design & Thermodynamics:** A sliding window scanner that computes hybridization free energy ($\Delta G^\circ_{37}$) of ASOs targeting the ISS-N1 intron 7 region using **Freier 1986 Nearest-Neighbor doublets**.
> 3.  **3D Docking & Metropolis Monte Carlo:** Converts ligands using Meeko, runs AutoDock Vina API on the PDB 8R62 structure (U1 snRNP complex), and feeds the pose into a **Metropolis Monte Carlo simulator** at $310\text{ K}$ for 1,000 steps using Lennard-Jones (12-6) and explicit hydrogen bonding. It filters out candidates that undergo pocket escape (monitored via RMSD drift).
> 4.  **Tox21 ML Model:** A Gradient Boosting Classifier trained on 1,999 compounds from the Tox21 database for Mitochondrial Membrane Potential Disruption (SR-MMP), achieving a cross-validated **ROC-AUC of 0.875** using Morgan Fingerprints (ECFP4) + physical descriptors.
> 5.  **5-Compartment PBPK ODE Engine:** Simulates local drug accumulation in the brain parenchyma and CSF. The PK parameters are fitted using SciPy's curve-fitting algorithms against published clinical data from Finkel et al. (2017) and Baranello et al. (2021) for Nusinersen and Risdiplam (RMSE: 0.33 and 0.28).
> 6.  **Monte Carlo Virtual Trial (10,000 patients):** Explores a "Full Cure" protocol combining a transient ASO bridge with permanent CRISPR Base Editing (CBE) and daily oral modulators. Splicing PD is coupled with a motor neuron decay ODE and mapped to Markov transition probability matrices for clinical milestones (sitting, walking).
> 
> ### Why we are posting:
> We packaged the entire pipeline, datasets, and a highly detailed **in-vitro validation wet-lab protocol** (qPCR Exon 6/8 primers, transfection in GM03813 SMA fibroblasts, Western blot protocols) into the repo.
> 
> We are looking for:
> -   **Wet-lab collaborators** interested in synthesizing and testing our de novo candidates.
> -   **Bioinformaticians** to review our MMC energy functions and expand the pipeline to other splicing-defective diseases like Duchenne Muscular Dystrophy (DMD).
> 
> Check out the repository and the full documentation here:
> 🔗 [INSERT GITHUB LINK HERE]
> 
> Feedback and PRs are highly welcome!
