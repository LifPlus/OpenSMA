# OpenSMA Twitter/X Promotion Thread Templates

This file contains optimized Twitter/X thread templates designed to announce the project to the global scientific, bioinformatics, and DeSci communities.

---

## THREAD 1: Scientific & DeSci Focused Thread (Recommended)

### Tweet 1: Introduction (Hook)
> 1/ Spinal Muscular Atrophy (SMA) treatments (Zolgensma, Spinraza) are among the most expensive drugs in the world. Can we democratize this process?
> 
> Today we are open-sourcing **OpenSMA**, an end-to-end computational drug discovery pipeline and systems biology simulator for SMA! 🧵👇
> 
> #DeSci #OpenScience #Bioinformatics #RareDiseases

### Tweet 2: ASO and Thermodynamics
> 2/ The pipeline screens ASOs targeting the pre-mRNA ISS-N1 intron 7 region using the Freier 1986 Nearest-Neighbor model.
> 
> While the FDA-approved Nusinersen control scores a $\Delta G^\circ_{37}$ of $-24.6\text{ kcal/mol}$, our top de novo designed ASO candidate exhibits a thermodynamic binding stability of **$-36.3\text{ kcal/mol}$**! 🧬

### Tweet 3: 3D Docking & Metropolis Monte Carlo
> 3/ We docked de novo small molecules into the PDB 8R62 U1-RNA pocket using AutoDock Vina (best binding energy: $-7.7\text{ kcal/mol}$).
> 
> To test pocket stability under physiological conditions ($310\text{ K}$), we ran 1000-step **Metropolis Monte Carlo** dynamics, tracking RMSD drift. 💻

### Tweet 4: Clinical PBPK ODE Engine
> 4/ Moving beyond basic PK models, we built a **5-compartment PBPK ODE** system simulating drug absorption, blood-brain barrier passage, and accumulation in the CSF and brain parenchyma.
> 
> Calibrated against Phase 3 trial clinical data (Finkel 2017 & Baranello 2021) with RMSE of $0.33$ and $0.28$.

### Tweet 5: 10,000 Patient Virtual Trial
> 5/ We simulated a multi-modal "Full Cure" protocol: bridging transient ASO dosage with permanent CRISPR Base Editing (CBE) genotoxicity checks.
> 
> 5-year Monte Carlo cohort results (10,000 virtual patients):
> 📈 Cumulative survival: `%98.7`
> 🚶 Probability of independent walking: `%50.2`

### Tweet 6: Call to Action & GitHub Repo
> 6/ The entire pipeline code, ML models, PBPK solvers, and a detailed **in-vitro wet-lab experimental validation protocol** (fibroblast transfection, Western Blot, qPCR primers) are MIT licensed.
> 
> Review the code, run the pipeline, and collaborate:
> 🔗 [INSERT YOUR GITHUB REPO LINK HERE]
> 
> #OpenSource #RareDisease #Biotech

---

## Visual Attachments Guide
*   **Attach to Tweet 1 or 5:** The generated `full_cure_simulation.png` plot.
*   **Attach to Tweet 4:** The multi-compartment `patient_simulation_plots.png` plot.
