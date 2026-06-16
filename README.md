# OpenSMA: Open-Source (DeSci) Spinal Muscular Atrophy Drug Discovery and Simulation Platform

OpenSMA is an open-source computational medicine and systems biology project focused on **Decentralized Science (DeSci)**. It is designed to screen and optimize novel therapeutic candidates (small molecules, ASOs, and CRISPR gRNAs) for Spinal Muscular Atrophy (SMA) and simulate their longitudinal pharmacokinetic and clinical response profiles.

---

## Repository Structure

The files in the repository are organized as follows:

```text
OpenSMA_OpenSource_Release/
├── LICENSE                     # MIT Open Source License
├── requirements.txt            # Python dependencies list
├── README.md                   # This overview file
├── data/                       # 3D Protein and genetic inputs
│   ├── smn2_receptor_real.pdb  # Raw PDB structure of the SMN2 splicing complex
│   ├── receptor_clean.pdb      # Cleaned 3D structure for docking and MMC
│   └── smn2_target_region.fasta # Target mRNA sequence for ASO scanning
├── docs/                       # Scientific monographs and lab protocols
│   ├── OpenSMA_Lab_Protocol.md # Wet-lab validation protocol
│   ├── OpenSMA_Master_Technical_Report.md # Full pipeline mathematical/biophysical report
│   ├── OpenSMA_Academic_Monograph.md # Academic monograph on SMA pathology and physics equations
│   ├── OpenSMA_Clinical_PBPK_Math.md # PBPK ODE derivation guide
│   ├── OpenSMA_Validation_Dossier.md # Machine learning toxicity and genomic off-target validation
│   ├── Production_and_Cost_Analysis.md # Manufacturing cost analysis report
│   ├── Project_Architecture_and_Logic.md # Pipeline workflow overview
│   ├── OpenSMA_Manuscript.md # Draft research manuscript
│   └── White_Paper_OpenSMA.md # OpenSMA executive white paper
├── src/                        # 12-Phase pipeline source code
│   ├── run_opensma_pipeline.py # Main end-to-end pipeline executor
│   ├── molecule_optimizer.py   # De novo small molecule generator (Genetic Algorithm)
│   ├── aso_designer.py         # Thermodynamic ASO design module
│   ├── crispr_designer.py      # CRISPR/CBE gRNA design module
│   ├── admet_screener.py       # ADMET filtration module
│   ├── aso_toxicity.py         # ASO off-target toxicity scanner
│   ├── docking_engine.py       # Open3DAlign structure alignment and Vina docking engine
│   ├── toxicity_model.py       # Tox21 ML mitochondrial toxicity classifier
│   ├── ode_calibration_and_screening.py # Literature-calibrated PBPK ODE engine
│   ├── final_scorer.py         # Multi-criteria scoring and candidate ranker
│   ├── patient_sim.py          # Multi-compartment PBPK patient simulator
│   ├── full_cure_sim.py        # 10,000 patient Monte Carlo virtual trial simulation
│   └── pocket_dynamics.py      # Metropolis Monte Carlo 3D pocket stability simulator
├── promotional_kit/            # Promotion and community sharing kit
│   ├── community_sharing_guide.md # Roadmap for sharing on DeSci and social networks
│   ├── desci_pitch_templates.md # Pitch copy for Discord and ResearchHub uploads
│   ├── linkedin_and_reddit_posts.md # Custom text templates for LinkedIn and Reddit
│   └── twitter_thread_templates.md # 6-tweet scientific communication thread template
└── results/                    # Directory where simulation outputs are saved
```

---

## 🛠️ Installation and Setup

### 1. Installing Dependencies
The platform relies on **RDKit** for chemical informatics, **Biopython** for sequence manipulation, and **SciPy**, **Scikit-learn**, **Pandas**, and **Matplotlib** for mathematical modeling, machine learning, and data visualization.

Create a Python virtual environment (virtualenv) and install the dependencies:

```bash
# Create a virtual environment
python3 -m venv sma_env
source sma_env/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### 2. Running the Pipeline
To execute the complete 12-phase drug discovery and clinical simulation pipeline, run:

```bash
python src/run_opensma_pipeline.py
```

The pipeline will:
1. Generate de novo small molecule candidates using the Genetic Algorithm.
2. Scan and design ASO sequences using thermodynamic models.
3. Align 3D conformations and dock compounds into the pre-mRNA pocket.
4. Filter out candidates failing ADMET, off-target, and machine learning safety thresholds.
5. Solve multi-compartment PBPK ODE systems for personalized virtual patients.
6. Run 10,000 patient Monte Carlo cohort simulations and visualize the longitudinal clinical outcomes.

---

## 🔬 Wet-Lab Experimental Validation

The repository includes a detailed wet-lab experimental blueprint: **[OpenSMA_Lab_Protocol.md](docs/OpenSMA_Lab_Protocol.md)**.
This protocol outlines the step-by-step methods (RT-qPCR splicing assays using Exon 6/8 flanking primers, Western Blotting, transfection in GM03813 SMA fibroblasts, and SMI-32 motor neuron survival assays) needed to synthesize and validate the designed candidates in-vitro.

---

## 📜 License and Contribution

This project is licensed under the **MIT License**. Scientists, pharmacologists, and researchers worldwide are free to download, modify, and extend this codebase for non-profit and academic drug discovery research.
