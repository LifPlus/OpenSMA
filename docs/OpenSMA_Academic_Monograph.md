# OpenSMA: Spinal Muscular Atrophy (SMA) Therapeutic Discovery & Systems Biology Simulation Academic Monograph

This scientific monograph serves as the primary academic reference document containing all biochemical, biophysical, and pharmacokinetic theories, equations, and literature sources underlying the OpenSMA platform.

---

## 1. Genetics and Splicing Pathophysiology of Spinal Muscular Atrophy (SMA)

Spinal Muscular Atrophy (SMA) is an autosomal recessive neuromuscular disorder characterized by the degeneration of spinal motor neurons leading to progressive muscle wasting. The molecular etiology of the disease is the homozygous deletion or loss-of-function mutations of the *SMN1* (Survival Motor Neuron 1) gene located on chromosome 5q13.2 (Lefebvre et al., 1995).

### 1.1. SMN1 and SMN2 Gene Structures and Transcription Homology
Due to an evolutionary duplication event, the human genome contains the *SMN2* (Survival Motor Neuron 2) gene, an almost identical copy located adjacent to the *SMN1* gene. These two genes share 99.9% sequence homology. The most critical functional difference is a single nucleotide transition in Exon 7 (position +6 of the coding region):

$$\text{SMN1 Exon 7 (+6): } 5'\text{-GGT TTC C C GA CAA- } 3' \implies \text{Normal Splicing}$$
$$\text{SMN2 Exon 7 (+6): } 5'\text{-GGT TTC T GA CAA- } 3' \implies \text{Exon Skipping (Delta7 SMN)}$$

Although this C $\to$ T transition is synonymous (both codons encode Glutamic Acid), it severely disrupts the pre-mRNA splicing process.

### 1.2. ESE and ESS Balance: SRSF1 and hnRNP A1/A2 Mechanisms
Splicing is the process where exon-intron boundaries on pre-mRNA are recognized and processed by the spliceosome complex. This mechanism depends on the delicate balance of protein-RNA interactions between exonic splicing enhancers (ESEs) and exonic splicing silencers (ESSs).
1.  **SRSF1 (SF2/ASF) Interaction:** In the *SMN1* gene, the +6 Cytosine (C) in Exon 7 forms the core of an ESE motif. This motif is recognized by the Serine/Arginine-Rich Splicing Factor 1 (SRSF1) protein. SRSF1 binding recruits the U1 snRNP (small nuclear ribonucleoprotein) complex to the 5' splice site (5'SS), directing Exon 7 inclusion into the mature mRNA.
2.  **hnRNP A1 Repression:** In the *SMN2* gene, the +6 Uracil (U) dramatically decreases the binding affinity of SRSF1. Simultaneously, it facilitates the cooperative binding of **hnRNP A1** (Heterogeneous Nuclear Ribonucleoprotein A1) to the adjacent **ISS-N1** (Intronic Splicing Silencer N1) region in Intron 7. hnRNP A1 induces pre-mRNA looping or sterically blocks spliceosome access, leading to Exon 7 exclusion (skipping).

Consequently, ~90% of *SMN2* transcription results in an unstable, rapidly degraded protein isoform lacking the C-terminus ($\Delta7$ SMN). Only ~10% results in functional full-length (FL) SMN protein, which is insufficient for motor neuron survival.

---

## 2. Scaffold-Based Genetic Algorithm (GA) Optimization of Small Molecules

The in-silico design of small molecule splicing modulators (e.g., Risdiplam derivatives) is executed using a directed Genetic Algorithm (GA) loop.

### 2.1. SMILES Representation and SMARTS Mutation Rules
Molecules are represented as Simplified Molecular Input Line Entry System (SMILES) strings. Mutation operations are performed using RDKit with predefined SMARTS reaction rules (Reaction SMARTS) to ensure chemical synthesizability and stability.

The mutation reaction matrix is defined below:

| Mutation Type | Reaction SMARTS Rule | Chemical Effect |
| :--- | :--- | :--- |
| **Fluorination** | `[c:1][H]>>[c:1]F` | Addition of Fluorine to aromatic carbons (Metabolic blockade) |
| **Methylation** | `[c:1][H]>>[c:1]C` | Increases lipophilicity and hydrophobic pocket fit |
| **Methoxylation** | `[c:1][H]>>[c:1]OC` | Modifies TPSA and hydrogen bond acceptor parameters |
| **Chlorination** | `[c:1][H]>>[c:1]Cl` | Enables potential halogen bonding interactions |
| **Aromatic Nitrogen Swap** | `[c:1][H]>>[n:1]` | Adjusts ring basicity and polarity |

### 2.2. Multi-Parameter Optimization (MPO) Fitness Function
The genetic algorithm is guided by a composite fitness function:

$$\text{Fitness} = (\text{QED} \times 100) + S_{\text{BBB\_bonus}} - S_{\text{penalties}} + S_{\text{3D\_stability}}$$

#### QED (Quantitative Estimate of Druglikeness) Calculation
Based on the Bickerton et al. (2012) (*Nature Chemistry*) model, QED is computed as the geometric mean of individual desirability functions ($d_i$) for eight key physicochemical properties (MW, LogP, TPSA, HBD, HBA, RotBonds, Aromatic Rings, Alerts):

$$\text{QED} = \left( \prod_{i=1}^{8} d_i(x_i) \right)^{\frac{1}{8}}$$

#### BBB (Blood-Brain Barrier) Penetration Criteria
For SMA therapeutics, CNS penetration is mandatory:
-   $2.0 \le \text{LogP} \le 5.0 \implies S_{\text{BBB\_bonus}} = S_{\text{BBB\_bonus}} + 20$
-   $\text{TPSA} < 90 \text{ Å}^2 \implies S_{\text{BBB\_bonus}} = S_{\text{BBB\_bonus}} + 20$
-   $\text{MW} < 500 \text{ Da} \implies S_{\text{BBB\_bonus}} = S_{\text{BBB\_bonus}} + 10$

#### Physicochemical Penalties
-   $\text{HBD} > 3 \implies S_{\text{penalties}} = S_{\text{penalties}} + 10$
-   $\text{HBA} > 7 \implies S_{\text{penalties}} = S_{\text{penalties}} + 10$
-   $\text{LogP} > 5.5 \implies S_{\text{penalties}} = S_{\text{penalties}} + 20$
-   $\text{MW} > 600 \text{ Da} \implies S_{\text{penalties}} = S_{\text{penalties}} + 30$

---

## 3. ASO Hybridization Thermodynamics and the Nearest-Neighbor Model

The physical binding stability of designed ASOs targeting the pre-mRNA ISS-N1 region is determined by the standard free energy of hybridization ($\Delta G^\circ_{37}$).

### 3.1. Freier 1986 Parameters Table
Free energy calculations utilize the RNA-RNA Nearest-Neighbor model defined by Freier et al. (1986) (*PNAS*). The standard free energy contributions ($\Delta G^\circ_{37}$, kcal/mol) of the 16 doublets at $37^\circ\text{C}$ are:

| Doublet (5' $\to$ 3') | $\Delta G^\circ_{37}$ (kcal/mol) | Doublet (5' $\to$ 3') | $\Delta G^\circ_{37}$ (kcal/mol) |
| :--- | :---: | :--- | :---: |
| **AA / UU** | -0.9 | **GA / CU** | -2.1 |
| **AU / UA** | -0.9 | **GG / CC** | -3.1 |
| **UA / AU** | -1.1 | **CG / GC** | -2.0 |
| **CA / GU** | -1.8 | **GC / CG** | -3.4 |
| **CC / GG** | -3.1 | **UG / AC** | -2.1 |
| **CU / GA** | -1.7 | **UU / AA** | -0.9 |

### 3.2. Total Free Energy and Initiation/Terminal Penalties
The total free energy of hybridization is computed by summing the doublet values and adding helix initiation and terminal AU pair penalties:

$$\Delta G^\circ_{37\text{, total}} = \Delta G^\circ_{\text{init}} + \sum_{i=1}^{L-1} \Delta G^\circ_{\text{doublet}}(i, i+1) + N_{\text{term\_AU}} \cdot \Delta G^\circ_{\text{term\_AU}}$$

Where:
-   $\Delta G^\circ_{\text{init}} = +3.4 \text{ kcal/mol}$ (Helix nucleation energy)
-   $\Delta G^\circ_{\text{term\_AU}} = +0.5 \text{ kcal/mol}$ (Terminal AU pair instability penalty)
-   $N_{\text{term\_AU}}$: Number of terminal A or U bases at the ends of the helix ($0, 1$, or $2$).

---

## 4. 3D Molecular Docking and Metropolis Monte Carlo Simulation

The binding and physical stability of candidate molecules within the experimental 3D NMR structure of the SMN2 pre-mRNA-U1 snRNP complex (PDB 8R62) are simulated at body temperature.

### 4.1. AutoDock Vina Scoring Function
AutoDock Vina utilizes an empirical, knowledge-based free energy scoring function:

$$c(t) = h(d(t)) \quad \text{where} \quad d(t) = r_i - r_j$$

Free energy of binding score ($\Delta G_{\text{bind}}$):

$$\Delta G_{\text{bind}} = w_1 \cdot \text{gauss}_1(d) + w_2 \cdot \text{gauss}_2(d) + w_3 \cdot \text{repulsion}(d) + w_4 \cdot \text{hydrophobic}(d) + w_5 \cdot \text{hbond}(d) + w_{\text{rot}} \cdot N_{\text{rot}}$$

Where:
-   **$\text{gauss}_1(d) = e^{-(d/0.5)^2}$** (Short-range attraction)
-   **$\text{gauss}_2(d) = e^{-((d-3.0)/2.0)^2}$** (Long-range attraction)
-   **$\text{repulsion}(d) = d^2$** (for $d < 0$; sterics repulsion)
-   **$\text{hydrophobic}(d)$:** Hydrophobic contact score between hydrophobic atoms.
-   **$\text{hbond}(d)$:** Hydrogen bond donor-acceptor geometric term.
-   **$N_{\text{rot}}$:** Number of active rotatable bonds (conformational entropy loss penalty).

### 4.2. Metropolis Monte Carlo Dynamics and Thermal Fluctuations
The mobility and stability of the ligand in the pocket are verified at body temperature ($T = 310\text{ K}$) using the Metropolis Monte Carlo (MMC) algorithm.
1.  **Potential Energy ($V$):** Computed using Lennard-Jones (12-6) potentials and hydrogen bonding attraction terms:
    $$V_{\text{LJ}} = \sum_{i \in \text{ligand}} \sum_{j \in \text{pocket}} 4\epsilon \left[ \left(\frac{\sigma}{r_{ij}}\right)^{12} - \left(\frac{\sigma}{r_{ij}}\right)^6 \right]$$
    Where $\sigma = 3.2\text{ Å}$ and $\epsilon = 0.15\text{ kcal/mol}$. A steric clash penalty of $+150\text{ kcal/mol}$ is applied if any two atoms are closer than $1.2\text{ Å}$.
2.  **Perturbation:** In each step, the ligand is subjected to a random translation ($\Delta x \in [-0.06, 0.06]\text{ Å}$) and three-axis rotation ($\Delta \theta \in [-0.025, 0.025]\text{ rad}$).
3.  **Metropolis Acceptance Criterion:** The energy change $\Delta E = E_{\text{new}} - E_{\text{old}}$ is calculated.
    -   If $\Delta E \le 0$, the move is **accepted** ($P_{\text{accept}} = 1.0$).
    -   If $\Delta E > 0$, the move is accepted with Boltzmann probability:
        $$P_{\text{accept}} = e^{-\frac{\Delta E}{k_B T}}$$
        Where $k_B = 0.001987\text{ kcal/(mol}\cdot\text{K)}$ and $T = 310\text{ K} \implies k_B T \approx 0.616\text{ kcal/mol}$.
4.  **RMSD (Root-Mean-Square Deviation) Tracking:** The displacement of the ligand from its docked starting position is tracked:
    $$\text{RMSD} = \sqrt{\frac{1}{N} \sum_{i=1}^{N} (x_i - x_{i,0})^2 + (y_i - y_{i,0})^2 + (z_i - z_{i,0})^2}$$
    Molecules maintaining an RMSD $< 2.2\text{ Å}$ after 1,000 steps are classified as *STABLE*, while those exceeding $4.0\text{ Å}$ are classified as *UNSTABLE* due to pocket escape.

---

## 5. Reference Scientific Papers and Literature Sources

1.  **Lefebvre, S., et al. (1995).** "Identification and characterization of a spinal muscular atrophy-determining gene." *Cell*, 80(1), 155-165.
2.  **Freier, S. M., et al. (1986).** "Improved free-energy parameters for predictions of RNA duplex stability." *Proceedings of the National Academy of Sciences*, 83(24), 9373-9377.
3.  **Wager, T. T., et al. (2010).** "Moving beyond Rules: The development of a central nervous system multiparameter optimization (CNS MPO) approach to enable alignment of druggability parameters for rational drug design." *ACS Chemical Neuroscience*, 1(6), 435-449.
4.  **Bickerton, G. R., et al. (2012).** "Quantifying the chemical beauty of drugs." *Nature Chemistry*, 4(2), 90-98.
5.  **Tran, Q. D., et al. (2020).** "Structural basis for targeted splicing modulation by small molecules." *Nature*, 581(7806), 105-110. (PDB 8R62 Modeling Reference).
6.  **Finkel, R. S., et al. (2017).** "Nusinersen versus Sham Control in Infantile-Onset Spinal Muscular Atrophy." *New England Journal of Medicine*, 377(18), 1723-1732. (ENDEAVOR/NURTURE Clinical Trial Reference).
7.  **Baranello, G., et al. (2021).** "Risdiplam in Type 1 Spinal Muscular Atrophy." *New England Journal of Medicine*, 384(10), 915-923. (FIREFISH Clinical Trial Reference).
