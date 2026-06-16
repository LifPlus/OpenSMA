# OpenSMA: Open-Source (DeSci) Wet-Lab Validation Protocol

This document serves as the step-by-step experimental validation protocol describing how to test, formulate, and validate the small molecule and ASO (Antisense Oligonucleotide) candidates designed by the OpenSMA platform in-vitro and in-vivo.

---

## 1. Experimental Design and Starting Materials

### 1.1. Cell Lines
*   **SMA Patient Fibroblast Cell Line:** GM03813 (Coriell Institute). Genotype: *SMN1* homozygous deletion, *SMN2* 2 copies (Type 1 SMA clinical phenotype).
*   **Control Cell Line:** GM03815 (Healthy carrier parent fibroblasts).
*   **General Expression Control:** HEK293T (ATCC CRL-11268).

### 1.2. Reagents and Chemicals
*   **Culture Medium:** DMEM (Dulbecco's Modified Eagle Medium - Gibco, Cat #11965092), 10% Fetal Bovine Serum (FBS - Gibco), 1% Penicillin-Streptomycin.
*   **Transfection Reagent:** Lipofectamine 3000 (Invitrogen) or MC3-based Lipid Nanoparticles (LNPs).
*   **ASO Chemistry:** 2'-O-methoxyethyl (2'-MOE) modified phosphorothioate (PS) backbone oligonucleotides (obtained from custom synthesis suppliers at PAGE or HPLC purity).
*   **Primary Antibodies:** Mouse anti-SMN (Clone 2B1, Sigma-Aldrich, Cat #S2944), Rabbit anti-Beta-Actin (Cell Signaling Technology, Cat #4970).
*   **SMI-32 Antibody:** Mouse anti-Neurofilament H Non-Phosphorylated (SMI-32, BioLegend, Cat #801701 - Motor neuron-specific marker).

---

## 2. Protocol A: ASO and LNP Transfection (In Vitro Splicing Modulation)

### 2.1. LNP Formulation (Microfluidic Chip Method)
1.  **Lipid Preparation:** Weigh the ionizable lipid DLin-MC3-DMA, DSPC, Cholesterol, and DMG-PEG2000 at a molar ratio of 50:10:38.5:1.5, respectively. Dissolve in anhydrous 100% Ethanol to obtain a total lipid concentration of 10 mM.
2.  **ASO Preparation:** Dissolve the synthesized ASO candidate in 50 mM Sodium Citrate buffer (pH 4.0) to achieve a 0.2 mg/mL nucleic acid solution.
3.  **Mixing (Formulation):** Impinge the aqueous phase and the organic ethanol phase in a microfluidic mixer at a flow rate ratio of 3:1 (Total flow rate = 12 mL/min).
4.  **Dialysis:** Dialyze the mixture against PBS (pH 7.4) buffer using a Slide-A-Lyzer dialysis cassette (10k MWCO) at 4°C for 16 hours to completely remove the ethanol.
5.  **Size Analysis:** Verify LNP size (target: 70-90 nm) and polydispersity index (PDI < 0.15) using a Dynamic Light Scattering (DLS) instrument.

### 2.2. Transfection of Fibroblast Cells
1.  Seed GM03813 patient fibroblasts in 6-well plates at $2 \times 10^5$ cells per well and incubate for 24 hours at 37°C in a 5% $CO_2$ incubator.
2.  Once cells reach 70-80% confluence, replace the growth medium with serum-free Opti-MEM (Gibco).
3.  Add ASO-LNP formulations (or naked ASOs complexed with Lipofectamine 3000) to wells at various concentrations (10 n M, 50 nM, 100 nM, 200 nM). Add empty LNPs (with Scrambled ASO) to control wells.
4.  After 4-6 hours, replace the transfection medium with standard DMEM containing 10% FBS. Incubate cells for 48 hours before analysis.

---

## 3. Protocol B: Analysis of Exon 7 Inclusion Rate via RT-qPCR

### 3.1. Total RNA Isolation and cDNA Synthesis
1.  Isolate total RNA from transfected cells using TRIzol (Invitrogen) reagent via standard chloroform phase separation.
2.  Verify RNA concentration and purity on a Nanodrop spectrophotometer ($A_{260}/A_{280} \ge 1.8$).
3.  Convert 1 µg of total RNA into cDNA using the High-Capacity cDNA Reverse Transcription Kit (Applied Biosystems) with random hexamer primers.

### 3.2. Semi-Quantitative PCR and qPCR Conditions
*   **RT-PCR Primers Determining Exon 7 Inclusion (Human SMN2 specific):**
    *   **Forward (Exon 6):** 5'-GCTATCATAATTTTGTTCATTTT-3'
    *   **Reverse (Exon 8):** 5'-CCAGATTCTCTTGATGATGC-3'
*   **PCR Reaction Mix:** 2 µL cDNA, 1.25 U Taq DNA Polymerase, 0.4 µM of each primer, 0.2 mM dNTPs, PCR Buffer (Total volume: 25 µL).
*   **PCR Thermal Cycler Protocol:**
    1.  Initial Denaturation: 95°C for 3 min
    2.  Cycling (30 repeats):
        *   Denaturation: 95°C for 30 sec
        *   Annealing: 55°C for 30 sec
        *   Extension: 72°C for 45 sec
    3.  Final Extension: 72°C for 5 min
4.  **Gel Electrophoresis:** Run PCR products on a 2% agarose gel in the presence of ethidium bromide.
    *   **Expected Band Sizes:** Full-length SMN2 transcript (incorporating Exon 7) = **263 bp**, skipped transcript ($\Delta7$-SMN2) = **209 bp**.
5.  Analyze band intensities using ImageJ software. Calculate the Exon 7 inclusion percentage using the following equation:
    $$\text{Exon 7 Inclusion \%} = \frac{\text{Intensity (263 bp Band)}}{\text{Intensity (263 bp Band)} + \text{Intensity (209 bp Band)}} \times 100$$

---

## 4. Protocol C: Quantification of SMN Protein Levels via Western Blot

### 4.1. Protein Isolation
1.  Lyse cells on ice using RIPA Lysis Buffer containing protease and phosphatase inhibitor cocktails.
2.  Centrifuge lysates at 14,000 rpm for 15 min at 4°C and collect the supernatant.
3.  Determine protein concentration using the BCA Protein Assay Kit (Thermo Fisher).

### 4.2. SDS-PAGE and Membrane Transfer
1.  Denature 30 µg of protein samples with Laemmli loading buffer at 95°C for 5 min.
2.  Load samples onto a 10% SDS-PAGE gel and run at 100V.
3.  Transfer proteins from the gel to a PVDF membrane (Bio-Rad) using the Wet Transfer method at 300 mA for 1.5 hours.

### 4.3. Antibody Incubation and Imaging
1.  Block the membrane in TBST (Tris-buffered saline, 0.1% Tween-20) containing 5% non-fat dry milk powder for 1 hour at room temperature.
2.  Incubate with primary antibodies (Mouse anti-SMN, 1:1000 and Rabbit anti-Beta-Actin, 1:5000) for 16 hours at 4°C.
3.  Wash the membrane 3 times for 10 min each with TBST.
4.  Incubate with HRP-conjugated secondary antibodies (Anti-mouse IgG-HRP and Anti-rabbit IgG-HRP, 1:10,000) for 1 hour at room temperature.
5.  Visualize bands using ECL Western Blotting Substrate (Pierce) on a chemiluminescence imaging system (ChemiDoc). Calculate the fold increase of SMN protein normalized against the Beta-Actin reference.

---

## 5. Protocol D: Motor Neuron Survival and SMI-32 Staining Analysis

Test the neuroprotective capability of the small molecules or ASOs in co-culture models:
1.  **Cell Isolation:** Isolate spinal cord motor neurons from embryonic day 14 (E14) mice and culture in Neurobasal Medium containing 2% B-27, 0.5 mM L-Glutamine, and 10 ng/mL BDNF, GDNF, and CNTF.
2.  **SMA Phenotype Induction:** Trigger degeneration by transfecting motor neurons with siRNA or exposing them to low-dose toxins, or adding conditioned medium obtained from patient fibroblasts.
3.  **Drug Treatment:** Add candidate small molecules (optimal dose range: 10 nM - 10 µM) or ASOs to the cells.
4.  **Immunofluorescence Staining (SMI-32):**
    *   Fix cells after 48 hours using 4% Paraformaldehyde (PFA).
    *   Permeabilize using Triton X-100 (0.1%).
    *   Incubate with **SMI-32 Primary Antibody** (1:1000) overnight at 4°C.
    *   Incubate with Alexa Fluor 488 conjugated secondary antibody (1:500) for 1 hour at room temperature.
    *   Stain nuclei using DAPI.
5.  **Survival Counting:** Count healthy SMI-32 positive motor neurons with long axons under a fluorescence microscope. Compare survival rates between untreated degenerative controls and rescued groups.

---

## 6. Protocol E: In Vivo Validation (SMA Delta7 Mouse Model)

### 6.1. Animal Model
*   **Model:** $Smn^{-/-}; SMN2^{+/+}; SMN\Delta7^{+/+}$ (Jackson Laboratory, Strain #005025). These mice recapitulate human Type 1 SMA, surviving ~13-15 days if untreated.

### 6.2. Dosing and Routes of Administration
*   **Small Molecules (Oral):** Starting from postnatal day 1 (P1), administer daily oral gavage doses of 1-10 mg/kg using a micro-syringe.
*   **ASO Treatment (Intrathecal Injection):**
    *   Perform intrathecal injection on postnatal day 2 (P2) under cryogenic anesthesia (light freezing on ice).
    *   Inject 2 µL of LNP-ASO formulation (dose: 10-20 µg) at the L3-L4 level using a 33G Hamilton micro-syringe.

### 6.3. Clinical Parameter Monitoring
1.  **Survival Analysis:** Record natural death days and construct Kaplan-Meier survival curves.
2.  **Body Weight Monitoring:** Weigh mice daily from P1 onwards to record development curves.
3.  **Motor Function Tests (Righting Reflex):** Place the mouse on its back. Measure the time in seconds it takes to flip back onto all four paws (Righting time). Healthy mice right themselves in <2 seconds. SMA mice show progressive delay or loss of righting ability.

---

## 7. Protocol F: Chemical Purity Verification via HPLC and LC-MS

Verify the purity of synthesized small molecules before carrying out in-vitro tests:
1.  **HPLC Conditions:**
    *   Column: C18 Reverse-Phase Column (4.6 mm x 150 mm, 5 µm).
    *   Mobile Phase A: Ultrapure water containing 0.1% Trifluoroacetic acid (TFA).
    *   Mobile Phase B: Acetonitrile containing 0.1% TFA.
    *   Gradient Run: Transition from 5% B to 95% B over 20 minutes; Flow Rate: 1.0 mL/min.
    *   Detector: UV (254 nm and 280 nm).
    *   **Acceptance Criterion:** The ratio of the main peak area to the total peak area must be $\ge 98.5\%$.
2.  **LC-MS Mass Analysis:** Verify that the mass-to-charge ($m/z$) ratio obtained via electrospray ionization mass spectrometry (ESI-MS) matches the theoretical molecular mass of the designed chemical formula.
