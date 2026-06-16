# OpenSMA LinkedIn ve Reddit Gönderi Şablonları

Bu dosya, projenizi profesyonel ağlarda ve akademik Reddit topluluklarında paylaşarak ilgi çekmek için özelleştirilmiş uzun format gönderi şablonlarını içerir.

---

## 1. LinkedIn Gönderi Şablonu (Profesyonel & DeSci Odaklı)

**Başlık:** SMA Hastalığına Karşı Açık Kaynaklı İlaç Keşfi ve Klinik Simülasyon: OpenSMA Platformu Yayında!

**Metin:**
> Spinal Müsküler Atrofi (SMA) gibi nadir ve tedavisi milyonlarca doları bulan hastalıklar için pre-klinik ilaç keşif süreçlerini nasıl hızlandırabilir ve demokratikleştirebiliriz?
> 
> Bu soruya yanıt olarak geliştirdiğimiz tamamen açık kaynaklı **OpenSMA** platformunu duyurmaktan mutluluk duyuyorum. Platform, de novo küçük molekül tasarımından başlayıp 5 yıllık klinik prognozu saniyeler içinde simüle eden 12 aşamalı hesaplamalı bir keşif motorudur.
> 
> **OpenSMA Neleri Bir Araya Getiriyor?**
> 🧬 **Smart ASO Tasarımı:** pre-mRNA hibridizasyon kararlılığını Freier 1986 Nearest-Neighbor termodinamik parametreleriyle tahmin eden gelişmiş RNA bağlama analitiği.
> 🧪 **3D Moleküler Yerleşim ve Dinamik:** U1 snRNP-RNA cebine (PDB 8R62) karşı AutoDock Vina yerleşimi ve vücut sıcaklığında ($310\text{ K}$) çalışan 1000 adımlık Metropolis Monte Carlo dinamik stabilitesi.
> 📊 **Kişiselleştirilmiş PBPK Simülasyonu:** 5-bölmeli (Gut, Plasma, Tissue, CSF, Brain) kütle transfer ODE sistemi. Modellerimiz Nusinersen ve Risdiplam Faz 3 klinik verileriyle kalibre edilmiştir (Uyum RMSE: $0.33$ ve $0.28$).
> 📈 **Full Cure Protokolü:** Kalıcı CRISPR Base Editing öncesinde motor nöron dejenerasyonunu durduracak "Bridging ASO" ve oral modülatör içeren kombinasyon terapisinin 10.000 sanal hasta üzerindeki uzun vadeli Monte Carlo klinik sonuçları.
> 
> **Neden Açık Kaynak ve DeSci (Merkeziyetsiz Bilim)?**
> İlaç keşif süreçlerinin kapalı kapılar ardında yapılması ve yüksek maliyetler, nadir hastalık araştırma hızını yavaşlatıyor. OpenSMA ile kodlarımızı, veri setlerimizi ve fibroblast hücre hatlarında denenebilecek detaylı **Wet-Lab Islak Laboratuvar Protokolünü** MIT lisansı ile paylaşıyoruz. Amacımız, deneysel biyologlar ve hesaplamalı bilimciler arasında köprü kurmaktır.
> 
> Bizi GitHub üzerinden takip edebilir, projeye katkı sunabilir veya laboratuvarınızda test edebilirsiniz:
> 🔗 [GITHUB REPO LINKINIZ]
> 
> #DeSci #OpenScience #ComputationalBiology #Bioinformatics #SMATreatment #Biotech

---

## 2. Reddit r/bioinformatics Gönderi Şablonu (Teknik & Algoritmik)

**Başlık:** [Showcase] OpenSMA: An open-source 12-phase pipeline for SMA drug discovery - Generative GA, RNA Nearest-Neighbor Thermodynamics, Metropolis Monte Carlo, and 5-Compartment PBPK ODEs.

**Metin:**
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
> 🔗 [LINK TO GITHUB REPO]
> 
> Feedback and PRs are highly welcome!
