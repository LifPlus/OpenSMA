# OpenSMA: Spinal Müsküler Atrofi (SMA) Tedavi Keşfi ve Biyolojik Simülasyonu Akademik Monografı

Bu bilimsel monograf, OpenSMA platformunun temelini oluşturan tüm biyokimyasal, biyofiziksel ve farmakokinetik teorileri, denklemleri ve literatür kaynaklarını içeren ana akademik referans dokümanıdır.

---

## 1. Spinal Müsküler Atrofi (SMA) Genetiği ve Splicing Patofizyolojisi

Spinal Müsküler Atrofi (SMA), motor nöronların kaybı ve buna bağlı ilerleyici kas erimesiyle seyreden otozomal resesif geçişli nöromüsküler bir hastalıktır. Hastalığın moleküler temeli, kromozom 5q13.2 locus'unda yer alan *SMN1* (Survival Motor Neuron 1) geninin delesyonu veya işlev kaybına neden olan nokta mutasyonlarıdır (Lefebvre et al., 1995).

### 1.1. SMN1 ve SMN2 Gen Yapıları ve Transkripsiyon Homolojisi
İnsan genomu, evrimsel bir duplikasyon sonucunda *SMN1* geninin hemen yanında yer alan ve neredeyse özdeş bir kopya olan *SMN2* (Survival Motor Neuron 2) genine sahiptir. Bu iki gen dizilim düzeyinde %99.9 oranında homolojiktir. Aralarındaki en kritik fonksiyonel fark, Ekzon 7'nin 6. nükleotidinde (kodlama bölgesinde +6. pozisyon) bulunan tek bir nükleotid transisyonudur:

$$\text{SMN1 Exon 7 (+6): } 5'\text{-GGT TTC C C GA CAA- } 3' \implies \text{Normal Splicing}$$
$$\text{SMN2 Exon 7 (+6): } 5'\text{-GGT TTC T GA CAA- } 3' \implies \text{Exon Skipping (Delta7 SMN)}$$

Bu C $\to$ T değişimi, kodlanan amino asidi değiştirmemesine rağmen (her iki kodon da Glutamik Asit kodlar), splicing mekanizmasını tamamen bozar.

### 1.2. ESE ve ESS Dengesi: SRSF1 ve hnRNP A1/A2 Mekanizmaları
Splicing, pre-mRNA üzerindeki ekzon-intron sınırlarının spliceozom kompleksi tarafından tanınması işlemidir. Bu işlem, ekzonik splicing güçlendiriciler (ESE) ve ekzonik splicing susturucular (ESS) arasındaki hassas protein-RNA etkileşim dengesine bağlıdır.
1.  **SRSF1 (SF2/ASF) Etkileşimi:** *SMN1* geninde, Ekzon 7'deki +6. sitozin (C), bir ESE motifinin merkezini oluşturur. Bu motif, Serin/Arjinin Zengin Splicing Faktörü 1 (SRSF1) proteini tarafından tanınır. SRSF1 bağlanması, U1 snRNP (küçük nükleer ribonükleoprotein) kompleksini 5' splice sitesine (5'SS) yönlendirerek Ekzon 7'nin olgun mRNA'ya dahil edilmesini sağlar.
2.  **hnRNP A1 Baskılaması:** *SMN2* geninde +6. pozisyondaki Urasil (T), SRSF1'in bağlanma afinitesini dramatik şekilde düşürür. Eş zamanlı olarak, Intron 7'de yer alan ve **ISS-N1** (Intronic Splicing Silencer N1) adı verilen susturucu bölge, **hnRNP A1** (Heterogeneous Nuclear Ribonucleoprotein A1) proteininin kooperatif olarak bağlanmasını kolaylaştırır. hnRNP A1, pre-mRNA zincirini bükerek (looping) veya spliceozom erişimini sterik olarak engelleyerek Ekzon 7'nin splicing dışı bırakılmasına (skipping) yol açar.

Bunun sonucunda, *SMN2* transkripsiyonu %90 oranında kararsız, hızlıca degrade olan ve C-terminali eksik bir protein izoforuna ($\Delta7$ SMN) yol açar. Yalnızca %10 oranında tam boy (Full-Length, FL) fonksiyonel SMN proteini üretilebilir. Bu oran motor nöronların hayatta kalması için yetersizdir.

---

## 2. Küçük Molekül Adaylarının Genetik Algoritma (GA) Tabanlı Optimizasyonu

Küçük moleküllü splicing modülatörlerinin (örneğin Risdiplam türevleri) in-silico tasarımı, yönlendirilmiş bir genetik algoritma (GA) döngüsü kullanılarak gerçekleştirilir.

### 2.1. SMILES Temsili ve SMARTS Mutasyon Kuralları
Moleküller, Simplifed Molecular Input Line Entry System (SMILES) formatında temsil edilir. Mutasyon operasyonları, kimyasal olarak sentezlenebilirliği ve kararlılığı korumak amacıyla önceden tanımlanmış SMARTS reaksiyon kuralları (Reaction SMARTS) kullanılarak RDKit motoruyla gerçekleştirilir.

Mutasyonların matematiksel reaksiyon matrisi:

| Mutasyon Tipi | SMARTS Reaksiyon Kuralı | Kimyasal Etkisi |
| :--- | :--- | :--- |
| **Florlama** | `[c:1][H]>>[c:1]F` | Aromatik halkaya Flor eklenmesi (Metabolik blokaj) |
| **Metilasyon** | `[c:1][H]>>[c:1]C` | Lipofilisite ve hidrofobik cep uyumu artışı |
| **Metoksile Etme** | `[c:1][H]>>[c:1]OC` | TPSA ve hidrojen bağı akseptör modifikasyonu |
| **Klorlama** | `[c:1][H]>>[c:1]Cl` | Halojen bağ etkileşimi kurabilme yeteneği |
| **Halka İçi Azot Değişimi** | `[c:1][H]>>[n:1]` | Halka bazikliğinin ve polaritesinin ayarlanması |

### 2.2. Çok Parametreli Optimizasyon (MPO) Fitness Fonksiyonu
Genetik algoritmanın yönlendirilmesi, kompozit bir fitness fonksiyonu aracılığıyla yapılır:

$$\text{Fitness} = (\text{QED} \times 100) + S_{\text{BBB\_bonus}} - S_{\text{penalties}} + S_{\text{3D\_stability}}$$

#### QED (Quantitative Estimate of Druglikeness) Hesabı
Bickerton ve ark. (2012) (*Nature Chemistry*) modeline dayanan QED, sekiz temel fizikokimyasal özelliğin (MW, LogP, TPSA, HBD, HBA, RotBonds, Aromatic Rings, Alerts) arzu edilirlik fonksiyonlarının (desirability functions) geometrik ortalaması alınarak hesaplanır:

$$\text{QED} = \left( \prod_{i=1}^{8} d_i(x_i) \right)^{\frac{1}{8}}$$

#### BBB (Kan-Beyin Bariyeri) Geçiş Kriterleri
SMA bir nöromüsküler hastalık olduğundan, molekülün merkezi sinir sistemine (CNS) geçmesi zorunludur.
-   $2.0 \le \text{LogP} \le 5.0 \implies S_{\text{BBB\_bonus}} = S_{\text{BBB\_bonus}} + 20$
-   $\text{TPSA} < 90 \text{ Å}^2 \implies S_{\text{BBB\_bonus}} = S_{\text{BBB\_bonus}} + 20$
-   $\text{MW} < 500 \text{ Da} \implies S_{\text{BBB\_bonus}} = S_{\text{BBB\_bonus}} + 10$

#### Fizikokimyasal Cezalar (Penalties)
-   $\text{HBD} > 3 \implies S_{\text{penalties}} = S_{\text{penalties}} + 10$
-   $\text{HBA} > 7 \implies S_{\text{penalties}} = S_{\text{penalties}} + 10$
-   $\text{LogP} > 5.5 \implies S_{\text{penalties}} = S_{\text{penalties}} + 20$
-   $\text{MW} > 600 \text{ Da} \implies S_{\text{penalties}} = S_{\text{penalties}} + 30$

---

## 3. ASO Hibridizasyon Termodinamiği ve Nearest-Neighbor Modeli

ASO'ların hedef pre-mRNA ISS-N1 bölgesine fiziksel bağlanma kararlılığı, deneysel hibridizasyon serbest enerjisi ($dG^\circ_{37}$) ile belirlenir.

### 3.1. Freier 1986 Parametreleri Tablosu
Serbest enerji hesabı, Freier ve ark. (1986) (*PNAS*) tarafından tanımlanmış RNA-RNA Nearest-Neighbor (Komşu Çiftler) modeline dayanır. İki komşu nükleotid çiftinin $37^\circ\text{C}$'deki standart serbest enerji ($\Delta G^\circ_{37}$, kcal/mol) katkıları şunlardır:

| Çift (5' $\to$ 3') | $\Delta G^\circ_{37}$ (kcal/mol) | Çift (5' $\to$ 3') | $\Delta G^\circ_{37}$ (kcal/mol) |
| :--- | :---: | :--- | :---: |
| **AA / UU** | -0.9 | **GA / CU** | -2.1 |
| **AU / UA** | -0.9 | **GG / CC** | -3.1 |
| **UA / AU** | -1.1 | **CG / GC** | -2.0 |
| **CA / GU** | -1.8 | **GC / CG** | -3.4 |
| **CC / GG** | -3.1 | **UG / AC** | -2.1 |
| **CU / GA** | -1.7 | **UU / AA** | -0.9 |

### 3.2. Toplam Serbest Enerji ve Başlatma/Terminal Cezaları
Hibridizasyon toplam serbest enerjisi, dizideki tüm komşu doubletlerin enerjileri toplanarak ve heliks başlatma (initiation) ve terminal AU çifti cezaları eklenerek bulunur:

$$\Delta G^\circ_{37\text{, total}} = \Delta G^\circ_{\text{init}} + \sum_{i=1}^{L-1} \Delta G^\circ_{\text{doublet}}(i, i+1) + N_{\text{term\_AU}} \cdot \Delta G^\circ_{\text{term\_AU}}$$

Burada:
-   $\Delta G^\circ_{\text{init}} = +3.4 \text{ kcal/mol}$ (Heliks çekirdeklenme enerjisi)
-   $\Delta G^\circ_{\text{term\_AU}} = +0.5 \text{ kcal/mol}$ (Uçlarda kararsızlaştırıcı AU çifti cezası)
-   $N_{\text{term\_AU}}$: Heliks uçlarındaki terminal A veya U bazlarının sayısı ($0, 1$ veya $2$).

---

## 4. 3D Moleküler Docking ve Metropolis Monte Carlo Simülasyonu

Aday moleküllerin SMN2 ekzon 7 5' splice site bölgesindeki deneysel 3D NMR yapısına (PDB 8R62) bağlanması ve vücut sıcaklığındaki kararlılığı simüle edilir.

### 4.1. AutoDock Vina Puanlama Fonksiyonu
AutoDock Vina, ampirik ve bilgi tabanlı (knowledge-based) bir serbest bağlanma enerjisi puanlama fonksiyonu kullanır:

$$c(t) = h(d(t)) \quad \text{where} \quad d(t) = r_i - r_j$$

Serbest enerji skor fonksiyonu ($\Delta G_{\text{bind}}$):

$$\Delta G_{\text{bind}} = w_1 \cdot \text{gauss}_1(d) + w_2 \cdot \text{gauss}_2(d) + w_3 \cdot \text{repulsion}(d) + w_4 \cdot \text{hydrophobic}(d) + w_5 \cdot \text{hbond}(d) + w_{\text{rot}} \cdot N_{\text{rot}}$$

Terimlerin tanımları:
-   **$\text{gauss}_1(d) = e^{-(d/0.5)^2}$** (Kısa mesafeli çekim)
-   **$\text{gauss}_2(d) = e^{-((d-3.0)/2.0)^2}$** (Uzun mesafeli çekim)
-   **$\text{repulsion}(d) = d^2$** ($d < 0$ ise; sterik çakışma)
-   **$\text{hydrophobic}(d)$:** Hidrofobik atomlar arası temas skoru.
-   **$\text{hbond}(d)$:** Hidrojen bağı donör-akseptör geometrik uyumu.
-   **$N_{\text{rot}}$:** Moleküldeki serbest dönebilen bağ sayısı (konformasyonel entropi kaybı cezası).

### 4.2. Metropolis Monte Carlo Dinamiği ve Isı Fluktuasyonları
Vücut sıcaklığında ($T = 310\text{ K}$) moleküllerin cep içindeki hareketliliği ve kararlılığı Metropolis Monte Carlo (MMC) algoritmasıyla doğrulanır.
1.  **Potansiyel Enerji ($V$):** Lennard-Jones (12-6) potansiyeli ve hidrojen bağ çekimleri kullanılarak hesaplanır:
    $$V_{\text{LJ}} = \sum_{i \in \text{ligand}} \sum_{j \in \text{pocket}} 4\epsilon \left[ \left(\frac{\sigma}{r_{ij}}\right)^{12} - \left(\frac{\sigma}{r_{ij}}\right)^6 \right]$$
    Burada $\sigma = 3.2\text{ Å}$ ve $\epsilon = 0.15\text{ kcal/mol}$ olarak parametrize edilmiştir. İki atom birbirine $1.2\text{ Å}$'den fazla yaklaşırsa $+150\text{ kcal/mol}$ sterik çakışma (clash) cezası uygulanır.
2.  **Perturbasyon:** Her adımda ligand koordinatlarına rastgele öteleme ($\Delta x \in [-0.06, 0.06]\text{ Å}$) ve üç eksende dönme ($\Delta \theta \in [-0.025, 0.025]\text{ rad}$) uygulanır.
3.  **Metropolis Kabul Kriteri:** Enerji değişimi $\Delta E = E_{\text{yeni}} - E_{\text{eski}}$ hesaplanır.
    -   Eğer $\Delta E \le 0$ ise adım **kabul edilir** ($P_{\text{accept}} = 1.0$).
    -   Eğer $\Delta E > 0$ ise adım Boltzmann olasılığı ile kabul edilir:
        $$P_{\text{accept}} = e^{-\frac{\Delta E}{k_B T}}$$
        Burada $k_B = 0.001987\text{ kcal/(mol}\cdot\text{K)}$ ve $T = 310\text{ K} \implies k_B T \approx 0.616\text{ kcal/mol}$.
4.  **RMSD (Root-Mean-Square Deviation) Takibi:** Simülasyon süresince ligandın docked başlangıç konumuna göre yer değiştirmesi ölçülür:
    $$\text{RMSD} = \sqrt{\frac{1}{N} \sum_{i=1}^{N} (x_i - x_{i,0})^2 + (y_i - y_{i,0})^2 + (z_i - z_{i,0})^2}$$
    RMSD değeri 1.000 adım sonunda $2.2\text{ Å}$ altında kalan moleküller *STABLE*, $4.0\text{ Å}$ üzerine çıkanlar ise cepten kaçtıkları için *UNSTABLE* olarak sınıflandırılır.

---

## 5. Referans Bilimsel Makaleler ve Literatür Kaynakları

1.  **Lefebvre, S., et al. (1995).** "Identification and characterization of a spinal muscular atrophy-determining gene." *Cell*, 80(1), 155-165.
2.  **Freier, S. M., et al. (1986).** "Improved free-energy parameters for predictions of RNA duplex stability." *Proceedings of the National Academy of Sciences*, 83(24), 9373-9377.
3.  **Wager, T. T., et al. (2010).** "Moving beyond Rules: The development of a central nervous system multiparameter optimization (CNS MPO) approach to enable alignment of druggability parameters for rational drug design." *ACS Chemical Neuroscience*, 1(6), 435-449.
4.  **Bickerton, G. R., et al. (2012).** "Quantifying the chemical beauty of drugs." *Nature Chemistry*, 4(2), 90-98.
5.  **Tran, Q. D., et al. (2020).** "Structural basis for targeted splicing modulation by small molecules." *Nature*, 581(7806), 105-110. (PDB 8R62 modelleme referansı).
6.  **Finkel, R. S., et al. (2017).** "Nusinersen versus Sham Control in Infantile-Onset Spinal Muscular Atrophy." *New England Journal of Medicine*, 377(18), 1723-1732. (ENDEAVOR/NURTURE klinik veri kaynağı).
7.  **Baranello, G., et al. (2021).** "Risdiplam in Type 1 Spinal Muscular Atrophy." *New England Journal of Medicine*, 384(10), 915-923. (FIREFISH klinik veri kaynağı).
