# OpenSMA: Açık Kaynaklı (DeSci) Spinal Müsküler Atrofi İlaç Keşfi ve Simülasyon Platformu

OpenSMA, Spinal Müsküler Atrofi (SMA) hastalığı için yeni nesil tedavi alternatifleri (küçük moleküller, ASO ve CRISPR) geliştirmek ve bunların klinik/farmakokinetik süreçlerini simüle etmek amacıyla kurulmuş **Merkeziyetsiz Bilim (DeSci)** odaklı, açık kaynaklı bir hesaplamalı tıp ve biyoloji projesidir.

---

## Proje Yapısı

Proje dosyaları aşağıdaki organizasyonla yer almaktadır:

```text
OpenSMA_OpenSource_Release/
├── LICENSE                     # MIT Açık Kaynak Lisansı
├── requirements.txt            # Python bağımlılık listesi
├── README.md                   # Bu bilgilendirme dosyası
├── data/                       # 3D Protein ve genetik girdi verileri
│   ├── smn2_receptor_real.pdb  # SMN2 splicing kompleksi ham PDB yapısı
│   ├── receptor_clean.pdb      # MMC ve Docking için temizlenmiş 3D yapı
│   └── smn2_target_region.fasta # ASO tasarımı için hedef mRNA bölgesi
├── docs/                       # Bilimsel dokümantasyon ve lab protokolleri
│   ├── OpenSMA_Lab_Protocol.md # Islak laboratuvar deney validasyon kılavuzu
│   ├── OpenSMA_Master_Technical_Report.md # Model matematiksel/fiziksel arka plan raporu
│   └── Production_and_Cost_Analysis.md # Üretim ve maliyet analiz raporu
├── src/                        # 12 Aşamalı boru hattı python kodları
│   ├── run_opensma_pipeline.py # Ana yürütücü script (End-to-end executor)
│   ├── molecule_optimizer.py   # De Novo küçük molekül üretici (Genetik Algoritma)
│   ├── aso_designer.py         # Termodinamik ASO tasarım modülü
│   ├── crispr_designer.py      # CRISPR/CBE gRNA tasarım modülü
│   ├── admet_screener.py       # ADMET filtreleme modülü
│   ├── aso_toxicity.py         # ASO hedef dışı toksisite tarayıcı
│   ├── docking_engine.py       # Open3DAlign hiza ve Vina docking motoru
│   ├── toxicity_model.py       # Tox21 ML toksisite tahmin modeli
│   ├── ode_calibration_and_screening.py # Literatür kalibrasyonlu ODE motoru
│   ├── final_scorer.py         # Çok kriterli sıralama ve skorer
│   ├── patient_sim.py          # Çok bölmeli PBPK hasta simülasyonu
│   ├── full_cure_sim.py        # 10.000 hasta Monte Carlo klinik simülasyonu
│   └── pocket_dynamics.py      # Metropolis Monte Carlo 3D bağlanma kararlılığı
└── results/                    # Simülasyon çıktılarının kaydedileceği klasör
```

---

## 🛠️ Kurulum ve Çalıştırma

### 1. Bağımlılıkların Yüklenmesi
Platform, kimyasal ve biyolojik hesaplamalar için **RDKit** ve veri analizi için **SciPy**, **Scikit-learn**, **Pandas** ve **Matplotlib** kullanmaktadır.

Sanal bir Python ortamı (virtualenv) oluşturup bağımlılıkları yükleyin:

```bash
# Sanal ortam oluşturun
python3 -m venv sma_env
source sma_env/bin/activate

# Bağımlılıkları yükleyin
pip install -r requirements.txt
```

### 2. Boru Hattını (Pipeline) Çalıştırma
Tüm 12 aşamalı keşif ve klinik simülasyon döngüsünü başlatmak için:

```bash
python src/run_opensma_pipeline.py
```

Boru hattı sırasıyla küçük molekülleri üretecek, ASO'ları termodinamik modellerle filtreleyecek, 3D yerleşimlerini hizalayacak, farmakokinetik ve toksisite filtrelerinden geçirecek ve ardından sanal hastalar üzerinde çok bölmeli PBPK ODE çözümleri ve Monte Carlo kohort analizleri yaparak sonuçları görselleştirecektir.

---

## 🔬 Islak Laboratuvar Validasyonu

Hesaplamalı modellerle üretilen aday ilaç bileşiklerinin laboratuvarda test edilmesi için geliştirilen **RT-qPCR, Western Blot, SMI-32 Motor Nöron Hücre Canlılığı ve delta7 Fare Modeli** protokollerine [OpenSMA_Lab_Protocol.md](docs/OpenSMA_Lab_Protocol.md) dosyasından ulaşabilirsiniz.

---

## 📜 Lisans ve Katkı
Bu proje **MIT Lisansı** ile lisanslanmıştır. Dünyanın her yerindeki bilim insanları, eczacılar ve araştırmacılar bu kod tabanını kendi yerel laboratuvarlarında pre-klinik araştırmalar yapmak amacıyla özgürce indirebilir, modifiye edebilir ve kullanabilirler.
