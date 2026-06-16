# OpenSMA: Açık Kaynaklı (DeSci) Islak Laboratuvar Validasyon Protokolü

Bu doküman, OpenSMA platformu tarafından tasarlanan küçük molekül ve ASO (Antisens Oligonükleotid) adaylarının laboratuvar ortamında (in vitro ve in vivo) nasıl test edileceğini, formüle edileceğini ve doğrulanacağını açıklayan adım adım deneysel validasyon protokolüdür.

---

## 1. Deneysel Tasarım ve Başlangıç Materyalleri

### 1.1. Hücre Hatları (Cell Lines)
*   **SMA Hasta Fibroblast Hücre Hattı:** GM03813 (Coriell Institute). Genotip: *SMN1* homozigot delesyon, *SMN2* 2 kopya (Tip 1 SMA klinik fenotipi).
*   **Kontrol Hücre Hattı:** GM03815 (Sağlıklı taşıyıcı anne/baba fibroblastları).
*   **Genel Ekspresyon Kontrolü:** HEK293T (ATCC CRL-11268).

### 1.2. Reaktifler ve Kimyasallar
*   **Kültür Medyumu:** DMEM (Dulbecco's Modified Eagle Medium - Gibco, Cat #11965092), %10 Fetal Bovine Serum (FBS - Gibco), %1 Penicillin-Streptomycin.
*   **Transfeksiyon Reaktifi:** Lipofectamine 3000 (Invitrogen) veya MC3 bazlı Lipit Nanopartiküller (LNP).
*   **ASO Kimyası:** 2'-O-metoksietil (2'-MOE) modifiyeli fosforotiyoat (PS) omurgalı oligonükleotidler (özel sentez firmalarından PAGE veya HPLC saflığında temin edilir).
*   **Birincil Antikorlar:** Mouse anti-SMN (Clone 2B1, Sigma-Aldrich, Cat #S2944), Rabbit anti-Beta-Actin (Cell Signaling Technology, Cat #4970).
*   **SMI-32 Antikoru:** Mouse anti-Neurofilament H Non-Phosphorylated (SMI-32, BioLegend, Cat #801701 - Motor nöron spesifik işaretleyici).

---

## 2. Protokol A: ASO ve LNP Transfeksiyonu (In Vitro Splicing Modülasyonu)

### 2.1. LNP Formulasyonu (Mikroakışkan Çip Yöntemi)
1.  **Lipit Hazırlığı:** İyonize edilebilir lipit DLin-MC3-DMA, DSPC, Kolesterol ve DMG-PEG2000 sırasıyla 50:10:38.5:1.5 molar oranında tartılır. Susuz %100 Etanol içinde çözülerek toplam lipit konsantrasyonu 10 mM yapılır.
2.  **ASO Hazırlığı:** Sentezlenen ASO adayı, 50 mM Sodyum Sitrat tamponunda (pH 4.0) çözülerek 0.2 mg/mL nükleik asit çözeltisi elde edilir.
3.  **Karıştırma (Formülasyon):** Mikroakışkan karıştırıcıda sulu faz ile organik etanol fazı 3:1 akış hızı oranında (Total akış hızı = 12 mL/dk) hızla çarpıştırılır.
4.  **Diyaliz:** Karışım, PBS (pH 7.4) tamponuna karşı Slide-A-Lyzer diyaliz kaseti (10k MWCO) kullanılarak 4°C'de 16 saat boyunca diyaliz edilir (etanolün uzaklaştırılması için).
5.  **Boyut Analizi:** DLS (Dinamik Işık Saçılımı) cihazında LNP boyutu (hedef: 70-90 nm) ve polidispersite indeksi (PDI < 0.15) doğrulanır.

### 2.2. Fibroblast Hücrelerine Transfeksiyon
1.  GM03813 hasta fibroblastları, 6 kuyucuklu plakalara kuyucuk başına $2 \times 10^5$ hücre olacak şekilde ekilir ve %5 $CO_2$, 37°C ortamında 24 saat inkübe edilir.
2.  Hücreler %70-80 yoğunluğa ulaştığında, besi ortamı serumsuz Opti-MEM (Gibco) ile değiştirilir.
3.  ASO-LNP formülasyonu (veya Lipofectamine 3000 ile komplekslenmiş çıplak ASO) kuyucuklara farklı konsantrasyonlarda (10 nM, 50 nM, 100 nM, 200 nM) eklenir. Kontrol kuyucuklarına sadece boş LNP (Scrambled ASO) verilir.
4.  4-6 saat sonra besi ortamı standart %10 FBS'li DMEM ile değiştirilir. Hücreler analiz için 48 saat boyunca inkübe edilir.

---

## 3. Protokol B: RT-qPCR ile Ekzon 7 Dahil Etme Oranının Analizi

### 3.1. Toplam RNA İzolasyonu ve cDNA Sentezi
1.  Transfekte edilen hücrelerden TRIzol (Invitrogen) reaktifi kullanılarak standart kloroform faz ayrımı yöntemiyle toplam RNA izole edilir.
2.  RNA konsantrasyonu ve saflığı Nanodrop spektrofotometresinde ($A_{260}/A_{280} \ge 1.8$) doğrulanır.
3.  1 µg toplam RNA, High-Capacity cDNA Reverse Transcription Kit (Applied Biosystems) kullanılarak rastgele hekzametrik primerlerle cDNA'ya dönüştürülür.

### 3.2. Yarı Kantitatif PCR ve qPCR Koşulları
*   **Ekzon 7 Dahil Edilme Durumunu Belirleyen RT-PCR Primerleri (Human SMN2 spesifik):**
    *   **Forward (Exon 6):** 5'-GCTATCATAATTTTGTTCATTTT-3'
    *   **Reverse (Exon 8):** 5'-CCAGATTCTCTTGATGATGC-3'
*   **PCR Reaksiyon Karışımı:** 2 µL cDNA, 1.25 U Taq DNA Polimeraz, 0.4 µM her bir primer, 0.2 mM dNTP, PCR Tamponu (Toplam hacim: 25 µL).
*   **PCR Termal Döngü Protokolü:**
    1.  Başlangıç Denatürasyonu: 95°C, 3 dk
    2.  Döngü (30 tekrar):
        *   Denatürasyon: 95°C, 30 sn
        *   Bağlanma (Annealing): 55°C, 30 sn
        *   Uzama (Extension): 72°C, 45 sn
    3.  Nihai Uzama: 72°C, 5 dk
4.  **Jel Elektroforezi:** PCR ürünleri %2'lik agaroz jelde ethidium bromide eşliğinde yürütülür.
    *   **Beklenen Bant Boyutları:** Ekzon 7 içeren tam boy transcript (Full-length SMN2) = **263 bp**, Ekzon 7 içermeyen transcript (Delta7-SMN2) = **209 bp**.
5.  Bant yoğunlukları ImageJ programı ile analiz edilir. Ekzon 7 dahil edilme yüzdesi aşağıdaki formülle hesaplanır:
    $$\text{Ekzon 7 Inclusion \%} = \frac{\text{Yoğunluk (263 bp Bandı)}}{\text{Yoğunluk (263 bp Bandı)} + \text{Yoğunluk (209 bp Bandı)}} \times 100$$

---

## 4. Protokol C: Western Blot ile SMN Protein Seviyesinin Ölçümü

### 4.1. Protein İzolasyonu
1.  Hücreler buz üzerinde RIPA Lysis Buffer (proteaz ve fosfataz inhibitör kokteyli içeren) ile lize edilir.
2.  Lizatlar 4°C'de 14.000 rpm'de 15 dk santrifüj edilerek süpernatant alınır.
3.  Protein konsantrasyonu BCA Protein Assay Kit (Thermo Fisher) ile tayin edilir.

### 4.2. SDS-PAGE ve Membrana Aktarım
1.  30 µg protein örneği Laemmli yükleme tamponu ile 95°C'de 5 dk denatüre edilir.
2.  Örnekler %10'luk SDS-PAGE jeline yüklenir ve 100V'ta yürütülür.
3.  Jeldeki proteinler PVDF membrana (Bio-Rad) ıslak aktarım (Wet Transfer) yöntemiyle 300 mA'de 1.5 saat boyunca transfer edilir.

### 4.3. Antikor İnkübasyonu ve Görüntüleme
1.  Membran, %5 yağsız süt tozu içeren TBST (Tris-buffered saline, %0.1 Tween-20) içinde oda sıcaklığında 1 saat bloke edilir.
2.  Birincil antikorlar (Mouse anti-SMN, 1:1000 ve Rabbit anti-Beta-Actin, 1:5000) ile 4°C'de 16 saat inkübe edilir.
3.  TBST ile 3 kez 10 dk yıkama yapılır.
4.  HRP-konjuge ikincil antikorlar (Anti-mouse IgG-HRP ve Anti-rabbit IgG-HRP, 1:10.000) ile oda sıcaklığında 1 saat inkübe edilir.
5.  ECL Western Blotting Substrate (Pierce) kullanılarak kemilüminesans görüntüleme cihazında (ChemiDoc) bantlar görüntülenir. Beta-Aktin referansına göre SMN protein artış katı hesaplanır.

---

## 5. Protokol D: Motor Nöron Hayatta Kalma ve SMI-32 Staining Analizi

Küçük moleküllerin veya ASO'ların dejenerasyonu durdurma yetisi ko-kültür modellerinde test edilir:
1.  **Hücre İzolasyonu:** Gebeliğin 14. günündeki (E14) fare embriyolarından spinal kord motor nöronları izole edilerek Neurobasal Medium (%2 B-27, 0.5 mM L-Glutamine, 10 ng/mL BDNF, GDNF ve CNTF içeren) ortamında kültürlenir.
2.  **SMA Koşullarının Taklidi:** Motor nöronlara siRNA veya düşük doz toksin verilerek dejenerasyon tetiklenir veya hasta fibroblastlarından elde edilen şartlandırılmış medyum (conditioned medium) eklenir.
3.  **İlaç Uygulaması:** Aday küçük moleküller (optimum doz: 10 nM - 10 µM) veya ASO'lar hücrelere eklenir.
4.  **İmmünofloresan Boyama (SMI-32):**
    *   Hücreler 48 saat sonra %4 Paraformaldehit (PFA) ile fikse edilir.
    *   Triton X-100 (%0.1) ile geçirgenleştirme yapılır.
    *   **SMI-32 Birincil Antikoru** (1:1000) ile 4°C'de gece boyu inkübe edilir.
    *   Alexa Fluor 488 konjuge ikincil antikor (1:500) ile oda sıcaklığında 1 saat inkübe edilir.
    *   DAPI ile çekirdek boyaması yapılır.
5.  **Hayatta Kalma Sayımı:** SMI-32 pozitif olan, uzun aksonlara sahip sağlıklı motor nöron sayıları floresan mikroskobu altında sayılır. Tedavi edilmeyen dejeneratif kontrol grubu ile kurtarılan grup karşılaştırılır.

---

## 6. Protokol E: In Vivo Validasyon (SMA Delta7 Fare Modeli)

### 6.1. Hayvan Modeli
*   **Model:** $Smn^{-/-}; SMN2^{+/+}; SMN\Delta7^{+/+}$ (Jackson Laboratory, Strain #005025). Bu fareler insan Tip 1 SMA kliniğini taklit eder ve tedavi edilmezlerse ortalama 13-15 gün yaşarlar.

### 6.2. Dosing ve Uygulama Yolları
*   **Küçük Moleküller (Oral):** Doğumdan sonraki 1. günden (P1) itibaren farelere günlük oral gavage yoluyla 1-10 mg/kg dozunda ilaç çözeltisi (mikro-şırınga ile ağızdan) verilir.
*   **ASO Tedavisi (İntratekal Enjeksiyon):** 
    *   Doğumdan sonraki 2. günde (P2) farelere kriyojenik anestezi (buz üzerinde hafif dondurma) altında intratekal (spinal kanala doğrudan) enjeksiyon yapılır.
    *   33G Hamilton mikro-şırıngası ile L3-L4 seviyesinden 2 µL hacimde LNP-ASO formülasyonu (doz: 10-20 µg) enjekte edilir.

### 6.3. Klinik İzlem Parametreleri
1.  **Sağkalım Analizi (Survival Curve):** Farelerin doğal ölüm günleri kaydedilerek Kaplan-Meier eğrileri oluşturulur.
2.  **Vücut Ağırlığı Gelişimi:** P1'den itibaren her gün hassas terazi ile tartılarak gelişim grafiği çıkarılır.
3.  **Motor Fonksiyon Testleri (Righting Reflex):** Fare sırt üstü bırakılır. Tekrar dört ayağı üzerine dönme süresi (Righting time) saniye cinsinden ölçülür. Sağlıklı farelerde bu süre <2 saniyedir. SMA farelerinde ise hastalık ilerledikçe bu süre uzar veya dönme yeteneği kaybolur.

---

## 7. Protokol F: HPLC ve LC-MS ile Kimyasal Saflık Doğrulaması

Sentezlenen veya izole edilen küçük moleküllerin saflığı in vitro testlerden önce kanıtlanmalıdır:
1.  **HPLC Koşulları:**
    *   Kolon: C18 Ters Faz (Reverse-Phase) Kolon (4.6 mm x 150 mm, 5 µm).
    *   Mobil Faz A: %0.1 Trifloroasetik asit (TFA) içeren ultra saf su.
    *   Mobil Faz B: %0.1 TFA içeren Asetonitril.
    *   Gradient Akış: %5 B'den %95 B'ye 20 dakikada geçiş, Akış Hızı: 1.0 mL/dk.
    *   Dedektör: UV (254 nm ve 280 nm).
    *   **Kabul Kriteri:** Ana pik alanının toplam pik alanına oranı $\ge \%98.5$ olmalıdır (kromatogram saflık doğrulaması).
2.  **LC-MS Kütle Analizi:** Moleküler iyonizasyon (ESI-MS) ile elde edilen kütle/yük ($m/z$) oranı, tasarlanan formülün teorik kütlesi ile tam uyum göstermelidir.
