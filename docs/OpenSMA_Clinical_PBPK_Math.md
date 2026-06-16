# OpenSMA: Klinik PBPK ve Sistem Biyolojisi Matematiksel Modelleme Rehberi

Bu rehber, OpenSMA platformunda kullanılan çok bölmeli PBPK (Fizyoloji Tabanlı Farmakokinetik) diferansiyel denklemlerinin, motor nöron korunum dinamiklerinin ve Monte Carlo klinik kohort simülasyonunun tüm matematiksel türetimlerini içerir.

---

## 1. 4-Bölmeli PBPK Diferansiyel Denklem Sisteminin Türetilmesi

İlacın oral veya intratekal yolla uygulanmasından itibaren santral sinir sistemine (CNS) dağılımı ve eliminasyonu, kütle korunumu (mass transport) yasalarına dayanan doğrusal olmayan diferansiyel denklemlerle (ODE) modellenir.

```
       [ Oral Dosing ]
              │
              ▼
        ┌───────────┐
        │  1. Gut   │ ──( ka )──► Plazma
        └───────────┘
                            ┌───────────┐
                            │ 2. Tissue │
                            └───────────┘
                                  ▲
                               ( k_dist )
                                  ▼
                        ┌───────────────────┐
                        │     3. Plasma     │ ──( k_elim )──► Eliminasyon
                        └───────────────────┘
                             ▲         │
                     ( k_out_bbb )  ( k_in_bbb )  [ Kan-Beyin Bariyeri ]
                             │         ▼
                        ┌───────────────────┐
                        │     4. Brain      │
                        └───────────────────┘
                             ▲         │
                      ( k_clear_csf ) ( k_diff_csf )
                             │         ▼
                        ┌───────────────────┐
                        │      5. CSF       │  ◄── [ Intrathecal Dosing ]
                        └───────────────────┘
```

### 1.1. Bölme Denklemleri
1.  **Mide/Bağırsak Bölmesi (Gut - $C_{\text{gut}}$):**
    Ağızdan alınan ilacın bağırsaktan emilerek plazmaya geçişi birinci derece kinetikle temsil edilir:
    $$\frac{dC_{\text{gut}}}{dt} = -k_a \cdot C_{\text{gut}}$$
    Burada $k_a$ emilim (absorption) hız sabitidir ($t^{-1}$).

2.  **Kan Plazması Bölmesi (Plasma - $C_{\text{plasma}}$):**
    Plazmadaki ilaç konsantrasyonu bağırsaktan gelen akış, periferal dokularla ($C_{\text{tissue}}$) pasif dağılım, böbreklerden atılım ve kan-beyin bariyeri (BBB) üzerinden beyin ile olan alışveriş dengesinden oluşur:
    $$\frac{dC_{\text{plasma}}}{dt} = \frac{k_a \cdot C_{\text{gut}} \cdot V_{\text{gut}}}{V_{\text{plasma}}} - k_{\text{elim}} \cdot C_{\text{plasma}} - k_{\text{dist}} \cdot (C_{\text{plasma}} - C_{\text{tissue}}) - k_{\text{in\_bbb}} \cdot C_{\text{plasma}} + k_{\text{out\_bbb}} \cdot \frac{V_{\text{brain}}}{V_{\text{plasma}}} \cdot C_{\text{brain}}$$
    Burada:
    -   $V_x$: İlgili bölmenin hacmi (L).
    -   $k_{\text{elim}}$: Eliminasyon hız sabiti.
    -   $k_{\text{dist}}$: Dokular arası dağılım hızı.
    -   $k_{\text{in\_bbb}}$ / $k_{\text{out\_bbb}}$: BBB giriş ve çıkış taşıma hızları.

3.  **Beyin Parankimi Bölmesi (Brain - $C_{\text{brain}}$):**
    İlacın hedef motor nöronların bulunduğu beyin parankimine geçişi ve buradan BOS'a (Beyin Omurilik Sıvısı) difüzyonu:
    $$\frac{dC_{\text{brain}}}{dt} = k_{\text{in\_bbb}} \cdot \frac{V_{\text{plasma}}}{V_{\text{brain}}} \cdot C_{\text{plasma}} - k_{\text{out\_bbb}} \cdot C_{\text{brain}} - k_{\text{diff\_csf}} \cdot (C_{\text{brain}} - C_{\text{csf}})$$

4.  **Beyin Omurilik Sıvısı Bölmesi (CSF - $C_{\text{csf}}$):**
    BOS konsantrasyonu intratekal (spinal iğne ile) ilaç enjeksiyonunun yapıldığı bölgedir. Beyin parankimi ile difüzyon ve BOS dolaşımı ile temizlenme (clearance) denklemi:
    $$\frac{dC_{\text{csf}}}{dt} = k_{\text{diff\_csf}} \cdot \frac{V_{\text{brain}}}{V_{\text{csf}}} \cdot (C_{\text{brain}} - C_{\text{csf}}) - k_{\text{clear\_csf}} \cdot C_{\text{csf}}$$

---

## 2. Farmakodinamik (PD) ve Splicing Kinetiği

Beyin parankimindeki veya omurilik sıvısındaki ilaç konsantrasyonu ($C_{\text{local}}$), pre-mRNA splicing mekanizmasını uyararak fonksiyonel SMN protein sentezini artırır.

### 2.1. Hill Denklemi ile Proteinin Uyarılması (Boost)
İlaç etkisine bağlı SMN protein artış hızı ($Boost_{\text{SMN}}$), klasik kooperatif bağlanmayı tanımlayan Hill denklemiyle hesaplanır:

$$\text{Boost}_{\text{SMN}}(t) = \frac{E_{\max} \cdot C_{\text{local}}(t)^n}{EC_{50}^n + C_{\text{local}}(t)^n}$$

Burada:
-   $E_{\max}$: Ulaşılabilecek maksimum artış katsayısı (% bazal değere göre).
-   $EC_{50}$: Maksimum etkinin yarısını sağlayan ilaç konsantrasyonu (afinite ölçüsü).
-   $n$: Hill katsayısı (kooperatiflik derecesi; platformda $n=1.2$ olarak kalibre edilmiştir).

### 2.2. Motor Nöron Kaybı Diferansiyel Denklemi
Tedavi edilmeyen Tip 1 SMA hastalarında motor nöron havuzu ($MN$) hızla erir. SMN protein restorasyonu ($smn(t)$), motor nöron dejenerasyon hızını yavaşlatır:

$$\text{protection}(t) = \min\left(1.0, \max\left(0.0, \frac{smn(t)}{100.0}\right)\right)$$

$$\frac{d[MN]}{dt} = -r_{\text{decay}} \cdot \left(1.0 - \text{protection}(t) \cdot \eta \right) \cdot [MN]$$

Burada:
-   $r_{\text{decay}}$: Tedavi edilmeyen hastadaki bazal motor nöron kayıp hızı ($0.05 \text{ /ay}$).
-   $\eta$: Tedavinin motor nöronları koruma verimliliği ($0.85$ ila $0.95$ arası).
-   $[MN]$: Kalan motor nöron yüzdesi (% başlangıç değerine göre).

---

## 3. Markov Modellemesi ve Monte Carlo Kohort Simülasyonu

Simüle edilen 10.000 sanal hasta üzerinde, her hastanın motor nöron canlılığına ($[MN](t)$) bağlı olarak klinik motor dönüm noktalarına ulaşma olasılıkları Markov geçiş olasılıkları ile hesaplanır.

### 3.1. Motor Milestone Olasılık Fonksiyonları
Motor fonksiyon kazanımları (bağımsız oturma, ayakta durma, yürüme), motor nöron havuzunun büyüklüğüne bağlı doğrusal olmayan olasılık fonksiyonlarıdır:

$$P(\text{Sit}(t)) = P(\text{Sit}_{\max}) \cdot \left[ \frac{[MN](t) - MN_{\text{sit\_threshold}}}{100.0 - MN_{\text{sit\_threshold}}} \right]^2 \quad (\text{Eğer } [MN](t) > MN_{\text{sit\_threshold}})$$

$$P(\text{Walk}(t)) = P(\text{Walk}_{\max}) \cdot \left[ \frac{[MN](t) - MN_{\text{walk\_threshold}}}{100.0 - MN_{\text{walk\_threshold}}} \right]^3 \quad (\text{Eğer } [MN](t) > MN_{\text{walk\_threshold}})$$

Burada:
-   $MN_{\text{sit\_threshold}} = 35.0 \implies$ Bağımsız oturabilmek için motor nöron havuzunun en az %35'inin korunmuş olması gerekir.
-   $MN_{\text{walk\_threshold}} = 60.0 \implies$ Yürüyebilmek için en az %60 koruma şarttır.
-   Üsler ($2$ ve $3$), motor fonksiyon koordinasyonu için gereken sinaptik olgunlaşmanın doğrusal olmayan karmaşıklığını modellemek için konulmuştur.

### 3.2. Kaplan-Meier Sağkalım (Survival) Fonksiyonu
Monte Carlo döngüsünde, her sanal hastanın hayatta kalma durumu, motor nöron kaybının kritik eşik olan %20'nin altına düşüp düşmediği ile sorgulanır. Her ay için kohorttaki sağkalım olasılığı ($S(t)$) hesaplanır:

$$S(t) = \prod_{i=1}^{t} \left( 1 - \frac{d_i}{n_i} \right)$$

Burada:
-   $d_i$: $i$ ayında ölen (motor nöronu %20'nin altına düşen) sanal hasta sayısı.
-   $n_i$: $i$ ayının başındaki canlı sanal hasta sayısı.

Bu matematiksel altyapı, OpenSMA'nın 5 yıllık klinik kohort sonuçlarındaki yürüme ve hayatta kalma eğrilerinin bilimsel temelini oluşturur.
