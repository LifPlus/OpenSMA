# OpenSMA Twitter/X Tanıtım Zinciri (Thread) Şablonları

Bu dosya, projenizi Twitter/X platformunda paylaşarak en yüksek etkileşimi ve bilimsel takibi (özellikle DeSci topluluğundan) elde etmeniz için tasarlanmış tweet zinciri şablonlarını içerir.

---

## ŞABLON 1: Bilimsel ve DeSci Odaklı Zincir (Önerilen)

### Tweet 1: Giriş (Hook)
> 1/ Spinal Müsküler Atrofi (SMA) tedavileri (Zolgensma, Spinraza) dünyanın en pahalı ilaçları arasında. Peki bu süreci demokratikleştirebilir miyiz?
> 
> Tamamen açık kaynaklı, 12 aşamalı hesaplamalı ilaç keşfi ve klinik PBPK simülasyon platformu **OpenSMA**'yı yayına aldık! 🧵👇
> 
> #DeSci #OpenScience #Bioinformatics #RareDiseases

### Tweet 2: ASO ve Termodinamik Tasarım
> 2/ Platform, pre-mRNA ISS-N1 cebini hedefleyen ASO'ları Freier 1986 Nearest-Neighbor modeliyle tarıyor.
> 
> FDA onaylı Nusinersen referansı $-24.6\text{ kcal/mol}$ bağlanma enerjisi verirken, de novo tasarladığımız aday **$-36.3\text{ kcal/mol}$** serbest enerji ile kararlılık sınırlarını zorluyor. 🧬

### Tweet 3: 3D Docking ve Metropolis Monte Carlo Cebin Testi
> 3/ Tasarlanan küçük molekülleri PDB 8R62 U1-RNA cebinde AutoDock Vina ile dock ettik (Afinite: $-7.7\text{ kcal/mol}$).
> 
> Ardından $310\text{ K}$ (vücut sıcaklığı) altında 1000 adımlık **Metropolis Monte Carlo** simülasyonuyla termal cep kararlılığını test ettik. Dinamik kaçış analizi hazır! 💻

### Tweet 4: Klinik PBPK Modelleme
> 4/ Basit PK/PD modellerinin ötesine geçerek; bağırsak, kan plazması, kan-beyin bariyeri (BBB) ve BOS (omurilik sıvısı) geçişlerini içeren **5-bölmeli PBPK ODE** sistemi kurguladık.
> 
> Modeli Finkel 2017 ve Baranello 2021 klinik verileriyle kalibre ettik (Uyum RMSE: $0.33$ ve $0.28$).

### Tweet 5: 10.000 Sanal Hasta Simülasyonu (Monte Carlo)
> 5/ Kalıcı Base Editing (CBE) gen tedavisi öncesi motor nöronları koruyacak "ASO Köprüsü" içeren **Full Cure Kombinasyon Protokolünü** test ettik.
> 
> 5 yıllık Monte Carlo klinik kohort simülasyon sonuçları:
> 📈 Kalıcı hayatta kalma: `%98.7`
> 🚶 Bağımsız yürüme olasılığı: `%50.2`

### Tweet 6: Çağrı ve Kodlar (Call to Action)
> 6/ Projenin tüm kodları, veri setleri, matematiksel denklemleri ve en önemlisi fibroblast testleri için hazırladığımız **Islak Laboratuvar (Wet-Lab) Protokolü** açık kaynaklıdır!
> 
> Katkı sunmak, verileri doğrulamak ve laboratuvarda denemek için: 
> 🔗 [GITHUB LINKINIZI BURAYA EKLEYİN]
> 
> #OpenSource #RareDisease #Biotech

---

## Görsel İpuçları (Tweetlere Eklenecek Resimler)
*   **Tweet 1 veya 5'e:** Platformun ürettiği `full_cure_simulation.png` grafiğini ekleyin.
*   **Tweet 4'e:** `patient_simulation_plots.png` grafiğini ekleyin.
