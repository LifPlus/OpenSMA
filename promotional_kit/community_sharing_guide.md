# OpenSMA Açık Kaynak Yayılım ve Tanıtım Kılavuzu

Bu kılavuz, projenizi GitHub'a yükledikten sonra küresel bilimsel ve yazılımsal topluluklara ulaştırmak için adım adım yapmanız gerekenleri içeren bir yol haritasıdır.

---

## ADIM 1: GitHub Reposunu Hazırlama ve Canlıya Alma
1.  **GitHub Hesabınızda Yeni Bir Repo Oluşturun:**
    *   İsim: `OpenSMA` veya `OpenSMA-Drug-Discovery`
    *   Açıklama: `An open-source, end-to-end computational drug discovery pipeline and systems biology simulator for Spinal Muscular Atrophy.`
    *   Erişim: **Public (Kamuya Açık)**.
2.  **Dosyaları Yükleyin:**
    *   `OpenSMA_OpenSource_Release` klasörünün içindeki tüm dosyaları (kodlar, veri setleri, `docs/` ve `promotional_kit/` klasörleri dahil) repoya yükleyin.
3.  **Etiketleri (Topics) Seçin:**
    *   Repo ayarlarından şu etiketleri ekleyin: `desci`, `bioinformatics`, `drug-discovery`, `python`, `sma`, `rdkit`, `pharmacokinetics`, `molecular-dynamics`.

---

## ADIM 2: Sosyal Medyada Paylaşım (Görünürlük)
1.  **Twitter/X Zinciri Yayınlayın:**
    *   `twitter_thread_templates.md` dosyasındaki Şablon 1'i kopyalayın.
    *   Tweetlere `full_cure_simulation.png` ve `patient_simulation_plots.png` grafiklerini ekleyin.
    *   Özellikle `#DeSci`, `#Bioinformatics` ve `#RareDiseases` hashtag'lerini kullanın.
2.  **LinkedIn'de Profesyonel Gönderi Paylaşın:**
    *   `linkedin_and_reddit_posts.md` dosyasındaki LinkedIn şablonunu kullanın.
    *   Kişisel ağınızdaki biyoteknoloji, eczacılık ve genetik araştırmacısı tanıdıklarınızı etiketleyin.

---

## ADIM 3: Bilimsel ve DeSci Forumları
1.  **Reddit Paylaşımları:**
    *   **r/bioinformatics:** `linkedin_and_reddit_posts.md` dosyasındaki teknik Reddit şablonunu paylaşın. Bu toplulukta kodun temizliği ve algoritmalar çok soru alacaktır, gelen yorumlara teknik olarak cevap verin.
    *   **r/desci:** DeSci topluluğuna hitap eden bir özet paylaşın.
2.  **ResearchHub Yüklemesi:**
    *   [researchhub.com](https://www.researchhub.com) adresine gidin.
    *   "Upload" butonuna tıklayarak `OpenSMA_Master_Technical_Report.md` veya `OpenSMA_Academic_Monograph.md` dosyasını belge olarak yükleyin.
    *   Açıklamaya `desci_pitch_templates.md` içindeki ResearchHub özetini ekleyin.
3.  **LabDAO / Molecule Discord Sunucuları:**
    *   Toplulukların Discord davet linklerini bularak sunucularına katılın (LabDAO: `discord.gg/labdao`, Molecule: `molecule.to` üzerinden yönlendirmeler).
    *   `desci_pitch_templates.md` dosyasındaki Discord tanıtım metnini ilgili `#collaboration` veya `#project-ideas` kanallarına gönderin.

---

## ADIM 4: Geri Bildirimleri Takip Etme ve Güncelleme
*   **Issues (Sorunlar):** Kullanıcıların veya araştırmacıların GitHub üzerinde açacağı "Issue" bildirimlerini takip edin. Kodun çalışması veya parametrelerle ilgili sorular gelebilir.
*   **PRs (Pull Requests):** Diğer biyoinformatikçiler kodunuzu optimize etmek veya DMD (Duchenne) gibi başka bir hastalığa uyarlamak için kod eklemek isteyebilir (PR gönderebilirler). Bunları inceleyerek birleştirin.
*   **Wet-Lab İletişimleri:** Wet-lab deneysel protokolümüzü görüp kendi laboratuvarında test etmek isteyen akademisyenler çıkacaktır. Onlarla işbirliği kurarak in-vitro test sonuçlarını (RT-qPCR jel görüntüleri, WB band yoğunlukları) repoya deneysel doğrulamalar (`experimental_validation/` klasörü) olarak ekleyin. Bu, projenizi sıradan bir yazılımdan **uluslararası akredite bir bilimsel yayına** dönüştürecektir.
