import random
import numpy as np
import matplotlib.pyplot as plt
import time

E_BASLANGIC = random.uniform(2000000.0, 2000100.0)# eV (2 MeV - Fisyon Nötronu)
E_TERMAL = 0.025         # eV (Termal Enerji Sınırı)
SIMULASYON_NOTRON_SAYISI = 100000  # Her malzeme için simüle edilecek nötron sayısıe

class Atom:
    """
    Bir atom türünü tanımlar.
    A: Kütle Numarası
    sigma_s: Elastik Saçılma Tesir Kesiti (barn) - Ortalama/Efektif
    sigma_a: Absorpsiyon Tesir Kesiti (barn) - Ortalama/Efektif
    """
    def __init__(self, sembol, A, sigma_s, sigma_a):
        self.sembol = sembol
        self.A = A
        self.sigma_s = sigma_s
        self.sigma_a = sigma_a
        # Alfa parametresi (Çarpışma mekaniği için kilit parametre)
        # alpha = ((A-1)/(A+1))^2
        self.alpha = ((self.A - 1) / (self.A + 1)) ** 2

class Malzeme:
    """
    Moderatör malzemesini tanımlar (Örn: H2O, Grafit).
    bilesenler: {'AtomNesnesi': adet} sözlüğü (Örn: {H: 2, O: 1})
    """
    def __init__(self, isim, bilesenler, yogunluk):
        self.isim = isim
        self.bilesenler = bilesenler # Sözlük formatında: {atom_nesnesi: adet}
        self.yogunluk = yogunluk # g/cm3 (İleri hesaplamalar için eklendi)
        
        # Makroskobik saçılma olasılıklarını hesapla
        self.toplam_sigma_s = 0
        self.toplam_sigma_a = 0
        self.atom_olasiliklari = []
        
        for atom, adet in self.bilesenler.items():
            self.toplam_sigma_s += atom.sigma_s * adet
            self.toplam_sigma_a += atom.sigma_a * adet
            
        # Hangi atomla çarpışılacağını belirlemek için ağırlıklar
        simdiki_agirlik = 0
        for atom, adet in self.bilesenler.items():
            pay = (atom.sigma_s * adet) / self.toplam_sigma_s
            self.atom_olasiliklari.append((simdiki_agirlik + pay, atom))
            simdiki_agirlik += pay

    def hedef_atom_sec(self):
        """Çarpışma yapılacak atomu olasılığa göre seçer."""
        r = random.random()
        for esik, atom in self.atom_olasiliklari:
            if r <= esik:
                return atom
        return self.atom_olasiliklari[-1][1] # Güvenlik için sonuncuyu döndür

def carpisma_simulasyonu(notron_enerjisi, hedef_atom):
    """
    Elastik saçılma kinematiğini uygular.
    Geriye çarpışma sonrası yeni enerjiyi döndürür.
    """
    alpha = hedef_atom.alpha
    
    # Kütle Merkezi (CM) sisteminde saçılma açısı (varsayım: izotropik)
    # cos_theta -1 ile +1 arasında rastgele seçilir.
    cos_theta = 2 * random.random() - 1 
    
    # Enerji düşüş formülü:
    # E' = 0.5 * E * [ (1+alpha) + (1-alpha)*cos_theta ]
    yeni_enerji = 0.5 * notron_enerjisi * ((1 + alpha) + (1 - alpha) * cos_theta)
    
    return yeni_enerji

def moderasyon_testi(malzeme):
    """
    Belirli bir malzeme için N adet nötronu simüle eder.
    """
    carpisma_sayilari = []
    absorbe_sayisi = 0
    termalize_olanlar = 0
    
    print(f"--- {malzeme.isim} Simülasyonu Başlıyor ---")
    
    for i in range(SIMULASYON_NOTRON_SAYISI):
        E = E_BASLANGIC
        carpisma = 0
        aktif = True
        
        while E > E_TERMAL and aktif:
            carpisma += 1
            
            # 1. Adım: Absorpsiyon Kontrolü
            # Sigma_a / (Sigma_a + Sigma_s) olasılığı ile nötron yutulur
            P_abs = malzeme.toplam_sigma_a / (malzeme.toplam_sigma_a + malzeme.toplam_sigma_s)
            
            if random.random() < P_abs:
                aktif = False
                absorbe_sayisi += 1
            else:
                # 2. Adım: Saçılma (Hangi atomla çarpışacağız?)
                hedef = malzeme.hedef_atom_sec()
                
                # 3. Adım: Yeni Enerji Hesabı
                E = carpisma_simulasyonu(E, hedef)
        
        if aktif: # Nötron ölmedi ve termal enerjiye indi
            carpisma_sayilari.append(carpisma)
            termalize_olanlar += 1
            
    return carpisma_sayilari, termalize_olanlar, absorbe_sayisi

# --- 2. STANDART MALZEME TANIMLARI ---
H = Atom("Hidrojen", A=1, sigma_s=20.0, sigma_a=0.332) 
D = Atom("Döteryum", A=2, sigma_s=3.4, sigma_a=0.0005)
C = Atom("Karbon", A=12, sigma_s=4.8, sigma_a=0.0035)
O = Atom("Oksijen", A=16, sigma_s=3.8, sigma_a=0.0002)

# Malzemeler
su = Malzeme("Hafif Su (H2O)", {H: 2, O: 1}, yogunluk=1.0)
agir_su = Malzeme("Ağır Su (D2O)", {D: 2, O: 1}, yogunluk=1.1)
grafit = Malzeme("Grafit (C)", {C: 1}, yogunluk=1.7)
polietilen = Malzeme("Polietilen (CH2)", {C: 1, H: 2}, yogunluk=0.95)

# --- 3. KULLANICI TANIMLI YENİ MALZEME GİRİŞİ ---
def kullanici_malzeme_olustur():
    print("\n--- YENİ MODERATÖR TASARIMI (AR-GE MODU) ---")
    isim = input("Malzeme Adı (örn: YeniModeratör): ")
    
    # Basitleştirilmiş: Tek tip atomdan oluştuğunu varsayalım veya ortalama A girelim
    try:
        A_user = float(input("Ortalama Kütle Numarası (A) [örn: Li-7 için 7]: "))
        Sig_s_user = float(input("Ortalama Saçılma Tesir Kesiti (barn) [örn: 1.4]: "))
        Sig_a_user = float(input("Ortalama Absorpsiyon Tesir Kesiti (barn) [örn: 0.04]: "))
        Rho_user = float(input("Yoğunluk (g/cm3) [örn: 1.8]: "))
        
        # Yeni atom ve malzeme oluştur
        UserAtom = Atom("UserAtom", A=A_user, sigma_s=Sig_s_user, sigma_a=Sig_a_user)
        YeniMalzeme = Malzeme(isim, {UserAtom: 1}, yogunluk=Rho_user)
        return YeniMalzeme
    except ValueError:
        print("Hata: Lütfen sayısal değerler giriniz.")
        return None

# --- 4. ANA PROGRAM AKIŞI ---

malzeme_listesi = [su, agir_su, grafit, polietilen]

# Kullanıcıya sorma
print("Standart malzemelere ek olarak kendi malzemenizi test etmek ister misiniz? (E/H)")
cevap = input().upper()
if cevap == 'E':
    ozel_malzeme = kullanici_malzeme_olustur()
    if ozel_malzeme:
        malzeme_listesi.append(ozel_malzeme)

sonuclar = {}

print(f"\nSimülasyon Başlıyor... ({SIMULASYON_NOTRON_SAYISI} nötron/malzeme)\n")
for mat in malzeme_listesi:
    start_time = time.time()
    carpisma_listesi, termal, absorbe = moderasyon_testi(mat)
    end_time = time.time()
    
    ortalama_carpisma = np.mean(carpisma_listesi) if carpisma_listesi else 0
    verim = (termal / SIMULASYON_NOTRON_SAYISI) * 100
    
    sonuclar[mat.isim] = {
        "Ortalama Çarpışma": ortalama_carpisma,
        "Termalizasyon Verimi (%)": verim,
        "Absorbe Olan": absorbe,
        "Data": carpisma_listesi
    }
    
    print(f"--> {mat.isim} Tamamlandı. Ort. Çarpışma: {ortalama_carpisma:.1f}, Süre: {end_time-start_time:.2f}s")

# --- 5. SONUÇLARIN GÖRSELLEŞTİRİLMESİ ---

isimler = list(sonuclar.keys())
ort_carpisma_degerleri = [sonuclar[k]["Ortalama Çarpışma"] for k in isimler]
verim_degerleri = [sonuclar[k]["Termalizasyon Verimi (%)"] for k in isimler]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Grafik 1: Ortalama Çarpışma Sayısı (Düşük olması iyi yavaşlatıcı demek)
bars = ax1.bar(isimler, ort_carpisma_degerleri, color=['blue', 'cyan', 'gray', 'green', 'purple'][:len(isimler)])
ax1.set_title(f'2 MeV -> 0.025 eV İçin Ortalama Çarpışma Sayısı\n(Daha az = Daha Hızlı Yavaşlatma)')
ax1.set_ylabel('Çarpışma Sayısı')
ax1.grid(axis='y', linestyle='--', alpha=0.7)

# Değerleri barların üzerine yaz
for bar in bars:
    yval = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2, yval + 1, f'{yval:.1f}', ha='center', va='bottom')

# Grafik 2: Moderasyon Verimliliği (Absorbe olmadan termalize olma oranı)
bars2 = ax2.bar(isimler, verim_degerleri, color=['blue', 'cyan', 'gray', 'green', 'purple'][:len(isimler)])
ax2.set_title('Termalizasyon Başarısı (Absorpsiyondan Kaçış %)')
ax2.set_ylabel('Verim (%)')
ax2.set_ylim(0, 105)
ax2.grid(axis='y', linestyle='--', alpha=0.7)

for bar in bars2:
    yval = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2, yval + 1, f'%{yval:.1f}', ha='center', va='bottom')

plt.suptitle('ModSim: Moderatör Malzemeleri Karşılaştırmalı Analizi', fontsize=16)
plt.tight_layout()
plt.show()

# İstatistiksel Döküm
print("\n--- DETAYLI ANALİZ RAPORU ---")
print(f"{'Malzeme':<20} | {'Ort. Çarpışma':<15} | {'Termalizasyon %':<15} | {'Yorum'}")
print("-" * 80)
for isim, veri in sonuclar.items():
    c = veri["Ortalama Çarpışma"]
    v = veri["Termalizasyon Verimi (%)"]
    yorum = "Yüksek Verim" if v > 95 else "Absorpsiyon Riski"
    if c > 100: yorum += ", Hacim Gerektirir"
    else: yorum += ", Kompakt"
    print(f"{isim:<20} | {c:<15.1f} | %{v:<14.1f} | {yorum}")