import random
import numpy as np
import matplotlib.pyplot as plt
import time

# --- 1. FİZİKSEL SABİTLER ---
AVOGADRO = 0.6022        # 10^24 cinsinden
E_BASLANGIC = 2000000.0  # eV
E_TERMAL = 0.025         # eV
SIMULASYON_NOTRON_SAYISI = 10000 # Her malzeme için deneme sayısı

class Atom:
    def __init__(self, sembol, A, sigma_s, sigma_a):
        self.sembol = sembol
        self.A = A
        self.sigma_s = sigma_s # barn
        self.sigma_a = sigma_a # barn
        self.alpha = ((self.A - 1) / (self.A + 1)) ** 2

class Malzeme:
    def __init__(self, isim, bilesenler, yogunluk):
        self.isim = isim
        self.bilesenler = bilesenler
        self.yogunluk = yogunluk # g/cm3
        
        # 1. Mol Ağırlığı
        self.mol_agirligi = sum([atom.A * adet for atom, adet in bilesenler.items()])
        
        # 2. Moleküler Sayı Yoğunluğu (N_mol)
        self.N_molekul = (self.yogunluk * AVOGADRO) / self.mol_agirligi
        
        # 3. Makroskobik Tesir Kesitleri Hesabı
        self.Sigma_s_toplam = 0.0
        self.Sigma_a_toplam = 0.0
        self.atom_secim_listesi = [] 
        
        print(f"\n--- {self.isim} Verileri İşleniyor ---")
        print(f"  > Yoğunluk: {self.yogunluk} g/cm3, Mol Ağırlığı: {self.mol_agirligi:.2f} g/mol")
        
        for atom, adet in self.bilesenler.items():
            N_i = adet * self.N_molekul # Atomik yoğunluk
            Sigma_s_i = N_i * atom.sigma_s
            Sigma_a_i = N_i * atom.sigma_a
            
            self.Sigma_s_toplam += Sigma_s_i
            self.Sigma_a_toplam += Sigma_a_i
            
            # Seçim için listeye ekle
            self.atom_secim_listesi.append([atom, Sigma_s_i])
            print(f"  > {atom.sembol}: N={N_i:.4f}, Sig_s={Sigma_s_i:.4f} cm^-1")

        self.Sigma_tot = self.Sigma_s_toplam + self.Sigma_a_toplam
        self.P_absorpsiyon = self.Sigma_a_toplam / self.Sigma_tot
        
        # 4. Seçim Olasılıkları (Kümülatif)
        kümülatif = 0
        self.hedef_atom_araliklari = []
        for item in self.atom_secim_listesi:
            pay = item[1] / self.Sigma_s_toplam
            kümülatif += pay
            self.hedef_atom_araliklari.append((kümülatif, item[0]))

    def hedef_atom_sec(self):
        r = random.random()
        for esik, atom in self.hedef_atom_araliklari:
            if r <= esik: return atom
        return self.hedef_atom_araliklari[-1][1]

# --- SİMÜLASYON FONKSİYONLARI ---
def carpisma_kinematigi(enerji, atom):
    cos_theta = 2 * random.random() - 1
    alpha = atom.alpha
    return 0.5 * enerji * ((1 + alpha) + (1 - alpha) * cos_theta)

def moderasyon_simulasyonu(malzeme):
    print(f"\n>>> {malzeme.isim} Simülasyonu Başladı...")
    carpisma_verileri = []
    absorbe_sayisi = 0
    termalize_sayisi = 0
    
    start = time.time()
    for i in range(SIMULASYON_NOTRON_SAYISI):
        E = E_BASLANGIC
        carpisma = 0
        aktif = True
        while E > E_TERMAL and aktif:
            if random.random() < malzeme.P_absorpsiyon:
                aktif = False
                absorbe_sayisi += 1
            else:
                hedef = malzeme.hedef_atom_sec()
                carpisma += 1
                E = carpisma_kinematigi(E, hedef)
        
        if aktif:
            carpisma_verileri.append(carpisma)
            termalize_sayisi += 1
    
    sure = time.time() - start
    return carpisma_verileri, termalize_sayisi, absorbe_sayisi, sure

def analiz_modu():
    while True:
        print("\n" + "#"*40)
        print("      ANALİZ MODU: MALZEME PERFORMANS İNCELEMESİ      ")
        print("#"*40)
        print("Hangi parametrenin etkisini incelemek istiyorsunuz?")
        print("[1] YOĞUNLUK Etkisi -> Ortalama Serbest Yol")
        print("[2] ABSORPSİYON KESİTİ Etkisi -> Verim %")
        print("[3] SAÇILMA KESİTİ Etkisi -> Verim %")
        print("[4] Ana Menüye Dön")
            
        secim = input("\nSeçiminiz (1-4): ").strip()
        if secim not in ['1','2','3','4']:
            print("\n[HATA] Geçersiz seçim!")
            return
        elif secim == '4':
            break
        elif secim == '1':  
            analiz_yogunluk_etkisi()
        elif secim == '2':
            analiz_absorpsiyon_etkisi()
        elif secim == '3':
            analiz_sacilma_etkisi()        
def analiz_yogunluk_etkisi():
    """
    Hipotez: Yoğunluk arttıkça atomlar sıkılaşır, nötronun çarpışmadan aldığı yol (MFP) kısalır.
    """
    print("Analiz 1: Yoğunluk vs. Ortalama Serbest Yol (MFP)")
    print("Sabit: Su benzeri malzeme (H2O)")
    print("Değişken: Yoğunluk (0.1 g/cm3'ten 3.0 g/cm3'e kadar)\n")
    
    yogunluklar = np.linspace(0.1, 3.0, 30) #30 farklı nokta
    mfp_degerleri = []
    
    
    for rho in yogunluklar:
        # Geçici Malzeme Oluştur
        mat = Malzeme(f"Rho-{rho:.1f}", {H: 2, O: 1}, yogunluk=rho)
        # Ortalama Serbest Yol (Mean Free Path) = 1 / Sigma_total
        mfp = 1.0 / mat.Sigma_tot
        mfp_degerleri.append(mfp)
    
    # Grafik
    plt.figure(figsize=(10, 6))
    plt.plot(yogunluklar, mfp_degerleri, color='purple', linewidth=2.5)
    plt.title("Yoğunluğun Nötron Yoluna Etkisi", fontsize=14)
    plt.xlabel("Malzeme Yoğunluğu (g/cm³)", fontsize=12)
    plt.ylabel("Ortalama Serbest Yol (cm)", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.text(yogunluklar[-1], mfp_degerleri[-1], "Daha Yoğun = Daha Kısa Yol", ha='right', va='bottom')
    plt.show()

def analiz_absorpsiyon_etkisi():
    """
    Hipotez: Sigma_a arttıkça nötronlar termalize olamadan yutulur (Verim düşer).
    """
    print("Analiz 2: Absorpsiyon Kesiti (Sigma_a) vs. Verim")
    print("Sabit: Hidrojen benzeri saçılma (20 barn), Yoğunluk 1.0")
    print("Değişken: Sigma_a (0.01 barn'dan 2.0 barn'a kadar)\n")
    
    sig_a_degerleri = np.linspace(0.01, 2.0, 20)
    verimler = []
    
    print("Hesaplanıyor")
    
    for sig_a in sig_a_degerleri:
        # Sanal atom oluştur: Saçılması H gibi (20), Absorpsiyonu değişken
        atom_test = Atom("Test", A=1, sigma_s=20.0, sigma_a=sig_a)
        mat = Malzeme("TestMat", {atom_test: 1}, yogunluk=1.0)
        
        # Hızlı Simülasyon (N=1000 yeterli trend için)
        _, termal, _, _ = moderasyon_simulasyonu_hizli(mat, 1000) 
        verim = (termal / 1000) * 100
        verimler.append(verim)
        
    # Grafik
    plt.figure(figsize=(10, 6))
    plt.plot(sig_a_degerleri, verimler, color='red', marker='o', linestyle='-')
    plt.title("Absorpsiyon Kesitinin Moderasyon Verimine Etkisi", fontsize=14)
    plt.xlabel("Mikroskobik Absorpsiyon Kesiti(barn)", fontsize=12)
    plt.ylabel("Termalizasyon Verimi (%)", fontsize=12)
    plt.grid(True)
    plt.fill_between(sig_a_degerleri, verimler, color='red', alpha=0.1)
    plt.show()

def analiz_sacilma_etkisi():
    """
    Hipotez: Saçılma kesiti (Sigma_s) arttıkça, absorpsiyon olasılığı (P_abs) düşer, verim artar.
    Formül: P_scat = Sigma_s / (Sigma_s + Sigma_a)
    """
    print("Analiz 3: Saçılma Kesiti (Sigma_s) vs. Verim")
    print("Sabit: Absorpsiyon (0.33 barn - Su gibi), Yoğunluk 1.0")
    print("Değişken: Sigma_s (1 barn'dan 50 barn'a kadar)\n")
    
    sig_s_degerleri = np.linspace(1.0, 50.0, 25)
    verimler = []
    
    print("Hesaplanıyor")
    
    for sig_s in sig_s_degerleri:
        # Sanal atom: Absorpsiyon sabit, Saçılma değişken
        atom_test = Atom("TestScat", A=1, sigma_s=sig_s, sigma_a=0.332)
        mat = Malzeme("TestMat", {atom_test: 1}, yogunluk=1.0)
        
        # Hızlı Simülasyon
        _, termal, _, _ = moderasyon_simulasyonu_hizli(mat, 1000)
        verim = (termal / 1000) * 100
        verimler.append(verim)

    # Grafik
    plt.figure(figsize=(10, 6))
    plt.plot(sig_s_degerleri, verimler, color='blue', marker='s', linestyle='-')
    plt.title("Saçılma Kesitinin Moderasyon Verimine Etkisi", fontsize=14)
    plt.xlabel("Mikroskobik Saçılma Kesiti(barn)", fontsize=12)
    plt.ylabel("Termalizasyon Verimi (%)", fontsize=12)
    plt.grid(True)
    plt.show()
# analiz menüsü için hızlı simülasyon
def moderasyon_simulasyonu_hizli(malzeme, n_sayisi):
    termalize = 0
    absorbe = 0
    for _ in range(n_sayisi):
        E = 2000000.0
        aktif = True
        while E > 0.025 and aktif:
            if random.random() < malzeme.P_absorpsiyon:
                aktif = False
                absorbe += 1
            else:
                atom = malzeme.hedef_atom_sec()
                alpha = atom.alpha
                cos_theta = 2 * random.random() - 1
                E = 0.5 * E * ((1 + alpha) + (1 - alpha) * cos_theta)
        if aktif: termalize += 1
    return [], termalize, absorbe, 0
# Yeni madde ekleme
def yeni_malzeme_tasarla():
    print("\n" + "#"*40)
    print("      AR-GE MODU: YENİ MODERATÖR TASARIMI      ")
    print("#"*40)
    
    try:
        isim = input("Yeni Malzemenin Adı (örn: Lityum Hidrür): ")
        yogunluk = float(input(f"{isim} için Yoğunluk girin (g/cm3): "))
        
        bilesenler = {}
        print("\n--- Atomları Tanımlayın ---")
        
        while True:
            atom_sembol = input("\nAtom Sembolü (örn: Li, Be, F): ")
            A = float(input(f"  > {atom_sembol} Kütle Numarası (A): "))
            sig_s = float(input(f"  > {atom_sembol} Saçılma Kesiti (sigma_s, barn): "))
            sig_a = float(input(f"  > {atom_sembol} Absorpsiyon Kesiti (sigma_a, barn): "))
            adet = int(input(f"  > Formülde kaç adet {atom_sembol} var?: "))
            
            yeni_atom = Atom(atom_sembol, A, sig_s, sig_a)
            bilesenler[yeni_atom] = adet
            
            devam = input("  > Başka atom ekleyecek misiniz? (E/H): ").strip().upper()
            if devam != 'E':
                break
        
        print(f"\n[BAŞARILI] {isim} sisteme eklendi!")
        return Malzeme(isim, bilesenler, yogunluk)
        
    except ValueError:
        print("\n[HATA] Lütfen sayısal değerleri doğru giriniz!")
        return None

# --- ANA PROGRAM AKIŞI ---

# Standart Atomlar
H = Atom("H", 1, 20.0, 0.332)
D = Atom("D", 2, 3.4, 0.0005)
C = Atom("C", 12, 4.8, 0.0035)
O = Atom("O", 16, 3.8, 0.0002)

# Standart Malzemeler
standart_malzemeler = [
    Malzeme("Hafif Su", {H: 2, O: 1}, yogunluk=1.0),
    Malzeme("Ağır Su", {D: 2, O: 1}, yogunluk=1.105),
    Malzeme("Grafit", {C: 1}, yogunluk=1.70),
    Malzeme("Polietilen", {C: 1, H: 2}, yogunluk=0.95)
]

aktif_malzemeler = standart_malzemeler.copy()

while True:
    print("\n" + "="*50)
    print("      MODSIM: NÖTRON SİMÜLATÖRÜ      ")
    print("="*50)
    print("[1] Standart Karşılaştırma Testi (Sütun Grafikleri)")
    print("[2] Yeni Malzeme Ekle (Ar-Ge Modu)")
    print("[3] Parametrik Analiz Laboratuvarı (Çizgi Grafikleri)")
    print("[4] Çıkış")
    
    ana_secim = input("\nNe yapmak istersiniz? (1-4): ").strip()
    
    if ana_secim == '1':
    
        sonuclar = {}
        print("\nSimülasyon Başlıyor...")
        for mat in aktif_malzemeler:
            data, termal, absorbe, sure = moderasyon_simulasyonu(mat)
            ort = np.mean(data) if data else 0
            verim = (termal / SIMULASYON_NOTRON_SAYISI) * 100
            sonuclar[mat.isim] = {"Ort": ort, "Verim": verim}
            print(f"--> {mat.isim}: Ort. Çarpışma={ort:.1f}, Verim=%{verim:.1f}")
            
        # Grafik
        isimler = list(sonuclar.keys())
        degerler = [sonuclar[k]["Ort"] for k in isimler]
        verimler = [sonuclar[k]["Verim"] for k in isimler]
        
        
        renkler = ["#80d91a", "#d60e68", "#b61010", '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Grafik 1
        ax1.bar(isimler, degerler, color=renkler[:len(isimler)])
        ax1.set_title("Yavaşlatma Performansı (Düşük = İyi)")
        ax1.set_ylabel("Ortalama Çarpışma Sayısı")
        for i, v in enumerate(degerler):
            ax1.text(i, v + 1, f"{v:.1f}", ha='center')

        # Grafik 2
        ax2.bar(isimler, verimler, color=renkler[:len(isimler)])
        ax2.set_title("Nötron Ekonomisi (Yüksek = İyi)")
        ax2.set_ylabel("Termalizasyon Verimi (%)")
        
        
        if verimler:
            ax2.set_ylim(min(verimler) - 5, 105)
        else:
            ax2.set_ylim(0, 105)
            
        for i, v in enumerate(verimler):
            ax2.text(i, v + 0.2, f"%{v:.1f}", ha='center')

        plt.suptitle(f"Moderatör Karşılaştırma Analizi (N={SIMULASYON_NOTRON_SAYISI})", fontsize=14)
        plt.tight_layout()
        plt.show()

    elif ana_secim == '2':
        yeni = yeni_malzeme_tasarla()
        if yeni:
            aktif_malzemeler.append(yeni)
            print("Malzeme listeye eklendi! (1. seçenekte görebilirsiniz)")

    elif ana_secim == '3':
        analiz_modu()

    elif ana_secim == '4':
        print("Çıkış yapılıyor")
        break