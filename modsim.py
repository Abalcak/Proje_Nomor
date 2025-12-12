import random
import numpy as np
import matplotlib.pyplot as plt
import time

# --- 1. FİZİKSEL SABİTLER ---
AVOGADRO = 0.6022        # 10^24 cinsinden (Barn birimiyle uyumlu)
E_BASLANGIC = 2000000.0  # eV
E_TERMAL = 0.025         # eV
SIMULASYON_NOTRON_SAYISI = 100000 # Her malzeme için deneme sayısı

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

# --- YENİ MADDE EKLEME SİHİRBAZI ---
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

# 1. Standart Atomlar
H = Atom("H", 1, 20.0, 0.332)
D = Atom("D", 2, 3.4, 0.0005)
C = Atom("C", 12, 4.8, 0.0035)
O = Atom("O", 16, 3.8, 0.0002)

# 2. Standart Malzemeler
malzeme_listesi = [
    Malzeme("Hafif Su", {H: 2, O: 1}, yogunluk=1.0),
    Malzeme("Ağır Su", {D: 2, O: 1}, yogunluk=1.105),
    Malzeme("Grafit", {C: 1}, yogunluk=1.70),
    Malzeme("Polietilen", {C: 1, H: 2}, yogunluk=0.95),
    Malzeme("MESİTİLEN", {C: 9,H:12}, yogunluk=0.864)                                                                        
    
]

# 3. KULLANICIYA SOR: Yeni madde eklensin mi?
print("\n" + "="*50)
print("MODSIM: EĞİTİM AMAÇLI NÖTRON SİMÜLASYONU")
print("="*50)

secim = input("Kendi tasarladığınız bir malzemeyi test etmek ister misiniz? (E/H): ").strip().upper()

if secim == 'E':
    ozel_malzeme = yeni_malzeme_tasarla()
    if ozel_malzeme is not None:
        malzeme_listesi.append(ozel_malzeme)

# 4. Simülasyonu Çalıştır
sonuclar = {}

print("\n" + "="*50)
print(f"SİMÜLASYON BAŞLIYOR (Her malzeme için {SIMULASYON_NOTRON_SAYISI} Nötron)")
print("="*50)

for mat in malzeme_listesi:
    data, termal, absorbe, sure = moderasyon_simulasyonu(mat)
    
    ort = np.mean(data) if data else 0
    verim = (termal / SIMULASYON_NOTRON_SAYISI) * 100
    
    sonuclar[mat.isim] = {"Ort": ort, "Verim": verim}
    print(f"--> SONUÇ: Ort. Çarpışma: {ort:.1f} | Verim: %{verim:.1f} | Süre: {sure:.2f}s")

# 5. GRAFİK
isimler = list(sonuclar.keys())
degerler = [sonuclar[k]["Ort"] for k in isimler]
verimler = [sonuclar[k]["Verim"] for k in isimler]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Renk paleti (Standartlar + Yeni eklenen için Kırmızı)
renkler = ['#1f77b4', '#17becf', '#7f7f7f', '#2ca02c', '#d62728', '#9467bd']

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
ax2.set_ylim(min(verimler)-5, 105) 
for i, v in enumerate(verimler):
    ax2.text(i, v + 0.2, f"%{v:.1f}", ha='center')

plt.suptitle("ModSim: Standart ve Yenilikçi Moderatör Analizi", fontsize=16)
plt.tight_layout()
plt.show()