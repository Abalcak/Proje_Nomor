import streamlit as st
import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time

# --- SAYFA YAPILANDIRMASI ---
st.set_page_config(
    page_title="MODSIM: Nötron Simülatörü",
    page_icon="none",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- 1. FİZİKSEL SABİTLER ---
AVOGADRO = 0.6022        # 10^24 cinsinden
E_BASLANGIC = 2000000.0  # eV
E_TERMAL = 0.025         # eV
# Simülasyon nötron sayısı artık kullanıcı tarafından değiştirilebilir

# --- SINIF TANIMLARI ---
# Not: modsim.py'dan alınan temel mantık
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
        
        for atom, adet in self.bilesenler.items():
            N_i = adet * self.N_molekul # Atomik yoğunluk
            Sigma_s_i = N_i * atom.sigma_s
            Sigma_a_i = N_i * atom.sigma_a
            
            self.Sigma_s_toplam += Sigma_s_i
            self.Sigma_a_toplam += Sigma_a_i
            
            # Seçim için listeye ekle
            self.atom_secim_listesi.append([atom, Sigma_s_i])

        self.Sigma_tot = self.Sigma_s_toplam + self.Sigma_a_toplam
        # Sıfıra bölme hatası önlemi
        if self.Sigma_tot > 0:
            self.P_absorpsiyon = self.Sigma_a_toplam / self.Sigma_tot
        else:
            self.P_absorpsiyon = 0
            
        # 4. Seçim Olasılıkları (Kümülatif)
        kümülatif = 0
        self.hedef_atom_araliklari = []
        if self.Sigma_s_toplam > 0:
            for item in self.atom_secim_listesi:
                pay = item[1] / self.Sigma_s_toplam
                kümülatif += pay
                self.hedef_atom_araliklari.append((kümülatif, item[0]))
        else:
            # Eğer hiç saçılma yoksa, varsayılan olarak ilk atomu al (Hata önleme)
            if self.atom_secim_listesi:
                 self.hedef_atom_araliklari.append((1.0, self.atom_secim_listesi[0][0]))

    def hedef_atom_sec(self):
        r = random.random()
        for esik, atom in self.hedef_atom_araliklari:
            if r <= esik: return atom
        if self.hedef_atom_araliklari:
            return self.hedef_atom_araliklari[-1][1]
        return None

# --- SİMÜLASYON FONKSİYONLARI ---
def carpisma_kinematigi(enerji, atom):
    cos_theta = 2 * random.random() - 1
    alpha = atom.alpha
    return 0.5 * enerji * ((1 + alpha) + (1 - alpha) * cos_theta)

def moderasyon_simulasyonu(malzeme, notron_sayisi, ilerleme_cubugu=None):
    if ilerleme_cubugu:
        ilerleme_cubugu.progress(0, text=f"{malzeme.isim} simüle ediliyor...")
        
    carpisma_verileri = []
    absorbe_sayisi = 0
    termalize_sayisi = 0
    ornek_iz = [] 
    
    start = time.time()
    
    # İlerleme çubuğu güncelleme adımı
    guncelleme_adim = max(1, notron_sayisi // 20)
    
    for i in range(notron_sayisi):
        E = E_BASLANGIC
        carpisma = 0
        aktif = True
        
        bu_ornek_kayit = (i == 0)
        if bu_ornek_kayit:
            ornek_iz.append(E) 
            
        while E > E_TERMAL and aktif:
            if random.random() < malzeme.P_absorpsiyon:
                aktif = False
                absorbe_sayisi += 1
            else:
                hedef = malzeme.hedef_atom_sec()
                if hedef:
                    carpisma += 1
                    E = carpisma_kinematigi(E, hedef)
                    if bu_ornek_kayit:
                        ornek_iz.append(E)
                else:
                    # Hedef atom yoksa döngüyü kır (Hata durumu)
                    aktif = False
        
        if aktif:
            carpisma_verileri.append(carpisma)
            termalize_sayisi += 1
            
        if ilerleme_cubugu and (i + 1) % guncelleme_adim == 0:
            oran = (i + 1) / notron_sayisi
            ilerleme_cubugu.progress(oran, text=f"{malzeme.isim} simüle ediliyor... (%{int(oran*100)})")
            
    sure = time.time() - start
    if ilerleme_cubugu:
        ilerleme_cubugu.empty()
        
    return carpisma_verileri, termalize_sayisi, absorbe_sayisi, sure, ornek_iz

# --- VERİ HAZIRLIĞI ---
if 'atomlar' not in st.session_state:
    st.session_state.atomlar = {
        "H": Atom("H", 1, 20.0, 0.332),
        "D": Atom("D", 2, 3.4, 0.0005),
        "C": Atom("C", 12, 4.8, 0.0035),
        "O": Atom("O", 16, 3.8, 0.0002)
    }

if 'malzemeler' not in st.session_state:
    H = st.session_state.atomlar["H"]
    D = st.session_state.atomlar["D"]
    C = st.session_state.atomlar["C"]
    O = st.session_state.atomlar["O"]
    
    st.session_state.malzemeler = [
        Malzeme("Hafif Su", {H: 2, O: 1}, yogunluk=1.0),
        Malzeme("Ağır Su", {D: 2, O: 1}, yogunluk=1.105),
        Malzeme("Grafit", {C: 1}, yogunluk=1.70),
        Malzeme("Polietilen", {C: 1, H: 2}, yogunluk=0.95)
    ]

# --- ARAYÜZ ---
st.title("MODSIM: Nötron Moderasyon Simülatörü")
st.markdown("""
Bu uygulama, farklı malzemelerin nötron yavaşlatma (moderasyon) performanslarını
Monte Carlo yöntemi ile simüle eder.
""")

# Yan Menü
with st.sidebar:
    st.header("⚙️ Kontrol Paneli")
    mod = st.radio("Çalışma Modu", ["Standart Test & Karşılaştırma", "Yeni Malzeme Ekle", "Parametrik Analiz"])

    
    st.subheader("Simülasyon Ayarları")
    sim_notron_sayisi = st.slider("Nötron Sayısı (Her malzeme için)", 100, 20000, 1000, step=100)
    
    st.info(f"Başlangıç Enerjisi: {E_BASLANGIC/1e6} MeV\nTermal Sınır: {E_TERMAL} eV")

# --- MOD 1: STANDART TEST ---
if mod == "Standart Test & Karşılaştırma":
    st.subheader("📊 Malzeme Performans Karşılaştırması")
    
    col1, col2 = st.columns([3, 1])
    with col1:
        st.write("Aşağıdaki listedeki malzemeler simüle edilecektir:")
    with col2:
        if st.button("Simülasyonu Başlat", type="primary"):
            st.session_state.run_simulation = True
            
    # Malzeme seçimi/görüntüleme
    selected_materials = st.multiselect(
        "Dahil Edilecek Malzemeler",
        [m.isim for m in st.session_state.malzemeler],
        default=[m.isim for m in st.session_state.malzemeler]
    )
    
    active_mats = [m for m in st.session_state.malzemeler if m.isim in selected_materials]

    if st.session_state.get('run_simulation', False) and active_mats:
        sonuclar = {}
        ornek_izler = {}
        progress_bar = st.progress(0, text="Hazırlanıyor...")
        
        start_total = time.time()
        
        tabs = st.tabs(["Özet Tablo", "Detaylı Grafikler"])
        
        for idx, mat in enumerate(active_mats):
            data, termal, absorbe, sure, iz = moderasyon_simulasyonu(mat, sim_notron_sayisi, progress_bar)
            
            ort = np.mean(data) if data else 0
            verim = (termal / sim_notron_sayisi) * 100
            
            sonuclar[mat.isim] = {
                "Ort. Çarpışma": ort,
                "Verim (%)": verim,
                "Absorbe (%)": (absorbe / sim_notron_sayisi) * 100,
                "Süre (sn)": sure
            }
            ornek_izler[mat.isim] = iz
        
        sim_sure = time.time() - start_total
        st.success(f"Simülasyon tamamlandı! Toplam Süre: {sim_sure:.2f} saniye")
        st.session_state.run_simulation = False # Reset

        with tabs[0]:
            df_sonuclar = pd.DataFrame.from_dict(sonuclar, orient='index')
            st.dataframe(df_sonuclar)

        with tabs[1]:
            isimler = list(sonuclar.keys())
            degerler = [sonuclar[k]["Ort. Çarpışma"] for k in isimler]
            verimler = [sonuclar[k]["Verim (%)"] for k in isimler]
            renkler = plt.cm.tab10(np.arange(len(isimler)))

            # Grafik Alanı
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
            
            # Grafik 1
            ax1.bar(isimler, degerler, color=renkler)
            ax1.set_title("Ortalama Çarpışma Sayısı (Daha düşük = Daha hızlı yavaşlatma)")
            ax1.set_ylabel("Adet")
            ax1.grid(axis='y', alpha=0.3)
            
            # Grafik 2
            bars = ax2.bar(isimler, verimler, color=renkler)
            ax2.set_title("Termalizasyon Verimi (%)")
            ax2.set_ylabel("%")
            ax2.set_ylim(0, 105)
            ax2.grid(axis='y', alpha=0.3)
            
            st.pyplot(fig)
            
            # Enerji Kaybı Grafiği
            st.markdown("#### 📉 Tekil Nötron Enerji Kaybı (Logaritmik)")
            fig2, ax3 = plt.subplots(figsize=(10, 4))
            for i, isim in enumerate(isimler):
                iz = ornek_izler[isim]
                ax3.plot(iz, label=isim, linewidth=2)
            
            ax3.set_yscale('log')
            ax3.set_title("Enerji Düşüş Seyri")
            ax3.set_xlabel("Çarpışma Sayısı")
            ax3.set_ylabel("Enerji (eV)")
            ax3.axhline(y=0.025, color='red', linestyle=':', label='Termal Sınır')
            ax3.legend()
            ax3.grid(True, which="both", alpha=0.3)
            st.pyplot(fig2)

# --- MOD 2: AR-GE MODU ---
elif mod == "Yeni Malzeme Ekle":
    st.subheader("Yeni Moderatör Malzemesi Ekle")
    
    with st.form("yeni_malzeme_form"):
        col1, col2 = st.columns(2)
        with col1:
            yeni_isim = st.text_input("Malzeme Adı", placeholder="Örn: Lityum Hidrür")
            yeni_yogunluk = st.number_input("Yoğunluk (g/cm³)", min_value=0.1, value=1.0, step=0.1)
        
        st.markdown("##### Bileşen Atomlar")
        
        # Basitlik için en fazla 2 farklı atom türü eklemeye izin verelim şimdilik
        # Daha dinamik bir yapı için session_state kullanılabilir ama form içinde zor.
        
        c1, c2, c3, c4 = st.columns(4)
        with c1:
            atom1_sembol = st.selectbox("1. Atom", list(st.session_state.atomlar.keys()), key="a1")
        with c2:
            atom1_adet = st.number_input("Adet", 1, 100, 1, key="ad1")
        with c3:
            atom2_sembol = st.selectbox("2. Atom (Opsiyonel)", ["Yok"] + list(st.session_state.atomlar.keys()), key="a2")
        with c4:
            atom2_adet = st.number_input("Adet", 0, 100, 0, key="ad2")

        sunmitted = st.form_submit_button("Malzemeyi Ekle")
        
        if sunmitted:
            if not yeni_isim:
                st.error("Lütfen bir malzeme adı girin.")
            else:
                bilesenler = {}
                # Atom 1
                atom_obj1 = st.session_state.atomlar[atom1_sembol]
                bilesenler[atom_obj1] = atom1_adet
                
                # Atom 2
                if atom2_sembol != "Yok" and atom2_adet > 0:
                    atom_obj2 = st.session_state.atomlar[atom2_sembol]
                    bilesenler[atom_obj2] = atom2_adet
                
                yeni_mat = Malzeme(yeni_isim, bilesenler, yeni_yogunluk)
                st.session_state.malzemeler.append(yeni_mat)
                st.success(f"✅ {yeni_isim} başarıyla eklendi! 'Standart Test' modunda kullanabilirsiniz.")

    st.info("💡 Not: Şu anda sadece sistemde tanımlı atomları (H, D, C, O) kullanabilirsiniz.")

# --- MOD 3: PARAMETRİK ANALİZ ---
elif mod == "Parametrik Analiz":
    st.subheader("Parametrik Etki Analizi")
    
    analiz_tipi = st.selectbox("Analiz Parametresi", [
        "Yoğunluk Etkisi (Ortalama Serbest Yol)",
        "Absorpsiyon Kesiti Etkisi (Verim)",
        "Saçılma Kesiti Etkisi (Verim)"
    ])
    
    if analiz_tipi == "Yoğunluk Etkisi (Ortalama Serbest Yol)":
        st.markdown("**Hipotez:** Yoğunluk arttıkça atomlar sıkışır, nötronun serbest yolu kısalır.")
        
        rho_min, rho_max = st.slider("Yoğunluk Aralığı (g/cm³)", 0.1, 5.0, (0.1, 3.0))
        
        if st.button("Analizi Çalıştır"):
            H = st.session_state.atomlar["H"]
            O = st.session_state.atomlar["O"]
            yogunluklar = np.linspace(rho_min, rho_max, 30)
            mfp_degerleri = []
            
            for rho in yogunluklar:
                mat = Malzeme(f"Rho-{rho:.1f}", {H: 2, O: 1}, yogunluk=rho)
                if mat.Sigma_tot > 0:
                    mfp = 1.0 / mat.Sigma_tot
                else:
                    mfp = 0
                mfp_degerleri.append(mfp)
            
            fig, ax = plt.subplots()
            ax.plot(yogunluklar, mfp_degerleri, color='purple', linewidth=2.5)
            ax.set_title("Yoğunluk vs Ortalama Serbest Yol")
            ax.set_xlabel("Yoğunluk (g/cm³)")
            ax.set_ylabel("MFP (cm)")
            ax.grid(True)
            st.pyplot(fig)
            
    elif analiz_tipi == "Absorpsiyon Kesiti Etkisi (Verim)":
        st.markdown("**Hipotez:** Absorpsiyon kesiti arttıkça verim düşer.")
        if st.button("Analizi Çalıştır"):
            sig_a_degerleri = np.linspace(0.01, 2.0, 15)
            verimler = []
            
            bar = st.progress(0)
            for i, sig_a in enumerate(sig_a_degerleri):
                atom_test = Atom("Test", A=1, sigma_s=20.0, sigma_a=sig_a)
                mat = Malzeme("TestMat", {atom_test: 1}, yogunluk=1.0)
                # Hızlı simülasyon (az nötronlu)
                _, termal, _, _, _ = moderasyon_simulasyonu(mat, 500, None)
                verimler.append((termal/500)*100)
                bar.progress((i+1)/len(sig_a_degerleri))
            
            bar.empty()
            fig, ax = plt.subplots()
            ax.plot(sig_a_degerleri, verimler, 'r-o')
            ax.set_title("Absorpsiyon Kesiti vs Verim")
            ax.set_xlabel("Sigma_a (barn)")
            ax.set_ylabel("Verim (%)")
            ax.grid(True)
            st.pyplot(fig)

    elif analiz_tipi == "Saçılma Kesiti Etkisi (Verim)":
        st.markdown("**Hipotez:** Saçılma kesiti arttıkça verim artar.")
        if st.button("Analizi Çalıştır"):
            sig_s_degerleri = np.linspace(1.0, 50.0, 15)
            verimler = []
            
            bar = st.progress(0)
            for i, sig_s in enumerate(sig_s_degerleri):
                atom_test = Atom("Test", A=1, sigma_s=sig_s, sigma_a=0.332)
                mat = Malzeme("TestMat", {atom_test: 1}, yogunluk=1.0)
                _, termal, _, _, _ = moderasyon_simulasyonu(mat, 500, None)
                verimler.append((termal/500)*100)
                bar.progress((i+1)/len(sig_s_degerleri))
            
            bar.empty()
            fig, ax = plt.subplots()
            ax.plot(sig_s_degerleri, verimler, 'b-s')
            ax.set_title("Saçılma Kesiti vs Verim")
            ax.set_xlabel("Sigma_s (barn)")
            ax.set_ylabel("Verim (%)")
            ax.grid(True)
            st.pyplot(fig)
