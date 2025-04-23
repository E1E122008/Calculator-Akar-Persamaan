# Metode Setengah Interval Calculator

Aplikasi web untuk menghitung akar persamaan menggunakan metode setengah interval (interval halving method).

## Fitur

- Input fungsi matematika dalam format Python (contoh: x**2 - 4)
- Input interval awal (x₀, x₁)
- Input toleransi (ε)
- Tampilan hasil perhitungan
- Tabel iterasi
- Visualisasi grafik fungsi dan iterasi

## Persyaratan

- Python 3.7+
- Flask
- NumPy
- SymPy

## Instalasi

1. Clone repository ini
2. Buat virtual environment (opsional tapi direkomendasikan):
   ```bash
   python -m venv venv
   source venv/bin/activate  # Linux/Mac
   venv\Scripts\activate     # Windows
   ```
3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Menjalankan Aplikasi

1. Jalankan server Flask:
   ```bash
   python app.py
   ```
2. Buka browser dan akses `http://localhost:5000`

## Cara Penggunaan

1. Masukkan fungsi matematika dalam format Python (contoh: x**2 - 4)
2. Masukkan nilai interval awal x₀ dan x₁
3. Masukkan nilai toleransi ε (default: 0.0001)
4. Klik tombol "Hitung"
5. Hasil akan ditampilkan dalam bentuk:
   - Nilai akar yang ditemukan
   - Tabel iterasi
   - Grafik fungsi

## Contoh Input

- Fungsi: x**2 - 4
- x₀: 0
- x₁: 3
- ε: 0.0001

Hasil yang diharapkan: x ≈ 2 (akar positif dari x² - 4 = 0) 