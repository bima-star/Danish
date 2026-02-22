# Dasar Analisis DEG menggunakan Bahasa Pemograman R dan Aplikasi GO dan KEGG pada dataset GSE48018

Pendahuluan

GSE48018 merupakan dataset microarray ekspresi gen yang berasal dari dalam sel darah perifer sebelum dan pada tiga titik waktu setelah pemberian vaksin influenza trivalen pada subjek pria, dan untuk menghubungkannya dengan respons antibodi terhadap vaksin. GSe48018 mencangkup 431 array mencangkup sampel sebelum pemberian vaksin influenza trivalen yang diberi label hari ke-O sebanyak 111 sampel. setelah vaksin, terdapat tiga tik waktu yaitu Sampel Hari ke-1: 110 Sampel Hari ke-3: 101 Sampel Hari ke-14: 109 Sampel. Platform yang digunakan adalah GPL6947 Ilumina HumanHT-12 V3.0 expression beadchip. Tujuan analis yaitu Mengidentifikasi gen yang mengalami perubahan ekspresi signifikan (Differentially Expressed Genes / DEG) antara individu yang divaksinasi dibandingkan dengan kontrol (baseline).

Metode

Alat yang digunakan adalah GEO2R, g:Profiler, KEGG MAPPER, bahasa pemograman R dn bahan yang digunakan adalah GSE48018 dari NCBI. parameter digunakan adalah pembagian kelompok pada dataset GSE48018 untuk group 1 adalah sampel Baseline (sebelum diberi vaksin pada hari ke-O) dan group 2 adalah sampel vaksin (setelah pemberian vaksin pada hari 1, 3 dan 14). kriteria signifikan yaitu Adjusted P-value < 0.05,Log2 Fold Change >1 atau <1. Preprocessing adalah Kontrol kualitas data, transformasi log2, dan inspeksi normalisasi. Analisis Statistik Menggunakan package limma untuk pemodelan linear. visualisasi menggunakan:

Volcano Plot: Melihat signifikansi statistik vs fold change.

Heatmap: Pola klastering 50 gen teratas.

UMAP: Reduksi dimensi untuk melihat sebaran sampel.

Hasil

Jumlah gen yang signifikan yaitu 263 gen, dengan gen up-regulated yaitu 81 dan downregulated 182. Kemudian pada hasil Go up-regulated tidak terdapat go term yang signifikan, pada down-regulated pada domain BP respon pada stres, CC pada sitosol dan MF protein binding. kemudian hasil KEGG jalur yang sering muncul adalah Herpes simplex virus 1 infection dan Influenza A dengan jumlah hits masing masing 10 pada hasil Go domain BP dengan term ID GO:0006950 (SUSD6) berfungsi dalam penekanan aktivitas dan kematian sel, domain CC pada term id GO:0005829 (SH2D4A) berperan dalam perkembangan sel T. pada domain MF (DRC12) mengatur pergeseran mikrotubulus pada aksonem bergerak.

![](images/e7e2eada8a51ca8c8e27c0e048ae86075052f568ceea0871e5684dd0b952ff6c.jpg)

![](images/1dcb2eb8d821c10f24e11b1b3609e607228478c1e145af118c0a709ae987945d.jpg)  
Gambar 1. Volcano plot pada GEO dataset GSE48018   
Gambar 2. 50 Differentially Expressed Genes (DEGs) teratas pada GEO dataset GSE48018

![](images/1a5a61f1a5ff9ab5802d22cf18efab85e1ea5cddc7cbad44dd54c7c5e7b0e15c.jpg)

![](images/1e0ef7bab4f6c8a44d2150c9c5b455b9f25bc6cdd666f802ecafcb7c4aa80665.jpg)

Gambar 3. Hasil Gene Ontology (DEGs) teratas pada GEO dataset GSE48018 down-regulated

![](images/687e89575caecce3fcb7588d309d900f8fdc25737e379dd4f3840c1295b518dc.jpg)  
Gambar 4. Hasil Gene Ontology (DEGs) teratas pada GEO dataset GSE48018 Up-regulated

![](images/e04c59055f9983f3385cb02ee5a049b91a01db00973f4bea17405e2e6c524b76.jpg)

![](images/6ace1ee1e0cc73b827ddf549dca6038291d57a28b25739791f8f56bc87f58586.jpg)  
Gambar 5. Jalur Herpes simplex virus 1 infection pada KEGG GEO dataset GSE48018

![](images/379f57b41ce2adbcd97b41073915aa5a8d364ee9e197c759e4bd6a98e625ada6.jpg)  
Gambar 6. Jalur Influenza A pada KEGG GEO dataset GSE48018

# Kesimpulan

Analisis GO menunjukkan dominansi proses respon terhadap stress, protein binding dan sitosol, sedangkan analisis KEGG menguatkannya dengan jalur Influenza A dan Herpes simplex virus 1 infectio
