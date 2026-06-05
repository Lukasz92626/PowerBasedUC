# Power-Based Unit Commitment — implementacja CPLEX/C++
 
Implementacja w C++ ścisłego sformułowania MIP przedstawionego w artykule:

> **Morales-España, G., Gentile, C., & Ramos, A. (2015)**  
> [*Tight MIP formulations of the power-based unit commitment problem.*](https://link.springer.com/article/10.1007/s00291-015-0400-4?utm_source=chatgpt.com#auth-Germ_n-Morales_Espa_a-Aff1)
> OR Spectrum, 37(4), 929–950. DOI: [10.1007/s00291-015-0400-4]

Projekt przekłada sformułowanie matematyczne z powyższego artykułu na działający solver CPLEX/C++. Numeracja równań, definicje zmiennych i ograniczenia są zgodne z oryginalną publikacją.
 
---
 
## Opis
 
Implementacja rozwiązuje problem harmonogramowania generatorów cieplnych w krótkim horyzoncie czasowym. Kluczową cechą podejścia jest **wyraźne rozróżnienie między mocą a energią**, co pozwala uniknąć problemu niedopuszczalnej dostawy energii występującego w tradycyjnych sformułowaniach opartych na energii. Model jawnie uwzględnia trajektorie mocy podczas rozruchu i odstawienia.
 
Główne cechy:
- Opis powłoki wypukłej (convex hull) dla jednostek **szybko** i **wolnorozruchowych**
- **Trajektorie mocy** podczas rozruchu i odstawienia dla jednostek wolnorozruchowych (liniowa rampa)
- **Zdolności rozruchowe i odstawienne** (SU/SD) dla jednostek szybkorozruchowych
- Ograniczenia **minimalnego czasu pracy i postoju**
- **Limity generacji** z korektą zdolności SU/SD
- **Ograniczenia rampowania** ze zintegrowaną rezerwą w górę i w dół
- **Energia** obliczana metodą trapezów z całkowitej mocy wyjściowej
- Systemowe **bilansowanie mocy** i wymagania dotyczące **rezerwy**
---
 
## Sformułowanie problemu
 
Celem jest minimalizacja całkowitych kosztów operacyjnych w horyzoncie $T$ okresów i $G$ generatorów:
 
$$\min \sum_{g,t} \left[ C^{LV}_g \, e_{g,t} + C^{NL}_g \, u_{g,t} + C^{SU*}_g \, v_{g,t} + C^{SD*}_g \, w_{g,t} \right]$$
 
gdzie $C^{SU\*}_g = C^{SU}_g + C^{NL}_g \cdot \text{SUD}_g$ oraz $C^{SD\*}_g = C^{SD}_g + C^{NL}_g \cdot \text{SDD}_g$ uwzględniają koszt bezobciążeniowy w trakcie rozruchu/odstawienia.
 
### Zmienne decyzyjne
 
| Zmienna | Typ | Opis |
|---|---|---|
| $u_{g,t}$ | Binarna | 1 jeśli jednostka $g$ jest **online** w okresie $t$ |
| $v_{g,t}$ | Binarna | 1 jeśli jednostka $g$ **startuje** w okresie $t$ |
| $w_{g,t}$ | Binarna | 1 jeśli jednostka $g$ **odstawia** w okresie $t$ |
| $p_{g,t}$ | Ciągła | Moc **powyżej minimum** |
| $\bar{p}_{g,t}$ | Ciągła | **Całkowita** moc wyjściowa (z trajektoriami SU/SD) |
| $e_{g,t}$ | Ciągła | **Energia** wyprodukowana w okresie $t$ (metoda trapezów) |
| $r^+_{g,t}$, $r^-_{g,t}$ | Ciągła | Rezerwa **w górę / w dół** |
 
### Ograniczenia (numeracja równań z artykułu)
 
| Równ. | Opis |
|---|---|
| (4)–(6) | Limity generacji z uwzględnieniem zdolności SU/SD |
| (7) | Logika rozruchu/odstawienia: $u_t - u_{t-1} = v_t - w_t$ |
| (8) | Minimalny **czas pracy** (krocząca suma rozruchów) |
| (9) | Minimalny **czas postoju** (krocząca suma odstawień) |
| (12) | Całkowita moc dla jednostek **wolnorozruchowych** |
| (14) | Całkowita moc dla jednostek **szybkorozruchowych** |
| (15) / (31) | Energia metodą trapezów: $e_t = (\bar{p}_{t-1} + \bar{p}_t)/2$ |
| (19) | **Bilans mocy** systemu |
| (20)–(21) | Wymagania dotyczące **rezerwy** |
| (27)–(28) | **Rampowanie** ze zintegrowaną rezerwą |
 
---
 
## Struktura projektu
 
```
PowerBasedUC/
├── PowerBasedUC.cpp   # Model: dane, zmienne, ograniczenia, solver, wyniki
├── CMakeLists.txt     # Konfiguracja CMake (Windows + CPLEX Community)
├── .gitignore
└── README.md
```
 
Cała logika zawarta jest w `PowerBasedUC.cpp`:
 
| Komponent | Opis |
|---|---|
| Struktura `ThermalFleet` | Kontener danych wejściowych |
| `buildFleet()` | Wbudowana instancja testowa (3 generatory, 6 okresów) |
| `makeStartupTraj()` / `makeShutdownTraj()` | Generatory liniowych trajektorii dla jednostek wolnorozruchowych |
| `main()` | Budowa modelu CPLEX, rozwiązywanie, wyświetlanie wyników |
 
---
 
## Wymagania
 
- **IBM ILOG CPLEX** (projekt skonfigurowany pod CPLEX Studio Community 2212)
- **CMake** 3.8+
- Kompilator **C++20** (zalecany MSVC na Windows)
---
 
## Budowanie
 
### Windows (Visual Studio + CMake)
 
1. Należy zainstalować [IBM ILOG CPLEX Studio Community Edition](https://www.ibm.com/products/ilog-cplex-optimization-studio) pod domyślną ścieżką:
   ```
   C:\Program Files\IBM\ILOG\CPLEX_Studio_Community2212
   ```
 
2. Następnie trzeba sklonować repozytorium:
   ```bat
   git clone https://github.com/Lukasz92626/PowerBasedUC.git
   cd PowerBasedUC
   ```
 
3. Kolejnym krokiem jest skonfigurowanie i zbudowanie za pomocą CMake:
   ```bat
   cmake -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build --config Release
   ```
 
   Jeśli CPLEX jest zainstalowany w innym miejscu, podaj ścieżkę jawnie:
   ```bat
   cmake -B build -DCPLEX_ROOT_PATH="D:/IBM/CPLEX_Studio_Community2212"
   cmake --build build --config Release
   ```
 
4. Uruchom solver:
   ```bat
   build\Release\PowerBasedUC.exe
   ```
 
---
 
## Instancja testowa
 
Wbudowana instancja (`buildFleet()`) zawiera **3 generatory** i horyzont **6 godzin**:
 
| Gen | Typ | $\underline{P}$ | $\overline{P}$ | SUD / SDD | $T^U$ / $T^D$ | $C^{LV}$ |
|---|---|---|---|---|---|---|
| 0 | Szybkorozruchowy | 120 MW | 350 MW | 1 / 1 h | 2 / 2 h | 18 $/MWh |
| 1 | Szybkorozruchowy | 80 MW | 200 MW | 1 / 1 h | 2 / 1 h | 24 $/MWh |
| 2 | Wolnorozruchowy | 40 MW | 120 MW | 2 / 2 h | 1 / 1 h | 40 $/MWh |
 
Profil obciążenia: `[180, 420, 520, 470, 260, 150]` MW  
Rezerwa w górę/dół: `[40, 50, 60, 50, 30, 20]` MW
 
---
 
## Format wyjścia
 
```
Objective = 38042.50
 
============================
Generator 0
============================
t   ON  SU  SD  p_abs  p_tot  r+   r-   E
0   0   0   0   0.00   0.00   0.00 0.00 0.00
1   1   1   0   110.00 230.00 40.00 ...
...
```
 
Legenda kolumn:
 
| Kolumna | Znaczenie |
|---|---|
| `ON` | $u_{g,t}$ — jednostka online |
| `SU` / `SD` | $v_{g,t}$ / $w_{g,t}$ — rozruch / odstawienie |
| `p_abs` | $p_{g,t}$ — moc powyżej minimum |
| `p_tot` | $\bar{p}_{g,t}$ — całkowita moc wyjściowa |
| `r+` / `r-` | $r^+_{g,t}$ / $r^-_{g,t}$ — rezerwa w górę / w dół |
| `E` | $e_{g,t}$ — energia (MWh) |
 
---
 
## Rozszerzanie modelu
 
Aby użyć własnych danych, zastąp ciało funkcji `buildFleet()`. Kluczowe pola:
 
```cpp
fleet.isQuickStart    // true = szybkorozruchowy (SU/SD), false = wolnorozruchowy (trajektorie)
fleet.startupTraj     // PSU_{g,i} — dla szybkorozruchowych zostaw {}
fleet.shutdownTraj    // PSD_{g,i} — dla szybkorozruchowych zostaw {}
fleet.initialStatus   // > 0: godziny online,  < 0: godziny offline,  0: wyłączony
fleet.load            // zapotrzebowanie systemu w każdym okresie [MW]
fleet.reserveUp       // wymagana rezerwa w górę w każdym okresie [MW]
fleet.reserveDown     // wymagana rezerwa w dół w każdym okresie [MW]
```
 
Dla jednostek wolnorozruchowych funkcje `makeStartupTraj(Pmin, SUD)` i `makeShutdownTraj(Pmin, SDD)` generują liniową rampę od 0 do $\underline{P}$ i z powrotem.
 
---
 
## Literatura
 
1. Morales-España, G., Gentile, C., & Ramos, A. (2015). Tight MIP formulations of the power-based unit commitment problem. *OR Spectrum*, 37(4), 929–950.
2. Morales-España, G., Latorre, J. M., & Ramos, A. (2013). Tight and compact MILP formulation of start-up and shut-down ramping in unit commitment. *IEEE Trans. Power Syst.*, 28(2), 1288–1296.
3. Rajan, D., & Takriti, S. (2005). Minimum up/down polytopes of the unit commitment problem with start-up costs. *IBM Research Report RC23628*.
