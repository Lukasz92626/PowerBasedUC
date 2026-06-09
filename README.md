# Power-Based Unit Commitment — implementacja CPLEX/C++

Implementacja w C++ ścisłego sformułowania MIP przedstawionego w artykule:

> **Germán Morales-España, Claudio Gentile & Andres Ramos (2015)**
> [*Tight MIP formulations of the power-based unit commitment problem.*](https://link.springer.com/article/10.1007/s00291-015-0400-4)
> OR Spectrum, 37(4), 929–950. DOI: [10.1007/s00291-015-0400-4](https://doi.org/10.1007/s00291-015-0400-4)

Projekt przekłada sformułowanie matematyczne z powyższego artykułu na działający solver CPLEX/C++. Numeracja równań, definicje zmiennych i ograniczenia są zgodne z oryginalną publikacją.

---

## Spis treści

1. [Problem Unit Commitment](#1-problem-unit-commitment)
2. [Moc a energia — kluczowe rozróżnienie](#2-moc-a-energia--kluczowe-rozróżnienie)
3. [Opis implementacji](#3-opis-implementacji)
4. [Sformułowanie matematyczne](#4-sformułowanie-matematyczne)
5. [Powłoka wypukła i przewaga nad tradycyjnymi sformułowaniami](#5-powłoka-wypukła-i-przewaga-nad-tradycyjnymi-sformułowaniami)
6. [Struktura projektu](#6-struktura-projektu)
7. [Wymagania](#7-wymagania)
8. [Budowanie](#8-budowanie)
9. [Instancja testowa](#9-instancja-testowa)
10. [Format wyjścia](#10-format-wyjścia)
11. [Rozszerzanie modelu](#11-rozszerzanie-modelu)

---

## 1. Problem Unit Commitment

**Unit Commitment (UC)** to krótkoterminowy problem harmonogramowania jednostek wytwórczych, rozwiązywany codziennie przez operatorów systemów elektroenergetycznych oraz uczestników rynku energii. Horyzont planowania wynosi typowo od jednej doby do jednego tygodnia, z rozdzielczością godzinową.

Formalnie UC można zapisać jako:

$$\min_{x,\, p} \quad b^\top x + c^\top p$$

$$\text{s.t.} \quad Fx \leq f, \quad x \in \{0,1\}^n$$

$$H p \leq h$$

$$Ax + Bp \leq g$$

gdzie $x$ to wektor binarnych decyzji zaangażowania (on/off, rozruch, odstawienie), a $p$ to wektor ciągłych decyzji dyspozycji mocy. Funkcja celu minimalizuje sumę kosztów zaangażowania $b^\top x$ (koszty bezobciążeniowe, rozruchu i odstawienia) oraz kosztów produkcji $c^\top p$.

Problem UC należy do klasy **NP-trudnych** — liczba możliwych kombinacji zaangażowania rośnie wykładniczo z liczbą jednostek i liczbą okresów. W praktyce rozwiązuje się go metodą **Mixed-Integer Programming (MIP)**, korzystając z solvera takiego jak IBM ILOG CPLEX. Kluczowym czynnikiem wpływającym na czas obliczeń jest **ścisłość sformułowania** — im bliższe optimum LP-relaksacji jest optimum MIP, tym mniej węzłów branch-and-bound trzeba zbadać.

---

## 2. Moc a energia — kluczowe rozróżnienie

Tradycyjne sformułowania UC opierają się na **harmonogramowaniu energii**: zmienna decyzyjna reprezentuje energię wyprodukowaną w danej godzinie, traktowaną jako stały poziom mocy przez cały okres. Takie podejście prowadzi do **niedopuszczalnych harmonogramów**, których nie da się fizycznie zrealizować.

### Przykład niedopuszczalności

Rozważmy jednostkę z $\underline{P} = 100\ \text{MW}$, $\overline{P} = 300\ \text{MW}$ i prędkością rampowania $200\ \text{MW/h}$.

**Tradycyjny harmonogram energii** (niedopuszczalny fizycznie):

```
Godzina 1:  100 MW  (jednostka przy minimum)
Godzina 2:  300 MW  (skok do maksimum — wymaga natychmiastowej zmiany)
```

Jednostka wytwarza 100 MW na końcu godziny 1. Żeby wyprodukować 300 MW w godzinie 2, musiałaby osiągnąć maksimum **natychmiast** na początku godziny 2 — co wymagałoby nieskończonej prędkości rampowania.

**Sformułowanie oparte na mocy** (poprawne):

```
Koniec godziny 1:  100 MW
Koniec godziny 2:  300 MW  ← osiągalne dopiero pod koniec godziny 2
Energia w godz. 2: (100 + 300) / 2 = 200 MWh  ← metoda trapezów
```

W podejściu proponowanym przez Morales-España et al. zmienne decyzyjne reprezentują **moc na końcu każdego okresu**, a energia jest obliczana jako pole trapezu między dwoma kolejnymi wartościami mocy:

$$e_{g,t} = \frac{\bar{p}_{g,t-1} + \bar{p}_{g,t}}{2}$$

Pozwala to dokładnie modelować fizyczne ograniczenia rampowania i trajektorie mocy podczas rozruchu i odstawienia, których tradycyjne podejście po prostu ignoruje.

### Trajektorie rozruchu i odstawienia

Dodatkowym problemem tradycyjnych sformułowań jest założenie, że jednostki **zaczynają i kończą pracę dokładnie przy mocy minimalnej** $\underline{P}$. W rzeczywistości proces rozruchu i odstawienia przebiega według z góry określonej trajektorii mocy:

- **Jednostki wolnorozruchowe** — wymagają więcej niż jednego okresu na dojście od 0 do $\underline{P}$ (rozruch) lub zejście z $\underline{P}$ do 0 (odstawienie). Trajektoria $P^{SU}_i$ / $P^{SD}_i$ jest predefiniowana.
- **Jednostki szybkorozruchowe** — mogą osiągnąć dowolną wartość między $\underline{P}$ a $SU$ w ciągu jednego okresu. Parametry $SU$ i $SD$ definiują zdolność rozruchową i odstawienną ($SU, SD \geq \underline{P}$).

Energia wyprodukowana podczas rozruchu/odstawienia jest realnie obecna w systemie, ale tradycyjne sformułowania jej nie alokują — co powoduje błędy w rozliczeniach rynkowych i wymaga większych rezerw operacyjnych.

---

## 3. Opis implementacji

Implementacja rozwiązuje problem harmonogramowania generatorów cieplnych w krótkim horyzoncie czasowym. Kluczową cechą podejścia jest **wyraźne rozróżnienie między mocą a energią**, co pozwala uniknąć problemu niedopuszczalnej dostawy energii. Model jawnie uwzględnia trajektorie mocy podczas rozruchu i odstawienia.

Główne cechy:

- Opis **powłoki wypukłej** (convex hull) dla jednostek szybko- i wolnorozruchowych
- **Trajektorie mocy** podczas rozruchu i odstawienia dla jednostek wolnorozruchowych (liniowa rampa)
- **Zdolności rozruchowe i odstawienne** (SU/SD) dla jednostek szybkorozruchowych
- Ograniczenia **minimalnego czasu pracy i postoju**
- **Limity generacji** z korektą zdolności SU/SD
- **Ograniczenia rampowania** ze zintegrowaną rezerwą w górę i w dół
- **Energia** obliczana metodą trapezów z całkowitej mocy wyjściowej
- Systemowe **bilansowanie mocy** i wymagania dotyczące **rezerwy**

---

## 4. Sformułowanie matematyczne

### Funkcja celu

Minimalizacja całkowitych kosztów operacyjnych w horyzoncie $T$ okresów i $G$ generatorów (równanie (18) z artykułu):

$$\min \sum_{g \in G} \sum_{t \in T} \left[ C^{LV}_g \, e_{g,t} + C^{NL}_g \, u_{g,t} + C^{SU*}_g \, v_{g,t} + C^{SD*}_g \, w_{g,t} \right]$$

Efektywne koszty rozruchu i odstawienia uwzględniają koszt bezobciążeniowy w trakcie tych procesów (równania (18a)–(18b)):

$$C^{SU*}_g = C^{SU}_g + C^{NL}_g \cdot \text{SUD}_g, \qquad C^{SD*}_g = C^{SD}_g + C^{NL}_g \cdot \text{SDD}_g$$

### Zmienne decyzyjne

| Zmienna | Typ | Opis |
|---|---|---|
| $u_{g,t}$ | Binarna | 1 jeśli jednostka $g$ jest **online** (powyżej $\underline{P}$) w okresie $t$ |
| $v_{g,t}$ | Binarna | 1 jeśli jednostka $g$ **startuje** w okresie $t$ |
| $w_{g,t}$ | Binarna | 1 jeśli jednostka $g$ **odstawia** w okresie $t$ |
| $p_{g,t}$ | Ciągła ≥ 0 | Moc **powyżej minimum** $\underline{P}$ |
| $\bar{p}_{g,t}$ | Ciągła ≥ 0 | **Całkowita** moc wyjściowa (z trajektoriami SU/SD) |
| $e_{g,t}$ | Ciągła ≥ 0 | **Energia** wyprodukowana w okresie $t$ |
| $r^+_{g,t}$, $r^-_{g,t}$ | Ciągła ≥ 0 | Rezerwa **w górę / w dół** |

### Ograniczenia

#### Logika zaangażowania (równanie (7))

$$u_{g,t} - u_{g,t-1} = v_{g,t} - w_{g,t} \qquad \forall g,\, t$$

#### Minimalny czas pracy — convex hull (równanie (8))

$$\sum_{i=t-T^U_g+1}^{t} v_{g,i} \leq u_{g,t} \qquad \forall g,\, t \in [T^U_g,\, T]$$

#### Minimalny czas postoju — convex hull (równanie (9))

$$\sum_{i=t-T^D_g+1}^{t} w_{g,i} \leq 1 - u_{g,t} \qquad \forall g,\, t \in [T^D_g,\, T]$$

#### Limity generacji z uwzględnieniem SU/SD (równania (4)–(5))

$$p_{g,t} + r^+_{g,t} \leq (\overline{P}_g - \underline{P}_g)\, u_{g,t} - (\overline{P}_g - SD_g)\, w_{g,t+1} + (SU_g - \underline{P}_g)\, v_{g,t+1} \qquad \forall g,\, t < T$$

$$p_{g,T} + r^+_{g,T} \leq (\overline{P}_g - \underline{P}_g)\, u_{g,T} \qquad \forall g$$

Kluczowa własność: jeśli jednostka odstawia w następnym okresie ($w_{t+1}=1$), dostępny zakres mocy zmniejsza się, bo jednostka musi zejść do poziomu $SD_g < \overline{P}_g$. Analogicznie, jeśli startuje ($v_{t+1}=1$), może tymczasowo osiągnąć $SU_g > \underline{P}_g$.

#### Całkowita moc — jednostki szybkorozruchowe (równanie (14))

$$\bar{p}_{g,t} = \underline{P}_g \left( u_{g,t} + v_{g,t+1} \right) + p_{g,t} \qquad \forall g \in G^Q,\, t$$

#### Całkowita moc — jednostki wolnorozruchowe (równanie (12))

$$\bar{p}_{g,t} = \underline{P}_g \left( u_{g,t} + v_{g,t+1} \right) + p_{g,t} + \sum_{i=1}^{SUD_g} P^{SU}_{g,i}\, v_{g,\,t-i+SUD_g+2} + \sum_{i=2}^{SDD_g+1} P^{SD}_{g,i}\, w_{g,\,t-i+2} \qquad \forall g \in G^S,\, t$$

Sumy po $i$ dodają wkład trajektorii rozruchu i odstawienia w okresy, w których jednostka jeszcze nie jest (lub już nie jest) w stanie online.

#### Energia — metoda trapezów (równania (15), (31))

$$e_{g,t} = \frac{\bar{p}_{g,t-1} + \bar{p}_{g,t}}{2} \qquad \forall g,\, t$$

#### Rampowanie z rezerwą (równania (27)–(28))

$$p_{g,t} + r^+_{g,t} - p_{g,t-1} \leq RU_g \cdot u_{g,t-1} + (SU_g - \underline{P}_g)\, v_{g,t} \qquad \forall g,\, t$$

$$p_{g,t-1} - \left( p_{g,t} - r^-_{g,t} \right) \leq RD_g \cdot u_{g,t} + (SD_g - \underline{P}_g)\, w_{g,t} \qquad \forall g,\, t$$

Ponieważ $p_{g,t}$ mierzy moc **powyżej** $\underline{P}_g$, ograniczenia rampowania nie wymagają dodatkowych zmiennych binarnych ani parametrów big-M — w odróżnieniu od tradycyjnych sformułowań.

#### Bilans mocy systemu (równanie (19))

$$\sum_{g \in G} \bar{p}_{g,t} \geq D_t \qquad \forall t$$

#### Wymagania dotyczące rezerwy (równania (20)–(21))

$$\sum_{g \in G} r^+_{g,t} \geq D^+_t, \qquad \sum_{g \in G} r^-_{g,t} \geq D^-_t \qquad \forall t$$

---

## 5. Powłoka wypukła i przewaga nad tradycyjnymi sformułowaniami

### Convex hull

Morales-España et al. dowodzą, że układ ograniczeń (4)–(15) stanowi **opis powłoki wypukłej** zbioru dopuszczalnych rozwiązań dla pojedynczej jednostki (Twierdzenia 1–5 w artykule). Oznacza to, że LP-relaksacja tego sformułowania daje automatycznie całkowite rozwiązanie dla problemu self-UC — MIP rozwiązywany jest w czasie LP, bez konieczności branch-and-bound.

Formalnie, dla zbioru:

$$D_T = \{(u, v, p) \in \mathbb{R}^{3T-1}_+ \mid (u,v,p)\ \text{spełnia}\ (4)\text{–}(11)\}$$

zachodzi $C_T = D_T$, gdzie $C_T$ to powłoka wypukła punktów całkowitoliczbowych. Wynik ten rozszerza wcześniejsze prace dotyczące convex hull dla minimalnych czasów pracy/postoju (Rajan & Takriti 2005) o ograniczenia generacji i zdolności SU/SD.

### Porównanie z tradycyjnymi sformułowaniami

W literaturze funkcjonują dwa powszechnie stosowane sformułowania energetyczne:

- **1bin** (Carrión & Arroyo 2006) — jedno binarne na jednostkę i okres; rozruch/odstawienie wyrażone przez zaangażowanie
- **3bin** (FERC 2012) — trzy binarne zmienne ($u$, $v$, $w$) z convex hull dla minimalnych czasów pracy/postoju wg Rajan & Takriti

Wyniki obliczeniowe z artykułu (Table 4 i Table 6) pokazują dramatyczną przewagę sformułowania **Pw** (power-based):

| Sformułowanie | IntGap (self-UC) | RootIGap (IEEE 118) | Czas MIP (64 dni) |
|---|---|---|---|
| **Pw** (niniejsza implementacja) | **0 %** | **0.089 %** | LP time |
| 3bin | 0.88 % | 0.375 % | 12.01 s |
| 1bin | 2.57 % | 1.052 % | 13.79 s |

Dla problemu self-UC sformułowanie Pw osiąga IntGap = 0 — MIP jest rozwiązywany jako LP bez branch-and-bound. Dla sieci IEEE 118-bus z 54 jednostkami i 60-godzinnym horyzontem, Pw jako jedyne rozwiązuje wszystkie przypadki w limicie czasu (10 000 s), podczas gdy 3bin i 1bin nie są w stanie ukończyć obliczeń.

Przewaga wynika z dwóch własności jednocześnie:
- **Ścisłość** — mniejszy RootIGap oznacza, że solver potrzebuje mniej węzłów B&B
- **Kompaktowość** — Pw ma mniej ograniczeń niż 3bin i 1bin (patrz Table 5 w artykule)

---

## 6. Struktura projektu

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
| Struktura `ThermalFleet` | Kontener danych wejściowych dla całej floty |
| `buildFleet()` | Wbudowana instancja testowa (3 generatory, 6 okresów) |
| `makeStartupTraj()` / `makeShutdownTraj()` | Generatory liniowych trajektorii dla jednostek wolnorozruchowych |
| `main()` | Budowa modelu CPLEX, rozwiązywanie, wyświetlanie wyników |

---

## 7. Wymagania

- **IBM ILOG CPLEX** (projekt skonfigurowany pod CPLEX Studio Community 2212)
- **CMake** 3.8+
- Kompilator **C++20** (zalecany MSVC na Windows)

---

## 8. Budowanie

### Windows (Visual Studio + CMake)

1. Zainstaluj [IBM ILOG CPLEX Studio Community Edition](https://www.ibm.com/products/ilog-cplex-optimization-studio) pod domyślną ścieżką:
   ```
   C:\Program Files\IBM\ILOG\CPLEX_Studio_Community2212
   ```

2. Sklonuj repozytorium:
   ```bat
   git clone https://github.com/Lukasz92626/PowerBasedUC.git
   cd PowerBasedUC
   ```

3. Skonfiguruj i zbuduj za pomocą CMake:
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

## 9. Instancja testowa

Wbudowana instancja (`buildFleet()`) zawiera **3 generatory** i horyzont **6 godzin**:

| Gen | Typ | $\underline{P}$ | $\overline{P}$ | SUD / SDD | $T^U$ / $T^D$ | $C^{LV}$ | $C^{NL}$ |
|---|---|---|---|---|---|---|---|
| 0 | Szybkorozruchowy | 120 MW | 350 MW | 1 / 1 h | 2 / 2 h | 18 $/MWh | 120 $/h |
| 1 | Szybkorozruchowy | 80 MW | 200 MW | 1 / 1 h | 2 / 1 h | 24 $/MWh | 80 $/h |
| 2 | Wolnorozruchowy | 40 MW | 120 MW | 2 / 2 h | 1 / 1 h | 40 $/MWh | 30 $/h |

Profil obciążenia: `[180, 420, 520, 470, 260, 150]` MW
Rezerwa w górę/dół: `[40, 50, 60, 50, 30, 20]` MW

---

## 10. Format wyjścia

```
========================================================================
  SOLVER SUMMARY
========================================================================
  Optimal objective :  47320.00 $
  Best bound        :  47320.00 $
  MIP gap           :  0.00 %
  Solve time        :  0.02 s
  B&B nodes         :  0
  Constraints       :  195
  Variables         :  144
========================================================================

========================================================================
  GENERATOR 0  [quick-start]   P_min=120.00 MW  P_max=350.00 MW  c_var=18.00 $/MWh
========================================================================
t   status         p_abs[MW] p_tot[MW]   r+[MW]   r-[MW]    E[MWh]
------------------------------------------------------------------------
0   ON+SU              60.00    180.00     0.00    60.00     90.00
1   ON                180.00    300.00     0.00    50.00    240.00
...
------------------------------------------------------------------------
  Cost breakdown:  no-load=720.00  startup=720.00  shutdown=0.00  variable=23940.00  total=25380.00 $
========================================================================

========================================================================
  SYSTEM BALANCE CHECK
========================================================================
t   demand[MW]    gen[MW]   slack    r+req   sum_r+    r-req   sum_r-
------------------------------------------------------------------------
0       180.00     220.00      OK    40.00    40.00    40.00    60.00
...
------------------------------------------------------------------------
  Overall: ALL CONSTRAINTS SATISFIED
========================================================================
```

Legenda kolumn tabeli generatora:

| Kolumna | Zmienna | Znaczenie |
|---|---|---|
| `status` | $u_{g,t}$, $v_{g,t}$, $w_{g,t}$ | Stan jednostki: `ON`, `OFF`, `ON+SU`, `ON+SD` |
| `p_abs` | $p_{g,t}$ | Moc powyżej minimum $\underline{P}$ [MW] |
| `p_tot` | $\bar{p}_{g,t}$ | Całkowita moc wyjściowa (z trajektoriami) [MW] |
| `r+` / `r-` | $r^+_{g,t}$ / $r^-_{g,t}$ | Rezerwa w górę / w dół [MW] |
| `E` | $e_{g,t}$ | Energia wyprodukowana w okresie $t$ [MWh] |
