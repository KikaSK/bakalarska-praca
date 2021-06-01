# Bakalárska práca
## Návod na inštaláciu

Algoritmus je implementovaný v jazyku C++, používala som prostredie Visual Studio Code, 
voľne dostupné na [stiahnutie](https://code.visualstudio.com/).
Ako prvé je potrebné si nainštalovať knižnicu `GiNaC`. Návod na inštaláciu môžete nájsť
[TU](https://www.ginac.de/Download.html). Odporúčam [verziu 1.8.0.](https://www.ginac.de/ginac-1.8.0.tar.bz2). Na funkčnosť GiNaC-u je potrebné stiahnutie aj knižnice [CLN](https://www.ginac.de/CLN/), preto je potrebné pozrieť si návod na inštaláciu GiNaC.

## Návod pre OS Linux na stiahnutue a spustenie programu

Samotný algoritmus môžete nájsť na githube na [tomto](https://github.com/KikaSK/bakalarska-praca/tree/refactor) odkaze. Je potrebné si naklonovať repozitár:

1. V priečinku kam chcete repozitár naklonovať otvorte terminál.
2. Do terminálu postupne skopírujte nasledujúce príkazy:
    - `git clone https://github.com/KikaSK/bakalarska-praca`
    - `cd bakalarska-praca`
3. Do terminálu môžete postupne napísať príkazy:
    - `cd Kód/bakalarka`
    - `make`

Po poslednom príkaze sa spustia predvolené vstupné súbory, 
v priečinku "bakalarska-praca/Kód/bakalarka/outputs/" sa zobrazia vytvorené .obj súbory, ktoré si môžete pozrieť napríklad na [tejto](https://3dviewer.net/) stránke.

## Návod na spustenie iných ako predvolených súborov

Pre spustenie iných ako predvolených vstupných súborov je najpríjemnejší spôsob zatvoriť terminál, a otvoriť celý priečinok `bakalarka` vo Visual Studio Code (alebo podobnom prostredí).

V priečinku `bakalarka/inputs` sa nachádzajú vstupné súbory:
1. V podpriečinku `finite_surfaces` sú všetky končné, uzavreté plochy, ktoré sme prezentovali v BP aj s rôznymi dĺžkami hrán.
2. V podpričinku `infinite_surfaces` sú všetky nekonečné ohraničené plochy, ktoré sme prezentovali v BP aj s rôznymi dĺžkami hrán.
3. V podpriečinku `cut_surfaces` sú niektoré konečné aj nekonečné plochy s rôznym ohraničením.

V konkrétnom priečinku nájdeme vstupné súbory s názvom `input0`, `input1`, ...

#### Funkcie spúšťajúce vstupné súbory

V prostredí si otvoríme súbor `main.cpp` a približne v 3/4 súboru sa nachádza funkcia `main()`, v nej sú 
napísané 2 riadky spúšťajúce vstupné súbory aj s komentármi.
Prvý riadok vyzerá nasledovne:
    `run_input(0, "/finite_surfaces/sphere", "my_run_input");`

Funkcia `run_input` má 3 vstupné parametre:
- Druhý vstupný parameter udáva priečinok v ktorom bude funkcia hľadať vstupný súbor.
- Prvý parameter označuje číslo vstupného súboru, ktorý spustí. Čím menšie číslo vstupného súboru, tým väčšia je veľkosť hrany, teda tým kratšie trvá výpočet. Toto platí pre všetky priečinky okrem `cut_surfaces`. Odporúčam teda najprv skúšať spúšťať súbory `input0`. Algoritmus nie je optimalizovaný na rýchlosť, teda pre veľmi malé veľkosti hrán môže počítať pomerne dlho. 
- Tretí parameter určuje predponu, ktorou nazve výstupný súbor.

Druhý riadok vyzerá nasledovne:
    `run_all(0, 2, "/infinite_surfaces/hyperboloid", "my_run_all");`

Funkcia `run_all` má 4 vstupné parametre:
- Tretí vstupný parameter udáva priečinok v ktorom bude funkcia hľadať vstupný súbor.
- Prvý vstupný parameter označuje začiatočné číslo vstupného súboru, ktorý spustí. 
- Druhý vstupný parameter označuje koncové číslo vstupného súboru, ktorý spustí.
- Štvrtý parameter určuje predponu, ktorou nazve výstupné súbory.

Funkcia run_all slúži na spustenie viacerých súborov za sebou. Medzi začiatočným a koncovým 
číslom spustí všetky súbory. Teda ak je začiatočné číslo 0 a koncové 2, tak spustí súbory
`input0`, `input1`, `input2`.

Tieto riadky teda môžete modifikovať, pridávať a odoberať podľa potreby a spúšťať tak ktorýkoľvek zo vstupných súborov v priečinku `bakalarka/input`. Na spustenie je potrebné otvoriť v prostredí terminál a napísať príkaz `make`.

Po spustení konkrétneho vstupu sa v priečinku `bakalarka/output` vytvorí výstupný `.obj` súbor,
tento súbor sa počas behu programu po každých 50 krokoch prepisuje, môžeme teda sledovať aj
progres algoritmu. Výstupný súbor si môžeme pozrieť napr. na [tejto](https://3dviewer.net/) stránke.

## Návod na spustenie merania kritérií kvality

Ak by ste si chceli spustiť aj meranie kritérii kvality, je potrebné si otvoriť súbor
`measure.cpp` a opäť nájsť úplne naspodku súboru funkciu `main()`.

V priečinku `bakalarka/measure` sa nachádzajú priečinky `inputs`, `outputs`
a `measure_data`. 
- V priečinku `inputs` sa nachádzajú vstupné súbory.
- V priečinku `outputs` sa nachádzajú príslušné výstupné .obj modely.
- V priečinku `measure_data` sa nachádzajú namerané dáta.

Príkazy `run_input` a `run_all` fungujú podobne ako pri tvorbe modelov. V komentároch v main() 
je popísaný presný postup.

Meranie sa spúšťa v termináli príkazom `make run_measure` a vygeneruje výstupné súbory v priečinku `measure/measure_data`. V tomto priečinku sa už nachádzajú nejaké namerané dáta, ak však meranie spustíte, novo-namerané súbory sa prepíšu.
