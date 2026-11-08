# Etap 20 — semantycznie równoważne strategie aktualizacji certyfikatów

## Co dokładnie porównano
- Eager: po każdej zmianie aktualizuje ścieżkę od węzła do korzenia.
- Lazy: oznacza przodków jako brudnych i oblicza ich stany dopiero przed kolejnym zapytaniem.
- Full: przelicza wszystkie węzły przy każdym odczycie certyfikatu.

W każdej strategii użyto **tych samych** drzew, dwustanowych relacji zgodności, wektorów wkładów, strumieni aktualizacji i harmonogramu odczytu certyfikatów. Sprawdzono wyniki po każdym odczycie względem pełnego DP. 18 uruchomień (2 rozmiary × 3 częstotliwości × 3 powtórzenia), 120 aktualizacji na przebieg, 24 zapytania.

## Ważne ograniczenia
- Drzewa i dozwolone relacje stanów są **syntetyczne**. Dane nie dowodzą biologicznych zależności anotacji.
- `query_every_updates` oznacza odczyt po każdej 1., 5. albo 20. aktualizacji.
- Czasy obejmują pełen przebieg 120 aktualizacji oraz wymaganą liczbę odczytów. Nie porównywać tych wartości bez uwzględnienia częstotliwości odczytów.
- To prototyp Python/NumPy, a nie porównanie z natywną implementacją DBSP/F-IVM.
- Porównanie pełnego DP może istotnie faworyzować metody cache; najważniejsze jest rozróżnienie eager vs lazy.
- Obserwowane różnice czasów zależą od obciążenia maszyny i wymagają odtworzenia na sprzęcie docelowym.

## Wyniki
Patrz `summary.csv` i `fig20_read_frequency.png`. Wariant lazy jest wolniejszy od eager przy odczycie po każdej zmianie, a szybszy gdy wiele zmian poprzedza odczyt.

## Reprodukcja
`python -m pip install numpy matplotlib`
`python run_stage20.py`

## Konsekwencja dla ICDE
Należy jawnie modelować nie tylko koszt pojedynczej aktualizacji, ale mieszane obciążenie (updates + certificate reads). Dokładne DP i znaczniki dirty są klasycznymi technikami: nie deklarować ich jako samodzielnej nowości. Do mocnego wkładu potrzebne są dane o realnych sekwencjach wersji referencji i równoważne baseline'y z aktualizacją widoków.
