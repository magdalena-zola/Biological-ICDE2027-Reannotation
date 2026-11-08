# Publikacja repozytorium biologicznego ICDE 2027

Projekt zawiera artykuł LaTeX i PDF, kod, małe wejścia QIIME 2 pierwszego badania, zamrożone wyniki dwóch eksperymentów i testy kontrolne. Drugi zewnętrzny zbiór QIIME 2 jest wskazany adresem wraz z sumą SHA-256 — nie jest kopiowany w repozytorium.

1. Zweryfikuj prawa do redystrybucji spakowanych plików wejściowych (szczególnie `work/inputs/*.qza`, `*.qzv` oraz pliku z połączonymi klasyfikacjami). Zatwierdź prawa autorów do kodu i wybór licencji. Zobacz `RIGHTS_AND_RELEASE.md`.
2. Na Windows zainstaluj Git oraz GitHub CLI `gh`, następnie wykonaj `gh auth login`.
3. Rozpakuj paczkę `Biological_ICDE2027_GITHUB_PUBLISH_Windows.zip`, a potem uruchom:

   ```powershell
   powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\PUBLIKUJ_BIOLOGICZNY_ICDE2027.ps1
   ```

4. Skrypt wymaga jawnego potwierdzenia prawa do publikacji i stworzy **nowe** publiczne repozytorium pod podaną nazwą — nie zmieni repozytorium SynOD ani TideMark.
5. Otwórz powstały adres bez zalogowania, sprawdź workflow GitHub Actions i dopiero wtedy podaj rzeczywisty URL w CMT.

**Uwaga:** repozytorium nie zostało publicznie opublikowane w chwili wygenerowania tej paczki; skrypt nie zastępuje wymaganego potwierdzenia licencji/danych. W przypadku EAB konieczne jest zapewnienie pełnej odtwarzalności wszystkich zadeklarowanych wyników bez wyjątków.
