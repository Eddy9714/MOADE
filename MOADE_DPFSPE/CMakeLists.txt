# CMakeList.txt: progetto CMake per MOADE_DPFSPE. Includere l'origine e definire
# qui la logica specifica del progetto.
#
cmake_minimum_required (VERSION 3.8)

# Aggiungere l'origine all'eseguibile di questo progetto.
add_executable (MOADE_DPFSPE "MOADE_DPFSPE.cpp" "MOADE_DPFSPE.h" "Globali.cpp" "Globali.h" "Istanza.cpp" "Istanza.h" "Modulo.cpp" "Modulo.h" "Permutazione.cpp" "Permutazione.h" "PermutazioneI.cpp" "PermutazioneI.h" "PermutazioneST.cpp" "PermutazioneST.h" "PermutazioneT.h" "PermutazioneT.cpp" "Main.cpp" "Main.h" "GruppoPDZN.cpp" "GruppoPDZN.h" "IndiciRandom.cpp" "IndiciRandom.h")

# TODO: Aggiungere i test e, se necessario, installare le destinazioni.
