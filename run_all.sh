#!/bin/sh
# Reproduces every computation of the paper.  Total time about 20 minutes on two cores
# (the exhaustive search for six speeds up to 110 is the slow part).
set -e
make
echo "== Theorem 3.7: the three critical subtori of (R/Z)^6";      python3 classify6.py
echo "== Theorem 3.8: the two critical subtori of (R/Z)^5";        python3 lr5.py
echo "== Theorem 4.2: one-very-fast-runner tori, six speeds";      ./fastrunner 30000
echo "== Theorem 4.2: one-very-fast-runner tori, five speeds";     ./fastrunner5 20000
echo "== Proposition 4.1 (numerical check)";                       ./prejump
echo "== Section 5: the relative spectrum of U_C (Theorem 5.6)";   python3 uc_theorem.py
echo "== consistency check of (5.2) against exact ML up to 300";   ./ucdump 300 > ucdump300.txt && python3 jk_uc2.py 300
echo "== direct scan of U_C up to 300";                            ./ucscan 300
echo "== exhaustive searches (Section 6)";
./survey 200 4 > survey_n4_200.txt; tail -2 survey_n4_200.txt
./survey 130 5 > survey_n5_130.txt; tail -2 survey_n5_130.txt
./survey 110 6 > survey_n6_110.txt; tail -2 survey_n6_110.txt
./survey 50 7  > survey_n7_50.txt;  tail -2 survey_n7_50.txt
./survey 40 8  > survey_n8_40.txt;  tail -2 survey_n8_40.txt
echo "done"
