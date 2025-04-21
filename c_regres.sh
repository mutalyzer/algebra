set -o errexit

valgrind --error-exitcode=1 ./a.out TATCCGTCCG GAAAAATCCG
valgrind --error-exitcode=1 ./a.out GCCTGCCTAT GCCTAAAAAG
valgrind --error-exitcode=1 ./a.out TGCCTTA TAGCAGCCC
valgrind --error-exitcode=1 ./a.out CCC CACCCACA
valgrind --error-exitcode=1 ./a.out TTGGGAA CCGAATA
valgrind --error-exitcode=1 ./a.out GACAT ACACTCCCAT
valgrind --error-exitcode=1 ./a.out GAG AGGATAAGT
valgrind --error-exitcode=1 ./a.out CACGCTCGT ACACT
valgrind --error-exitcode=1 ./a.out GTTTA TCTTTCTGC
valgrind --error-exitcode=1 ./a.out TCAACAAAGG CAA
valgrind --error-exitcode=1 ./a.out CTTATAAT CCTACCG
