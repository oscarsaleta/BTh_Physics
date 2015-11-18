# Cal definir:
# nom : nom del fitxer
# k : valor inicial de l'índex a l'animació
# kmax : valor final de l'índex a l'animació
# p (s) : temps de pausa entre frames

print "k=", k
splot nom using 1:2:3 index k w l
pause p
k=k+1
if (k<=kmax) reread
