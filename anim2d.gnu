# Cal definir:
# k : valor inicial de l'índex a l'animació
# dk : increment en k
# kmax : valor final de l'índex a l'animació
# p (s) : temps de pausa entre frames
# nom : nom del fitxer amb wf
print "k=", k
splot nom u 1:2:3 index k w l
pause p
k=k+dk
if (k<kmax) reread
