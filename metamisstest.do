cd C:\Users\rmjwiww\Documents
clear
input study r1 f1 m1 r2 f2 m2
1 10 80 10 10 80 0
2 20 70 10 10 50 40
end
gen logimor=_n

metamiss r1 f1 m1 r2 f2 m2, logimor(1) nograph
metamiss r1 f1 m1 r2 f2 m2, logimor(2) nograph
metamiss r1 f1 m1 r2 f2 m2, logimor(logimor) nograph

metamiss2 r1 f1 m1 r2 f2 m2, impmean(1) metanopt(nograph)
metamiss2 r1 f1 m1 r2 f2 m2, impmean(2) metanopt(nograph)
metamiss2 r1 f1 m1 r2 f2 m2, impmean(logimor) metanopt(nograph)
