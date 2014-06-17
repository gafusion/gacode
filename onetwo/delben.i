c
c --- INCLUDE file delben.i
c
c --- NOTE that common block /DELBEN/ is now properly data-aligned.
c               PLEASE maintain this when you make modifications. ... JF
c
      common /delben/ qben(kjm1),
     .                rben(kjm1), bq(kjm1), cq(kjm1), dq(kjm1), qs,
     .                nben
c
