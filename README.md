# PSMMA
These routines are used for testing Example 1 in our TPDS paper.
"A parallel structured divide-and-conquer algorithm for symmetric tridiagonal eigenvalue problems".

SRRSC is implemented in Cauchylowrank.f90.


Its usage is quite easy. The testing routine is psmma_testcauchy.f90. You can change
1) n = 8192, ...  to choose different sizes of matrices
2) nb = 64, ...   to choose different N_B
3) Redist=.true.  to test PSMMA_WRedist
         =.false. do not redistribute matrix A from BCDD to BDD


Right now these routine always compute low-rank approximation, and
you can turn it off by modifying the lines 246-247 of pscauchy_compute.f90
to
!         call is_intersect( LengthC,CIndex,LengthR,RIndex,islowrank )
         islowrank=.false.

It will turn off the low-rank approximation and then you can test PSMMA_Nlowrank.


The compilation is easy and you only need to modify SLmake.inc and Makefile.
Two examples of SLmake.inc are included, one is used on my Laptop and another is used
on Tianhe2.
~            
