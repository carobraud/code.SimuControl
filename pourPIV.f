C --> NX     : NOMBRE DE POINTS SUIVANT LA COORDONNEE X
C --> NY     : NOMBRE DE POINTS SUIVANT LA COORDONNEE Y
C --> NK     : NOMBRE DE POINTS POUR LA MATRICE D INFLUENCE
C --> NPE    : NOMBRE DE PROCESSEURS POUR MPI
C --> NZP    : NOMBRE DE POINTS SUIVANT LA COORDONNE Z
C --> NZ_MPI : NOMBRE DE POINTS SUIVANT LA COORDONNE Z 
C --> NZ     : NOMBRE DE MODES
C              SUR CHAQUE PROCESSEUR POUR MPI
C --> NA     : NOMBRE DE POINTS POUR LE CALCUL DES TERMES
C              NON-LINEAIRES : NA=NZP       -> sans antialiasing
C                              NA=NZP+NZP/2 -> avec antialiasing
C
C Attention : Recompiler le fichier subDNS.f !!!
c
      parameter (npe=1)
      parameter (nx=1024,ny=97,nzp=4)
      parameter (nz_mpi=nzp/npe,nz=nz_mpi/2)
      parameter (na=nzp)
      parameter (nk=2*(nx-1)+2*(ny-1))
c
      parameter (nin=75,nout=277)
c      parameter (nin=75,nout=400)
      parameter (nwal=nout-nin-2)
