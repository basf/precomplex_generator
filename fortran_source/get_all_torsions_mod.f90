! This file is part of the precomplex_generator
!
! Copyright (C) BASF SE 2022
!
! The precomplex_generator is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! The precomplex_generator is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with the precomplex_generator. If not, see <https://www.gnu.org/licenses/>.


MODULE fortran_utils

   IMPLICIT NONE

   REAL, PARAMETER :: pi=4.0*atan(1.0)

CONTAINS

   SUBROUTINE remove_blank (line)

      IMPLICIT NONE

      CHARACTER (len=*), INTENT(OUT)  :: line
      INTEGER                         :: n

      !line=trim(adjustl(line))

      DO n=1,len(line)
         IF (line(:1) == ' ') THEN
            line=line(2:)
         ELSE
            EXIT
         END IF
      END DO

   END SUBROUTINE

   INTEGER FUNCTION chemical_number(atnam)

      IMPLICIT NONE

      CHARACTER(2), INTENT(IN) :: atnam

      INTEGER               :: i

      CHARACTER *2 element
      DIMENSION element(109)

      DATA element /'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si',&
         'P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu',&
         'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru',&
         'Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr',&
         'Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W',&
         'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',&
         'Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf',&
         'Db','Sg','Bh','Hs','Mt'/

      l1: DO i=1,109
         IF (element(i) == atnam) THEN
            chemical_number=i
            EXIT l1
         END IF
      END DO l1

   END FUNCTION chemical_number

END MODULE


Program get_torsions
   IMPLICIT NONE

END PROGRAM get_torsions


SUBROUTINE get_tors(datFile)

   ! ----------------------------------------------------------------------------
   !  Function that reads in all torsions from the given input file (datfile)   !
   ! ----------------------------------------------------------------------------

   USE fortran_utils, ONLY : chemical_number, remove_blank
   IMPLICIT NONE
   INTEGER                        :: i,j,io_error,istat
   INTEGER                        :: natom, ntors
   INTEGER, ALLOCATABLE           :: con(:,:), tors(:,:), atomic_number(:)
   CHARACTER(90)                  :: line, line2, line3, datFile
   CHARACTER(2), ALLOCATABLE      :: atnam(:)

!---------------------------------------------------------------------------------------------------

   OPEN (86,file=datFile,iostat=io_error,position='REWIND',status='old')

   natom=0
   DO
      READ(86,*,iostat=io_error) line
      IF (io_error /= 0) EXIT
      natom=natom+1
   END DO

   REWIND(86)

   ALLOCATE(atnam(natom),stat=istat);           IF(istat/=0) STOP 'allocation error: atnam'
   ALLOCATE(con(natom,natom+1),stat=istat);     IF(istat/=0) STOP 'allocation error: con'
   ALLOCATE(tors(natom,3),stat=istat);          IF(istat/=0) STOP 'allocation error: tors'
   ALLOCATE(atomic_number(natom),stat=istat);   IF(istat/=0) STOP 'allocation error: atomic_number'

   con=0

   DO i=1,natom
      READ(86,'(a)') line
      CALL remove_blank(line)
      line=trim(line(index(line,' '):)) ! delete i
      CALL remove_blank(line)
      line2=trim(line(:index(line,' ')))
      READ(line2,*) atnam(i)
      atomic_number(i)=chemical_number(atnam(i))
      line=trim(line(index(line,' '):)) ! delete atnam(i)
      READ(line,*) con(i,1)
      DO j=2,con(i,1)+1
         CALL remove_blank(line)
         line=trim(line(index(line,' '):)) ! delete con(i,j)
         READ(line,*) con(i,j)
      END DO
   END DO

   CLOSE(86)


!---------------------------------------------------------------------------------------------------


   ! Call "get_all_torsions" to obtain all required torsional angles and write them to tors.dat

   OPEN(99,file='tors.dat')
   CLOSE(99, status='delete')
   OPEN(98,file='tors.dat',iostat=io_error,position='REWIND',status='new')


   tors=0
   ntors=0

   call get_all_torsions(natom, atomic_number, con, ntors, tors)

   WRITE(98,'(a)') '$tors'
   DO i=1,ntors
      WRITE(line2,*) tors(i,1)
      CALL remove_blank(line2)
      WRITE(line3,*) tors(i,2)
      CALL remove_blank(line3)
      WRITE(line,*) tors(i,3)
      CALL remove_blank(line)
      WRITE(98,'(a)') trim(line2(:index(line2,' ')))//' '//trim(line3(:index(line3,' ')))//' '&
      & //trim(line(:index(line,' ')))
   END DO
   WRITE(98,'(a)') '$end'
   CLOSE(98)

END SUBROUTINE get_tors


SUBROUTINE get_all_torsions (natom, atomic_number, con, ntors, tors)

! ----------------------------------------------------------------------------
!  From all theretically possible torisons, select only the "required" ones  !
! ----------------------------------------------------------------------------


   IMPLICIT NONE

   INTEGER, INTENT(IN)      :: natom, atomic_number(natom), con(natom,natom+1)
   INTEGER, INTENT(OUT)     :: ntors, tors(natom,3)
   INTEGER                  :: i, j, k, l, m, o, p
   LOGICAL                  :: check, check2, check3, ring

   ntors=0
   DO i=1,natom                   ! i = Atom number of the 1st atom of the torsion (definition)
      DO k=2,con(i,1)+1
         j=con(i,k)               ! j = Atom number of the 2nd atom of the torsion (definition)
         IF (j > i) THEN          ! Avoid duplicate torsional angles

            check=.FALSE.


            ! Do not rotate, if atom i or atom j has only 1 neighbor atom
            IF ((con(i,1) > 1) .AND. (con(j,1) > 1)) THEN
               check=.TRUE.
            END IF


            ! Do not rotate, if the bond is within a ring.
            IF (check) THEN
               call check_ring(i,j,natom,con,ring)
               IF (ring) THEN
                  check=.FALSE.
               END IF
            END IF


            ! C#C-bonds do not need to get rotated
            IF (check) THEN
               IF ((atomic_number(i) == 6) .AND. (con(i,1) == 2)) THEN
                  IF (con(i,2) == j) THEN
                     m=con(i,3)
                  ELSE
                     m=con(i,2)
                  END IF
                  IF ((atomic_number(m) == 6) .AND. (con(m,1) == 2)) THEN
                     check=.FALSE.
                  END IF
               END IF
            END IF

            IF (check) THEN
               IF ((atomic_number(j) == 6) .AND. (con(j,1) == 2)) THEN
                  IF (con(j,2) == i) THEN
                     m=con(j,3)
                  ELSE
                     m=con(j,2)
                  END IF
                  IF ((atomic_number(m) == 6) .AND. (con(m,1) == 2)) THEN
                     check=.FALSE.
                  END IF
               END IF
            END IF


            ! R-CX does not need to get rotated (linear C only has 2 neighbors)
            IF ((check) .AND. (atomic_number(i) == 6)) THEN
               IF (con(i,1) == 2) THEN
                  IF (con(i,2) == j) THEN
                     m=con(i,3)
                  ELSE
                     m=con(i,2)
                  END IF
                  IF (con(m,1) == 1) THEN
                     check=.FALSE.
                  END IF
               END IF
            END IF

            IF ((check) .AND. (atomic_number(j) == 6)) THEN
               IF (con(j,1) == 2) THEN
                  IF (con(j,2) == i) THEN
                     m=con(j,3)
                  ELSE
                     m=con(j,2)
                  END IF
                  IF (con(m,1) == 1) THEN
                     check=.FALSE.
                  END IF
               END IF
            END IF


            ! C(i)=N(j)=N does not need to get rotated
            IF (check) THEN
               IF ((atomic_number(i) == 6) .AND. (con(i,1) == 3 )) THEN
                  IF ((atomic_number(j) == 7) .AND. (con(j,1) == 2 )) THEN
                     IF (con(j,2) == i) THEN
                        p=con(j,3)
                     ELSE
                        p=con(j,2)
                     END IF
                     IF ((atomic_number(p) == 7) .AND. (con(p,1) == 1 )) THEN
                        check=.FALSE.
                     END IF
                  END IF
               END IF

               IF ((atomic_number(j) == 6) .AND. (con(j,1) == 3 )) THEN
                  IF ((atomic_number(i) == 7) .AND. (con(i,1) == 2 )) THEN
                     IF (con(i,2) == j) THEN
                        p=con(i,3)
                     ELSE
                        p=con(i,2)
                     END IF
                     IF ((atomic_number(p) == 7) .AND. (con(p,1) == 1 )) THEN
                        check=.FALSE.
                     END IF
                  END IF
               END IF
            END IF


            ! Rotation by 120 degree (check2=.FALSE.) or 180 degree (check2=.TRUE.)?
            check2=.FALSE.
            IF (check) THEN

               ! N=N only needs to get rotated by 180 degree
               IF ((atomic_number(i) == 7) .AND. (con(i,1) == 2 )) THEN
                  IF ((atomic_number(j) == 7) .AND. (con(j,1) == 2 )) THEN
                     check2=.TRUE.
                  END IF
               END IF


               ! C=C=C only needs to get rotated by 180 degree
               ! and C=C=CH2 not at all
               IF ((atomic_number(i) == 6) .AND. (con(i,1) == 3 )) THEN
                  IF ((atomic_number(j) == 6) .AND. (con(j,1) == 2 )) THEN
                     check2=.TRUE.
                     IF (con(j,2) == i) THEN
                        p=con(j,3)
                     ELSE
                        p=con(j,2)
                     END IF
                     IF ((atomic_number(p) == 6) .AND. (con(p,1) == 3 )) THEN
                        check3=.TRUE.
                        IF (con(p,2) == j) THEN
                           o=con(p,3)
                        ELSE
                           o=con(p,2)
                        END IF
                        DO l=1,con(p,1)
                           IF (con(p,l+1) /= j) THEN
                              IF (con(con(p,l+1),1) /= 1) THEN
                                 check3=.FALSE.
                              ELSEIF (atomic_number(con(p,l+1)) /= atomic_number(o)) THEN
                                 check3=.FALSE.
                              END IF
                           END IF
                        END DO
                     END IF
                  END IF
               END IF

               IF ((atomic_number(j) == 6) .AND. (con(j,1) == 3 )) THEN
                  IF ((atomic_number(i) == 6) .AND. (con(i,1) == 2 )) THEN
                     check2=.TRUE.
                     IF (con(i,2) == j) THEN
                        p=con(i,3)
                     ELSE
                        p=con(i,2)
                     END IF
                     IF ((atomic_number(p) == 6) .AND. (con(p,1) == 3 )) THEN
                        check3=.TRUE.
                        IF (con(p,2) == i) THEN
                           o=con(p,3)
                        ELSE
                           o=con(p,2)
                        END IF
                        DO l=1,con(p,1)
                           IF (con(p,l+1) /= i) THEN
                              IF (con(con(p,l+1),1) /= 1) THEN
                                 check3=.FALSE.
                              ELSEIF (atomic_number(con(p,l+1)) /= atomic_number(o)) THEN
                                 check3=.FALSE.
                              END IF
                           END IF
                        END DO
                     END IF
                  END IF
               END IF


               IF ((con(i,1) == 3) .AND. (atomic_number(i) == 6) .AND. (con(j,1) == 3) .AND. (atomic_number(j) == 6)) THEN
                  check2=.TRUE. !rotate by 180 degree (C=C double bond)
               END IF


               ! For C=C=C only one of the C=C-bonds must get rotated
               IF (check2) THEN
                  IF ((atomic_number(i) == 6) .AND. (con(i,1) == 2)) THEN
                     DO l=1,ntors
                        IF (tors(l,3) == 180 ) THEN
                           IF ((i == tors(l,1)) .OR. (i == tors(l,2))) THEN
                              check=.FALSE.
                           END IF
                        END IF
                     END DO
                  END IF
                  IF (( atomic_number(j) == 6) .AND. (con(j,1) == 2)) THEN
                     DO l=1,ntors
                        IF (tors(l,3) == 180 ) THEN
                           IF ((j == tors(l,1)) .OR. (j == tors(l,2))) THEN
                              check=.FALSE.
                           END IF
                        END IF
                     END DO
                  END IF
               END IF

               ! For C(aromat)-C(sp2) rotate by 120 degree
               IF (check2) THEN
                  IF ((atomic_number(i) == 6) .AND. (con(i,1) == 3)) THEN
                     IF ((atomic_number(j) == 6) .AND. (con(j,1) == 3)) THEN
                        IF (con(i,1) == 3) THEN
                           ring=.FALSE.
                           IF (con(j,2) /= i) THEN
                              m=con(j,2)
                           ELSE
                              m=con(j,3)
                           END IF
                           call check_ring(j,m,natom,con,ring)
                           IF (ring) THEN
                              check2=.FALSE. ! We check if both centers are planar in the python code
                           END IF
                           IF (con(i,2) /= j) THEN
                              m=con(i,2)
                           ELSE
                              m=con(i,3)
                           END IF
                           call check_ring(i,m,natom,con,ring)
                           IF (ring) THEN
                              check2=.FALSE. ! We check if both centers are planar in the python code
                           END IF
                        END IF
                     END IF
                  END IF
               END IF
            END IF

            IF (check) THEN
               ntors=ntors+1
               tors(ntors,1)=i
               tors(ntors,2)=j
               IF (check2) THEN
                  tors(ntors,3)=180
               ELSE
                  tors(ntors,3)=120
               END IF
            END IF
         END IF
      END DO
   END DO

END SUBROUTINE

SUBROUTINE check_ring(atom_1,atom_2,natom,con,ring)

   ! ----------------------------------------------------------------------------
   !  Detects, IF atom 1 and atom 2 are connected via a ring.                   !
   ! ----------------------------------------------------------------------------

   IMPLICIT NONE

   INTEGER, INTENT(IN)           :: atom_1,atom_2,natom,con(natom,natom+1)
   LOGICAL, INTENT(OUT)          :: ring

   INTEGER                       :: i,j,n,n_neighbors,istat
   INTEGER, ALLOCATABLE          :: neighbor_list(:)
   LOGICAL                       :: check


   ALLOCATE(neighbor_list(natom),stat=istat);  IF(istat/=0) STOP 'allocation error: rot'
   neighbor_list=0

   n_neighbors=1
   neighbor_list(n_neighbors)=atom_1
   ring=.FALSE.

   DO i=1,con(atom_1,1)
      IF (con(atom_1,i+1) /= atom_2) THEN
         n_neighbors=n_neighbors+1
         neighbor_list(n_neighbors)=con(atom_1,i+1)
      END IF
   END DO

   n=1
   l1: DO WHILE (n < n_neighbors)
      n=n+1
      DO i=1,con(neighbor_list(n),1)
         check=.TRUE.
         l2:       DO j=1,n_neighbors
            IF (con(neighbor_list(n),i+1) == neighbor_list(j)) THEN
               check=.FALSE.    ! Atom con(neighbor_list(n),i+1) is already included in neighbor_list
               EXIT l2
            END IF
         END DO l2

         IF (check) THEN
            IF (con(neighbor_list(n),i+1) /= atom_2) THEN
               n_neighbors=n_neighbors+1
               neighbor_list(n_neighbors)=con(neighbor_list(n),i+1)
            ELSE
               ring=.TRUE.
               !WRITE(*,*) 'WARNING: You have asked for a torsion inside a ring:',atom_1, atom_2
               EXIT l1
            END IF
         END IF
      END DO
   END DO l1

   DEALLOCATE(neighbor_list)

END SUBROUTINE
