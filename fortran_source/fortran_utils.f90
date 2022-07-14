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

REAL(8) FUNCTION norm(v_in)

  INTEGER                         :: i
  REAL(8), INTENT(IN)             :: v_in(3)
  REAL(8)                         :: temp(3)

  DO i = 1,3
      temp(i) = v_in(i)**2
  END DO

  norm = sqrt(sum(temp(:)))

END FUNCTION norm

SUBROUTINE normalize(v_inout)

  INTEGER                         :: i
  REAL(8), INTENT(INOUT)          :: v_inout(3)
  REAL(8)                         :: temp(3), len

  DO i = 1,3
      temp(i) = v_inout(i)**2
  END DO

  len = sqrt(sum(temp(:)))

  DO i = 1,3
      v_inout(i) = v_inout(i) / len
  END DO

END SUBROUTINE normalize

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

REAL(8) FUNCTION atomic_mass(atnam)

  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: atnam

  REAL(8) element
  DIMENSION element(109)

  DATA element /1.0079, 4.0026, &
                6.941, 9.0122, 10.811, 12.011, 14.007, 15.999, 18.988, 20.180, &
                22.990, 24.305, 26.982, 29.086, 30.974, 32.065, 35.453, 39.948, &
                39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, &
                  69,723, 72.64, 74.922, 78.96, 79.904, 83.798, &
                85.468, 87.62, 88.906, 91.224, 92.906, 95.96, 98.91, 101.07, 102.91, 106.42, 107.87, 112.41, &
                  114.82, 118.71, 121.76, 127.6, 126.9, 131.29, &
                132.91, 137.33, 138.91, 140.12, 140.91, 144.24, 146.9, 150.36, 151.96, 157.25, 158.93, 162.5, 164.93, 167.26, &
                  168.93, 173.05, 174.97, 178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59, &
                  204.38, 207.2, 208.98, 209.98, 210.0, 222.0, &
                223, 226.03, 227.0, 232.04, 231.04, 238.03, 237.05, 244.10, 243.10, 247.10, 251.10, 254.10, 257.10, 258.0, 259.0, &
                  260.0, 261.0, 262.0, 263.0, 262.0, 265.0, 266.0 /

  atomic_mass = element(atnam)

END FUNCTION atomic_mass

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

CHARACTER(2) FUNCTION chemical_symbol(index)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: index

  CHARACTER(2)        :: elements(109)
  DATA elements / &
          'H',  'He', 'Li', 'Be',         &
          'B',  'C',  'N',  'O',  'F',    &
          'Ne', 'Na', 'Mg', 'Al', 'Si',   &
          'P',  'S',  'Cl', 'Ar', 'K',    &
          'Ca', 'Sc', 'Ti', 'V',  'Cr',   &
          'Mn', 'Fe', 'Co', 'Ni', 'Cu',   &
          'Zn', 'Ga', 'Ge', 'As', 'Se',   &
          'Br', 'Kr', 'Rb', 'Sr', 'Y',    &
          'Zr', 'Nb', 'Mo', 'Tc', 'Ru',   &
          'Rh', 'Pd', 'Ag', 'Cd', 'In',   &
          'Sn', 'Sb', 'Te', 'I',  'Xe',   &
          'Cs', 'Ba', 'La', 'Ce', 'Pr',   &
          'Nd', 'Pm', 'Sm', 'Eu', 'Gd',   &
          'Tb', 'Dy', 'Ho', 'Er', 'Tm',   &
          'Yb', 'Lu', 'Hf', 'Ta', 'W',    &
          'Re', 'Os', 'Ir', 'Pt', 'Au',   &
          'Hg', 'Tl', 'Pb', 'Bi', 'Po',   &
          'At', 'Rn', 'Fr', 'Ra', 'Ac',   &
          'Th', 'Pa', 'U',  'Np', 'Pu',   &
          'Am', 'Cm', 'Bk', 'Cf', 'Es',   &
          'Fm', 'Md', 'No', 'Lr', 'Rf',   &
          'Db', 'Sg', 'Bh', 'Hs', 'Mt'    /

  chemical_symbol = elements(index)

END FUNCTION chemical_symbol

SUBROUTINE rotation(coord,phi,vec)

  IMPLICIT NONE
  REAL(8), intent(inout)     ::  coord(3)
  REAL(8), INTENT(IN)        ::  phi,vec(3)
  REAL(8), ALLOCATABLE       ::  old(:)

  ALLOCATE(old(3))

  old=coord

  coord(1)=old(1)*(cos(phi)+(vec(1)**2)*(1-cos(phi)))+old(2)*(vec(1)*vec(2)*(1-cos(phi))+vec(3)*sin(phi))+ &
          old(3)*(vec(1)*vec(3)*(1-cos(phi))-vec(2)*sin(phi))
  coord(2)=old(1)*(vec(1)*vec(2)*(1-cos(phi))-vec(3)*(sin(phi)))+old(2)*(cos(phi)+(vec(2)**2)*(1-cos(phi)))+ &
          old(3)*(vec(2)*vec(3)*(1-cos(phi))+vec(1)*sin(phi))
  coord(3)=old(1)*(vec(1)*vec(3)*(1-cos(phi))+vec(2)*sin(phi))+old(2)*(vec(2)*vec(3)*(1-cos(phi))-vec(1)*sin(phi))+ &
          old(3)*(cos(phi)+(vec(3)**2)*(1-cos(phi)))

  DEALLOCATE(old)

END SUBROUTINE rotation

subroutine write_movie(natom, atomic_number, coord, i_conf)

  implicit none
  integer, parameter                              :: dp = kind(1.0d0)
  integer, intent(in)                             :: natom            
  integer, intent(in)                             :: atomic_number(natom)
  integer                                         :: j, io_error
  real(kind=dp), intent(in), dimension(natom, 3)  :: coord             
  integer, intent(in)                             :: i_conf
  logical                                         :: exist
  
  inquire(file="movie.xyz", exist=exist)
  if (exist) then
    open(94, file="movie.xyz", status="old", position="append", action="write")
  else
    open(94, file="movie.xyz", status="new", action="write")
  end if
  
  write(94,*) natom
  write(94,*) real(i_conf)
  do j=1,natom
     write(94,'(a3,3f14.7)') chemical_symbol(atomic_number(j)), coord(j,:)
  end do
  close(94)

end subroutine write_movie

SUBROUTINE write_molecule(natom, atomic_number, coord, which_conf)

  IMPLICIT NONE
  INTEGER, parameter                              :: dp = kind(1.0)
  INTEGER, INTENT(IN)                             :: natom
  INTEGER, INTENT(IN)                             :: atomic_number(natom)
  INTEGER                                         :: j
  CHARACTER(15)                                   :: filename
  REAL(8), INTENT(IN), DIMENSION(natom, 3)        :: coord
  INTEGER, INTENT(IN)                             :: which_conf
  
  WRITE(filename,'(a4,i0.7,a4)') 'conf', which_conf,'.xyz'
  OPEN(92,file=filename)

  WRITE(92,*) natom
  WRITE(92,*) ' '
  DO j=1,natom
      WRITE(92,'(a3,3f14.7)') chemical_symbol(atomic_number(j)), coord(j,:)
  END DO
  CLOSE(92)

END SUBROUTINE write_molecule

REAL(8) FUNCTION cylin_func(natom, coord, adds)

    IMPLICIT NONE

    INTEGER, INTENT(IN)             :: adds(2), natom
    REAL(8), INTENT(IN)             :: coord(natom,3)
    REAL(8)                         :: v1(3), vr1(3), vr2(3), vat(3)
    REAL(8)                         :: dist_cyl, cp(3), ch1, ch2, min_dist_cyl
    INTEGER                         :: i
    LOGICAL                         :: found

    vr1(:) = coord(adds(1),:)
    vr2(:) = coord(adds(2),:)
    v1(:)  = vr2(:) - vr1(:)
    ch1 = 0; ch2 = 0; min_dist_cyl = 999; found = .False.

    ! Calculate radius of smallest cylinder around add-vector, defined by add-atoms and closest disturbing atom

    DO i = 1,natom
        IF (i /= adds(1) .and. i /= adds(2)) THEN
            vat(:) = coord(i,:)
            ch1 = dotp((vat-vr1),(vr2-vr1))
            ch2 = dotp((vat-vr2),(vr1-vr2))
            IF (ch1 > 0 .and. ch2 < 0) THEN
                call crossp((vat-vr1), (vr2-vr1), cp)
                dist_cyl = norm(cp) / norm(vr2-vr1)
                IF (dist_cyl <= min_dist_cyl .and. dist_cyl > 1E-3) THEN
                    min_dist_cyl = dist_cyl
                    found = .True.
                END IF
            END IF
        END IF
    END DO

    IF (found) THEN
        cylin_func = min_dist_cyl
    ELSE
        cylin_func = 1000
    END IF

END FUNCTION cylin_func

SUBROUTINE cop(natom, coord, atomic_number, con, add_is, cop_out, neighbours, mass)

    INTEGER, INTENT(IN)                             :: add_is, con(natom, natom+1), natom, atomic_number(natom)
    REAL(8), INTENT(IN)                             :: coord(natom, 3)
    REAL(8)                                         :: nb(6,4), avg(3), full_mass, m_avg(3)
    REAL(8), INTENT(INOUT)                          :: cop_out(3)
    INTEGER                                         :: nn, i, j
    LOGICAL                                         :: neighbours, mass
    
    nn=0; nb=0; avg = 0; cop_out = 0; full_mass = 0; m_avg = 0
    
    IF (neighbours) THEN

      DO j = 1,6
        IF ((con(add_is,j+1) > 0)) THEN
          nb(j,1) = con(add_is,j+1)
          IF (.NOT. all(abs(coord(int(nb(j,1)),1:3)) < 1E-3)) THEN
            
            nb(j,2:4) = coord(int(nb(j,1)),1:3)
            avg(:) = avg(:) + nb(j,2:4)
            m_avg(:) = m_avg(:) + nb(j,2:4) * atomic_mass(atomic_number(int(nb(j,1))))
            full_mass = full_mass + atomic_mass(atomic_number(int(nb(j,1))))
            nn = nn+1

          END IF
        END IF
      END DO
        
      IF (mass) THEN
        cop_out(:) = m_avg(:)/full_mass
      ELSE IF (.NOT. mass) THEN
        cop_out(:) = avg(:)/nn
      END IF

    ELSE IF (.NOT. neighbours) THEN

      DO i = 1, natom
          avg(:) = avg(:) + coord(i,:)
          m_avg = m_avg(:) + coord(i,:) * atomic_mass(atomic_number(i))
          full_mass = full_mass + atomic_mass(atomic_number(i))
      END DO

      IF (mass) THEN
        cop_out(:) = m_avg(:)/full_mass
      ELSE IF (.NOT. mass) THEN
        cop_out(:) = avg(:)/natom
      END IF

  END IF

END SUBROUTINE cop

SUBROUTINE dihed(p0, p1, p2, p3, dih)

    IMPLICIT NONE

    REAL(8), DIMENSION(3)           :: p0, p1, p2, p3, cr, b0, b1, b2, v, w
    REAL(8)                         :: x, y, dih, temp

    dih = 0.1

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    call normalize(b1)

    v = b0 - dotp(b0, b1) * b1
    w = b2 - dotp(b2, b1) * b1
    x = dotp(v, w)

    call crossp(b1, v, cr)

    y = dotp(cr, w)
    temp = abs(atan2(y,x))
    dih = temp * 180.0/pi

END SUBROUTINE dihed

REAL(8) FUNCTION dotp(v1, v2)

    INTEGER                         :: i
    REAL(8)                         :: temp(3)
    REAL(8), INTENT(IN)             :: v1(3), v2(3)

    DO i = 1,3
        temp(i) = v1(i)*v2(i)
    END DO

    dotp = sum(temp(:))

END FUNCTION dotp

SUBROUTINE crossp(vin1, vin2, vout)

    IMPLICIT NONE

    REAL(8), INTENT(IN)                  :: vin1(3), vin2(3)
    REAL(8), INTENT(INOUT)               :: vout(3)

    vout(1) = vin1(2)*vin2(3) - vin1(3)*vin2(2)
    vout(2) = vin1(3)*vin2(1) - vin1(1)*vin2(3)
    vout(3) = vin1(1)*vin2(2) - vin1(2)*vin2(1)

END SUBROUTINE crossp

SUBROUTINE anglefunc(v1, v2, x)

    IMPLICIT NONE
    REAL(8), INTENT(IN)                 :: v1(3), v2(3)
    REAL(8), INTENT(OUT)                :: x
    REAL(8)                             :: dotprod, before, denom

    x=0

    dotprod = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
    denom = (norm(v1) * norm(v2))
    
    before = dotprod/denom
    
    IF (before > 1.0) THEN
        before = 1.0 - 1E-6
    ELSE IF (before < -1.0) THEN
        before = -1.0 + 1E-6
    END IF

    x = (abs(acos(before)) * 180)/pi

END SUBROUTINE anglefunc

SUBROUTINE kov_rad(at_1,at_2,rad_out)

  IMPLICIT NONE
  INTEGER                             ::  i
  CHARACTER(2), INTENT(IN)            ::  at_1, at_2
  CHARACTER(2)                        ::  nam
  REAL(8), INTENT(OUT)                ::  rad_out
  REAL(8)                             ::  rad, rad_1, rad_2

  DIMENSION nam(49), rad(49)

  DATA nam /'H','Li','Be','B','C','N','O','F','Na','Mg','Al','Si','P','S','Cl',&
          & 'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',&
          & 'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I'/

  ! Radii of elements found at http://de.wikipedia.org/wiki/Kovalenter_Radius
  ! Ionic radii for alkaline and alkaline earth metals + Al found in ref: https://doi.org/10.1107/S0567739476001551

  DATA rad /0.32,0.76,0.45,0.82,0.77,0.71,0.73,0.71,1.02,0.72,0.803,1.11,1.06,1.02,0.99,&
          & 1.38,1.00,1.44,1.36,1.25,1.27,1.39,1.25,1.26,1.21,1.38,1.31,1.26,1.22,1.21,1.16,1.14,&
          & 1.66,1.32,1.62,1.48,1.37,1.45,1.31,1.26,1.35,1.31,1.53,1.48,1.44,1.41,1.38,1.35,1.33/

l11: DO i=1,49
     IF (at_1 == nam(i)) THEN
        rad_1=rad(i)
        EXIT l11
     END IF
   END DO l11
l12: DO i=1,49
     IF (at_2 == nam(i)) THEN
        rad_2=rad(i)
        EXIT l12
     END IF
   END DO l12
   rad_out=1.3*(rad_1+rad_2)

END SUBROUTINE kov_rad

SUBROUTINE kov_rad_num(at_num_1,at_num_2,rad_out)

  IMPLICIT NONE
  INTEGER                             ::  i
  INTEGER, INTENT(IN)                 ::  at_num_1,at_num_2
  CHARACTER(2)                        ::  at_1, at_2
  CHARACTER(2)                        ::  nam
  REAL(8), INTENT(OUT)                ::  rad_out
  REAL(8)                             ::  rad, rad_1, rad_2

  DIMENSION nam(49), rad(49)

  DATA nam /'H','Li','Be','B','C','N','O','F','Na','Mg','Al','Si','P','S','Cl',&
          & 'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',&
          & 'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I'/

  ! Radii of elements found at http://de.wikipedia.org/wiki/Kovalenter_Radius
  ! Ionic radii for alkaline and alkaline earth metals + Al found in ref: https://doi.org/10.1107/S0567739476001551

  DATA rad /0.32,0.76,0.45,0.82,0.77,0.71,0.73,0.71,1.02,0.72,0.535,1.11,1.06,1.02,0.99,&
          & 1.38,1.00,1.44,1.36,1.25,1.27,1.39,1.25,1.26,1.21,1.38,1.31,1.26,1.22,1.21,1.16,1.14,&
          & 1.66,1.32,1.62,1.48,1.37,1.45,1.31,1.26,1.35,1.31,1.53,1.48,1.44,1.41,1.38,1.35,1.33/

  at_1=chemical_symbol(at_num_1)
  at_2=chemical_symbol(at_num_2)

  l11: DO i=1,49
      IF (at_1 == nam(i)) THEN
          rad_1=rad(i)
          EXIT l11
      END IF
  END DO l11
  l12: DO i=1,49
      IF (at_2 == nam(i)) THEN
          rad_2=rad(i)
          EXIT l12
      END IF
  END DO l12
  rad_out=1.3*(rad_1+rad_2)

END SUBROUTINE kov_rad_num

SUBROUTINE reac_vec(natom, coord, atomic_number, con, add_is, rvec, break_bond, planar, stereo)

    ! ----------------------------------------------------------------------
    !  Function that determines the reaction vector(s) for a reactive atom !
    ! ----------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER, INTENT(IN)                             :: natom,add_is,con(natom,natom+1),break_bond(2), atomic_number(natom)
    LOGICAL, INTENT(OUT)                            :: planar
    REAL(8), INTENT(IN)                             :: coord(natom, 3)
    INTEGER                                         :: i, j, nn, id, br_atom1, br_atom2, istat
    REAL(8)                                         :: nb(6,4), ang, dih, br_angle, len_break
    REAL(8)                                         :: cop_n(3), cop_all(3), copm_n(3), copm_all(3), br_vec(3), rvec(3)
    REAL(8)                                         :: tet_vec(3), rvec_temp(3), copm(3), v_plane(3), p(3), d
    REAL(8), DIMENSION(3)                           :: v1, v2, v3, v4, va, vb, vc
    LOGICAL                                         :: found
    INTEGER                                         :: stereo

    nn = 0; nb = 1000; dih = 0.1; cop_n = 0; cop_all = 0; va=0; vb=0; vc=0; rvec=0; br_angle = 0.1; id = 0
    planar = .FALSE.; ang=0.1; br_atom1 = 0; br_atom2 = 0; len_break = 0; found=.FALSE.

    id = 0
    DO j = 1,6
        IF (con(add_is,j+1) > 0) THEN
            nn=nn+1
            nb(j,1) = con(add_is,j+1)
            nb(j,2:4) = coord(nint(nb(j,1)),1:3)
            ! Shape of nb is number of rows = number of neighbors, 
            ! In each row, the columns are = ID of the neighbor atom(s) and its coordinates
        END IF
    END DO

    call cop(natom, coord, atomic_number, con, add_is, cop_all, .FALSE., .FALSE.)
    call cop(natom, coord, atomic_number, con, add_is, cop_n, .TRUE., .FALSE.)
    call cop(natom, coord, atomic_number, con, add_is, copm_all, .FALSE., .TRUE.)
    call cop(natom, coord, atomic_number, con, add_is, copm_n, .TRUE., .TRUE.)

    copm(:) = 0.5*cop_all(:) + 0.5*copm_all(:)

    v1 = coord(add_is,:)-nb(1,2:4)
    v2 = coord(add_is,:)-nb(2,2:4)
    v3 = coord(add_is,:)-nb(3,2:4)
    v4 = coord(add_is,:)-nb(4,2:4)

    IF (abs(nn - 4) <= 1E-3) THEN

        call anglefunc(v1, v2 ,ang)
        call dihed(nb(1,2:4), nb(2,2:4), nb(3,2:4), nb(4,2:4), dih)

        IF (((50 <= ang) .AND. (ang <= 130)) .AND. ((55 <= dih) .AND. (dih <= 85))) THEN

            IF (((break_bond(1) > 0) .AND. (abs(break_bond(1)) < 1000)) .AND. &
                ((break_bond(2) > 0) .AND. (abs(break_bond(2)) < 1000))) THEN

                br_atom1 = break_bond(1); br_atom2 = break_bond(2)
                br_vec = coord(br_atom1,:) - coord(br_atom2,:)
                len_break = br_vec(1)**2 + br_vec(2)**2 + br_vec(3)**2

                IF (len_break > 0.1) THEN
                    DO i=1,4
                        rvec_temp(:) = coord(add_is,:)-nb(i,2:4)
                        CALL anglefunc(rvec_temp, br_vec, br_angle)
                        IF (abs(br_angle) < 5.0 .OR. abs(abs(br_angle) - 180) < 5) THEN
                          tet_vec(:)=rvec_temp(:)
                          found = .TRUE.
                        END IF
                    END DO

                    IF (found) THEN
                        rvec(:) = tet_vec(:)
                    ELSE
                        rvec_temp(:) = coord(add_is,:)-nb(i,2:4)
                    END IF
                END IF
            ELSE
                rvec = nb(1,2:4)-coord(add_is,:)
            END IF
        ELSE
            rvec = nb(1,2:4)-coord(add_is,:)
        END IF


    ELSE IF (abs(nn - 3) <= 1E-3) THEN

        call anglefunc(v1, v2 ,ang)
        call dihed(nb(1,2:4),coord(add_is,:),nb(2,2:4),nb(3,2:4), dih)

        IF (abs(dih - 180d0) < 1E-6) THEN
            dih = 180d0 - 1E-6
        END IF
        
        IF (((100 <= ang) .AND. (ang <= 140)) .AND. (((160 <= dih) .AND. (dih <= 190)) .OR. (abs(dih) < 10))) THEN

        IF (stereo == 0) THEN
            call crossp(v1, v2, rvec)
            planar = .True.
        ELSE IF (stereo == 1) THEN
            call crossp(v2, v1, rvec)
            planar = .True.
        END IF


        ELSE IF (((100 <= ang) .AND. (ang <= 120)) .AND. ((110 <= dih) .AND. (dih < 160))) THEN

            rvec = coord(add_is,:) - cop_n

        END IF


    ELSE IF (abs(nn - 2) <= 1E-3) THEN


      call anglefunc(v1, v2 ,ang)
      IF ((170 <= ang) .AND. (ang <= 190)) THEN

          v_plane = nb(1,2:4) - nb(2,2:4)
          copm(:)=copm(:)-coord(add_is,:)
          d = dotp(copm(:), v_plane(:)) / norm(v_plane(:))
          call normalize(v_plane(:))
          p(:) = d*v_plane(:)
          rvec(:) = -copm(:) - p(:)

      ELSE
          rvec = coord(add_is,:) - cop_n
      END IF


    ELSE IF (abs(nn - 1) <= 1E-3) THEN
        rvec = coord(add_is,:) - cop_n

    ELSE IF (abs(nn) < 1E-3) THEN
        rvec = (/ 0, 0, 1 /)
    END IF


    call normalize(rvec)

    IF (.not. rvec(1) == rvec(1)) THEN

        id = 0; nn=0
        DO j = 1,6
            IF (con(add_is,j+1) > 0) THEN
                nn=nn+1
                nb(j,1) = con(add_is,j+1)
                nb(j,2:4) = coord(nint(nb(j,1)),1:3)
                ! Shape of nb is number of rows = Anzahl Nachbarn,
                ! in jeder Reihe die Spalten = ID des Nachbarn und dessen Koordinaten
            END IF
        END DO
        nn=nn-1

        call cop(natom, coord, atomic_number, con, add_is, cop_all, .FALSE., .FALSE.)
        call cop(natom, coord, atomic_number, con, add_is, cop_n, .TRUE., .FALSE.)
        call cop(natom, coord, atomic_number, con, add_is, copm_all, .FALSE., .TRUE.)
        call cop(natom, coord, atomic_number, con, add_is, copm_n, .TRUE., .TRUE.)

        copm(:) = 0.5*cop_all(:) + 0.5*copm_all(:)

        v1 = coord(add_is,:)-nb(1,2:4)
        v2 = coord(add_is,:)-nb(2,2:4)
        v3 = coord(add_is,:)-nb(3,2:4)
        v4 = coord(add_is,:)-nb(4,2:4)

        IF (abs(nn - 4) <= 1E-3) THEN

            call anglefunc(v1, v2 ,ang)
            call dihed(nb(1,2:4), nb(2,2:4), nb(3,2:4), nb(4,2:4), dih)

            IF (((50 <= ang) .AND. (ang <= 130)) .AND. ((55 <= dih) .AND. (dih <= 85))) THEN

                IF (((break_bond(1) > 0) .AND. (abs(break_bond(1)) < 1000)) .AND. &
                        ((break_bond(2) > 0) .AND. (abs(break_bond(2)) < 1000))) THEN

                    br_atom1 = break_bond(1); br_atom2 = break_bond(2)
                    br_vec = coord(br_atom1,:) - coord(br_atom2,:)
                    len_break = br_vec(1)**2 + br_vec(2)**2 + br_vec(3)**2

                    IF (len_break > 0.1) THEN
                        DO i=1,4
                            rvec_temp(:) = coord(add_is,:)-nb(i,2:4)
                            CALL anglefunc(rvec_temp, br_vec, br_angle)
                            IF (abs(br_angle) < 5.0 .OR. abs(abs(br_angle) - 180) < 5) THEN
                                tet_vec(:)=rvec_temp(:)
                                found = .TRUE.
                            END IF
                        END DO

                        IF (found) THEN
                            rvec(:) = tet_vec(:)
                        ELSE
                            rvec_temp(:) = coord(add_is,:)-nb(i,2:4)
                        END IF
                    END IF
                ELSE
                    rvec = nb(1,2:4)-coord(add_is,:)
                END IF
            ELSE
                rvec = nb(1,2:4)-coord(add_is,:)
            END IF


        ELSE IF (abs(nn - 3) <= 1E-3) THEN

            call anglefunc(v1, v2 ,ang)
            call dihed(nb(1,2:4),coord(add_is,:),nb(2,2:4),nb(3,2:4), dih)

            IF (abs(dih - 180d0) < 1E-6) THEN
                dih = 180d0 - 1E-6
            END IF

            IF (((100 <= ang) .AND. (ang <= 140)) .AND. (((160 <= dih) .AND. (dih <= 190)) .OR. (abs(dih) < 10))) THEN

                call crossp(v1, v2, rvec)
                planar = .True.

            ELSE IF (((100 <= ang) .AND. (ang <= 120)) .AND. ((110 <= dih) .AND. (dih < 160))) THEN

                rvec = coord(add_is,:) - cop_n

            END IF


        ELSE IF (abs(nn - 2) <= 1E-3) THEN


            call anglefunc(v1, v2 ,ang)
            IF ((170 <= ang) .AND. (ang <= 190)) THEN

                v_plane = nb(1,2:4) - nb(2,2:4)
                copm(:)=copm(:)-coord(add_is,:)
                d = dotp(copm(:), v_plane(:)) / norm(v_plane(:))
                call normalize(v_plane(:))
                p(:) = d*v_plane(:)
                rvec(:) = -copm(:) - p(:)

            ELSE
                rvec = coord(add_is,:) - cop_n
            END IF

        ELSE IF (abs(nn - 1) <= 1E-3) THEN
            rvec = coord(add_is,:) - cop_n

        ELSE IF (abs(nn) < 1E-3) THEN
            rvec = (/ 0, 0, 1 /)
        END IF

        call normalize(rvec)

    END IF

END SUBROUTINE

SUBROUTINE ranking_func(this_case, ranking_values, counter, ranking_values_add, coord, fortran_weights, con, &
                        ranking_values_avg, ranking_values_min, ranking_values_max, max_ranking, final_coord, replaced,&
                        break_bond, add_i, atomic_number, sadd_thresh, rot_bonds_neigh, rot_bonds_path, intramol_small_ring)

    ! ---------------------------------------------------------------------------------
    !  Calculation and evaluation of a "ranking function" based on the ranking_values 
    ! (properties from which we assume that they are important for a suited precomplex!
    ! ---------------------------------------------------------------------------------

    IMPLICIT NONE

    REAL(8), INTENT(IN)             :: ranking_values(:), fortran_weights(:), sadd_thresh
    INTEGER, INTENT(IN)             :: counter, con(:,:), break_bond(:), add_i(:), atomic_number(:)
    INTEGER, INTENT(IN)             :: rot_bonds_neigh, rot_bonds_path
    CHARACTER(90), INTENT(IN)       :: this_case
    LOGICAL, INTENT(INOUT)          :: replaced

    REAL(8), INTENT(INOUT)          :: ranking_values_add(:), ranking_values_avg(:), ranking_values_min(:)
    REAL(8), INTENT(INOUT)          :: ranking_values_max(:)
    REAL(8), INTENT(INOUT)          :: max_ranking, coord(:,:), final_coord(:,:)
    REAL(8)                         :: final_values(7), denom, overall_ranking, fin_val_avg, rvec_angle, temp_val
    REAL(8)                         :: rvec1(3), rvec2(3), sign, abs_fin, sum_kov_rad
    INTEGER                         :: i, j, k, num_rank, temp(2), natom
    LOGICAL                         :: exist, planar, discard
    LOGICAL                         :: intramol_small_ring

    replaced = .FALSE.; planar = .FALSE.; rvec1 = 0; rvec2 = 0; discard = .FALSE.
    temp = shape(coord)
    natom = temp(1)


    ! VALUE NUMBER 1 : dmin

    final_values(1) = ranking_values(1)-0.832

    IF (final_values(1) > 1) THEN
        final_values(1) = sqrt(final_values(1))
    END IF
    final_values(1) = 0.5*final_values(1)


    ! VALUE NUMBER 2

    ranking_values_add(2) = ranking_values_add(2) + ranking_values(2)
    ranking_values_avg(2) = ranking_values_add(2) / counter
    final_values(2) = 3*(1/(1+exp(-(ranking_values(2)/ranking_values_avg(2) - 1)*20)))
    final_values(2) = (1d0 + 0.5 * (log(real(rot_bonds_neigh+rot_bonds_path)**2))) * final_values(2)

    ! VALUE NUMBER 7

    final_values(7) = 0

    IF (this_case == 'oneadd') THEN

        final_values(3:7) = 0

        overall_ranking = fortran_weights(1) * final_values(1) + &
                          fortran_weights(2) * final_values(2)


    ELSE IF (this_case == 'twoadds' .OR. this_case =='intramol' .OR. this_case =='shuttle') THEN

        ! VALUE NUMBER 3 : sadd

        call kov_rad_num(add_i(1), add_i(2), sum_kov_rad)

        IF (ranking_values(3) > sadd_thresh) THEN

            ! Highest prio threshold = user input
            discard = .TRUE.
            ! write(*,*) "discarded due to sadd thresh"

        ELSE
            IF (intramol_small_ring) THEN
                IF (ranking_values(3) < 0.7*sum_kov_rad) THEN
                    discard = .TRUE.
                    ! write(*,*) "discarded due to sadd thresh < 1.1*sum_kov_rad"
                ELSE IF (ranking_values(3) < 1.5*sum_kov_rad) THEN       
                    final_values(3) = 0
                ELSE
                    final_values(3) = 1+(5/(1+exp(2*ranking_values(3)-6)))
                END IF
            ELSE
                IF (ranking_values(3) < 1.1*sum_kov_rad) THEN
                    discard = .TRUE.
                    ! write(*,*) "discarded due to sadd thresh < 1.1*sum_kov_rad"
                ELSE IF (ranking_values(3) < 1.5*sum_kov_rad) THEN      
                    final_values(3) = 0
                ELSE
                    final_values(3) = 1+(5/(1+exp(2*ranking_values(3)-6)))
                END IF
            END IF

            IF (rot_bonds_path < 3) THEN
                final_values(3) = 3*final_values(3)
            ELSE IF (rot_bonds_path == 3) THEN
                final_values(3) = 2*final_values(3)
            ELSE IF ((rot_bonds_path > 3) .AND. (rot_bonds_path <= 6)) THEN
                final_values(3) = final_values(3)
            ELSE IF (rot_bonds_path > 6) THEN
                final_values(3) = (1 + 0.5*(log(real(rot_bonds_path)))) * final_values(3)
            END IF

        END IF

!       VALUE NUMBER 4 : rvec alignment

        final_values(4) = 5 - 5*ranking_values(4)
        

!       VALUE NUMBER 5: Dummy atoms!

        IF (ranking_values(5) > 0.0) THEN  
            final_values(5) = 0
            discard = .True.
            ! write(*,*) "discarded due to dummy atoms", ranking_values(5)
        ELSE
            final_values(5) = 20 * (1/(1+exp(6*ranking_values(5)))) - 10
        END IF


        ! VALUE NUMBER 6

        IF (ranking_values(6) > 100) THEN
            final_values(6) = 0
        ELSE
            final_values(6) = 8/(1+exp(-2*ranking_values(6))) - 4
        END IF


        ! Overall Ranking
        IF (discard) THEN
            overall_ranking = 0
        ELSE
            overall_ranking = fortran_weights(1) * final_values(1) + &
                              fortran_weights(2) * final_values(2) + &
                              fortran_weights(3) * final_values(3) + &
                              fortran_weights(4) * final_values(4) + &
                              fortran_weights(5) * final_values(5) + &
                              fortran_weights(6) * final_values(6)
        END IF
    END IF

    IF (.NOT. discard) THEN
        !write(*,*) overall_ranking, max_ranking
        IF (overall_ranking > max_ranking) THEN
            max_ranking = overall_ranking
            replaced = .TRUE.
            final_coord(:,:) = coord(:,:)
            ! new
            inquire(file="final_values.txt", exist=exist)
            if (exist) then
                open(73, file="final_values.txt", status="old", position="append", action="write")
            else
                open(73, file="final_values.txt", status="new", action="write")
                WRITE(73,'(6A20)',advance="no") "max_min","avg_dist","add_dist","angle_sub","dummy","cylin"
                WRITE(73,'(A20)') "overall"
            end if

            WRITE(73,'(7F20.6)') final_values(1), final_values(2), final_values(3), final_values(4), &
                                final_values(5), final_values(6), overall_ranking


            close(73)
        END IF
    END IF

END SUBROUTINE ranking_func


SUBROUTINE ranking_value_sub(coord, con, atomic_numbers, fortran_weights, this_case, counter, all_bonds, sadd, break_bond, &
                             ranking_values, rot_pun, stereo1, stereo2)

    ! ---------------------------------------------------------------------------------------------------------------------
    !  Calculation of the ranking_values = properties from which we assume that they are important for a suited precomplex!
    ! ---------------------------------------------------------------------------------------------------------------------

    IMPLICIT NONE

    REAL(8), INTENT(IN)             :: coord(:,:), fortran_weights(:), all_bonds
    INTEGER, INTENT(IN)             :: atomic_numbers(:), sadd(:), con(:,:), break_bond(:), rot_pun
    CHARACTER(90), INTENT(IN)       :: this_case
    INTEGER, INTENT(INOUT)          :: counter
    REAL(8), INTENT(OUT)            :: ranking_values(7)

    INTEGER                         :: natoms, i, j, k, at_num_sum, istat, factor
    REAL(8)                         :: at_dist, avg_dist, max_min_weight, avg_dist_weight, max_min, add_dist, rvec1(3)
    REAL(8)                         :: rvec2(3), numrot, angle_sub
    REAL(8)                         :: temp_val, angle_diff, all_dist, kov_dist, kov_dist_sum, dummy1(3), dummy2(3)
    REAL(8)                         :: dummy, cylin
    REAL(8)                         :: max_min_temp, radius
    REAL(8), ALLOCATABLE            :: dists(:)
    LOGICAL                         :: planar1, planar2
    INTEGER                         :: stereo1, stereo2

    rvec1=0; rvec2=0;
    max_min = 0; avg_dist = 0; add_dist = 0; angle_diff = 0; dummy = 0; cylin = 0; numrot = 0
    max_min_temp = 1000

    at_dist = 0; all_dist = 0; at_num_sum = 0; add_dist = 0; k = 1; natoms = size(atomic_numbers); planar1 = .FALSE.
    kov_dist_sum = 0; dummy = 0; cylin = 0; max_min = 0; planar2 = .FALSE.

    DO i = 1, natoms
      DO j = 1, natoms
        IF (j > i) THEN

          CALL kov_rad_num(atomic_numbers(i), atomic_numbers(j), kov_dist)
          kov_dist_sum = kov_dist_sum + kov_dist
          at_dist=at_dist + sqrt((coord(i,1)-coord(j,1))**2 + (coord(i,2)-coord(j,2))**2 + (coord(i,3)-coord(j,3))**2)

          max_min=sqrt((coord(i,1)-coord(j,1))**2 + (coord(i,2)-coord(j,2))**2 + (coord(i,3)-coord(j,3))**2)

          IF ((max_min > kov_dist) .AND. (max_min < max_min_temp)) THEN
            max_min_temp = max_min
          END IF

        END IF
      END DO
    END DO

    max_min = max_min_temp
    avg_dist = (at_dist - all_bonds)/kov_dist_sum

    !!! Number of Rotations

!    numrot = -0.2*(log10(real(rot_pun)) / log10(20d0))
    numrot = 0

    IF (this_case /= 'oneadd') THEN

        !!! Second add distance !!!

        add_dist = sqrt((coord(sadd(1),1)-coord(sadd(2),1))**2 + (coord(sadd(1),2)-coord(sadd(2),2))**2 + &
                        (coord(sadd(1),3)-coord(sadd(2),3))**2)

        !!! Reaction vector alignment !!!
        CALL reac_vec(natoms, coord, atomic_numbers, con, sadd(1), rvec1, break_bond, planar1, stereo1)
        CALL reac_vec(natoms, coord, atomic_numbers, con, sadd(2), rvec2, break_bond, planar2, stereo2)

        temp_val = rvec1(1)*rvec2(1)+rvec1(2)*rvec2(2)+rvec1(3)*rvec2(3)

        IF (temp_val >= 1) THEN         !!! same direction
            angle_diff = 0d0
        ELSE IF (temp_val <= -1) THEN   !!! opposite direction
            angle_diff = 180d0
        ELSE
            angle_diff = (180/pi)*acos(temp_val)
        END IF

        IF ((120 < angle_diff) .and. (angle_diff < 240)) THEN   !.OR. &
           angle_sub = abs(180-angle_diff)/60
        ELSE IF ((-60 < angle_diff) .and. (angle_diff < 60)) THEN
            IF (planar1 .OR. planar2) then                      
                angle_sub = abs(angle_diff)/60                  
            ELSE                                                
                angle_sub = 1d0                                 
            END IF                                             
        ELSE
            angle_sub = 1d0
        END IF

        IF (angle_sub /= angle_sub) THEN
                STOP 'angle_sub equals NaN'
                angle_sub = 0d0
        END IF

        !!! Dummy atom distance !!!

        dummy1 = coord(sadd(1),:) + 0.5*rvec1(:)
        dummy2 = coord(sadd(2),:) + 0.5*rvec2(:)

        dummy = (norm(dummy1-dummy2)/norm(coord(sadd(1),:)-coord(sadd(2),:)) - 1)

        IF (dummy /= dummy) THEN
                !!! STOP 'dummy equals NaN'
                dummy = 0d0
        END IF

        !!! Minimal distance from reaction vector !!!

        cylin = cylin_func(natoms, coord, sadd)

        !!! Assignment of calculated properties to ranking values

        ranking_values(1) = max_min
        ranking_values(2) = avg_dist
        ranking_values(3) = add_dist
        ranking_values(4) = angle_sub
        ranking_values(5) = dummy
        ranking_values(6) = cylin
        ranking_values(7) = numrot

    ELSE IF (this_case == 'oneadd') THEN

        ranking_values(1) = max_min
        ranking_values(2) = avg_dist
        ranking_values(3) = 0
        ranking_values(4) = 0
        ranking_values(5) = 0
        ranking_values(6) = 0
        ranking_values(7) = numrot
    END IF

    counter = counter + 1

END SUBROUTINE ranking_value_sub

END MODULE
