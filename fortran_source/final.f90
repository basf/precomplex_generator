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

PROGRAM gen_conformers_all
  IMPLICIT NONE

END PROGRAM gen_conformers_all


SUBROUTINE file_input(xyzFile, datFile, this_case, fortran_weights, max_mem, rot_bonds_neigh, rot_bonds_path, stereo1, &
stereo2, intramol_small_ring)

    ! Input-Parameters from Python code!

    USE fortran_utils, ONLY : remove_blank, chemical_number

    IMPLICIT NONE
    INTEGER                        :: natom,n_tor,n_del_bonds,i,j,io_error,istat,start
    INTEGER, INTENT(IN)            :: rot_bonds_neigh, rot_bonds_path
    INTEGER, ALLOCATABLE           :: con(:,:), torsions_i(:,:),atomic_number(:),del_bonds(:,:)
    INTEGER                        :: add_i(2), temp1, temp2, break_bond(2), n_break
    CHARACTER(90)                  :: line, xyzFile, datFile, this_case
    CHARACTER(2)                   :: atnam
    REAL(8), ALLOCATABLE           :: coord(:,:),torsions_a(:)
    REAL(8), INTENT(IN)            :: fortran_weights(7), max_mem
    REAL(8)                        :: sadd_thresh
    INTEGER, INTENT(IN)            :: stereo1, stereo2
    LOGICAL, INTENT(IN)            :: intramol_small_ring

    !---------------------------------------------------------------------------------------------------
    temp1 = 0; temp2 = 0

    OPEN (99,file=xyzFile,iostat=io_error,position='REWIND',status='old')
    READ(99,*) natom

    ALLOCATE(atomic_number(natom),stat=istat);   IF(istat/=0) STOP 'allocation error: atomic_number'
    ALLOCATE(coord(natom,3),stat=istat);         IF(istat/=0) STOP 'allocation error: coord'
    ALLOCATE(con(natom,natom+1),stat=istat);     IF(istat/=0) STOP 'allocation error: con'

    READ(99,'(a)') line                            ! Read atom symbols and coordinates
    DO i=1,natom
        READ(99,'(a)') line ! '(a)' READs whole line in
        CALL remove_blank(line)
        READ(line,*) atnam, coord(i,1), coord(i,2), coord(i,3)
        atomic_number(i)=chemical_number(atnam)
    END DO

    con=0
    OPEN (86,file=datFile,iostat=io_error,position='REWIND',status='old')
    DO i=1,natom
        READ(86,'(a)') line
        CALL remove_blank(line)
        line=trim(line(index(line,' '):)) ! delete i 
        CALL remove_blank(line)
        line=trim(line(index(line,' '):)) ! delete atnam(i)
        READ(line,*) con(i,1) 
        DO j=2,con(i,1)+1
            CALL remove_blank(line)
            line=trim(line(index(line,' '):)) ! delete con(i,j)
            READ(line,*) con(i,j)
        END DO
    END DO

    CLOSE(86)

    ! Do bonds exist which shall be ignored?
    REWIND(99)
    READ(99,*,iostat=io_error) line
    CALL remove_blank(line)
    DO WHILE (line /= '$del_bond' .AND. io_error >= 0)
        READ(99,*,iostat=io_error) line
        CALL remove_blank(line)
    END DO
    n_del_bonds=0
    IF ( io_error >= 0) THEN
        READ(99,*,iostat=io_error) line
        IF ( io_error < 0 ) STOP 'ERROR: Delete bond specification does not end with $end.'
        CALL remove_blank(line)
        DO WHILE (line /= '$end' .AND. io_error >= 0)
            n_del_bonds=n_del_bonds+1
            READ(99,*) line
            CALL remove_blank(line)
        END DO
    END IF

    ALLOCATE(del_bonds(n_del_bonds,2),stat=istat); IF(istat/=0) STOP 'allocation error: del_bonds'

    del_bonds=0
    IF ( io_error >= 0) THEN
        REWIND(99)
        READ(99,*) line
        CALL remove_blank(line)
        DO WHILE (line /= '$del_bond')
            READ(99,*) line
            CALL remove_blank(line)
        END DO

        DO i=1,n_del_bonds
            READ(99,*) del_bonds(i,1), del_bonds(i,2)
        END DO
    END IF

    ! Read in bonds which might get broken
    REWIND(99)
    READ(99,*,iostat=io_error) line
    CALL remove_blank(line)
    DO WHILE (line /= '$breaks' .AND. io_error >= 0)
        READ(99,*,iostat=io_error) line
        CALL remove_blank(line)
    END DO
    
    IF ( io_error >= 0) THEN
        READ(99,*,iostat=io_error) line
        IF ( io_error < 0 ) STOP 'ERROR: Break bond specification does not end with $end.'
        CALL remove_blank(line)
        DO WHILE (line /= '$end' .AND. io_error >= 0)
            n_break=n_break+1
            READ(99,*) line
            CALL remove_blank(line)
        END DO
    END IF

    break_bond=0

    IF ( io_error >= 0) THEN
        REWIND(99)
        READ(99,*) line
        CALL remove_blank(line)
        DO WHILE (line /= '$breaks')
            READ(99,*) line
            CALL remove_blank(line)
        END DO

        READ(99,*) temp1, temp2
        IF (((abs(temp1) < 1000) .AND. (temp1 > 0)) .AND. ((abs(temp2) < 1000) .AND. (temp2 > 0))) THEN
            break_bond(1) = temp1
            break_bond(2) = temp2
        END IF
    END IF

    ! Read in torsional angles used for rotational sampling
    REWIND(99)
    READ(99,*,iostat=io_error) line
    CALL remove_blank(line)
    DO WHILE (line /= '$tors' .AND. io_error >= 0)
        READ(99,*,iostat=io_error) line
        CALL remove_blank(line)
    END DO
    IF ( io_error < 0) STOP 'ERROR: $tors was not found in input file.'

    n_tor=0
    READ(99,*,iostat=io_error) line
    CALL remove_blank(line)
    DO WHILE (line /= '$end' .AND. io_error >= 0)
        READ(99,*) line
        CALL remove_blank(line)
        n_tor=n_tor+1
    END DO
    IF ( io_error < 0 ) STOP 'ERROR: Torsion specification does not end with $end.'

    IF (this_case /= 'intramol') THEN
        n_tor=n_tor+1
        start=2
    ELSE
        start=1
    END IF

    ALLOCATE(torsions_a(n_tor),stat=istat);   IF(istat/=0) STOP 'allocation error: torsions_a'
    ALLOCATE(torsions_i(n_tor,2),stat=istat); IF(istat/=0) STOP 'allocation error: torsions_i'
    torsions_a=0
    torsions_i=0

    REWIND(99)
    READ(99,*) line
    CALL remove_blank(line)
    DO WHILE (line /= '$tors')
        READ(99,*) line
        CALL remove_blank(line)
    END DO

    DO i=start,n_tor
        READ(99,*) torsions_i(i,1), torsions_i(i,2), torsions_a(i)
    END DO

    IF (this_case == 'intramol') THEN

        ! Read in second/intramolecular add bond

        REWIND(99)
        READ(99,*,iostat=io_error) line
        CALL remove_blank(line)
        DO WHILE (line /= '$add' .AND. io_error >= 0)
            READ(99,*,iostat=io_error) line
            CALL remove_blank(line)
        END DO
        IF ( io_error < 0) STOP 'ERROR: $add was not found in input file.'

        READ(99,*,iostat=io_error) line
        CALL remove_blank(line)
        DO WHILE (line /= '$end' .AND. io_error >= 0)
            READ(99,*) line
            CALL remove_blank(line)
        END DO
        IF ( io_error < 0 ) STOP 'ERROR: Second Add specification does not end with $end.'

        REWIND(99)
        READ(99,*) line
        CALL remove_blank(line)
        DO WHILE (line /= '$add')
            READ(99,*) line
            CALL remove_blank(line)
        END DO
        READ(99,*) add_i(1), add_i(2), sadd_thresh

    ELSE IF ((this_case == 'twoadds') .OR. (this_case == 'shuttle')) THEN

        ! Read in first add bond, around which the precomplex generator shall rotate

        REWIND(99)
        READ(99,*,iostat=io_error) line
        CALL remove_blank(line)
        DO WHILE (line /= '$fadd' .AND. io_error >= 0)
            READ(99,*,iostat=io_error) line
            CALL remove_blank(line)
        END DO
        IF ( io_error < 0) STOP 'ERROR: $fadd was not found in input file.'

        READ(99,*,iostat=io_error) line
        CALL remove_blank(line)
        DO WHILE (line /= '$end' .AND. io_error >= 0)
            READ(99,*) line
            CALL remove_blank(line)
        END DO
        IF ( io_error < 0 ) STOP 'ERROR: First Add specification does not end with $end.'

        REWIND(99)
        READ(99,*) line
        CALL remove_blank(line)
        DO WHILE (line /= '$fadd')
            READ(99,*) line
            CALL remove_blank(line)
        END DO
        READ(99,*) torsions_i(1,1), torsions_i(1,2), torsions_a(1)

        ! Read in second/intramolecular add bond

        REWIND(99)
        READ(99,*,iostat=io_error) line
        CALL remove_blank(line)
        DO WHILE (line /= '$add' .AND. io_error >= 0)
            READ(99,*,iostat=io_error) line
            CALL remove_blank(line)
            !WRITE(*,*) 'tors: ', line
        END DO
        IF ( io_error < 0) STOP 'ERROR: $add was not found in input file.'

        READ(99,*,iostat=io_error) line
        CALL remove_blank(line)
        DO WHILE (line /= '$end' .AND. io_error >= 0)
            READ(99,*) line
            CALL remove_blank(line)
        END DO
        IF ( io_error < 0 ) STOP 'ERROR: Second Add specification does not end with $end.'

        REWIND(99)
        READ(99,*) line
        CALL remove_blank(line)
        DO WHILE (line /= '$add')
            READ(99,*) line
            CALL remove_blank(line)
        END DO
        READ(99,*) add_i(1), add_i(2), sadd_thresh

    ELSE IF (this_case == 'oneadd') THEN

        ! Read in "add bond", around which the precomplex generator shall rotate

        REWIND(99)
        READ(99,*,iostat=io_error) line
        CALL remove_blank(line)
        DO WHILE (line /= '$fadd' .AND. io_error >= 0)
            READ(99,*,iostat=io_error) line
            CALL remove_blank(line)
        END DO
        IF ( io_error < 0) STOP 'ERROR: $fadd was not found in input file.'

        READ(99,*,iostat=io_error) line
        CALL remove_blank(line)
        DO WHILE (line /= '$end' .AND. io_error >= 0)
            READ(99,*) line
            CALL remove_blank(line)
        END DO
        IF ( io_error < 0 ) STOP 'ERROR: First Add specification does not end with $end.'

        REWIND(99)
        READ(99,*) line
        CALL remove_blank(line)
        DO WHILE (line /= '$fadd')
            READ(99,*) line
            CALL remove_blank(line)
        END DO
        READ(99,*) torsions_i(1,1), torsions_i(1,2), torsions_a(1)
        add_i = 1
        sadd_thresh = 6d0

    END IF

    CLOSE(99)

    CALL gen_conformers(natom,atomic_number,coord,con,n_tor,torsions_i,torsions_a,n_del_bonds,del_bonds,add_i,&
            this_case,start,break_bond, fortran_weights, max_mem, sadd_thresh, rot_bonds_neigh, rot_bonds_path, &
            stereo1, stereo2, intramol_small_ring)

END SUBROUTINE file_input


SUBROUTINE gen_conformers(natom,atomic_number,coord,con,n_torsions,torsions_i,torsions_a,n_del_bonds,del_bonds,add_i,&
                          this_case,start,break_bond, fortran_weights, max_mem, sadd_thresh, rot_bonds_neigh, &
                          rot_bonds_path, stereo1, stereo2, intramol_small_ring)

    ! -------------------------------------------------------------------------------------------------------------
    !  Function that generates and evaluates conformations of the initial structure(s) and provides the precomplex 
    !  structure (which is the "best suited" input structure for e.g. SE-MGSM according to the given parametrizations)
    ! -------------------------------------------------------------------------------------------------------------

    USE fortran_utils, ONLY: pi, ranking_value_sub, rotation, kov_rad_num, write_molecule, ranking_func, write_movie, reac_vec
    IMPLICIT NONE

    INTEGER, INTENT(IN)               :: natom, add_i(2),con(natom, natom + 1),start,n_torsions
    INTEGER, INTENT(IN)               :: rot_bonds_neigh, rot_bonds_path
    INTEGER, INTENT(IN)               :: atomic_number(natom),n_del_bonds,del_bonds(n_del_bonds,2)
    INTEGER, INTENT(IN)               :: break_bond(2), torsions_i(n_torsions, 2)
    REAL(8), INTENT(IN)               :: coord(natom, 3),torsions_a(n_torsions), fortran_weights(7), max_mem
    REAL(8), INTENT(IN)               :: sadd_thresh

    CHARACTER(90), INTENT(IN)         :: this_case
    CHARACTER(90)                     :: strategy
    CHARACTER(50)                     :: spacing
    INTEGER                           :: n_tor,n_conf,n_conf_old,i_conf,j_w,i,j,k,l,n,nt,istat,counter,stereo1,stereo2
    INTEGER                           :: grid,conf_num, write_conf, angle_count(50), angle_id, rot_pun
    INTEGER, ALLOCATABLE              :: con_new(:,:),rot(:),natom_f(:),natom_t(:)
    INTEGER, ALLOCATABLE              :: natom_f_new(:),natom_t_new(:)
    INTEGER, ALLOCATABLE              :: n_w(:)
    REAL(8), ALLOCATABLE              :: coord_new(:,:),vec(:),coord_tmp(:,:), final_coord(:,:)
    REAL(8), ALLOCATABLE              :: w(:),dw(:),dw_new(:),O1(:,:),O2(:,:),O3(:,:),O4(:,:)
    REAL(8), ALLOCATABLE              :: mag(:),vec1(:),vec2(:), angle_max(:)
    REAL(8)                           :: dist,dist_new,phi,radius,min_dist,dist_sum
    REAL(8)                           :: max_min_dist, max_dist_sum, dummy1(3), dummy2(3)
    REAL(8)                           :: diff,min_dumm, cl_rad, ranking_values(7)
    REAL(8)                           :: min_clash,all_bonds,dist_diff,rvec_angle,min_add_dist
    REAL(8)                           :: max_ranking,ranking,cylin_out, temp_val
    REAL(8)                           :: angles_size, angle_sub, all_dist, all_dist_bond
    REAL(8)                           :: ranking_values_avg(7), ranking_values_add(7), temp, rand
    REAL(8)                           :: ranking_values_min(7), ranking_values_max(7)
    INTEGER (KIND=2), ALLOCATABLE     :: angles(:,:)
    INTEGER (KIND=2)                  :: angle_incr(50)
    LOGICAL                           :: ex, planar, too_large, bonds_calc, replaced, exist
    LOGICAL, ALLOCATABLE              :: tor(:,:), tor_new(:,:), ring(:)
    LOGICAL                           :: intramol_small_ring


    integer :: values(1:8), kk
    integer, dimension(:), allocatable :: seed
    real(8) :: r

    ALLOCATE(vec(3),stat=istat);      IF(istat/=0) STOP 'allocation error: vec'
    ALLOCATE(vec1(3),stat=istat);     IF(istat/=0) STOP 'allocation error: vec1'
    ALLOCATE(vec2(3),stat=istat);     IF(istat/=0) STOP 'allocation error: vec2'



    write_conf = 0; grid = 1; ranking_values = 0; angle_id = 0
    ranking_values_add = 0; ranking_values_avg = 0; replaced = .FALSE.; angle_count(:) = 0
    ranking_values_min(:) = 1000; ranking_values_max(:) = 0; max_ranking = 0; angle_incr(:) = 0

    !---------------------------------------------------------------------------------------------------

        ALLOCATE(rot(natom),stat=istat);             IF(istat/=0) STOP 'allocation error: rot'
        ALLOCATE(coord_tmp(natom+1,3),stat=istat);     IF(istat/=0) STOP 'allocation error: coord_tmp'
        ALLOCATE(con_new(natom,natom+1),stat=istat); IF(istat/=0) STOP 'allocation error: con_new'


        con_new=0

        n_tor = n_torsions

        ALLOCATE(tor(n_tor,natom),stat=istat);   IF(istat/=0) STOP 'allocation error: tor'
        ALLOCATE(ring(n_tor),stat=istat);        IF(istat/=0) STOP 'allocation error: ring'
        ALLOCATE(natom_f(n_tor),stat=istat);     IF(istat/=0) STOP 'allocation error: natom_f'
        ALLOCATE(natom_t(n_tor),stat=istat);     IF(istat/=0) STOP 'allocation error: natom_t'
        ALLOCATE(dw(n_tor),stat=istat);          IF(istat/=0) STOP 'allocation error: dw'

        natom_f = torsions_i(:, 1) ! Vector which contains the index of the first atom of the $tors section
        natom_t = torsions_i(:, 2) ! Vector which contains the index of the second atom of the $tors section
        dw = torsions_a            ! Vector with the float numbers of the torsional angles


        DO i=start,n_tor    ! Index i = 1 is reserved for first add
            ex=.TRUE.
            DO j=2, con(natom_f(i), 1) + 1                      ! Are natom_f and natom_t connected by a bond?
                IF (con(natom_f(i), j) == natom_t(i)) THEN
                    ex=.FALSE.
                    EXIT
                ENDIF
            END DO
            IF (ex) THEN
                WRITE(*,*) 'ERROR: The atoms ',natom_f(i),' and ',natom_t(i),' are not connected via a direct bond!'
                WRITE(*,*) 'Please select only directly bonded atoms here.'
                STOP
            END IF
        END DO

        ! Check for torsions inside a ring
        tor=.FALSE.
        ring=.TRUE.
        DO k=1,n_tor
            nt=1                            ! Which atoms have to be rotated? -> tor
            rot=0
            rot(nt)=natom_t(k)
            tor(k,rot(nt))=.TRUE.
            DO i=2,con(natom_t(k),1)+1
                IF (con(natom_t(k),i) /= natom_f(k)) THEN
                    nt=nt+1
                    rot(nt)=con(natom_t(k),i)
                    tor(k,rot(nt))=.TRUE.
                END IF
            END DO

            n=1
            l1:  DO WHILE (n < nt)
                n=n+1
                DO i=2,con(rot(n),1)+1
                    ex=.TRUE.
                    DO j=1,nt
                        IF (con(rot(n),i) == rot(j)) THEN
                            ex=.FALSE.
                            EXIT
                        END IF
                    END DO
                    IF (ex) THEN
                        IF (con(rot(n),i) /= natom_f(k)) THEN
                            nt=nt+1
                            rot(nt)=con(rot(n),i)
                            tor(k,rot(nt))=.TRUE.
                        ELSE
                            ring(k)=.FALSE.
                            WRITE(*,*) 'WARNING: You have asked for a torsion inside a ring!'
                            STOP 'Torsion inside ring'
                            EXIT l1
                        END IF
                    END IF
                END DO
            END DO l1
            tor(k,rot(1))=.FALSE.
            nt=nt-1
        END DO

        n=0
        DO k=1,n_tor
            IF (ring(k)) THEN
                n=n+1
            END IF
        END DO

        ALLOCATE(tor_new(n,natom),stat=istat);          IF(istat/=0) STOP 'allocation error: tor_new'
        ALLOCATE(natom_f_new(n),stat=istat);            IF(istat/=0) STOP 'allocation error: natom_f_new'
        ALLOCATE(natom_t_new(n),stat=istat);            IF(istat/=0) STOP 'allocation error: natom_t_new'
        ALLOCATE(dw_new(n),stat=istat);                 IF(istat/=0) STOP 'allocation error: dw_new'

        n=0
        DO k=1,n_tor
            IF (ring(k)) THEN
                n=n+1
                tor_new(n,:)=tor(k,:)
                natom_f_new(n)=natom_f(k)
                natom_t_new(n)=natom_t(k)
                dw_new(n)=dw(k)
            END IF
        END DO

        n_tor=n

        DEALLOCATE(tor, natom_f, natom_t, dw)
        ALLOCATE(tor(n_tor,natom),stat=istat);              IF(istat/=0) STOP 'allocation error: tor'
        ALLOCATE(natom_f(n_tor),stat=istat);                IF(istat/=0) STOP 'allocation error: natom_f'
        ALLOCATE(natom_t(n_tor),stat=istat);                IF(istat/=0) STOP 'allocation error: natom_t'
        ALLOCATE(w(n_tor),stat=istat);                      IF(istat/=0) STOP 'allocation error: w'
        ALLOCATE(dw(n_tor),stat=istat);                     IF(istat/=0) STOP 'allocation error: dw'
        ALLOCATE(n_w(n_tor),stat=istat);                    IF(istat/=0) STOP 'allocation error: n_w'
        ALLOCATE(mag(n_tor),stat=istat);                    IF(istat/=0) STOP 'allocation error: mag'
        ALLOCATE(coord_new(natom,3),stat=istat); IF(istat/=0)         STOP 'allocation error: coord_new'
        ALLOCATE(O1(natom,3),stat=istat); IF(istat/=0)                STOP 'allocation error: O1'
        ALLOCATE(O2(natom,3),stat=istat); IF(istat/=0)                STOP 'allocation error: O2'
        ALLOCATE(O3(natom,3),stat=istat); IF(istat/=0)                STOP 'allocation error: O3'
        ALLOCATE(O4(natom,3),stat=istat); IF(istat/=0)                STOP 'allocation error: O4'
        ALLOCATE(angle_max(n_tor),stat=istat); IF(istat/=0)           STOP 'allocation error: angle_max'

        O1 = 0; O2 = 0; O4 = 0; O3 = 0; O4 = 0; min_add_dist = 999; max_ranking = 0; ranking = 0; angle_sub = 1
        planar = .FALSE.; min_dist = 999; max_min_dist = 0; max_dist_sum = 0; min_dumm = 999; min_clash=0; cl_rad = 0.75
        too_large=.FALSE.; temp = 1; rand = 0

        tor=tor_new
        natom_f=natom_f_new
        natom_t=natom_t_new
        dw=dw_new

        DEALLOCATE(tor_new, natom_f_new, natom_t_new, dw_new)

        IF (n_tor > 0) THEN
            n_conf=1
            DO i=1, n_tor
                n_conf = n_conf * ceiling(360/dw(i))
            END DO

            angles_size = (real(n_conf)*real(n_tor)*2d0) / (1E9)

            IF (angles_size > 12d0) THEN
                WRITE(*,'(A45, F6.3, A4)') 'RAM usage for angles-Matrix too large: ', angles_size, ' GB.'
                WRITE(*,'(A18, F6.3, A4)') 'Available memory: ', 1.2*max_mem, ' GB.'
                WRITE(*,'(A45, F6.3, A3)') 'Will take initial structure.'
                too_large = .TRUE.
            END IF
        ELSE
            n_conf=1   ! No conformers will be generated, only the input structure survives.
            ALLOCATE(angles(n_conf,1),stat=istat); IF(istat/=0) STOP 'First allocation error: angles'
        END IF

        IF (.NOT. too_large) THEN
            IF (n_tor > 0.1) THEN
                ALLOCATE(angles(n_conf,n_tor),stat=istat); IF(istat/=0) STOP 'Third allocation error: angles'
            END IF
        ELSE
            n_conf=1
            ALLOCATE(angles(n_conf,1),stat=istat); IF(istat/=0) STOP 'Fourth allocation error: angles'
        END IF

        ! Construct a list of all permutation of all possible values for all torsion angles

        IF (n_conf > 1) THEN
            n_conf=1
            angles(n_conf,:)=0

            n_w(:)=0
            DO j_w=1,n_tor
                n_conf_old=n_conf
                DO i=1,n_conf_old
                    w(j_w)=0
                    DO WHILE (w(j_w) < 360-dw(j_w))
                        w(j_w)=w(j_w)+dw(j_w)
                        n_conf=n_conf+1
                        DO j=1,n_tor
                            IF (j /= j_w) THEN
                                angles(n_conf,j)=angles(i,j)
                            ELSE
                                angles(n_conf,j)=int(w(j_w))
                            END IF
                        END DO
                    END DO
                END DO
            END DO

            ! Create conformers
            DO j_w=1,n_tor
                vec(:)=coord(natom_t(j_w),:)-coord(natom_f(j_w),:)
                mag(j_w)=sqrt((vec(1)**2)+(vec(2)**2)+(vec(3)**2))
            END DO

            IF (n_conf <= nint(1/log10(natom+1d0) * 5E4)) THEN
                grid = 1
            ELSE
                temp = n_conf/grid

                DO WHILE (temp > nint(1/log10(natom+1d0) * 5E4))
                    IF (mod(n_conf, 3) < 1E-3) THEN
                        grid = grid * 3
                    ELSE IF (mod(n_conf, 2) < 1E-3) THEN
                        grid = grid * 2
                    END IF
                    temp = n_conf/grid
                END DO
            END IF

            ! WRITE(*,'(A8, I20, A8, I8, A8, I12)') 'CONF', n_conf, 'GRID', grid, 'NCALC', n_conf/grid

            ! EXTRA-OUTPUT for conformer sampling

            inquire(file="conf_info.txt", exist=exist)
            if (exist) then
                close(74, status='delete')
            end if

            OPEN(74, file="conf_info.txt", status="new", action="write")
            WRITE(74,'(3I20.1)') n_conf, grid, n_conf/grid

            ALLOCATE(final_coord(natom,3),stat=istat); IF(istat/=0)         STOP 'allocation error: final_coord'

            i_conf=1; counter = 1; all_bonds=0; bonds_calc = .FALSE.; final_coord(:,:) = 0

main_loop: DO k=1,n_conf,grid

                IF (n_tor > 0) THEN
                    coord_tmp(:,:)=coord(:,:)
                    DO j_w=1,n_tor
                        phi=angles(k,j_w)*pi/180.0
                        vec(:)=coord_tmp(natom_t(j_w),:)-coord_tmp(natom_f(j_w),:)
                        vec=vec/mag(j_w)
                        DO i=1,natom
                            coord_new(i,:)=coord_tmp(i,:)-coord_tmp(natom_t(j_w),:)
                            IF (tor(j_w,i)) THEN
                                IF (phi /= 0) THEN
                                    call rotation(coord_new(i,:),phi,vec)
                                END IF
                            END IF
                        END DO
                        coord_tmp(:,:)=coord_new(:,:)
                    END DO
                ELSE
                    coord_new(:,:)=coord(:,:)
                END IF

                con_new=0
                ex=.TRUE.

        l2:     DO i=1,natom
                    n=1
                    DO j=1,natom
                        IF (i /= j) THEN

                            CALL kov_rad_num(atomic_number(i), atomic_number(j), radius)

                            dist=(coord_new(i,1)-coord_new(j,1))**2
                            dist=dist+((coord_new(i,2)-coord_new(j,2))**2)
                            dist=dist+((coord_new(i,3)-coord_new(j,3))**2)
                            dist=sqrt(dist)

                            ! check if newly generated conformation is "valid"/reasonable

                            IF (dist < radius .AND. dist > 0.0 ) THEN

                                n=n+1
                                con_new(i,1)=con_new(i,1)+1
                                con_new(i,n)=j

                                IF (.NOT.(ANY(con(i, 2:) == con_new(i, n)))) THEN
                                    ex=.FALSE.                                 ! Steric clash
                                ELSE
                                    IF ((j>i) .AND. (.NOT. bonds_calc)) THEN   ! do bonds_calc only once
                                        all_bonds = all_bonds + dist
                                    END IF
                                END IF

                                IF (.NOT. ex) THEN
                                    IF (n_del_bonds > 0) THEN
                                        DO l=1,n_del_bonds
                                            IF (i == del_bonds(l,1)) THEN
                                                IF (j == del_bonds(l,2)) THEN
                                                    ex=.TRUE.
                                                END IF
                                            END IF
                                            IF (i == del_bonds(l,2)) THEN
                                                IF (j == del_bonds(l,1)) THEN
                                                    ex=.TRUE.
                                                END IF

                                            END IF
                                        END DO
                                    END IF
                                END IF

                                IF (.NOT. ex) THEN
                                    EXIT l2
                                END IF
                            END IF
                        END IF
                    END DO
                END DO l2

                bonds_calc = .TRUE.

                IF (ex) THEN

                    !!! OBSOLETE - Calculation of number of tors for rot-punishment
                    DO angle_id = 1, n_tor
                        IF (angles(k, angle_id) > 1) THEN
                            angle_incr(angle_id) = 1
                        ELSE IF (angles(k, angle_id) < 1) THEN
                            angle_incr(angle_id) = 0
                        END IF
                    END DO
                    rot_pun = sum(angle_incr)
                    !!!

                    CALL ranking_value_sub(coord_new(:,:), con(:,:), atomic_number(:), fortran_weights(:), &
                                           this_case, counter, all_bonds, add_i, break_bond, ranking_values(:), rot_pun, &
                                           stereo1, stereo2)

                    IF (i_conf == 1) THEN
                        ranking_values_add(:) = ranking_values(:)
                        ranking_values_avg(:) = ranking_values(:)
                    END IF

                    CALL ranking_func(this_case, ranking_values(:), i_conf, ranking_values_add(:), coord_new(:,:), &
                          fortran_weights(:), con(:,:),  ranking_values_avg(:), ranking_values_min(:), &
                          ranking_values_max(:), max_ranking, final_coord(:,:), replaced, break_bond(:),&
                          add_i(:), atomic_number(:), sadd_thresh, rot_bonds_neigh, rot_bonds_path, intramol_small_ring)

                    IF (replaced) THEN
                        CALL write_movie(natom, atomic_number(:), final_coord(:,:), i_conf) 
                    END IF

                    i_conf=i_conf+1
                END IF

            END DO main_loop


            IF (max_ranking > 1E-3) THEN
                CALL write_molecule(natom, atomic_number(:), final_coord(:,:), nint(max_ranking*100))

                ! FURTHER OUTPUT FOR ANALYSIS !
                inquire(file="best_rank.txt", exist=exist)
                if (exist) then
                    open(79, file="best_rank.txt", status="old", position="append", action="write")
                else
                    open(79, file="best_rank.txt", status="new", action="write")
                end if

                close(79)

            END IF

            DEALLOCATE(con_new,rot,natom_f,natom_t,n_w,coord_new)
            DEALLOCATE(vec,coord_tmp,w,dw,angles,vec1,vec2,O1,O2,O3,O4)

        ELSE
            CALL write_molecule(natom, atomic_number(:), coord(:,:), 1)
        END IF

END SUBROUTINE gen_conformers
