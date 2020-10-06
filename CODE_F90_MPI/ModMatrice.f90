module ModMatrice
  use ModMesh
  use MPI
  implicit none

  integer,parameter::nbr_point_par_ele=3,nombre_node=10181,nombre_ddl=2
  



contains


  !JACOBIENNE DE LA TRANSFORMATION
  subroutine JACOBIENNE(T_node,T_conn,jac,i,nbr_tri)
    integer,intent(in)::i,nbr_tri
    TYPE(tab_node),intent(in)::T_node
    integer,dimension(1:nbr_tri,1:5),intent(in)::T_conn
    real(rki),dimension(1:2,1:2),intent(out)::jac
    integer::P1,P2,P3
    !i : itérant sur la liste des connectivité des élé triangle
    !P numero du noeud de l'élément
    jac=0.d0
    P1=T_conn(i,2);P2=T_conn(i,3);P3=T_conn(i,4)

    jac(1,1)=T_node%c_node(P2,1)-T_node%c_node(P1,1);jac(1,2)=T_node%c_node(P2,2)-T_node%c_node(P1,2)

    jac(2,1)=T_node%c_node(P3,1)-T_node%c_node(P1,1);jac(2,2)=T_node%c_node(P3,2)-T_node%c_node(P1,2)


  end subroutine JACOBIENNE

  !DETERMINANT DE J
  subroutine DETERMINANT_JAC(jac,det)
    real(rki),dimension(1:2,1:2),intent(in)::jac
    real(rki),intent(out)::det
    det=0.
    det=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
  end subroutine DETERMINANT_JAC

  !INVERSE DE LA JACOBIENNE:
  subroutine INVERSE_JAC(jac,det,inv_jac)
    real(rki),dimension(1:2,1:2),intent(in)::jac
    real(rki),dimension(1:2,1:2),intent(out)::inv_jac
    real(rki),intent(in)::det
    real(rki)::inv_det
    inv_jac=0.
    !call DETERMINANT_JAC(jac,det)
    inv_det=1/det

    inv_jac(1,1)=jac(2,2);inv_jac(2,1)=-jac(2,1)
    inv_jac(1,2)=-jac(1,2);inv_jac(2,2)=jac(1,1)

    inv_jac=inv_det*inv_jac
  end subroutine INVERSE_JAC

  !MATRICE MODIFIEE DE INV_J POUR CORRESPONDRE A LA FORMULATION DU PB
  subroutine Q_INV_J(inv_jac,Q)
    real(rki),dimension(1:2,1:2),intent(in)::inv_jac
    real(rki),dimension(1:3,1:4),intent(out)::Q
    Q=0.
    Q(1,1)=inv_jac(1,1);Q(1,2)=inv_jac(1,2);Q(1,3:4)=0.
    Q(2,1:2)=0.;Q(2,3)=inv_jac(2,1);Q(2,4)=inv_jac(2,2)
    Q(3,1)=inv_jac(2,1);Q(3,2)=inv_jac(2,2);Q(3,3)=inv_jac(1,1);Q(3,4)=inv_jac(1,2)
  end subroutine Q_INV_J

  !RETOURNE B i.e Q.DN avec DN derivé fct de forme dans base ele ref
  subroutine GRAD_N(Q,DN,B)
    real(rki),dimension(1:3,1:4),intent(in)::Q
    real(rki),dimension(1:4,1:6),intent(in)::DN
    real(rki),dimension(1:3,1:6),intent(inout)::B !B:=OPERATEUR DERIVATION dans la base de ref :epsilon=B.dpl(élé)
    integer::i,j,k
    B=0.
    DO i=1,3
       DO j=1,6
          DO k=1,4
             B(i,j)=B(i,j)+Q(i,k)*DN(k,j)
          END DO
       END DO
    END DO
    ! ou B=MATMUL(Q,DN)
  end subroutine GRAD_N

  !MATRICE DE RIGIDITE ELEMENTAIRE DE l'ELEMENT CONSIDERE:
  subroutine MATRICE_RIGIDITE_ELEMENT(B,HOOK,det,K_ele)
    real(rki),dimension(1:3,1:6),intent(in)::B
    real(rki),dimension(1:3,1:3),intent(in)::HOOK
    real(rki),intent(in)::det
    real(rki),dimension(1:6,1:6),intent(out)::K_ele
    real(rki),dimension(1:6,1:3)::inter
    real(rki),dimension(1:6,1:6)::trans
    real(rki)::ep
    integer::i,j,k
    ep=1
    K_ele=0.
    inter=0.
    !inter=tranpose(B)*HOOK:
    DO j=1,6
       DO i=1,3
          DO k=1,3
             inter(j,i)=inter(j,i)+B(k,j)*HOOK(k,i)
          END DO
       END DO
    END DO

    !SINON FAIRE inter=MATMUL(TRANSPOSE(B),HOOK)
    !k_ele=inter*B:
    DO i=1,6
       DO j=1,6
          DO k=1,3
             K_ele(i,j)=K_ele(i,j)+inter(i,k)*B(k,j)
          END DO
       END DO
    END DO
    !SINON FAIRE K_ele=MATMUL(inter,B)

    K_ele=0.5*ABS(det)*K_ele!*ep
!!$trans=TRANSPOSE(K_ele)
!!$PRINT*,' '
!!$PRINT*,'K_ele'
!!$DO i=1,6
!!$   PRINT*,K_ele(i,:)!-trans(i,:)
!!$END DO

  end subroutine MATRICE_RIGIDITE_ELEMENT

  SUBROUTINE DPL_IMPO(K_ele,T_conn,T_node,nbr_tri,i,condi_x,condi_y)
    real(rki),dimension(1:6,1:6),intent(inout)::K_ele
    integer,intent(in)::i,nbr_tri
    TYPE(tab_node),intent(in)::T_node
    integer,dimension(1:nbr_tri,1:5),intent(in)::T_conn
    real(rki),intent(in)::condi_x,condi_y
    integer::P1,P2,P3,k1,k,k2
    real(rki)::a,b
    !i : itérant sur la liste des connectivité des élé triangle
    !P numero du noeud de l'élément
    P1=T_conn(i,2);P2=T_conn(i,3);P3=T_conn(i,4)
    DO k=2,4
       IF((T_node%c_node(T_conn(i,k),1)==condi_x))THEN !dpl bloqué en x=condi_x, i.e dpl(x)=0:
          k1=2*k-3 !bijection
          a=K_ele(k1,k1);K_ele(k1,:)=0.0d0;K_ele(:,k1)=0.0d0;K_ele(k1,k1)=1.0d0!a*nbr_tri
          PRINT*,'CONDI_X',T_node%c_node(T_conn(i,k),1)
       ELSE IF((T_node%c_node(T_conn(i,k),2)==condi_y))THEN !dpl bloqué en y=condi_y, i.e dpl(y)=0:
          k2=2*k-2 !bijection
          b=K_ele(k2,k2);K_ele(k2,:)=0.0d0;k_ele(:,k2)=0.0d0;K_ele(k2,k2)=1.0d0!b*nbr_tri
          PRINT*,'CONDI_Y',T_node%c_node(T_conn(i,k),2)
       END IF
    END DO
  END SUBROUTINE DPL_IMPO

  !SUBROUTINE DE LOCALISATION:
  subroutine LOC_i(n,NDDL,nbr_tri,T_conn,i,localisation)
    integer,intent(in)::n,NDDL,nbr_tri
    integer,dimension(1:nbr_tri,1:5),intent(in)::T_conn
    integer,intent(in)::i
    integer,dimension(1:n*NDDL),intent(out)::localisation
    integer::j,k,l
    localisation=0
    k=1

    DO WHILE(k<n*NDDL+1)
       DO j=2,4 !1,n de 2à4 car stocke colonne 2à4
          DO l=1,NDDL
             localisation(k)=(T_conn(i,j)-1)*NDDL+l
             
             k=k+1
          END DO
       END DO
    END DO
  end subroutine LOC_i


  !PRE-TRI
  SUBROUTINE PRE_TRI(T1,T2,T3,NBR_C,SIZE,rang,L1,L2,L3,RC1,RC2,RC3,status,stateinfo)
    integer,dimension(MPI_STATUS_SIZE),intent(inout)::status
    INTEGER,intent(in)::T1,T2,T3,NBR_C,SIZE,rang !L1=TAB_CL(1,:), SIZE = la taille de L1
    INTEGER,intent(inout)::stateinfo
    INTEGER,dimension(1:SIZE),intent(in)::L1,L2  !L2=TAB_CL(2,:)
    REAL(rki),dimension(1:SIZE),intent(in)::L3   !L3=TAB_V(:)
    INTEGER,dimension(1:4)::PIV1,PIV2,other
    INTEGER::cpt,t,i,j,he,TAILLE,GAP
    INTEGER,dimension(:),allocatable::EX1,EX2
    real(rki),dimension(:),allocatable::EX3
    INTEGER,dimension(:),allocatable,intent(out)::RC1,RC2
    real(rki),dimension(:),allocatable,intent(out)::RC3
    INTEGER::TAG
    CHARACTER(len=3)::rank
    PIV1(1)=0;PIV1(2)=T1;PIV1(3)=T2;PIV1(4)=T3
    PIV2(1)=T1;PIV2(2)=T2;PIV2(3)=T3;PIV2(4)=NBR_C

    write(rank,fmt='(i3.3)')rang


    DO t=3,0,-1!0,3 SI ON FAIT LA BOUCLE DE t=0,3 IL Y A UN BUG: PIV1(4)==0 ALORS QU'IL VAUT T3
       
       cpt=0
       DO j=1,SIZE
             IF(L1(j)<=PIV2(t+1))THEN
                IF(L1(j)>PIV1(t+1))THEN
                   
                   cpt=cpt+1
                END IF
             END IF
       END DO

       IF(cpt>0)THEN
          ALLOCATE (EX1(1:cpt),EX2(1:cpt),EX3(1:cpt));TAG=1
          
       ELSE
          TAG=0
       END IF
       
       cpt=0
       TAILLE=0;other=0
       IF(rang/=t)THEN
          IF(TAG==1)THEN
             DO i=1,SIZE
                IF((L1(i)>PIV1(t+1)).AND.(L1(i)<=PIV2(t+1)))THEN
                   cpt=cpt+1
                   
                   EX1(cpt)=L1(i);EX2(cpt)=L2(i);EX3(cpt)=L3(i)
                END IF
             END DO
          END IF
          
          call MPI_SEND(cpt,1,MPI_INTEGER,t,10+rang+t,MPI_COMM_WORLD,stateinfo)
          
          IF(TAG==1)THEN
             CALL MPI_SEND(EX1(1:cpt),cpt,MPI_INTEGER,t,30+rang+t,MPI_COMM_WORLD,stateinfo)
             CALL MPI_SEND(EX2(1:cpt),cpt,MPI_INTEGER,t,60+rang+t,MPI_COMM_WORLD,stateinfo)
             CALL MPI_SEND(EX3(1:cpt),cpt,MPI_DOUBLE_PRECISION,t,90+rang+t,MPI_COMM_WORLD,stateinfo)
             DEALLOCATE(EX1,EX2,EX3)
          END IF
          
       ELSE
          IF(TAG==1)THEN
             DO i=1,SIZE
                IF((L1(i)>PIV1(t+1)).AND.(L1(i)<=PIV2(t+1)))THEN
                   cpt=cpt+1
                   
                   EX1(cpt)=L1(i);EX2(cpt)=L2(i);EX3(cpt)=L3(i)
                END IF
             END DO
          END IF
          other(rang)=cpt
          TAILLE=cpt
          DO he=0,3
             IF(he/=t)THEN
                CALL MPI_RECV(other(he),1,MPI_INTEGER,he,10+he+t,MPI_COMM_WORLD,status,stateinfo)
                TAILLE=TAILLE+other(he)
                
             END IF
          END DO
          ALLOCATE(RC1(1:TAILLE),RC2(1:TAILLE),RC3(1:TAILLE))
          IF(TAG==1)THEN
             DO i=1,other(rang)
                RC1(i)=EX1(i);RC2(i)=EX2(i);RC3(i)=EX3(i)
             END DO
             DEALLOCATE(EX1,EX2,EX3)
          END IF
          GAP=other(rang)
          
          DO he=0,3
             IF(other(he)>0.AND.he/=t)THEN!t ou rang)then
                CALL MPI_RECV(RC1(GAP+1:GAP+other(he)),other(he),MPI_INTEGER,he,30+he+t,MPI_COMM_WORLD,status,stateinfo)
                CALL MPI_RECV(RC2(GAP+1:GAP+other(he)),other(he),MPI_INTEGER,he,60+he+t,MPI_COMM_WORLD,status,stateinfo)
                CALL MPI_RECV(RC3(GAP+1:GAP+other(he)),other(he),MPI_DOUBLE_PRECISION,he,90+he+t,MPI_COMM_WORLD,status,stateinfo)
                GAP=GAP+other(he)
             END IF
          END DO
          
       END IF


END DO
PRINT*,'RANG',rang,'FIN BOUCLE t, t=',t

END SUBROUTINE PRE_TRI



!POUR LE FORMAT SPARSE:
  SUBROUTINE MAKE_A_LISTE1(rang,NNZ_LOW_PART,CHARGE_ME,ite,it1,k_ite,n_NDDL,LOC,K_ele,TAB_CL,TAB_V)
    integer,intent(in)::rang,NNZ_LOW_PART,CHARGE_ME,ite,it1,n_NDDL
    !itération actuelle
    !nbr non zero dans la partie triangulaire inferieur de K_ele
    !nbr d'élément traité par chaque proc
    integer,dimension(1:n_NDDL),intent(in)::LOC
    real(rki),dimension(1:n_NDDL,1:n_NDDL),intent(in)::K_ele
    integer,intent(inout)::k_ite
    integer,dimension(1:NNZ_LOW_PART*CHARGE_ME,1:2),intent(out)::TAB_CL
    !STOCK cows rows where NNZ in K_ele_exp
    real(rki),dimension(1:NNZ_LOW_PART*CHARGE_ME),intent(out)::TAB_V
    !STOCK VAL nnz de K_ele_exp partie trian inf
    integer::i,j
    character(len=3)::rank
    WRITE(rank,fmt='(i3.3)')rang
    !integer::lenght
    !lenght=NNZ_LOW_PART*CHARGE_ME
    !TAB_CL=0;TAB_V=0.d0
    IF(ite==it1)THEN
       k_ite=1
    END IF
    !***POUR "simuler" un format CCS SYMETRIQUE
    !k_ite=1
   
    DO j=1,n_NDDL
       DO i=1,n_NDDL !j,n_NDDL
          !IF(LOC(i)>=LOC(j))THEN
             IF(K_ele(i,j)/=0.d0)THEN
                TAB_CL(k_ite,1)=LOC(j);TAB_CL(k_ite,2)=LOC(i)
                TAB_V(k_ite)=K_ele(i,j)
                k_ite=k_ite+1
             ELSE
                TAB_CL(k_ite,1)=0;TAB_CL(k_ite,2)=0
                TAB_V(k_ite)=0.d0
                k_ite=k_ite+1
             END IF
          !END IF
       END DO
    END DO

  END SUBROUTINE MAKE_A_LISTE1
  


  subroutine EXPANSER_RIGIDITE_LOCALE(K_ele,localisation,n,NDDL,Nd,K_ele_expansee)
    integer,intent(in)::n,NDDL,Nd
    integer,dimension(1:n*NDDL),intent(in)::localisation
    real(rki),dimension(1:n*NDDL,1:n*NDDL),intent(in)::K_ele
    real(rki),dimension(1:Nd*NDDL,1:Nd*NDDL),intent(out)::K_ele_expansee
    integer::k,l
    PRINT*,'DANS EXPANSER RIGIDITE LOC'
    K_ele_expansee=0.d0
    DO k=1,n*NDDL
       print*,'ex_ri_loc',' k=',k
       DO l=1,n*NDDL
          !print*,'l=',l
          K_ele_expansee(localisation(k),localisation(l))=K_ele(k,l)
       END DO
    END DO
    PRINT*,'FIN EXP RIG LOC'
  end subroutine EXPANSER_RIGIDITE_LOCALE
  !******************************************************************************************************


  subroutine FORCE_LOCALE_SUIVANT_Y(nume,P,longueur,n,Nd,NDDL,nbr_tri,T_node,T_conn,F_LOCALE,MODEL)
    real(rki),intent(in)::P,longueur
    integer,intent(in)::nume,n,Nd,NDDL,nbr_tri,MODEL
    TYPE(tab_node),intent(in)::T_node
    integer,dimension(1:nbr_tri,1:5),intent(in)::T_conn
    real(rki),dimension(1:n*NDDL),intent(out)::F_LOCALE
    !integer::i,k,j,l
    integer::P1,P2,P3
    real(rki)::X1,X2,X3,Y1,Y2,Y3,ep
    ep=0.00003
    F_LOCALE=0.d0
    P1=T_conn(nume,2);P2=T_conn(nume,3);P3=T_conn(nume,4)
    X1=T_node%c_node(P1,1);X2=T_node%c_node(P2,1);X3=T_node%c_node(P3,1)
    Y1=T_node%c_node(P1,2);Y2=T_node%c_node(P2,2);Y3=T_node%c_node(P3,2)
    

    SELECT CASE(MODEL)
       CASE(0)
          IF(Y1==0. .AND. Y2==0.)THEN
             F_LOCALE(2)=-P*abs(X1-X2)*0.5;F_LOCALE(4)=F_LOCALE(2)
          ELSE IF(Y1==0. .AND. Y3==0.)THEN
             F_LOCALE(2)=-P*abs(X1-X3)*0.5;F_LOCALE(6)=F_LOCALE(2)
          ELSE IF(Y2==0. .AND. Y3==0.)THEN
             F_LOCALE(4)=-P*abs(X2-X3)*0.5;F_LOCALE(6)=F_LOCALE(4)

          ELSE IF(Y1==longueur .AND. Y2==longueur)THEN
             F_LOCALE(2)=P*abs(X1-X2)*0.5;F_LOCALE(4)=F_LOCALE(2)
          ELSE IF(Y1==longueur .AND. Y3==longueur)THEN
             F_LOCALE(2)=P*abs(X1-X3)*0.5;F_LOCALE(6)=F_LOCALE(2)
          ELSE IF(Y2==longueur .AND. Y3==longueur)THEN
             F_LOCALE(4)=P*abs(X2-X3)*0.5;F_LOCALE(6)=F_LOCALE(4)
          END IF
       CASE(1)
          IF(Y1==longueur .AND. Y2==longueur)THEN
             F_LOCALE(2)=P*abs(X1-X2)*0.5;F_LOCALE(4)=F_LOCALE(2)
          ELSE IF(Y1==longueur .AND. Y3==longueur)THEN
             F_LOCALE(2)=P*abs(X1-X3)*0.5;F_LOCALE(6)=F_LOCALE(2)
          ELSE IF(Y2==longueur .AND. Y3==longueur)THEN
             F_LOCALE(4)=P*abs(X2-X3)*0.5;F_LOCALE(6)=F_LOCALE(4)
          END IF
    END SELECT
  end subroutine FORCE_LOCALE_SUIVANT_Y


subroutine FORCE_LOCALE_SUIVANT_X(nume,P,largeur,n,Nd,NDDL,nbr_tri,T_node,T_conn,F_LOCALE)
    real(rki),intent(in)::P,largeur
    integer,intent(in)::nume,n,Nd,NDDL,nbr_tri
    TYPE(tab_node),intent(in)::T_node
    integer,dimension(1:nbr_tri,1:5),intent(in)::T_conn
    real(rki),dimension(1:n*NDDL),intent(out)::F_LOCALE
    !integer::i,k,j,l
    integer::P1,P2,P3
    real(rki)::X1,X2,X3,Y1,Y2,Y3

    F_LOCALE=0.0d0
    P1=T_conn(nume,2);P2=T_conn(nume,3);P3=T_conn(nume,4)
    X1=T_node%c_node(P1,1);X2=T_node%c_node(P2,1);X3=T_node%c_node(P3,1)
    Y1=T_node%c_node(P1,2);Y2=T_node%c_node(P2,2);Y3=T_node%c_node(P3,2)
    !PRINT*,'Y ',Y1,' ',Y2,' ',Y3,'longueur',longueur,' P',P1,' ',P2,' ',P3
    !P=E_L*facteur_de_dpl !159*10**9*0.001

!!$    IF(X1==0.0d0 .AND. X2==0.0d0)THEN
!!$       !PRINT*,'$$$$$$$$$$$$$$$$$$$$$1',-P*abs(Y1-Y2)*0.5,' ',abs(Y1-Y2)
!!$       F_LOCALE(1)=-P*abs(Y1-Y2)*0.5;F_LOCALE(3)=F_LOCALE(1)
!!$    ELSE IF(X1==0.0d0 .AND. X3==0.0d0)THEN
!!$       !PRINT*,'$$$$$$$$$$$$$$2'
!!$       F_LOCALE(1)=-P*abs(Y1-Y3)*0.5;F_LOCALE(5)=F_LOCALE(1)
!!$    ELSE IF(X2==0.0d0 .AND. X3==0.0d0)THEN
!!$       !PRINT*,'$$$$$$$$$$$$$$$3'
!!$       F_LOCALE(3)=-P*abs(Y2-Y3)*0.5;F_LOCALE(5)=F_LOCALE(3)

    IF(ABS(X1-largeur)<=0.0000001 .AND. ABS(X2-largeur)<=0.0000001)THEN!(X1==largeur .AND. X2==largeur)THEN
       !PRINT*,'$$$$$$$$$$*******************$$$$$$$$$$$$$$$$$$$$$$$4',abs(Y1-Y2)
       F_LOCALE(1)=P*abs(Y1-Y2)*0.5;F_LOCALE(3)=F_LOCALE(1)
    ELSE IF(ABS(X1-largeur)<=0.0000001 .AND. ABS(X3-largeur)<=0.0000001)THEN!(X1==largeur .AND. X3==largeur)THEN
       !PRINT*,'*********5'
       F_LOCALE(1)=P*abs(Y1-Y3)*0.5;F_LOCALE(5)=F_LOCALE(1)
    ELSE IF(ABS(X2-largeur)<=0.0000001 .AND. ABS(X3-largeur)<=0.0000001)THEN!(X2==largeur .AND. X3==largeur)THEN
       !PRINT*,'********6'
       F_LOCALE(3)=P*abs(Y2-Y3)*0.5;F_LOCALE(5)=F_LOCALE(3)
    END IF

  end subroutine FORCE_LOCALE_SUIVANT_X









  subroutine FORCE_LOCALE_SUIVANT_Y_EXPANSEE(n,NDDL,Nd,localisation,F_LOCALE,F_EXPANSEE)
    integer,intent(in)::n,NDDL,Nd
    integer,dimension(1:n*NDDL),intent(in)::localisation
    real(rki),dimension(1:n*NDDL),intent(in)::F_LOCALE
    real(rki),dimension(1:Nd*NDDL),intent(out)::F_EXPANSEE
    integer::k,rank
    
    F_EXPANSEE=0.d0
    DO k=1,n*NDDL
       F_EXPANSEE(localisation(k))=F_LOCALE(k)
    END DO
    
  end subroutine FORCE_LOCALE_SUIVANT_Y_EXPANSEE

  subroutine node_loc(n,NDDL,nbr_tri,ite,rang,localisation,T_conn)
    integer,intent(in)::n,NDDL,nbr_tri,ite,rang
    integer,dimension(1:n*NDDL),intent(in)::localisation
    integer,dimension(1:nbr_tri,1:5),intent(in)::T_conn
    integer::k,j
    character(len=3)::rank
    write(rank,fmt='(i3.3)')rang
    OPEN(10+rang,file='NODE_LOC'//trim(adjustl(rank))//'.dat',position='append')
    DO j=2,4
       WRITE(10+rang,*)T_conn(ite,1),T_conn(ite,j),localisation(1:6)
       
    END DO
    CLOSE(10+rang)
  end subroutine node_loc



  SUBROUTINE WR_K(userchoice,Nd,NDDL,NNZ)
    integer,intent(inout)::userchoice
    integer,intent(in)::Nd,NDDL,NNZ
    !character(len=21),parameter::R_F='CLV_PC001_RANK000.dat'
    character(len=10),parameter::R_F='CCS000.dat'
    character(len=10),parameter::W_R='K_VIEW.vtk'
    integer::i,j,a,b,d,err
    real(rki)::c
    real(rki),dimension(:,:),allocatable::K
    
    K=0.d0
    SELECT CASE(userchoice)
    CASE(1)
       err=0
       ALLOCATE(K(1:Nd*NDDL,1:ND*NDDL))
       !allocate(K(1:Nd*NDDL))
       K=0.d0
       OPEN(30,file=W_R)
       write(30,fmt='(1A26)') '# vtk DataFile Version 3.0'
       write(30,fmt='(1A4)') 'cell'
       write(30,fmt='(1A5)') 'ASCII'
       write(30,fmt='(1A25)') 'DATASET STRUCTURED_POINTS'
       write(30,fmt='(1A11,1I5,1A1,1I5,1A2)') 'DIMENSIONS ', Nd*NDDL ,' ' , Nd*NDDL, ' 1'
       write(30,fmt='(1A12)') 'ORIGIN 0 0 0'
       write(30,fmt='(1A13)') 'SPACING 1 1 1'
       write(30,fmt='(1A11,1I9)') 'POINT_DATA ' , (Nd*NDDL)**2
       write(30,fmt='(1A18)') 'SCALARS cell float'
       write(30,fmt='(1A20)') 'LOOKUP_TABLE default'

       i=1

       OPEN(40,file=R_F)
       i=1
       DO WHILE(i<NNZ+1)
          READ(40,*)a,b,c,d
          PRINT*,'i=',i,' ',a,b,c,d
          K(b,a)=c
          i=i+1
          PRINT*,'i=',i
       END DO
       DO i=1,Nd*NDDL
          WRITE(30,*)K(i,:)
       END DO


       CLOSE(40)
       CLOSE(30)

       deallocate(K)
       userchoice=0
       if(err<0)then
          print*,'fin de fichier '
       else
          print*,'erreur de lecture'
       end if
    CASE DEFAULT
       userchoice=0
    END SELECT
  END SUBROUTINE WR_K





  subroutine SAVE_K_me(Nd,NDDL,K_GLOBALE,rang,num_pc,stateinfo,status)
    integer,dimension(MPI_STATUS_SIZE),intent(inout)::status
    integer,intent(inout)::stateinfo
    integer,intent(in)::Nd,NDDL,rang,num_pc
    real(rki),dimension(1:Nd*NDDL,1:Nd*NDDL),intent(in)::K_GLOBALE
    character(len=3)::rank,pc
    integer::i,j
    write(rank,fmt='(i3.3)')rang
    write(pc,fmt='(i3.3)')num_pc
    open(20+rang,file='PC_'//trim(adjustl(pc))//'_pour_'//'K_'//trim(adjustl(rank))//'.dat')

    DO j=1,Nd*NDDL
       DO i=1,ND*NDDL
          write(20+rang,*)K_GLOBALE(i,j)
       END DO
    END DO
    close(20+rang)
  end subroutine SAVE_K_ME



  subroutine VIEW_K_GLOBALE(Nd,NDDL,K_GLOBALE)

    integer,intent(in)::Nd,NDDL
    real(rki),dimension(1:Nd*NDDL,1:Nd*NDDL),intent(in)::K_GLOBALE
    real(rki),dimension(1:Nd*NDDL)::A
    integer::i,j
    open(20,file='view_K_globale.vtk')
    write(20,*)'# vtk DataFile Version 3.0'
    write(20,*)'cell'
    write(20,*)'ASCII'
    write(20,*)'DATASET STRUCTURED_POINTS'
    write(20,*)'DIMENSIONS ',Nd*NDDL,' ',Nd*NDDL,' 1'
    write(20,*)'ORIGIN 0 0 0 '
    write(20,*)'SPACING 1 1 1'
    write(20,*)'POINT_DATA ',(Nd*NDDL)**2
    write(20,*)'SCALARS cell float'
    write(20,*)'LOOKUP_TABLE default'
    DO j=1,Nd*NDDL
       A(:)=K_GLOBALE(:,j)
       !DO i=1,ND*NDDL
       !IF(K_GLOBALE(i,j)==0.d0)THEN
       write(20,*)A(:)
       !END DO
    END DO
    close(20)
  end subroutine VIEW_K_GLOBALE


end module ModMatrice






