module ModPGC
  use mpi
  use ModMesh
  implicit none

contains

  subroutine PGC(rang,it1,itN,NNZ_P,kmax,TOL,stateinfo,status,X,VAL,J_IND,I_cpt,b,Nd,NDDL,RES,H_NM)
    integer,intent(in)::rang,it1,itN,NNZ_P,Nd,NDDL,kmax
    real(rki),intent(in)::TOL
    CHARACTER(LEN=*),intent(in)::H_NM
    !NNZ_P := nbr de non zero par processeur
    integer,intent(inout)::stateinfo
    integer,dimension(MPI_STATUS_SIZE),intent(inout)::status
    real(rki),dimension(it1:itN),intent(inout)::X
    real(rki),dimension(it1:itN),intent(out)::RES
    real(rki),dimension(1:NNZ_P),intent(in)::VAL !ele non nul
    integer,dimension(1:NNZ_P),intent(in)::J_ind !num des cols ou ele non nul
    integer,dimension(1:itN-it1+2),intent(in)::I_cpt
    real(rki),dimension(it1:itN),intent(in)::b
    real(rki),dimension(it1:itN)::r,z
    real(rki),dimension(1:Nd*NDDL)::p,Ap
    real(rki)::beta,alpha,gamma
    real(rki)::beta_P,alpha_p
    real(rki)::drz,dApp,drz_old !<r|z>, <Ap|p>
    real(rki)::drz_P,dApp_P
    integer,dimension(1:2,1:4)::CHARGE
    integer::k
    IF(rang==0)THEN
       OPEN(20,file='./RES/'//trim(H_NM)//'res.txt')
    END IF
    !INI:
    beta=0;beta_P=0;k=0
    CALL COMM_CHARGE(rang,stateinfo,status,it1,itN,CHARGE)
    PRINT*,'FIN COMM_CHARGE'
    !INITIALISATION DE X à Xo=0.0d0:
    X=0.0d0
    !ro=b-A*X et beta =norm2(ro):
    r(it1:itN)=b(it1:itN)
    beta_P=DOT_PRODUCT(r(it1:itN),r(it1:itN))
    PRINt*,'INIT beta_P'
    CALL REDUCTION_REEL(rang,stateinfo,status,beta,beta_P)
    beta=sqrt(beta)
    PRINT*,'REDUCTION INITIALE DE beta'
    !zo=INV_DIAG(A).ro: (precondi jacobi)
    CALL PRE_DIAG(rang,it1,itN,NNZ_P,Nd,NDDL,VAL,J_IND,I_cpt,z(it1:itN),r(it1:itN))
    !INI p à z:
    p(it1:itN)=z(it1:itN)
    PRINT*,'INIT zo et po'
    IF(rang==0)THEN
       WRITE(20,*)k,beta
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)
    PRINT*,'DEBUT BOUCLE PGC'
    DO WHILE((k<kmax+1) .AND. (beta>TOL))
       k=k+1
       !CALCUL DE Ap_k: Ap=A.p
                  !SEND&RECV POUR QUE TOUS LES PROCS connaissent P(1:Nd*NDDL):
       CALL SEND_RECV_VECT(rang,CHARGE,Nd,NDDL,stateinfo,status,p)
       CALL SP_PMV(rang,it1,itN,NNZ_P,Nd,NDDL,VAL,J_IND,I_cpt,p,Ap)
       !CALCUL de <Ap|p>:
       dApp_P=DOT_PRODUCT(Ap(it1:itN),p(it1:itN))
       CALL REDUCTION_REEL(rang,stateinfo,status,dApp,dApp_P)
       !CALCUL de <r|z>:
       drz_P=DOT_PRODUCT(r(it1:itN),z(it1:itN))
       CALL REDUCTION_REEL(rang,stateinfo,status,drz,drz_P)
       !on doit stocker drz dans dzr_old:
       drz_old=drz
       !CALCUL de alpha_k: alpha=<r|z>/<Ap|p>:
       alpha=drz/dApp
       !CALCUL DE X_k+1:=X_k+alpha_k.p_k:
       X(it1:itN)=X(it1:itN)+alpha*p(it1:itN)
       !CALCUL DE r_k+1=r_k-alpha_k*Ap_k:
       r(it1:itN)=r(it1:itN)-alpha*Ap(it1:itN)
       !CALCUL DE z_k+1=INV_DIAG(A)*r:
       CALL PRE_DIAG(rang,it1,itN,NNZ_P,Nd,NDDL,VAL,J_IND,I_cpt,z(it1:itN),r(it1:itN))
       !RECALCUL DE <r|z>:
       drz_P=DOT_PRODUCT(r(it1:itN),z(it1:itN))
       CALL REDUCTION_REEL(rang,stateinfo,status,drz,drz_P)
       !CALCUL DE GAMMA:=drz/drz_old
       gamma=drz/drz_old
       !CALCUL DE P_K+1=Z_K+GAMMA*P_K:
       p(it1:itN)=z(it1:itN)+gamma*p(it1:itN)
       !CALCUL DE LA NORME DU RESIDU BETA:
       beta_P=DOT_PRODUCT(r(it1:itN),r(it1:itN))
       CALL REDUCTION_REEL(rang,stateinfo,status,beta,beta_P)
       beta=sqrt(beta)
       IF(rang==0)THEN
          WRITE(20,*)k,beta
          PRINT*,'k=',k,' beta=',beta,'TOL',TOL
       END IF
    END DO
    IF(rang==0)THEN
       CLOSE(20)
    END IF
    RES=r
  end subroutine PGC


  subroutine COMM_CHARGE(rang,stateinfo,status,it1,itN,CHARGE)
    integer,intent(in)::rang,it1,itN
    integer,intent(inout)::stateinfo
    integer,dimension(MPI_STATUS_SIZE),intent(inout)::status
    integer,dimension(1:2,1:4),intent(out)::CHARGE
    integer::i,he
    CHARGE(1,rang+1)=it1;CHARGE(2,rang+1)=itN

    DO he=0,3
       IF(rang/=he)THEN
          CALL MPI_SEND(CHARGE(1:2,rang+1),2,MPI_INTEGER,he,10+rang,MPI_COMM_WORLD,stateinfo)
       ELSE
          DO i=0,3
             IF(i/=rang)THEN
                CALL MPI_RECV(charge(1:2,i+1),2,MPI_INTEGER,i,10+i,MPI_COMM_WORLD,status,stateinfo)
             END IF
          END DO
       END IF
    END DO
  end subroutine COMM_CHARGE
  
  subroutine SEND_RECV_VECT(rang,CHARGE,Nd,NDDL,stateinfo,status,V_P)
    integer,intent(in)::rang,Nd,NDDL
    integer,dimension(1:2,1:4),intent(in)::CHARGE
    integer,intent(inout)::stateinfo
    integer,dimension(MPI_STATUS_SIZE),intent(inout)::status
    real(rki),dimension(1:Nd*NDDL),intent(inout)::V_P
    !real(rki),dimension(1:Nd*NDDL),intent(out)::V
    integer::i,he,G,D
    
     DO he=0,3
        IF(rang/=he)THEN
           G=CHARGE(1,rang+1);D=CHARGE(2,rang+1)
          CALL MPI_SEND(V_P(G:D),D-G+1,MPI_DOUBLE_PRECISION,he,10+rang,MPI_COMM_WORLD,stateinfo)
       ELSE
          DO i=0,3
             IF(i/=rang)THEN
                G=CHARGE(1,i+1);D=CHARGE(2,i+1)
                CALL MPI_RECV(V_P(G:D),D-G+1,MPI_DOUBLE_PRECISION,i,10+i,MPI_COMM_WORLD,status,stateinfo)
             END IF
          END DO
       END IF
    END DO
  end subroutine SEND_RECV_VECT
    
  subroutine REDUCTION_REEL(rang,stateinfo,status,R,R_P)
    integer,intent(in)::rang
    integer,intent(inout)::stateinfo
    integer,dimension(MPI_STATUS_SIZE),intent(inout)::status
    real(rki),intent(in)::R_P
    real(rki),intent(inout)::R
    real(rki)::reduc
    integer::he
    !*** SEND R_P:
    IF(rang/=0)THEN
       CALL MPI_SEND(R_P,1,MPI_DOUBLE_PRECISION,0,10+rang,MPI_COMM_WORLD,stateinfo)
    ELSE
       R=R_P
       do he=1,3
          CALL MPI_RECV(reduc,1,MPI_DOUBLE_PRECISION,he,10+he,MPI_COMM_WORLD,status,stateinfo)
          R=R+reduc
       end do
    END IF
    !*** SEND R:
    IF(rang==0)THEN
       do he=1,3
          CALL MPI_SEND(R,1,MPI_DOUBLE_PRECISION,he,10+he,MPI_COMM_WORLD,stateinfo)
       end do
    ELSE
       CALL MPI_RECV(R,1,MPI_DOUBLE_PRECISION,0,10+rang,MPI_COMM_WORLD,status,stateinfo)
    END IF
    
  end subroutine REDUCTION_REEL

  subroutine PRE_DIAG(rang,it1,itN,NNZ_P,Nd,NDDL,VAL,J_IND,I_cpt,z,r)
    integer,intent(in)::rang,it1,itN,NNZ_P,Nd,NDDL
    real(rki),dimension(1:NNZ_P),intent(in)::VAL
    integer,dimension(1:NNZ_P),intent(in)::J_IND
    integer,dimension(1:itN-it1+2),intent(in)::I_cpt
    real(rki),dimension(it1:itN),intent(out)::z !z=INV_DIAG(A).r
    real(rki),dimension(it1:itN),intent(in)::r !le residu
    !INV_DIAG := l'inverse de DIAG(A) pour le preconditionneur, avec A.X=b
    integer::i,k,a,C,d,b,NNZ_R
    !NNZ_P := nbr de non zero par processeur
    !NNZ_R := nbr de non zero par ligne par proc
    d=0;a=1;b=0;C=0
    DO i=it1,itN
       NNZ_R=I_cpt(a+1)-I_cpt(a)
       IF(NNZ_R/=0)THEN
          DO k=1,NNZ_R
             b=k+d
             C=J_IND(b)
             if(C==i)then
                z(i)=r(C)/VAL(b)
             end if
          END DO
          d=d+NNZ_R
       END IF
       a=a+1
    END DO
    
  end subroutine PRE_DIAG

  !PRODUIT MATRICE VECTEUR:
  subroutine SP_PMV(rang,it1,itN,NNZ_P,Nd,NDDL,VAL,J_IND,I_cpt,VECT,PMV)
    integer,intent(in)::rang,it1,itN,NNZ_P,Nd,NDDL
    real(rki),dimension(1:NNZ_P),intent(in)::VAL
    integer,dimension(1:NNZ_P),intent(in)::J_IND
    integer,dimension(1:itN-it1+2),intent(in)::I_cpt
    real(rki),dimension(1:Nd*NDDL),intent(in)::VECT
    real(rki),dimension(1:Nd*NDDL),intent(out)::PMV
    integer::i,j,k,a,b,d,C,NNZ_R
    d=0;a=1;b=0;C=0;PMV=0.0d0
    DO i=it1,itN
       NNZ_R=I_cpt(a+1)-I_cpt(a)
       IF(NNZ_R/=0)THEN
          do k=1,NNZ_R
             b=d+k
             C=J_IND(b)
             PMV(i)=PMV(i)+VAL(b)*VECT(C)
          end do
          d=d+NNZ_R
       END IF
       a=a+1
    END DO
  end subroutine SP_PMV

  !FORMAT SPARSE
  subroutine CCS(S_SIZE,Nd,NDDL,CCS_it1,CCS_itN,S_CL,S_V,VAL,J_IND,I_cpt)
    integer,intent(in)::S_SIZE,Nd,NDDL,CCS_it1,CCS_itN
    integer,dimension(1:S_SIZE,1:2),intent(in)::S_CL
    real(rki),dimension(1:S_SIZE),intent(in)::S_V
    real(rki),dimension(1:S_SIZE),intent(out)::VAL
    integer,dimension(1:S_SIZE),intent(out)::J_IND
    integer,dimension(1:(CCS_itN-CCS_it1)+2),intent(out)::I_cpt
    !integer,dimension(1:(Nd*NDDL)+1),intent(out)::I_cpt
    integer::i,j,k,NNZ_R,cpt
    I_cpt(1)=1;cpt=0;k=0
    DO i=CCS_it1,CCS_itN
       NNZ_R=0;k=k+1
       DO j=1,S_SIZE
          IF(S_CL(j,1)==i)THEN
             NNZ_R=NNZ_R+1;cpt=cpt+1
             VAL(cpt)=S_V(j);J_IND(cpt)=S_CL(j,2)
          END IF
       END DO
       I_cpt(k+1)=NNZ_R+I_cpt(k)
    END DO
  end subroutine CCS
  

  subroutine SEND_RECV_F_GLO(rang,CHARGE,it1,itN,Nd,NDDL,stateinfo,status,V_P,F_GP)
    integer,intent(in)::rang,Nd,NDDL,it1,itN
    integer,dimension(1:2,1:4),intent(in)::CHARGE
    integer,intent(inout)::stateinfo
    integer,dimension(MPI_STATUS_SIZE),intent(inout)::status
    real(rki),dimension(it1:itN),intent(inout)::V_P
    real(rki),dimension(1:Nd*NDDL),intent(inout)::F_GP !F_GLOBAL PAR PROC
    real(rki),dimension(:),allocatable::V_reduc
    integer::i,he,G,D
    ALLOCATE(V_reduc(it1:itN))
    V_P(it1:itN)=F_GP(it1:itN)
    DO he=0,3
       IF(rang/=he)THEN
          G=CHARGE(1,he+1);D=CHARGE(2,he+1)
          CALL MPI_SEND(F_GP(G:D),D-G+1,MPI_DOUBLE_PRECISION,he,10+rang,MPI_COMM_WORLD,stateinfo)
       ELSE
          DO i=0,3
             IF(i/=rang)THEN
                !G=CHARGE(1,rang+1);D=CHARGE(2,rang+1)
                !CALL MPI_RECV(V_reduc(G:D),D-G+1,MPI_DOUBLE_PRECISION,i,10+i,MPI_COMM_WORLD,status,stateinfo)
                CALL MPI_RECV(V_reduc(it1:itN),itN-it1+1,MPI_DOUBLE_PRECISION,i,10+i,MPI_COMM_WORLD,status,stateinfo)
                V_P(it1:itN)=V_P(it1:itN)+V_reduc(it1:itN)
             END IF
          END DO
       END IF
    END DO
    DEALLOCATE(V_reduc)
  end subroutine SEND_RECV_F_GLO

  subroutine SEND_RECV_SOL(rang,CHARGE,CCS_it1,CCS_itN,Nd,NDDL,stateinfo,status,X,SOL_X,SOL_Y)
    integer,intent(in)::rang,Nd,NDDL,CCS_it1,CCS_itN
    integer,dimension(1:2,1:4),intent(in)::CHARGE
    integer,intent(inout)::stateinfo
    integer,dimension(MPI_STATUS_SIZE),intent(inout)::status
    real(rki),dimension(CCS_it1:CCS_itN),intent(in)::X
    real(rki),dimension(1:Nd),intent(inout)::SOL_X,SOL_Y !DPL_X DPL_Y PAR PROC
    real(rki),dimension(:),allocatable::SOL_Xreduc1,SOL_Yreduc1,SOL_Xreduc2,SOL_Yreduc2
    integer::i,he,G,D,k,SIZE_X,a,b
    a=0;b=0
    ALLOCATE(SOL_Xreduc1(1:Nd),SOL_Yreduc1(1:Nd),SOL_Xreduc2(1:Nd),SOL_Yreduc2(1:Nd))
    SIZE_X=CCS_itN-CCS_it1+1
    SOL_X=0.0d0;SOL_Y=0.0d0
    SOL_Xreduc1=0.0d0;SOL_Yreduc1=0.0d0;SOL_Xreduc2=0.0d0;SOL_Yreduc2=0.0d0
    DO k=CCS_it1,CCS_itN
       IF(mod(k,2)==0)THEN
          SOL_Y(k/2)=X(k)
       ELSE
          a=(k-mod(k,2))/2+mod(k,2)
          SOL_X(a)=X(k)
       END IF
    END DO
    SOL_Xreduc2=SOL_X
    SOL_Yreduc2=SOL_Y
    DO he=0,3
       IF(rang/=he)THEN
          !G=CHARGE(1,rang+1);D=CHARGE(2,rang+1)
          CALL MPI_SEND(SOL_Xreduc2(1:Nd),Nd,MPI_DOUBLE_PRECISION,he,10+rang,MPI_COMM_WORLD,stateinfo)
          CALL MPI_SEND(SOL_Yreduc2(1:Nd),Nd,MPI_DOUBLE_PRECISION,he,10+rang,MPI_COMM_WORLD,stateinfo)
       ELSE
          DO i=0,3
             IF(i/=rang)THEN
                !G=CHARGE(1,i+1);D=CHARGE(2,i+1)
                CALL MPI_RECV(SOL_Xreduc1(1:Nd),Nd,MPI_DOUBLE_PRECISION,i,10+i,MPI_COMM_WORLD,status,stateinfo)
                CALL MPI_RECV(SOL_Yreduc1(1:Nd),Nd,MPI_DOUBLE_PRECISION,i,10+i,MPI_COMM_WORLD,status,stateinfo)
                SOL_X=SOL_X+SOL_Xreduc1
                SOL_Y=SOL_Y+SOL_Yreduc1
             END IF
          END DO
       END IF
    END DO
    DEALLOCATE(SOL_Xreduc1,SOL_Yreduc1,SOL_Xreduc2,SOL_Yreduc2)
  end subroutine SEND_RECV_SOL


  SUBROUTINE FIRST_REDUCTION_SIGMA_OU_DEF(rang,Nd,COEF_POND,SIGMA_GLO,DUMMY1,DUMMY2,stateinfo,status)
!OPERATION DE REDUCTION FACON TRI PAIR-IMPAIR:
    integer,intent(inout)::stateinfo
    integer,dimension(MPI_STATUS_SIZE),intent(inout)::status
    INTEGER,intent(in)::rang,Nd
    INTEGER,dimension(1:Nd),intent(inout)::COEF_POND
    INTEGER,dimension(:),allocatable,intent(inout)::DUMMY2
    REAL(rki),dimension(1:Nd),intent(inout)::SIGMA_GLO
    REAL(rki),dimension(:),allocatable,intent(inout)::DUMMY1
     IF(rang==1)THEN
        CALL MPI_SEND(COEF_POND(1:Nd),Nd,MPI_INTEGER,0,30,MPI_COMM_WORLD,stateinfo)
        CALL MPI_SEND(SIGMA_GLO(1:Nd),Nd,MPI_DOUBLE_PRECISION,0,31,MPI_COMM_WORLD,stateinfo)
     ELSE IF(rang==0)THEN
        ALLOCATE(DUMMY1(1:Nd),DUMMY2(1:Nd))
        CALL MPI_RECV(DUMMY2(1:Nd),Nd,MPI_INTEGER,1,30,MPI_COMM_WORLD,status,stateinfo)
        CALL MPI_RECV(DUMMY1(1:Nd),Nd,MPI_DOUBLE_PRECISION,1,31,MPI_COMM_WORLD,status,stateinfo)
        COEF_POND=COEF_POND+DUMMY2
        SIGMA_GLO=SIGMA_GLO+DUMMY1
        DEALLOCATE(DUMMY1,DUMMY2)
     ELSE IF(rang==2)THEN
        CALL MPI_SEND(COEF_POND(1:Nd),Nd,MPI_INTEGER,3,40,MPI_COMM_WORLD,stateinfo)
        CALL MPI_SEND(SIGMA_GLO(1:Nd),Nd,MPI_DOUBLE_PRECISION,3,41,MPI_COMM_WORLD,stateinfo)
     ELSE IF(rang==3)THEN
        ALLOCATE(DUMMY1(1:Nd),DUMMY2(1:Nd))
        CALL MPI_RECV(DUMMY2(1:Nd),Nd,MPI_INTEGER,2,40,MPI_COMM_WORLD,status,stateinfo)
        CALL MPI_RECV(DUMMY1(1:Nd),Nd,MPI_DOUBLE_PRECISION,2,41,MPI_COMM_WORLD,status,stateinfo)
        COEF_POND=COEF_POND+DUMMY2
        SIGMA_GLO=SIGMA_GLO+DUMMY1
        DEALLOCATE(DUMMY1,DUMMY2)
     END IF
   END SUBROUTINE FIRST_REDUCTION_SIGMA_OU_DEF




end module ModPGC
