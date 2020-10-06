program main
  use ModMesh
  use ModMatrice
  use ModParra
  use ModTri
  use mpi
  use ModPGC
  implicit none

  !mpirun -n 4 --mca pml ob1 ./a.out

  integer::rang,Np,stateinfo,me
  integer,dimension(MPI_STATUS_SIZE)::status
  integer::num_pc,indice,userchoice
  !POUR LE MAILLAGE:
  INTEGER::PX,PY,P1,P2,P3
  INTEGER::MARK1,MARK2 ! cf ModMesh.f90 pour info
  !n:nbr de noeuds par élé, ici 3
  !NDDL:nbr de degre de liberté aux noeuds, ici 2
  !Nd:nbr de noeuds total
  !taille matrice rigidité locale: n.NDDL par n.NDDL
  !taille matrice globale et expansee: Nd.NDDL par Nd.NDDL

  !POUR DEFINIR LES NOMS DES FICHIERS
  character(len=20)::mesh_file
  character(len=3)::rank,N_PC,u_dir
  character(len=4)::INPUT_DATA
  character(len=25)::DEF_NODE
  character(len=25)::DEF_X,DEF_Y,RES_X,RES_Y
  !*** POUR MEF
  integer::nbr_node,nbr_tri,Nd
  INTEGER::it1,itn,CHARGE_TOT,CHARGE_ME,k_ite
  integer::i,j,k,l,he,compteur_ite,NNZ,NNZ1,NNZ2,NNZ0,NNZ3,TRUE_NNZ,NZ
  integer,parameter::n=3,NDDL=2
  integer::n_NDDL
  !T_node TABLEAU DE CORD DES NOEUDS          
  TYPE(tab_node)::T_node
  !ALLOCATE(T_node%n_node(1:nbr_node),T_node%c_node(1:nbr_node,1:3))
  integer,dimension(:,:),allocatable::T_conn

  !*** POUR MEF
  real(rki),dimension(1:2,1:2)::jac,inv_jac
  real(rki),dimension(1:3,1:4)::Q
  real(rki),dimension(1:4,1:6)::DN !DERIVE FCT DE FORME DANS BASE ELE REF
  real(rki),dimension(1:3,1:6)::B
  real(rki),dimension(1:3,1:3)::HOOK
  real(rki),dimension(1:6,1:6)::K_ele
  real(rki),dimension(1:6,1:3)::inter
  integer,dimension(1:n*NDDL)::localisation
  integer,dimension(1:n*NDDL,1:3)::liste

  !CONDITION DPL BLOQUE:
  real(rki)::condi_x,condi_y
  INTEGER::MODEL !POUR SAVOIR SI 1/4 DE PLAQUE OU ENTIERE
  !*** POUR TRI
  INTEGER::T1,T2,T3,Z !Z nbr ZERO
  INTEGER::TAILLE
  INTEGER,dimension(:),allocatable::RC1,RC2
  real(rki),dimension(:),allocatable::RC3

  integer,dimension(:,:),allocatable::TAB_CL,T_CL,T01_CL,T23_CL,L01_CL,L23_CL,T03_CL,L03_CL,S_CL
  real(rki),dimension(:),allocatable::TAB_V,T_V,T01_V,T23_V,L01_V,L23_V,T03_V,L03_V,S_V

  !*** POUR MEF
  real(rki)::det
  
  REAL(rki)::longueur,largeur,RAY 
  
  REAL(rki)::P,P_MAX
  INTEGER::N_P
  INTEGER::N_STEP,STEP
  real(rki),dimension(1:n*NDDL)::F_LOCALE
  real(rki),dimension(:),allocatable::F_LOCALE_EXPANSEE,F_GLOBAL,F
  

  !*** POUR PGC ET FORMAT CCS:
  real(rki)::TOL
  integer::CCS_it1,CCS_itN,AA,BB,CC,DD,KMAX
  real(rki),dimension(:),allocatable::VAL
  integer,dimension(:),allocatable::J_IND
  integer,dimension(:),allocatable::I_cpt
  real(rki),dimension(:),allocatable::X,SOLX,SOLY,RES,RESX,RESY
  integer,dimension(1:2,1:4)::TAB_IT,TAB_CHARGE

  !*** POUR CALCUL DE SIGMA ET DEFORMATION:
  INTEGER::DIR,DIR_ALL,TRACE
  REAL(rki),dimension(1:6)::DPL
  REAL(rki),dimension(1:3)::SIGMA_LOC,DEF_LOC
  REAL(rki),dimension(1:3,1:6)::AB
  REAL(rki),dimension(:),allocatable::SIGMA_GLO,DUMMY1,DEF_GLO,DUMMY3
  INTEGER,dimension(:),allocatable::COEF_POND,DUMMY2,COEF_POND_2,DUMMY4 !coeff de pondération
  !*** POUR LES MODULE E_XX,E_YYou G_XY:
  !REAL(rki),dimension(:),allocatable::MOD_ING
  REAL(rki)::MOY_MOD_ING,SOM_SIG,SOM_DEF,X_TRACE,Y_TRACE !moyenne et point de niveau pour relever S_YY sur la ligne
  INTEGER::CPT_TRACE,H_A,H_B
  CHARACTER(len=20)::HOST_NM
  !*** CPU TIME
  real(rki)::start,end
  
  !******************************************************
  !******************************************************
  CALL MPI_INIT(stateinfo) !DEBUT REGION //
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,stateinfo)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,Np,stateinfo)
  PRINT*,'AFFICHAGE NBR PROC:'
  IF(rang==0)THEN
     PRINT*,'**********************************'
     PRINT*,'NOMBRE DE PROCESSUS ACTIF:',Np
     PRINT*,'**********************************'
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)
  PRINT*,'|'
  PRINT*,'PROCESSUS #',rang
  PRINT*,'|'
  IF(rang==0)THEN
     PRINT*,'**********************************'
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)
  !******************************************************
  indice=0
  !INI
  n_NDDL=n*NDDL

  !condi_dpl:
  !condi_x=0.0125d0;condi_y=0.08d0
  !condi_x=0.036d0;condi_y=0.0d0!0.2305d0
  !INITIALAISATION DN:
  DN=0.0d0
  DN(1,1)=-1.0d0;DN(1,3)=1.0d0;DN(2,1)=-1.0d0;DN(2,5)=1.0d0;DN(3,2)=-1.0d0;DN(3,4)=1.0d0
  DN(4,2)=-1.0d0;DN(4,6)=1.0d0

  !INITIALISATION DE HOOK: on switch c11 avec c22 car traction sens y
  HOOK=0.d0

  !OTCH2:
  HOOK(1,1)=43.0630d0;HOOK(2,2)=HOOK(1,1)
  HOOK(1,2)=13.4340d0;HOOK(2,1)=HOOK(1,2)
  HOOK(3,3)=14.8140d0
  HOOK=HOOK*10**9
  HOOK(1,3)=0.0d0;HOOK(3,1)=0.0d0;HOOK(2,3)=0.0d0;HOOK(3,2)=0.0d0

 

  !mesh_file='HOLE.o.msh'
  write(rank,fmt='(i3.3)')rang
  !mesh_file='HOLE_'//trim(adjustl(rank))//'.o.msh'

  !*** LECTURE DE DATA.txt POUR LES PARAMETRES DU CODE:
  INPUT_DATA='DATA'!.txt'
  CALL READ_DATA(rang,INPUT_DATA,mesh_file,nbr_node,nbr_tri,MARK1,MARK2,KMAX,TOL,P_MAX,STEP,N_STEP &
       ,MODEL,longueur,largeur,RAY,condi_x,condi_y)!,DIR)

  SELECT CASE(MODEL)
     CASE(1)
        X_TRACE=0.0d0;Y_TRACE=0.0d0
     CASE(0)
        X_TRACE=largeur/2.0d0;Y_TRACE=longueur/2.0d0   
  END SELECT

  num_pc=1
  write(N_PC,fmt='(i3.3)')num_pc
  Nd=nbr_node
  CHARGE_TOT=nbr_tri
  
!BOUCLE POUR VARTIATION DE CHARGE ET EVENTUELLLEMENT SI L'ON VEUT DISTRIBUER LA BOUCLE SUR PLUSIEURS PC:
  CALL hostnm(HOST_NM)!LE NOM DU PC
IF(STEP==1)THEN

   IF(HOST_NM=='foissac')THEN
      CALL CHARGE(rang,5,H_A,H_B,N_STEP,0)!!H_A=1;H_B=10
   ELSE IF(HOST_NM=='enlene')THEN
      CALL CHARGE(rang,5,H_A,H_B,N_STEP,0)
   ELSE IF(HOST_NM=='combarelles')THEN
     CALL CHARGE(rang,5,H_A,H_B,N_STEP,0)! H_A=21;H_B=30
   ELSE IF(HOST_NM=='cosquer')THEN
      CALL CHARGE(rang,5,H_A,H_B,N_STEP,0)!H_A=31;H_B=40
   ELSE IF(HOST_NM=='cussac')THEN
      CALL CHARGE(rang,5,H_A,H_B,N_STEP,0)!H_A=41;H_B=50
   END IF
ELSE IF(STEP==0)THEN
   HOST_NM='DELL'
   H_A=1;H_B=1;N_STEP=1
END IF
P=0.0d0
  DO N_P=H_A,H_B
     P=P+N_P*(P_MAX/N_STEP)!N_P*POURCENTAGE_P*P_MAX

  ALLOCATE(T_conn(1:nbr_tri,1:5))

  ALLOCATE(F_LOCALE_EXPANSEE(1:Nd*NDDL),F_GLOBAL(1:Nd*NDDL))



     F_GLOBAL=0.d0
 




  CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)



  PRINT*,'READ MESH'
  call read_mesh(mesh_file,nbr_node,nbr_tri,T_node,T_conn,rang,MARK1,MARK2)



  CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)
  !*** DEFINITION DES CHARGES PAR PROCESSEURS
  CALL CHARGE(rang,Np,it1,itN,CHARGE_TOT,0)!CHARGE_TOT=4970,enplus)
  CHARGE_ME=itN-it1+1
  CALL COMM_CHARGE(rang,stateinfo,status,it1,itN,TAB_CHARGE)
  ALLOCATE(TAB_CL(1:36*CHARGE_ME,1:2),TAB_V(1:36*CHARGE_ME))
  

  CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)



  

  compteur_ite=1
  PRINT*,'DEBUT BOUCLE MEF'
  start=MPI_Wtime()
  DO i=it1,itN !1,nbr_tri
     PRINT*,'RANG',rang,'i=',i
     
     
     call JACOBIENNE(T_node,T_conn,jac,i,nbr_tri)
     
     call DETERMINANT_JAC(jac,det)
     
     call INVERSE_JAC(jac,det,inv_jac)
     
     call Q_INV_J(inv_jac,Q)
     
     call GRAD_N(Q,DN,B)
    
     call MATRICE_RIGIDITE_ELEMENT(B,HOOK,det,K_ele)
     !pour les dpl bloqué:
     IF(MODEL==1)THEN !MODEL 1/4 DE PLAQUE
        call DPL_IMPO(K_ele,T_conn,T_node,nbr_tri,i,condi_x,condi_y)
     END IF
     

     call LOC_i(n,NDDL,nbr_tri,T_conn,i,localisation)
     
     !**************************************************************************************
     CALL MAKE_A_LISTE1(rang,36,CHARGE_ME,i,it1,k_ite,n_NDDL,localisation,K_ele,TAB_CL,TAB_V) !POUR FORMAT CREUX
     !21 pour low tri part but use ful matrix i.e 6*6=36
     !*************************************************************************************

     call FORCE_LOCALE_SUIVANT_Y(i,P,longueur,n,Nd,NDDL,nbr_tri,T_node,T_conn,F_LOCALE,MODEL)
     
     !*********************************************************************************************
     call FORCE_LOCALE_SUIVANT_Y_EXPANSEE(n,NDDL,Nd,localisation,F_LOCALE,F_LOCALE_EXPANSEE)
     F_GLOBAL=F_GLOBAL+F_LOCALE_EXPANSEE
     !*********************************************************************************************
     

     compteur_ite=compteur_ite+1
  END DO !FIN BOUCLE i
  end=MPI_Wtime()
  OPEN(60+rang,file='./TPS_CPU/'//trim(HOST_NM)//'temps_cpu.txt',position='append')
  write(60+rang,*)'RANG=',rang,'MEF_1',end-start
  CLOSE(60+rang)
  DEALLOCATE(F_LOCALE_EXPANSEE)



  CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)
  PRINT*,'FIN BOUCLE'


  !*******************PRE_TRI*********************************
  !*** POUR TRI PAR ORDRE CROISSANT SUR LES COLS DE K_GLOBALE
  !*** EN VUE DU FORMAT CCS
  T1=((Nd*NDDL)-mod((Nd*NDDL),4))/4;T2=2*T1;T3=3*T1
  PRINT*,size(TAB_V),mod(size(TAB_V),4),T1,T2,T3,Nd*NDDL
  start=MPI_Wtime()
  CALL PRE_TRI(T1,T2,T3,Nd*NDDL,size(TAB_V),rang,TAB_CL(:,1),TAB_CL(:,2),TAB_V,RC1,RC2,RC3,status,stateinfo)
  end=MPI_Wtime()
  OPEN(60+rang,file='./TPS_CPU/'//trim(HOST_NM)//'temps_cpu.txt',position='append')
  write(60+rang,*)'RANG=',rang,'PRE_TRI',end-start
  CLOSE(60+rang)
  CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)
  DEALLOCATE(TAB_CL,TAB_V)
  IF(size(RC1)>0)THEN
     
     ALLOCATE(TAB_CL(1:size(RC1),1:2),TAB_V(1:size(RC1)))
     
     TAB_CL(:,1)=RC1;TAB_CL(:,2)=RC2;TAB_V=RC3
     PRINT*,TAB_CL(1:3,1)
  END IF
  DEALLOCATE(RC1,RC2,RC3)
  PRINT*,'FIN PRE TRI'
  !*********************TRI*************************************
  PRINT*,'DEBUT TRI:'
  start=MPI_Wtime()
  IF(SIZE(TAB_V)>0)THEN!(SIZE(RC1)>0)THEN
     
    CALL  DOUBLON(size(TAB_V),rang,TAB_CL,TAB_V,Z)
    
    CALL TRI_CR(size(TAB_V),Z,rang,TAB_CL,TAB_V,S_CL,S_V,NNZ)
    
  END IF
  end=MPI_Wtime()
  OPEN(60+rang,file='./TPS_CPU/'//trim(HOST_NM)//'temps_cpu.txt',position='append')
  write(60+rang,*)'RANG=',rang,'TRI',end-start
  CLOSE(60+rang)
  DEALLOCATE(TAB_CL,TAB_V)
   !*******************FIN TRI*********************************

  CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)
  PRINT*,'RANG',rang,' PREMIER TRI DONE'


  !*** ON REMPLIT LE FORMAT CCS (A est SYM DONC CCS=CCR)
  CCS_it1=minval(S_CL(:,1));CCS_itN=maxval(S_CL(:,1))
  ALLOCATE(VAL(1:size(S_V)),J_IND(1:size(S_V)),I_cpt(1:(CCS_itN-CCS_it1)+2))!I_cpt(1:(Nd*NDDL)+1))
  start=MPI_Wtime()
  CALL CCS(size(S_V),Nd,NDDL,CCS_it1,CCS_itN,S_CL,S_V,VAL,J_IND,I_cpt)
  end=MPI_Wtime()
  OPEN(60+rang,file='./TPS_CPU/'//trim(HOST_NM)//'temps_cpu.txt',position='append')
  write(60+rang,*)'RANG=',rang,'CCS',end-start,' NNZ=',size(VAL)
  CLOSE(60+rang)
  !***POUR VERIFICATION DECOMMENTER:

!!$  open(30+rang,file='SPARSE'//trim(adjustl(rank))//'.dat')
!!$  BB=CCS_it1;CC=0
!!$  DO i=1,size(I_cpt)-1!NNZ!21*CHARGE_ME
!!$     AA=I_cpt(i+1)-I_cpt(i)
!!$     IF(AA/=0)THEN
!!$        DO j=1,AA
!!$           write(30+rang,*)BB,J_IND(j+CC),VAL(j+CC)!,TRUE_NNZ!,NNZ!21*CHARGE_ME
!!$        END DO
!!$     END IF
!!$     BB=BB+1;CC=CC+AA
!!$  END DO
!!$  close(30+rang)


  
  DEALLOCATE (S_CL,S_V)
  CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)



  PRINT*,'rang=',rang,'it1 itN=',it1,itN
  !*** SWITCH DE LA CHARGE PAR PROCS POUR EN VUE DU SOLVEUR:
  it1=CCS_it1;itN=CCS_itN
  PRINT*,'rang=',rang,'CCS_it1 CCS_itN=',CCS_it1,CCS_itN
  
  !*** CELA S'APPARENTE A UNE OPERATION DE REDUCTION PAR PARTIE 
  !*** SUR F_GLOBAL POUR QUE TOUT LES PROCS AIENT L'ADDITION DE F_GLOBAL SUR LEUR it1:itN
  !*** STOCKE DANS F:
  ALLOCATE(F(it1:itN))
  CALL COMM_CHARGE(rang,stateinfo,status,it1,itN,TAB_IT)
  CALL SEND_RECV_F_GLO(rang,TAB_IT,it1,itN,Nd,NDDL,stateinfo,status,F,F_GLOBAL)
  CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)
  DEALLOCATE(F_GLOBAL)

  
  !*** ALLOC SOL INIT ET PGC:
  ALLOCATE(X(it1:itN),RES(it1:itN))
  !TOL=0.000001
  PRINT*,'DEBUT PGC PRECONDITIONEUR DIAGONAL'
  start=MPI_Wtime()
  CALL PGC(rang,it1,itN,size(VAL),KMAX,TOL,stateinfo,status,X,VAL,J_IND,I_cpt,F,Nd,NDDL,RES,HOST_NM)
  end=MPI_Wtime()
  OPEN(60+rang,file='./TPS_CPU/'//trim(HOST_NM)//'temps_cpu.txt',position='append')
  write(60+rang,*)'RANG=',rang,'PGC+WR_res',end-start,' NNZ=',size(VAL),'it1',it1,'itN',itN
  CLOSE(60+rang)
  


  ALLOCATE(SOLX(1:Nd),SOLY(1:Nd))
  ALLOCATE(RESX(1:Nd),RESY(1:Nd))
  CALL SEND_RECV_SOL(rang,TAB_IT,CCS_it1,CCS_itN,Nd,NDDL,stateinfo,status,X,SOLX,SOLY)
  CALL SEND_RECV_SOL(rang,TAB_IT,CCS_it1,CCS_itN,Nd,NDDL,stateinfo,status,RES,RESX,RESY)
  DEALLOCATE(X,F,RES)
  PRINT*,'rang=',rang,' SOL_X(1:5)',SOLX(1:5)
  !*** MODIFICATION DU FICHIER DE MAILLAGE POUR VOIR LA DEFORMEE:
  DEF_NODE='./DPL/'//trim(HOST_NM)//'NODE.txt'
  DEF_X='./DPL/'//trim(HOST_NM)//'DPL_X.txt'
  DEF_Y='./DPL/'//trim(HOST_NM)//'DPL_Y.txt'
  RES_X='./RES/'//trim(HOST_NM)//'RES_X.txt'
  RES_Y='./RES/'//trim(HOST_NM)//'RES_Y.txt'
  CALL DEFORMEE_DPL(rang,mesh_file,DEF_NODE,DEF_X,DEF_Y,Nd,SOLX,SOLY)
  CALL MAP_RESIDUS(rang,RES_X,RES_Y,Nd,RESX,RESY)
  !***
  DEALLOCATE(VAl,J_IND,I_cpt)
  
  DEALLOCATE(RESX,RESY)
  !*** CALCUL DES DEFORMATION ET DE SIGMA:
  DO DIR_ALL=1,3
     
     DIR=DIR_ALL
     write(u_dir,fmt='(i3.3)')DIR
     ALLOCATE(DEF_GLO(1:Nd),COEF_POND_2(1:Nd))
     DEF_GLO=0.d0;COEF_POND_2=0
     ALLOCATE(SIGMA_GLO(1:Nd),COEF_POND(1:Nd))
     SIGMA_GLO=0.d0;COEF_POND=0
     !DIR=2



     compteur_ite=1
     PRINT*,'DEBUT BOUCLE POUR CALCUL DE SIGMA'
     start=MPI_WTIME()
     DO i=TAB_CHARGE(1,rang+1),TAB_CHARGE(2,rang+1) 
        !PRINT*,'RANG',rang,'i=',i
        call JACOBIENNE(T_node,T_conn,jac,i,nbr_tri)

        call DETERMINANT_JAC(jac,det)


        call INVERSE_JAC(jac,det,inv_jac)


        call Q_INV_J(inv_jac,Q)


        call GRAD_N(Q,DN,B)


        
        !FORMATION DE DPL
        T1=T_conn(i,2);T2=T_conn(i,3);T3=T_conn(i,4)

        DPL(1)=SOLX(T1);DPL(2)=SOLY(T1)
        DPL(3)=SOLX(T2);DPL(4)=SOLY(T2)
        DPL(5)=SOLX(T3);DPL(6)=SOLY(T3)
        !PRINT*,'rang',rang,'T1',T1,'SOLX(T1)',SOLX(T1),'DPL(1)',DPL(1)
        
        !CALCUL DE DEFORMATION: DEF_LOC:
        DEF_LOC=MATMUL(B,DPL)
        !CALCUL DE SIGMA_LOC:
        SIGMA_LOC=MATMUL(HOOK,DEF_LOC)
        

        !REDUCTION/EXTRAPOLATION AUX NOEUDS DE DEFORMATION: SIGMA_GLO:
        DEF_GLO(T1)=DEF_GLO(T1)+DEF_LOC(DIR);COEF_POND_2(T1)=COEF_POND_2(T1)+1
        DEF_GLO(T2)=DEF_GLO(T2)+DEF_LOC(DIR);COEF_POND_2(T2)=COEF_POND_2(T2)+1
        DEF_GLO(T3)=DEF_GLO(T3)+DEF_LOC(DIR);COEF_POND_2(T3)=COEF_POND_2(T3)+1
        !REDUCTION/EXTRAPOLATION AUX NOEUDS DE SIGMA: SIGMA_GLO:
        SIGMA_GLO(T1)=SIGMA_GLO(T1)+SIGMA_LOC(DIR);COEF_POND(T1)=COEF_POND(T1)+1
        SIGMA_GLO(T2)=SIGMA_GLO(T2)+SIGMA_LOC(DIR);COEF_POND(T2)=COEF_POND(T2)+1
        SIGMA_GLO(T3)=SIGMA_GLO(T3)+SIGMA_LOC(DIR);COEF_POND(T3)=COEF_POND(T3)+1



        compteur_ite=compteur_ite+1
     END DO !FIN BOUCLE i

     CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)
     IF(DIR_ALL==3)THEN
        DEALLOCATE(SOLX,SOLY)
     END IF


     !OPERATION DE REDUCTION FACON TRI PAIR-IMPAIR:
     CALL FIRST_REDUCTION_SIGMA_OU_DEF(rang,Nd,COEF_POND,SIGMA_GLO,DUMMY1,DUMMY2,stateinfo,status)
     CALL FIRST_REDUCTION_SIGMA_OU_DEF(rang,Nd,COEF_POND_2,DEF_GLO,DUMMY3,DUMMY4,stateinfo,status)


     CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)

!LAST_REDUCTION POUR DEFORMATION:
     IF(rang==3)THEN
        CALL MPI_SEND(COEF_POND_2(1:Nd),Nd,MPI_INTEGER,0,50,MPI_COMM_WORLD,stateinfo)
        CALL MPI_SEND(DEF_GLO(1:Nd),Nd,MPI_DOUBLE_PRECISION,0,51,MPI_COMM_WORLD,stateinfo)
        DEALLOCATE(DEF_GLO,COEF_POND_2)
     ELSE IF(rang==0)THEN
        ALLOCATE(DUMMY3(1:Nd),DUMMY4(1:Nd))
        CALL MPI_RECV(DUMMY4(1:Nd),Nd,MPI_INTEGER,3,50,MPI_COMM_WORLD,status,stateinfo)
        CALL MPI_RECV(DUMMY3(1:Nd),Nd,MPI_DOUBLE_PRECISION,3,51,MPI_COMM_WORLD,status,stateinfo)
        COEF_POND_2=COEF_POND_2+DUMMY4
        DEF_GLO=DEF_GLO+DUMMY3
        DEALLOCATE(DUMMY3,DUMMY4)
        IF(DIR_ALL==1)THEN
           OPEN(60,file='./DEF/'//trim(HOST_NM)//'def_xx.txt')
        ELSE IF(DIR_ALL==2)THEN
           OPEN(60,file='./DEF/'//trim(HOST_NM)//'def_yy.txt')
        ELSE IF(DIR_ALL==3)THEN
           OPEN(60,file='./DEF/'//trim(HOST_NM)//'def_xy.txt')
        END IF
        DO i=1,Nd
           DEF_GLO(i)=DEF_GLO(i)/COEF_POND_2(i)
           WRITE(60,*)i,DEF_GLO(i)
        END DO
        CLOSE(60)
        IF(DIR_ALL==2)THEN
           open(90,file='./DEF/'//trim(HOST_NM)//'TRACE_DEF_YY.txt')
           CPT_TRACE=0
           SOM_DEF=0.0d0
           DO i=1,Nd
              TRACE=T_node%n_node(i)
              !PRINT*,'TRACE=',TRACE,' T_node%c_node(TRACE)=',T_node%c_node(TRACE,2)
              IF(T_node%c_node(TRACE,2)==Y_TRACE .AND. T_node%c_node(TRACE,1)>=X_TRACE )THEN
                 write(90,*)T_node%c_node(TRACE,1)-RAY,DEF_GLO(TRACE)
                 !PRINT*,T_node%c_node(TRACE,1),SIGMA_GLO(TRACE)
              END IF
              SOM_DEF=SOM_DEF+DEF_GLO(TRACE)
              CPT_TRACE=CPT_TRACE+1
           END DO
           SOM_DEF=SOM_DEF/CPT_TRACE
           close(90)
        END IF
        DEALLOCATE(COEF_POND_2)!DEALLOCATE(DEF_GLO,COEF_POND_2)
     ELSE
        DEALLOCATE(DEF_GLO,COEF_POND_2)
     END IF



!LAST_REDUCTION POUR SIGMA:
     IF(rang==3)THEN
        CALL MPI_SEND(COEF_POND(1:Nd),Nd,MPI_INTEGER,0,50,MPI_COMM_WORLD,stateinfo)
        CALL MPI_SEND(SIGMA_GLO(1:Nd),Nd,MPI_DOUBLE_PRECISION,0,51,MPI_COMM_WORLD,stateinfo)
        DEALLOCATE(SIGMA_GLO,COEF_POND)
     ELSE IF(rang==0)THEN
        ALLOCATE(DUMMY1(1:Nd),DUMMY2(1:Nd))
        CALL MPI_RECV(DUMMY2(1:Nd),Nd,MPI_INTEGER,3,50,MPI_COMM_WORLD,status,stateinfo)
        CALL MPI_RECV(DUMMY1(1:Nd),Nd,MPI_DOUBLE_PRECISION,3,51,MPI_COMM_WORLD,status,stateinfo)
        COEF_POND=COEF_POND+DUMMY2
        SIGMA_GLO=SIGMA_GLO+DUMMY1
        DEALLOCATE(DUMMY1,DUMMY2)
        IF(DIR_ALL==1)THEN
           OPEN(60,file='./SIGMA/'//trim(HOST_NM)//'sigma_xx.txt')
        ELSE IF(DIR_ALL==2)THEN
           OPEN(60,file='./SIGMA/'//trim(HOST_NM)//'sigma_yy.txt')
        ELSE IF(DIR_ALL==3)THEN
           OPEN(60,file='./SIGMA/'//trim(HOST_NM)//'sigma_xy.txt')
        END IF
        DO i=1,Nd
           SIGMA_GLO(i)=SIGMA_GLO(i)/COEF_POND(i)
           WRITE(60,*)i,SIGMA_GLO(i)
        END DO
        CLOSE(60)
        IF(DIR_ALL==2)THEN
           open(90,file='./SIGMA/'//trim(HOST_NM)//'TRACE_SIGMA_YY.txt')
           open(95,file='./MODULE_EVO/'//trim(HOST_NM)//'ELASTICITE_LIN.txt',position='append')
           SOM_SIG=0.d0
           CPT_TRACE=0
           DO i=1,Nd
              TRACE=T_node%n_node(i)
              !PRINT*,'TRACE=',TRACE,' T_node%c_node(TRACE)=',T_node%c_node(TRACE,2)
              IF(T_node%c_node(TRACE,2)==Y_TRACE .AND. T_node%c_node(TRACE,1)>=X_TRACE )THEN
                 write(90,*)T_node%c_node(TRACE,1)-RAY,SIGMA_GLO(TRACE)
                 !PRINT*,T_node%c_node(TRACE,1),SIGMA_GLO(TRACE)
              END IF
              SOM_SIG=SOM_SIG+SIGMA_GLO(TRACE)
              CPT_TRACE=CPT_TRACE+1
           END DO
           close(90)
           SOM_SIG=SOM_SIG/CPT_TRACE
           WRITE(95,*)SOM_DEF,SOM_SIG,P
           close(95)
           
        END IF
        DEALLOCATE(COEF_POND)!DEALLOCATE(SIGMA_GLO,COEF_POND)
     ELSE
        DEALLOCATE(SIGMA_GLO,COEF_POND)
     END IF



     IF(rang==0)THEN
        !allocate(MOD_ING(1:Nd))
        MOY_MOD_ING=0.0d0
        DO i=1,Nd
           MOY_MOD_ING=MOY_MOD_ING+SIGMA_GLO(i)/DEF_GLO(i)
           !MOD_ING(i)=SIGMA_GLO(i)/DEF_GLO(i)
           !MOY_MOD_ING=MOY_MOD_ING+MOD_ING(i)
        END DO
        MOY_MOD_ING=MOY_MOD_ING/Nd
        OPEN(10+DIR_ALL,file='./MODULE_EVO/'//trim(HOST_NM)//'EVO_MOD_ING_4AVERAGE'//trim(adjustl(u_dir))//'.txt',POSITION='APPEND')
        WRITE(10+DIR_ALL,*)P,MOY_MOD_ING
        CLOSE(10+DIR_ALL)
        DEALLOCATE(SIGMA_GLO,DEF_GLO)
     END IF



  end=MPI_WTIME()
  OPEN(70+rang,file='./TPS_CPU/'//trim(HOST_NM)//'temps_cpu.txt',position='append')
  !write(70+rang,*)'RANG=',rang,'MEF2+SIGMA',end-start
  IF(DIR_ALL==1)THEN
           write(70+rang,*)'RANG=',rang,'MEF2+SIGMA_XX',end-start
        ELSE IF(DIR_ALL==2)THEN
           write(70+rang,*)'RANG=',rang,'MEF2+SIGMA_YY',end-start
        ELSE IF(DIR_ALL==3)THEN
           write(70+rang,*)'RANG=',rang,'MEF2+SIGMA_XY',end-start
        END IF
  CLOSE(70+rang)


END DO
  CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)




  DEALLOCATE(T_node%n_node,T_node%c_node)
  DEALLOCATE(T_conn)
  PRINT*,'DEALLOCATE T_node et T_conn'

  CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)
END DO


  CALL MPI_BARRIER(MPI_COMM_WORLD,stateinfo)
  CALL MPI_FINALIZE(stateinfo)

end program main



