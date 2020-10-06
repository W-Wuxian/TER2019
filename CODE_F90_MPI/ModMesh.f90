Module ModMesh

  implicit none
  integer,parameter::rki=8
  
  TYPE tab_node
     integer,dimension(:),allocatable::n_node
     real(rki),dimension(:,:),allocatable::c_node !1::x, 2::y, 3::z
  END TYPE tab_node

!LES 2 TYPES SUIVANT SONT POUR STOCKER LES NODES QUI SONT SUR LES LIGNES EN SUD ET NORD
!C'EST POUR LES FORCES
  TYPE info_pts
     integer::num_node
     real(rki)::x_node,y_node
  END TYPE info_pts

  TYPE TAB_INFO
     TYPE(info_pts),dimension(1:6)::liste
  END TYPE TAB_INFO
contains
  

  subroutine read_mesh(mesh_file,nbr_node,nbr_tri,T_node,T_conn,rang,MARK1,MARK2)
    character(len=*),intent(in)::mesh_file
    character(len=25)::CHEMIN
    integer,intent(in)::nbr_node,nbr_tri,rang,MARK1,MARK2
    !MARK1 1ere ligne où les éléments sont triangles, MARK2  dernière ligne où les élé sont triangles
    integer::k,j,b,c,n1,n2,n3,n4,n5,n6,n7,n8,T
    real::a                          !a,b,c dummy variables
    real(rki)::x,y,z                      !x,y,z coordonnées nodes
    integer::num_node                !numero du node
 
    TYPE(tab_node),intent(out)::T_node
    integer,dimension(1:nbr_tri,1:5),intent(out)::T_conn
    ALLOCATE(T_node%n_node(1:nbr_node),T_node%c_node(1:nbr_node,1:3))
    CHEMIN='./MESH/'//mesh_file
    open(10,file=CHEMIN)!mesh_file)
    k=1
    DO WHILE(k<6)
       READ(10,*)
       k=k+1
    END DO
    PRINT*,'k=',k
    T=0
    DO WHILE(k<nbr_node+6)
       read(10,*)num_node,x,y,z
       T=T+1
       T_node%n_node(num_node)=num_node
       T_node%c_node(num_node,1)=x
       T_node%c_node(num_node,2)=y
       T_node%c_node(num_node,3)=z
       IF(RANG==0)THEN
          PRINT*,'k=',k,'  ',T_node%n_node(num_node),T_node%c_node(num_node,:)
       END IF
       !print*,'k=',k,'  ',num_node,x,y,z
       k=k+1
    END DO
    PRINT*,k
    !k=10187
    read(10,*) !lit $EndNodes
    k=k+1 
    PRINT*,k
    read(10,*) !lit $Elements
    k=k+1 
    PRINT*,k
    read(10,*) !lit NBR d'ele total
    k=k+1
    PRINT*,k,'avvv'
    DO WHILE(k<MARK1)
       read(10,*)
       k=k+1
    END DO

    
    j=1
    DO WHILE(k<MARK2+1)
       k=k+1 !on lit cette ligne
       read(10,*)n1,n2,n3,n4,n5,n6,n7,n8
       !PRINT*,k,'  ',n1,n2,n3,n4,n5,n6,n7,n8
       T_conn(j,1)=n1  !num ele
       T_conn(j,2)=n6;T_conn(j,3)=n7;T_conn(j,4)=n8 !num sommets
       T_conn(j,5)=n5 !num elementary tag
       IF(rang==0)THEN
          PRINT*,j,' ',T_conn(j,1),T_conn(j,2),T_conn(j,3),T_conn(j,4),T_conn(j,5)
       END IF
       j=j+1
       
    END DO

    close(10)
    PRINT*,'NBR NODES= ',nbr_node,' T=',T,' j=',j-1
  end subroutine read_mesh

  SUBROUTINE READ_DATA(rang,file_name,mesh_name,Nd,N_TRI,FIRST_ROW,LAST_ROW,KMAX,TOL&
       ,PMAX,STEP,N_STEP,MODEL,LONG,LARG,RAY,condi_x,condi_y)!,DIR)
    INTEGER,intent(in)::rang
    CHARACTER(LEN=*),intent(in)::file_name
    
    CHARACTER(LEN=*),intent(out)::mesh_name
    INTEGER,intent(out)::Nd,N_TRI,FIRST_ROW,LAST_ROW,KMAX!,DIR
    REAL(rki),intent(out)::TOL
    REAL(rki),intent(out)::PMAX
    INTEGER,intent(out)::STEP,N_STEP,MODEL
    REAL(rki),intent(out)::LONG,LARG,RAY,condi_x,condi_y
    OPEN(10+rang,file=file_name)
    read(10+rang,*)
    read(10+rang,*)mesh_name
    read(10+rang,*)
    read(10+rang,*)Nd
    read(10+rang,*)
    read(10+rang,*)N_TRI
    read(10+rang,*)
    read(10+rang,*)FIRST_ROW
    read(10+rang,*)
    read(10+rang,*)LAST_ROW
    read(10+rang,*)
    read(10+rang,*)KMAX
    read(10+rang,*)
    read(10+rang,*)TOL
    !read(10+rang,*)
    !read(10+rang,*)DIR
    read(10+rang,*)
    read(10+rang,*)PMAX
    read(10+rang,*)
    read(10+rang,*)STEP
    read(10+rang,*)
    read(10+rang,*)N_STEP
    read(10+rang,*)
    read(10+rang,*)MODEL
    read(10+rang,*)
    read(10+rang,*)LONG
    read(10+rang,*)
    read(10+rang,*)LARG
    read(10+rang,*)
    read(10+rang,*)RAY
    read(10+rang,*)
    read(10+rang,*)condi_x
    read(10+rang,*)
    read(10+rang,*)condi_y
    CLOSE(10+rang)
  END SUBROUTINE READ_DATA
    
  SUBROUTINE DEFORMEE_DPL(rang,INPUT_MESH,OUTPUT_DEF,OUTPUT_X,OUTPUT_Y,Nd,SOLX,SOLY)
    INTEGER,intent(in)::rang,Nd
    CHARACTER(LEN=*),intent(in)::INPUT_MESH
    CHARACTER(LEN=*),intent(out)::OUTPUT_DEF,OUTPUT_X,OUTPUT_Y
    REAL(rki),DIMENSION(1:Nd),intent(in)::SOLX,SOLY
    integer::k,a
    real(rki)::x,y,z
    character(len=25)::CHEMIN
    CHEMIN='./MESH/'//INPUT_MESH
    IF(rang==0)THEN
    OPEN(10,file=CHEMIN)!INPUT_MESH)
    OPEN(20,file=OUTPUT_DEF)
    READ(10,*);READ(10,*);READ(10,*);READ(10,*);READ(10,*)
    DO k=1,Nd
       READ(10,*)a,x,y,z
       WRITE(20,*)a,x+SOLX(k),y+SOLY(k),z
    END DO
    CLOSE(20)
    CLOSE(10)
    ELSE IF(rang==1)THEN
       OPEN(30,file=OUTPUT_X)
       DO k=1,Nd
          WRITE(30,*)k,SOLX(k)
       END DO
       CLOSE(30)
       ELSE IF(rang==2)THEN
          OPEN(40,file=OUTPUT_Y)
          DO k=1,Nd
             WRITE(40,*)k,SOLY(k)
          END DO
          CLOSE(40)
    END IF
  END SUBROUTINE DEFORMEE_DPL




SUBROUTINE MAP_RESIDUS(rang,OUTPUT_X,OUTPUT_Y,Nd,SOLX,SOLY)
    INTEGER,intent(in)::rang,Nd
    
    CHARACTER(LEN=*),intent(out)::OUTPUT_X,OUTPUT_Y
    REAL(rki),DIMENSION(1:Nd),intent(in)::SOLX,SOLY
    integer::k
    real(rki)::x,y,z
    
    IF(rang==0)THEN
       OPEN(30,file=OUTPUT_X)
       DO k=1,Nd
          WRITE(30,*)k,SOLX(k)
       END DO
       CLOSE(30)
       ELSE IF(rang==3)THEN
          OPEN(40,file=OUTPUT_Y)
          DO k=1,Nd
             WRITE(40,*)k,SOLY(k)
          END DO
          CLOSE(40)
    END IF
  END SUBROUTINE MAP_RESIDUS


end Module ModMesh
