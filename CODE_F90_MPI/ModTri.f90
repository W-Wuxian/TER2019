module ModTri
  use ModMesh
  use mpi
  implicit none

contains

  SUBROUTINE DOUBLON(SIZE,rang,T_CL,T_V,Z)
    INTEGER,intent(in)::SIZE,rang
    INTEGER,dimension(1:SIZE,1:2),intent(inout)::T_CL
    REAL(rki),dimension(1:SIZE),intent(inout)::T_V
    INTEGER,intent(out)::Z !NBR DE ZERO
    INTEGER::i,j
    !ON SUPPOSE QUE T_V(:)/=0
    DO j=1,SIZE-1
       IF(T_V(j)/=0)THEN
          DO i=j+1,SIZE
             IF(T_V(i)/=0)THEN
                !PRINT*,'RANG',rang,j,i,'DOUBLON'
                IF((T_CL(j,1)==T_CL(i,1)).AND.(T_CL(j,2)==T_CL(i,2)))THEN
                   T_V(j)=T_V(j)+T_V(i);T_V(i)=0
                END IF
             END IF
          END DO
       END IF
    END DO
    Z=0
    DO i=1,SIZE
       IF(T_V(i)==0.d0)THEN
          Z=Z+1
       END IF
    END DO
  END SUBROUTINE DOUBLON

  SUBROUTINE TRI_CR(SIZE,Z,rang,T_CL,T_V,S_CL,S_V,NNZ)!TRI OVER ALL T_CL(:,1) ET T_CL(:,2)
    INTEGER,intent(in)::SIZE,Z,rang
    INTEGER,dimension(1:SIZE,1:2),intent(in)::T_CL
    INTEGER,dimension(1:SIZE)::T_TRI
    REAL(rki),dimension(1:SIZE),intent(in)::T_V
    INTEGER,intent(out)::NNZ
    INTEGER,dimension(:,:),intent(out),allocatable::S_CL
    REAL(rki),dimension(:),intent(out),allocatable::S_V
    INTEGER::i,j,DIM,pt
    DIM=SIZE-Z
    ALLOCATE(S_CL(1:DIM,1:2),S_V(1:DIM))
    T_TRI=0;S_CL=0;S_V=0.d0
    DO j=1,SIZE
       !PRINT*,'RANG',RANG,'j=',j
       IF(T_V(j)/=0)THEN
          !PRINT*,'RANG',rang,'IF'
          DO i=1,SIZE
             !PRINT*,'RANG',rang,'i=',i
             IF(T_V(i)/=0)THEN
                !PRINT*,'RANG',rang,j,i,'TRI'
                IF(T_Cl(j,1)==T_CL(i,1).AND.T_CL(j,2)>=T_CL(i,2))THEN
                   T_TRI(j)=T_TRI(j)+1
                ELSE IF(T_CL(j,1)>T_CL(i,1))THEN
                   T_TRI(j)=T_TRI(j)+1
                END IF
             END IF
          END DO
       END IF
    END DO

    DO i=1,SIZE
       IF(T_TRI(i)/=0)THEN
          pt=T_TRI(i)
          S_CL(pt,1)=T_CL(i,1);S_CL(pt,2)=T_CL(i,2);S_V(pt)=T_V(i);
       END IF
    END DO
    NNZ=0
    DO i=1,DIM
       IF(S_V(i)/=0)THEN
          NNZ=NNZ+1
       END IF
    END DO

    PRINT*,'RANG',rang,'DIM',DIM,'NNZ',NNZ

  END SUBROUTINE TRI_CR



  


end module ModTri
