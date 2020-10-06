program poi
implicit none
integer::i,k
REAL(kind=8)::DX,DY,SUM_DX,SUM_DY,POISSON
POISSON=0.0d0;SUM_DX=0.0d0;SUM_DY=0.0d0
OPEN(10,file='DELLdef_xx.txt')
DO i=1,27430
   READ(10,*)k,DX
   SUM_DX=SUM_DX+DX
END DO
CLOSE(10)
OPEN(11,file='DELLdef_yy.txt')
DO i=1,27430
   READ(11,*)k,DY
   SUM_DY=SUM_DY+DY
END DO
CLOSE(11)
POISSON=-SUM_DX/SUM_DY
PRINT*,'POI=',POISSON


end program poi
