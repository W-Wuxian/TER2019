MODULE ModParra
use mpi

implicit none


CONTAINS

subroutine CHARGE(rang,Np,it1,itN,CHARGE_TOT,enplus)
    integer,intent(in)::rang,Np,CHARGE_TOT,enplus
    integer,intent(out)::it1,itN
    real::coeff_repartition

    coeff_repartition=CHARGE_TOT/Np
    SELECT CASE(Np)
    CASE(1)
       it1=1
       itN=CHARGE_TOT
    CASE DEFAULT 
       if(rang<mod(CHARGE_TOT,Np))then
          it1=rang*(coeff_repartition+1)+1
          it1=it1+enplus
          itN=(rang+1)*(coeff_repartition+1)
          itN=itN+enplus
          print*,'RANGRANG',rang,' coeff rep=',coeff_repartition,' mod(chargetot,Np)', mod(CHARGE_TOT,Np)

       else
          it1=1+mod(CHARGE_TOT,Np)+rang*coeff_repartition
          it1=it1+enplus
          itN=it1+coeff_repartition-1
          !itN=itN
       end if
    END SELECT
    print*,'je suis proc',rang,'je vais de',it1,'à',itN
  end subroutine CHARGE

END MODULE ModParra
