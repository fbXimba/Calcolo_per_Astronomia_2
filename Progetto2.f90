!gfortran -fbounds-check -o Progetto2.exe Progetto2.f90
!./Progetto2.exe
MODULE var
 IMPLICIT NONE
 
 INTEGER, PARAMETER::nd=6
 REAL*8::m(nd,nd)
 REAL*8,DIMENSION(nd)::c,nx1,nx2,nx
 REAL*8::Y1,y_1,Y2,y_2
 REAL*8,PARAMETER::mp=1.66d-24,k=1.38064852d-16   !c.g.s.
 !alcune variabili ridichiarate per intent

 CONTAINS 

END MODULE var

!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------

MODULE MOD_L
 USE var
 IMPLICIT NONE
 
 CONTAINS

SUBROUTINE sistema(T,m,c,y_)
 IMPLICIT NONE 
 REAL*8,INTENT(INOUT)::m(nd,nd),c(nd)
 REAL*8,INTENT(IN)::T,y_
 REAL*8::d1,d2,d3
  
  d1=max(b1(t),a1(T))
  d2=max(b2(T),(a2(t)+a3(T)))
  d3=max(b3(T),a4(T))
  
  !scrivo la matrice : il sistema
  !le prime 3 righe sono /d per avere valori migliori per gauss (jordan)
  !per evitare denominatori a 0: 1=1, 2=3, 3=4, 4=2, 5=6, 6=5
  
  !definisco m
  m=0.d0
  !equazione invariata
  m(1,1)=b1(T)/d1
  m(1,2)=-a1(T)/d1
  !equazione 2 = ex 4
  m(2,1)=1.d0
  m(2,2)=1.d0
  !equazione 3 = ex 2
  m(3,3)=b2(T)/d2
  m(3,4)=-(a2(T)+a3(T))/d2
  !equazione 4 = ex 3
  m(4,4)=-b3(T)/d3
  m(4,5)=a4(T)/d3
  !equazione 5 = ex 6
  m(5,3)=1.d0
  m(5,4)=1.d0
  m(5,5)=1.d0
  !equazione 6 = ex 5
  m(6,2)=-1.d0
  m(6,4)=-1.d0
  m(6,5)=-2.d0
  m(6,6)=1.d0
 
 !definisco c 
  c=0.d0
  c(2)=1.d0
  c(5)=y_

END SUBROUTINE sistema

!------------------------------------------------------------------------------------------------------------
SUBROUTINE jordan(a,c,x) 
  IMPLICIT NONE
  REAL*8,INTENT(INOUT):: a(nd,nd), c(nd)   !modifica la matrice
  REAL*8,INTENT(OUT):: x(nd)
  REAL*8::fakt,aux
  INTEGER::i,j,q
  
  DO i=1,nd ! elimino la variabile i
    aux=a(i,i)
    DO j=1,nd ! scelgo la riga su cui eliminare
       a(i,j)=a(i,j)/aux               
    END DO
    c(i)=c(i)/aux
     
    DO j=1,nd ! scelgo la riga su cui eliminare
      IF(i/=j) THEN 
        fakt=a(j,i)/a(i,i)      
        DO q=1,nd ! agisco su tutti i termini della riga
          a(j,q)=a(j,q)-a(i,q)*fakt
        END DO
        c(j)=c(j)-fakt*c(i)
      END IF
    END DO
  END DO  
  x=c
   
END SUBROUTINE jordan

!------------------------------------------------------------------------------------------------------------
SUBROUTINE cooling(Ltot,L_v,nx,T)
 IMPLICIT NONE
 REAL*8::lex,lion,lrec,ldrec,lff
 REAL*8,DIMENSION(6),INTENT(IN)::nx
 REAL*8,INTENT(IN)::T
 REAL*8,INTENT(OUT)::Ltot
 REAL*8,DIMENSION(5),INTENT(OUT)::L_v
 REAL*8::a1,a2,b1,b2,b3,c1,c2,c3,gff
 !gli n sono /nH ----> tutte le componenti sono /nH**2, nH=1 per i primi due step quindi non cambiano in modulo
 
 a1=nx(1)*7.50d-19*(exp(-118348.0d0/t))/(1.d0+(t/1.d5)**(0.5d0))
 a2=nx(4)*5.54d-17*(t**(-0.397d0))*(exp(-473638.0d0/t))/(1.d0+(t/1.d5)**(0.5d0))
 lex=nx(6)*(a1+a2)
 L_v(1)=lex
 
 b1=nx(1)*1.27d-21*(t**(0.5d0))*(exp(-157809.1d0/t))/(1.d0+(t/1.d5)**(0.5d0))
 b2=nx(3)*9.38d-22*(t**(0.5d0))*(exp(-285335.4d0/t))/(1.d0+(t/1.d5)**(0.5d0))
 b3=nx(4)*4.95d-22*(t**(0.5d0))*(exp(-631515.0d0/t))/(1.d0+(t/1.d5)**(0.5d0))
 lion=nx(6)*(b1+b2+b3)
 L_v(2)=lion
 
 c1=nx(2)*8.70d-27*(t**(0.5d0))*((t/1.d3)**(-0.2d0))/(1.d0+(t/1.d6)**(0.7d0))
 c2=nx(4)*1.55d-26*(t**(0.3647d0))
 c3=nx(5)*3.48d-26*(t**(0.5d0))*((t/1.d3)**(-0.2d0))/(1.d0+(t/1.d6)**(0.7d0))
 lrec=nx(6)*(c1+c2+c3)
 L_v(3)=lrec
 
 ldrec=nx(4)*nx(6)*1.24d-13*(t**(-1.5d0))*(exp(-470000.0d0/t))*(1.d0+0.3d0*(exp(-94000.0d0/t)))
 L_v(4)=ldrec
 
 gff=1.1d0+0.34d0*(exp(-((5.5d0-log10(t))**2.d0)/3.d0))
 lff=(nx(2)+nx(4)+4.d0*nx(5))*nx(6)*1.42d-27*gff*(t**(0.5d0))
 L_v(5)=lff
 
 Ltot=lex+lion+lrec+ldrec+lff
 
END SUBROUTINE cooling

!------------------------------------------------------------------------------------------------------------
SUBROUTINE s_tc(ne,nH,X,L,T,mu,d,u,t_cool)  !ne=ne/nH
 IMPLICIT NONE
 REAL*8,INTENT(IN)::ne,nH,X,T,L
 REAL*8,INTENT(OUT)::mu,d,u,t_cool
 REAL*8::y_,t_cool0
 
 y_=(1.d0-X)/(4.d0*X)
 
 mu=(1.d0+4.d0*y_)/(1.d0+y_+ne)   !mu=(1+4y)/(1+y+ne/nH)
 
 d=nH*mp/X
 
 !ua=(1.d0/(5.d0/3.d0-1.d0))*((k*Tmax)/(mu*mp))
 u=(1.d0/(5.d0/3.d0-1.d0))*((k*T)/(mu*mp))
 
 !t_coola=ua/(L/d)5
 t_cool=u/((L*nH*nH)/d)

END SUBROUTINE s_tc

!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------

!costanti dipendenti dalla temperatura t a1, a2, a3, a4, b1, b2, b3
!a=alpha H+,He+,d,He++ riocmbinazione

  FUNCTION a1(T)
    IMPLICIT NONE
    REAL*8,INTENT(IN)::T
    REAL*8::a1
    
    a1=8.4d-11*(T**(-0.5d0))*((T/1.d3)**(-0.2d0))/(1+(T/1.d6)**(0.7d0))
    
  END FUNCTION a1
  
  FUNCTION a2(T)
    IMPLICIT NONE
    REAL*8,INTENT(IN)::T
    REAL*8::a2
    
    a2=1.5d-10*(T**(-0.6353d0))
    
  END FUNCTION a2

  FUNCTION a3(T)
    IMPLICIT NONE
    REAL*8,INTENT(IN)::T
    REAL*8::a3
    
    a3=1.9d-3*(T**(-1.5d0))*(exp(-470000.0d0/T))*(1+0.3d0*(exp(-94000.0d0/T)))
    
  END FUNCTION a3

  FUNCTION a4(T)
    IMPLICIT NONE
    REAL*8,INTENT(IN)::T
    rEAL*8::a4
    
    a4=3.36d-10*(T**(-0.5d0))*((T/1.d3)**(-0.2d0)/(1+(T/1.d6)**(0.7d0)))  
    
  END FUNCTION a4

!b=GAMMA eH0,eHe0,eHe+ ionizzazione
  FUNCTION b1(T)
    IMPLICIT NONE
    REAL*8,INTENT(IN)::T
    REAL*8::b1
    
    b1=5.85d-11*(T**0.5d0)*(exp(-157809.1d0/T))/(1+(T/1.d5)**(0.5d0))
    
  END FUNCTION b1

  FUNCTION b2(T)
    IMPLICIT NONE
    REAL*8,INTENT(IN)::T
    REAL*8::b2
    
    b2=2.38d-11*(T**0.5d0)*(exp(-285335.4d0/T))/(1+(T/1.d5)**(0.5d0))
    
  END FUNCTION b2

  FUNCTION b3(T)
    IMPLICIT NONE
    REAL*8,INTENT(IN)::T
    REAL*8::b3
    
    b3=5.68d-12*(T**0.5d0)*(exp(-631515.0d0/T))/(1+(T/1.d5)**(0.5d0))
    
  END FUNCTION b3

END MODULE MOD_L

!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------

MODULE x_T
 USE MOD_L
 IMPLICIT NONE
 
 CONTAINS
 
  SUBROUTINE lambda_T(T,L)
   IMPLICIT NONE
   REAL*8,INTENT(IN)::T
   REAL*8,INTENT(OUT)::L
   REAL*8::y_,X
   REAL*8,DIMENSION(5)::o5
   X=0.76d0 !usiamo solo questa X nello step 4/5
   y_=(1.d0-X)/(4.d0*X)       !y_=Y/(4.d0*(1.d0-Y))
  
   CALL sistema(T,m,c,y_)
   CALL jordan(m,c,nx)
   CALL cooling(L,o5,nx,T)

  END SUBROUTINE lambda_T

!------------------------------------------------------------------------------------------------------------
  SUBROUTINE mu_T(T,mu)
   IMPLICIT NONE
   REAL*8,INTENT(IN)::T
   REAL*8,INTENT(OUT)::mu
   REAL*8::y_,X
   X=0.76d0 !usiamo solo questa X quando serve
   y_=(1.d0-X)/(4.d0*X)       !y_=Y/(4.d0*(1.d0-Y))
    
   CALL sistema(T,m,c,y_)
   CALL jordan(m,c,nx)

   mu=(1.d0+4.d0*y_)/(1.d0+y_+nx(6))

  END SUBROUTINE mu_T

END MODULE x_T

!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
MODULE secante
 USE x_T
 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE raphson(x0,u)  !T=x0
  IMPLICIT NONE
  REAL*8,EXTERNAL:: fun, derivata
  REAL*8,INTENT(INOUT)::x0
  REAL*8,INTENT(IN)::u
  REAL*8,PARAMETER::toll=1.d-12
  INTEGER,PARAMETER::max_iter=200
  REAL*8::xnew,shift
  INTEGER::iter
  
  iter=0
  DO
    iter=iter+1
    IF(iter>max_iter) THEN
      EXIT
    END IF
    xnew=x0-fun(x0,u)/derivata(x0,u)
    shift=ABS(x0-xnew)
    IF((shift<toll).or.(iter==max_iter)) THEN
      EXIT
    END IF
    x0=xnew
  END DO
  
  IF(x0<1.d4)THEN
    x0=1.d4 
  END IF
 
 END SUBROUTINE raphson


END MODULE secante

!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
MODULE s4
 USE secante
 IMPLICIT NONE
 
 CONTAINS

  SUBROUTINE dxdy(T,L0,u0,y,f)  !RHS    !y=ua
   REAL*8,INTENT(IN)::y,L0,u0
   REAL*8,INTENT(OUT)::f
   REAL*8,INTENT(INOUT)::T
   REAL*8::L
    
   !SECANTE u--->T
   CALL raphson(T,(y*u0))
   
   CALL lambda_T(T,L)
  
  !Ho L come L/nH**2 ma (L/nH**2)/(L0/nH**2)=L/L0 per qualsiasi valore di nH costante
   f=-L/L0

  END SUBROUTINE dxdy
  
!------------------------------------------------------------------------------------------------------------
  SUBROUTINE heun(T,L0,u0,h,yold,ynew)  !ODE_SOLVER,RK2   
   REAL*8, INTENT(IN)::h ,L0,u0
   REAL*8,INTENT(IN)::yold
   REAL*8,INTENT(OUT)::ynew
   REAL*8,INTENT(INOUT)::T
   REAL*8::f,fnew

   CALL dxdy(T,L0,u0,yold,f)
   
   ynew=yold+h*f   
         
   CALL dxdy(T,L0,u0,ynew, fnew)
  
   ynew=yold+ 0.5d0*h*(f+fnew) 
       
  END SUBROUTINE heun
  
!------------------------------------------------------------------------------------------------------------  
  SUBROUTINE ode(xmin,xmax,L0,tc,u0,a,b,x,y)    !a=y=ua, b=T
   REAL*8,DIMENSION(0:1),INTENT(OUT)::x
   REAL*8,DIMENSION(0:1),INTENT(OUT)::y
   REAL*8,INTENT(IN)::xmin,xmax,L0,tc,u0
   INTEGER::i
   REAL*8::h
   REAL*8,INTENT(OUT)::a,b
   
   h=((xmax-xmin)/tc)

   y(0)=a
   
   DO i=0,1
     x(i)=xmin/tc+i*h
   END DO 

   CALL heun(b,L0,u0,h, y(0), y(1)) 
   
  END SUBROUTINE ode
  
END MODULE s4

!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
PROGRAM PROGETTO2
 USE s4
 IMPLICIT NONE
 INTEGER::j,q,n
 REAL*8::T,interv,h,nH
 REAL*8, PARAMETER::X1=0.76d0,X2=1.d0 
 REAL*8,DIMENSION(5)::L_v1,L_v2
 REAL*8::L1,L2,L
 REAL*8::d,mu,u,t_cool
 REAL*8::L0,tc,t1,t2,time,ui,u0,ua,uii
 REAL*8,DIMENSION(2)::T_v,time_v,ui_v
 REAL*8,DIMENSION(3)::nH_v
 
 !y_=nHe/nH=Y/(4.d0*(1.d0-Y)), quando nH=1--->y_=nHe 
 Y1=1.d0-X1
 y_1=Y1/(4.d0*(1.d0-Y1))
 Y2=1.d0-X2                   != 0           
 y_2=Y2/(4.d0*(1.d0-Y2))      != 0
 
 !file di output con i risultati delle n (per X=1 le concentrazioni di He saranno 0)
 OPEN(23, file="dati1_1.dat")
 WRITE(23,*)"#   N         T          nH0/nH        nH1/nH       nHe0/nHe      nHe1/nHe      nHe2/nHe        ne", &
 "         X=0.76  Y=0.24  nH=1  nHe=7.8947E-002"
 
 OPEN(24, file="dati1_2.dat")
 WRITE(24,*)"#   N         T          nH0/nH        nH1/nH        nHe0          nHe1          nHe2           ne", &
 "         X=1   Y=0   nH=1   nHe=0"
 
 !file di output per i coefficienti di cooling
 OPEN(25, file="dati2_1.dat")
 WRITE(25,*)"#   N        T        Λtot          Λex           Λion          Λrec          Λdrec         Λff", &
 "         X=0.76    (Λ/nH^2, nH=1)"
 
 OPEN(26, file="dati2_2.dat")
 WRITE(26,*)"#   N        T        Λtot          Λex           Λion          Λrec          Λdrec         Λff", &
 "         X=1     (Λ/nH^2, nH=1)"
 
 !file di output per lo step 4 e 5
 OPEN(30, file="dati4_1.dat")
 WRITE(30,*)"#   tempo                     T                        u                  nH=0.1   T6"
 
 OPEN(31, file="dati4_2.dat")
 WRITE(31,*)"#   tempo                     T                        u                  nH=1.0   T6"
 
 OPEN(32, file="dati4_3.dat")
 WRITE(32,*)"#   tempo                     T                        u                  nH=10.0  T6"
 
 OPEN(33, file="dati5_0.dat")
 WRITE(33,*)"#   tempo                     T                        u                  nH=1.0   T7"


 !200 T equidistanti in logaritmo
 interv=(8.d0-4.d0)/199.d0
 DO q=1,200
 
   T=10**(4.d0+(q-1.d0)*interv)
   PRINT*,"Punto",q,"con T=",T
      
   !richesta 1:calcolo degli stati di ionizzazione (n/nH) con nH=1.d0 
   !X=0.76
   CALL sistema(T,m,c,y_1)
   CALL jordan(m,c,nx1)
  
   !X=1
   CALL sistema(T,m,c,y_2)  !sovrascrivo m e c  
   CALL jordan(m,c,nx2)
  
  
   !scrivo i risultati sui file di output 1
   717 FORMAT(I6,F14.2,6(1pe14.4))
  
   WRITE(23,717) q,T,nx1(1),nx1(2),(nx1(3)/y_1),(nx1(4)/y_1),(nx1(5)/y_1),nx1(6)   !nHe=y_1  (in modulo perchè nH=1)
  
   WRITE(24,717) q,T,nx2(1),nx2(2),nx2(3),nx2(4),nx2(5),nx2(6)      !nHe=y_2=0  (in modulo perchè nH=1)


   !richiesta 2: calcolo delle componenti del cooling e somma totale
   !X=0.76
   CALL cooling(L1,L_v1,nx1,T)
  
   !X=1
   CALL cooling(L2,L_v2,nx2,T)
  
   !scrivo i risultati sui file di output 2
   727 FORMAT(I6,F10.2,6(1pe14.4))
  
   WRITE(25,727) q,(log10(T)),((log10(L1))),(log10(L_v1(1))),(log10(L_v1(2))),(log10(L_v1(3))),(log10(L_v1(4))),(log10(L_v1(5)))
  
   WRITE(26,727) q,(log10(T)),((log10(L2))),(log10(L_v2(1))),(log10(L_v2(2))),&
   (log10(L_v2(3))),(log10(L_v2(4))),(log10(L_v2(5)))
  
 END DO

 CLOSE(23)
 CLOSE(24)
 CLOSE(25)
 CLOSE(26)
 
 PRINT*," "
 PRINT*,"I risultati sono salvati nei file: dati1_1.dat, dati1_2.dat, dati2_1.dat, dati2_2.dat"
 PRINT*," "

 !le 3 nH, X=X1=0.76
 nH_v=(/1.d-1,1.d0,1.d1/)
 T_v=(/1.d6,1.d7/)
 
 
 DO q=1,2 !T 6,7
  DO j=1,3 !nH 0.1,1.0,10
    
    IF((q==2).and.(j/=2)) THEN
     CONTINUE   !per T7 non devo fare nH=0.1,10.0
    ELSE        ! T6 con nH=0.1,1.0,10.0 e T7 con nH=1.0
      
      T=T_v(q)
        
      nH=nH_v(j)
        
      !calcolo le condizione per la dimensionalizzazione ed il ciclo a T6
      CALL lambda_T(T,L)
      L0=L
    
      CALL s_tc(nx(6),nH,X1,L,T,mu,d,u,t_cool)  !nh=1
      tc=t_cool
      u0=u 
    
      t1=0.d0
      t2=0.d0
      PRINT*,"Calcoli per: T=",T,"nH=",nH
      
      DO WHILE(t2<10.d0*tc)
        
        !X=0.76
        CALL sistema(T,m,c,y_1)
        CALL jordan(m,c,nx)  
        CALL cooling(L,L_v1,nx,T)
        CALL s_tc(nx(6),nH,X1,L,T,mu,d,u,t_cool)
        
        !time-step, intervallo 
        h=1.d-3*t_cool
        t2=t1+h
        
        !adimensionalizzo u perchè la ode è adimensionalizzata
        ua=u/u0
        
        !risoluzione di du/dt adimensionalzzata
        CALL ode(t1,t2,L0,tc,u0,ua,T,time_v,ui_v)
        ui=ui_v(2)     !adimensionale
        time=time_v(2)*tc
        
        !ridimensionalizzo la u per raphson
        uii=ui*u0  
        
        !metodo della secante per stimare T dalla nuova u  
        CALL raphson(T,uii)
        
        !scrivo i risultati
        n=27+j+2*q
        WRITE(n,*)(log10(time/tc)),(log10(T)),(ui)
        
        t1=t2
      END DO
    
    END IF 
  END DO  
 END DO  
 
 CLOSE(30)
 CLOSE(31)
 CLOSE(32)
 CLOSE(34)
 
 PRINT*," "
 PRINT*,"I risultati sono salvati nei file: dati4_1.dat, dati4_2.dat, dati4_3.dat, dati5_0.dat"
  
END PROGRAM PROGETTO2

!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
REAL*8 FUNCTION fun(T,u)
  USE x_T
  IMPLICIT NONE
  REAL*8::T,u
  REAL*8::mu
  CALL mu_T(T,mu)
  fun=T-(5.d0/3.d0-1.d0)*u*mu*mp/k
END FUNCTION fun

!------------------------------------------------------------------------------------------------------------
REAL*8 FUNCTION derivata(T,u)
  IMPLICIT NONE
  REAL*8, PARAMETER:: h=1.d-4
  REAL*8, EXTERNAL:: fun
  REAL*8:: T,u
  derivata=(fun(T+h,u)-fun(T,u))/h
END FUNCTION derivata

!------------------------------------------------------------------------------------------------------------
!--------------------------------------------the-end---------------------------------------------------------
