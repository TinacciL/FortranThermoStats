  program FortranThermoStats
  use parametri  
   
  implicit none 
 
  real(8) :: x(n),y(n),z(n), freq(n), theta(n), fvib(n),theta_tot, fp_vib, Ain, Bin, Cin, Iin, &
             fp_rot, mass(n), fp_tot, fp_el, fp_nuc, fp_trl,                                   &
             ix, iy, iz, ixz, ixy, iyz, Inerzia(3,3), Inerzia1(3,3), M, bx, by, bz, pi, ht, A, &
             S_tot, S_vib, S_trl, S_rot, G, V, P, E_tot, E_vib, E_trl, E_rot, Cv, Cp, T, mol,  &
             E_totj, E_vibj, E_trlj, E_rotj, S_totj, S_vibj, S_trlj, S_rotj, E_vibs, E_vibsj,  &
             Cvr
  character(2) :: nome(n)
  integer :: k, i, nato, j, sig, sim, nmv

  !---caso in cui voglia leggere e stampare su file-------------------------------------------
  !open (unit = 1, file = 'risultati.dat', status = 'replace', action = 'write')
  !open (unit = 2, file = 'input', status = 'old', action = 'read')  
  !-------------------------------------------------------------------------------------------
  
  pi = acos(-1.d0)
  ht =  h/(2.d0*pi) 
  
  !---letture dati ingresso-------------------------------------------------------------------
  read (*,*) !commento input
  read (*,*) T
  read (*,*) P
  read (*,*) V
  read (*,*) mol
  read (*,*) nato
  read (*,*) sim
  read (*,*) sig

  do i=1, nato 
     read (*,*) nome(i), x(i), y(i), z(i)
  enddo

   k=0

   do 
      k= k+1
      read(*,*,end=10) freq(k)
   enddo
  !---fine input-----------------------------------------------------------------------------
   
  !----associazione vettore nome a vettore masse e calcolo massa molecola M------------------
10 do i=1, nato
      do j=1, 118
         if (nome(i) == element(j)) then
            mass(i)  = weights(j) 
         endif
      enddo
   enddo   

   do i=1, nato
      M= M + mass(i)
   enddo   
  !-------------------------------------------------------------------------------------------  
   
  !Parametri sistema da legge gas ideali
  if (P == 0) then
     P = (mol * Av *kb * T) / (V * 101325)
  endif
  if (V == 0) then 
     V = (mol * Av *kb * T) / (P * 101325)
  endif
  if (mol == 0) then
     mol = (P * V * 101325) / (Av *kb * T)
  endif

  !assegnare modi nomali di vibrazione
  if (sim == 0) then
     nmv = 3 * nato - 5 
  else if (sim == 1) then
     nmv = 3 * nato - 6
  endif    
  
  !---------------Inerzia-------------------------------------------------------------- 
  !calcolo baricentro e cambio sistema riferimento 
   do i=1, nato
      bx = bx + mass(i)*x(i)/M
      by = by + mass(i)*y(i)/M
      bz = bz + mass(i)*z(i)/M
   enddo
          
     x(i) = x(i) - bx
     y(i) = y(i) - bx
     z(i) = z(i) - bx
         
  !costruzione elementi tensore inerzia e diagonalizzazione
    do i=1, nato
       ix = ix + mass(i) * ((y(i)**2) + (z(i)**2))
       iy = iy + mass(i) * ((z(i)**2) + (x(i)**2))
       iz = iz + mass(i) * ((y(i)**2) + (x(i)**2))
       ixz = ixz + mass(i) * x(i) * z(i)
       ixy = ixy + mass(i) * x(i) * y(i)
       iyz = iyz +  mass(i) * y(i) * z(i)
    enddo   
    
     Inerzia(1,1) = ix
     Inerzia(1,2) = ixy
     Inerzia(1,3) = ixz
     Inerzia(2,1) = ixy
     Inerzia(2,2) = iy
     Inerzia(2,3) = iyz
     Inerzia(3,1) = ixz
     Inerzia(3,2) = iyz
     Inerzia(3,3) = iz

    call jacobi (Inerzia, Inerzia1, abserr, 3)
     
   ! assi principali di inerzia in  kg e metri da uma e amstrong 
     Ain = Inerzia(1,1) * 1.660539d-47 
     Bin = Inerzia(2,2) * 1.660539d-47 
     Cin = Inerzia(3,3) * 1.660539d-47  
   !---------------------------------------------------------------------------------------------------  
     
   !---Funzione partizione-----------------------------------------------------------------------------     
   !NOTA: vib= vibrazionale, trl=traslazionale, el=eletronica, nuc=nucleare,rot= rotazionale  tot=totale

   !vibrazionale e temperature vibrazionali
   fp_vib = 1.d0

   do k=1, nmv 
      theta(k)  =  freq(k) * c * 100* h / Kb
      fvib(k)   =  1 / (1 - exp(-theta(k) / T))
      fp_vib    =  fp_vib * fvib(k)
      theta_tot = theta_tot + theta(k)
   enddo
  
     
   !rotazionale lineare o non 
   if (sim == 0) then                                     
      I = Ain + Bin + Cin
      fp_rot = 2 * T * I *kb / (sig * h**2)                                   
   else if (sim == 1) then
     fp_rot = ( ( ( (2*kb*T) / ( (ht)**2) ) )**(1.5d0)) *  ( (pi*Ain*Bin*Cin)**(0.5d0) ) /sig
   endif     

 
     fp_trl = V /  h**3  *  (2 * pi * kb * T * M * 1.660539e-27)**(1.5d0) / Av

     fp_el = 1

     fp_nuc = 1
          

     fp_tot = fp_vib + fp_rot + fp_el + fp_nuc + fp_trl
   !----fine funzione partizione-----------------------------------------------------------------------  
     
   !----energia in termodinamica classica--------------------------------------------------------------
   !NOTA: suffisso -j indica  joule/mol, senza kcal/mol

     !energia rotazionale per molecola lineare o non in joule/mol    
   if (sim == 0) then     
      E_rotj  = (R * T) / mol                                   
   else if (sim == 1) then
      E_rotj  = (1.5d0 * R * T) / mol                     
   endif

   E_rot   = E_rotj * 0.2388459 *0.001                    
   E_trlj  = (1.5d0 * R * T) /mol                          
   E_trl   = E_trlj * 0.2388459 *0.001                    
   E_vibj  = (R * T) / mol                               
   E_vib   = E_vibj * 0.2388459 * 0.001                   

   E_totj  = (E_vibj + E_rotj + E_trlj)                    
   E_tot   = (E_vib + E_rot + E_trl)               
   !---------------------------------------------------------------------------------------------------   
     
   !grandezze termodinamiche
   
    S_vibj = Av * kb * log(fp_vib) + E_vibj / T
    S_vib  = S_vibj * 0.2388459
    S_rotj = Av * kb * log(fp_rot) + E_rotj / T
    S_rot  = S_rotj * 0.2388459
    S_trlj = Av * kb * log(fp_trl) + E_trlj / T
    S_trl  = S_trlj * 0.2388459
    S_tot  = S_vib + S_rot + S_trl
    S_totj = S_vibj + S_rotj + S_trlj     

   !Calcolo Cv 
    do i=1, nmv
       Cv = Cv + R * (theta(i) / T)**2.d0 * exp(theta(i) / T ) / ( exp(theta(i) / T) - 1)**2.d0
    enddo

    Cvr = 1.5d0 * R
    Cv  = ( Cv + 1.5d0 * R + Cvr) * 0.2388459


    !----Dati in OutPut---------------------------------------------------------------------------------
    write (*,*) 'Temperatura (K):              ', T
    write (*,*) 'Volume sistema (m^3):         ', V
    write (*,*) 'Pressione Sistema (atm):      ', P
    write (*,*) 'Moli Sistema':                ', mol
    write (*,*) 'Numero di Atomi:              ', nato
    write (*,*) 'Numero Modi Normali:          ', k-1
    write (*,*) 'Assi di Simmetria Molecola:   ', sig
    
    if (sim == 0) then
       write (*,*) 'Simmetria Molecola:           ', 'Lineare'
    else if (sim == 1) then
       write (*,*) 'Simmetria Molecola:           ', 'Non Lineare'
    endif
    
    if (k-1 /= nmv) then
       write (*,*) 'ERROR: numero modi normali vibrazione diverso da frequenze in lettura'
       stop
    endif
    
    write (*,*) '                        Energia (kcal*mol^-1)      Energia (J*mol^-1) '
    write (*,*) 'Energia Vibrazionale:  ', E_vib,               E_vibj
    write (*,*) 'Energia Rotazionale:   ', E_rot,               E_rotj
    write (*,*) 'Energia Traslazionale: ', E_trl,               E_trlj
    write (*,*) 'Energia Totale Sistema:', E_tot,               E_totj

   do i=1, nmv
      write(*,*) 'Frequenza (cm^-1):', freq(i), 'Temperatura Vibrazionale (K):', theta(i)
   enddo
  
   do i=1, nato
      write (*,*) 'Atomo, Massa (u) e Coordinate (A):', nome(i), mass(i),  x(i), y(i), z(i)
   enddo
  
    write(*,*) 'Massa Molecola (u):', M     
    write(*,*) 'Baricentro (A):    '
    write(*,*)  bx, by, bz    
    write(*,*) 'Tensore Inerzia (u*A^2):'
    write(*,*)  Inerzia(1,1), Inerzia(1,2), Inerzia(1,3)
    write(*,*)  Inerzia(2,1), Inerzia(2,2), Inerzia(2,3)
    write(*,*)  Inerzia(3,1), Inerzia(3,2), Inerzia(3,3)    
    write(*,*) 'Assi Principali di Inerzia (Kg*m^2)'
    Write(*,*)  Ain
    Write(*,*)  Bin
    Write(*,*)  Cin
    write(*,*) 'Funzione Partizione Rotazionale:  ', fp_rot
    write(*,*) 'Funzione Partizione Vibrazionale: ', fp_vib
    write(*,*) 'Funzione Partizione Traslazionale:', fp_trl
    write(*,*) 'Funzione Partizione Elettronica:  ', fp_el
    write(*,*) 'Funzione Partizione Nucleare:     ', fp_nuc
    write(*,*) 'Funzione Partizione Totale:       ', fp_tot     
    write(*,*) '                          Entropia (cal*mol^-1*K^-1) Entropia (J*mol^-1*K^-1)'    
    write(*,*) 'Entropia Vibrazionale:  ', S_vib,               S_vibj
    write(*,*) 'Entropia Rotazionale:   ', S_rot,               S_rotj
    write(*,*) 'Entropia Traslazionale: ', S_trl,               S_trlj
    write(*,*) 'Entropia Totale Sistema:', S_tot,               S_totj
    write(*,*) 'Cv:',  Cv
    
end program FortranThermoStats




subroutine Jacobi(a,x,abserr,n)
!===========================================================
! Evaluate eigenvalues and eigenvectors
! of a real symmetric matrix a(n,n): a*x = lambda*x 
! method: Jacoby method for symmetric matrices 
! Alex G. (December 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - number of equations
! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! output ...
! a(i,i) - eigenvalues
! x(i,j) - eigenvectors
! comments ...
!===========================================================
implicit none
integer i, j, k
double precision  b2, bar
double precision beta, coeff, c, s, cs, sc
real(8), intent(inout) :: a(n,n), x(n,n)
real(8), intent(in) :: abserr
integer, intent(in) :: n

x = 0.0
do i=1,n
  x(i,i) = 1.0
end do

! find the sum of all off-diagonal elements (squared)
b2 = 0.0
do i=1,n
  do j=1,n
    if (i.ne.j) b2 = b2 + a(i,j)**2
  end do
end do

if (b2 <= abserr) return

! average for off-diagonal elements /2
bar = 0.5*b2/float(n*n)

do while (b2.gt.abserr)
  do i=1,n-1
    do j=i+1,n
      if (a(j,i)**2 <= bar) cycle  ! do not touch small elements
      b2 = b2 - 2.0*a(j,i)**2
      bar = 0.5*b2/float(n*n)
! calculate coefficient c and s for Givens matrix
      beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
      coeff = 0.5*beta/sqrt(1.0+beta**2)
      s = sqrt(max(0.5+coeff,0.0))
      c = sqrt(max(0.5-coeff,0.0))
! recalculate rows i and j
      do k=1,n
        cs =  c*a(i,k)+s*a(j,k)
        sc = -s*a(i,k)+c*a(j,k)
        a(i,k) = cs
        a(j,k) = sc
      end do
! new matrix a_{k+1} from a_{k}, and eigenvectors 
      do k=1,n
        cs =  c*a(k,i)+s*a(k,j)
        sc = -s*a(k,i)+c*a(k,j)
        a(k,i) = cs
        a(k,j) = sc
        cs =  c*x(k,i)+s*x(k,j)
        sc = -s*x(k,i)+c*x(k,j)
        x(k,i) = cs
        x(k,j) = sc
      end do
    end do
  end do
end do
return
end subroutine Jacobi
  
