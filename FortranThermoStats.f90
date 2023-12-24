program FortranThermoStats
  use parameters

  implicit none

  real(8) :: x(n), y(n), z(n), freq(n), theta(n), fvib(n), theta_tot, fp_vib, Ain, Bin, Cin, Iin, &
             fp_rot, mass(n), fp_tot, fp_el, fp_nuc, fp_trl,                                   &
             ix, iy, iz, ixz, ixy, iyz, Inerzia(3,3), Inerzia1(3,3), M, bx, by, bz, pi, ht, A, &
             S_tot, S_vib, S_trl, S_rot, G, V, P, E_tot, E_vib, E_trl, E_rot, Cv, Cp, T, mol,  &
             E_totj, E_vibj, E_trlj, E_rotj, S_totj, S_vibj, S_trlj, S_rotj, E_vibs, E_vibsj,  &
             Cvr
  character(2) :: nome(n)
  integer :: k, i, nato, j, sig, sim, nmv
  
  !---case where reading and printing to file is desired-------------------------------------------
  open (unit = 1, file = 'risultati.dat', status = 'replace', action = 'write')
  open (unit = 2, file = 'input.txt', status = 'old', action = 'read')  
  !-------------------------------------------------------------------------------------------

  pi = acos(-1.d0)
  ht =  h/(2.d0*pi) 

  !---reading input data-------------------------------------------------------------------
  
  read (2,*) !comment input
  read (2,*) T
  read (2,*) P
  read (2,*) V
  read (2,*) mol
  read (2,*) nato
  read (2,*) sim
  read (2,*) sig
  
  do i=1, nato 
     read (2,*) nome(i), x(i), y(i), z(i)
  enddo

   k=0

   do 
      k= k+1
      read(2,*,end=10) freq(k)
   enddo

  !---end of input-----------------------------------------------------------------------------
   
  !----associating atom names with mass vector and calculating molecular mass M------------------
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
   
  !Parameters for ideal gas law
  if (P == 0) then
     P = (mol * Av *kb * T) / (V * 101325)
  endif
  if (V == 0) then 
     V = (mol * Av *kb * T) / (P * 101325)
  endif
  if (mol == 0) then
     mol = (P * V * 101325) / (Av *kb * T)
  endif
 
  !assigning normal vibrational modes
  if (sim == 0) then
     nmv = 3 * nato - 5 
  else if (sim == 1) then
     nmv = 3 * nato - 6
  endif    
  
  !---------------Inertia-------------------------------------------------------------- 
  !calculating center of mass and changing reference frame 
   do i=1, nato
      bx = bx + mass(i)*x(i)/M
      by = by + mass(i)*y(i)/M
      bz = bz + mass(i)*z(i)/M
   enddo
          
     x(i) = x(i) - bx
     y(i) = y(i) - bx
     z(i) = z(i) - bx
         
  !constructing inertia tensor elements and diagonalizing
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
     
   ! principal axes of inertia in kg and meters from uma and angstrom 
     Ain = Inerzia(1,1) * 1.660539d-47 
     Bin = Inerzia(2,2) * 1.660539d-47 
     Cin = Inerzia(3,3) * 1.660539d-47  
   !---------------------------------------------------------------------------------------------------  
     
   !---Partition Function-----------------------------------------------------------------------------     
   !NOTE: vib= vibrational, trl=translational, el=electronic, nuc=nuclear, rot= rotational  tot=total

   !vibrational and vibrational temperatures
   fp_vib = 1.d0

   do k=1, nmv 
      theta(k)  =  freq(k) * c * 100* h / Kb
      fvib(k)   =  1 / (1 - exp(-theta(k) / T))
      fp_vib    =  fp_vib * fvib(k)
      theta_tot = theta_tot + theta(k)
   enddo
  
   !linear or non-linear rotational
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
   !----end Partition Function-----------------------------------------------------------------------  
     
   !----energy in classical thermodynamics--------------------------------------------------------------
   !NOTE: suffix -j indicates  joule/mol, without kcal/mol

   !rotational energy for linear or non-linear molecule in joule/mol    
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
     
   !thermodynamic quantities
   
    S_vibj = Av * kb * log(fp_vib) + E_vibj / T
    S_vib  = S_vibj * 0.2388459
    S_rotj = Av * kb * log(fp_rot) + E_rotj / T
    S_rot  = S_rotj * 0.2388459
    S_trlj = Av * kb * log(fp_trl) + E_trlj / T
    S_trl  = S_trlj * 0.2388459
    S_tot  = S_vib + S_rot + S_trl
    S_totj = S_vibj + S_rotj + S_trlj     

   !calculating Cv 
    do i=1, nmv
       Cv = Cv + R * (theta(i) / T)**2.d0 * exp(theta(i) / T ) / ( exp(theta(i) / T) - 1)**2.d0
    enddo

    Cvr = 1.5d0 * R
    Cv  = ( Cv + 1.5d0 * R + Cvr) * 0.2388459


    !----Output Data---------------------------------------------------------------------------------
    write (1,*) 'Temperature (K):              ', T
    write (1,*) 'System Volume (m^3):          ', V
    write (1,*) 'System Pressure (atm):        ', P
    write (1,*) 'System Moles:                 ', mol
    write (1,*) 'Number of Atoms:              ', nato
    write (1,*) 'Number of Normal Modes:       ', k-1
    write (1,*) 'Molecule Symmetry Axes:       ', sig
    
    if (sim == 0) then
       write (1,*) 'Molecule Symmetry:            ', 'Linear'
    else if (sim == 1) then
       write (1,*) 'Molecule Symmetry:            ', 'Nonlinear'
    endif
    
    if (k-1 /= nmv) then
       write (*,*) 'ERROR: Number of normal vibrational modes different from read frequencies'
       stop
    endif
    
    write (1,*) '                        Energy (kcal*mol^-1)      Energy (J*mol^-1) '
    write (1,*) 'Vibrational Energy:    ', E_vib,               E_vibj
    write (1,*) 'Rotational Energy:     ', E_rot,               E_rotj
    write (1,*) 'Translational Energy:  ', E_trl,               E_trlj
    write (1,*) 'Total System Energy:   ', E_tot,               E_totj

   do i=1, nmv
      write (1,*) 'Frequency (cm^-1):', freq(i), 'Vibrational Temperature (K):', theta(i)
   enddo
  
   do i=1, nato
      write (1,*) 'Atom, Mass (u), and Coordinates (A):', nome(i), mass(i),  x(i), y(i), z(i)
   enddo
  
    write (1,*) 'Molecule Mass (u):', M     
    write (1,*) 'Center of Mass (A):'
    write (1,*)  bx, by, bz    
    write (1,*) 'Inertia Tensor (u*A^2):'
    write (1,*)  Inerzia(1,1), Inerzia(1,2), Inerzia(1,3)
    write (1,*)  Inerzia(2,1), Inerzia(2,2), Inerzia(2,3)
    write (1,*)  Inerzia(3,1), Inerzia(3,2), Inerzia(3,3)    
    write (1,*) 'Principal Axes of Inertia (Kg*m^2)'
    Write (1,*)  Ain
    Write (1,*)  Bin
    Write (1,*)  Cin
    write (1,*) 'Rotational Partition Function:  ', fp_rot
    write (1,*) 'Vibrational Partition Function: ', fp_vib
    write (1,*) 'Translational Partition Function:', fp_trl
    write (1,*) 'Electronic Partition Function:  ', fp_el
    write (1,*) 'Nuclear Partition Function:     ', fp_nuc
    write (1,*) 'Total Partition Function:       ', fp_tot     
    write (1,*) '                        Entropy (cal*mol^-1*K^-1) Entropy (J*mol^-1*K^-1)'    
    write (1,*) 'Vibrational Entropy:  ', S_vib,               S_vibj
    write (1,*) 'Rotational Entropy:   ', S_rot,               S_rotj
    write (1,*) 'Translational Entropy: ', S_trl,               S_trlj
    write (1,*) 'Total System Entropy: ', S_tot,               S_totj
    write (1,*) 'Cv:',  Cv
    
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