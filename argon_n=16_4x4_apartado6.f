      program argon
      implicit none
      
      real*8, dimension(:,:), allocatable :: r
      real*8, dimension(:,:), allocatable :: r_aux
      real*8, dimension(:,:), allocatable :: v
      real*8, dimension(:,:), allocatable :: acelx
      real*8, dimension(:,:), allocatable :: acely
      real*8, dimension(:,:), allocatable :: a1
      real*8, dimension(:,:), allocatable :: omega
      real*8  d, aux_a, rad, dt, t!, dran_u, ang_aux
      integer i,j,n,p,q
      real*8 pdec
      real*8 paux(2)
      !Variables para c lculos postenriores
      real*8 Ec,Ep,Et,Temp
      integer cont
      real*8, dimension(:,:), allocatable :: radio
      real*8 posini(2)
      logical lejos
      
      Temp=0.d0
      cont=0
      t=0.d0
      lejos=.FALSE.
      
      !call dran_ini(783247)
      
      dt=0.002
      n=16

      allocate(r(n,2))
      allocate(r_aux(n,2))
      allocate(v(n,2))
      allocate(acelx(n,n))
      allocate(acely(n,n))
      allocate(a1(n,2))
      allocate(omega(n,2))
      allocate(radio(n,n))

      d=4.d0 !Dimensiones (dxd) de mi cristal bidimensional
              !Debe ser un n£mero natural

      !Abro archivos
      
      open(4,file='acel.txt')
      open(2,file='argon.txt')
      
      open(5,file='energia.txt')
      open(1,file='prueba.txt')
      
      open(6,file='vel_t0.txt')
      open(7,file='vel_t1_1.txt')
      
      open(8,file='Temperatura.txt')
      
      open(9,file="fluctuacion_pos.txt")

      !Inicializo a cero la aceleraciones para luego sumarlas
      do i=1,n
         do j=1,n
          acelx(i,j)=0.d0  !Aceleraci¢n de i ejercida por j
          acely(i,j)=0.d0
          radio(i,j)=0.d0
         enddo
      enddo

      do i=1,n
         do j=1,2
          r_aux(i,j)=0.d0
          r_aux(i,j)=0.d0
         enddo
      enddo
      
      !Inicializo posiciones iniciales random uniforme
      r(1,1)=0.5d0
      r(1,2)=0.5d0
      r(2,1)=1.5d0
      r(2,2)=0.5d0
      r(3,1)=2.5d0
      r(3,2)=0.5d0
      r(4,1)=3.5d0
      r(4,2)=0.5d0
      r(5,1)=0.5d0
      r(5,2)=1.5d0
      r(6,1)=1.5d0
      r(6,2)=1.5d0
      r(7,1)=2.5d0
      r(7,2)=1.5d0
      r(8,1)=3.5d0
      r(8,2)=1.5d0
      r(9,1)=0.5d0
      r(9,2)=2.5d0
      r(10,1)=1.5d0
      r(10,2)=2.5d0
      r(11,1)=2.5d0
      r(11,2)=2.5d0
      r(12,1)=3.5d0
      r(12,2)=2.5d0
      r(13,1)=0.5d0
      r(13,2)=3.5d0
      r(14,1)=1.5d0
      r(14,2)=3.5d0
      r(15,1)=2.5d0
      r(15,2)=3.5d0
      r(16,1)=3.5d0
      r(16,2)=3.5d0
      
      posini(1)=r(6,1)
      posini(2)=r(6,2)

      write(2,*) ((r(i,j), j=1,2), i=1,n)
      write(2,*)
      write(2,*)

      !Inicializo las velocidades iniciales a modulo 1 y  ngulo random

      do i=1,n
         do j=1,2
           v(i,j)=0.d0
         enddo
      enddo
      !Histograma de velocidades en t=0
      do i=1,n
        write(6,*) (v(i,1)**2.d0+v(i,2)**2.d0)**0.5d0, v(i,1), v(i,2)
      enddo

      
      !Aceleraciones iniciales
      
      do i=1,n-1
         do j=i+1,n
            !Cond contorno en fuerza
            do q=1,2
               if ((r(j,q)-r(i,q)).gt.(d/2.d0)) then
                  paux(q)=r(j,q)-d
               else if ((r(j,q)-r(i,q)).lt.(-d/2.d0)) then
                  paux(q)=r(j,q)+d
               else
                  paux(q)=r(j,q)
               endif
            enddo
            !Sigo
            rad=((paux(1)-r(i,1))**2.d0+(paux(2)-r(i,2))**2.d0)**0.5d0
            radio(i,j)=rad
            if (3.ge.rad) then
              aux_a=24.d0*(2.d0*rad**(-13.d0)-rad**(-7.d0))
              acelx(i,j)=-aux_a*(paux(1)-r(i,1))/rad
              acely(i,j)=-aux_a*(paux(2)-r(i,2))/rad
            else
              acelx(i,j)=0.d0
              acely(i,j)=0.d0
            end if
            acelx(j,i)=-acelx(i,j)
            acely(j,i)=-acely(i,j)
            radio(j,i)=radio(i,j)
         enddo
      enddo
      !Aceleraci¢n total de cada part¡cula
      do i=1,n
        do j=1,n
           a1(i,1)=a1(i,1)+acelx(i,j)
           a1(i,2)=a1(i,2)+acely(i,j)
        enddo
      enddo

      
      !ENERGÖAS:
                !Inicializo
      Ec=0.d0
      Ep=0.d0
      Et=0.d0
                !Energ¡a cin‚tica:
      do i=1,n
         Ec=Ec+0.5d0*(v(i,1)**2.d0+v(i,2)**2.d0)
         do j=i+1,n
           Ep=Ep+4.d0*(radio(i,j)**(-12.d0)-radio(i,j)**(-6.d0))
         enddo
      enddo
      Et=Ec+Ep
      write(5,*) t, Ec, Ep, Et


!!!!!!!ALGORITMO PARA CµLCULO MECµNICA  -------------------------------
      do p=1,30000
        !Calculo omega
        omega=v+dt*0.5d0*a1
        !Calculo posici¢n siguiente
        r_aux=r+dt*omega
        !Cond peri¢dicas en posici¢n
        do i=1,n
          do j=1,2
          pdec=(r_aux(i,j)/d)-dfloat(int(r_aux(i,j)/d))
          if (r_aux(i,j).gt.d) then
           r(i,j)=pdec*d
          else if (r_aux(i,j).lt.0.d0) then
           r(i,j)=(pdec+1.d0)*d
          else
           r(i,j)=r_aux(i,j)
          endif
          enddo
        enddo

!!!!!!!Calculo aceleraci¢n siguiente
      do i=1,n
        do j=1,n
           acelx(i,j)=0.d0
           acely(i,j)=0.d0
        enddo
      enddo

      do i=1,n-1
         do j=i+1,n
            !Cond contorno en fuerza
            do q=1,2
               if ((r(j,q)-r(i,q)).gt.(d/2.d0)) then
                  paux(q)=r(j,q)-d
               else if ((r(j,q)-r(i,q)).lt.(-d/2.d0)) then
                  paux(q)=r(j,q)+d
               else
                  paux(q)=r(j,q)
               endif
            enddo
            !Sigo
            rad=((paux(1)-r(i,1))**2.d0+(paux(2)-r(i,2))**2.d0)**0.5d0
            radio(i,j)=rad
            if (3.ge.rad) then
              aux_a=24.d0*(2.d0*rad**(-13.d0)-rad**(-7.d0))
              acelx(i,j)=-aux_a*(paux(1)-r(i,1))/rad
              acely(i,j)=-aux_a*(paux(2)-r(i,2))/rad
            else
              acelx(i,j)=0.d0
              acely(i,j)=0.d0
            end if
            acelx(j,i)=-acelx(i,j)
            acely(j,i)=-acely(i,j)
            radio(j,i)=radio(i,j)
         enddo
      enddo
      !Aceleraci¢n total de cada part¡cula
      do i=1,n
        do q=1,2
           a1(i,q)=0.d0
        enddo
      enddo
      do i=1,n
        do j=1,n
           a1(i,1)=a1(i,1)+acelx(i,j)
           a1(i,2)=a1(i,2)+acely(i,j)
        enddo
      enddo
      
        !Calculo velocidad siguiente
        v=omega+dt*0.5d0*a1
           !Reescalo la velocidad den tiempos indicados
           if((t.le.20).and.(t+dt.ge.20)) then
            v=v*1.5d0
           elseif((t.le.30).and.(t+dt.ge.30)) then
            v=v*1.5d0
           elseif((t.le.35).and.(t+dt.ge.35)) then
            v=v*1.5d0
           elseif((t.le.45).and.(t+dt.ge.45)) then
            v=v*1.5d0
           endif

        !Calculo todo lo demÂ s
        t=t+dt

        !Escribo posiciones en fichero

        write(2,*) ((r(i,j), j=1,2), i=1,n)

        write(2,*)
        write(2,*)
        
        !Escribo fluctuaciones en posici¢n de part¡cula 6
        
        write(9,*) t, (r(6,1)-posini(1))**2.d0+(r(6,2)-posini(2))**2.d0
        
        !Escribo aceleraciones en fichero
        write(4,*) ((a1(i,j), j=1,2), i=1,n)

        !Energ¡as en cada paso
        Ec=0.d0
        Ep=0.d0
        Et=0.d0
                !Energ¡a cin‚tica:
      do i=1,n
         Ec=Ec+0.5d0*(v(i,1)**2.d0+v(i,2)**2.d0)
         do j=i+1,n
           Ep=Ep+4.d0*(radio(i,j)**(-12.d0)-radio(i,j)**(-6.d0))
         enddo
      enddo
      
      Et=Ec+Ep
      write(5,*) t, Ec, Ep, Et
      
         !Temperatura
         if((t.ge.50.d0).and.(t.le.100.d0)) then
           cont=cont+1
           do i=1,n
             Temp=Temp+0.5d0*(v(i,1)**2.d0+v(i,2)**2.d0)/dfloat(n)
             write(7,*) (v(i,1)**2.d0+v(i,2)**2.d0)**0.5d0,v(i,1),v(i,2)
           enddo
         endif
      enddo
      
      Temp=Temp/dfloat(cont)
      write(8,*) "Valor medio de temperatura entre t=50-100s", Temp

      close(2)
      close(4)
      close(5)
      close(6)
      close(7)
      close(8)
      close(9)


      stop
      end

      !include 'dranxor2_new.f'
