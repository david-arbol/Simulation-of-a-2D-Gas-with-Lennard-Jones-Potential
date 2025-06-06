      ! Dynamics and Temperature starting with 80 particles in random positions, ensuring that they are not too close to each other.
      program argon
      implicit none
      
      real*8, dimension(:,:), allocatable :: r
      real*8, dimension(:,:), allocatable :: r_aux
      real*8, dimension(:,:), allocatable :: v
      real*8, dimension(:,:), allocatable :: acelx
      real*8, dimension(:,:), allocatable :: acely
      real*8, dimension(:,:), allocatable :: a1
      real*8, dimension(:,:), allocatable :: omega
      real*8  d, aux_a, rad, dt, t, dran_u, ang_aux
      integer i,j,n,p,q
      real*8 pdec
      real*8 paux(2)
      !Variables para cculos postenriores
      real*8 Ec,Ep,Et,Temp
      integer cont
      real*8, dimension(:,:), allocatable :: radio
      logical lejos
      
      Temp=0.d0
      cont=0
      t=0.d0
      lejos=.FALSE.
      
      call dran_ini(783247)
      
      dt=0.002
      n=80

      allocate(r(n,2))
      allocate(r_aux(n,2))
      allocate(v(n,2))
      allocate(acelx(n,n))
      allocate(acely(n,n))
      allocate(a1(n,2))
      allocate(omega(n,2))
      allocate(radio(n,n))
      

      
      d=20.d0 !Dimensiones (dxd) de mi cristal bidimensional
              !Debe ser un nσero natural
              

      !Abro archivos
      
      open(4,file='acel.txt')
      open(2,file='argon.txt')
      
      open(5,file='energia.txt')
      open(1,file='prueba.txt')
      
      open(6,file='vel_t0.txt')
      open(7,file='vel_t1_1.txt')
      
      open(8,file='Temperatura.txt')

      !Inicializo a cero la aceleraciones para luego sumarlas
      do i=1,n
         do j=1,n
          acelx(i,j)=0.d0  !Aceleraci▋ de i ejercida por j
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
      do i=1,n
        lejos=.FALSE.
        DOWHILE (lejos.eqv..FALSE.)
        lejos=.TRUE.
        do j=1,2
           r(i,j)=dran_u()*d
        enddo
        do p=1,i-1
          do q=1,2
            if ((r(p,q)-r(i,q)).gt.(d/2.d0)) then
               paux(q)=r(p,q)-d
            else if ((r(p,q)-r(i,q)).lt.(-d/2.d0)) then
               paux(q)=r(p,q)+d
            else
               paux(q)=r(p,q)
            endif
          enddo
       radio(i,p)=((paux(1)-r(i,1))**2.d0+(paux(2)-r(i,2))**2.d0)**0.5d0
          if (radio(i,p).le.1.1d0) then
             lejos=lejos.and..FALSE.
          endif
        enddo
        ENDDO
      enddo

      write(2,*) ((r(i,j), j=1,2), i=1,n)
      write(2,*)
      write(2,*)

      !Inicializo las velocidades iniciales a modulo 1 y gulo random

      do i=1,n
           ang_aux=dran_u()*4.d0*datan(1.d0)*2.d0
           v(i,1)=4.d0*cos(ang_aux)
           v(i,2)=4.d0*sin(ang_aux)
      enddo
      
      !do i=1,n
      !     v(i,1)=dran_u()
      !     v(i,2)=0.d0
      !enddo
      
      
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
      !Aceleraci▋ total de cada part︷ula
      do i=1,n
        do j=1,n
           a1(i,1)=a1(i,1)+acelx(i,j)
           a1(i,2)=a1(i,2)+acely(i,j)
        enddo
      enddo

      
      !ENERG淲S:
                !Inicializo
      Ec=0.d0
      Ep=0.d0
      Et=0.d0
                !Energ｛ cinica:
      do i=1,n
         Ec=Ec+0.5d0*(v(i,1)**2.d0+v(i,2)**2.d0)
         do j=i+1,n
           Ep=Ep+4.d0*(radio(i,j)**(-12.d0)-radio(i,j)**(-6.d0))
         enddo
      enddo
      Et=Ec+Ep
      write(5,*) t, Ec, Ep, Et


!!!!!!!ALGORITMO PARA C無CULO MEC煮ICA  -------------------------------
      do p=1,75000
        !Calculo omega
        omega=v+dt*0.5d0*a1
        !Calculo posici▋ siguiente
        r_aux=r+dt*omega
        !Cond peri▃icas en posici▋
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

!!!!!!!Calculo aceleraci▋ siguiente
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
      !Aceleraci▋ total de cada part︷ula
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

        !Calculo todo lo dem�s
        t=t+dt

        !Escribo posiciones en fichero

        write(2,*) ((r(i,j), j=1,2), i=1,n)

        write(2,*)
        write(2,*)

        !Escribo aceleraciones en fichero
        write(4,*) ((a1(i,j), j=1,2), i=1,n)

        !Energ｛s en cada paso
        Ec=0.d0
        Ep=0.d0
        Et=0.d0
                !Energ｛ cinica:
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


      stop
      end


      include 'dranxor2_new.f'
