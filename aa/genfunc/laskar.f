	program laskar
	implicit none
	integer npt,nfrq,maxline,nline,kplot,kount,i,j,n,kstop,ithet,nthet,iphi,nphi
	parameter(npt=2048,nfrq=npt/2,maxline=5,nthet=10,nphi=10)
	real*4 pi,torad,y(6),dydx(6),yscal(6),f(3),t,t0,eps,h0,hdid,hnext,
     >		deltat,e1,e2,theta,phi,tmax,val
	parameter (pi=3.14159,torad=pi/180.)
	real*4 q1,q2,rc2,xout(3,npt),vout(3,npt),xx(npt),yy(npt),zz(npt)
	common/Psidat/q1,q2,rc2
	external derivs,rfun
	real*4 data(2*nfrq+2),cfreq(4,maxline),freqs(3)
	real*4 rfun,r,ej0,sthet,cthet,sphi,cphi,Psi
	common/angles/ej0,sthet,cthet,sphi,cphi

	data yscal/6*1.e0/,h0/1.e0/,kplot/10/,eps/1.e-7/
c
	open(13,file='laskar.dat',status='unknown')

	rc2=0.1
	q1=.9
	q2=.7
	ej0=0.5
c	write(*,*)' Enter theta, phi, tmax '
c	read(*,*) theta,phi,tmax
c	theta=theta*torad
c	phi=phi*torad
	tmax=500.
	deltat=tmax/float(2*nfrq-1)

	do ithet=1,nthet
c		ithet=1
c		iphi=1
		theta=acos(ithet/float(nthet+1))
		print*,ithet
		do iphi=1,nphi
			phi=.5*pi*iphi/float(nphi+1)
			sthet=sin(theta)
			cthet=cos(theta)
			sphi=sin(phi)
			cphi=cos(phi)
			call root(.001,3.,1.e-4,rfun,val,r)

			y(1)=r*sthet*cphi
			y(2)=r*sthet*sphi
			y(3)=r*cthet
			y(4)=0.
			y(5)=0.
			y(6)=0.

			t=0.
			t0=t
			kount=1
			xout(1,1)=y(1)
			xout(2,1)=y(2)
			xout(3,1)=y(3)
			call energy(y,e1)
c			write(*,'(a,3f8.3,i6)')' Energy:', e1,ej0
			call derivs(t,y,dydx)
20			continue
			call bsstep(y,dydx,6,t,h0,eps,yscal,hdid,hnext,derivs)
			call derivs(t,y,dydx)
			call nfillup(y,dydx,6,t,hdid,xout,vout,3,npt,kount,t0,deltat)
			h0=hnext

			if(t.lt.tmax.and.kount.lt.npt)  go to 20
c			write(*,'(4f8.3,i6)') t,y(1),y(2),y(3),kount

			call energy(y,e2)
			e2=e2-e1
			write(*,'(a,3(1pe11.3))')' E0, E, dE =', ej0,e1,e2
			do n=1,kount
				xx(n)=xout(1,n)
				yy(n)=xout(2,n)
				zz(n)=xout(3,n)
			enddo
c			call plots(kount,xx,zz,-2.,2.,-2.,2.,'x',1,'y',1,-1.01,kplot)
			kplot=1
			call realft(xx,nfrq,1)
			call newsperg(xx,deltat,nfrq,cfreq,maxline,nline)
			freqs(1)=cfreq(1,1)
			if(freqs(1).eq.0.) freqs(1)=cfreq(1,2)
c			write(*,*)'     Freq        Amp      Phase    Tres/Tres0'
c			do i=1,nline
c				write(*,'(4f11.4)') (cfreq(j,i),j=1,4)
c			enddo
			call realft(yy,nfrq,1)
			call newsperg(yy,deltat,nfrq,cfreq,maxline,nline)
			freqs(2)=cfreq(1,1)
			if(freqs(2).eq.0.) freqs(2)=cfreq(1,2)
c			write(*,*)'     Freq        Amp      Phase    Tres/Tres0'
c			do i=1,nline
c				write(*,'(4f11.4)') (cfreq(j,i),j=1,4)
c			enddo
			call realft(zz,nfrq,1)
			call newsperg(zz,deltat,nfrq,cfreq,maxline,nline)
			freqs(3)=cfreq(1,1)
			if(freqs(3).eq.0.) freqs(3)=cfreq(1,2)
c			write(*,*)'     Freq        Amp      Phase    Tres/Tres0'
c			do i=1,nline
c				write(*,'(4f11.4)') (cfreq(j,i),j=1,4)
c			enddo
c			call grend(0)
			call compress(13,freqs,3)
		enddo
	enddo
	stop
	end

	function rfun(r)
	implicit none
	real*4 rfun,r,ej0,sthet,cthet,sphi,cphi,Psi
	common/angles/ej0,sthet,cthet,sphi,cphi
	rfun=ej0-Psi(r*sthet*cphi,r*sthet*sphi,r*cthet)
	return
	end

	subroutine force(y,f)
	implicit none
	real*4 q1,q2,rc2,y(6),f(3),r2
	common/Psidat/q1,q2,rc2
	r2=rc2+y(1)**2 +(y(2)/q1)**2+(y(3)/q2)**2
	f(1)=-y(1)/r2
	f(2)=-y(2)/(q1*q1*r2)
	f(3)=-y(3)/(q2*q2*r2)
	return
	end

	function Psi(x,y,z)
	implicit none
	real*4 Psi,x,y,z
	real*4 q1,q2,rc2,r2
	common/Psidat/q1,q2,rc2
	r2=rc2+x**2+(y/q1)**2+(z/q2)**2
	Psi=.5*log(r2)
	return
	end

	subroutine derivs(t,y,dydt)
	implicit real*4 (a-h,o-z)
	dimension y(6),dydt(6),f(3)
	call force(y,f)
	dydt(1)=y(4)
	dydt(2)=y(5)
	dydt(3)=y(6)
	dydt(4)=f(1)
	dydt(5)=f(2)
	dydt(6)=f(3)
	return
	end

	subroutine energy(y,ej)
	implicit none
	real*4 y(6),ej,Psi
	ej=Psi(y(1),y(2),y(3))+.5*(y(4)**2+y(5)**2+y(6)**2)
	return
	end

	subroutine newsperg(data,deltat,nfrq,cfreq,maxline,nline)
	parameter(pi=3.1415926,tpi=2.*pi,acc=.001)
	real*4 cfreq(4,maxline)
	complex data(nfrq+1),z(2049),zs,sp,s,sm

	data(1)=cmplx(real(data(1)),0.)
	data(nfrq+1)=cmplx(real(data(nfrq+1)),0.)
	tpion=2.*pi/float(2*nfrq-1)
	prefac=pi/float(2*nfrq-1)
	nline=0

    1	continue
c							Calc 2nd diff
	do i=2,nfrq
		z(i)=prefac*(data(i+1)+data(i-1)-2.*data(i))
	enddo

   	z(1)=2.*prefac*(data(2)-data(1))

c							Locate peak
	tres=0.
	biggest=0.
	do i=1,nfrq-1
		this=z(i)*conjg(z(i))
		tres=tres+this
		if(this.gt.biggest) then
			ibig=i
			biggest=this
		endif
	enddo
   	tres=sqrt(tres)
	if(nline.eq.0) then
		tres0=tres
	else
		cfreq(4,nline)=tres/tres0
		if(tres/tres0.lt.acc.or.nline.ge.maxline) return
	endif
    	nline=nline+1
c						Find frequency
	if(ibig.eq.1) then
		freq=0.
		zs=-z(1)/tpi
		amp=real(zs)
		phase=0.
	else
		sp=z(ibig+1)
		s=z(ibig)
		sm=z(ibig-1)

		xm=-real((s+2.*sm)/(s-sm))
		xp=real((s+2.*sp)/(s-sp))
		x=.5*(xp+xm)

		freq=(float(ibig-1)-x)*tpion/deltat
	
c						Find amplitude
		zs=z(ibig)*.5*(x+1.)*x*(x-1.)
		if(abs(x).lt.0.03) then
		s=z(ibig)*.5*(x+1.)*(x-1.)/pi
		else
			s=z(ibig)*.5*(x+1.)*x*(x-1.)/sin(pi*x)
		endif

		amp=sqrt(s*conjg(s))
		phase=-atan2(aimag(s),real(s))+pi*x
	endif
	cfreq(1,nline)=freq
	cfreq(2,nline)=amp
	cfreq(3,nline)=phase

	call subtract(data,nfrq,deltat,freq,zs)
c	call nsubtract(data,nfrq,freq,zalpha,x,ibig,zs)
	go to 1

	end

	subroutine subtract(data,nfrq,deltat,freq,zs)
	parameter(pi=3.1415926,tpi=2.*pi)
	complex data(nfrq+1),zs

	T=float(2*nfrq)*deltat

	if(freq.gt.0.)  then
		do i=1,nfrq
			domeg=float(i-1)*tpi/T-freq
			sins=-sin(.5*domeg*deltat)
			data(i)=data(i)-zs/sins
		enddo
   	else
		data(1)=data(1)-zs*float(2*nfrq)
	endif

   	return
	end

	subroutine nsubtract(data,nfrq,freq,zalpha,x0,ibig,zs0)
	implicit none
	integer nfrq,ibig,i
	real*4 pi,freq,x0,x,sins
	parameter(pi=3.1415926)
	complex data(nfrq+1),zs0,zs,zalpha

	if(freq.gt.0.)  then
		zs=zs0
		do i=ibig,nfrq
			x=x0-float(i-ibig)
			sins=sin(pi*x/float(2*nfrq))
			data(i)=data(i)-zs/sins
			zs=-zs/zalpha
		enddo
		zs=-zs0*zalpha
		do i=ibig-1,1,-1
			x=x0-float(i-ibig)
			sins=sin(pi*x/float(2*nfrq))
			data(i)=data(i)-zs/sins
			zs=-zs*zalpha
		enddo
   	else
		data(1)=data(1)-zs0*float(2*nfrq)
	endif

   	return
	end

	subroutine nfillup(y2,yd2,nv,t2,hdid,xout,vout,nc,npt,kount,t0,deltat)
	implicit real*4 (a-h,o-z)
	common/bbs1/y1(10),yd1(10)
	dimension y2(nv),yd2(nv)
	real*4 xout(nc,npt),vout(nc,npt)
	save
c
	t1=t2-hdid
	t3=0.
	if(hdid.le.0.) stop' hdid <=0 in fillup'
	do 20  i=1,nc
		a0=(y2(i)-(y1(i)+hdid*(y1(nc+i)+hdid*.5*yd1(nc+i))))/hdid**3
		a1=(yd2(i)-(yd1(i)+hdid*(yd1(nc+i)+hdid*3.*a0)))/hdid**3
		a2=(yd2(nc+i)-(yd1(nc+i)+6.*hdid*(a0+hdid*a1)))/hdid**3
		t=t0
		k=kount
5		t=t+deltat
		if(t.lt.t2)  then
			if(k.eq.npt) go to 20
			k=k+1
			d1=t-t1
			d2=t-t2
			xout(i,k)=y1(i)+d1*(yd1(i)+d1*(.5*yd1(nc+i)+d1*(a0+d2*(a1+.5*d2*a2))))
     		vout(i,k)=yd1(i)+d1*(yd1(nc+i)+d1*(3.*(a0+d2*(a1+.5*d2*a2)))+d1*(a1+d2*a2))
   		go to 5
	endif
20	continue
	kount=k
	t0=t-deltat
	return
	end

	include '\\u\\f\\spress\\fourier.f'
	include '\\u\\f\\spress\\realft.f'