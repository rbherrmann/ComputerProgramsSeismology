	uaz=270
	uth=54.7
	vaz=30
	vth=54.7
	waz=150
	wth=54.7
        dr=3.1415927/180.
        sr6=sqrt(6.)
c  X=E
	a21= sr6*sin(uaz*dr)*sin(uth*dr)
	a22= sr6*sin(vaz*dr)*sin(vth*dr)
	a23= sr6*sin(waz*dr)*sin(wth*dr)
        write(6,*)a21,a22,a23
c  Y=N
	a11= sr6*cos(uaz*dr)*sin(uth*dr)
	a12= sr6*cos(vaz*dr)*sin(vth*dr)
	a13= sr6*cos(waz*dr)*sin(wth*dr)
        write(6,*)a11,a12,a13
c  Z=Z
        a31 = sr6*cos(uth*dr)
        a32 = sr6*cos(vth*dr)
        a33 = sr6*cos(wth*dr)
	write(6,*)a31,a32,a33
	end
