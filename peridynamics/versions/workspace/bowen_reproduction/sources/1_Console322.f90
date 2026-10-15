program main
implicit none
integer ndivx,ndivy,totnode,nt,maxfam,nnum,cnode,i,j,tt,nbnd,ntotnode,totint,totright,nnum1,totelement,ccnode,totnode1,totele
parameter(ndivx=200)
parameter(ndivy=100)
parameter(nbnd=3)
parameter(ntotnode=(ndivx+20)*ndivy) 
parameter(nt=10000)
parameter(maxfam=100)
real*8 length,width,thick,dx,delta,dens,emod,pratio,area,vol,bc,emod1,A,aaa,bbb,bc1,C,S
real*8 sedload1,sedload2,dt,totime,ctime,idist,fac,radij,nlength,dforce1,dforce2,pi,theta 
real*8 scr,scx,scy,coordx,coordy,cn,cn1,cn2,pforce11(ntotnode,2),pforce22(ntotnode,2),dispdispdisp(ntotnode)
real*8 coord(ntotnode,2),pforce(ntotnode,2),pforceold(ntotnode,2),bforce(ntotnode,2),stendens(ntotnode,2)
real*8 fncst(ntotnode,2),fncstold(ntotnode,2),disp(ntotnode,2),vel(ntotnode,2),velhalfold(ntotnode,2),velhalf(ntotnode,2)
real*8 acc(ntotnode,2),massvec(ntotnode,2),ske(4,4),coord1(ntotnode,2),wy(4),skewy(4),appres
integer numfam(ntotnode,1),pointfam(ntotnode,1),nodefam(10000000,1),point(ntotnode),nodes(ntotnode,2)
integer hh(ntotnode),m,n,hhh(ntotnode),nodenode(ntotnode,49),nnum2,node(ntotnode,4),cnode1,cnode2,cnode3,cnode4
real*8 coord2(ntotnode,2),aaaa,bbbb,coord3(ntotnode,2,49),coordcoord(ntotnode,2,49),dispdisp(ntotnode,2,49)
real*8 N1(ntotnode,49),N2(ntotnode,49),N3(ntotnode,49),N4(ntotnode,49),pforce33(ntotnode,2)
pi=dacos(-1.0d0)

do i=1,ntotnode 
	coord(i,1)=0.0d0
	coord(i,2)=0.0d0
	numfam(i,1)=0
	pointfam(i,1)=0
	pforce(i,1)=0.0d0
	pforce(i,2)=0.0d0
	pforce11(i,1)=0.0d0
	pforce11(i,2)=0.0d0
	pforce22(i,1)=0.0d0
	pforce22(i,2)=0.0d0
	pforce33(i,1)=0.0d0
	pforce33(i,2)=0.0d0
	pforceold(i,1)=0.0d0
	pforceold(i,2)=0.0d0
	bforce(i,1)=0.0d0
	bforce(i,2)=0.0d0
	stendens(i,1)=0.0d0
	stendens(i,2)=0.0d0
	fncst(i,1)=1.0d0
	fncstold(i,1)=1.0d0  
	fncst(i,2)=1.0d0
	fncstold(i,2)=1.0d0  
	disp(i,1)=0.0d0
	disp(i,2)=0.0d0
	vel(i,1)=0.0d0
	vel(i,2)=0.0d0
	velhalfold(i,1)=0.0d0
	velhalfold(i,2)=0.0d0
	velhalf(i,1)=0.0d0
	velhalf(i,2)=0.0d0
	acc(i,1)=0.0d0
	acc(i,2)=0.0d0
	massvec(i,1)=0.0d0
	massvec(i,2)=0.0d0
    point(i)=0
enddo

do i=1,10000000
	nodefam(i,1)=0
enddo

length=200.0d-3
width=100.0d-3
dx=length/ndivx
thick=dx
delta=3.015d0*dx
dens=7850.0d0
emod=2.5d9
emod1=2.0d11
A=pi*(3.0d-3)**2/4.0d0
pratio=1.0d0/3.0d0
area=dx*dx
vol=area*dx
bc=9.0d0*emod/(pi*thick*(delta**3))
bc1=9.0d0*emod1/(pi*thick*(delta**3))
sedload1=9.0d0/16.0d0*emod*1.0d-6
sedload2=9.0d0/16.0d0*emod*1.0d-6 
dt=1.0d0
totime=nt*dt
ctime=0.0d0
idist=0.0d0
fac=0.0d0
radij=dx/2.0d0
nnum=0
nnum1=0
nnum2=0
cnode=0
nlength=0.0d0
dforce1=0.0d0
dforce2=0.0d0
appres=2.0d7
C=69.0d0/sqrt(7162.0d0)
S=49.0d0/sqrt(7162.0d0)
ske(1,1)=emod1*A/(sqrt(2.0d0)*dx)*C*C
ske(1,2)=emod1*A/(sqrt(2.0d0)*dx)*C*S
ske(1,3)=-1.0d0*emod1*A/(sqrt(2.0d0)*dx)*C*C
ske(1,4)=-1.0d0*emod1*A/(sqrt(2.0d0)*dx)*C*S
ske(2,1)=ske(1,2)
ske(2,2)=emod1*A/(sqrt(2.0d0)*dx)*S*S
ske(2,3)=-1.0d0*emod1*A/(sqrt(2.0d0)*dx)*C*S
ske(2,4)=-1.0d0*emod1*A/(sqrt(2.0d0)*dx)*S*S
ske(3,1)=ske(1,3)
ske(3,2)=ske(2,3)
ske(3,3)=emod1*A/(sqrt(2.0d0)*dx)*C*C
ske(3,4)=emod1*A/(sqrt(2.0d0)*dx)*C*S
ske(4,1)=ske(1,4)
ske(4,2)=ske(2,4)
ske(4,3)=ske(3,4)
ske(4,4)=emod1*A/(sqrt(2.0d0)*dx)*S*S
    
do i=1,201
    do j=1,101
        coordx=-1.0d0/2.0d0*length+(i-1)*dx
        coordy=-1.0d0/2.0d0*width+(j-1)*dx
        nnum=nnum+1
        point(nnum)=0
        coord(nnum,1)=coordx
        coord(nnum,2)=coordy
    enddo
enddo
totint=nnum    

do i=1,3
    do j=1,101
        nnum=nnum+1
        point(nnum)=0
        coord(nnum,1)=1.0d0/2.0d0*length+i*dx
        coord(nnum,2)=-1.0d0/2.0d0*width+(j-1)*dx
    enddo
enddo    
totright=nnum

do i=1,3
    do j=1,101
        nnum=nnum+1
        point(nnum)=0
        coord(nnum,1)=-1.0d0/2.0d0*length-i*dx
        coord(nnum,2)=-1.0d0/2.0d0*width+(j-1)*dx
    enddo
enddo    
totnode1=nnum

do i=1,206
    do j=1,100
        nnum2=nnum2+1
        coord2(nnum2,1)=-1.0d0/2.0d0*length+(i-1)*dx-3.0d0*dx+(dx/2.0d0)
        coord2(nnum2,2)=-1.0d0/2.0d0*width+(j-1)*dx+(dx/2.0d0)
    enddo
enddo
totele=nnum2

do i=1,totele
    do j=1,totnode1
        aaa=coord(j,1)-coord2(i,1)
        bbb=coord(j,2)-coord2(i,2)
        if((abs(aaa+(dx/2.0d0)).le.(1.0d-6)).and.(abs(bbb+(dx/2.0d0)).le.(1.0d-6)))then
            node(i,1)=j
        elseif((abs(aaa-(dx/2.0d0)).le.(1.0d-6)).and.(abs(bbb+(dx/2.0d0)).le.(1.0d-6)))then
            node(i,2)=j
        elseif((abs(aaa-(dx/2.0d0)).le.(1.0d-6)).and.(abs(bbb-(dx/2.0d0)).le.(1.0d-6)))then     
            node(i,3)=j
        elseif((abs(aaa+(dx/2.0d0)).le.(1.0d-6)).and.(abs(bbb-(dx/2.0d0)).le.(1.0d-6)))then
            node(i,4)=j
        endif
    enddo
enddo

do i=1,100
    nnum=nnum+1
    point(nnum)=1
    coord(nnum,1)=-70.0d0*dx+(i-1)*140.0d0/100.0d0*dx+(dx/2.0d0)*140.0d0/100.0d0
    coord(nnum,2)=-50.0d0*dx+(i-1)*dx+(dx/2.0d0)
enddo
totnode=nnum

do i=1,99
    nnum1=nnum1+1
    coord1(nnum1,1)=-70.0d0*dx+i*140.0d0/100.0d0*dx
    coord1(nnum1,2)=-50.0d0*dx+i*dx
enddo
totelement=nnum1

do i=1,totelement
    do j=(totnode1+1),totnode
        aaa=coord(j,1)-coord1(i,1)
        bbb=coord(j,2)-coord1(i,2)
        if((abs(aaa+140.0d0/100.0d0*dx/2.0d0).le.(1.0d-6)).and.(abs(bbb+dx/2.0d0).le.(1.0d-6)))then
            nodes(i,1)=j
        elseif((abs(aaa-140.0d0/100.0d0*dx/2.0d0).le.(1.0d-6)).and.(abs(bbb-dx/2.0d0).le.(1.0d-6)))then 
            nodes(i,2)=j
        endif
    enddo
enddo
    
do i=(totnode1+1),totnode
    hh(i)=0
    do m=1,7
        do n=1,7
            hh(i)=hh(i)+1
            coord3(i,1,hh(i))=coord(i,1)-3.0d0*dx+(m-1)*dx
            coord3(i,2,hh(i))=coord(i,2)-3.0d0*dx+(n-1)*dx
        enddo
    enddo
enddo
          
do i=(totnode1+1),totnode
    hhh(i)=0
    do j=1,hh(i)
        idist=dsqrt((coord3(i,1,j)-coord(i,1))**2+(coord3(i,2,j)-coord(i,2))**2)
        if((idist.le.delta).and.(idist.ge.(1.0d-6)))then
            if((coord3(i,1,j).ge.(-103.0d0*dx)).and.(coord3(i,1,j).le.(103.0d0*dx)))then
                if((coord3(i,2,j).ge.(-50.0d0*dx)).and.(coord3(i,2,j).le.(50.0d0*dx)))then
                    hhh(i)=hhh(i)+1
                    coordcoord(i,1,hhh(i))=coord3(i,1,j)
                    coordcoord(i,2,hhh(i))=coord3(i,2,j)
                endif
            endif
        endif
    enddo
enddo
        
do i=(totnode1+1),totnode
    do j=1,hhh(i)
        do m=1,totele
            cnode=node(m,1)
            ccnode=node(m,3)
            aaa=(coordcoord(i,1,j)-coord(cnode,1))*(coordcoord(i,1,j)-coord(ccnode,1))
            bbb=(coordcoord(i,2,j)-coord(cnode,2))*(coordcoord(i,2,j)-coord(ccnode,2))
            if((aaa<=0.0d0).and.(bbb<=0.0d0))then
                nodenode(i,j)=m
                aaaa=(coordcoord(i,1,j)-coord(cnode,1))-(dx/2.0d0)
                bbbb=(coordcoord(i,2,j)-coord(cnode,2))-(dx/2.0d0)
                N1(i,j)=(1-aaaa/(dx/2.0d0))*(1-bbbb/(dx/2.0d0))/4.0d0
                N2(i,j)=(1+aaaa/(dx/2.0d0))*(1-bbbb/(dx/2.0d0))/4.0d0
                N3(i,j)=(1+aaaa/(dx/2.0d0))*(1+bbbb/(dx/2.0d0))/4.0d0
                N4(i,j)=(1-aaaa/(dx/2.0d0))*(1+bbbb/(dx/2.0d0))/4.0d0
            endif
        enddo
    enddo
enddo
   
do i=1,totnode1
    if(i.eq.1)then 
        pointfam(i,1)=1
    else
        pointfam(i,1)=pointfam(i-1,1)+numfam(i-1,1)
    endif
    do j=1,totnode1
        idist=dsqrt((coord(j,1)-coord(i,1))**2+(coord(j,2)-coord(i,2))**2)
        if(i.ne.j)then
            if(idist.le.delta)then
                numfam(i,1)=numfam(i,1)+1
                nodefam(pointfam(i,1)+numfam(i,1)-1,1)=j
            endif
        endif
    enddo
enddo    
    
do i=1,totnode1
    disp(i,1)=0.001d0*coord(i,1)
    disp(i,2)=0.0d0
enddo
do i=1,totnode1
    stendens(i,1)=0.0d0
    do j=1,numfam(i,1)
        cnode=nodefam(pointfam(i,1)+j-1,1)
        idist=dsqrt((coord(cnode,1)-coord(i,1))**2+(coord(cnode,2)-coord(i,2))**2)
        nlength=dsqrt((coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))**2+(coord(cnode,2)+disp(cnode,2)-coord(i,2)-disp(i,2))**2)
        if(idist.le.delta-radij)then
            fac=1.0d0
        elseif(idist.le.delta+radij)then
            fac=(delta+radij-idist)/(2.0d0*radij)
        else
            fac=0.0d0
        endif           
		stendens(i,1)=stendens(i,1)+0.5d0*0.5d0*bc*((nlength-idist)/idist)**2*idist*vol*fac 
    enddo
    fncst(i,1)=sedload1/stendens(i,1)
enddo
    
do i=1,totnode1
    disp(i,1)=0.0d0
    disp(i,2)=0.001d0*coord(i,2)
enddo
do i=1,totnode1
    stendens(i,2)=0.0d0
    do j=1,numfam(i,1)
        cnode=nodefam(pointfam(i,1)+j-1,1)
        idist=dsqrt((coord(cnode,1)-coord(i,1))**2+(coord(cnode,2)-coord(i,2))**2)
        nlength=dsqrt((coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))**2+(coord(cnode,2)+disp(cnode,2)-coord(i,2)-disp(i,2))**2)
        if(idist.le.delta-radij)then
            fac=1.0d0
        elseif(idist.le.delta+radij)then
            fac=(delta+radij-idist)/(2.0d0*radij)
        else
            fac=0.0d0
        endif        
        stendens(i,2)=stendens(i,2)+0.5d0*0.5d0*bc*((nlength-idist)/idist)**2*idist*vol*fac  
    enddo
    fncst(i,2)=sedload2/stendens(i,2)
enddo    

do i=1,totnode
    vel(i,1)=0.0d0
    disp(i,1)=0.0d0 
    vel(i,2)=0.0d0
    disp(i,2)=0.0d0 
enddo
    
do i=1,totnode
    if(point(i)==0)then
        massvec(i,1)=0.25d0*dt*dt*(pi*(delta)**2*thick)*bc/dx*30.0d0
        massvec(i,2)=0.25d0*dt*dt*(pi*(delta)**2*thick)*bc/dx*30.0d0
    else
        massvec(i,1)=0.25d0*dt*dt*(pi*(delta)**2*thick)*bc1/dx*30.0d0
        massvec(i,2)=0.25d0*dt*dt*(pi*(delta)**2*thick)*bc1/dx*30.0d0   
    endif
enddo        
    
do tt = 1,nt
    write(*,*) "tt=",tt

    do i=(totint+1),totright
        disp(i,1)=(2.0d-7)*tt*dt
        disp(i,2)=0.0d0
    enddo
    do i=(totright+1),totnode1
        disp(i,1)=(-2.0d-7)*tt*dt
        disp(i,2)=0.0d0
    enddo
    
    do i=1,totnode1
        pforce11(i,1)=0.0d0
        pforce11(i,2)=0.0d0
        do j=1,numfam(i,1)            
            cnode=nodefam(pointfam(i,1)+j-1,1)
            idist=dsqrt((coord(cnode,1)-coord(i,1))**2+(coord(cnode,2)-coord(i,2))**2)
            nlength=dsqrt((coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))**2+(coord(cnode,2)+disp(cnode,2)-coord(i,2)-disp(i,2))**2)
            if(idist.le.delta-radij)then
                fac=1.0d0
            elseif(idist.le.delta+radij)then
                fac=(delta+radij-idist)/(2.0d0*radij)
            else
                fac=0.0d0
            endif
            if(dabs(coord(cnode,2)-coord(i,2)).le.1.0d-10)then
                theta=0.0d0
            elseif(dabs(coord(cnode,1)-coord(i,1)).le.1.0d-10)then
                theta=90.0d0*pi/180.0d0
            else
                theta=datan(dabs(coord(cnode,2)-coord(i,2))/dabs(coord(cnode,1)-coord(i,1)))
            endif
            scx=(fncst(i,1)+fncst(cnode,1))/2.0d0
            scy=(fncst(i,2)+fncst(cnode,2))/2.0d0
            scr=1.0d0/(((dcos(theta))**2.0d0/(scx)**2.0d0)+((dsin(theta))**2.0d0/(scy)**2.0d0))
            scr=dsqrt(scr)
            dforce1=bc*(nlength-idist)/idist*vol*scr*fac*(coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))/nlength 
            dforce2=bc*(nlength-idist)/idist*vol*scr*fac*(coord(cnode,2)+disp(cnode,2)-coord(i,2)-disp(i,2))/nlength      
            pforce11(i,1)=pforce11(i,1)+dforce1      
            pforce11(i,2)=pforce11(i,2)+dforce2  
        enddo
    enddo    
    
    do i=1,totelement
        cnode=nodes(i,1)
        ccnode=nodes(i,2)
        pforce22(cnode,1)=0.0d0
        pforce22(cnode,2)=0.0d0
        pforce22(ccnode,1)=0.0d0
        pforce22(ccnode,2)=0.0d0
    enddo
    
    do i=1,totelement
        cnode=nodes(i,1)
        ccnode=nodes(i,2)
        wy(1)=disp(cnode,1)
        wy(2)=disp(cnode,2)
        wy(3)=disp(ccnode,1)
        wy(4)=disp(ccnode,2)
        do m=1,4
            skewy(m)=0.0d0
            do n=1,4
                skewy(m)=skewy(m)+ske(m,n)*wy(n)
            enddo
        enddo
        pforce22(cnode,1)=pforce22(cnode,1)+skewy(1)*(-1.0d0)/vol
        pforce22(cnode,2)=pforce22(cnode,2)+skewy(2)*(-1.0d0)/vol
        pforce22(ccnode,1)=pforce22(ccnode,1)+skewy(3)*(-1.0d0)/vol
        pforce22(ccnode,2)=pforce22(ccnode,2)+skewy(4)*(-1.0d0)/vol
    enddo    
    
    do i=1,totnode1
        pforce33(i,1)=0.0d0
        pforce33(i,2)=0.0d0
    enddo
    
    do i=(totnode1+1),totnode
        pforce33(i,1)=0.0d0
        pforce33(i,2)=0.0d0
        do j=1,hhh(i)
            cnode=nodenode(i,j)
            cnode1=node(cnode,1)
            cnode2=node(cnode,2)
            cnode3=node(cnode,3)
            cnode4=node(cnode,4)
            dispdisp(i,1,j)=N1(i,j)*disp(cnode1,1)+N2(i,j)*disp(cnode2,1)+N3(i,j)*disp(cnode3,1)+N4(i,j)*disp(cnode4,1)
            dispdisp(i,2,j)=N1(i,j)*disp(cnode1,2)+N2(i,j)*disp(cnode2,2)+N3(i,j)*disp(cnode3,2)+N4(i,j)*disp(cnode4,2)
            idist=dsqrt((coordcoord(i,1,j)-coord(i,1))**2+(coordcoord(i,2,j)-coord(i,2))**2)
            nlength=dsqrt((coordcoord(i,1,j)+dispdisp(i,1,j)-coord(i,1)-disp(i,1))**2+(coordcoord(i,2,j)+dispdisp(i,2,j)-coord(i,2)-disp(i,2))**2)
            dforce1=bc*(nlength-idist)/idist*vol*(coordcoord(i,1,j)+dispdisp(i,1,j)-coord(i,1)-disp(i,1))/nlength 
            dforce2=bc*(nlength-idist)/idist*vol*(coordcoord(i,2,j)+dispdisp(i,2,j)-coord(i,2)-disp(i,2))/nlength
            pforce33(i,1)=pforce33(i,1)+dforce1    
            pforce33(i,2)=pforce33(i,2)+dforce2     
            pforce33(cnode1,1)=pforce33(cnode1,1)+N1(i,j)*(-1.0d0)*dforce1
            pforce33(cnode1,2)=pforce33(cnode1,2)+N1(i,j)*(-1.0d0)*dforce2
            pforce33(cnode2,1)=pforce33(cnode2,1)+N2(i,j)*(-1.0d0)*dforce1
            pforce33(cnode2,2)=pforce33(cnode2,2)+N2(i,j)*(-1.0d0)*dforce2
            pforce33(cnode3,1)=pforce33(cnode3,1)+N3(i,j)*(-1.0d0)*dforce1
            pforce33(cnode3,2)=pforce33(cnode3,2)+N3(i,j)*(-1.0d0)*dforce2
            pforce33(cnode4,1)=pforce33(cnode4,1)+N4(i,j)*(-1.0d0)*dforce1
            pforce33(cnode4,2)=pforce33(cnode4,2)+N4(i,j)*(-1.0d0)*dforce2
        enddo
    enddo
            
    do i=1,totnode
        pforce(i,1)=pforce11(i,1)+pforce22(i,1)+pforce33(i,1)
        pforce(i,2)=pforce11(i,2)+pforce22(i,2)+pforce33(i,2)
    enddo    

    cn=0.0d0
	cn1=0.0d0
	cn2=0.0d0
	do i=1,totnode
        if(velhalfold(i,1).ne.0.0d0)then
            cn1=cn1-disp(i,1)*disp(i,1)*((pforce(i,1)+bforce(i,1))/massvec(i,1)-(pforceold(i,1)+bforce(i,1))/massvec(i,1))/(dt*velhalfold(i,1))
        endif
        if(velhalfold(i,2).ne.0.0d0)then
            cn1=cn1-disp(i,2)*disp(i,2)*((pforce(i,2)+bforce(i,2))/massvec(i,2)-(pforceold(i,2)+bforce(i,2))/massvec(i,2))/(dt*velhalfold(i,2))
        endif
		cn2=cn2+disp(i,1)*disp(i,1)
		cn2=cn2+disp(i,2)*disp(i,2)
    enddo

	if(cn2.ne.0.0d0)then
        if((cn1/cn2)>0.0d0)then 
            cn=2.0d0*dsqrt(cn1/cn2)
        else
            cn=0.0d0
        endif
    else
        cn=0.0d0
    endif
	if(cn>2.0d0)then
		cn=1.9d0
	endif

	do i=1,totnode
		if(tt.eq.1)then
            velhalf(i,1)=1.0d0*dt/massvec(i,1)*(pforce(i,1)+bforce(i,1))/2.0d0	
            velhalf(i,2)=1.0d0*dt/massvec(i,2)*(pforce(i,2)+bforce(i,2))/2.0d0
        else	
            velhalf(i,1)=((2.0d0-cn*dt)*velhalfold(i,1)+2.0d0*dt/massvec(i,1)*(pforce(i,1)+bforce(i,1)))/(2.0d0+cn*dt)
            velhalf(i,2)=((2.0d0-cn*dt)*velhalfold(i,2)+2.0d0*dt/massvec(i,2)*(pforce(i,2)+bforce(i,2)))/(2.0d0+cn*dt)
        endif 
        vel(i,1)=0.5d0*(velhalfold(i,1)+velhalf(i,1))	
        vel(i,2)=0.5d0*(velhalfold(i,2)+velhalf(i,2))	
        disp(i,1)=disp(i,1)+velhalf(i,1)*dt	
        disp(i,2)=disp(i,2)+velhalf(i,2)*dt	   
        velhalfold(i,1)=velhalf(i,1)
        velhalfold(i,2)=velhalf(i,2)
		pforceold(i,1)=pforce(i,1)
		pforceold(i,2)=pforce(i,2)
    enddo    
    
    do i=1,totnode
        dispdispdisp(i)=dsqrt(disp(i,1)**2+disp(i,2)**2)
    enddo
    
    if(tt.eq.nt)then
        open(26,file="coord_disp_pd_10000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),disp(i,1),disp(i,2),dispdispdisp(i)
            enddo
        close(26)
        open(26,file="coord_10000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),disp(i,1),disp(i,2),dispdispdisp(i)
                endif
            enddo
        close(26)
        open(26,file="xx.txt")
            do i=1,totint
                if(abs(coord(i,2))<(1.0d-6))then
                    write(26,222) coord(i,1),coord(i,2),disp(i,1),disp(i,2),dispdispdisp(i)
                endif
            enddo
        close(26)
        open(26,file="yy.txt")
            do i=1,totint
                if(abs(coord(i,1))<(1.0d-6))then
                    write(26,222) coord(i,1),coord(i,2),disp(i,1),disp(i,2),dispdispdisp(i)
                endif
            enddo
        close(26)
    elseif(tt.eq.8000)then
        open(26,file="coord_disp_pd_8000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),disp(i,1),disp(i,2),dispdispdisp(i)
            enddo
        close(26)
        open(26,file="coord_8000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),disp(i,1),disp(i,2),dispdispdisp(i)
                endif
            enddo
        close(26)
    elseif(tt.eq.5000)then
        open(26,file="coord_disp_pd_5000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),disp(i,1),disp(i,2),dispdispdisp(i)
            enddo
        close(26)
        open(26,file="coord_5000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),disp(i,1),disp(i,2),dispdispdisp(i)
                endif
            enddo
        close(26)
    elseif(tt.eq.4000)then
        open(26,file="coord_disp_pd_4000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),disp(i,1),disp(i,2),dispdispdisp(i)
            enddo
        close(26)
        open(26,file="coord_4000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),disp(i,1),disp(i,2),dispdispdisp(i)
                endif
            enddo
        close(26)
    elseif(tt.eq.2000)then
        open(26,file="coord_disp_pd_2000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),disp(i,1),disp(i,2),dispdispdisp(i)
            enddo
        close(26)
        open(26,file="coord_2000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),disp(i,1),disp(i,2),dispdispdisp(i)
                endif
            enddo
        close(26)
    elseif(tt.eq.6000)then
        open(26,file="coord_disp_pd_6000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),disp(i,1),disp(i,2),dispdispdisp(i)
            enddo
        close(26)
        open(26,file="coord_6000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),disp(i,1),disp(i,2),dispdispdisp(i)
                endif
            enddo
        close(26)
    endif
enddo 
222 format(e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5)
end program main