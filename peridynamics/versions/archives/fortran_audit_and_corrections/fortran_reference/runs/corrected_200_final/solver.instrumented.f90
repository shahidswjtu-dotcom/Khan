program main
implicit none
integer first_yield_step, yielded_rods
real*8 max_us,max_uc,max_vs,max_vc,net_internal(3),net_interface(3),node_volume
integer i,j,nnum,totnode1,nnum1,totelement11,totelement1,ntotnode,k,m,totnode2,totelement2,totnode3,totelement,totelement3
integer totint,totele,nnum2,tottop,totbottom,totnode,cnode,tt,nt,ccnode,n,cnode1,cnode2,cnode3,cnode4,cnode5,cnode6,cnode7,cnode8
parameter(ntotnode=52000) 
parameter(nt=100000)
integer hh(ntotnode),nodes(ntotnode,2),point(ntotnode),hhh(ntotnode),node(ntotnode,8),nodenode(ntotnode,350)
integer numfam(ntotnode,1),pointfam(ntotnode,1),nodefam(20000000,1)
real*8 coord(ntotnode,3),dx,coord1(ntotnode,3),aaa,bbb,ccc,coord3(ntotnode,3,350),coordcoordcoord(ntotnode,3,200),dx1
real*8 coordcoord(ntotnode,3),aaaa,bbbb,cccc,N1(ntotnode,200),N2(ntotnode,200),N3(ntotnode,200),N4(ntotnode,200)
real*8 N5(ntotnode,200),N6(ntotnode,200),N7(ntotnode,200),N8(ntotnode,200),disp(ntotnode,3),stendens(ntotnode,3)
real*8 idist,nlength,fncst(ntotnode,3),sedload1,sedload2,sedload3,fac,radij,delta,vel(ntotnode,3),massvec(ntotnode,3)
real*8 bc,bc1,thick,dt,pi,theta,phi,KE1(2,2),KE2(2,2),KE3(2,2),KE4(2,2),emod,emod1,A1,A2,A3,A4,L1,L2,L3,L4,velhalfold(ntotnode,3),velhalf(ntotnode,3)
real*8 pforceold(ntotnode,3),bforce(ntotnode,3),vol,dforce1,dforce2,dforce3,scr,scx,scy,scz,cn,cn1,cn2,wy(6),skewy(6),scr0,scr1
real*8 pforce11(ntotnode,3),pforce22(ntotnode,3),pforce33(ntotnode,3),pforce(ntotnode,3),dispdisp(ntotnode,3,200),dispdispdisp(ntotnode)
real*8 ske1(6,6),ske2(6,6),ske3(6,6),ske4(6,6),pratio,area,dmg(ntotnode,1),dmgpar1,dmgpar2,fail(ntotnode,200),failfail(ntotnode,200)
real*8 T1(2,6),T1T(6,2),T1TKE(6,2),T2(2,6),T2T(6,2),T2TKE(6,2),T3(2,6),T3T(6,2),T3TKE(6,2),T4(2,6),T4T(6,2),T4TKE(6,2)
real*8 s0,sc,sy,sssss,loadloadx,loadloady,loadloadz,loadxx(100000),loadyy(100000),loadzz(100000),loadloadzong(100000),dispdispld(100000)
real*8 ddt,pendiam,penrad,penccx,penccy,penccz,penicx,penicy,penicz,ppenicx,ppenicy,ppenicz,penalx,penaly,penalz,penux,penuy,penuz
real*8 pentx,penty,pentz,pentfx,pentfy,pentfz,newpx,newpy,newpz,dpenx,dpeny,dpenz,dpend,dpux,dpuy,dpuz,dprx,dpry,dprz,dpcx,dpcy,dpcz
real*8 olddispx,olddispy,olddispz,oldvelx,oldvely,oldvelz,dpenfx,dpenfy,dpenfz,olddisp(ntotnode,3),dens,ppenccx,ppenccy,ppenccz,vol1
real*8 ppenux,ppenuy,ppenuz,pnewpx,pnewpy,pnewpz,pdpenx,pdpeny,pdpenz,pdpend,forceld(100000,3),aaaaaa,bbbbbb,dengs,acc(ntotnode,3),oldvel(ntotnode,3)
pi=dacos(-1.0d0)

do i=1,ntotnode 
    coord(i,1)=0.0d0
    coord(i,2)=0.0d0
    coord(i,3)=0.0d0
    coord1(i,1)=0.0d0
    coord1(i,2)=0.0d0
    coord1(i,3)=0.0d0
    pforce(i,1)=0.0d0
    pforce(i,2)=0.0d0
    pforce(i,3)=0.0d0
    pforce11(i,1)=0.0d0
    pforce11(i,2)=0.0d0
    pforce11(i,3)=0.0d0
    pforce22(i,1)=0.0d0
    pforce22(i,2)=0.0d0
    pforce22(i,3)=0.0d0
    pforce33(i,1)=0.0d0
    pforce33(i,2)=0.0d0
    pforce33(i,3)=0.0d0
    pforceold(i,1)=0.0d0
    pforceold(i,2)=0.0d0
    pforceold(i,3)=0.0d0
    bforce(i,1)=0.0d0
    bforce(i,2)=0.0d0
    bforce(i,3)=0.0d0
    stendens(i,1)=0.0d0
    stendens(i,2)=0.0d0
    stendens(i,3)=0.0d0
    fncst(i,1)=1.0d0 
    fncst(i,2)=1.0d0
    fncst(i,3)=1.0d0
    disp(i,1)=0.0d0
    disp(i,2)=0.0d0
    disp(i,3)=0.0d0
    vel(i,1)=0.0d0
    vel(i,2)=0.0d0
    vel(i,3)=0.0d0
    velhalfold(i,1)=0.0d0
    velhalfold(i,2)=0.0d0
    velhalfold(i,3)=0.0d0
    velhalf(i,1)=0.0d0
    velhalf(i,2)=0.0d0
    velhalf(i,3)=0.0d0
    massvec(i,1)=0.0d0
    massvec(i,2)=0.0d0
    massvec(i,3)=0.0d0
    point(i)=0
    numfam(i,1)=0
    pointfam(i,1)=0
    do j=1,200
        fail(i,j)=1.0
        failfail(i,j)=1.0
    enddo
    dmg(i,1)=0.0d0
enddo

!scr0=-0.03d0*100000000000.0d0
s0=5.00456d-5
sc=1.25114d-3
sy=2.0d-3
dx=19.0d-3
dx1=19.0d-3
nnum=0
nnum1=0
nnum2=0
totnode=0
totnode1=0
delta=3.015d0*dx
emod=3.05d10
emod1=2.08d11
pratio=1.0d0/4.0d0
area=dx*dx
vol=area*dx
bc=12.0d0*emod/(pi*(delta**4))
bc1=12.0d0*emod1/(pi*(delta**4))
sedload1=0.6d0*emod*1.0d-6    
sedload2=0.6d0*emod*1.0d-6
sedload3=0.6d0*emod*1.0d-6
idist=0.0d0
fac=0.0d0
radij=dx/2.0d0
cnode=0
ccnode=0
nlength=0.0d0
dforce1=0.0d0
dforce2=0.0d0
dforce3=0.0d0
A1=pi*(26.0d-3)**2/4.0d0
A2=pi*(10.0d-3)**2/4.0d0
A3=pi*(8.0d-3)**2/4.0d0
A4=pi*(8.0d-3)**2/4.0d0
L1=dx
L2=dx
L3=10.0d-3
L4=10.0d-3
dt=3.0d-6
do i=1,20000000
    nodefam(i,1)=0
enddo
    

dens=2500.0d0
dengs=7850.0d0
pendiam=0.01d0 
penrad=pendiam/2.0d0
!penmass=1.57d0
!penvx=0.0d0 
!penvy=-32.0d0 
!penvz=0.0d0 
penccx=0.0d0 
penccy=0.0d0 
penccz=0.0d0
penicx=-10.5d0*dx
penicy=0.0d0
penicz=9.0d0*dx+penrad-0.999999999999999999999999999d0*dx/2.0d0  

ppenicx=10.5d0*dx
ppenicy=0.0d0
ppenicz=9.0d0*dx+penrad-0.999999999999999999999999999d0*dx/2.0d0  

!penax=0.0d0 
!penay=0.0d0 
!penaz=0.0d0 
penalx=0.0d0 
penaly=0.0d0 
penalz=0.0d0 
penux=0.0d0 
penuy=0.0d0 
penuz=0.0d0 

ppenux=0.0d0 
ppenuy=0.0d0 
ppenuz=0.0d0 

!penwx=0.0d0 
!penwy=0.0d0 
!penwz=0.0d0 
pentx=0.0d0 
penty=0.0d0 
pentz=0.0d0
pentfx=0.0d0
pentfy=0.0d0 
pentfz=0.0d0 
!pentmx=0.0d0 
!pentmy=0.0d0 
!pentmz=0.0d0 
newpx=0.0d0 
newpy=0.0d0 
newpz=0.0d0 
dpenx=0.0d0 
dpeny=0.0d0 
dpenz=0.0d0 
dpend=0.0d0
dpux=0.0d0 
dpuy=0.0d0 
dpuz=0.0d0 
dprx=0.0d0 
dpry=0.0d0 
dprz=0.0d0 
dpcx=0.0d0 
dpcy=0.0d0 
dpcz=0.0d0 
olddispx=0.0d0 
olddispy=0.0d0 
olddispz=0.0d0 
oldvelx=0.0d0 
oldvely=0.0d0 
oldvelz=0.0d0
dpenfx=0.0d0 
dpenfy=0.0d0 
dpenfz=0.0d0 
!dpenmx=0.0d0 
!dpenmy=0.0d0 
!dpenmz=0.0d0 
ddt=3.0d-6
vol1=dx*(530.66d-6)
cn=3.0d6
    
do i=1,181
    nnum=nnum+1
    coord(nnum,1)=-90.0d0*dx1+(i-1)*dx1
    coord(nnum,2)=-0.05d0
    coord(nnum,3)=-0.11d0
enddo
do i=1,181
    nnum=nnum+1
    coord(nnum,1)=-90.0d0*dx1+(i-1)*dx1
    coord(nnum,2)=0.05d0
    coord(nnum,3)=-0.11d0
enddo
totnode1=nnum


do i=1,180
    nnum1=nnum1+1
    coord1(nnum1,1)=-90.0d0*dx1+(i-1)*dx1+(dx1/2.0d0)
    coord1(nnum1,2)=-0.05d0
    coord1(nnum1,3)=-0.11d0
enddo
do i=1,180
    nnum1=nnum1+1
    coord1(nnum1,1)=-90.0d0*dx1+(i-1)*dx1+(dx1/2.0d0)
    coord1(nnum1,2)=0.05d0
    coord1(nnum1,3)=-0.11d0
enddo
totelement1=nnum1
totelement=nnum1

do i=1,totelement1
    do j=1,totnode1
        aaa=coord(j,1)-coord1(i,1)
        bbb=coord(j,2)-coord1(i,2)
        ccc=coord(j,3)-coord1(i,3)
        if((abs(aaa+(dx1/2.0d0)).le.(1.0d-6)).and.(abs(bbb).le.(1.0d-6)).and.(abs(ccc).le.(1.0d-6)))then
            nodes(i,1)=j
        elseif((abs(aaa-(dx1/2.0d0)).le.(1.0d-6)).and.(abs(bbb).le.(1.0d-6)).and.(abs(ccc).le.(1.0d-6)))then
            nodes(i,2)=j
        endif   
    enddo
enddo
    
do i=1,totnode1
    hh(i)=0
    do m=1,7
        do n=1,7
            do k=1,7
                hh(i)=hh(i)+1
                coord3(i,1,hh(i))=coord(i,1)-3.0d0*dx+(m-1)*dx
                coord3(i,2,hh(i))=coord(i,2)-3.0d0*dx+(n-1)*dx
                coord3(i,3,hh(i))=coord(i,3)-3.0d0*dx+(k-1)*dx
            enddo
        enddo
    enddo
enddo    
    
do i=1,totnode1
    hhh(i)=0
    do j=1,hh(i)
        idist=dsqrt((coord3(i,1,j)-coord(i,1))**2+(coord3(i,2,j)-coord(i,2))**2+(coord3(i,3,j)-coord(i,3))**2)
        if((idist.le.delta).and.(idist.ge.(1.0d-6)))then
            if((coord3(i,1,j).ge.(-90.5d0*dx)).and.(coord3(i,1,j).le.(90.5d0*dx)))then
                if((coord3(i,2,j).ge.(-5.0d0*dx)).and.(coord3(i,2,j).le.(5.0d0*dx)))then
                    if((coord3(i,3,j).ge.(-8.5d0*dx)).and.(coord3(i,3,j).le.(8.5d0*dx)))then
                        hhh(i)=hhh(i)+1
                        coordcoordcoord(i,1,hhh(i))=coord3(i,1,j)
                        coordcoordcoord(i,2,hhh(i))=coord3(i,2,j)
                        coordcoordcoord(i,3,hhh(i))=coord3(i,3,j)
                    endif
                endif
            endif
        endif
    enddo
enddo    
    
do k=1,18
    do i=1,182
        do j=1,11
            nnum=nnum+1
            coord(nnum,1)=-91.0d0*dx+(dx/2.0d0)+(i-1)*dx
            coord(nnum,2)=-5.0d0*dx+(j-1)*dx
            coord(nnum,3)=-9.0d0*dx+(dx/2.0d0)+(k-1)*dx
        enddo
    enddo
enddo
totint=nnum    

open(26,file="eee.txt")
    do i=(1+totnode1),totint
        write(26,222) coord(i,1),coord(i,2),coord(i,3)
    enddo
close(26)

do k=1,17
    do i=1,181
        do j=1,10
            nnum2=nnum2+1
            coordcoord(nnum2,1)=-91.0d0*dx+i*dx
            coordcoord(nnum2,2)=-5.0d0*dx+(dx/2.0d0)+(j-1)*dx
            coordcoord(nnum2,3)=-9.0d0*dx+k*dx
        enddo
    enddo
enddo
totele=nnum2

open(26,file="ccc.txt")
    do i=1,totele
        write(26,222) coordcoord(i,1),coordcoord(i,2),coordcoord(i,3)
    enddo
close(26)

do i=1,totele
    do j=(totnode1+1),totint
        aaa=coord(j,1)-coordcoord(i,1)
        bbb=coord(j,2)-coordcoord(i,2)
        ccc=coord(j,3)-coordcoord(i,3)
        if((abs(aaa+(dx/2.0d0)).le.(1.0d-6)).and.(abs(bbb+(dx/2.0d0)).le.(1.0d-6)).and.(abs(ccc+(dx/2.0d0)).le.(1.0d-6)))then
            node(i,1)=j
        elseif((abs(aaa-(dx/2.0d0)).le.(1.0d-6)).and.(abs(bbb+(dx/2.0d0)).le.(1.0d-6)).and.(abs(ccc+(dx/2.0d0)).le.(1.0d-6)))then
            node(i,2)=j
        elseif((abs(aaa-(dx/2.0d0)).le.(1.0d-6)).and.(abs(bbb-(dx/2.0d0)).le.(1.0d-6)).and.(abs(ccc+(dx/2.0d0)).le.(1.0d-6)))then     
            node(i,3)=j
        elseif((abs(aaa+(dx/2.0d0)).le.(1.0d-6)).and.(abs(bbb-(dx/2.0d0)).le.(1.0d-6)).and.(abs(ccc+(dx/2.0d0)).le.(1.0d-6)))then
            node(i,4)=j
        elseif((abs(aaa+(dx/2.0d0)).le.(1.0d-6)).and.(abs(bbb+(dx/2.0d0)).le.(1.0d-6)).and.(abs(ccc-(dx/2.0d0)).le.(1.0d-6)))then
            node(i,5)=j
        elseif((abs(aaa-(dx/2.0d0)).le.(1.0d-6)).and.(abs(bbb+(dx/2.0d0)).le.(1.0d-6)).and.(abs(ccc-(dx/2.0d0)).le.(1.0d-6)))then
            node(i,6)=j
        elseif((abs(aaa-(dx/2.0d0)).le.(1.0d-6)).and.(abs(bbb-(dx/2.0d0)).le.(1.0d-6)).and.(abs(ccc-(dx/2.0d0)).le.(1.0d-6)))then     
            node(i,7)=j
        elseif((abs(aaa+(dx/2.0d0)).le.(1.0d-6)).and.(abs(bbb-(dx/2.0d0)).le.(1.0d-6)).and.(abs(ccc-(dx/2.0d0)).le.(1.0d-6)))then
            node(i,8)=j
        endif
    enddo
enddo

open(26,file="ddd.txt")
    do i=1,totele
        write(26,333) i,node(i,1),node(i,2),node(i,3),node(i,4),node(i,5),node(i,6),node(i,7),node(i,8)
    enddo
close(26)

!do i=1,3
!    do j=1,11
!        do k=1,3
!            nnum=nnum+1
!            coord(nnum,1)=-13.0d0*dx+(dx/2.0d0)+(i-1)*dx
!            coord(nnum,2)=-5.0d0*dx+(j-1)*dx
!            coord(nnum,3)=8.5d0*dx+k*dx
!        enddo
!    enddo
!enddo  
!do i=1,3
!    do j=1,11
!        do k=1,3
!            nnum=nnum+1
!            coord(nnum,1)=13.0d0*dx-(dx/2.0d0)-(i-1)*dx
!            coord(nnum,2)=-5.0d0*dx+(j-1)*dx
!            coord(nnum,3)=8.5d0*dx+k*dx
!        enddo
!    enddo
!enddo            
tottop=nnum    

!do i=1,4
!    do j=1,11
!        do k=1,3
!            nnum=nnum+1
!            coord(nnum,1)=-85.0d0*dx+(dx/2.0d0)+(i-1)*dx
!            coord(nnum,2)=-5.0d0*dx+(j-1)*dx
!            coord(nnum,3)=-9.0d0*dx-(dx/2.0d0)-(k-1)*dx
!        enddo
!    enddo
!enddo
!do i=1,4
!    do j=1,11
!        do k=1,3
!            nnum=nnum+1
!            coord(nnum,1)=85.0d0*dx-(dx/2.0d0)-(i-1)*dx
!            coord(nnum,2)=-5.0d0*dx+(j-1)*dx
!            coord(nnum,3)=-9.0d0*dx-(dx/2.0d0)-(k-1)*dx
!        enddo
!    enddo
!enddo            
totbottom=nnum
totnode=nnum

write(*,*) "totnode=",totnode

do i=1,totnode
    point(i)=0
enddo
do i=1,totnode1
    point(i)=1
enddo

open(26,file="aaa.txt")
    do i=1,totnode
        if(point(i)==1)then
            write(26,222) coord(i,1),coord(i,2),coord(i,3)
        endif
    enddo
close(26)

open(26,file="bbb.txt")
    do i=1,totnode
        if(point(i)==0)then
            write(26,222) coord(i,1),coord(i,2),coord(i,3)
        endif
    enddo
close(26)



do i=1,totnode1
    do j=1,hhh(i)
        do m=1,totele
            cnode=node(m,1)
            ccnode=node(m,7)
            aaa=(coordcoordcoord(i,1,j)-coord(cnode,1))*(coordcoordcoord(i,1,j)-coord(ccnode,1))
            bbb=(coordcoordcoord(i,2,j)-coord(cnode,2))*(coordcoordcoord(i,2,j)-coord(ccnode,2))
            ccc=(coordcoordcoord(i,3,j)-coord(cnode,3))*(coordcoordcoord(i,3,j)-coord(ccnode,3))
            if((aaa<=0.0d0).and.(bbb<=0.0d0).and.(ccc<=0.0d0))then
                nodenode(i,j)=m
                aaaa=(coordcoordcoord(i,1,j)-coord(cnode,1))-(dx/2.0d0)
                bbbb=(coordcoordcoord(i,2,j)-coord(cnode,2))-(dx/2.0d0)
                cccc=(coordcoordcoord(i,3,j)-coord(cnode,3))-(dx/2.0d0)
                N1(i,j)=(1-aaaa/(dx/2.0d0))*(1-bbbb/(dx/2.0d0))*(1-cccc/(dx/2.0d0))/8.0d0
                N2(i,j)=(1+aaaa/(dx/2.0d0))*(1-bbbb/(dx/2.0d0))*(1-cccc/(dx/2.0d0))/8.0d0
                N3(i,j)=(1+aaaa/(dx/2.0d0))*(1+bbbb/(dx/2.0d0))*(1-cccc/(dx/2.0d0))/8.0d0
                N4(i,j)=(1-aaaa/(dx/2.0d0))*(1+bbbb/(dx/2.0d0))*(1-cccc/(dx/2.0d0))/8.0d0
                N5(i,j)=(1-aaaa/(dx/2.0d0))*(1-bbbb/(dx/2.0d0))*(1+cccc/(dx/2.0d0))/8.0d0
                N6(i,j)=(1+aaaa/(dx/2.0d0))*(1-bbbb/(dx/2.0d0))*(1+cccc/(dx/2.0d0))/8.0d0
                N7(i,j)=(1+aaaa/(dx/2.0d0))*(1+bbbb/(dx/2.0d0))*(1+cccc/(dx/2.0d0))/8.0d0
                N8(i,j)=(1-aaaa/(dx/2.0d0))*(1+bbbb/(dx/2.0d0))*(1+cccc/(dx/2.0d0))/8.0d0
            endif
        enddo
    enddo
enddo
     
do i=(totnode1+1),totnode
    if(i.eq.(totnode1+1))then 
        pointfam(i,1)=1
    else
        pointfam(i,1)=pointfam(i-1,1)+numfam(i-1,1)
    endif
    do j=(totnode1+1),totnode
        idist=dsqrt((coord(j,1)-coord(i,1))**2+(coord(j,2)-coord(i,2))**2+(coord(j,3)-coord(i,3))**2)
        if(i.ne.j)then
            if(idist.le.delta)then
                numfam(i,1)=numfam(i,1)+1
                nodefam(pointfam(i,1)+numfam(i,1)-1,1)=j
            endif
        endif
    enddo
enddo         
    
do i=(totnode1+1),totnode
    disp(i,1)=0.001d0*coord(i,1)
    disp(i,2)=0.0d0
    disp(i,3)=0.0d0
enddo
do i=(totnode1+1),totnode
    stendens(i,1)=0.0d0
    do j=1,numfam(i,1)
        cnode=nodefam(pointfam(i,1)+j-1,1)
        idist=dsqrt((coord(cnode,1)-coord(i,1))**2+(coord(cnode,2)-coord(i,2))**2+(coord(cnode,3)-coord(i,3))**2)
        nlength=dsqrt((coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))**2+(coord(cnode,2)+disp(cnode,2)-coord(i,2)-disp(i,2))**2+(coord(cnode,3)+disp(cnode,3)-coord(i,3)-disp(i,3))**2)
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
    
do i=(totnode1+1),totnode
    disp(i,1)=0.0d0
    disp(i,2)=0.001d0*coord(i,2)
    disp(i,3)=0.0d0
enddo
do i=(totnode1+1),totnode
    stendens(i,2)=0.0d0
    do j=1,numfam(i,1)
        cnode=nodefam(pointfam(i,1)+j-1,1)
        idist=dsqrt((coord(cnode,1)-coord(i,1))**2+(coord(cnode,2)-coord(i,2))**2+(coord(cnode,3)-coord(i,3))**2)
        nlength=dsqrt((coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))**2+(coord(cnode,2)+disp(cnode,2)-coord(i,2)-disp(i,2))**2+(coord(cnode,3)+disp(cnode,3)-coord(i,3)-disp(i,3))**2)
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

do i=(totnode1+1),totnode
    disp(i,1)=0.0d0
    disp(i,2)=0.0d0
    disp(i,3)=0.001d0*coord(i,3)
enddo
do i=(totnode1+1),totnode
    stendens(i,3)=0.0d0
    do j=1,numfam(i,1)
        cnode=nodefam(pointfam(i,1)+j-1,1)
        idist=dsqrt((coord(cnode,1)-coord(i,1))**2+(coord(cnode,2)-coord(i,2))**2+(coord(cnode,3)-coord(i,3))**2)
        nlength=dsqrt((coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))**2+(coord(cnode,2)+disp(cnode,2)-coord(i,2)-disp(i,2))**2+(coord(cnode,3)+disp(cnode,3)-coord(i,3)-disp(i,3))**2)
        if(idist.le.delta-radij)then
            fac=1.0d0
        elseif(idist.le.delta+radij)then
            fac=(delta+radij-idist)/(2.0d0*radij)
        else
            fac=0.0d0
        endif        
        stendens(i,3)=stendens(i,3)+0.5d0*0.5d0*bc*((nlength-idist)/idist)**2*idist*vol*fac  
    enddo
    fncst(i,3)=sedload2/stendens(i,3)
enddo        
    
do i=1,totnode
    vel(i,1)=0.0d0
    disp(i,1)=0.0d0 
    vel(i,2)=0.0d0
    disp(i,2)=0.0d0 
    vel(i,3)=0.0d0
    disp(i,3)=0.0d0 
enddo
    
do i=1,totnode
    if(point(i)==0)then
       massvec(i,1)=0.25d0*dt*dt*((4.0d0/3.0d0)*pi*(delta)**3)*bc/dx*30.0d0
       massvec(i,2)=0.25d0*dt*dt*((4.0d0/3.0d0)*pi*(delta)**3)*bc/dx*30.0d0
       massvec(i,3)=0.25d0*dt*dt*((4.0d0/3.0d0)*pi*(delta)**3)*bc/dx*30.0d0
    else
        massvec(i,1)=0.25d0*dt*dt*((4.0d0/3.0d0)*pi*(delta)**3)*bc1/dx*30.0d0
        massvec(i,2)=0.25d0*dt*dt*((4.0d0/3.0d0)*pi*(delta)**3)*bc1/dx*30.0d0
        massvec(i,3)=0.25d0*dt*dt*((4.0d0/3.0d0)*pi*(delta)**3)*bc1/dx*30.0d0
    endif    
enddo    
    
KE1(1,1)=emod1*A1/L1
KE1(1,2)=emod1*A1/L1*(-1.0d0)
KE1(2,1)=emod1*A1/L1*(-1.0d0)
KE1(2,2)=emod1*A1/L1    


do i=1,2
    do j=1,6
        T1(i,j)=0.0d0
        T1T(j,i)=0.0d0
    enddo
enddo
T1(1,1)=1.0d0
T1(2,4)=1.0d0
T1T(1,1)=1.0d0
T1T(4,2)=1.0d0
do i=1,6
    do j=1,2
        T1TKE(i,j)=0.0d0
        do m=1,2
            T1TKE(i,j)=T1TKE(i,j)+T1T(i,m)*KE1(m,j)
        enddo
    enddo
enddo
do i=1,6
    do j=1,6
        ske1(i,j)=0.0d0
        do m=1,2
            ske1(i,j)=ske1(i,j)+T1TKE(i,m)*T1(m,j)
        enddo
    enddo
enddo    
    
open(91,file="history.csv",status="replace")
write(91,'(A)') "step,time_s,probe_z_m,contact_x_N,contact_z_N,prescribed_z_next_m"
first_yield_step=0
open(93,file="diagnostics.csv",status="replace")
write(93,'(A)') "step,time_s,max_steel_u_m,max_concrete_u_m,max_steel_speed_m_s,max_concrete_speed_m_s,yielded_rods,net_internal_x_N,net_internal_y_N,net_internal_z_N,net_interface_x_N,net_interface_y_N,net_interface_z_N"
open(94,file="first_yield.txt",status="replace")
do tt=1,200
    write(*,*) "tt=",tt

 !   if(tt<=10000)then
 !       do i=(totint+1),tottop
 !           disp(i,1)=0.0d0
 !           disp(i,2)=0.0d0
 !           disp(i,3)=((-8.0d-11)*tt)*tt*dt
 !       enddo
 !   else
 !       do i=(totint+1),tottop
 !           disp(i,1)=0.0d0
 !           disp(i,2)=0.0d0
 !           disp(i,3)=(-8.0d-3)+(-8.0d-7)*(tt-10000)*dt
 !       enddo
 !   endif
    
!    do i=(totint+1),tottop
!        disp(i,1)=0.0d0
!        disp(i,2)=0.0d0
!        if(coord(i,3)>=(170.0d-3))then
!            disp(i,1)=0.0d0
!            disp(i,2)=0.0d0
!            disp(i,3)=(-1.0d-7)*tt*dt
!        endif
!    enddo

    do i=(totnode1+1),totnode
        pforce11(i,1)=0.0d0
        pforce11(i,2)=0.0d0
        pforce11(i,3)=0.0d0
        dmgpar1=0.0d0
        dmgpar2=0.0d0
        do j=1,numfam(i,1)            
            cnode=nodefam(pointfam(i,1)+j-1,1)
            idist=dsqrt((coord(cnode,1)-coord(i,1))**2+(coord(cnode,2)-coord(i,2))**2+(coord(cnode,3)-coord(i,3))**2)
            nlength=dsqrt((coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))**2+(coord(cnode,2)+disp(cnode,2)-coord(i,2)-disp(i,2))**2+(coord(cnode,3)+disp(cnode,3)-coord(i,3)-disp(i,3))**2)
            if(idist.le.delta-radij)then
                fac=1.0d0
            elseif(idist.le.delta+radij)then
                fac=(delta+radij-idist)/(2.0d0*radij)
            else
                fac=0.0d0
            endif
            if(dabs(coord(cnode,3)-coord(i,3))<=1.0d-10)then
                if(dabs(coord(cnode,2)-coord(i,2))<=1.0d-10)then
		            theta=0.0d0
		        elseif(dabs(coord(cnode,1)-coord(i,1))<=1.0d-10)then
		            theta=90.0d0*pi/180.0d0
		        else
		            theta=datan2(dabs(coord(cnode,2)-coord(i,2)),dabs(coord(cnode,1)-coord(i,1)))
                endif
                phi=90.0d0*pi/180.0d0
		        scx=(fncst(i,1)+fncst(cnode,1))/2.0d0
		        scy=(fncst(i,2)+fncst(cnode,2))/2.0d0
		        scz=(fncst(i,3)+fncst(cnode,3))/2.0d0
		        scr=1.0d0/(((dcos(theta)*dsin(phi))**2/(scx)**2)+((dsin(theta)*dsin(phi))**2/(scy)**2)+((dcos(phi))**2/(scz)**2))
		        scr = dsqrt(scr)
	        elseif(dabs(coord(cnode,1)-coord(i,1))<=1.0d-10.and.dabs(coord(cnode,2)-coord(i,2))<=1.0d-10)then
	            scz=(fncst(i,3)+fncst(cnode,3))/2.0d0
	            scr=scz
            else
                theta=datan2(dabs(coord(cnode,2)-coord(i,2)),dabs(coord(cnode,1)-coord(i,1)))
		        phi=dacos(dabs(coord(cnode,3)-coord(i,3))/idist)
		        scx=(fncst(i,1)+fncst(cnode,1))/2.0d0
		        scy=(fncst(i,2)+fncst(cnode,2))/2.0d0
		        scz=(fncst(i,3)+fncst(cnode,3))/2.0d0
		        scr=1.0d0/(((dcos(theta)*dsin(phi))**2/(scx)**2)+((dsin(theta)*dsin(phi))**2/(scy)**2)+((dcos(phi))**2/(scz)**2))
		        scr=dsqrt(scr)
            endif       
            sssss=(nlength-idist)/idist
            if(sssss<=s0)then
                dforce1=bc*sssss*vol*scr*fac*(coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))/nlength 
                dforce2=bc*sssss*vol*scr*fac*(coord(cnode,2)+disp(cnode,2)-coord(i,2)-disp(i,2))/nlength   
                dforce3=bc*sssss*vol*scr*fac*(coord(cnode,3)+disp(cnode,3)-coord(i,3)-disp(i,3))/nlength  
            elseif((sssss>s0).and.(sssss<=sc))then
                dforce1=bc*s0*(sc-sssss)/(sc-s0)*vol*scr*fac*(coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))/nlength 
                dforce2=bc*s0*(sc-sssss)/(sc-s0)*vol*scr*fac*(coord(cnode,2)+disp(cnode,2)-coord(i,2)-disp(i,2))/nlength   
                dforce3=bc*s0*(sc-sssss)/(sc-s0)*vol*scr*fac*(coord(cnode,3)+disp(cnode,3)-coord(i,3)-disp(i,3))/nlength  
            elseif(sssss>sc)then
                dforce1=0.0d0
                dforce2=0.0d0
                dforce3=0.0d0
            endif
            pforce11(i,1)=pforce11(i,1)+dforce1      
            pforce11(i,2)=pforce11(i,2)+dforce2  
            pforce11(i,3)=pforce11(i,3)+dforce3  
 !          if(coord(i,3)<=(0.0d0))then
 !              if(abs(coord(i,1))<=(40.0d0*dx))then
                  if(((nlength-idist)/idist)>sc)then
  !                  if(((nlength-idist)/idist)>scr1)then
                        fail(i,j)=0
                  endif
 !               endif
 !           endif
            dmgpar1=dmgpar1+fail(i,j)*vol*fac
            dmgpar2=dmgpar2+vol*fac
        enddo
        dmg(i,1)=1.0d0-dmgpar1/dmgpar2
    enddo    
    
    
    loadloadx=0.0d0
    loadloady=0.0d0
    loadloadz=0.0d0
    do i=(totint+1),tottop
        loadloadx=loadloadx+pforce(i,1)*vol
        loadloady=loadloady+pforce(i,2)*vol
        loadloadz=loadloadz+pforce(i,3)*vol
    enddo
    loadxx(tt)=loadloadx
    loadyy(tt)=loadloady
    loadzz(tt)=loadloadz
    loadloadzong(tt)=dsqrt(loadloadx**2+loadloady**2+loadloadz**2)
    
    do i=(totnode1+1),totint
        if((abs((coord(i,1)-0.5d0*dx)).le.(1.0d-6)).and.(abs(coord(i,2)).le.(1.0d-6)).and.(abs(coord(i,3)+0.5d0*dx).le.(1.0d-6)))then
            dispdispld(tt)=disp(i,3)
        endif
    enddo
    
    do i=1,totelement
        cnode=nodes(i,1)
        ccnode=nodes(i,2)
        pforce22(cnode,1)=0.0d0
        pforce22(cnode,2)=0.0d0
        pforce22(cnode,3)=0.0d0
        pforce22(ccnode,1)=0.0d0
        pforce22(ccnode,2)=0.0d0
        pforce22(ccnode,3)=0.0d0
    enddo
    
    do i=1,totelement1
        cnode=nodes(i,1)
        ccnode=nodes(i,2)
        idist=dsqrt((coord(ccnode,1)-coord(cnode,1))**2+(coord(ccnode,2)-coord(cnode,2))**2+(coord(ccnode,3)-coord(cnode,3))**2)
        nlength=dsqrt((coord(ccnode,1)+disp(ccnode,1)-coord(cnode,1)-disp(cnode,1))**2+(coord(ccnode,2)+disp(ccnode,2)-coord(cnode,2)-disp(cnode,2))**2+(coord(ccnode,3)+disp(ccnode,3)-coord(cnode,3)-disp(cnode,3))**2)
        sssss=(nlength-idist)/idist
        wy(1)=disp(cnode,1)
        wy(2)=disp(cnode,2)
        wy(3)=disp(cnode,3)
        wy(4)=disp(ccnode,1)
        wy(5)=disp(ccnode,2)
        wy(6)=disp(ccnode,3)        
        do m=1,6
            skewy(m)=0.0d0
            do n=1,6
                skewy(m)=skewy(m)+ske1(m,n)*wy(n)
            enddo
        enddo
        if(sssss>sy.and.first_yield_step==0)then
            first_yield_step=tt
            write(94,*) tt,i,cnode,ccnode,sssss,disp(cnode,:),disp(ccnode,:)
            flush(94)
        endif
        if(sssss<=sy)then
            pforce22(cnode,1)=pforce22(cnode,1)+skewy(1)*(-1.0d0)/vol1
            pforce22(cnode,2)=pforce22(cnode,2)+skewy(2)*(-1.0d0)/vol1
            pforce22(cnode,3)=pforce22(cnode,3)+skewy(3)*(-1.0d0)/vol1
            pforce22(ccnode,1)=pforce22(ccnode,1)+skewy(4)*(-1.0d0)/vol1
            pforce22(ccnode,2)=pforce22(ccnode,2)+skewy(5)*(-1.0d0)/vol1
            pforce22(ccnode,3)=pforce22(ccnode,3)+skewy(6)*(-1.0d0)/vol1
        else
            pforce22(cnode,1)=pforce22(cnode,1)+416.0d0*530.66d0/vol1
            pforce22(cnode,2)=pforce22(cnode,2)+0.0d0
            pforce22(cnode,3)=pforce22(cnode,3)+0.0d0
            pforce22(ccnode,1)=pforce22(ccnode,1)-416.0d0*530.66d0/vol1
            pforce22(ccnode,2)=pforce22(ccnode,2)+0.0d0
            pforce22(ccnode,3)=pforce22(ccnode,3)+0.0d0
        endif
    enddo    
    
    do i=1,totnode
        pforce33(i,1)=0.0d0
        pforce33(i,2)=0.0d0
        pforce33(i,3)=0.0d0
    enddo
    
    do i=1,totnode1
        pforce33(i,1)=0.0d0
        pforce33(i,2)=0.0d0
        pforce33(i,3)=0.0d0
        dmgpar1=0.0d0
        dmgpar2=0.0d0
        do j=1,hhh(i)
            cnode=nodenode(i,j)
            cnode1=node(cnode,1)
            cnode2=node(cnode,2)
            cnode3=node(cnode,3)
            cnode4=node(cnode,4)
            cnode5=node(cnode,5)
            cnode6=node(cnode,6)
            cnode7=node(cnode,7)
            cnode8=node(cnode,8)
            dispdisp(i,1,j)=N1(i,j)*disp(cnode1,1)+N2(i,j)*disp(cnode2,1)+N3(i,j)*disp(cnode3,1)+N4(i,j)*disp(cnode4,1)+N5(i,j)*disp(cnode5,1)+N6(i,j)*disp(cnode6,1)+N7(i,j)*disp(cnode7,1)+N8(i,j)*disp(cnode8,1)
            dispdisp(i,2,j)=N1(i,j)*disp(cnode1,2)+N2(i,j)*disp(cnode2,2)+N3(i,j)*disp(cnode3,2)+N4(i,j)*disp(cnode4,2)+N5(i,j)*disp(cnode5,2)+N6(i,j)*disp(cnode6,2)+N7(i,j)*disp(cnode7,2)+N8(i,j)*disp(cnode8,2)
            dispdisp(i,3,j)=N1(i,j)*disp(cnode1,3)+N2(i,j)*disp(cnode2,3)+N3(i,j)*disp(cnode3,3)+N4(i,j)*disp(cnode4,3)+N5(i,j)*disp(cnode5,3)+N6(i,j)*disp(cnode6,3)+N7(i,j)*disp(cnode7,3)+N8(i,j)*disp(cnode8,3)
            idist=dsqrt((coordcoordcoord(i,1,j)-coord(i,1))**2+(coordcoordcoord(i,2,j)-coord(i,2))**2+(coordcoordcoord(i,3,j)-coord(i,3))**2)
            nlength=dsqrt((coordcoordcoord(i,1,j)+dispdisp(i,1,j)-coord(i,1)-disp(i,1))**2+(coordcoordcoord(i,2,j)+dispdisp(i,2,j)-coord(i,2)-disp(i,2))**2+(coordcoordcoord(i,3,j)+dispdisp(i,3,j)-coord(i,3)-disp(i,3))**2)
            sssss=(nlength-idist)/idist
            if(sssss<=s0)then            
                dforce1=bc*sssss*vol*(coordcoordcoord(i,1,j)+dispdisp(i,1,j)-coord(i,1)-disp(i,1))/nlength 
                dforce2=bc*sssss*vol*(coordcoordcoord(i,2,j)+dispdisp(i,2,j)-coord(i,2)-disp(i,2))/nlength
                dforce3=bc*sssss*vol*(coordcoordcoord(i,3,j)+dispdisp(i,3,j)-coord(i,3)-disp(i,3))/nlength
            elseif((sssss>s0).and.(sssss<=sc))then
                dforce1=bc*s0*(sc-sssss)/(sc-s0)*vol*(coordcoordcoord(i,1,j)+dispdisp(i,1,j)-coord(i,1)-disp(i,1))/nlength 
                dforce2=bc*s0*(sc-sssss)/(sc-s0)*vol*(coordcoordcoord(i,2,j)+dispdisp(i,2,j)-coord(i,2)-disp(i,2))/nlength
                dforce3=bc*s0*(sc-sssss)/(sc-s0)*vol*(coordcoordcoord(i,3,j)+dispdisp(i,3,j)-coord(i,3)-disp(i,3))/nlength
            elseif(sssss>sc)then
                dforce1=0.0d0
                dforce2=0.0d0
                dforce3=0.0d0
            endif
            pforce33(i,1)=pforce33(i,1)+dforce1    
            pforce33(i,2)=pforce33(i,2)+dforce2     
            pforce33(i,3)=pforce33(i,3)+dforce3   
            pforce33(cnode1,1)=pforce33(cnode1,1)+N1(i,j)*(-1.0d0)*dforce1*(vol1/vol)
            pforce33(cnode1,2)=pforce33(cnode1,2)+N1(i,j)*(-1.0d0)*dforce2*(vol1/vol)
            pforce33(cnode1,3)=pforce33(cnode1,3)+N1(i,j)*(-1.0d0)*dforce3*(vol1/vol)
            pforce33(cnode2,1)=pforce33(cnode2,1)+N2(i,j)*(-1.0d0)*dforce1*(vol1/vol)
            pforce33(cnode2,2)=pforce33(cnode2,2)+N2(i,j)*(-1.0d0)*dforce2*(vol1/vol)
            pforce33(cnode2,3)=pforce33(cnode2,3)+N2(i,j)*(-1.0d0)*dforce3*(vol1/vol)
            pforce33(cnode3,1)=pforce33(cnode3,1)+N3(i,j)*(-1.0d0)*dforce1*(vol1/vol)
            pforce33(cnode3,2)=pforce33(cnode3,2)+N3(i,j)*(-1.0d0)*dforce2*(vol1/vol)
            pforce33(cnode3,3)=pforce33(cnode3,3)+N3(i,j)*(-1.0d0)*dforce3*(vol1/vol)
            pforce33(cnode4,1)=pforce33(cnode4,1)+N4(i,j)*(-1.0d0)*dforce1*(vol1/vol)
            pforce33(cnode4,2)=pforce33(cnode4,2)+N4(i,j)*(-1.0d0)*dforce2*(vol1/vol)
            pforce33(cnode4,3)=pforce33(cnode4,3)+N4(i,j)*(-1.0d0)*dforce3*(vol1/vol)
            pforce33(cnode5,1)=pforce33(cnode5,1)+N5(i,j)*(-1.0d0)*dforce1*(vol1/vol)
            pforce33(cnode5,2)=pforce33(cnode5,2)+N5(i,j)*(-1.0d0)*dforce2*(vol1/vol)
            pforce33(cnode5,3)=pforce33(cnode5,3)+N5(i,j)*(-1.0d0)*dforce3*(vol1/vol)
            pforce33(cnode6,1)=pforce33(cnode6,1)+N6(i,j)*(-1.0d0)*dforce1*(vol1/vol)
            pforce33(cnode6,2)=pforce33(cnode6,2)+N6(i,j)*(-1.0d0)*dforce2*(vol1/vol)
            pforce33(cnode6,3)=pforce33(cnode6,3)+N6(i,j)*(-1.0d0)*dforce3*(vol1/vol)
            pforce33(cnode7,1)=pforce33(cnode7,1)+N7(i,j)*(-1.0d0)*dforce1*(vol1/vol)
            pforce33(cnode7,2)=pforce33(cnode7,2)+N7(i,j)*(-1.0d0)*dforce2*(vol1/vol)
            pforce33(cnode7,3)=pforce33(cnode7,3)+N7(i,j)*(-1.0d0)*dforce3*(vol1/vol)
            pforce33(cnode8,1)=pforce33(cnode8,1)+N8(i,j)*(-1.0d0)*dforce1*(vol1/vol)
            pforce33(cnode8,2)=pforce33(cnode8,2)+N8(i,j)*(-1.0d0)*dforce2*(vol1/vol)
            pforce33(cnode8,3)=pforce33(cnode8,3)+N8(i,j)*(-1.0d0)*dforce3*(vol1/vol)
      !     if(coord(i,3)<=(0.0d0))then
                if(((nlength-idist)/idist)>sc)then
      !          if(((nlength-idist)/idist)>scr1)then
                    failfail(i,j)=0
                endif
      !     endif
            dmgpar1=dmgpar1+failfail(i,j)*vol
            dmgpar2=dmgpar2+vol
        enddo
        dmg(i,1)=1.0d0-dmgpar1/dmgpar2
    enddo
            
    do i=1,totnode
        pforce(i,1)=pforce11(i,1)+pforce22(i,1)+pforce33(i,1)
        pforce(i,2)=pforce11(i,2)+pforce22(i,2)+pforce33(i,2)
        pforce(i,3)=pforce11(i,3)+pforce22(i,3)+pforce33(i,3)
    enddo    

!    cn=0.0d0
!	cn1=0.0d0
!	cn2=0.0d0
!	do i=1,totnode
!        if(velhalfold(i,1).ne.0.0d0)then
!            cn1=cn1-disp(i,1)*disp(i,1)*((pforce(i,1)+bforce(i,1))/massvec(i,1)-(pforceold(i,1)+bforce(i,1))/massvec(i,1))/(dt*velhalfold(i,1))
!        endif
!        if(velhalfold(i,2).ne.0.0d0)then
!            cn1=cn1-disp(i,2)*disp(i,2)*((pforce(i,2)+bforce(i,2))/massvec(i,2)-(pforceold(i,2)+bforce(i,2))/massvec(i,2))/(dt*velhalfold(i,2))
!        endif
!        if(velhalfold(i,3).ne.0.0d0)then
!            cn1=cn1-disp(i,3)*disp(i,3)*((pforce(i,3)+bforce(i,3))/massvec(i,3)-(pforceold(i,3)+bforce(i,3))/massvec(i,3))/(dt*velhalfold(i,3))
!        endif
!		cn2=cn2+disp(i,1)*disp(i,1)
!		cn2=cn2+disp(i,2)*disp(i,2)
!		cn2=cn2+disp(i,3)*disp(i,3)
 !   enddo

!	if(cn2.ne.0.0d0)then
!        if((cn1/cn2)>0.0d0)then 
!            cn=2.0d0*dsqrt(cn1/cn2)
!        else
!            cn=0.0d0
!        endif
!    else
!        cn=0.0d0
!    endif
!	if(cn>2.0d0)then
!		cn=1.9d0
!	endif

!	do i=1,tottop
!		if(tt.eq.1)then
 !           velhalf(i,1)=1.0d0*dt/massvec(i,1)*(pforce(i,1)+bforce(i,1))/2.0d0	
 !           velhalf(i,2)=1.0d0*dt/massvec(i,2)*(pforce(i,2)+bforce(i,2))/2.0d0
 !           velhalf(i,3)=1.0d0*dt/massvec(i,3)*(pforce(i,3)+bforce(i,3))/2.0d0
 !       else	
 !           velhalf(i,1)=((2.0d0-cn*dt)*velhalfold(i,1)+2.0d0*dt/massvec(i,1)*(pforce(i,1)+bforce(i,1)))/(2.0d0+cn*dt)
 !           velhalf(i,2)=((2.0d0-cn*dt)*velhalfold(i,2)+2.0d0*dt/massvec(i,2)*(pforce(i,2)+bforce(i,2)))/(2.0d0+cn*dt)
 !           velhalf(i,3)=((2.0d0-cn*dt)*velhalfold(i,3)+2.0d0*dt/massvec(i,3)*(pforce(i,3)+bforce(i,3)))/(2.0d0+cn*dt)
 !       endif 
 !       vel(i,1)=0.5d0*(velhalfold(i,1)+velhalf(i,1))	
 !       vel(i,2)=0.5d0*(velhalfold(i,2)+velhalf(i,2))	
 !       vel(i,3)=0.5d0*(velhalfold(i,3)+velhalf(i,3))	
 !       disp(i,1)=disp(i,1)+velhalf(i,1)*dt	
 !       disp(i,2)=disp(i,2)+velhalf(i,2)*dt	   
 !       disp(i,3)=disp(i,3)+velhalf(i,3)*dt	 
 !       velhalfold(i,1)=velhalf(i,1)
 !       velhalfold(i,2)=velhalf(i,2)
 !       velhalfold(i,3)=velhalf(i,3)
!		pforceold(i,1)=pforce(i,1)
!		pforceold(i,2)=pforce(i,2)
!		pforceold(i,3)=pforce(i,3)
!    enddo    
    
    
    
net_internal=0.0d0
net_interface=0.0d0
do i=1,totnode
    node_volume=vol
    if(point(i)==1)node_volume=vol1
    net_internal=net_internal+pforce(i,:)*node_volume
    net_interface=net_interface+pforce33(i,:)*node_volume
enddo
 do i=1,tottop
    if(point(i)==0)then
        acc(i,1)=(pforce(i,1)+bforce(i,1)-cn*vel(i,1))/dens
        acc(i,2)=(pforce(i,2)+bforce(i,2)-cn*vel(i,2))/dens 
        acc(i,3)=(pforce(i,3)+bforce(i,3)-cn*vel(i,3))/dens 
        oldvel(i,1)=vel(i,1)
        oldvel(i,2)=vel(i,2)
        oldvel(i,3)=vel(i,3)
        vel(i,1)=vel(i,1)+acc(i,1)*dt
        vel(i,2)=vel(i,2)+acc(i,2)*dt
        vel(i,3)=vel(i,3)+acc(i,3)*dt
        olddisp(i,1)=disp(i,1)
        olddisp(i,2)=disp(i,2)
        olddisp(i,3)=disp(i,3)
        disp(i,1)=disp(i,1)+vel(i,1)*dt
        disp(i,2)=disp(i,2)+vel(i,2)*dt
        disp(i,3)=disp(i,3)+vel(i,3)*dt
    else
        acc(i,1)=(pforce(i,1)+bforce(i,1)-cn*vel(i,1))/dengs
        acc(i,2)=(pforce(i,2)+bforce(i,2)-cn*vel(i,2))/dengs 
        acc(i,3)=(pforce(i,3)+bforce(i,3)-cn*vel(i,3))/dengs 
        oldvel(i,1)=vel(i,1)
        oldvel(i,2)=vel(i,2)
        oldvel(i,3)=vel(i,3)
        vel(i,1)=vel(i,1)+acc(i,1)*dt
        vel(i,2)=vel(i,2)+acc(i,2)*dt
        vel(i,3)=vel(i,3)+acc(i,3)*dt
        olddisp(i,1)=disp(i,1)
        olddisp(i,2)=disp(i,2)
        olddisp(i,3)=disp(i,3)
        disp(i,1)=disp(i,1)+vel(i,1)*dt
        disp(i,2)=disp(i,2)+vel(i,2)*dt
        disp(i,3)=disp(i,3)+vel(i,3)*dt
    endif    
 enddo   
 
 
    do i=1,totnode
        if((abs(coord(i,1)+10.5d0*dx).le.(1.0d-6)).and.(abs(coord(i,3)-8.5d0*dx).le.(1.0d-6)))then
            disp(i,1)=0.0d0
            vel(i,1)=0.0d0
        endif
    enddo
    do i=1,totnode
        if((abs(coord(i,1)-10.5d0*dx).le.(1.0d-6)).and.(abs(coord(i,3)-8.5d0*dx).le.(1.0d-6)))then
            disp(i,1)=0.0d0
            vel(i,1)=0.0d0
        endif
    enddo
 
   
   do i=1,totnode
        if((abs(coord(i,1)+81.5d0*dx).le.(1.0d-6)).and.(abs(coord(i,3)+8.5d0*dx).le.(1.0d-6)))then
            disp(i,3)=0.0d0
            vel(i,3)=0.0d0
        endif
    enddo
    do i=1,totnode
        if((abs(coord(i,1)-81.5d0*dx).le.(1.0d-6)).and.(abs(coord(i,3)+8.5d0*dx).le.(1.0d-6)))then
            disp(i,3)=0.0d0
            vel(i,3)=0.0d0
        endif
    enddo
    

    
    do i=1,totnode
        dispdispdisp(i)=dsqrt(disp(i,1)**2+disp(i,2)**2+disp(i,3)**2)
    enddo
    
    
    
    
    
    
    
    
    
    
    
    
    pentfx=0.0d0
    pentfy=0.0d0
    pentfz=0.0d0
    penccx=penicx+penux
    penccy=0.0d0
    penccz=penicz+penuz
    do i=1,totnode
        newpx=coord(i,1)+disp(i,1)
        newpy=0.0d0
        newpz=coord(i,3)+disp(i,3)
        dpenx=newpx-penccx
        dpeny=0.0d0
        dpenz=newpz-penccz
        dpend=dsqrt(dpenx*dpenx+dpeny*dpeny+dpenz*dpenz)
        if(dpend.lt.penrad)then
            dpux=dpenx/dpend
            dpuy=dpeny/dpend
            dpuz=dpenz/dpend
            dprx=dpux*penrad
            dpry=dpuy*penrad
            dprz=dpuz*penrad            
            dpcx=penccx+dprx 
            dpcy=0.0d0
            dpcz=penccz+dprz
            olddispx=disp(i,1)
            olddispy=0.0d0
            olddispz=disp(i,3)
            oldvelx=vel(i,1)
            oldvely=0.0d0
            oldvelz=vel(i,3)
            disp(i,1)=dpcx-coord(i,1)
            disp(i,3)=dpcz-coord(i,3)
            vel(i,1)=(disp(i,1)-olddisp(i,1))/ddt
            vel(i,3)=(disp(i,3)-olddisp(i,3))/ddt
            if(point(i)==0)then
            dpenfx=-1.0d0*dens*(vel(i,1)-oldvelx)/ddt*vol
            else
                dpenfx=-1.0d0*dengs*(vel(i,1)-oldvelx)/ddt*vol1
            endif
            dpenfy=0.0d0
            if(point(i)==0)then
            dpenfz=-1.0d0*dens*(vel(i,3)-oldvelz)/ddt*vol
            else
                dpenfz=-1.0d0*dengs*(vel(i,3)-oldvelz)/ddt*vol1
            endif
            pentfx=pentfx+dpenfx
            pentfy=0.0d0
            pentfz=pentfz+dpenfz    
        endif
    enddo
    
    ppenccx=ppenicx+ppenux
    ppenccy=0.0d0
    ppenccz=ppenicz+ppenuz
    do i=1,totnode
        pnewpx=coord(i,1)+disp(i,1)
        pnewpy=0.0d0
        pnewpz=coord(i,3)+disp(i,3)
        pdpenx=pnewpx-ppenccx
        pdpeny=0.0d0
        pdpenz=pnewpz-ppenccz
        pdpend=dsqrt(pdpenx*pdpenx+pdpeny*pdpeny+pdpenz*pdpenz)
        if(pdpend.lt.penrad)then
            dpux=pdpenx/pdpend
            dpuy=pdpeny/pdpend
            dpuz=pdpenz/pdpend
            dprx=dpux*penrad
            dpry=dpuy*penrad
            dprz=dpuz*penrad            
            dpcx=ppenccx+dprx 
            dpcy=0.0d0
            dpcz=ppenccz+dprz
            olddispx=disp(i,1)
            olddispy=0.0d0
            olddispz=disp(i,3)
            oldvelx=vel(i,1)
            oldvely=0.0d0
            oldvelz=vel(i,3)
            disp(i,1)=dpcx-coord(i,1)
            disp(i,3)=dpcz-coord(i,3)
            vel(i,1)=(disp(i,1)-olddisp(i,1))/ddt
            vel(i,3)=(disp(i,3)-olddisp(i,3))/ddt
            if(point(i)==0)then
            dpenfx=-1.0d0*dens*(vel(i,1)-oldvelx)/ddt*vol
            else
                dpenfx=-1.0d0*dengs*(vel(i,1)-oldvelx)/ddt*vol1
            endif
            dpenfy=0.0d0
            if(point(i)==0)then
            dpenfz=-1.0d0*dens*(vel(i,3)-oldvelz)/ddt*vol
            else
                dpenfz=-1.0d0*dengs*(vel(i,3)-oldvelz)/ddt*vol1
            endif
            pentfx=pentfx+dpenfx
            pentfy=0.0d0
            pentfz=pentfz+dpenfz    
        endif
    enddo
    
    forceld(tt,1)=pentfx
    forceld(tt,2)=pentfy
    forceld(tt,3)=pentfz
    
    
    aaaaaa=1.0d0*tt/100000.0d0
    bbbbbb=(1.0d-2)*(aaaaaa**3)*(10.0d0-15.0d0*aaaaaa+6.0d0*(aaaaaa**2))
    
    penux=0.0d0
    penuy=0.0d0
    penuz=(-1.0d0)*bbbbbb
    ppenux=0.0d0
    ppenuy=0.0d0
    ppenuz=(-1.0d0)*bbbbbb 
    
    
    
    
    
    
    
    
    
    
write(91,'(I0,5(A,ES25.16E3))') tt,",",tt*dt,",",dispdispld(tt),",",forceld(tt,1),",",forceld(tt,3),",",penuz
if(mod(tt,100)==0.or.tt==first_yield_step)then
    max_us=0.0d0
    max_uc=0.0d0
    max_vs=0.0d0
    max_vc=0.0d0
    do i=1,totnode
        if(point(i)==1)then
            max_us=max(max_us,dsqrt(sum(disp(i,:)**2)))
            max_vs=max(max_vs,dsqrt(sum(vel(i,:)**2)))
        else
            max_uc=max(max_uc,dsqrt(sum(disp(i,:)**2)))
            max_vc=max(max_vc,dsqrt(sum(vel(i,:)**2)))
        endif
    enddo
    yielded_rods=0
    do i=1,totelement1
        cnode=nodes(i,1)
        ccnode=nodes(i,2)
        idist=dsqrt(sum((coord(ccnode,:)-coord(cnode,:))**2))
        nlength=dsqrt(sum((coord(ccnode,:)+disp(ccnode,:)-coord(cnode,:)-disp(cnode,:))**2))
        if((nlength-idist)/idist>sy)yielded_rods=yielded_rods+1
    enddo
    write(93,'(I0,5(A,ES25.16E3),A,I0,6(A,ES25.16E3))') tt,",",tt*dt,",",max_us,",",max_uc,",",max_vs,",",max_vc,",",yielded_rods,",",net_internal(1),",",net_internal(2),",",net_internal(3),",",net_interface(1),",",net_interface(2),",",net_interface(3)
    flush(93)
endif
   if(tt.eq.nt)then
        open(26,file="coord_disp_pd_100000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_100000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_100000.txt")
		    do i=1,nt
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_100000.txt")
		    do i=1,nt
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.40000)then
        open(26,file="coord_disp_pd_40000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_40000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_40000.txt")
		    do i=1,40000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_40000.txt")
		    do i=1,40000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.50000)then
        open(26,file="coord_disp_pd_50000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_50000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_50000.txt")
		    do i=1,50000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_50000.txt")
		    do i=1,50000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.60000)then
        open(26,file="coord_disp_pd_60000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_60000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_60000.txt")
		    do i=1,60000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_60000.txt")
		    do i=1,60000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.70000)then
        open(26,file="coord_disp_pd_70000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_70000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_70000.txt")
		    do i=1,70000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_70000.txt")
		    do i=1,70000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.80000)then
        open(26,file="coord_disp_pd_80000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_80000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_80000.txt")
		    do i=1,80000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_80000.txt")
		    do i=1,80000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.90000)then
        open(26,file="coord_disp_pd_90000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_90000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_90000.txt")
		    do i=1,90000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_90000.txt")
		    do i=1,90000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.20000)then
        open(26,file="coord_disp_pd_20000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_20000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_20000.txt")
		    do i=1,20000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_20000.txt")
		    do i=1,20000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.25000)then
        open(26,file="coord_disp_pd_25000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_25000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_25000.txt")
		    do i=1,25000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_25000.txt")
		    do i=1,25000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
        
        
        
        
        
        
        
    elseif(tt.eq.65000)then
        open(26,file="coord_disp_pd_65000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_65000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_65000.txt")
		    do i=1,65000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_65000.txt")
		    do i=1,65000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
        
        
        
        
        
        
        
    elseif(tt.eq.61000)then
        open(26,file="coord_disp_pd_61000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_61000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_61000.txt")
		    do i=1,61000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_61000.txt")
		    do i=1,61000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
        
        
        
        
    elseif(tt.eq.62000)then
        open(26,file="coord_disp_pd_62000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_62000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_62000.txt")
		    do i=1,62000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_62000.txt")
		    do i=1,62000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
        
        
        
    elseif(tt.eq.63000)then
        open(26,file="coord_disp_pd_63000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_63000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_63000.txt")
		    do i=1,63000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_63000.txt")
		    do i=1,63000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
        
        
        
        
        
        
        
    elseif(tt.eq.64000)then
        open(26,file="coord_disp_pd_64000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_64000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_64000.txt")
		    do i=1,65000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_64000.txt")
		    do i=1,64000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
        
        
    elseif(tt.eq.66000)then
        open(26,file="coord_disp_pd_66000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_66000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_66000.txt")
		    do i=1,66000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_66000.txt")
		    do i=1,66000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
        
        
    elseif(tt.eq.67000)then
        open(26,file="coord_disp_pd_67000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_67000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_67000.txt")
		    do i=1,67000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_67000.txt")
		    do i=1,67000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
        
        
        
    elseif(tt.eq.68000)then
        open(26,file="coord_disp_pd_68000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_68000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_68000.txt")
		    do i=1,68000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_68000.txt")
		    do i=1,68000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
        
        
    elseif(tt.eq.69000)then
        open(26,file="coord_disp_pd_69000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_69000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_69000.txt")
		    do i=1,65000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_69000.txt")
		    do i=1,69000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
        
        



    elseif(tt.eq.55000)then
        open(26,file="coord_disp_pd_55000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_55000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_55000.txt")
		    do i=1,55000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_55000.txt")
		    do i=1,55000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)





    elseif(tt.eq.30000)then
        open(26,file="coord_disp_pd_30000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_30000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_30000.txt")
		    do i=1,30000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_30000.txt")
		    do i=1,30000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.35000)then
        open(26,file="coord_disp_pd_35000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_35000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_35000.txt")
		    do i=1,35000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_35000.txt")
		    do i=1,35000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.10000)then
        open(26,file="coord_disp_pd_10000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_10000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_10000.txt")
		    do i=1,10000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_10000.txt")
		    do i=1,10000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.16000)then
        open(26,file="coord_disp_pd_16000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_16000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_16000.txt")
		    do i=1,16000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_16000.txt")
		    do i=1,16000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.14000)then
        open(26,file="coord_disp_pd_14000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_14000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_14000.txt")
		    do i=1,14000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_14000.txt")
		    do i=1,14000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.18000)then
        open(26,file="coord_disp_pd_18000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_18000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_18000.txt")
		    do i=1,18000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_18000.txt")
		    do i=1,18000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.12000)then
        open(26,file="coord_disp_pd_12000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_12000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_12000.txt")
		    do i=1,12000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_12000.txt")
		    do i=1,12000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.8000)then
        open(26,file="coord_disp_pd_8000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_8000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_8000.txt")
		    do i=1,8000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_8000.txt")
		    do i=1,8000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.6000)then
        open(26,file="coord_disp_pd_6000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_6000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_6000.txt")
		    do i=1,6000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_6000.txt")
		    do i=1,6000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.4000)then
        open(26,file="coord_disp_pd_4000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_4000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_4000.txt")
		    do i=1,4000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_4000.txt")
		    do i=1,4000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.200)then
        open(26,file="coord_disp_pd_200.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_200.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_200.txt")
		    do i=1,200
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_200.txt")
		    do i=1,200
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.2000)then
        open(26,file="coord_disp_pd_2000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_2000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_2000.txt")
		    do i=1,2000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_2000.txt")
		    do i=1,2000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    elseif(tt.eq.5000)then
        open(26,file="coord_disp_pd_5000.txt")
            do i=1,totint
                write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
            enddo
        close(26)
        open(26,file="coord_5000.txt")
            do i=1,totint
                if(point(i)==1)then
                    write(26,222) coord(i,1),coord(i,2),coord(i,3),disp(i,1),disp(i,2),disp(i,3),dispdispdisp(i),dmg(i,1)
                endif
            enddo
        close(26)
		open(26,file="load_disp_5000.txt")
		    do i=1,5000
			    write(26,111) dispdispld(i),loadxx(i),loadyy(i),loadzz(i),loadloadzong(i)
            enddo
		close(26)
		open(26,file="load_displd_5000.txt")
		    do i=1,5000
			    write(26,111) dispdispld(i),forceld(i,1),forceld(i,3)
            enddo
		close(26)
    endif
enddo 
111 format(e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5)
222 format(e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5,3x,e12.5)
333 format(I12,3x,I12,3x,I12,3x,I12,3x,I12,3x,I12,3x,I12,3x,I12,3x,I12,3x,I12,3x,I12,3x,I12,3x)
close(93)
close(94)
close(91)
open(92,file="final_state.txt",status="replace")
do i=1,totnode
write(92,'(I0,10(1X,ES25.16E3))') i,coord(i,:),disp(i,:),vel(i,:),dmg(i,1)
enddo
close(92)
end program main    
