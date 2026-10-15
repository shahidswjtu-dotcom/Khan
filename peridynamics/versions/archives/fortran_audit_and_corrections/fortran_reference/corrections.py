"""Explicitly selected source corrections; original source is never modified."""
import re

AVAILABLE = {"steel_components", "interface_balance", "support_velocity", "contact_mass", "contact_tangent"}


def apply_corrections(source, selected):
    if selected:
        source = source.replace("\r\n", "\n")
    unknown = set(selected) - AVAILABLE
    if unknown:
        raise ValueError("Unknown corrections: " + repr(unknown))
    if "steel_components" in selected:
        for component in (2, 3):
            old = f"pforce22(ccnode,{component})=pforce22(ccnode,1)+0.0d0"
            new = f"pforce22(ccnode,{component})=pforce22(ccnode,{component})+0.0d0"
            assert source.count(old) == 1
            source = source.replace(old, new)
    if "interface_balance" in selected:
        # Preserve steel force density (and its constitutive calibration).
        # Concrete force density = -N_a * f_s * V_s/V_c, so physical
        # forces and moments cancel for partition-of-unity interpolation.
        pattern = r"(N[1-8]\(i,j\)\*\(-1\.0d0\)\*dforce[123])"
        source, count = re.subn(pattern, r"\1*(vol1/vol)", source)
        assert count == 24
    if "support_velocity" in selected:
        # Only the original four displacement constraints; do not alter
        # support locations or add anchoring components.
        start = source.index("    do i=1,totnode\n        if((abs(coord(i,1)+10.5d0*dx)")
        end = source.index("        dispdispdisp(i)=", start)
        block = source[start:end]
        for c in (1, 3):
            old = f"            disp(i,{c})=0.0d0"
            assert block.count(old) == 2
            block = block.replace(old, old + f"\n            vel(i,{c})=0.0d0")
        source = source[:start] + block + source[end:]
    if "contact_mass" in selected:
        # Reaction must be minus the contacted node's actual momentum change.
        start = source.index("    pentfx=0.0d0", source.index("do tt=1,"))
        block = source[start:]
        for c, name, velocity, prior in [(1, "dpenfx", "vel(i,1)", "oldvelx"),
                                          (3, "dpenfz", "vel(i,3)", "oldvelz")]:
            old = f"            {name}=-1.0d0*dens*({velocity}-{prior})/ddt*vol"
            new = f"            if(point(i)==0)then\n{old}\n            else\n" + \
                  f"                {name}=-1.0d0*dengs*({velocity}-{prior})/ddt*vol1\n            endif"
            assert block.count(old) == 2
            block = block.replace(old, new)
        source = source[:start] + block
    if "contact_tangent" in selected:
        # The declared cylinder's normal lies in x-z. Frictionless projection
        # has no y correction or y impulse; preserve its tangential motion.
        start = source.index("    pentfx=0.0d0", source.index("do tt=1,"))
        end = source.index("    forceld(tt,1)=", start)
        block = source[start:end]
        for old in ["            disp(i,2)=0.0d0\n", "            vel(i,2)=0.0d0\n"]:
            assert block.count(old) == 2
            block = block.replace(old, "")
        source = source[:start] + block + source[end:]
    return source


def add_diagnostics(source):
    source = source.replace("\r\n", "\n")
    declarations = """integer first_yield_step, yielded_rods
real*8 max_us,max_uc,max_vs,max_vc,net_internal(3),net_interface(3),node_volume
"""
    source = source.replace("implicit none\n", "implicit none\n" + declarations, 1)
    source = source.replace("do tt=1,", 'first_yield_step=0\n'
        'open(93,file="diagnostics.csv",status="replace")\n'
        'write(93,\'(A)\') "step,time_s,max_steel_u_m,max_concrete_u_m,max_steel_speed_m_s,max_concrete_speed_m_s,yielded_rods,net_internal_x_N,net_internal_y_N,net_internal_z_N,net_interface_x_N,net_interface_y_N,net_interface_z_N"\n'
        'open(94,file="first_yield.txt",status="replace")\n'
        "do tt=1,", 1)
    old = "        if(sssss<=sy)then"
    assert source.count(old) == 1
    source = source.replace(old, """        if(sssss>sy.and.first_yield_step==0)then
            first_yield_step=tt
            write(94,*) tt,i,cnode,ccnode,sssss,disp(cnode,:),disp(ccnode,:)
            flush(94)
        endif
""" + old)
    # Record force sums before integration/contact and before constraints.
    old = " do i=1,tottop\n    if(point(i)==0)then"
    assert source.count(old) == 1
    source = source.replace(old, """net_internal=0.0d0
net_interface=0.0d0
do i=1,totnode
    node_volume=vol
    if(point(i)==1)node_volume=vol1
    net_internal=net_internal+pforce(i,:)*node_volume
    net_interface=net_interface+pforce33(i,:)*node_volume
enddo
""" + old)
    old = "   if(tt.eq.nt)then"
    diag = """if(mod(tt,100)==0.or.tt==first_yield_step)then
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
"""
    assert source.count(old) == 1
    source = source.replace(old, diag + old)
    source = source.replace("close(91)", "close(93)\nclose(94)\nclose(91)", 1)
    return source
