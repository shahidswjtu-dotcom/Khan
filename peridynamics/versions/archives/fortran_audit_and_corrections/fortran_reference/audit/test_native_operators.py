"""Native compiled regression tests using exact operator blocks from the supplied source.

Run with ``python audit/test_native_operators.py`` or pytest. The full original
source and corrected generated source supply the blocks; no rewritten Python
force operator is used. Only tiny harness meshes/inputs are synthetic.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from corrections import apply_corrections

V = 19e-3*530.66e-6
VC = (19e-3)**3
CAP = 416*530.66
SOURCE = (ROOT/'Console2.original.f90').read_text().replace('\r\n','\n')

def extract_stiffness(source):
    start = source.index('KE1(1,1)=emod1*A1/L1')
    end = source.index('do tt=1,nt',start)
    return source[start:end]

def extract_rod(source):
    marker = '    do i=1,totelement1\n        cnode=nodes(i,1)\n        ccnode=nodes(i,2)\n        idist='
    start = source.index(marker)
    end = source.index('    do i=1,totnode\n        pforce33(i,1)=0.0d0',start)
    return source[start:end]

def extract_interface(source):
    start = source.index('            pforce33(i,1)=pforce33(i,1)+dforce1')
    last = '            pforce33(cnode8,3)='
    end = source.index('\n',source.index(last,start))
    return source[start:end]

def compile_native(program,name,folder):
    compiler = shutil.which('gfortran') or str(ROOT/'toolchain/usr/bin/gfortran-14')
    source = folder/(name+'.f90');source.write_text(program)
    executable = folder/name
    flags = ['-std=legacy','-ffree-line-length-none','-O0','-g','-fcheck=all',
             '-fbacktrace','-ffpe-trap=invalid,overflow,zero']
    if str(ROOT/'toolchain') in compiler:
        flags += ['-B/usr/lib/gcc/x86_64-linux-gnu/14/',
                  '-L'+str(ROOT/'toolchain/usr/lib/x86_64-linux-gnu'),
                  '-Wl,-rpath,'+str(ROOT/'toolchain/usr/lib/x86_64-linux-gnu')]
    subprocess.run([compiler,*flags,str(source),'-o',str(executable)],check=True,capture_output=True,text=True)
    return executable

def rod_program(source):
    return '''program native_rod
implicit none
integer i,j,m,n,cnode,ccnode,totelement1,nodes(1,2)
real*8 coord(2,3),disp(2,3),pforce22(2,3),idist,nlength,sssss,wy(6),skewy(6),sy,vol1
real*8 KE1(2,2),ske1(6,6),T1(2,6),T1T(6,2),T1TKE(6,2),emod1,A1,L1,pi,strain,seed
character(64) arg
call get_command_argument(1,arg)
read(arg,*) strain
call get_command_argument(2,arg)
read(arg,*) seed
pi=acos(-1.0d0)
emod1=2.08d11
A1=pi*(26.0d-3)**2/4.0d0
L1=19.0d-3
vol1=L1*(530.66d-6)
sy=2.0d-3
coord=0.0d0
coord(2,1)=L1
disp=0.0d0
disp(2,1)=strain*L1
pforce22=0.0d0
pforce22(2,2)=seed/vol1
pforce22(2,3)=2.0d0*seed/vol1
totelement1=1
nodes(1,1)=1
nodes(1,2)=2
''' + extract_stiffness(source)+extract_rod(source)+'''
do i=1,2
write(*,'(3(ES26.17E3,1X))') pforce22(i,:)*vol1
enddo
end program native_rod
'''

def interface_program(source):
    return '''program native_interface
implicit none
integer i,j,cnode1,cnode2,cnode3,cnode4,cnode5,cnode6,cnode7,cnode8,m
real*8 pforce33(9,3),N1(1,1),N2(1,1),N3(1,1),N4(1,1),N5(1,1),N6(1,1),N7(1,1),N8(1,1)
real*8 dforce1,dforce2,dforce3,vol,vol1,nodevol
i=1
j=1
cnode1=2
cnode2=3
cnode3=4
cnode4=5
cnode5=6
cnode6=7
cnode7=8
cnode8=9
vol=(19.0d-3)**3
vol1=19.0d-3*530.66d-6
! These positive trilinear weights sum to one. The transfer block itself is original.
N1=.01d0
N2=.03d0
N3=.06d0
N4=.10d0
N5=.15d0
N6=.18d0
N7=.22d0
N8=.25d0
pforce33=0.0d0
dforce1=1.0d9
dforce2=-2.0d9
dforce3=3.0d9
''' + extract_interface(source)+'''
do m=1,9
nodevol=vol
if(m==1)nodevol=vol1
write(*,'(3(ES26.17E3,1X))')pforce33(m,:)*nodevol
enddo
end program native_interface
'''

def extract_contact_reaction(source):
    start = source.index('            vel(i,3)=(disp(i,3)-olddisp(i,3))/ddt')
    start = source.index('\n',start)+1
    end = source.index('            pentfx=pentfx+dpenfx',start)
    return source[start:end]

def contact_program(source):
    return """program native_contact
implicit none
integer i,point(1)
real*8 dens,dengs,vol,vol1,ddt,vel(1,3),oldvelx,oldvelz,dpenfx,dpenfy,dpenfz
character(32) arg
call get_command_argument(1,arg)
read(arg,*)point(1)
i=1
dens=2500.0d0
dengs=7850.0d0
vol=(19.0d-3)**3
vol1=19.0d-3*530.66d-6
ddt=3.0d-6
vel(1,1)=1.0d0
vel(1,2)=0.0d0
vel(1,3)=2.0d0
oldvelx=0.2d0
oldvelz=0.3d0
"""+extract_contact_reaction(source)+"""
write(*,'(3(ES26.17E3,1X))')dpenfx,dpenfy,dpenfz
end program native_contact
"""

def extract_support(source):
    start=source.index('    do i=1,totnode\n        if((abs(coord(i,1)+10.5d0*dx)')
    end=source.index('        dispdispdisp(i)=',start)
    return source[start:end].rsplit('    do i=1,totnode',1)[0]

def support_program(source):
    return """program native_support
implicit none
integer i,totnode
real*8 dx,coord(4,3),disp(4,3),vel(4,3)
dx=19.0d-3
totnode=4
coord=0.0d0
coord(1,1)=-10.5d0*dx
coord(2,1)=10.5d0*dx
coord(3,1)=-81.5d0*dx
coord(4,1)=81.5d0*dx
coord(1,3)=8.5d0*dx
coord(2,3)=8.5d0*dx
coord(3,3)=-8.5d0*dx
coord(4,3)=-8.5d0*dx
disp=1.0d0
vel=2.0d0
"""+extract_support(source)+"""
do i=1,totnode
write(*,'(6(ES26.17E3,1X))')disp(i,:),vel(i,:)
enddo
end program native_support
"""

def extract_first_contact_loop(source):
    start=source.index('    do i=1,totnode\n        newpx=coord(i,1)+disp(i,1)',source.index('do tt=1,nt'))
    end=source.index('    ppenccx=ppenicx+ppenux',start)
    return source[start:end]

def full_contact_program(source):
    return """program native_full_contact
implicit none
integer i,totnode,point(1)
real*8 coord(1,3),disp(1,3),olddisp(1,3),vel(1,3)
real*8 dens,dengs,vol,vol1,ddt,penrad,penccx,penccy,penccz
real*8 newpx,newpy,newpz,dpenx,dpeny,dpenz,dpend,dpux,dpuy,dpuz,dprx,dpry,dprz,dpcx,dpcy,dpcz
real*8 olddispx,olddispy,olddispz,oldvelx,oldvely,oldvelz,dpenfx,dpenfy,dpenfz,pentfx,pentfy,pentfz
character(32) arg
call get_command_argument(1,arg)
read(arg,*)point(1)
totnode=1
dens=2500.0d0
dengs=7850.0d0
vol=(19.0d-3)**3
vol1=19.0d-3*530.66d-6
ddt=3.0d-6
penrad=0.005d0
penccx=0.0d0
penccy=0.0d0
penccz=0.0d0
coord=0.0d0
coord(1,2)=0.013d0
coord(1,3)=0.001d0
disp=0.0d0
disp(1,1)=0.0001d0
disp(1,2)=0.007d0
disp(1,3)=0.0001d0
olddisp=0.0d0
olddisp(1,2)=0.007d0
vel(1,1)=0.3d0
vel(1,2)=0.4d0
vel(1,3)=0.6d0
pentfx=0.0d0
pentfy=0.0d0
pentfz=0.0d0
"""+extract_first_contact_loop(source)+"""
write(*,'(9(ES26.17E3,1X))')disp(1,:),vel(1,:),pentfx,pentfy,pentfz
end program native_full_contact
"""

class Operators:
    def __init__(self,folder):
        self.folder=Path(folder)
        self.rod_original=compile_native(rod_program(SOURCE),'rod_original',self.folder)
        self.rod_corrected=compile_native(rod_program(apply_corrections(SOURCE,['steel_components'])),'rod_corrected',self.folder)
        self.interface_original=compile_native(interface_program(SOURCE),'interface_original',self.folder)
        self.interface_corrected=compile_native(interface_program(apply_corrections(SOURCE,['interface_balance'])),'interface_corrected',self.folder)
        self.contact_original=compile_native(contact_program(SOURCE),'contact_original',self.folder)
        self.contact_corrected=compile_native(contact_program(apply_corrections(SOURCE,['contact_mass'])),'contact_corrected',self.folder)
        self.support_original=compile_native(support_program(SOURCE),'support_original',self.folder)
        self.support_corrected=compile_native(support_program(apply_corrections(SOURCE,['support_velocity'])),'support_corrected',self.folder)
        self.fullcontact_original=compile_native(full_contact_program(SOURCE),'fullcontact_original',self.folder)
        self.fullcontact_corrected=compile_native(full_contact_program(apply_corrections(SOURCE,['contact_mass','contact_tangent'])),'fullcontact_corrected',self.folder)
    def rod(self,corrected,strain,seed=0):
        executable=self.rod_corrected if corrected else self.rod_original
        output=subprocess.check_output([str(executable),str(strain),str(seed)],text=True)
        return np.fromstring(output,sep=' ').reshape(2,3)
    def interface(self,corrected):
        executable=self.interface_corrected if corrected else self.interface_original
        output=subprocess.check_output([str(executable)],text=True)
        return np.fromstring(output,sep=' ').reshape(9,3)

    def contact(self,corrected,point):
        executable=self.contact_corrected if corrected else self.contact_original
        output=subprocess.check_output([str(executable),str(point)],text=True)
        return np.fromstring(output,sep=' ')
    def support(self,corrected):
        executable=self.support_corrected if corrected else self.support_original
        output=subprocess.check_output([str(executable)],text=True)
        return np.fromstring(output,sep=' ').reshape(4,6)

    def full_contact(self,corrected,point):
        executable=self.fullcontact_corrected if corrected else self.fullcontact_original
        output=subprocess.check_output([str(executable),str(point)],text=True)
        return np.fromstring(output,sep=' ')

@pytest.fixture(scope='module')
def operators(tmp_path_factory):
    return Operators(tmp_path_factory.mktemp('native_fortran_operators'))

def test_native_yielded_rod_internal_force_balance(operators):
    before=operators.rod(False,.003)
    after=operators.rod(True,.003)
    np.testing.assert_allclose(before.sum(axis=0),[0,-CAP,-CAP],rtol=1e-13,atol=1e-8)
    np.testing.assert_allclose(after.sum(axis=0),0,atol=1e-8)
    np.testing.assert_allclose(after[:,0],before[:,0],atol=0)

def test_native_corrected_rod_preserves_existing_transverse_accumulator(operators):
    seed=17.
    before=operators.rod(False,.003,seed)
    after=operators.rod(True,.003,seed)
    np.testing.assert_allclose(after[1,1:],[seed,2*seed],rtol=1e-14)
    assert before[1,1] < -2e5 and before[1,2] < -2e5

def test_native_pre_yield_and_compression_are_unchanged(operators):
    for strain in [.001,-.001]:
        before=operators.rod(False,strain,17)
        after=operators.rod(True,strain,17)
        np.testing.assert_array_equal(before,after)

def test_native_interface_action_reaction_and_virtual_work(operators):
    before=operators.interface(False)
    after=operators.interface(True)
    density=np.array([1e9,-2e9,3e9])
    np.testing.assert_allclose(before.sum(axis=0),density*(V-VC),rtol=1e-13)
    np.testing.assert_allclose(after.sum(axis=0),0,atol=1e-9)
    np.testing.assert_array_equal(after[0],before[0])
    weights=np.array([.01,.03,.06,.10,.15,.18,.22,.25])
    velocities=np.random.default_rng(491).normal(size=(9,3))
    expected=np.dot(after[0],velocities[0]-weights@velocities[1:])
    power=float(np.sum(after*velocities))
    np.testing.assert_allclose(power,expected,rtol=1e-13)

def test_native_contact_reaction_uses_contacted_node_mass(operators):
    concrete=operators.contact(False,0)
    np.testing.assert_array_equal(concrete,operators.contact(True,0))
    expected_concrete=-2500*VC*np.array([.8,0,1.7])/3e-6
    expected_steel=-7850*V*np.array([.8,0,1.7])/3e-6
    np.testing.assert_allclose(concrete,expected_concrete,rtol=1e-14)
    np.testing.assert_allclose(operators.contact(True,1),expected_steel,rtol=1e-14)
    np.testing.assert_array_equal(operators.contact(False,1),concrete)

def test_native_support_constraints_reset_only_constrained_velocity(operators):
    before=operators.support(False)
    after=operators.support(True)
    np.testing.assert_array_equal(before[:,:3],after[:,:3])
    np.testing.assert_array_equal(before[:,3:],2)
    expected=np.full((4,3),2.)
    expected[:2,0]=0
    expected[2:,2]=0
    np.testing.assert_array_equal(after[:,3:],expected)

def test_native_frictionless_cylinder_preserves_tangential_motion_and_mass_reaction(operators):
    for point in [0,1]:
        before=operators.full_contact(False,point)
        after=operators.full_contact(True,point)
        assert before[1]==0 and before[4]==0
        assert after[1]==.007 and after[4]==.4
        np.testing.assert_array_equal(before[[0,2,3,5]],after[[0,2,3,5]])
        projected_xz=np.array([after[0],.001+after[2]])
        assert np.isclose(np.linalg.norm(projected_xz),.005,rtol=1e-14)
        actual_mass=7850*V if point else 2500*VC
        expected=-actual_mass*(after[[3,5]]-np.array([.3,.6]))/3e-6
        np.testing.assert_allclose(after[[6,8]],expected,rtol=1e-13)
        assert after[7]==0

if __name__=='__main__':
    folder=Path(__file__).parent/'native_operator_build'
    folder.mkdir(exist_ok=True)
    op=Operators(folder)
    tests=[test_native_yielded_rod_internal_force_balance,
           test_native_corrected_rod_preserves_existing_transverse_accumulator,
           test_native_pre_yield_and_compression_are_unchanged,
           test_native_interface_action_reaction_and_virtual_work,
           test_native_contact_reaction_uses_contacted_node_mass,
           test_native_support_constraints_reset_only_constrained_velocity,
           test_native_frictionless_cylinder_preserves_tangential_motion_and_mass_reaction]
    for test in tests:test(op)
    record={'original_bytes_sha256':hashlib.sha256((ROOT/'Console2.original.f90').read_bytes()).hexdigest(),
            'normalized_source_sha256':hashlib.sha256(SOURCE.encode()).hexdigest(),
            'passed_native_compiled_tests':[t.__name__ for t in tests],
            'original_yielded_rod_physical_force_N':op.rod(False,.003).tolist(),
            'corrected_yielded_rod_physical_force_N':op.rod(True,.003).tolist(),
            'original_interface_net_force_N':op.interface(False).sum(axis=0).tolist(),
            'corrected_interface_net_force_N':op.interface(True).sum(axis=0).tolist(),
            'original_contact_steel_reaction_N':op.contact(False,1).tolist(),
            'corrected_contact_steel_reaction_N':op.contact(True,1).tolist(),
            'original_full_contact_steel_disp_vel_reaction':op.full_contact(False,1).tolist(),
            'corrected_full_contact_steel_disp_vel_reaction':op.full_contact(True,1).tolist(),
            'build_directory':str(folder)}
    (Path(__file__).parent/'native_operator_verification.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))
