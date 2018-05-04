!!written by xiong wei
!!revised by Jiao shuping
module global
implicit none
!============================ 用户可能需要更新的量======================
integer,parameter::no_frames=1700! 轨迹文件的总帧数， 包括第一帧
integer,parameter::fileunit=20
real(8),parameter::equitime=1400.0 ! 系统平衡时间， 单位ps
real(8),parameter::timestep=0.001 !单位ps
real(8),parameter::PI=3.14159265
character(14),save::inistru_filename= 'water-graphene'
character(27)::strj_filename= 'atomacc490.lammpstrj'
integer::col,row !density_first_water 画密度云图网格
real(8),parameter:: grid=0.3 !Angstroms
integer,parameter::splits=36 !-90—90 度分多少份统计水分子个数
integer::count_no(splits)
real(8),parameter::kedu=180.0/splits
real(8)::no_hbonds_all ! 第一层水中平均每个水分子具有的氢键数目
real(8)::no_hbonds_in !D和A 水分子均在第一层中
real(8)::no_hbonds_out !D和A 水分子至少有一个不在第一层中
!RDFs
real(8),parameter::delta=0.05 !单位A
real(8),parameter::cutoff=10.0
real(8),allocatable::g_global(:)
integer::ceng
real(8)::rho !RDF 中平均密度
real(8),parameter::z_lo=-9.39877
real(8),parameter::z_hi=-2.8346
real(8),parameter::CC=1.42*1.0
!======================================================================
real(8)::q1
real(8)::q2
real(8)::qcos
real(8)::qsin
real(8)::waters ! 第一层水的数目
real(8)::struc_factor1,struc_factor2,struc_factor3
integer,parameter::n_q=45
integer,parameter::n_thelta=72
real(8)::struc_factor(0:n_q,0:n_thelta)
real(8)::struc_ave(0:n_q)
type atom_type
!!!!
integer::molecule_no ! 分子编号
integer::species ! 原子种类
real(8)::mass !质量
real(8)::charge !电荷
real(8)::position(3) ! 原子坐标
end type atom_type
type water_type
integer::oxygen_no ! 某水分子中氧原子编号
integer::hydrogen_no1 ! 氢原子编号
integer::hydrogen_no2 ! 氢原子编号
end type water_type
integer,save::no_atoms,no_bonds,no_angles ! 系统中原子、键、键角数目
integer,save::no_waters ! 水分子个数， 本程序中等于键角数目
integer,save::no_atom_types,no_bond_types,no_angle_types ! 原子种类的数目
real(8),save::time
real(8),save::xlo,xhi,ylo,yhi,zlo,zhi ! 模拟盒子大小
type(atom_type),allocatable,save::atom(:) ! 每个原子
type(water_type),allocatable,save::water(:) ! 每个水分子
real(8),allocatable,save::masses(:) ! 不同种类原子的质量
integer,allocatable,save::density(:,:)
end module
! 主程序， 后处理分析的入口。通过输入lammps 的初始构型文件初始化原子系统
!使用atom type 的数据结构存储所有原子信息
program water_structure
use global
implicit none
integer::frame,i,j
integer::eff_frames=0 ! 用于统计的有效帧数
character(len=50)::filename
open(99,file='test.txt',status='unknown')
write(99,*)'oooo'
call initial()
write(99,*)'oooo'
! 更新每一帧的原子坐标
open(fileunit,file=strj_filename,status='unknown')
frame=1
do while(frame<=no_frames)
call update_position()
write(*,*)'oooo'
if(time>equitime) then
call countwater()
write(*,*)'oooo'
eff_frames=eff_frames+1
call RDFs()
do i=0,n_q
do j=0,n_thelta
!!!!!
!1/A,here checked right 20110325
qcos=(1.0+DBLE(i)/DBLE(n_q)*(6.0-1.0))* &
cos(DBLE(j)/DBLE(n_thelta)*2.0*pi)
qsin=(1.0+DBLE(i)/DBLE(n_q)*(6.0-1.0))* &
sin(DBLE(j)/DBLE(n_thelta)*2.0*pi)
call Struc(i,j)
enddo
enddo
call density_first_water()
endif
frame=frame+1
enddo
open(333,file='goo_global.txt',status='unknown')
rho=rho/eff_frames
write(*,*) eff_frames,rho*(xhi-xlo)*(yhi-ylo)*eff_frames
write(333,*) rho/(z_hi-z_lo)*1000.0/33.3 ! 用块体水密度归一化
do i=1,ceng
write(333,*) (i-0.5)*delta,g_global(i)/eff_frames/rho
enddo
open(222,file='struc_factor.txt',status='unknown')
open(333,file='struc_factor_matlab.txt',status='unknown')
open(555,file='struc_ave.txt',status='unknown')
struc_factor=struc_factor/eff_frames
!用gnuplot 画图
do i=0,n_q
do j=0,n_thelta
write(222,*) 10.0*(1.0+DBLE(i)/DBLE(n_q)*(6.0-1.0)), &
DBLE(j)/DBLE(n_thelta)*360.0,Struc_Factor(i,j)
struc_ave(i)=struc_ave(i)+Struc_Factor(i,j)
enddo
write(222,*)
struc_ave(i)=struc_ave(i)/DBLE(n_thelta+1)
write(555,*) 10.0*(1.0+DBLE(i)/DBLE(n_q)*(6.0-1.0)),struc_ave(i)
enddo
!用matlab 画图
do i=n_q,0,-1
write(333,"(<n_thelta+1>f12.4)") Struc_Factor(i,:)
enddo
close(333)
! 进行统计分析
open(50,FILE='density_first_water.plt',STATUS='unknown')
write(50,*) 'TITLE = ”density_first_water 2D Plot”'
write(50,*) 'VARIABLES = ”X”, ”Y”, ”density”'
write(50,*) 'ZONE T=”BIG ZONE”, I=’,row,’ J=’,col,’F=POINT'
!!!!!!!
do i=1,col
do j=1,row
! write(50,’(2i5,2f9.4,i5)’) &
!i,j,real(i)/10.0-0.05,real(j)/10.0-0.05,density(i,j)
write(50,'(3f9.4)') (real(i)-0.5)*grid, &
(real(j)-0.5)*grid,REAL(density(i,j))/REAL(eff_frames)
enddo
enddo
close(50)
!
end program water_structure
subroutine initial()
use global
implicit none
integer::atom_no ! 原子编号
integer::water_no ! 水分子编号
integer::i,temp_int
integer::step ! 输出的步数
LOGICAL(4):: results
open(10,file=inistru_filename,status='unknown')
read(10,*)
read(10,*)
read(10,*) no_atoms
read(10,*) no_bonds
read(10,*) no_angles
no_waters=no_angles ! 水分子个数， 本程序中等于键角数目
read(10,*)
read(10,*) no_atom_types
read(10,*) no_bond_types
read(10,*) no_angle_types
read(10,*)
read(10,*) xlo,xhi
read(10,*) ylo,yhi
read(10,*) zlo,zhi
read(10,*)
read(10,*)
read(10,*)
col=floor((xhi-xlo)/grid)+1
row=floor((yhi-ylo)/grid)+1
write(*,*) "meshnum"
write(*,*) col,row
allocate(masses(no_atom_types))
allocate(atom(no_atoms))
allocate(water(no_waters))
!!!!!!
allocate(density(col,row))
ceng=ceiling(cutoff/delta) !RDFs 中的层数
allocate(g_global(ceng))
g_global=0.0
rho=0.0
density=0
i=1
do while(i<=no_atom_types)
read(10,*) temp_int,masses(i)
i=i+1
enddo
!read(10,*)
!read(10,*)
!read(10,*)
!i=1
!do while(i<=no_bond_types)
!read(10,*)
!i=i+1
!enddo
!read(10,*)
!read(10,*)
!read(10,*)
!i=1
!do while(i<=no_angle_types)
!read(10,*)
!i=i+1
!enddo
read(10,*)
read(10,*)
read(10,*)
i=1
do while(i<=no_atoms)
read(10,*) atom_no,atom(atom_no)%molecule_no,atom(atom_no)%species,&
atom(atom_no)%charge,atom(atom_no)%position
atom(atom_no)%mass=masses(atom(atom_no)%species)
if(atom_no/=i) print *,'numbers of atoms are disordered!'
i=i+1
enddo
read(10,*)
read(10,*)
read(10,*)
i=1
do while(i<=no_bonds)
read(10,*)
!!!!!!
i=i+1
enddo
read(10,*)
read(10,*)
read(10,*)
i=1
do while(i<=no_waters)
read(10,*) water_no,temp_int,water(water_no)%hydrogen_no1,&
water(water_no)%oxygen_no,water(water_no)%hydrogen_no2
if(water_no/=i) print *,'numbers of waters(angles) are disordered!'
i=i+1
enddo
no_hbonds_all=0.0 ! 第一层水中平均每个水分子具有的氢键数目
no_hbonds_in=0.0 ! D和A 水分子均在第一层中
no_hbonds_out=0.0
q1=1.5*CC
q2=0.5*SQRT(3.0)*CC
!q1=0.0
!q2=SQRT(3.0)*CC
struc_factor1=0.0
struc_factor2=0.0
struc_factor3=0.0
struc_factor=0.0
struc_ave=0.0
close(10)
end subroutine
subroutine update_position()
use global
implicit none
integer::i,temp_int
integer::atom_no ! 原子编号
integer::step ! 输出的步数
read(fileunit,*)
read(fileunit,*) step
read(fileunit,*)
read(fileunit,*) ! 该行是原子个数， 不用更新
read(fileunit,*)
read(fileunit,*) xlo,xhi
read(fileunit,*) ylo,yhi
read(fileunit,*) zlo,zhi
read(fileunit,*)
time=timestep*step !单位ps, 根据lammps 中对单位的定义可能不同
!!!!!
i=1
do while(i<=no_atoms)
! 更新原子坐标
read(fileunit,*) atom_no,temp_int,atom(atom_no)%position(:)
! 由于轨迹文件的特点， 需要化为绝对坐标
atom(atom_no)%position(1)=xlo+atom(atom_no)%position(1)*(xhi-xlo)
atom(atom_no)%position(2)=ylo+atom(atom_no)%position(2)*(yhi-ylo)
atom(atom_no)%position(3)=zlo+atom(atom_no)%position(3)*(zhi-zlo)
i=i+1
enddo
end subroutine
! 通过分析lammps 的轨迹输出文件， 统计靠近壁面的第一层水分子的数目
! 第一层水分子的范围划定是经验性质的， 即距离壁面3.0<d<=3.5 Angstroms
subroutine countwater()
use global
implicit none
integer::i,o_no,n
real(8)::x,y,z
n=0
do i=1,no_waters
o_no=water(i)%oxygen_no
if(atom(o_no)%position(3)>z_lo .and. &
atom(o_no)%position(3)<=z_hi) then
n=n+1
endif
if(atom(o_no)%species==2) then
else
write(*,*) atom(o_no)%species
endif
enddo
waters=n
end subroutine
! 通过分析lammps 的轨迹文件， 统计靠近壁面的第一层水分子（ 其实是氧原子） 的
! 二维面内径向分布函数以及第一层水的结构因子
! 第一层水分子的范围划定是经验性质的， 即距离壁面3.0<d<=3.5 Angstroms
subroutine RDFs()
use global
implicit none
integer::i,j,o_no,n
integer,allocatable::onos(:)
!!!!!!
real(8)::dx,dy,r,factor,x,y,z
real(8),allocatable::g(:)
real(8)::sf1_cos
real(8)::sf1_sin
real(8)::sf2_cos
real(8)::sf2_sin
real(8)::sf3_cos
real(8)::sf3_sin
n=waters
allocate(onos(n))
allocate(g(ceng))
g=0.0
onos=0
sf1_cos=0.0
sf1_sin=0.0
sf2_cos=0.0
sf2_sin=0.0
sf3_cos=0.0
sf3_sin=0.0
j=1
do i=1,no_waters
o_no=water(i)%oxygen_no
if(atom(o_no)%position(3)>z_lo .and. &
atom(o_no)%position(3)<=z_hi) then
onos(j)=o_no
j=j+1
endif
enddo
dx=xhi-xlo
dy=yhi-ylo
do i=1,n
!先求structure factor
x=atom(onos(i))%position(1)
y=atom(onos(i))%position(2)
z=atom(onos(i))%position(3)
sf1_cos=sf1_cos+cos((x/q1+y/q2)*pi)
sf1_sin=sf1_sin+sin((x/q1+y/q2)*pi)
sf2_cos=sf2_cos+cos((x/q1-y/q2)*pi)
sf2_sin=sf2_sin+sin((x/q1-y/q2)*pi)
sf3_cos=sf3_cos+cos((2.0*x/q1)*pi)
sf3_sin=sf3_sin+sin((2.0*x/q1)*pi)
!!!!
! 以下求RDF
do j=1,n
! 原区域
r=sqrt((atom(onos(i))%position(1)-atom(onos(j))%position(1))**2 &
+(atom(onos(i))%position(2)-atom(onos(j))%position(2))**2+ &
(atom(onos(i))%position(3)-atom(onos(j))%position(3))**2)
if(r>delta/10.0 .and. r<=cutoff) then
g(ceiling(r/delta))=g(ceiling(r/delta))+1.0
endif
!-0
r=sqrt((atom(onos(i))%position(1)-atom(onos(j))%position(1)+dx)**2 &
+(atom(onos(i))%position(2)-atom(onos(j))%position(2))**2+ &
(atom(onos(i))%position(3)-atom(onos(j))%position(3))**2)
if(r>delta/10.0 .and. r<=cutoff) then
g(ceiling(r/delta))=g(ceiling(r/delta))+1.0
endif
!+0
r=sqrt((atom(onos(i))%position(1)-atom(onos(j))%position(1)-dx)**2 &
+(atom(onos(i))%position(2)-atom(onos(j))%position(2))**2+ &
(atom(onos(i))%position(3)-atom(onos(j))%position(3))**2)
if(r>delta/10.0 .and. r<=cutoff) then
g(ceiling(r/delta))=g(ceiling(r/delta))+1.0
endif
!0-
r=sqrt((atom(onos(i))%position(1)-atom(onos(j))%position(1))**2 &
+(atom(onos(i))%position(2)-atom(onos(j))%position(2)+dy)**2+ &
(atom(onos(i))%position(3)-atom(onos(j))%position(3))**2)
if(r>delta/10.0 .and. r<=cutoff) then
g(ceiling(r/delta))=g(ceiling(r/delta))+1.0
endif
!0+
r=sqrt((atom(onos(i))%position(1)-atom(onos(j))%position(1))**2 &
+(atom(onos(i))%position(2)-atom(onos(j))%position(2)-dy)**2+ &
(atom(onos(i))%position(3)-atom(onos(j))%position(3))**2)
if(r>delta/10.0 .and. r<=cutoff) then
g(ceiling(r/delta))=g(ceiling(r/delta))+1.0
endif
!-+
r=sqrt((atom(onos(i))%position(1)-atom(onos(j))%position(1)+dx)**2 &
+(atom(onos(i))%position(2)-atom(onos(j))%position(2)-dy)**2+ &
(atom(onos(i))%position(3)-atom(onos(j))%position(3))**2)
if(r>delta/10.0 .and. r<=cutoff) then
g(ceiling(r/delta))=g(ceiling(r/delta))+1.0
!!!!!!
endif
!--
r=sqrt((atom(onos(i))%position(1)-atom(onos(j))%position(1)-dx)**2 &
+(atom(onos(i))%position(2)-atom(onos(j))%position(2)+dy)**2+ &
(atom(onos(i))%position(3)-atom(onos(j))%position(3))**2)
if(r>delta/10.0 .and. r<=cutoff) then
g(ceiling(r/delta))=g(ceiling(r/delta))+1.0
endif
!+-
r=sqrt((atom(onos(i))%position(1)-atom(onos(j))%position(1)+dx)**2 &
+(atom(onos(i))%position(2)-atom(onos(j))%position(2)+dy)**2+ &
(atom(onos(i))%position(3)-atom(onos(j))%position(3))**2)
if(r>delta/10.0 .and. r<=cutoff) then
g(ceiling(r/delta))=g(ceiling(r/delta))+1.0
endif
!++
r=sqrt((atom(onos(i))%position(1)-atom(onos(j))%position(1)-dx)**2 &
+(atom(onos(i))%position(2)-atom(onos(j))%position(2)-dy)**2+ &
(atom(onos(i))%position(3)-atom(onos(j))%position(3))**2)
if(r>delta/10.0 .and. r<=cutoff) then
g(ceiling(r/delta))=g(ceiling(r/delta))+1.0
endif
enddo
enddo
struc_factor1=struc_factor1+(sf1_cos**2+sf1_sin**2)/n
struc_factor2=struc_factor2+(sf2_cos**2+sf2_sin**2)/n
struc_factor3=struc_factor3+(sf3_cos**2+sf3_sin**2)/n
rho=rho+n/dx/dy
do i=1,ceng
factor=2.0*pi*delta*(i-0.5)*delta
g(i)=g(i)/n/factor
g_global(i)=g_global(i)+g(i)
enddo
deallocate(onos)
deallocate(g)
end subroutine
! 通过分析lammps 的轨迹输出文件， 统计靠近壁面的第一层水的structure factor
! 第一层水分子的范围划定是经验性质的， 即距离壁面3.0<d<=3.5 Angstroms
!!!!!
subroutine Struc(k,l)
use global
implicit none
integer,allocatable::onos(:)
integer::i,j,o_no,n,k,l
real(8)::x,y,z
real(8)::sf_cos
real(8)::sf_sin
n=waters
allocate(onos(n))
sf_cos=0.0
sf_sin=0.0
j=1
do i=1,no_waters
o_no=water(i)%oxygen_no
if(atom(o_no)%position(3)>z_lo .and. &
atom(o_no)%position(3)<=z_hi) then
onos(j)=o_no
j=j+1
endif
enddo
do i=1,n
!求structure factor
x=atom(onos(i))%position(1)
y=atom(onos(i))%position(2)
z=atom(onos(i))%position(3)
sf_cos=sf_cos+cos(x*qcos+y*qsin)
sf_sin=sf_sin+sin(x*qcos+y*qsin)
enddo
struc_factor(k,l)=struc_factor(k,l)+(sf_cos**2+sf_sin**2)/n
end subroutine
! 通过分析lammps 的轨迹文件， 统计靠近壁面的第一层水分子的密度分布画出云图
! 第一层水分子的范围划定是经验性质的， 即距离壁面小于5.0 Angstroms
subroutine density_first_water()
use global
implicit none
integer::i,j,ii,jj,o_no
real(8)::xo,yo,zo
i=1
!!!!!!
do while(i<=no_waters)
o_no=water(i)%oxygen_no
xo=atom(o_no)%position(1)
yo=atom(o_no)%position(2)
zo=atom(o_no)%position(3)
if(zo>z_lo .and. zo<=z_hi) then
if(xo<xlo) then
ii=((xhi+xo-2.0*xlo)/(xhi-xlo)*col)+1
elseif(xo>xhi) then
ii=((xo-xhi)/(xhi-xlo)*col)+1
elseif(xo==xlo) then
ii=1
elseif(xo==xhi) then
ii=col
else
ii=((xo-xlo)/(xhi-xlo)*col)+1
endif
if(yo<ylo) then
jj=((yhi+yo-2.0*ylo)/(yhi-ylo)*row)+1
elseif(yo>yhi) then
jj=((yo-yhi)/(yhi-ylo)*row)+1
elseif(yo==ylo) then
jj=1
elseif(yo==yhi) then
jj=row
else
jj=((yo-ylo)/(yhi-ylo)*row)+1
endif
density(ii,jj)=density(ii,jj)+1
endif
i=i+1
enddo
end subroutine