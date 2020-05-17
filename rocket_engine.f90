module fun
  implicit none
contains
  !paylord ratio
  function payload_r(dv,vj,fi) result(pr)
    real,intent(in)::dv,vj,fi
    real::pr
    pr = (exp(-dv/vj)-fi)/(1-fi)
  end function payload_r

  !exhaust velocity
  function exhaust_v(ita,gamma,R,M,Tf,pa,po) result(vj)
    real,intent(in)::ita,gamma,R,M,Tf,pa,po
    real::vj,medium,gt
    gt =gamma/(gamma-1)
    medium = (1 -(pa/po)**(1/gt))
    vj=sqrt(2*ita* gt * R * Tf * medium / M )
  end function exhaust_v

  !construction ratio
  function const_r(fis,rho,rhos) result(fi)
    real,intent(in)::fis,rho,rhos
    real::fi
    fi=(1/fis -1)*rho/rhos +1
    fi=1/fi
  end function const_r

  !temperature in combustion chamber
  function temp(no,hi,ho,nio,gamma,R,fini,To,MR,Mi) result(Tf)
    real,intent(in)::gamma, R, To, MR
    integer,intent(in)::fini
    real,intent(in)::no(fini),hi(2),ho(fini),nio(2),Mi(2)
    real::Tf,Q,k
    integer::i

    !reactants h
    Q=(1-nio(1))*hi(1)
    Q=(MR*Mi(1)/Mi(2) - nio(2))*hi(2)+Q

    !products h
     do i = 1, fini
        Q=Q-no(fini)*ho(fini)
     end do

     !solve Tf
     k=(sum(no)+sum(nio))*gamma*R/(gamma-1)
     Tf=Q/k +To
     if(Tf > 2500)then
        Tf = (2*Q+To)/(0.0004*Q+k)
     else
        Tf=Tf
     end if
   end function temp
   
   !dencity of propellant
   function dencity(rhoi,MR) result(rho)
     real,intent(in)::MR
     real,intent(in)::rhoi(2)
     real::rho
     rho=(1+MR)/(1/rhoi(1)+MR/rhoi(2))
   end function dencity

   !exhaust molecule
   function molecule(nio,no,Mi,Mo,fini) result(M)
     integer,intent(in)::fini
     real,intent(in)::nio(2),no(fini),Mi(2),Mo(fini)
     real::M,k
     integer::i
     k=nio(1)*Mi(1)+nio(2)*Mi(2)
     do i=1, fini
        k=k+no(i)*Mo(i)
     end do
     M=k/(sum(no)+sum(nio))
   end function molecule

   !specific thrust
   function spe_thr(vj,g) result(isp)
     real,intent(in)::vj,g
     real::isp
     isp=vj/g
   end function spe_thr

   !area of throat
   function ar_th(F,cf,po) result(at)
     real,intent(in)::F,cf,po
     real::at
     at=F/(cf*po)
   end function ar_th

   !thrust
   function thr(ml,pr,aci) result(F)
     real,intent(in)::ml,pr,aci
     real::F
     F=ml*aci/pr
   end function thr

   !constant of thrust
   function con_thr(gamma, pa, po) result(cf)
     real,intent(in)::gamma, pa, po
     real::cf
     cf=sqrt(2*gamma**2/(gamma-1) * (2/(gamma-1))**((gamma+1)/(gamma-1))*(1 -(pa/po)**((gamma-1)/gamma)))
   end function con_thr

   !length of combustion chamber
   function len_com(ar,ts,M,gamma,R,Tf) result(lenc)
     real,intent(in)::ar,ts,M,gamma,R,Tf
     real::lenc,medium,gt
     gt=((gamma+1)/2)**((gamma+1)/2*(gamma-1))
     medium = sqrt(M*gt/(gamma*R*Tf))
     lenc=ts/(medium*ar)
   end function len_com

   !diameter of combustion chamber
   function dia_com(ar,at) result(cd)
     real,intent(in)::ar,at
     real::cd,pi
     pi=3.14159265
     cd=sqrt((4*ar*at)/pi)
   end function dia_com

   !mass of propellant
   function mas_pro(pr,fi,ml) result(mp)
     real,intent(in)::pr,fi,ml
     real::mp
     mp=(1/pr-1)*(1-fi)*ml
   end function mas_pro

   !length of tank
   function len_tan(mp,rho,dtan) result(ltan)
     real,intent(in)::mp,rho,dtan
     real::ltan,pi
     pi=3.14159265
     ltan=4*mp/(rho*pi*dtan*dtan*1000)
     !単位換算で1000かける
   end function len_tan

   !consider about chemical things
   function chemino(ci,co,MR,Mi,fini) result(no)
     integer,intent(in)::fini
     real,intent(in)::MR,Mi(2),ci(2),co(fini)
     real::nio(2),no(fini)
     real::nox
     integer::i
     nox=MR*Mi(1)/Mi(2)
     if(ci(2)/ci(1)>=nox)then
        !fuel rich
        nio(1)=1-ci(1)/ci(2) * nox
        nio(2)=0
        do i = 1, fini
           no(i) = co(i)/ci(2) * nox
        end do
     else
        !fuel lean
        nio(1)=0
        nio(2)=nox-ci(2)/ci(1)
        do i =1, fini
           no(i) = co(i)/ci(1)
        end do
     end if
   end function chemino
   
   function cheminio(ci,co,MR,Mi,fini) result(nio)
     integer,intent(in)::fini
     real,intent(in)::MR,Mi(2),ci(2),co(fini)
     real::nio(2),no(fini)
     real::nox
     integer::i
     nox=MR*Mi(1)/Mi(2)
     if(ci(2)/ci(1)>=nox)then
        !fuel rich
        nio(1)=1-ci(1)/ci(2) * nox
        nio(2)=0
        do i = 1, fini
           no(i) = co(i)/ci(2) * nox
        end do
     else
        !fuel lean
        nio(1)=0
        nio(2)=nox-ci(2)/ci(1)
        do i =1, fini
           no(i) = co(i)/ci(1)
        end do
     end if
   end function cheminio
  end module fun
 

 program engine
   use fun
   implicit none
   integer::fini
   integer::j
   real::g,R,pa,ita,gamma,fis,rhos,To,MR,ci(2)
   real::ml,mp,pr,dv,vj,fi,rho,po,Tf,isp,F,at,ar,cf,aci,ts,lenc,cd,ltan,dtan,M
   real,allocatable::no(:),ho(:),Mo(:),co(:)
   real::nio(2),hi(2),Mi(2),rhoi(2)

   !初期条件
   g=9.81
   R=8.314
   pa=10000
   ita=0.96
   gamma=1.4
   fis=0.1
   rhos=0.3
   To=300
   ml=2000
   dv=7000
   aci=1.3*g
   ts=0.010
!反応物質について
   !燃料密度[g/cm^3]
   rhoi(1)=0.071
   !燃料の質量数[kg/mol]
   Mi(1)=0.002
   !燃料のエンタルピー[J/mol]
   hi(1)=0
   !燃料の完全反応の時の係数(整数）
   ci(1)=2
   
   !酸化剤密度[g/com^3]
   rhoi(2)=1.140
   !酸化剤の質量数[kg/mol]
   Mi(2)=0.032
   !酸化剤のエンタルピー[J/mol]
   hi(2)=0
   !酸化剤の完全反応の時の係数(整数）
   ci(2)=1
   
   !排気ガスの数（燃料のあまりを除く）（finiの個数分だけMo,hi,coは増やせます
   fini=1

   allocate(no(fini),ho(fini),Mo(fini),co(fini))
   !排気ガスの質量数[kg/mol]
   Mo(1)=0.018
   !酸化剤のエンタルピー[J/mol]
   ho(1)=-241900
   !酸化剤の完全反応の時の係数(整数）
   co(1)=2

 !ロケットそのものについて
   !燃焼室とスロート面積の比
   ar=5
   !燃焼室圧力
   po=12000000
   !タンク直径
   dtan = 2.0

   !計算開始
   do j = 1, 100
     MR = 0.1* j
     no = chemino(ci,co,MR,Mi,fini)
     nio =  cheminio(ci,co,MR,Mi,fini)
     rho = dencity(rhoi,MR)
     fi = const_r(fis,rho,rhos)
     Tf = temp(no,hi,ho,nio,gamma,R,fini,To,MR,Mi)
     M = molecule(nio,no,Mi,Mo,fini)
     vj = exhaust_v(ita,gamma,R,M,Tf,pa,po)
     pr = payload_r(dv,vj,fi)
     isp = spe_thr(vj,g)
     F = thr(ml,pr,aci)
     cf = con_thr(gamma, pa, po)
     at = ar_th(F,cf,po)
     lenc = len_com(ar,ts,M,gamma,R,Tf)
     cd = dia_com(ar,at)
     mp = mas_pro(pr,fi,ml)
     ltan = len_tan(mp,rho,dtan)
     print*, pr, po, MR, Tf, isp, cd, lenc, dtan, ltan
  end do
 end program engine
 
      
      
      
