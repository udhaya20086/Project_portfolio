usemathjax: true
---
title: Welcome to my blog
---

Hello

The objective of this project is to solve 1D supersonic nozzle flow equations in conservative and non-conservative forms using the Macormack Method, to perform grid dependence test, to find the number of cycles required to be run to reach convergence, to compare the normalized mass flow obtained from both the forms and find out which form of nozzle flow equation is faster.

The 1D supersonic nozzle flow equations in Conservative form are given below:

1) Continuity equation:

$$(del(rho'A'))/(del(t'))+(del(rho'A'V'))/(del(x'))=0$$

2) Momentum equation:

`(del(rho'A'V'))/(del(t'))+(del(rho'A'V'^2+(1/gamma)p'A'))/(del(x'))=(1/gamma)p'(del(A'))/(del(x'))`

3) Energy equation:

`(del(rho'(((e')/(gamma-1))+(gamma/2v'^2))A'))/(del(t'))+(del(rho'(((e')/(gamma-1))+(gamma/2v'^2))v'A'+p'A'v'))/(del(x'))=0`

The 1D supersonic nozzle flow equations in Non-conservative form are given below:

1)Continuity equation:

`(del(rho'))/(del(t'))=-rho'((delV')/(delx'))-rho'v'(del(lnA'))/(del(x'))-v'(del(rho'))/(del(x'))`

2)Momentum equation:

`(del(V'))/(del(t'))=-v'(del(v'))/(del(x'))-(1/gamma)((delT')/(delx')+(T')/(rho')(del(rho'))/(del(x')))`

3) Energy equation:

`(del(T'))/(del(t'))=-v'(del(T'))/(del(x'))-(gamma-1)T'((del(V'))/(del(x'))+v'(del(lnA'))/(del(x')))`

In the programs below CFL number is taken as 0.5, and the value of dt is obtained using the formulae `dt=CFLdx/(v+a)` where a=`a=sqrt(t)` here. From the formula of dt we can notice that dt value changes at each time step, this should not be allowed as it causes inaccuracy in solution in order to avoid this the lowest of all of dt obtained at each grid point is taken as the value of dt. 0.59 is chosen as initial mass flow rate as it is close to analytical value of steady state mass flow and number of grid points n=31 is chosen for convenience.

The program for solving 1D supersonic nozzle flow equations in conservative and non-conservative forms using the Macormack Method is attached below:

The main program:

clear all
close all
clc

n=input('Enter the number of grids spaces required ideally 30,40,50 or 60:')+1;
nt=1400;
x=linspace(0,3,n);

tic;
[rho1,v1,t1,m1]=nozzle_conservative_form(n,nt);
time_of_conservative_form=toc

tic;
[rho2,v2,t2,m2]=nozzle_non_conservative_form(n,nt);
time_of_non_conservative_form=toc

conservative_rho_initial=rho1(1)
non_conservative_rho_initial=rho2(1)

conservative_rho_final=rho1(n)
non_conservative_rho_final=rho2(n)

conservative_velocity_initial=v1(1)
non_conservative_velocity_initial=v2(1)

conservative_velocity_final=v1(n)
non_conservative_velocity_final=v2(n)

conservative_temperature_initial=t1(1)
non_conservative_temperature_initial=t2(1)

conservative_temperature_final=t1(n)
non_conservative_temperature_final=t2(n)

figure(3)
plot(x,m1,'b','linewidth',1.5)
hold on
plot(x,m2,'r','linewidth',1.5)
hold on
plot(x,0.578,'k','linewidth',1.5)
xlabel('Non dimensionalised distance along nozzle-->')
ylabel('Non dimensionalised mass flow rate-->')
legend('Conservative form','Non-conservative form','Analytical')
title(' Mass flow rate non-dimensionalised comparison between conservative and non-conservative forms')
Function to solve conservative form:

function [rho,v,t,m]=nozzle_conservative_form(n,nt)

  x=linspace(0,3,n);
  dx=x(2)-x(1);

  % Initial profiles, for conservtive form should be taken more accurately as conservtive form is more sensitive
  for i=1:n
    if(x(i)>=0&&x(i)<=0.5)
      rho(i)=1;
      t(i)=1;
    elseif(x(i)>=0.5&&x(i)<=1.5)
      rho(i)=1-0.366*(x(i)-0.5);
      t(i)=1-0.167*(x(i)-0.5);
    elseif(x(i)>=1.5&&x(i)<=3.1)
      rho(i)=0.634-0.3879*(x(i)-1.5);
      t(i)=0.833-0.3507*(x(i)-1.5);
    end
  end

  a=1+2.2*(x-1.5).^2;
  v=0.59./(rho.*a); 
  %0.59 is chosen as initial mass flow rate as it is close to analytical value of steady state mass flow
  gamma=1.4;

  %In conservtive form continuity, momentum and energy equations are taken in generic form for simplicity
  %Continuity equation in generic form is: doU1/dot'=-doF1/dox'
  %Momentum equation in generic form is: doU2/dot'=-doF2/dox'+J2
  %Continuity equation in generic form is: doU3/dot'=-doF3/dox'
  %The primitiv valriables in F1,F2 and F3 are converted into terms of U1,U2 and U3 so as to maintain purity which leads to higher stability

  p=rho.*t;
  U1=rho.*a;
  U2=rho.*a.*v;
  U3=U1.*((t./(gamma-1))+(gamma/2).*(v.^2));

  %Timestep is given by dt=CFL*dx/(v+a), here a =speed of sound which is equal to t^0.5.
  %As we can see from above timestep formulae the time step changes for each grid so it does not give accurate solution unless steady state is attained.
  %So, to avoid the above problem the least value of all of the time steps is taken as dt
  %CFL number here is taken as 0.5 but it can be anything less than 1 here 

  dta=(0.5.*(dx./((t.^0.5)+v)));
  dt=min(dta)

  for k=1:nt
    U1_old=U1;
    U2_old=U2;
    U3_old=U3;
    
    %Flux terms 
    F1=U2;
    F2=((U2.^2)./U1)+(((gamma-1)/gamma)*(U3-((gamma/2)*((U2.^2)./U1))));
    F3=((gamma*U2.*U3)./U1)-(gamma*(gamma-1)*(U2.^3))./(2*U1.^2);
    
    %Predictor method
    for i=2:n-1
      dF1dx(i)=(F1(i+1)-F1(i))/dx;
      dF2dx(i)=(F2(i+1)-F2(i))/dx;
      dF3dx(i)=(F3(i+1)-F3(i))/dx;
      dadx_p(i)=(a(i+1)-a(i))/dx;
      J2_p(i)=(1/gamma)*(rho(i)*t(i))*dadx_p(i);
      
      %Continuity equation
      dU1dt_p(i)=-dF1dx(i);
      
      %Momentum equation
      dU2dt_p(i)=-dF2dx(i)+J2_p(i);
      
      %Energy equation
      dU3dt_p(i)=-dF3dx(i);
      
      %Update
      U1(i)=U1(i)+dU1dt_p(i)*dt;
      U2(i)=U2(i)+dU2dt_p(i)*dt;
      U3(i)=U3(i)+dU3dt_p(i)*dt;
      
    end
    
    %Flux terms for corrector method
    F1=U2;
    F2=((U2.^2)./U1)+(((gamma-1)/gamma)*(U3-((gamma/2)*((U2.^2)./U1))));
    F3=((gamma*U2.*U3)./U1)-(gamma*(gamma-1)*(U2.^3))./(2*U1.^2);
    
    %Corrector method
    for i=2:n-1
      dF1dx(i)=(F1(i)-F1(i-1))/dx;
      dF2dx(i)=(F2(i)-F2(i-1))/dx;
      dF3dx(i)=(F3(i)-F3(i-1))/dx;
      dadx_c(i)=(a(i)-a(i-1))/dx;
      J2_c(i)=(1/gamma)*(rho(i)*t(i))*dadx_c(i);
      
      %Continuity equation
      dU1dt_c(i)=-dF1dx(i);
      
      %Momentum equation
      dU2dt_c(i)=-dF2dx(i)+J2_c(i);
      
      %Energy equation
      dU3dt_c(i)=-dF3dx(i);
         
    end
    
    %Averaging
    dU1dt=0.5*(dU1dt_p+dU1dt_c);
    dU2dt=0.5*(dU2dt_p+dU2dt_c);
    dU3dt=0.5*(dU3dt_p+dU3dt_c);
    
    for j=2:n-1
      U1(j)=U1_old(j)+dU1dt(j)*dt;
      U2(j)=U2_old(j)+dU2dt(j)*dt;
      U3(j)=U3_old(j)+dU3dt(j)*dt;
    end
    
    %Inlet boundary conditions
    U1(1)=rho(1)*a(1);
    %U1 is a constant here 
    U2(1)=2*U2(2)-U2(3);
    v(1)=U2(1)/U1(1);
    U3(1)=U1(1)*((t(1)/(gamma-1))+((gamma/2)*v(1)^2));
    
    %Outlet boundary conditions
    U1(n)=2*U1(n-1)-U1(n-2);
    U2(n)=2*U2(n-1)-U2(n-2);
    U3(n)=2*U3(n-1)-U3(n-2);
    
    rho=U1./a;
    v=U2./U1;
    t=(gamma-1)*((U3./U1)-((gamma*v.^2)/2));
    m=rho.*v.*a;
    mach=v./((t).^0.5);
    p=rho.*t;
    
    if (rem(k,200)==0)
      figure(1)
      plot(x,m,'linewidth',1.5)
      xlabel('Non dimensionalised distance along nozzle-->')
      ylabel('Non dimensionalised mass flow rate-->')
      legend('200th time step','400th time step','600th time step','800th time step','1000th time step','1200th time step','1400th time step')
      title('Mass flow rate non-dimensionalised at different time steps comparison for conservative form')
      hold on
    end
    
  end
  rho
end
Function to solve non-conservative form:

function [rho,v,t,m]=nozzle_non_conservative_form(n,nt)

  %input
  x=linspace(0,3,n);
  dx=x(2)-x(1);

  %calculate initial profiles
  rho=1-0.3146*x;
  t=1-0.2314*x;
  v=(0.1+1.09*x).*t.^0.5;
  a=1+2.2*(x-1.5).^2;
  gamma=1.4;

  %Time steps
  dta=(0.5.*(dx./((t.^0.5)+v)));
  dt=min(dta)

  %Outer time loop
  for k=1:nt
    
    rho_old=rho;
    v_old=v;
    t_old=t;
    
    %Predictor method
      
    for j=2:n-1
      
      dvdx=(v(j+1)-v(j))/dx;
      drhodx=(rho(j+1)-rho(j))/dx;
      dtdx=(t(j+1)-t(j))/dx;
      dlogadx=(log(a(j+1))-log(a(j)))/dx;
      
      %Continuity equation
      drhodt_p(j)=-rho(j)*dvdx-rho(j)*v(j)*dlogadx-v(j)*drhodx;
      
      %Momentum equation
      dvdt_p(j)=-v(j)*dvdx-(1/gamma)*(dtdx+(t(j)/rho(j))*drhodx);
      
      %Energy equation
      dtdt_p(j)=-v(j)*dtdx-(gamma-1)*t(j)*(dvdx+v(j)*dlogadx);
      
      %Solution update
      rho(j)=rho(j)+drhodt_p(j)*dt;
      v(j)=v(j)+dvdt_p(j)*dt;
      t(j)=t(j)+dtdt_p(j)*dt;
    
    end

    %Corrector method
      
    for j=2:n-1
      
      dvdx=(v(j)-v(j-1))/dx;
      drhodx=(rho(j)-rho(j-1))/dx;
      dtdx=(t(j)-t(j-1))/dx;
      dlogadx=(log(a(j))-log(a(j-1)))/dx;
      
      %Continuity equation
      drhodt_c(j)=-rho(j)*dvdx-rho(j)*v(j)*dlogadx-v(j)*drhodx;
      
      %Momentum equation
      dvdt_c(j)=-v(j)*dvdx-(1/gamma)*(dtdx+(t(j)/rho(j))*drhodx);
      
      %Energy equation
      dtdt_c(j)=-v(j)*dtdx-(gamma-1)*t(j)*(dvdx+v(j)*dlogadx);
    
    end

    drhodt=0.5*(drhodt_p+drhodt_c);
    dvdt=0.5*(dvdt_p+dvdt_c);
    dtdt=0.5*(dtdt_p+dtdt_c);
    
    %Final solution update
    for i=2:n-1
      rho(i)=rho_old(i)+drhodt(i)*dt;
      v(i)=v_old(i)+dvdt(i)*dt;
      t(i)=t_old(i)+dtdt(i)*dt;
    end
    
    %Applying boundary conditions
    
    %Inlet
    v(1)=2*v(2)-v(3);
    
    %Outlet
    v(n)=2*v(n-1)-v(n-2);
    rho(n)=2*rho(n-1)-rho(n-2);
    t(n)=2*t(n-1)-t(n-2);
    
    m=rho.*v.*a;
    mach=v./((t).^0.5);
    p=rho.*t;
    
    if (rem(k,200)==0)
      figure(2)
      plot(x,m,'linewidth',1.5)
      xlabel('Non dimensionalised distance along nozzle-->')
      ylabel('Non dimensionalised mass flow rate-->')
      legend('200th time step','400th time step','600th time step','800th time step','1000th time step','1200th time step','1400th time step')
      title('Mass flow rate non-dimensionalised at different time steps comparison for non-conservative form')
      hold on
    end
    
  end
    rho
end
The mass flow rates obtained at different time steps for conservative form is attached below:



From the above figure we can see that the from 800th step the solution is stabilizing, so 800 steps are required of the obtained dt which is 0.020363 for n=31 for the solution to converge in conservative form.

The mass flow rates obtained at different time steps for non-conservative form is given below:



Similarly from the above figure we can see that the solution is stabilizing at around 700th time step, so 700 time steps of the obtained dt value which is 0.020134 for n=31 is required for the solution to converge in non-conservative form. 

The attached figure below compares the the mass flow rates obtained in both conservative, non-conservative and analytical solution. This figure can be used for validation of the data obtained with John D. Anderson book chapter 7.



The below table gives the comparison of different values at number of grid points 31 and 41 for grid dependence test:

 

Conservative form n=31

Non-conservative form n=31

Conservative form n=41

Non-conservative form n=41

Density(Initial)

1

1

1

1

Density(Final)

0.0527

0.0528

0.0526

0.0527

Velocity(Initial)

0.0983

0.0987

0.0981

0.0983

Velocity(Final)

1.8662

1.8618

1.8648

1.8623

Temperature(Initial)

1

1

1

1

Temperature(Final)

0.30524

0.30828

0.3063

0.3080

From the above we can easily see that the increasing the number of grid points from 31 to 41 does not impact the values obtained, thus the solution obtained at n=31 is grid independent.

A table is attached below which compares the time taken by both the forms at n=31 and n=41:

Time taken by each form

Conservative form

Non-conservative form

n=31

17.511

17.898

n=41

19.360

22.133

From the above table we can easily see conservative form is faster than non-conservative form.
