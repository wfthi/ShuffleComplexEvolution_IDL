; Usage
;
; The main routine is sce_ua
;
; 22/3/2018 add the Shekel test function
;
; IDL > .r ShuffleComplexEvolution.pro
;
;       ftol=1e-3
;       sce_ua,ftol,'rosenbrock',[-5.,-2.],[5.,8.],10000,seed,parameters,merit,nevolution_steps=20,ncomplexes=2,extinction=1,max_nevolution_steps=20
;
;       where rosenbrock is the name of an IDL function that computes
;       the cost function to be minimized.
;       10000 is the maximum number of function calls
;
; test functions:
; IDL> test_sce
;
; ----------------------------------------------------------------------------
; https://en.wikipedia.org/wiki/Test_functions_for_optimization
pro test_sce,parameters,merit

ftol = 1e-5
print
print,' ------------------------------------'
print,'Test: Shekel m=10'
sce_ua,ftol,'shekel10',[0.,0.,0.,0.],[10.,10.,10.,10.],10000,seed,parameters,merit,ncomplexes=4,nevolution_steps=20,extinction=1
print,'Tolerance =',ftol
print,'global minimum at (4,4,4,4) = -10.5364'
print,'With extinction'

ftol=1e-4
print
print,' ------------------------------------'
print,'Test: Himmelblau'
sce_ua,ftol,'himmelblau',[-5.,-5.],[5.,5.],10000,seed,parameters,merit,ncomplexes=4,nevolution_steps=20
print,'Tolerance =',ftol
print,'Global minimum at (3.,2.)=0.'
print,'+3 others mimina'

ftol=1e-3
print
print,' ------------------------------------'
print,'Test: Rosenbrock'
sce_ua,ftol,'rosenbrock',[-5.,-2.],[5.,8.],10000,seed,parameters,merit,nevolution_steps=20,ncomplexes=2
sce_ua,ftol,'rosenbrock',[-5.,-2.],[5.,8.],10000,seed,parameters,merit,nevolution_steps=20,ncomplexes=2,extinction=1,max_nevolution_steps=20
print,'Tolerance =',ftol
print,'Global minimum at (1.,1.)=0.'

ftol=1e-3
print
print,' ------------------------------------'
print,'Goldstein Price'
print,'Bound X1=[-2,2], X2=[-2,2]'
print,'Global Optimum: 3.0,(0.0,-1.0)'
sce_ua,ftol,'goldstein_price',[-2.,-2.],[2.,2.],1000,seed,parameters,merit,ncomplexes=2,nevolution_steps=20
sce_ua,ftol,'goldstein_price',[-2.,-2.],[2.,2.],1000,seed,parameters,merit,ncomplexes=2,nevolution_steps=20,extinction=1
;stop

ftol=1e-4
print
print,' ------------------------------------'
print,'The Rastrigin Function has many hills and valleys'
print,'Bound: X1=[-1,1], X2=[-1,1]'
print,'Global Optimum: -2, (0,0)'
sce_ua,ftol,'rastrigin',[-1.,-1.],[1.,1.],10000,seed,parameters,merit,ncomplexes=2,min_ncomplexes=2
sce_ua,ftol,'rastrigin',[-1.,-1.],[1.,1.],10000,seed,parameters,merit,ncomplexes=2,extinction=1,min_ncomplexes=2
;stop
;
ftol=1e-4
print
print,' ------------------------------------'
print,'Test: 6 hump camelback'
print,'Bound: X1=[-3,3], X2=[-2,2]'
print,'True Optima: -1.031628453489877, (-0.08983,0.7126), (0.08983,-0.7126)'

sce_ua,ftol,'camelback',[-3.,-2.],[3.,2.],10000,seed,parameters,merit,nevolution_steps=10,ncomplexes=4
print
print,'Test: 6 hump camelback'
print,'Extinction option on'
sce_ua,ftol,'camelback',[-3.,-2.],[3.,2.],10000,seed,parameters,merit,nevolution_steps=10,ncomplexes=4,extinction=1,min_ncomplexes=1
;stop
;
ftol=1e-4
print
print,' ------------------------------------'
print,'Test: Griewank 2d'
print,'minimum at origin =0'
sce_ua,ftol,'griewank_2d',[-600.,-600.],[600.,600.],10000,seed,parameters,merit,nevolution_steps=20

ftol=1e-5
print
print,' ------------------------------------'
print,'Test: Griewank 10d'
print,'minimum at origin =0'
bl=replicate(-400.0d0,10)
bu=replicate(400.0d0,10)
sce_ua,ftol,'griewank_10d',bl,bu,10000,seed,parameters,merit,nevolution_steps=20,ncomplexes=2,extinction=1


ftol=1e-3
print
print,' ------------------------------------'
print,'Test: noisy quartic 30d'
print,'minimum =< 30'
n = 30
bl=replicate(-1.28d0,n)
bu=replicate(1.28d0,n)
sce_ua,ftol,'noisy_quartic_30d',bl,bu,15000,seed,parameters,merit,nevolution_steps=20,ncomplexes=6


ftol=1e-4
print
print,' ------------------------------------'
print,'Test: sphere'
print,'minimum at (0.,0.,0.)=0.'
bl=[-5.12d0,-5.12d0,-5.12d0]
bu=[5.12d0,5.12d0,5.12d0]
sce_ua,ftol,'sphere',bl,bu,10000,seed,parameters,merit,nevolution_steps=5,ncomplexes=2,extinction=1


ftol=1e-2
print
print,' ------------------------------------'
print,'Test: step'
print,'minimum at (-5.,-5.,-5.,-5.,-5.)=0.'
bl=[-5.12d0,-5.12d0,-5.12d0,-5.12d0,-5.12d0]
bu=[5.12d0,5.12d0,5.12d0,5.12d0,5.12d0]
sce_ua,ftol,'step',bl,bu,10000,seed,parameters,merit,nevolution_steps=5,ncomplexes=6,nelements_complex=5,extinction=4,max_nevolution_steps=50
    
ftol=1e-4
print
print,' ------------------------------------'
print,'Test: Ackley function'
; search domain -15=< x_i =< 30
print,'global minimum at (0.,0.,)=0.'
; several local minima
bl=[-15.d0,-15.d0]
bu=[30.d0,30.d0]
sce_ua,ftol,'ackley',bl,bu,10000,seed,parameters,merit,nevolution_steps=5,ncomplexes=2,nelements_complex=5,extinction=1,max_nevolution_steps=20

ftol=1e-4
print
print,' ------------------------------------'
print,'Test: Beale function'
print,'Bounds -4.5 =< x_i =< 4.5'
print,'Global minimum at (3.,0.5)=0.'
bl=[-4.5d0,-4.5d0]
bu=[4.5d0,4.5d0]
sce_ua,ftol,'beale',bl,bu,10000,seed,parameters,merit,nevolution_steps=5,ncomplexes=2,nelements_complex=5,extinction=1,max_nevolution_steps=20

ftol=1e-4
print
print,' ------------------------------------'
print,'Test: Booth function'
print,'bounds -10 =< x_i =< 10'
print,'several local minima'
print,' Global minimum at (1.,3.)=0.'
bl=[-10d0,-10d0]
bu=[10d0,10d0]
sce_ua,ftol,'booth',bl,bu,10000,seed,parameters,merit,nevolution_steps=5,ncomplexes=2,nelements_complex=5,extinction=1,max_nevolution_steps=20

ftol=1e-5
print
print,' ------------------------------------'
print,'Test: Schwefel function'
print,'bounds -500=< x_i =< 500 i=1,2'
print,'several local minima'
print,'Global minima at (420.9687d0,420.9687d0)=0.'
bl=[-500d0,-500d0]
bu=[500d0,500d0]
sce_ua,ftol,'schwefel',bl,bu,10000,seed,parameters,merit,nevolution_steps=10,ncomplexes=4,nelements_complex=5,max_nevolution_steps=30

ftol=1e-4
print
print,' ------------------------------------'
print,'Test: Goldstein Price function'
print,'bounds -1 =< x_i =< 2 '
print,'several local minima'
print,'Global minimum at   (0.,-1.)=3.0'
bl=replicate(-1.0d0,2)
bu=replicate(2.0d0,2)
sce_ua,ftol,'goldstein',bl,bu,10000,seed,parameters,merit,nevolution_steps=5,ncomplexes=2,nelements_complex=5,extinction=1,max_nevolution_steps=30

ftol=1e-4
print
print,' ------------------------------------'
print,'Test: Shubert function'
print,'bounds -10. =< x_i =< 10.'
print,'several local minima'
print,'Global minimum -186.7309'
;
bl=replicate(-10.0d0,2)
bu=replicate(10.0d0,2)
sce_ua,ftol,'shubert',bl,bu,10000,seed,parameters,merit,nevolution_steps=5,ncomplexes=2,nelements_complex=5,extinction=2,max_nevolution_steps=30

ftol=1e-4
print
print,' ------------------------------------'
print,'Easom function'
print,'bounds -10. < x_i < 10.'
print,'several local minima'
print,'Global minimum (pi,pi)=-1.'
;
bl=replicate(-10.0d0,2)
bu=replicate(10.0d0,2)
sce_ua,ftol,'easom',bl,bu,10000,seed,parameters,merit,nevolution_steps=5,ncomplexes=2,nelements_complex=5,max_nevolution_steps=30

ftol=1e-4
print
print,' ------------------------------------'
print,'f2_test function, n=10'
print,'bounds 0. < x_i < 1.'
print,'several local minima'
print,'Global minimum (0.4,0.4,...,0.4)=0.'
;
bl=replicate(0.0d0,10)
bu=replicate(1.0d0,10)
sce_ua,ftol,'f2_test',bl,bu,15000,seed,parameters,merit,ncomplexes=2,extinction=2

ftol=1e-4
print
print,' ------------------------------------'
print,'Test: Zakharov function, n=2'
print,'bounds -5. < x_i < 10.'
print,'no local minimum'
print,'Global minimum (0.,0.,...,0.)=0.'
;
bl=replicate(-5.0d0,2)
bu=replicate(1.0d1,2)
sce_ua,ftol,'zakharov',bl,bu,10000,seed,parameters,merit,nevolution_steps=5,ncomplexes=2,nelements_complex=5,max_nevolution_steps=20

ftol=1e-4
print
print,' ------------------------------------'
print,'Test: Branin function, n=2'
print,'bounds -5. < x_1 < 10., 0.<x_2<15  '
print,'no local minimum'
print,'Global minimum at (-pi,12.275), (pi,2.275), (9.42478,2.475)'
;
bl=[-5.0d0,0.0d0]
bu=[10.0d0,15.0d0]
sce_ua,ftol,'branin',bl,bu,10000,seed,parameters,merit,nevolution_steps=5,ncomplexes=2,nelements_complex=5,max_nevolution_steps=20

ftol=1e-4
print
print,' ------------------------------------'
print,'Test: f1_test function, n=2'
print,'bounds -1. < x_1 < 1.'
print,'manu local minima'
print,'Global minimum at (0.,0.)=-2.'
;
bl=[-1.0d0,-1.0d0]
bu=[1.0d0,1.0d0]
sce_ua,ftol,'f1_test',bl,bu,10000,seed,parameters,merit,nevolution_steps=5,ncomplexes=2,nelements_complex=5,max_nevolution_steps=20

ftol=1e-4
print
print,' ------------------------------------'
print,'Test: Zimmermann''s problem, n=2'
print,'bounds x_i >0. i=1,2'
print,'local minima?' 
print,'Global minimum at (7.,2.)=0.'
;
bl=[0.0d0,0.0d0] 
bu=[100.0d0,100.0d0] 
sce_ua,ftol,'zimmermann',bl,bu,10000,seed,parameters,merit,nevolution_steps=5,ncomplexes=10,nelements_complex=5,max_nevolution_steps=20,extinction=1

ftol=1e-4
print
print,' ------------------------------------'
print,'Test: Corana''s parabola, n=4'
print,'bounds -1000. < x_i < 1000. i=0,1,2,3'
print,'many local minima'
print,'Global minimum at abs(x_j)<0.05'
;
bl=replicate(-1d3,4)
bu=replicate(1d3,4) 
sce_ua,ftol,'corana',bl,bu,10000,seed,parameters,merit,nevolution_steps=5,ncomplexes=2,nelements_complex=5,max_nevolution_steps=20

end

; -------------------------------------------------------------------------------------
; Test Amoeba a Simplex locally-convergent optimizer

pro test_amoeba

ftol=1e-4
initial=10.*randomu(seed,2)-5.0
print,'initial values =',initial
print
print,'The standard Nelder-Meade Simplex method (called here amoeba) can fail.'
print,'It suceeds for the Rosenbrock (without local minima) function.'
print
result=amoeba(ftol,function_name='himmelblau',p0=initial,scale=[5.,5.],ncalls=ncalls)
print,'Amoeba Himmelblau function'
print,'Global minimum at (3.,2.)=0.'
print,'parameters =',result
print,'function =',call_function('himmelblau',result)
print,'number of function calls=',ncalls
print

result=amoeba(ftol,function_name='rosenbrock',p0=initial,scale=[5.,8.],ncalls=ncalls)
print,'Amoeba Rosenbrock function'
print,'Global minimum at (1.,1.)=0.'
print,'parameters =',result
print,'function =',call_function('rosenbrock',result)
print,'number of function calls=',ncalls
print

ftol=1e-6
result=amoeba(ftol,function_name='rastrigin',p0=initial,scale=[5.,5.],ncalls=ncalls)
print,'Amoeba Rastrigin function'
print,'Global Optimum: -2, (0,0)'
print,'parameters =',result
print,'function =',call_function('rastrigin',result)
print,'number of function calls=',ncalls

end

; ----------------------------------------------------------------------------
; Test functions

function himmelblau,x

; This is the Himmelblau Function (see textbook by himmelblau)
; Bound: X1=[-5,5], X2=[-5,5]
; Global Optimum: 0 (3,2)
; https://en.wikipedia.org/wiki/Himmelblau%27s_function  

f1 = (x(0)*x(0)+x(1)-11.0d0)
f2 = (x(0)+x(1)*x(1)-7.0d0)

return,f1*f1+f2*f2

end

function rosenbrock,x

; The Rosenbrock Function is unimodal
; Bound: X1=[-5,5], X2=[-2,8] 
; Global Optimum: 0,at (1,1) is not easy to find because it is situated in a
; valley with a flat bottom.
; https://en.wikipedia.org/wiki/Rosenbrock_function  

x1=x(0)
x2=x(1)
a = 100.0d0
f=a*(x2-x1*x1)^2+(1.0-x1)^2

return,f
 
end

function goldstein_price,x

; This is the Goldstein-Price Function
; Bound X1=[-2,2], X2=[-2,2]
; Global Optimum: 3.0,(0.0,-1.0)
;      
x1 = x(0)
x2 = x(1)
u1 = (x1 + x2 + 1.0d0)^2
u2 = 19.0d0 - 14.0d0*x1 + 3.0d0*x1^2 - 14.0d0*x2 + 6.0d0*x1*x2 +3.0d0*x2^2
u3 = (2.0d0*x1 - 3.0d0*x2)^2
u4 = 18.0d0 - 32.0d0*x1 + 12.0d0*x1^2 + 48.0d0*x2 -36.0d0*x1*x2 + 27.0d0*x2^2
u5 = u1 * u2
u6 = u3 * u4
f = (1.0d0 + u5) * (30.0d0 + u6)

return, f

end

function camelback,x

;  This is the Six-hump Camelback Function.
;  Bound: X1=[-3,3], X2=[-2,2]
;  True Optima: -1.031628453489877, (-0.08983,0.7126), (0.08983,-0.7126)
;  also called sixmin 

x1 = x(0)
x2 = x(1)
f = (4.0d0-2.1d0*x1^2+(x1^4.)/3.)*x1^2+x1*x2+(-4.0d0+4.0d0*x2^2)*x2^2
return,f

end

function rastrigin,x

;
;  The Rastrigin Function has many hills and valleys
;  Bound: X1=[-1,1], X2=[-1,1]
;  Global Optimum: -2, (0,0)
 
x1 = x(0)
x2 = x(1)
f = x1^2 + x2^2 - cos(18.0d0*x1) - cos(18.0d0*x2)

return,f

end

function two_peak_trap,x

; The global maximum of this function is x=0

if x ge 0.0d0 and x lt 15.0d0 then f=160.0d0*(15.0d0-x)/15.0d0
if x ge 15.0d0 and x le 20.0d0 then f=200.0d0*(x-15.0d0)/5.0d0

return,f

end

function griewank_2d,x

;  This is the Griewank Function (2-D or 10-D)
;  Bound: X(i)=[-600,600], for i=1,2,...,10
;  Global Optimum: 0, at origin

d = 200.0

u1 = 0.0
u2 = 1.0
for j =0,1 do begin
    u1 = u1 + x(j)*x(j)/d
    u2 = u2 * cos(x(j)/sqrt(double(j+1)))
endfor
f = u1 - u2 + 1.0d0

return,f

end

function griewank_10d,x

; Griwwabgk's function 10 dimensions
; bounds -400. < x_j < 400. j=0,...,9
; Global minimum at origin = 0.

n=10
d = 4000.0

u1 = 0.0
u2 = 1.0
for j =0,n-1 do begin
    u1 = u1 + x(j)*x(j)/d
    u2 = u2 * cos(x(j)/sqrt(double(j+1)))
endfor
f = u1 - u2 + 1.0d0

return,f

end


function step,x
; Storn & Price 1997
; 5 numbers

f=25.d0+total(fix(x))

return,f
end


function sphere,x

; bound xj = [-5.12, 5.12]
; minimum at (0.,0.,0.) fmin=0

f = x(0)*x(0)+x(1)*x(1)+x(2)*x(2)

return,f 
end

function ackley,x

; Ackley function.
; search domain -15=< x_i =< 30
; global minimum at (0,0,...)=0.
; several local minima
; The number of variables n should be adjusted below.
; The default value of n =2.
; 
n = 2;
a = 20.0d0; 
b = 0.2d0; 
c = 2.0d0*!pi;
s1 = 0.0d0; 
s2 = 0.0d0;
for i=0,n-1 do begin
   s1 = s1+x(i)^2;
   s2 = s2+cos(c*x(i));
endfor
f = -a*exp(-b*sqrt(1.0d0/double(n)*s1))-exp(1.0d0/double(n)*s2)+a+exp(1.0d0);
return,f
end

function beale,x
; 
; Beale function.
; The number of variables n = 2.
; Bounds -4.5 =< x_i =< 4.5
; Global minimum at (3.,0.5)=0.
f = (1.5d0-x(0)*(1.0d0-x(1)))^2+(2.25d0-x(0)*(1.0d0-x(1)^2))^2+(2.625d0-x(0)*(1-x(1)^3))^2
return,f
end

function booth,x
; 
; Booth function 
; bounds -10 =< x_i =< 10
; several local minima
; Global minimum at (1.,3.)=0.
; The number of variables n = 2.
;
f  = (x(0)+2.0d0*x(1)-7.0d0)^2+(2.0d0*x(0)+x(1)-5.0d0)^2
return,f
end

function schwefel,x

; Schwefel function
; bounds -500=< x_i =< 500 i=1,2,3,...,n
; several local minima
; Global minima at 420.9687d0*(1.,1.,...,1.)=0.
; The number of variables n should be adjusted below.
; The default value of n = 2.
; 
n = 2;
s = total(-x*sin(sqrt(abs(x))))
f = 418.9829d0*double(n)+s
;y=[420.9687d0,420.9687d0]
;print,418.9829d0*double(n)+total(-y*sin(sqrt(abs(y))))
return,f
end

function michalewicz,x

; Michalewicz function 
; bounds 0 =< x_i =< pi 
; several local minima
; Global minimum at 
;  n=2 -1.8013
;  n=5  -4.4687658
;  n=10 -9.66015
; The number of variables n should be adjusted below.
; The default value of n =2.
; 
;y=[2.20d0,1.570d0,1.29d0,1.92d0,1.72d0,1.57d0,1.45d0,1.76d0,1.66d0,1.57d0]

n = 10
m = 10.0d0
s = 0.0d0 
for i = 0,n-1 do begin
    s = s+sin(x(i))*(sin((double(i)*x(i)^2.0d0)/!dpi))^(2.0d0*m);
endfor
f = -s
return,f
end

function goldstein,x
 
; Goldstein and Price function 
; bounds -2 =< x_i =< 2 i=1,2
; several local minima
; Global minimum (0.,1.)=3.0d0
; The number of variables n = 2.
; 
a=1.0d0+(x(0)+x(1)+1.0d0)^2.0d0*$
        (19.0d0-14.0d0*x(0)$
         +3.0d0*x(0)^2.0d0-14.0d0*x(1)$
         +6.0d0*x(0)*x(1)+3.0d0*x(1)^2.0d0)

b=30.0d0+(2.0d0*x(0)-3.0d0*x(1))^2.0d0*$
        (18.0d0-32.0d0*x(0)$
         +12.0d0*x(0)^2.0d0$
         +48.0d0*x(1)$
         -36.0d0*x(0)*x(1)$
         +27.0d0*x(1)^2.0d0);

f = a*b

return,f
end

function easom,x

; Easom function
; n = 2
; bounds -10< x_i < 10 i=1,2
; several local minima
; Global minimum at (pi,pi)=-1.

f = -cos(x(0))*cos(x(1))*exp(-(x(0)-!dpi)^2-(x(1)-!dpi)^2)

return,f
end


function branin,x
;
; Branin RCOS function
; n = 2
; bounds -5<x_1<10, 0<x_2<15
; no local minimum
; Global minima at (-pi,12.275), (pi,2.275), (9.42478,2.475)
f = (x(1)-(5.1d0/(4.0d0*!dpi^2))*x(0)*x(0)+5.0d0*x(0)/!dpi-6.0d0)^2 $
    + 10.0d0*(1.0d0-1.0d0/(8.0d0*!dpi))*cos(x(0))+10.0d0
return,f
end

function f1_test,x
;
; n = 2
; bounds -1< x_i < 1 i=1,2
; many local minima
; Global minimum at (0.,0.)=-2
f = x(0)^2+x(1)^2-cos(18.0d0*x(0))-cos(18.0d0*x(1))

return,f

end


function f2_test,x

; n=10
; many local minima
; bounds 0< x_i < 1 i=1,2,...,10
; Global minimum at (0.4,0.4,...,0.4)=0.
;
a = 0.05d0
n=10
f=0.0d0
for i=0,n-1 do begin
f = f+min([abs(x(i)-0.2d0)+a,abs(x(i)-0.4d0),abs(x(i)-0.7d0)+a])
endfor
return,f
end


function shubert,x
;
; Shubert function
; bounds -10. =< x_i =< 10.
; several local minima
; Global minimum -186.7309
; The number of variables n =2.
; 
s1 = 0.0
s2 = 0.0
for i = 1,5 do begin 
    i1=double(i)
    i2=double(i+1)	
    s1 = s1 + i1*cos(i2*x(0)+i1)
    s2 = s2 + i1*cos(i2*x(1)+i1)
endfor
f = s1*s2
return,f
end

function zakharov,x

; Zakharov function
; n = 2,5,10
; bounds -5< x_i < 10 i=1,2,..., n
; no local minima
; Global minima at (0.,0.,...,0.)=0.
n=2
j=dindgen(n)+1
y = (total(0.50d*j*x))^2
f = total(x^2)+y*(1.0d0+y)

return,f
end

function zimmermann,x

; Zimmerman's problem
; n=2
; bounds x_i > 0 i=1,2

f= 9.0d0-x(0)-x(1)
p = (x(0)-3.0d0)^2+(x(1)-2.0d0)^2
q = x(0)*x(1)
if (p gt 16.0d0 or q gt 14.0d0) then f=100.0d0

return,f

end

function corana,x

; Corana's parabola
; n = 2
; bounds -1000. < < 1000.
; minimum at -0.05 < x_j < 0.05 j=0,1,2,3
; minimum = 0.

z = fix(abs(x/0.2d0)+0.49999d0)*(fix(x gt 0.) - fix(x lt 0.))*0.2d0
d = [1.0d0,1d3,1d1,1d2]
f = 0.0d0
for i=0,3 do begin
   if abs(x(i)-z(i)) lt 0.05d0 then begin
        signz=fix(z(i) gt 0.) - fix(z(i) lt 0.)
        f = f + (0.15d0*(z(i)-0.05d0*signz)^2*d(i))
    endif else begin
        f = f + d(i)*x(i)*x(i)
    endelse
endfor

return,f

end

function noisy_quartic_30d,x

; de Jong noisy quartic
; bounds -1.28 =< x_j -< 1.28

n=30
rn=randomu(seed)
f=0.0d0
for i=0,n-1 do begin
    f = f+ (double(i+1)*x(i)^4)*(1.d0+2.0d0*rn)
;    f = f+ (double(i+1)*x(i)^4)
endfor
;f=f+rn
return,f

end

function shekel10,x
; adapated from https://www.sfu.ca/~ssurjano Matlab version
; m local minima
; bounds 0. < (x1,x2,x3,x4) < 10. x=(x1,x2,x3,x4)
; global minimum at (4,4,4,4) f = -10.5364
  
m= 10
b = [1.,2.,2.,4.,4.,6.,3.,7.,5.,5.]*0.1

C = [[4.,1.,8.,6.,3.,2.,5.,8.,6.,7.],[4.,1.,8.,6.,7.,9.,3.,1.,2.,3.6],[4.,1.,8.,6.,3.,2.,5.,8.,6.,7.],[4.,1.,8.,6.,7.,9.,3.,1.,2.,3.6]]
outer = 0.
for ii = 0,m-1 do begin
   bi=b(ii)
   inner = 0.
   for jj = 0,3 do begin
      xj = x(jj)
      Cij = C(ii,jj)
      inner = inner + (xj-Cij)^2
   endfor
   outer = outer + 1./(inner+bi)
endfor
y = - outer
return,y

end
; ----------------------------------------------------------------------------
pro sce_simplex,ftol,func,bl,bu,max_func_calls,archival_parameters,archival_merit,seed,ncalls=ncalls,verbose=verbose

; NAME:
;	sce_simplex
;
; PURPOSE:
;	Multidimensional minimization of a function func(X), where
;	X is an N-dimensional vector, using the downhill sce_simplex
;	method of Nelder and Mead, 1965, Computer Journal, Vol 7, pp 308-313.
;
; CATEGORY:
;	Function minimization/maximization. Sce_Simplex method.
;
; CALLING SEQUENCE:
;	sce_simplex,Ftol,'function',...
;
; INPUT:
;    FTOL:  the fractional tolerance to be achieved in the function
;	value.  e.g. the fractional decrease in the function value in the
;	terminating step.  This should never be less than the
;	machine's single or double precision.
;
;    func: name of the function to be minimized/maximized
;
; OUTPUT:
;    archival_parameters
;    archival_merit
;
; OPTIONAL INPUTS:
;   ncalls:
;
; AUTHOR:
;      Wing-Fai Thi (May 2007)
;
; Reference:
;      Nelder & Mead, 1965, Computer Journal, Vol 7, pp 308-313.
;
; ----- Initialize parameters
  
nbu=n_elements(bu)
nbl=n_elements(bl)
if not(nbu eq nbl) then begin
    print,'Error in input boundaries'
    stop
endif
nparameters=nbu
bound = bu - bl

; ----- 

nvertices           = nparameters+1
vertex              = dblarr(nvertices,nparameters)
archival_parameters = dblarr(nvertices,nparameters)
merit               = dblarr(nvertices)
archival_merit      = dblarr(nvertices)
gnrng               = dblarr(nparameters)

; --- Generate initial vertices

ncalls  = 0

for i=0,nvertices-1 do begin
    vertex(i,*)=bl+randomu(seed,nparameters)*bound
    archival_parameters(i,*) = vertex(i,*)
    merit(i) = call_function(func,vertex(i,*))
    ncalls = ncalls + 1
end 
archival_merit = merit

; ----- Order the vertices

order   = sort(merit)
vertex  = vertex(order,*)  
merit   = merit(order)

; ----- Record the best and worst points

best_vertex  = vertex(0,*)
best_merit   = merit(0)
worst_vertex = vertex(nvertices-1,*) 
worst_merit  =  merit(nvertices-1)

; ----- Initial convergence criterium
    
success = 0 
nloop   = 0

; ----- Start search loop

while ncalls lt max_func_calls and  success eq 0 do begin

nloop = nloop + 1

; ----- Generate a new sce_simplex point

generate_offspring,func,vertex,merit,bl,bu, $
                   bl,bu,ncalls,new_vertex,new_merit, $
                   archival_parameters,archival_merit

; ----- Replace the worst point in Sce_Simplex with the new point:

vertex(nvertices-1,*) = new_vertex 
merit(nvertices-1)    = new_merit   

; ----- Order the sce_simplexes

order   = sort(merit)
vertex  = vertex(order,*)  
merit   = merit(order)

; ----- Record the best and worst points

best_vertex  = vertex(0,*)
best_merit   = merit(0)
worst_vertex = vertex(nvertices-1,*) 
worst_merit  =  merit(nvertices-1)

; ----- Computes the normalized geometric range of the parameters

for i=0,nparameters-1 do gnrng(i)=exp(mean(alog((max(vertex(*,i))-min(vertex(*,i)))/bound(i))))  

; ----- Test for success
; relative difference in cost function value between the best and
; the worst vertex is below the criterion ftol

d    = abs(worst_merit)+abs(best_merit)
crit = 2.0 * abs(worst_merit-best_merit)/d
;crit = min(gnrng)
; F-test method
; ftol = nparameters/(ndata-nparameters)*f_cvf(confidence,nparameters,ndata-nparameters)
;crit = abs((best_merit-worst_merit)/best_merit)

if crit lt ftol then success = 1

if keyword_set(verbose) then begin
    print,'Number of function calls =',ncalls,' best =',merit(0),' crit=',crit
endif
endwhile
; --------------------------------

; Display the results

if success eq 1 then begin
   print 
   print,'Convergence criterion :',crit
   print
   print,'number of function calls:',ncalls
   print
   print,'Results'
   print,'Parameters:'
   print,best_vertex
   print,'Evaluation function=',best_merit 
endif

if ncalls ge max_func_calls then begin
    print,'Optimization search terminated because the limit'
    print,'on the maximum number of trials'
    print, max_func_calls
    print,'has beem exceeded. Search was stopped at trial number:'
    print,ncalls
    print,'of the initial loop!'
;    stop 
endif
; -----
jsort = sort(archival_merit)
narchival = indgen(n_elements(jsort))
archival_parameters(narchival,*) = archival_parameters(jsort,*)
archival_merit                   = archival_merit(jsort)

;stop
end

; ----------------------------------------------------------------------------
pro sce_ua,ftol,func,bl,bu,max_func_calls,seed,archival_parameters $
          ,archival_merit $
          ,ncomplexes=ncomplexes $
          ,nelements_complex=nelements_complex $
          ,nevolution_steps=nevolution_steps $
          ,ncalls=ncalls $
          ,verbose=verbose $
          ,max_nshuffles=max_nshuffles $
          ,extinction=extinction $
          ,min_ncomplexes=min_ncomplexes $
          ,max_nevolution_steps=max_nevolution_steps

; Purpose:
;   Multidimensional minimization of a function func(x), where
;   x is an n-dimensional vector, using the Shuffled Complex
;   Evolution (SCE-UA) Optimization Method of Duan et al. with
;   some modifications.
;
; Description:
;   The SCE-UA method is an heuristic global optimization method that 
;   combines features from Genetic Algorithms and Sce_Simplex Algorithms.
;   There is a high probability to find a global minimum but this
;   has not been proven mathematically. 
;
;   The SCE-UA method starts with the initial selection of a
;   "population" of points distributed randomly or pseudo randomly
;   or quasi-randomly throughout the feasibe parameter space. 
;   The population is then partitioned into several (ncomplexes) "complexes",
;   each containing at least 2n+1 points, where n is the number of parameters 
;   to be constrained. 
;    Each complex evolves independently according to a
;   "reproduction" process that, in turn, uses the Simpex Method
;   (Nelder & Mead, 1965) but without the shrinking step, which is
;   replaced by a randomly generated point.  
;   At periodic stages (chosen by the user: nevolution_steps), 
;   the entire population is shuffled and points are reassigned to new complexes
;   formed so that the information gained by the previous complexes is shared. 
;   The best "population" points in each complex are grouped into a new
;   complex. The following best points are grouped into a second complex and
;   so forth. The idea is that the best point (the local minimum) in each complex 
;   will "bread" with other best points of an other complexes. The "offsprings" will have
;   characteristics (parameters) closer to the best of all (closer to that set of
;   parameters that give rise to the lower merit function).
;
;   The evolution and shuffling steps continue until the prescribed
;   convergence criteria are reached.
;
;   The method combines the advantages of deterministic (Sce_Simplex)
;   and stochastic (Genetic) methods. While the genetic algorithm
;   allows global optimization, the method is slow. On the hand, the
;   Simplex (Sce_simplex) method is rapid but prompt to find local minima.
;
;   In principle other local optimizer can be used instead of the Nelder-Mead
;   Sinmplex algorithm such as Sequencial Quadratic Programming (SQP).
;
;   The performance of a global optimization solver depends on two
;   characteristics, the effectiveness and the efficiency (Duan et
;   al. 1992)
;
;   Further improvements over the original algorithm includes the
;   inclusion of the concept of "extinction".
;   The principle is the decrease of number of complexes after a
;   certain number of generations, the worst points being eliminated.
;
;   The choice of the stop criterion ftol is important and difficult
;   to make.
;
; Reference:
;
; Duan, Q., A Global Optimization Strategy for Efficient and
;      Effective Calibration of Hydrologic Models, Ph.D.
;      dissertation, University of Arizona, Tucson, Arizona, 1991
;
; Duan, Q., V.K. Gupta, and S. Sorooshian, A Shuffled Complex
;      Evolution Approach for Effective and Efficient Global
;      Minimization, Journal of Optimization Theory and Its
;      Applications, Vol 61(3), 1993
;
; Duan, Q., S. Sorooshian, and V.K. Gupta, Effective and Efficient
;      Global Optimization for Conceptual Rainfall-Runoff Models,
;      Water Resources Research, Vol 28(4), pp. 1015-1031, 1992
;
; Duan, Q., Sorooshian S., & Gupta V. K, Optimal Use of the SCE-UA
;      Method for Calibrating Watershed Models, Journal of Hydrology, vol
;      158, 265-294, 1994
;
; Nelder & Mead, 1965, Computer Journal, Vol 7, pp 308-313.
;
;  
; Other methods: Simplex-Simulated annealing (SIMPSA) Cardoso et al. 1996
;                
; used in Thi et al. 2010 Monthly Notices of the Royal Astronomical Society, Volume 406, Issue 3, pp. 1409-1424
;  
; Input:
;
;   ftol tolerance
;   func name of the evaluation (merit, crirterion, cost, ...) function
;   bl   lower bounds to the parameters
;   bu   upper bounds to the parameters
;   max_func_calls maximum number of function calls
;
; Optional input (keyword):
;
;   ncomplexes : number of complexes (default 2)
;   nelements_complex : number of elements per complex
;   nevolution_steps : number of evolutionary steps
;   ncalls : maximum number of function calls
;   verbose
;
; Output:
;
;   archival_parameters  all parameters set
;   archival_merit       the merit of all parameter sets considered
;
; Author:
;
;   Wing-Fai Thi (SUPA, Institut for Astronomy, Royal Observatory Edinburgh)
;   wingfai.thi at google mail adress  
;
; History:
;  
;   24/04/2007  European Southern Observatory (Garching, Germany)
;               Version 1.0 
;
;   26/04/2007  Change the output to archival_parameters and
;               archival_merit take keep track of all function
;               evaluation for further parameter space statistical
;               analysis
;  
;   27/02/2008  change all simplex into sce_simplex to avoid name
;               space clash
;
;   29/02/2008  add extinction: decrease of the number of complexes with generations, 
;               the worst points being eliminated.
;
;   01/03/2018  public version  
;
; Licence: BSD
; -------------------------------------------------------------------------

; ----- Initialize SCE parameters
nb_minima=1
randomValue = randomu(seed)

nbu=n_elements(bu)
nbl=n_elements(bl)
if not(nbu eq nbl) then begin
    print,'Error in input boundaries'
    stop
endif
nparameters=nbu
gnrng=fltarr(nparameters)
bound = bu - bl
; 
if not keyword_set(max_nshuffle) then max_nshuffle=1000
nshuffle=0

; ----- Check the parameter bounds
for i=0,nparameters-1 do begin
    if bound(i) lt 0.0d0 then begin
        print,'Parameter #',i,'. Upper bound lower than lower bound. Please check!'
        stop
    endif
endfor  

; ----- Select the number of complexes ncomplexes, default = 2

if not keyword_set(min_ncomplexes) then min_ncomplexes = 2
if not keyword_set(ncomplexes) then ncomplexes = 2
if ncomplexes lt 1 then ncomplexes = 2
print,'Number of complexes =',ncomplexes

; ----- Select the number of elements in a complex, default = 2*nparameters+1

if not keyword_set(nelements_complex) then begin
    nelements_complex = 2*nparameters+1
endif    
if nelements_complex lt (2 * nparameters + 1) then begin
    nelements_complex = 2 * nparameters + 1
endif
print,'Number of elements in a complex =',nelements_complex

; ----- Total number of elements nsample

nsample = ncomplexes * nelements_complex

; ----- Number of evolution steps for each complex: nevolution_steps

if not keyword_set(max_nevolution_steps) then max_nevolution_steps = 50  ; arbitrary number
if not keyword_set(nevolution_steps) then begin
    nevolution_steps = nelements_complex
endif
if nevolution_steps lt nelements_complex then nevolution_steps = nelements_complex

; ----- Number of members in a sce_simplex

nelements_sce_simplex = nparameters+1
sce_simplex           = dblarr(nelements_sce_simplex,nparameters)
merit_sce_simplex     = dblarr(nelements_sce_simplex)

; ----- Generate sample
; Sample nsample points in the feasible parameter
; space and compute the criterion value (merit function) at each
; point.  In the absence of prior information, use a uniform probability
; distribution to generate a sample.

parameters          = dblarr(nsample,nparameters)
archival_parameters = dblarr(nsample,nparameters)
for i=0,nsample-1 do begin
    parameters(i,*)          = bl+randomu(seed,nparameters)*bound
    archival_parameters(i,*) = parameters(i,*)
endfor

nloop=0
ncalls=0
merit          = dblarr(nsample)
archival_merit = dblarr(nsample)
for i=0,nsample-1 do begin
    merit(i) = call_function(func,parameters(i,*))
    ncalls = ncalls + 1
end 
archival_merit = merit

; ---- Sort

idx        = sort(merit)
parameters = parameters(idx,*)
merit      = merit(idx)

; ----- Record the best and worst points

best_parameters  = parameters(0,*)
best_merit       = merit(0)
worst_parameters = parameters(nsample-1,*) 
worst_merit      = merit(nsample-1)

; ----- Define the complexes

cmplx       = dblarr(nelements_complex,nparameters)
merit_cmplx = dblarr(nelements_complex)

; ----- Assign a triangular probability distribution
    
m = double(nelements_complex)
proba=2.0d0*(m-indgen(nelements_complex))/m/(m+1.0d0)
proba_range=dblarr(nelements_complex+1)
for i=1,nelements_complex do proba_range(i)=total(proba(0:i-1))

; -----

bl_cmplx = dblarr(nparameters)
bu_cmplx = dblarr(nparameters)
xnstd = dblarr(nparameters)

; ----- Start search loop

success = 0  ; unsuccessful = 0

; -----
ncomplexes_tmp = ncomplexes

while ncalls lt max_func_calls and  success eq 0 and nshuffle lt max_nshuffle do begin

nloop    = nloop + 1
nshuffle = nshuffle + 1

; ----- Loop on complexes (sub-populations)

for icomp = 0, ncomplexes_tmp-1 do begin

   ; ----- Partition into complexes:
   ; Partition the nsample points into ncomplexes
   ; complexes, each containing nelements_complex points.  The complexes are
   ; partitioned in such a way that the first complex contains
   ; every ncomplexes*k+1 ranked point, the second complex contains
   ; every ncomplexes*k+2 ranked point, and so on, where k =
   ; 0,2,...,nelements_complex-1.

   ; print,'Population ',icomp+1,'  /',ncomplexes
   k1 = indgen(nelements_complex)
   k2 = k1*ncomplexes+icomp
   cmplx(k1,*) = parameters(k2,*)
   merit_cmplx = merit(k2)
   for i=0,nparameters-1 do begin
       ave_cmp = mean(cmplx(k1,i))
       sig_cmp = sigma(cmplx(k1,i))
       bl_cmplx(i) = max([ave_cmp - 2. * sig_cmp,bl(i)])
       bu_cmplx(i) = min([ave_cmp + 2. * sig_cmp,bu(i)])
   endfor
    
   ; Competitive Evolution of Sce_Simplexes
   ; Evolve the sub-population icomp for nevolution_steps        

   for evol_step = 0,nevolution_steps-1 do begin

   ; print, 'Generation ',evol_step+1,'   /',nevolution_steps 

   ; Select randomly the sce_simplex by
   ; sampling the complex according to 
   ; a linear probability distribution probability

      selected = replicate(-1,nelements_sce_simplex)
      nselected    = 0
      repeat begin
         rand=randomu(seed)
         w=where(rand gt proba_range)
         candidate=max(w)
         wchosen = where(candidate eq selected,nwchosen)
         if nwchosen eq 0 then begin ; select a candidate 
                                     ; that has not been chosen yet
            selected(nselected) = candidate
            nselected = nselected + 1
         endif
      endrep until (nselected eq nelements_sce_simplex)
           
      ; Order the sce_simplexes  

      merit_sce_simplex = merit_cmplx(selected)
      order = sort(merit_sce_simplex)
      selected_order = selected(order)
      sce_simplex(indgen(nelements_sce_simplex),*)=cmplx(selected_order,*)  
      merit_sce_simplex = merit_sce_simplex(order)
           
      ; Generate a new sce_simplex point

	if ncomplexes_tmp gt 1 then begin
         generate_offspring,func,sce_simplex,merit_sce_simplex,bl,bu, $
                         bl_cmplx,bu_cmplx,ncalls,sce_simplex_new,merit_new, $
                         archival_parameters,archival_merit
      endif else begin
         generate_offspring,func,sce_simplex,merit_sce_simplex,bl,bu, $
                         bl_cmplx,bu_cmplx,ncalls,sce_simplex_new,merit_new, $
                         archival_parameters,archival_merit; $
                         ;,/pure_simplex,/barycenter
	endelse

      ; Replace the worst point in Sce_Simplex with the new point:
      
      sce_simplex(nelements_sce_simplex-1,*)     = sce_simplex_new 
      merit_sce_simplex(nelements_sce_simplex-1) = merit_new
    
      ; Replace the sce_simplex into the complex;

      cmplx(selected_order,*)     = sce_simplex
      merit_cmplx(selected_order) = merit_sce_simplex
       
      ; Sort the complex
 
      order       = sort(merit_cmplx)
      cmplx(k1,*) = cmplx(order,*)
      merit_cmplx = merit_cmplx(order)

   ; end of inner loop for Competitive Evolution of Sce_Simplexes
   endfor
       
   ; Replace the complex back into the population;
   parameters(k2,*) = cmplx(k1,*)
   merit(k2) = merit_cmplx

; end of loop on Complex Evolution
endfor  

; ----- Shuffle/Rank points
; Sort the nsample points in order of increasing
; criterion value so that the first point represents the point
; with the smallest criterion value (best) and the last point
; represents the point with the largest criterion value (worst).

idx        = sort(merit)
parameters = parameters(idx,*)
merit      = merit(idx)

; ----- Record the best and worst points

nsample = nelements_complex*ncomplexes_tmp

best_parameters  = parameters(0,*)
best_merit       = merit(0)
worst_parameters = parameters(nsample-1,*) 
worst_merit      = merit(nsample-1)

;
; ----- Compute the standard deviation for each parameter xnstd
;       and the normalized geometric range of the parameters gnrng

for i=0,nparameters-1 do begin
    xnstd(i)=sigma(parameters(0:nsample-1,i))
    gnrng(i)=exp(mean(alog((max(parameters(0:nsample-1,i))-min(parameters(0:nsample-1,i)))/bound(i))))  
endfor

; ---- Define the criterion for finding the optimium
; There are alternative stopping criterion
;
; theses criteria may not be appropriate for the step function test
crit = min(gnrng)
crit = median(gnrng) ; seems to work better with step, not with Schwefel, Corana

;
; F-test method to adapt the paremeters
; ftol = nparameters/(ndata-nparameters)*f_cvf(confidence,nparameters,ndata-nparameters)

if nshuffle gt 1 then begin
;   crit = abs(previous_best_merit-best_merit)
   if crit gt 0.0d0 and crit lt ftol then success = 1
endif

previous_best_merit=best_merit

; ------- Extinction
; If this option is chosen, the complex with the worst points is eliminated after "extinction" shuffling 
; steps until min_ncomplexes are left.

if keyword_set(extinction) and (ncomplexes_tmp gt min_ncomplexes) then begin
	        if((nshuffle mod extinction) eq 0) then ncomplexes_tmp = ncomplexes_tmp -1
endif
; -------- Increase the number of evolutionary steps if there is only one complex left

if (ncomplexes_tmp eq 1) then nevolution_steps = max_nevolution_steps

; ----- Verbose  

if keyword_set(verbose) then begin
print
print,'shuffle=',nshuffle
print,'number of complexes=',ncomplexes_tmp
print,'Number of function calls =',ncalls,' best =',best_merit,' crit =',crit
print,'best parameters =',best_parameters
endif
; end of the search 

endwhile
; --------------------------------

best_parameters  = parameters(0,*)
best_merit       = merit(0)

; ----- Possibility of multi minima
;       1 march 2008

nbest = ncomplexes_tmp
for i=1,nbest do begin
    test=(abs(merit(i,*)-merit(0,*))/merit(0,*))
    if (test(0) lt 1d-3) then begin
       ; Possible multi global minima
         nb_minima=nb_minima+1
    endif
endfor    


; ----------------------------------
; Display the results

;print
;print,'Random generator seed:',seed

if success eq 1 then begin
   print 
   print,'Convergence criterion :',crit
   print
   print,'number of function calls:',ncalls
   print
   print,'Results'
   print,'Parameters:'
   print,best_parameters
   print,'Evaluation function=',best_merit 
   if(nb_minima gt 1) then begin
       print,'Possibility of',nb_minima,' global+local minima'
   endif    
endif
; -----
if ncalls ge max_func_calls then begin
    print,'Optimization search terminated because the limit'
    print,'on the maximum number of trials'
    print, max_func_calls
    print,'has beem exceeded. Search was stopped at trial number:'
    print,ncalls
    print,'of the initial loop!'
;    stop 
endif
; -----
jsort = sort(archival_merit)
narchival = indgen(n_elements(jsort))
archival_parameters(narchival,*) = archival_parameters(jsort,*)
archival_merit                   = archival_merit(jsort)

;stop
end

; ----------------------------------------------------------------------
pro generate_offspring,func,s,sf,bl,bu,bl_cmplx,bu_cmplx,ncalls,snew,fnew $
                      ,archival_parameters,archival_merit $
                      ,alpha=alpha,beta=beta,gamma=gamma,barycenter=barycenter $
                      ,expansion=expansion,pure_simplex=pure_simplex

;  This is the subroutine for generating a new point in a sce_simplex
;
;   s(.,.) = the sorted sce_simplex in order of increasing function values
;    sf(.) = function values in increasing order
;    bl(.) = the lower bound of the parameters
;    bu(.) = the upper bound of the parameters
;    cmplx
;    ncalls = number of function calls
;
;    nsce_simplex = number of members in a sce_simplex
; 
; LIST OF LOCAL VARIABLES
;   sb(.) = the best point of the sce_simplex
;   sw(.) = the worst point of the sce_simplex
;   fw = function value of the worst point
;   ce(.) = the centroid of the sce_simplex excluding wo
;   snew(.) = new point generated from the sce_simplex
;   fnew(.) = function values at the new point
;
; Author:
;
;   Wing-Fai Thi (SUPA, Institut for Astronomy, Royal Observatory Edinburgh)
;
; History
;   24/04/2007  European Southern Observatory (Garching, Germany)
;               Version 1.0 
;   26/04/2007  Change the output to archival_parameters and
;               archival_merit take keep track of all function
;               evaluation for further parameter space statistical analysis
;
;

dimension=size(s,/dimensions)
nsce_simplex= dimension(0)
nparameters=dimension(1)

; if not user-define, use the default Nelder-Mead coefficients

if not keyword_set(alpha) then alpha = 1.0
if not keyword_set(beta)  then beta  = 0.5
if not keyword_set(gamma) then gamma = 1.5

; Assign the best and worst points

sb = reform(s(0,*),nparameters)     & fb=sf(0)
sw = reform(s(nsce_simplex-1,*),nparameters) & fw=sf(nsce_simplex-1)

; Compute the centroid of the sce_simplex excluding the worst point
ce=dblarr(nparameters)

if not keyword_set(barycenter) then begin
    for i=0,nparameters-1 do ce(i)=mean(s(0:nsce_simplex-2,i))
endif else begin
; Weighted centroid
    weights = 1.0d0/sf(0:nsce_simplex-2)
    for i=0,nparameters-1 do ce(i)=total(weights*s(0:nsce_simplex-2,i))/total(weights)
endelse

fnew = fw + 100.0d0

; ----------------------
if keyword_set(pure_simplex) and keyword_set(expansion) then begin
   ; Attempt an expansion point, the expansion tends to overshoot

   snew = ce + gamma*(ce-sw)

   ; Compute the cost function at the expansion point if inside the bounds
   
   if (min(snew-bl) ge 0.0 and min(bu-snew) ge 0.0) then begin
     fnew = call_function(func,snew)
     archival_parameters = [archival_parameters,transpose(snew)]
     archival_merit = [archival_merit,fnew]
     ncalls = ncalls + 1
   endif

endif    
; ----------------------

; Attempt a reflection point if expansion point fails

if fnew gt fw then begin

; ------------------------------------------------------
    snew = ce + alpha*(ce-sw)

    ; Check if is outside the bounds:
 
    if (min(snew-bl) lt 0.0 or min(bu-snew) lt 0.0) then begin
      snew = bl_cmplx + randomu(seed,nparameters)*(bu_cmplx-bl_cmplx)
     ;print,'Reflection outside boundaries for the parameters!'
    endif 

    fnew = call_function(func,snew)
    archival_parameters = [archival_parameters,transpose(snew)]
    archival_merit = [archival_merit,fnew]
    ncalls = ncalls + 1

    ; Reflection failed; now attempt a contraction point
    if fnew gt fw then begin
       snew = sw + beta*(ce-sw)
       fnew = call_function(func,snew)
       archival_parameters = [archival_parameters,transpose(snew)]
       archival_merit = [archival_merit,fnew]
       ncalls = ncalls + 1
       ; Both reflection and contraction have failed, attempt a random point; 
       if fnew gt fw then begin
	    if keyword_set(pure_simplex) then begin
            ; inner contraction point (original Nelder & Mead)
            snew = sw - beta*(ce-sw)
          endif else begin    
          ; attempt a random point (original SCE) 
            snew = bl_cmplx + randomu(seed,nparameters)*(bu_cmplx-bl_cmplx) 
          ; attempt a random inner contraction point
;           snew = sw - beta*randomu(seed,nparameters)*(ce-sw)
          endelse
          ; -----
          archival_parameters = [archival_parameters,transpose(snew)]
          fnew = call_function(func,snew)
          archival_merit = [archival_merit,fnew]
;          print,'Reflection and Contraction failed! Choose a random point.'
          ncalls = ncalls + 1
       endif
   endif

endif
end
