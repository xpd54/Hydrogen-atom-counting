'''Copyright (C) 2014 Ravi Prakash ( vickyravi17@gmail.com )@ LNM Insittute Of Information Technology
All editing rights are reserved to the author and people related to this project.
Any academic or nonacademic use of this source code without informing author or related people is illegal.
Tested and Debugged on "Python 2.7.3 [GCC 4.6.3] on linux2" '''



''' 
	programm for solving Quartic Equation.
	It will read input from " inputSolutionForCubicSc.txt ".
	And write ans in " SolutionForCubicSc.txt "
'''

import math
R = 8.3144621; # R in J/K mol ...... last L is for long double
A0 = 0.02; # A0 in Pa.m^3/mol^2
a = -5060000; # a in m^3/mol
B0 = 20960000; # m^3/mol
b = -43590000; # m^3/mol
C = 5.04; #m^3K^3/mol

def complex_sqroot(argc):
	a=complex(argc)
	sqroot = math.sqrt((a.real**2)+(a.imag**2));
	y= math.sqrt((sqroot-a.real)/2.0);
	x= a.imag/(2.0*y);
	root1 = complex(x,y)
	root2 = complex(-x,-y)
	return (root1,root2)

def choose_numbers(a,b,c):
	if(a.real == 0 and a.imag == 0):
		return ( b , c )
	if(b.real == 0 and b.real == 0):
		return ( a , b )
	if(c.real == 0 and c.imag == 0):
		return ( a , b )
	if(a.imag!=0 and b.imag!=0):
		return ( a , b )
	if(b.imag!=0 and c.imag!=0):
		return ( b , c )
	if(c.imag !=0 and a.imag!=0):
		return ( c , a )
	else: # Error message just for information
		print("Error ********Something wrong with logic in quartic equation*******\n");
		print("Error *********** Program is terminating*****\n");
		exit(-1);
	
		

def get_root_of_cubic_equation(a,b,c,d):
	f = ( ( (3.0*c)/a ) - ((b**2)/(a**2)) ) / 3.0;
	g = ( ( (2.0*(b**3)) / (a**3) ) - ((9.0*b*c)/(a**2)) + ((27.0*d) / a ) ) / 27.0;
	h = ((g**2) / 4.0) + ((f**3) / 27.0);
	if(h < 0): # all three roots are real
		i = math.sqrt(((g**2) / 4.0) - h);
		j = i**(1.0/3.0);
		k = math.acos(-(g/(2.0*i)));
		l = -j;
		m = math.cos((k/3.0));
		n = math.sqrt(3.0) * math.sin(k/3.0);
		p = -(b/(3.0*a));
		x1 = ( 2.0 * j *math.cos(k/3.0) ) - (b/(3.0*a));
		x2 = ( l * (m+n) ) + p;
		x3 = ( l * (m-n) ) + p;
		first = complex(x1,0.0);   # all root are real so imaginary value are zero
		second = complex(x2,0.0);
		third = complex(x3,0.0)
		return (first,second,third);

	if(h > 0): # only one root is real
		r = -(g/2.0) + math.sqrt(h);
		s = r**(1.0/3.0);
		t = -(g/2.0) - math.sqrt(h);
		if(t<0):
			t = t * (-1);
			u = (-1)*t**(1.0/3.0);
		else:
			u = t**(1.0/3.0);
		x1 = complex( ( (s+u) - (b/(3.0*a)) ),0.0);
		x2 = complex(  ( -((s+u) / 2.0) - (b / (3.0*a ) ) ) , ( ((s-u) * math.sqrt(3.0))/2.0 ) );
		x3 = complex( (-((s+u) / 2.0) - (b / (3.0*a))) , ( -(((s-u) * math.sqrt(3.0))/2.0) ) );
		return (x1,x2,x3);


	if(h==0 and g==0 and f==0): # one special case when h==0 it will happen only when g and f will also be 0
		x1 = x2 = x3 = complex( ( -((d/a)**(1.0/3.0)) ) , 0.0);
		return (x1,x2,x3);
	else: # Error message just for information
		print("Error ********Something wrong with logic in cubic equation*******\n");
		print("Error *********** Program is terminating*****\n");
		exit(-1);
		
f = open('inputSolutionForCubicSc.txt','r');
constant = f.readline();
constant = float(constant);
constant = constant+2;
constant = -constant;
(x1,x2,x3)=get_root_of_cubic_equation(1,-3,2,constant);
print "1. " + str(x1)
print "2. " + str(x2)
print "3. " + str(x3)
print "********** Number Of Unit Cubes On A Given Edge ********"
print "********** j is Showing Imaginary Value*******"
option = int(raw_input());
if(option == 1):
	take = x1 ;
if(option == 2):
	take = x2 ;
if(option == 3):
	take = x3 ;
real_value = take.real;
f = open('SolutionForCubicSc.txt','w');
f.write(str(real_value));
f.close();
