'''Copyright (C) 2014 Ravi Prakash ( vickyravi17@gmail.com )@ LNM Insittute Of Information Technology
All editing rights are reserved to the author and people related to this project.
Any academic or nonacademic use of this source code without informing author or related people is illegal.
Tested and Debugged on "Python 2.7.3 [GCC 4.6.3] on linux2" '''



''' 
	programm for solving Quartic Equation.
	It will read input from " data.txt ".
	And write ans in " output.txt "
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
	

def get_root_of_quartic_equation( a, b, c, d, e):
	a = float(a)
	b = b / a;
	c = c / a;
	d = d / a;
	e = e / a;
	a = a / a;
	#print (a,b,c,d,e);
	f = c - ((3.0*(b**2))/8.0);
	g = d + ((b**3)/8.0) - ((b*c)/2.0);
	h = e - ( (3.0 * (b**4))/256.0 ) + (( (b**2) * c)/16.0) - ( (b*d)/4.0 );
	#print (f,g,h)
	q_a = a;
	q_b = (f/2.0);
	q_c = ((f*f - 4.0*h)/16.0);
	q_d = - (g*g) / 64.0;
	#print (q_a,q_b,q_c,q_d);
	(first,second,third) = get_root_of_cubic_equation(q_a,q_b,q_c,q_d);
	#print (first,second,third);
	(temp1,temp2) = choose_numbers(first,second,third);
	#print (temp1,temp2);
	(sq_p_first,sq_p_second) = complex_sqroot(temp1);
	(sq_q_first,sq_q_second) = complex_sqroot(temp2);
	if(sq_p_first.real > 0):
		p = sq_p_first;
	else:
		p = sq_p_second;
	if(sq_q_first.real > 0):
		q = sq_q_first;
	else:
		q = sq_q_second;
	temp = (p*q);
	r = complex( (-g / (8.0 * temp.real)) , 0.0 );
	s = complex( (b / (4.0 * a)) , 0.0);

	# x1 = p+q+r-s
	x1 = p+q+r-s;
	# x2 = p-q-r-s
	x2 = p-q-r-s;
	#x3 = -p+q-r-s
	x3 = -p+q-r-s;
	#x4 = -p-q+r-s
	x4 = -p-q+r-s;
	return (x1,x2,x3,x4);
#reading data from data.txt file written by C file

f = open('data.txt','r');
value = []
for val in f.read().split():
    value.append(float(val))
f.close()

(y1,y2,y3,y4) = get_root_of_quartic_equation(value[0],value[1],value[2],value[3],value[4]);
(x1,x2,x3,x4) = (y1,y2,y3,y4);
print "1. " + str(y1)
print "2. " + str(y2)
print "3. " + str(y3)
print "4. " + str(y4)
print "**********Select The Value for volume ********"
print "********** j is Showing Imaginary Value*******"
option = int(raw_input());
if(option == 1):
	take = x1 ;
if(option == 2):
	take = x2 ;
if(option == 3):
	take = x3 ;
if(option == 4):
	take = x4 ;
real_value = take.real;
f = open('output.txt','w');
f.write(str(real_value));
f.close();
print "Value You Have Selected "+str(real_value);
'''(x1,x2,x3,x4) = get_root_of_quartic_equation(7423.538334,-1937.269669,-40605172268.507228,-1769979459185018141.125000,705230359053.054166);
#(x1,x2,x3,x4) = get_root_of_quartic_equation(-20,5,17,-29,87);
#(x1,x2,x3) = get_root_of_cubic_equation( 1,-0.436718750000,1.113366088867,-0.028131544590 ); #cubic working fine
#(x1,x2,x3,x4) = get_root_of_quartic_equation(3,6,-123,-126,1080);
#(x1,x2,x3) = get_root_of_cubic_equation( 1,-21.25,20.25,0 )
print(x1,x2,x3,x4);'''
