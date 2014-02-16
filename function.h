/*Copyright (C) 2014 Ravi Prakash ( vickyravi17@gmail.com )@ LNM Insittute Of Information Technology
All editing rights are reserved to the author and people related to this project.
Any academic or nonacademic use of this source code without informing author or related people is illegal.
Tested and Debugged on "gcc version 4.6.3 (Ubuntu/Linaro 4.6.3-1ubuntu5)"*/



/*
	Due to precision error in C we are using python for calculating quartic equation.
	C code for this equation is written in this file and python is written in function.py .
	Both C file and python file are connected with file IO data.txt and output.txt
	
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define pi 3.14
const long double R = 8.3144621; // R in J/K mol ...... last L is for long double
const long double A0 = 0.02; // A0 in Pa.m^3/mol^2
const long double a = -5060000; // a in m^3/mol
const long double B0 = 20960000; // m^3/mol
const long double b = -43590000; // m^3/mol
const long double C = 5.04; //m^3K^3/mol

//............ complex structure .......//
struct complex_number
{
	long double real;
	long double imaginary;
};

struct pair_of_complex_numbers
{
	struct complex_number first;
	struct complex_number second;
};

struct three_complex_numbers
{
	struct complex_number first;
	struct complex_number second;
	struct complex_number third;
};

struct four_complex_numbers
{
	struct complex_number first;
	struct complex_number second;
	struct complex_number third;
	struct complex_number fourth;
};

struct complex_number complexAdd(struct complex_number a,struct complex_number b) // for a+b of complex numbers
{
	struct complex_number result;
	result.real = (a.real + b.real);
	result.imaginary = (a.imaginary + b.imaginary);
	return result;
}

struct complex_number complexSub(struct complex_number a,struct complex_number b) // for a-b of complex numbers
{
	struct complex_number result;
	result.real = a.real - b.real;
	result.imaginary = a.imaginary - b.imaginary;
	return result;
}

//************* complex number multiplication ***********//
	/*for a*b of complex numbers
	(5 + 6i) * (7 + 8i) This equals
	35 + 40i + 42i + 48i2 As we saw above,
	i2 = -1 so 48i2 = -48 So answer= -13 + 82i*/


struct complex_number complexMul(struct complex_number a,struct complex_number b)
{
	long double first = a.real * b.real;
	long double second = a.real * b.imaginary;
	long double third = a.imaginary * b.real;
	long double fourth = a.imaginary * b.imaginary;
	fourth = -fourth;
	struct complex_number result;
	result.real = first + fourth;
	result.imaginary = second + third;
	return result;
}


//************* complex number division *****************//

/*		(9 + 3i) ÷ (7 + 5i)
	Multiplying top and bottom by the conjugate:
	((9 + 3i)* (7 - 5i)) ÷ ((7 + 5i) * (7 - 5i))
	Which equals (78 - 24 i) ÷ 74
	Equals (78 ÷ 74) - (24i ÷ 74)
	Answer = 1.054054054054054 & -0.32432432432432434 i */

struct complex_number complexDivi(struct complex_number a,struct complex_number b)
{
	struct complex_number conjugate_b;
	conjugate_b.real = b.real;
	conjugate_b.imaginary = -b.imaginary;
	struct complex_number temp = complexMul(a,conjugate_b);
	struct complex_number conjugate_mulB = complexMul(b,conjugate_b);
	struct complex_number result;
	result.real = temp.real / conjugate_mulB.real;
	result.imaginary = temp.imaginary / conjugate_mulB.real;
	return result;
}


//**********complex sqroot *****************//
	/*Find the square root of 12 + 16i.
	r = Square Root (122 + 162)
	r = Square Root (144 + 256) = 20
	y = Square Root ((20-12)/2) = 2
	x = 16/(2*2) = 4
	root 1 = 4 + 2i
	root 2 = -4 - 2i*/

struct pair_of_complex_numbers complex_sqroot(struct complex_number a)// return pair of a complex number which has two root first and second
{
	long double sqroot;
	sqroot = sqrtl(powl(a.real,2.0L)+powl(a.imaginary,2.0L));
	long double y= sqrtl((sqroot-a.real)/2.0L);
	long double x= a.imaginary/(2.0L*y);
	struct complex_number root1;
	struct complex_number root2;
	root1.real = x;
	root1.imaginary = y;
	root2.real = -x;
	root2.imaginary = -y;
	struct pair_of_complex_numbers result;
	result.first = root1;
	result.second = root2;
	return result;
}

struct three_complex_numbers get_root_of_cubic_equation(long double a,long double b,long double c,long double d)
{
	//printf("%Lf %Lf %Lf %Lf\n",a,b,c,d);
	long double f = ( ( (3.0L*c)/a ) - (powl(b,2.0L)/powl(a,2.0L) ) ) / 3.0L;
	long double g = ( ( (2.0L*powl(b,3.0L)) / powl(a,3.0L) ) - ((9.0L*b*c)/powl(a,2.0L)) + ((27.0L*d) / a ) ) / 27.0L;
	long double h = (powl(g,2.0L) / 4.0L) + (powl(f,3.0L) / 27);
	//printf("%Lf %Lf %Lf\n",f,g,h);
	if(h < 0) // all three roots are real
	{
		long double i = sqrtl((powl(g,2.0L) / 4.0L) - h);
		long double j = cbrtl(i);
		long double k = acosl(-(g/(2.0L*i)));
		long double l = -j;
		long double m = cosl((k/3.0L));
		long double n = sqrtl(3.0L) * sinl(k/3.0L);
		long double p = -(b/(3.0L*a));
		long double x1 = ( 2.0L * j *cosl(k/3.0L) ) - (b/(3.0L*a));
		long double x2 = ( l * (m+n) ) + p;
		long double x3 = ( l * (m-n) ) + p;
		//printf("%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf\n",i,j,k,l,m,n,p,x1,x2,x3);
		struct three_complex_numbers result;
		result.first.real = x1;
		result.first.imaginary = 0;   // all root are real so imaginary value are zero
		result.second.real = x2;
		result.second.imaginary = 0;
		result.third.real = x3;
		result.third.imaginary = 0;
		return result;
	}

	if(h > 0) // only one root is real
	{
		long double r = -(g/2.0L) + sqrtl(h);
		long double s = cbrtl(r);
		long double t = -(g/2.0L) - sqrtl(h);
		long double u = cbrtl(t);
		long double x1 = (s+u) - (b/(3.0L*a));
		struct complex_number x2;
		x2.real = -((s+u) / 2.0L) - (b / (3.0L*a));
		x2.imaginary = ((s-u) * sqrtl(3.0L))/2.0L;
		struct complex_number x3;
		x3.real = -((s+u) / 2.0L) - (b / (3.0L*a));
		x3.imaginary = -(((s-u) * sqrtl(3.0L))/2.0L);
		//printf("%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf\n",f,g,h,r,s,t,u,x1,x2.real,x2.imaginary,x3.real,x3.imaginary);
		struct three_complex_numbers result;
		result.first.real = x1;
		result.first.imaginary = 0; // first root is real so imaginary value is zero
		result.second = x2;
		result.third = x3;
		return result;
	}

	if(h==0 && g==0 && f==0) // one special case when h==0 it will happen only when g and f will also be 0
	{
		long double x1,x2,x3;
		x1 = x2 = x3 = (cbrtl((d/a)) * (-1));
		struct three_complex_numbers result;
		result.third.real=result.second.real=result.first.real=x1; // all roots are equal and imaginary value are zero
		result.third.imaginary=result.second.imaginary=result.first.imaginary=0;
		return result;
	}
	else // Error message just for information
	{
		printf("Error ********Something wrong with logic in cubic equation*******\n");
		printf("Error *********** Program is terminating*****\n");
		exit(-1);
	}
}

struct four_complex_numbers get_root_of_quartic_equation(long double a,long double b,long double c,long double d,long double e)
{
	b = b / a;
	c = c / a;
	d = d / a;
	e = e / a;
	a = a / a;
	long double f = c - ((3.0L*powl(b,2.0L))/8.0L);
	long double g = d + (powl(b,3.0L)/8.0L) - ((b*c)/2.0L);
	long double h = e - ( (3.0L * powl(b,4.0L))/256.0L ) + ( powl(b,2.0L) * (c/16.0L)) - ( (b*d)/4.0L );
	//printf("%Lf %Lf %Lf\n",f,g,h);
	long double q_a = a;
	long double q_b = (f/2.0L);
	long double q_c = ((f*f - 4*h)/16);
	long double q_d = - g*g / 64.0L;
	//printf("%Lf\n",q_c);
	struct three_complex_numbers cubic = get_root_of_cubic_equation(q_a,q_b,q_c,q_d);
	//printf("%Lf %Lf %Lf %Lf %Lf %Lf \n",cubic.first.real,cubic.first.imaginary,cubic.second.real,cubic.second.imaginary,cubic.third.real,cubic.third.imaginary);
	struct complex_number p,q,r,s;
	if((cubic.first.imaginary==0.0L) && (cubic.second.imaginary!=0) && (cubic.third.imaginary !=0) )
	{
		struct pair_of_complex_numbers sq_p = complex_sqroot(cubic.second);
		struct pair_of_complex_numbers sq_q = complex_sqroot(cubic.third);
		if(sq_p.first.real > 0)
			p = sq_p.first;
		else
			p = sq_p.second;
		if(sq_q.first.real > 0)
			q = sq_q.first;
		else
			q = sq_q.second;
		struct complex_number temp;
		temp = complexMul(p,q);
		r.real = -g / (8.0L * temp.real);
		r.imaginary = 0;
		s.real = b / (4.0L * a);
		s.imaginary = 0;
		//printf("%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf\n",p.real,p.imaginary,q.real,q.imaginary,r.real,r.imaginary,s.real,s.imaginary);
	}
	else // Error message just for information
	{
		printf("Error ********Something wrong with logic in quartic equation*******\n");
		printf("Error *********** Program is terminating*****\n");
		exit(-1);
	}

	struct complex_number x1,x2,x3,x4;
	// x1 = p+q+r-s
	{
		x1 = complexAdd(p,q);
		x1 = complexAdd(x1,r);
		struct complex_number temp;
		temp.real = -s.real;
		temp.imaginary = -s.imaginary;
		x1 = complexAdd(x1,temp);
	}
	// x2 = p-q-r-s
	{
		struct complex_number temp;
		temp.real = -q.real;
		temp.imaginary = -q.imaginary;
		x2 = complexAdd(p,temp);
		temp.real = -r.real;
		temp.imaginary = -r.imaginary;
		x2 = complexAdd(x2,temp);
		temp.real = -s.real;
		temp.imaginary = -s.imaginary;
		x2 = complexAdd(x2,temp);
	}
	//x3 = -p+q-r-s
	{
		struct complex_number temp;
		temp.real = -p.real;
		temp.imaginary = -p.imaginary;
		x3 = complexAdd(temp,q);
		temp.real = -r.real;
		temp.imaginary = -r.imaginary;
		x3 = complexAdd(x3,temp);
		temp.real = -s.real;
		temp.imaginary = -s.imaginary;
		x3 = complexAdd(x3,temp);
	}
	//x4 = -p-q+r-s
	{
		struct complex_number temp;
		temp.real = -p.real;
		temp.imaginary = -p.imaginary;
		struct complex_number temp1;
		temp1.real = -q.real;
		temp1.imaginary = -q.imaginary;
		x4 = complexAdd(temp,temp1);
		x4 = complexAdd(x4,r);
		temp.real = -s.real;
		temp.imaginary = -s.imaginary;
		x4 = complexAdd(x4,temp);
	}
	struct four_complex_numbers result;
	result.first = x1;
	result.second = x2;
	result.third = x3;
	result.fourth = x4;
	//printf("%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf\n",result.first.real,result.first.imaginary,result.second.real,result.second.imaginary,result.third.real,result.third.imaginary,result.fourth.real,result.fourth.imaginary);
	return result;
}

struct return_of_Sc
{
	long double tank_volume_requared;
	long double no_of_unit_cell_required_to_design_the_tank;
};

struct return_of_Sc SC_arrangement(long double bubble_radius,long double no_of_bubble_in_tank)
{
	long double perct_of_touching_radius = 90.0;
	perct_of_touching_radius = perct_of_touching_radius / 100.0;
	long double unit_cube_edge_length;
	unit_cube_edge_length = (2.0*bubble_radius) + ((2.0*(1.0-perct_of_touching_radius))*bubble_radius); // in mm
	unit_cube_edge_length = unit_cube_edge_length / 1000.0; // in meter
	long double volume_of_unit_cube = powl(unit_cube_edge_length,3.0);
	long double no_of_bubble_per_unit = 1.0;
	long double tank_volume_requared = volume_of_unit_cube*no_of_bubble_in_tank;
	long double no_of_unit_cell_required_to_design_the_tank = tank_volume_requared / volume_of_unit_cube;
	struct return_of_Sc result;
	result.tank_volume_requared = tank_volume_requared;
	result.no_of_unit_cell_required_to_design_the_tank = no_of_unit_cell_required_to_design_the_tank;
	return result;
}

struct return_of_Bc
{
	long double tank_volume_requared;
	long double no_of_unit_cell_required_to_design_the_tank;
};

struct return_of_Bc BC_arrangement(long double bubble_radius,long double no_of_bubble_in_tank)
{
	long double unit_cube_digonal_length;
	unit_cube_digonal_length = (4.0*bubble_radius) * (1+(1-0.9)) ;
	long double unit_cube_edge_length = unit_cube_digonal_length / sqrtl(3.0);
	unit_cube_edge_length /=1000.0;
	long double volume_of_unit_cube = powl(unit_cube_edge_length,3.0);
	long double no_of_bubble_per_unit_cube = 2.0;
	long double tank_volume_requared = (volume_of_unit_cube * no_of_bubble_in_tank) / 2.0;
	long double no_of_unit_cell_requared_to_design_the_tank = tank_volume_requared / volume_of_unit_cube;
	struct return_of_Bc result;
	result.no_of_unit_cell_required_to_design_the_tank = no_of_unit_cell_requared_to_design_the_tank;
	result.tank_volume_requared = tank_volume_requared;
	return result;
}

struct return_of_Fc
{
	long double tank_volume_raquared;
	long double no_of_unit_cell_required_to_design_the_tank;
};

struct return_of_Fc FC_arrangement(long double bubble_radius,long double no_of_bubble_in_tank)
{
	long double unit_cube_digonal_length;
	unit_cube_digonal_length = (4.0*bubble_radius) + ((4.0*bubble_radius)*(1-0.9));
	long double unit_cube_edge_length = unit_cube_digonal_length / sqrtl(2.0);
	unit_cube_edge_length /= 1000.0;
	long double volume_of_unit_cube = powl(unit_cube_edge_length,3.0);
	long double no_bubble_per_unit_cube = 4.0;
	long double tank_volume_requared = ( volume_of_unit_cube * no_of_bubble_in_tank ) / 4.0;
	long double no_of_unit_cell_required_to_design_the_tank = tank_volume_requared / volume_of_unit_cube;
	struct return_of_Fc result;
	result.tank_volume_raquared = tank_volume_requared;
	result.no_of_unit_cell_required_to_design_the_tank = no_of_unit_cell_required_to_design_the_tank;
	return result;
}
