/*Copyright (C) 2014 Ravi Prakash ( vickyravi17@gmail.com )@ LNM Insittute Of Information Technology
All editing rights are reserved to the author and people related to this project.
Any academic or nonacademic use of this source code without informing author or related people is illegal.
Tested and Debugged on "gcc version 4.6.3 (Ubuntu/Linaro 4.6.3-1ubuntu5)"*/



/*

	To compile and run this code you need gcc and python install in you linux machine.
	This code won't work on windows cause it's uses linux system command line .
	To compile type "gcc main.c -lm" and to run " ./a.out "
	After that follow the terminal instruction .
	For better understanding open file "Project calcuations.xlsx".

*/


#include"function.h"
int main()
{
	long double MW = 2.0158L;// in g/mol
	printf("M.W in Kg/mol := %Lf\n",MW/1000); // in Kg/mol
	long double P; // pressure from user in kpa
	long double T; // Temperature from user in kelvin
	printf("Enter value of pressure\n");
	scanf("%Lf",&P);
	//P = 7423.538334;
	printf("Enter value of Temperature\n");
	scanf("%Lf",&T);
	//T = 233.0L;


	//****** for calculation of Quartic Equations  we need calculate coefficients of equation******//


	const long double Alpha = P;
	const long double Beta = -(R*T);
	const long double Gamma = (-(R*T*B0)) + ((R*C)/(powl(T,2))) +(A0);
	const long double Zeta = ((R*C*B0)/(powl(T,2)))+ (R*T*b*B0) - (a*A0);
	const long double Delta = (-(R*C*b*B0)/(powl(T,2)));

	// ********* these coefficients are defined as given value *******//
	// ******** these make an equation as      αV^4+βV^3+γV^2+δV+ξ= 0  coefficients are as Alpha Beta Gamma Zeta and Delta from right to left
	FILE *file;
	file = fopen("data.txt","w");
	fprintf(file,"%.30Lf\n%.30Lf\n%.30Lf\n%.30Lf\n%.30Lf",Alpha,Beta,Gamma,Zeta,Delta);
	//printf("%Lf %Lf %Lf %Lf %Lf\n",Alpha,Beta,Gamma,Zeta,Delta);
	fclose(file);
	// running python file and getting value form written file
	system("python function.py");
	// reading volume from file created and calculated by python
	long double volume;
	file = fopen("output.txt","r");
	if(file==NULL)
		printf("Error.....Something went wrong file not opened\n");
	else
		fscanf(file,"%Lf",&volume);
	printf("Volume You Have Selected = %Lf\n",volume);
	//long double density = (P*MW)/(R*T); // comment it to calculate by equation
	long double density = (MW/volume); // uncomment it to calculate by equation

	printf("Density is = %Lf\n",density);
	printf("Enter Fixed Mass Of H2 In Tank To Be Filled\n");
	long double fixed_mass_of_h2;
	// input for fixed mass of H2
	scanf("%Lf",&fixed_mass_of_h2);
	printf("Enter Bubble Radius(r)\n");
	long double r;
	// input for Bubble Radius
	scanf("%Lf",&r);
	long double bubble_volume;
	bubble_volume = (4.0L/3.0L) * pi * (r*r*r);
	//bubble_volume = (4.0L/3.0L) * M_PI * (r*r*r); // ask to instructor what to do with pi
	printf("Bubble Volume = %Lf mm^3\n",bubble_volume);
	bubble_volume = bubble_volume*pow(10,-9);
	printf("Bubble Volume = %.15LF m^3\n",bubble_volume);
	long double mass_of_hydrogen_per_bubble;
	mass_of_hydrogen_per_bubble = bubble_volume * density;
	printf("Mass Of Hydrogen Per Bubble = %.15Lf kg\n",mass_of_hydrogen_per_bubble);
	long double no_of_bubble_in_tank;
	no_of_bubble_in_tank = fixed_mass_of_h2 / mass_of_hydrogen_per_bubble;//did ceiling ask to instructor
	printf("No Of Bubble In Tank = %Lf\n",no_of_bubble_in_tank);
	long double total_bubble_volume = no_of_bubble_in_tank * bubble_volume;
	printf("Total Bubble Volume = %.15Lf m^3\n",total_bubble_volume);
	long double density_of_structure_material;
	printf("Enter Density Of Structure Material\n");
	// input for Density Of Structure Material
	scanf("%Lf",&density_of_structure_material);
	int user_arr;
	printf("*************Enter Option For SC BCC FCC***********\n");
	printf("1. SC\n2. BCC\n3. FCC\n");
	scanf("%d",&user_arr);
	if(user_arr == 1) // for SC
	{
		struct return_of_Sc Sc = SC_arrangement(r,no_of_bubble_in_tank);
		long double total_tank_volume = Sc.tank_volume_requared;// structure + bubble
		printf("Total Tank Volume = %Lf m^3\n",total_tank_volume);
		long double total_volume_of_structure;
		total_volume_of_structure = (total_tank_volume - total_bubble_volume);
		printf("Total Volume Of Structure = %.15Lf m^3\n",total_volume_of_structure);
		long double mass_of_structure = density_of_structure_material * total_volume_of_structure;
		printf("Mass Of Structure = %Lf kg\n",mass_of_structure);
		long double mass_of_h2_by_mass_of_structure = fixed_mass_of_h2 / mass_of_structure;
		mass_of_h2_by_mass_of_structure *= 100.0;
		printf("Mass Of H2 / Mass Of Structure = %Lf \n",mass_of_h2_by_mass_of_structure);


		//************************ Information about SC ******************
		printf("************************ Information about SC ******************\v");
		printf("Porosity (Ø)=   1-{V(Structure)/V(Tank)} = %Lf\n",1-(total_volume_of_structure / total_tank_volume));
		printf("Total Number Of Bubbles (N)= n^3-3(n^2)+2n -2 = %Lf",no_of_bubble_in_tank);
		long double no_of_unit_cubes_on_a_given_edge = 1790.0;
		printf("Where, n = Number Of Unit Cubes On A Given Edge = %Lf\n",no_of_unit_cubes_on_a_given_edge); // obtained from online solution
		long double total_no_of_unit_cube_in_tank_volume;
		printf("Thus, Total No. Of Unit Cubes In Tank Volume = %Lf\n",total_no_of_unit_cube_in_tank_volume = powl(no_of_unit_cubes_on_a_given_edge,3.0));
		long double diffrence_in_unit_cubes;
		printf("Difference In Unit Cubes = %Lf\n",diffrence_in_unit_cubes = powl(no_of_unit_cubes_on_a_given_edge,3.0) - Sc.no_of_unit_cell_required_to_design_the_tank);
		long double p_diff_value_ob;
		printf("Percentage Differnce w.r.t Value Obtained From Volume Calculation = %Lf\n",p_diff_value_ob = diffrence_in_unit_cubes/total_no_of_unit_cube_in_tank_volume);
		printf("\t%Lf\n",p_diff_value_ob*100.0);
	}

	if(user_arr == 2) // for BCC
	{
		struct return_of_Bc Bc = BC_arrangement(r,no_of_bubble_in_tank);
		long double total_tank_volume = Bc.tank_volume_requared;// structure + bubble
		printf("Total Tank Volume = %Lf m^3\n",total_tank_volume);
		long double total_volume_of_structure;
		total_volume_of_structure = (total_tank_volume - total_bubble_volume);
		printf("Total Volume Of Structure = %.15Lf m^3\n",total_volume_of_structure);
		long double mass_of_structure = density_of_structure_material * total_volume_of_structure;
		printf("Mass Of Structure = %Lf kg\n",mass_of_structure);
		long double mass_of_h2_by_mass_of_structure = fixed_mass_of_h2 / mass_of_structure;
		mass_of_h2_by_mass_of_structure *= 100.0;
		printf("Mass Of H2 / Mass Of Structure = %Lf \n",mass_of_h2_by_mass_of_structure);



		//************************ Information about BCC ******************
		printf("************************ Information about BCC ******************\v");
		printf("Porosity (Ø)=   1-{V(Structure)/V(Tank)} = %Lf\n",1-(total_volume_of_structure / total_tank_volume));
		printf("Total Number Of Bubbles (N)= n^3-3(n^2)+2n -2 = %Lf",no_of_bubble_in_tank);
		long double no_of_unit_cubes_on_a_given_edge = 1420.0;
		printf("Where, n = Number Of Unit Cubes On A Given Edge = %Lf\n",no_of_unit_cubes_on_a_given_edge); // obtained from online solution
		long double total_no_of_unit_cube_in_tank_volume;
		printf("Thus, Total No. Of Unit Cubes In Tank Volume = %Lf\n",total_no_of_unit_cube_in_tank_volume = powl(no_of_unit_cubes_on_a_given_edge,3.0));
		long double diffrence_in_unit_cubes;
		printf("Difference In Unit Cubes = %Lf\n",diffrence_in_unit_cubes = powl(no_of_unit_cubes_on_a_given_edge,3.0) - Bc.no_of_unit_cell_required_to_design_the_tank);
		long double p_diff_value_ob;
		printf("Percentage Differnce w.r.t Value Obtained From Volume Calculation = %Lf\n",p_diff_value_ob = diffrence_in_unit_cubes/total_no_of_unit_cube_in_tank_volume);
		printf("\t%Lf\n",p_diff_value_ob*100.0);

	}

	if(user_arr == 3) // for FCC
	{
		struct return_of_Fc Fc = FC_arrangement(r,no_of_bubble_in_tank);
		long double total_tank_volume = Fc.tank_volume_raquared;// structure + bubble
		printf("Total Tank Volume = %Lf m^3\n",total_tank_volume);
		long double total_volume_of_structure;
		total_volume_of_structure = (total_tank_volume - total_bubble_volume);
		printf("Total Volume Of Structure = %.15Lf m^3\n",total_volume_of_structure);
		long double mass_of_structure = density_of_structure_material * total_volume_of_structure;
		printf("Mass Of Structure = %Lf kg\n",mass_of_structure);
		long double mass_of_h2_by_mass_of_structure = fixed_mass_of_h2 / mass_of_structure;
		mass_of_h2_by_mass_of_structure *= 100.0;
		printf("Mass Of H2 / Mass Of Structure = %Lf \n",mass_of_h2_by_mass_of_structure);

		//************************ Information about BCC ******************
		printf("************************ Information about FCC ******************\v");
		printf("Porosity (Ø)=   1-{V(Structure)/V(Tank)} = %Lf\n",1-(total_volume_of_structure / total_tank_volume));
		printf("Total Number Of Bubbles (N)= n^3-3(n^2)+2n -2 = %Lf",no_of_bubble_in_tank);
		long double no_of_unit_cubes_on_a_given_edge = 1128.0;
		printf("Where, n = Number Of Unit Cubes On A Given Edge = %Lf\n",no_of_unit_cubes_on_a_given_edge); // obtained from online solution
		long double total_no_of_unit_cube_in_tank_volume;
		printf("Thus, Total No. Of Unit Cubes In Tank Volume = %Lf\n",total_no_of_unit_cube_in_tank_volume = powl(no_of_unit_cubes_on_a_given_edge,3.0));
		long double diffrence_in_unit_cubes;
		printf("Difference In Unit Cubes = %Lf\n",diffrence_in_unit_cubes = powl(no_of_unit_cubes_on_a_given_edge,3.0) - Fc.no_of_unit_cell_required_to_design_the_tank);
		long double p_diff_value_ob;
		printf("Percentage Differnce w.r.t Value Obtained From Volume Calculation = %Lf\n",p_diff_value_ob = diffrence_in_unit_cubes/total_no_of_unit_cube_in_tank_volume);
		printf("\t%Lf\n",p_diff_value_ob*100.0);

	}

}
