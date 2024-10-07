#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// VISHRUTI - THERMODYNAMICS
// #include <stdio.h>
// #include<math.h>
const double b = 2.897771955;
const double IDEAL_GAS_CONSTANT = 8.314;
const double PI = 3.14;
void heat_transfer();
void Thermo();
void classic_thermo();
void ktg();
void calorimetry();

void Thermo()
{   int n;
    printf("\nTHERMODYNAMICS!!");
    printf("\n1.Heat Transfer.");
    printf("\n2.Classic Thermodynamics & Thermal Properties");
    printf("\n3.Ideal Gas Laws");
    printf("\n4.Calorimetry");
    printf("\nEnter Chapter Number :");
    scanf("%d",&n);
    if (n==1)
    heat_transfer();
    else if(n==2)
    classic_thermo();
    else if (n==3)
    ktg();
    else if(n==4)
    calorimetry();
    return ;
}
void heat_transfer()
{   int sub_topic_1;
     printf("\n1.Conduction");
    printf("\n2.Thermal Resistance");
    printf("\n3.Emmisivity");
    printf("\n4.Laws.");
    printf("\nEnter Topic Number:");
    scanf("%d",&sub_topic_1);
    if(sub_topic_1==1)
    {
        double thermal_conductivity, cross_sectional_area, temperature_difference, thickness;


        printf("Enter thermal conductivity (k): ");
        scanf("%lf", &thermal_conductivity);


        printf("Enter cross-sectional area (A): ");
        scanf("%lf", &cross_sectional_area);


        printf("Enter temperature difference (dT): ");
        scanf("%lf", &temperature_difference);


        printf("Enter thickness (d): ");
        scanf("%lf", &thickness);


        double heat_flow_rate = (thermal_conductivity * cross_sectional_area * temperature_difference) / thickness;


        printf("Rate of heat flow (Q): %lf\n", heat_flow_rate);


    }
     else if(sub_topic_1==2)
    {
   
        double thickness, thermal_conductivity, cross_sectional_area;


        printf("Enter thickness (d): ");
        scanf("%lf", &thickness);


        printf("Enter thermal conductivity (k): ");
        scanf("%lf", &thermal_conductivity);


        printf("Enter cross-sectional area (A): ");
        scanf("%lf", &cross_sectional_area);


        double thermal_resistance = thickness / (thermal_conductivity * cross_sectional_area);


        printf("Thermal resistance (R): %lf\n", thermal_resistance);


    }
   
     else if(sub_topic_1==3)
      {
        double emitted_radiation, radiation_black_body;
        printf("Enter emitted radiation: ");
        scanf("%lf", &emitted_radiation);
        printf("Enter radiation emitted by a black body at the same temperature: ");
        scanf("%lf", &radiation_black_body);
        double emissivity = emitted_radiation / radiation_black_body;
        printf("Emissivity : %lf\n", emissivity);
      }
      else if (sub_topic_1==4)
      {
         double temperature;
         printf("Enter the absolute temperature of the black body (in Kelvin): ");
         scanf("%lf", &temperature);
        if (temperature <= 0)
        {
        printf("Invalid temperature. Please enter a positive temperature.\n");
        }
        double peak_wavelength = b / temperature;
         printf("Peak wavelength : %lf meters\n", peak_wavelength);


    }
     
    }


void classic_thermo()
{
    int sub_topic_2;
    printf("\n1.Temperature Conversion");
    printf("\n2.First Law of Thermodynamics(internal energy)");
    printf("\n3.Work Done by different process.");
    printf("\n4.Enthalpy");
    printf("\n5.Entropy");
    printf("\n6.Gibbs free energy");
    printf("\n7.Efficiency of cycles");
    printf("\nEnter the topic number:");
    scanf("%d",&sub_topic_2);
    if (sub_topic_2==2)
    {   double heat,work;
        printf("\nEnter the heat added to the system:");
        scanf("%lf",&heat);
        printf("\nWork done by system:");
        scanf("%lf",&work);
        double U=heat-work;
        printf("\nThe Internal energy is %lf",U);
    }






    if (sub_topic_2==3)
    {
      int sub_topic;  
      printf("\nChoose Process:");
      printf("\n1.Isothermal Irreversible");
      printf("\n2.Isothermal Reverssible");
      printf("\n3.Free Expansion");
      printf("\n4.Adiabatic process");
      printf("\n5.Isochoric Process");
      printf("\n6.Isobaric process");
      scanf("%d",sub_topic);
      if (sub_topic==1)
      {
        double external_pressure, volume_change;
        printf("Enter the external pressure (P_ext): ");
        scanf("%lf", &external_pressure);
        printf("Enter the change in volume (Delta V): ");
        scanf("%lf", &volume_change);
        double work_done = -external_pressure * volume_change;
        printf("Work done (W): %lf\n", work_done);
      }
      else if(sub_topic==3||sub_topic==5)
      {
        printf("\nWork Done=0");
      }
      else if (sub_topic==6)
      {
        double perssure,volume_in,volume_f;
        printf("\nConstant Pressure:");
        scanf("%lf",perssure);
        printf("\nInitial volume:");
        scanf("%lf",volume_in);
         printf("\nFinal volume:");
        scanf("%lf",volume_f);
        double volume_diff= volume_f-volume_in;
        double work_done= perssure*volume_diff;
        printf("\nWork Done=%lf",work_done);
      }
      else if (sub_topic==2)
      {
         
            double n, R, T, Vi, Vf;


       
            printf("Enter the number of moles of gas (n): ");
            scanf("%lf", &n);


            printf("Enter the absolute temperature (T) in Kelvin: ");
            scanf("%lf", &T);


            printf("Enter the initial volume (Vi): ");
            scanf("%lf", &Vi);


            printf("Enter the final volume (Vf): ");
            scanf("%lf", &Vf);


       
            R = IDEAL_GAS_CONSTANT;
            double work_done = -n * R * T * log(Vf / Vi);




            printf("Work done (W): %lf J\n", work_done);
      }
      else if (sub_topic==4)
      {
           
            double moles, Cv, delta_T, gamma;


            printf("Enter the number of moles of gas (n): ");
            scanf("%lf", &moles);


            printf("Enter the molar specific heat at constant volume (Cv): ");
            scanf("%lf", &Cv);


            printf("Enter the change in temperature (Delta T): ");
            scanf("%lf", &delta_T);
            printf("\nGamma:");
            scanf("%lf",&gamma);
            double work_done = -(moles * IDEAL_GAS_CONSTANT* delta_T) / (1.0 - gamma);
            printf("Work done (W): %lf J\n", work_done);


      }
     
     }
    if (sub_topic_2==4)
     {
            double internal_en,p,v,enthalpy;
            printf("\nINTERNAL ENERGY(U):");
            scanf("%lf",&internal_en);
            printf("\nPRESSURE:");
            scanf("%lf",&p);
            printf("\nVOLUME:");
            scanf("%lf",&v);
            enthalpy=(internal_en+(p*v));
            printf("\nENTHALPY=%lf",enthalpy);
     }
     if (sub_topic_2==5)
     {
            double heat, temp,entropy;
            printf("\nHat in rev process:");
            scanf("%lf",&heat);
            printf("\nTEMPERATURE:");
            scanf("%lf",&temp);
            entropy=heat/temp;
            printf("\nENTROPY=%lf",entropy);
     }
     if (sub_topic_2==6)
     {
         double enthalpy, temperature, entropy;
         printf("Enter the enthalpy (H): ");
        scanf("%lf", &enthalpy);


        printf("Enter the absolute temperature (T): ");
        scanf("%lf", &temperature);


        printf("Enter the entropy (S): ");
        scanf("%lf", &entropy);
        double gibbs_free_energy = enthalpy - (temperature * entropy);
        printf("Gibbs free energy (G): %lf J\n", gibbs_free_energy);


     }
     
     if (sub_topic_2==7)
     {
        double temperature_cold, temperature_hot;
        printf("Enter the absolute temperature of the cold reservoir (T_C): ");
        scanf("%lf", &temperature_cold);


        printf("Enter the absolute temperature of the hot reservoir (T_H): ");
        scanf("%lf", &temperature_hot);


       
        if (temperature_cold <= 0 || temperature_hot <= 0)
        {
            printf("Invalid temperatures. Please enter positive values for absolute temperatures.\n");
             
        }


   
        double efficiency = 1 - (temperature_cold / temperature_hot);
        printf("Carnot cycle efficiency (eta): %lf\n", efficiency);


     }
     if (sub_topic_2==1)
     {  float celsius,fahrenheit;
        printf("\n Enter the Temparature in Celsius : ");
        scanf("%f",&celsius);
        fahrenheit = (1.8 * celsius) + 32;
        printf("\n Temperature in Fahrenheit : %f ",fahrenheit);
        }
}
void ktg()
{   int sub_topic_3;
    printf("\n1.DIFFERENT SPEED");
    printf("\n2.KINETIC ENERGY CALCULATION");
    printf("\n3.SPECFIC HEAT CAPACITIES");
    printf("\n4.GAMMA VALUES");
    printf("\n5.MEAN FREE PATH");
    printf("\n6.PRESSURE EXTERTED BY GAS");
    scanf("%d",&sub_topic_3);
        if (sub_topic_3==1)
        {   int sub_topic_4;
            double T,M;
            printf("\nChoose the speed");
            printf("\n1.ROOT MEAN SQUARE");
            printf("\n2.AVERAGE SPEED");
            printf("\n3.MOST PROBABLE SPEED");
            scanf("%d",&sub_topic_4);
            if (sub_topic_4==1)
            {
                printf("\nTemperature:");
                scanf("%lf",&T);
                printf("\nMass");
                scanf("%lf",&M);
                double vrms=((3*IDEAL_GAS_CONSTANT*T)/(M));
                double Vrms=sqrt(vrms);
                printf("\nVrms=%lf",Vrms);
            }
            else if (sub_topic_4==2)
            {
                printf("\nTemperature:");
                scanf("%lf",&T);
                printf("\nMass");
                scanf("%lf",&M);
                double vav=((8*IDEAL_GAS_CONSTANT*T)/(PI*M));
                double Vav=sqrt(vav);
                printf("\nVrms=%lf",Vav);
            }
            else if (sub_topic_4==3)
            {
                printf("\nTemperature:");
                scanf("%lf",&T);
                printf("\nMass");
                scanf("%lf",&M);
                double vmp=((2*IDEAL_GAS_CONSTANT*T)/(M));
                double Vmp=sqrt(vmp);
                printf("\nVrms=%lf",Vmp);
            }
           
        }
       
        else if (sub_topic_3==2)
        {   int sub_topic_5;
            double T;
            printf("\nENTER TEMPERATURE:");
            scanf("%lf",&T);
            double KE= 1.5*IDEAL_GAS_CONSTANT*T;
            printf("KE=%lf",KE);


        }
        else if (sub_topic_3==3)
        {
            float DOF;
            float Cp,Cv;
            printf("\nEnter the Degree of Freedom of gas");
            scanf("%f",&DOF);
            Cp=(1+(DOF/2))*IDEAL_GAS_CONSTANT;
            Cv=Cp-IDEAL_GAS_CONSTANT;
            printf("\nValue of Cp :%f",Cp);
            printf("\nValue of Cv=%f",Cv);
\
        }
        else if (sub_topic_3==4)
        {  char RESPONSE;
           printf("\nMONOATOMIC GAS=5/3");
           printf("\nDIATOMIC GAS=7/5");
           printf("\nDo you want to calculate for POLYATOMIC GAS");
           fflush(stdin);
           scanf("%c",&RESPONSE);
           fflush(stdin);
           if (RESPONSE=='Y' || RESPONSE=='y')
           {    float DOF;
                printf("\nENTER DOF:");
                scanf("%f",&DOF);
                float gamma=((4+DOF)/(3+DOF));
                printf("The value of Gamma is: %f",gamma);
           }
           
        }
        else if (sub_topic_3==5)
        {
            float diameter,molecule_density;
            float mean_free;
            printf("\nEnter the diameter of molecule:");
            scanf("%f",&diameter);
            printf("\nMolecule density:");
            scanf("%f",&molecule_density);
            mean_free=(1/(1.4142*PI*molecule_density*diameter*diameter));
            printf("\nMEAN FREE PATH %f",mean_free);
        }
        else if (sub_topic_3==6)
        {
            double density,Vrms;
            double pressure;
            printf("\nDensity of molecule");
            scanf("%lf",&density);
            printf("\nVrms=");
            scanf("%lf",&Vrms);
            pressure=0.33333*density*Vrms*Vrms;
            printf("\nPressure exerted=%lf",pressure);
        }
       
       


    }
   
void calorimetry()
{   double M,S,dt,l;
    int formula;
    printf("\nChoose the formula");
    printf("\n1.Q=MSdt");
    printf("\n2.Q=ml");
    if (formula==1)
    {
            printf("\nEnter MASS(M):");
            scanf("%lf",&M);
            printf("\nSPECIFIC HEAT CAPACITY:");
            scanf("lf",&S);
            printf("\nCHANGE IN TEMPERATURE:");
            scanf("%lf",&dt);
            double heat_form1=M*S*dt;
            printf("\nHEAT =%lf",heat_form1);


    }
    else if (formula==2)
    {
            printf("\nENTER MASS");
            scanf("%lf",&M);
            printf("\nLATENT HEAT:");
            scanf("%lf",&l);
            double heat_form2=M*l;
            printf("\nHEAT:%lf",heat_form2);


    }
}






// Mithilesh - Gravitation:
// #include<stdio.h>
// #include<math.h>

void Gravitation()
{
    double F, m1, m2, x, G = 6.67408 * pow(10, -11);

    printf("Give the mass of first body(in KG):");
    scanf("%lf", &m1);
    printf("Give the mass of second body(in KG):");
    scanf("%lf", &m2);
    printf("Give the distance between centre of masses of two body(in meters):");
    scanf("%lf", &x);

    F = G * m1 * m2 / pow(x, 2);

    if (F < pow(10, -15))
        printf("\nThe Gravitation force between two given bodies is %.20lf Newton.", F);
    else if (F < pow(10, -10))
        printf("\nThe Gravitation force between two given bodies is %.15lf Newton.", F);
    else if (F < pow(10, -5))
        printf("\nThe Gravitation force between two given bodies is %.10lf Newton.", F);
    else if (F < pow(10, 0))
        printf("\nThe Gravitation force between two given bodies is %.5lf Newton.", F);
    else if (F < pow(10, 3))
        printf("\nThe Gravitation force between two given bodies is %.2lf Newton.", F);
    else if (F < pow(10, 6))
        printf("\nThe Gravitation force between two given bodies is %.2lfx10^3 Newton.", F / pow(10, 3));
    else if (F < pow(10, 9))
        printf("\nThe Gravitation force between two given bodies is %.2lfx10^6 Newton.", F / pow(10, 6));
    else
        printf("\nThe Gravitation force between two given bodies is %.0lfx10^6 Newton.", F / pow(10, 6));
}

void change_in_g()
{
    double G = 6.674 * pow(10, -11), M = 5.97219 * pow(10, 24), R = 6371 * pow(10, 3);
    double g = G * M / pow(R, 2), g_new;
    char choice;
here:
    printf("This program show you g at certian height and under certian depth and %% change in g from earth's surface.\nfor height:\th\nfor depth:\td");

    printf("\nGive your choice: ");
    scanf("%c", &choice);

    fflush(stdin);

    if (choice == 'h' || choice == 'H')
    {
        // for height
        double h; // in km

        printf("At what height you want g from eath's surface(in km): ");
        scanf("%lf", &h);

        g_new = g / pow((1 + 1000 * h / R), 2);
        printf("g = %.1lf ms^(-2)", g_new);
    }

    else if (choice == 'd' || choice == 'D')
    {
        // for depth
        double d; // in km

        printf("At what height you want g from eath's surface(in km): ");
        scanf("%lf", &d);

        g_new = g * (1 - 1000 * d / R);
        printf("g = %.1lf ms^(-2)", g_new);
    }

    else
    {
        printf("You give some wrong input.");
        goto here;
    }

    printf("\ng changed %lf %%", 100 * (g - g_new) / g);
}

// Center of Mass:
// #include <stdio.h>
// #include <stdlib.h>

void COM()
{
    int dim;
    int n;
    double *mass;
    double sum_mass = 0;

here:
    printf("In which dimesion do you want(1D/2D/3D):");
    scanf("%d", &dim);

    fflush(stdin);

    if (dim == 1)
    {
        double sum_x = 0;
        double *x;
        double com_x;

        printf("How many masses do you want:");
        scanf("%d", &n);

        mass = (double *)calloc(n, sizeof(double));
        x = (double *)calloc(n, sizeof(double));

        for (int i = 0; i < n; i++)
        {
            printf("\nGive mass of m%d (in kg):", i + 1);
            scanf("%lf", &mass[i]);

            printf("Give the coordinates of m%d\n", i + 1);
            printf("x%d(in m):", i + 1);
            scanf("%lf", &x[i]);

            sum_mass += mass[i];
            sum_x += x[i] * mass[i];
        }

        com_x = sum_x / sum_mass;

        printf("\nThe centre of mass of your system is:\n %lf m\n", com_x);
    }
    else if (dim == 2)
    {
        double sum_x = 0;
        double sum_y = 0;
        double *x;
        double *y;
        double com_x;
        double com_y;

        printf("How many masses do you want:");
        scanf("%d", &n);

        mass = (double *)calloc(n, sizeof(double));
        x = (double *)calloc(n, sizeof(double));
        y = (double *)calloc(n, sizeof(double));

        for (int i = 0; i < n; i++)
        {
            printf("\nGive mass of m%d (in kg):", i + 1);
            scanf("%lf", &mass[i]);

            printf("Give the coordinates of m%d\n", i + 1);
            printf("x%d(in m):", i + 1);
            scanf("%lf", &x[i]);
            printf("y%d(in m):", i + 1);
            scanf("%lf", &y[i]);

            sum_mass += mass[i];
            sum_x += x[i] * mass[i];
            sum_y += y[i] * mass[i];
        }

        com_x = sum_x / sum_mass;
        com_y = sum_y / sum_mass;

        printf("\nThe centre of mass of your system is:\n (%lf m, %lf m)\n", com_x, com_y);
    }
    else if (dim == 3)
    {
        double sum_x = 0;
        double sum_y = 0;
        double sum_z = 0;
        double *x;
        double *y;
        double *z;
        double com_x;
        double com_y;
        double com_z;

        printf("How many masses do you want:");
        scanf("%d", &n);

        mass = (double *)calloc(n, sizeof(double));
        x = (double *)calloc(n, sizeof(double));
        y = (double *)calloc(n, sizeof(double));
        z = (double *)calloc(n, sizeof(double));

        for (int i = 0; i < n; i++)
        {
            printf("\nGive mass of m%d (in kg):", i + 1);
            scanf("%lf", &mass[i]);

            printf("Give the coordinates of m%d\n", i + 1);
            printf("x%d(in m):", i + 1);
            scanf("%lf", &x[i]);
            printf("y%d(in m):", i + 1);
            scanf("%lf", &y[i]);
            printf("z%d(in m):", i + 1);
            scanf("%lf", &z[i]);

            sum_mass += mass[i];
            sum_x += x[i] * mass[i];
            sum_y += y[i] * mass[i];
            sum_z += z[i] * mass[i];
        }

        com_x = sum_x / sum_mass;
        com_y = sum_y / sum_mass;
        com_z = sum_z / sum_mass;

        printf("\nThe centre of mass of your system is:\n (%lf m, %lf m, %lf m)\n", com_x, com_y, com_z);
    }
    else
    {
        printf("You given something wrong input.\n");
        goto here;
    }
}

void COM_formula()
{
    /*
    half disc:4r/3pi
    annular half disc:4(r2^3-r1^3)/{3pi(r2^2-r1^2)}
    half ring:2r/pi
    triangular plate:h/3
    hollow hemisphere:r/2
    solid hemisphere:3r/8
    solid cone:h/4
    hollow cone:h/3
    */
    printf("\nR-radius\th-height\n");
    printf("\nHalf disc:\t\t (0,4R/3pi)[from centre]\nHalf ring:\t\t (0,2R/pi)[from centre]\nHollow hemisphere:\t (0,3R/8)[from centre]\nSolid hemisphere:\t (0,3R/8)[from centre]\nHollow cone:\t\t (0,h/3)[from centre of base]\nSolid cone:\t\t (0,h/4)[form centre of base]\nTriangualr plate:\t y=h/3[from base]");
    printf("\n\nr2-outer radius\tr1-inner radius");
    printf("\n\nAnnular half disc:\t(0,4(r2^3-r1^3)/{3pi(r2^2-r1^2)})[from centre]");
}

// Soham- WORK ,ENERGY,POWER:

// #include<stdio.h>
// #include<math.h>
// #include<stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

float work()
{
    float a;
    char c, b;
    float f, d;

    printf("give the value of force with unit ( n for newton,d for dyne)\n");
    scanf("%f%c", &f, &c);

    if (c == 'd')
    {
        f = f / pow(10, 5);
    }
    else
    {
    }

    printf("give the value of displacement with unit ( m for meter,k for kilometer,c for centimeter)\n");
    scanf("%f%c", &d, &b);

    if (b == 'k')
    {
        d = d * 1000;
    }
    else if (b == 'c')
    {
        d = d / 100;
    }
    else
    {
    }

    printf("give the value of angle in degree\n");
    scanf("%f", &a);

    a = a * (M_PI / 180.0);
    printf("work done is :%.10f joules\n", f * d * (cos(a)));
    return f * d * (cos(a));
}

void k_energy()
{
    char c, b;
    float m, v;

    printf("Give the value of mass with unit (g for gram, k for kilogram):\n");
    scanf("%f %c", &m, &c);

    if (c == 'g')
    {
        m = m / 1000;
    }

    printf("Give the value of velocity with unit (m for meter per sec, k for kilometer per sec, c for centimeter per sec):\n");
    scanf("%f %c", &v, &b);

    if (b == 'k')
    {
        v = v * 1000;
    }
    else if (b == 'c')
    {
        v = v / 100;
    }

    printf("kinetic energy is: %.9f joules\n", (1.0 / 2.0) * m * pow(v, 2));
}

void p_energy()
{
    char c, b;
    float m, h;
    float g = 9.8;

    printf("Give the value of mass with unit (g for gram, k for kilogram):\n");
    scanf("%f %c", &m, &c);

    if (c == 'g')
    {
        m = m / 1000;
    }

    printf("Give the value of height with unit (m for meter, k for kilometer, c for centimeter):\n");
    scanf("%f %c", &h, &b);

    if (b == 'k')
    {
        h = h * 1000;
    }
    else if (b == 'c')
    {
        h = h / 100;
    }

    printf("potential energy is: %.9f joules\n", m * g * h);
}

void power()
{
    float t;

    printf("give the value of time required to do work in sec\n");
    scanf("%f", &t);
    printf("power require to do this work is :%.10f", work() / t);
}

// HEAT CAPACITY AND HEAT TRANSFER

// #include<stdio.h>

float heat(char a)
{
    char b;
    float m, c, t;

    printf("Give the value of mass with unit (g for gram, k for kilogram):\n");
    scanf("%f %c", &m, &b);
    if (b == 'g')
    {
        m = m / 1000;
    }

    printf("Give the value of specific heat capacity:\n");
    scanf("%f", &c);

    if (a == 'c')
    {
        printf("%f", c * m);
        return 0;
    }

    printf("Give the value of change in temp:\n");
    scanf("%f", &t);

    printf("heat transfer is: %f joules\n", m * c * t);
    return 0;
}

// Shubham- one d,two d,nlm
// #include<stdio.h>
// #include<math.h>


void one_D()
{

    float initialPosition, velocity, time;

    // Input initial position, velocity, and time
    printf("Enter initial position: ");
    scanf("%f", &initialPosition);

    printf("Enter velocity: ");
    scanf("%f", &velocity);

    printf("Enter time: ");
    scanf("%f", &time);

    // Calculate final position using the motion equation: final_position = initial_position + (velocity * time)
    float finalPosition = initialPosition + (velocity * time);

    // Display the result
    printf("Final position after %.2f seconds: %.2f\n", time, finalPosition);
}
void two_D()
{
    float initialX, initialY, velocityX, velocityY, time;

    // Input initial position, velocity components, and time
    printf("Enter initial X position: ");
    scanf("%f", &initialX);

    printf("Enter initial Y position: ");
    scanf("%f", &initialY);

    printf("Enter X component of velocity: ");
    scanf("%f", &velocityX);

    printf("Enter Y component of velocity: ");
    scanf("%f", &velocityY);

    printf("Enter time: ");
    scanf("%f", &time);

    // Calculate final position using the motion equations: final_position_x = initial_position_x + (velocity_x * time)
    //                                                      final_position_y = initial_position_y + (velocity_y * time)
    float finalX = initialX + (velocityX * time);
    float finalY = initialY + (velocityY * time);

    // Display the result
    printf("Final position after %.2f seconds: X=%.2f, Y=%.2f\n", time, finalX, finalY);
}

void NLM()
{
    // Variables
    double mass, acceleration, force;

    // Input mass and acceleration from the user
    printf("Enter the mass of the object (in kg): ");
    scanf("%lf", &mass);

    printf("Enter the acceleration of the object (in m/s^2): ");
    scanf("%lf", &acceleration);

    // Calculate force using Newton's second law: F = m * a
    force = mass * acceleration;

    // Output the result
    printf("The force acting on the object is: %.2lf Newtons\n", force);
}

void shubham()
{

    char choice;

    printf("Give topic name(1D/2D/force): ");
    scanf("%c", &choice);

    fflush(stdin);

    if (choice == '1')
        one_D();
    else if (choice == '2')
        two_D();
    else if (choice == 'f')
        NLM();
}
// SHIVAM
// TOPIC NAME : SOLID AND FLUID MECHANICS

// #include <stdio.h>

// Function to calculate stress
float calculateStress(float force, float area)
{
    return force / area;
}

// Function to calculate strain
float calculateStrain(float changeInConfiguration, float originalConfiguration)
{
    return changeInConfiguration / originalConfiguration;
}

float calculateLongitudinalStrain(float changeInLength, float originalLength)
{
    return changeInLength / originalLength;
}

float calculatemodulasOfElasticity(float stress, float strain)
{
    return stress / strain;
}

float calculateVolumeStrain(float changeInVolume, float originalVolume)
{
    return changeInVolume / originalVolume;
}

float calculateYoungModulus(float normalStress, float longitudnalStrain)
{
    return normalStress / longitudnalStrain;
}

float calculateBulkModulus(float force, float area, float changeInVolume, float totalVolume)
{
    return -(force / area) / (changeInVolume / totalVolume);
}

float calculateCompressibility(float force, float area, float changeInVolume, float totalVolume)
{
    return (changeInVolume / totalVolume) / (force / area);
}

float calculateEnergyPerUnitVolume(float force, float area, float changeInLength, float originalLength)
{
    return 0.5 * (force * changeInLength) / (area * originalLength);
}

void shivam()
{
    int choice;
    float force, area, changeInLength, originalLength, changeInConfiguration, originalConfiguration, stress, strain, changeInVolume, originalVolume, normalStress, longitudnalStrain, totalVolume;

    printf("1. Calculate Stress\n");
    printf("2. Calculate Strain\n");
    printf("3. Calculate longitudinal Strain\n");
    printf("4. Calculate Modulus of Elasticity\n");
    printf("5. Calculate volume strain\n");
    printf("6. Calculate Young Modulus\n");
    printf("7. Calculate Bulk Modulus\n");
    printf("8. Calculate Compressibility\n");
    printf("9. Calculate Energy Per Unit Volume\n");
    printf("\n\n");
    printf("Select an option: ");
    scanf("%d", &choice);

    switch (choice)
    {
    case 1:
        printf("Enter force applied (in Newtons): ");
        scanf("%f", &force);
        printf("Enter cross-sectional area (in square meters): ");
        scanf("%f", &area);
        printf("Stress: %f N/m^2\n", calculateStress(force, area));
        break;

    case 2:
        printf("Enter change In Configuration(in meters): ");
        scanf("%f", &changeInConfiguration);
        printf("Enter original Configuration (in meters): ");
        scanf("%f", &originalConfiguration);
        printf("Strain: %f\n", calculateStrain(changeInConfiguration, originalConfiguration));
        break;

    case 3:
        printf("Enter change in length (in meters): ");
        scanf("%f", &changeInLength);
        printf("Enter original length (in meters): ");
        scanf("%f", &originalLength);
        printf("longitudinal strain : %f N/m^2\n", calculateLongitudinalStrain(changeInLength, originalLength));
        break;

    case 4:
        printf("Enter stress : ");
        scanf("%f", &stress);
        printf("Enter strain : ");
        scanf("%f", &strain);
        printf("Modulus of elasticity: %f\n", calculatemodulasOfElasticity(stress, strain));
        break;

    case 5:
        printf("Enter change In Volume : ");
        scanf("%f", &changeInVolume);
        printf("Enter original volume : ");
        scanf("%f", &originalVolume);
        printf("Volume strain : %f\n", calculateVolumeStrain(changeInVolume, originalVolume));
        break;

    case 6:
        printf("Enter normal Stress: ");
        scanf("%f", &normalStress);
        printf("Enter longitudnalStrain  : ");
        scanf("%f", &longitudnalStrain);
        printf("Young Modulus :%f", calculateYoungModulus(normalStress, longitudnalStrain));
        break;

    case 7:
        printf("Enter force : ");
        scanf("%f", &force);
        printf("Enter area: ");
        scanf("%f", &area);
        printf("Enter change in volume : ");
        scanf("%f", &changeInVolume);
        printf("Enter total volume : ");
        scanf("%f", &totalVolume);
        printf("bulk modulus : %f\n", calculateBulkModulus(force, area, changeInVolume, totalVolume));
        break;

    case 8:
        printf("Enter force : ");
        scanf("%f", &force);
        printf("Enter area : ");
        scanf("%f", &area);
        printf("Enter change in volume : ");
        scanf("%f", &changeInVolume);
        printf("Enter total volume: ");
        scanf("%f", &totalVolume);
        printf("compressibility: %f \n", calculateCompressibility(force, area, changeInVolume, totalVolume));
        break;

    case 9:
        printf("Enter force : ");
        scanf("%f", &force);
        printf("Enter area : ");
        scanf("%f", &area);
        printf("Enter change In Length : ");
        scanf("%f", &changeInLength);
        printf("Enter original Length: ");
        scanf("%f", &originalLength);
        printf("Energy Per Unit Volume: %f \n", calculateEnergyPerUnitVolume(force, area, changeInLength, originalLength));
        break;

    default:
        printf("Invalid choice\n");
    }

    
}

// 2.

// #include <stdio.h>
// #include <math.h>

// Function to calculate pressure
float calculatePressure(float force, float area)
{
    return (force / area);
}

// Function to calculate flow rate
float calculateFlowRate(float area, float velocity)
{
    return (area * velocity);
}

// Function to calculate velocity
float calculateVelocity(float flowRate, float area)
{
    return (flowRate / area);
}

float calculateContinuity(float area1, float velocity1, float area2, float velocity2)
{
    velocity2 = (area1 * velocity1) / area2;
    return (velocity2);
}

float calculateBernoulli(float pressure, float density, float velocity, float height, float constant)
{
    constant = pressure + 0.5 * density * velocity * velocity + density * 9.81 * height;
    return constant;
}

float calculateResistance(float viscosity, float length, float radius, float resistance)
{
    resistance = (56 * viscosity * length) / (22 * pow(radius, 4));
    return resistance;
}

void shivam2()
{
    int choice;
    float force, area, velocity, flowRate, area1, velocity1, area2, velocity2, pressure, density, height, constant, viscosity, length, radius, resistance;

    // Get user input
    printf("Choose a calculation:\n");
    printf("1. Pressure\n2. Velocity\n3. Flow Rate\n4. velocity2\n5. bernolli constant\n6. Resistance\n");
    printf("Enter your choice (1 to 6): ");
    scanf("%d", &choice);

    // Perform calculation based on user choice
    switch (choice)
    {
    case 1:
        printf("Enter force (N): ");
        scanf("%f", &force);

        printf("Enter area (m^2): ");
        scanf("%f", &area);

        // Calculate and display pressure
        printf("Pressure: %.2f Pa\n", calculatePressure(force, area));
        break;

    case 2:
        printf("Enter flow rate (m^3/s): ");
        scanf("%f", &flowRate);

        printf("Enter area (m^2): ");
        scanf("%f", &area);

        // Calculate and display velocity
        printf("Velocity: %.2f m/s\n", calculateVelocity(flowRate, area));
        break;

    case 3:
        printf("Enter area (m^2): ");
        scanf("%f", &area);

        printf("Enter velocity (m/s): ");
        scanf("%f", &velocity);

        // Calculate and display flow rate
        printf("Flow Rate: %.2f m^3/s\n", calculateFlowRate(area, velocity));
        break;

    case 4:
        printf("Enter area at section 1 (m^2): ");
        scanf("%f", &area1);

        printf("Enter velocity at section 1 (m/s): ");
        scanf("%f", &velocity1);

        printf("Enter area at section 2 (m^2): ");
        scanf("%f", &area2);

        printf("value of velocity2 : %.2f\n", calculateContinuity(area1, velocity1, area2, velocity2));
        break;

    case 5:
        printf("Enter pressure (Pa): ");
        scanf("%f", &pressure);

        printf("Enter density (kg/m^3): ");
        scanf("%f", &density);

        printf("Enter velocity (m/s): ");
        scanf("%f", &velocity);

        printf("Enter height (m): ");
        scanf("%f", &height);

        printf("Bernoulli constant : %.2f\n", calculateBernoulli(pressure, density, velocity, height, constant));
        break;

    case 6:
        printf("Enter viscosity of the fluid (Pa.s): ");
        scanf("%f", &viscosity);

        printf("Enter length of the pipe (m): ");
        scanf("%f", &length);

        printf("Enter radius of the pipe (m): ");
        scanf("%f", &radius);

        printf("resistance : %.4f\n", calculateResistance(viscosity, length, radius, resistance));
        break;

    default:
        printf("Invalid choice\n");
        
    }
}

// Mayank - Rotational Motion 🎉

// #include<stdio.h>
// #include<math.h>

float hollow_ring(float mass, float rad1)
{
    printf("Enter the mass and radius : ");
    scanf("%f %f", &mass, &rad1);
    float i = mass * rad1 * rad1;
    return i;
}

float thick_ring(float mass, float rad1, float rad2)
{
    printf("Enter the mass, inner radius and outer radius : ");
    scanf("%f %f %f", &mass, &rad1, &rad2);
    float i = (mass * ((rad1 * rad1) + (rad2 * rad2))) / 2;
    return i;
}

float spherical_shell(float mass, float rad1, float rad2)
{
    printf("Enter the mass, outer radius and inner radius : ");
    scanf("%f %f %f", &mass, &rad2, &rad1);
    float diff_5 = 0, diff_3 = 0;
    diff_5 = pow(rad2, 5) - pow(rad1, 5); // numerator difference
    diff_3 = pow(rad2, 3) - pow(rad1, 3); // denominator difference
    float i = mass * diff_5 / diff_3;
    return i;
}

float uniform_plate(float mass, float rad1)
{
    printf("Enter the mass and length : ");
    scanf("%f %f", &mass, &rad1);
    float i = mass * rad1 * rad1;
    return i;
}
float uniform_box(float mass, float rad1, float rad2)
{
    printf("Enter the mass, length and width : ");
    scanf("%f %f %f", &mass, &rad1, &rad2);
    float i = mass * ((rad1 * rad1) + (rad2 * rad2));
    return i;
}

void MoI()
{
    int n = 0;
    float mass, rad1, rad2, ans;
    printf("\n1.Thin ring/ hollow cylinder(Central)\n2.Thin ring(Diameter)\n3.Thick walled ring(Central)\n4.Uniform disc(Central)\n5.Uniform disc(Diameter)");
    printf("\n6.Thin walled hollow sphere(Central)\n7.Solid sphere(Central)\n8.Spherical shell(Central)\n9.Thin rod/ rectangular plate(Perpendicular-centre)");
    printf("\n10.Thin rod/ rectangular plate(Perpendicular-end)\n11.Uniform plate/box(Central)\n12.Solid cone(Central)\n13.Hollow cone(Central)\n\n");

    printf("Enter the number for which u want to find the moment of inertia : "); // here we have to give the code for the object whose calculations are to be done...
    scanf("%d", &n);

    if (n == 1)
    {
        ans = hollow_ring(mass, rad1);
        printf("\nMoment of Inertia = %.4f\n", ans);
        return;
    }
    else if (n == 2 || n == 4 || n == 13)
    {
        ans = hollow_ring(mass, rad1);
        printf("\nMoment of Inertia = %.4f\n", ans / 2);
        return;
    }
    else if (n == 3)
    {
        ans = thick_ring(mass, rad1, rad2);
        printf("\nMoment of Inertia = %.4f\n", ans);
        return;
    }
    else if (n == 5)
    {
        ans = hollow_ring(mass, rad1);
        printf("\nMoment of Inertia = %.4f\n", ans / 4);
        return;
    }
    else if (n == 6)
    {
        ans = hollow_ring(mass, rad1);
        printf("\nMoment of Inertia = %.4f\n", 2 * ans / 3);
        return;
    }
    else if (n == 7)
    {
        ans = hollow_ring(mass, rad1);
        printf("\nMoment of Inertia = %.4f\n", 2 * ans / 5);
        return;
    }
    else if (n == 8)
    {
        ans = spherical_shell(mass, rad1, rad2);
        printf("\nMoment of Inertia = %.4f\n", 2 * ans / 5);
        return;
    }
    else if (n == 9)
    {
        ans = uniform_plate(mass, rad1);
        printf("\nMoment of Inertia = %.4f\n", ans / 12);
        return;
    }
    else if (n == 10)
    {
        ans = uniform_plate(mass, rad1);
        printf("\nMoment of Inertia = %.4f\n", ans / 3);
        return;
    }
    else if (n == 11)
    {
        ans = uniform_box(mass, rad1, rad2);
        printf("\nMoment of Inertia = %.4f\n", ans / 12);
        return;
    }
    else if (n == 12)
    {
        ans = hollow_ring(mass, rad1);
        printf("\nMoment of Inertia = %.4f\n", 3 * ans / 10);
        return;
    }
}

// Mayank - Friction and Kinetics 🎉

// #include<stdio.h>
// #include<math.h>

float force_1(float u, float nf)
{
    printf("Enter the values of static coefficient and normal force : ");
    scanf("%f %f", &u, &nf);
    float i = u * nf;
    return i;
}

float force_2(float u, float mass)
{
    float g = 9.8;
    printf("Enter the values of static coefficient and mass : ");
    scanf("%f %f", &u, &mass);
    float i = u * mass * g;
    return i;
}

float force_3(float u, float mass, float angle)
{
    float g = 9.8;
    printf("Enter the values of static coefficient, mass, and angle of inclination(in deg) : ");
    scanf("%f %f %f", &u, &mass, &angle);
    float rad = angle * 3.14 / 180;
    float i = u * mass * g * cos(rad);
    return i;
}

void friction()
{
    int n = 0;
    float mass, nf, ff, angle, u;
    printf("\n1.Force of Friction \n2.Normal reaction on plane surface \n3.Normal reaction on inclined surface\n\n");
    printf("Choose a Calculation for Frictional force: ");
    scanf("%d", &n);
    if (n == 1)
    {
        ff = force_1(u, nf);
        printf("\nForce of friction = %.4f N\n", ff);
    }
    else if (n == 2)
    {
        ff = force_2(u, mass);
        printf("Force of friction = %.4f N\n", ff);
    }
    else if (n == 3)
    {
        ff = force_3(u, mass, angle);
        printf("Force of friction = %.4f N\n", ff);
    }
}

int main()
{

    while (1)
    {
        fflush(stdin);
        int choice;
    here:
        printf("\n\nWelcome to THERMOSPINSOLVER.\n\n");

        printf("If you want to close give '0'\n");

        printf("For chosing particualr topic give its serial number.\n");

        printf("1.Work\n2.Kinetic Energy\n3.Potential Energy\n4.Power\n5.Heat Capacity\n6.Heat Transfer\n");

        printf("7.Solid mechanics\n8.Fluid Mechanics\n");

        printf("9.Gravitation\n10.Change in g\n11.Center of mass(calculation)\n12.Center fo mass(formula)\n");

        printf("13.Moment of Inertia\n14.Friction\n");

        printf("15.Motion\n");

        printf("16.Thermodynamics\n");

        printf("Give your choice (serial number): ");
        scanf("%d", &choice);

        fflush(stdin);

        if (choice == 1)
            work();
        else if (choice == 2)
            k_energy();
        else if (choice == 3)
            p_energy();
        else if (choice == 5)
            heat('c');
        else if (choice == 4)
            power();
        else if (choice == 6)
            heat('a');
        else if (choice == 7)
            shivam();
        else if (choice == 8)
            shivam2();
        else if (choice == 9)
            Gravitation();
        else if (choice == 10)
            change_in_g();
        else if (choice == 11)
            COM();
        else if (choice == 12)
            COM_formula();
        else if (choice == 13)
            MoI();
        else if (choice == 14)
            friction();
        else if (choice == 15)
            shubham();
        else if (choice == 16)
            Thermo();
        else if (choice == 0)
            break;
        else
        {
            printf("You are not given any appropriate serial number.\n");
            goto here;
        }
    }
    return 0;
}
