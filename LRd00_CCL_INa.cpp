/* The Luo-Rudy Dynamic (LRd) Model of the Mammalian Ventricular Myocyte */
/* Gregory Faber */
/* This code requires a C++ compiler */
/* Detailed list of equations and model description are provided in */

/* Circ Res 1991;68:1501-1526 */
/* Circ Res 1994;74:1071-1096 */
/* Circ Res 1994;74:1097-1113 */
/* Circ Res 1995;77:140-152 */
/* Biophys J 1995;68:949-964 */
/* Cardiovasc Res 1997;35:256-272 */
/* Circulation 1999;99:2466-2474 */
/* Cardiovas Res 1999;42:530-542 */
/* Nature 1999;400:566-569 */
/* Circulation 2000;101:1192-1198 */
/* Biophy J 2000;78:2392-2404 */


/* IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USER SHOULD PACE THE MODEL UNTIL STEADY-STATE IS REACHED.
THIS REQUIRES APPROXIMATELY 5-20 MINUTES OF PACING DEPENDING ON THE RATE.*/
#include"stdafx.h"
#include <iostream> 
#include <iomanip> 
#include <math.h> 
#include <fstream> 
#include <stdlib.h> 
#include <stdio.h> 
#include <sys/timeb.h>

#define bcl 1000 // Basic Cycle Length (ms) 
#define beats 3 // Number of Beats 
//= == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
//Markov INa99 model
//@Colleen E.Clancy & Yoram Rudy(Natrue,1999;400:566-569)
//|------------Markov process of sodium channel------|
//|							  I_F<=>I_S				 |
//|		 	                // ||					 |
//|		 	 C_3<=>C_2<=>C_1<=>O_p					 |
//|--------------------------------------------------|
//@2017-12-1 Chen Xing-Ji added//
//= == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
//以下这些初值从哪来？
//double C_3;// Na C3 rate time constant
//double C_2;//Na C2 rate time constant
//double C_1; // Na C3 rate time constant
//double O_p;  //Na Op rate time constanttau
//double I_F;  //Na IF  rate time constant
double I_S;   //Na IS rate time constant

//double dC_3;    //Na C3 rate time constant
//double dC_2;    //Na C2 rate time constant
//double dC_1;    //Na C3 rate time constant
//double dO_p;    //Na Op rate time constanttau
//double dI_F;    //Na IF rate time constant
double dI_S;    //Na IS rate time constant

//double d2C_3;    //Na C3 rate time constant
//double d2C_2;    //Na C2 rate time constant
//double d2C_1;    //Na C3 rate time constant
//double d2O_p;    //Na Op rate time constanttau
//double d2I_F;    //Na IF rate time constant
double d2I_S;    //Na IS rate time constant
/*
double a11;//C_3->C_2;
double a12;//C_2->C_1;
double a13;//C_1->O_p
double b11;//C_2->C_3;
double b12;//C_1->C_2;
double b13;//O_p->C_1
double a2;//O_p->I_F
double b2;//I_F->O_p
double a3;//I_F->C_1;
double b3;//C_1->I_F;
double a4;//I_F->I_S
double b4;//I_S->I_F*/
double P1, P2;//total probability
double P_I;//, P_C;//zone total probability
double O_p1;    //Na Op rate time constanttau
double Ona;    //Na sodium channel open probability
double loctmp;
double dtold;
//= == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
//Markov INa02 model
//@Colleen E.Clancy & Yoram Rudy(Circulation.2002;105;1208-1213)
//|------------Markov process of sodium channel------|
//|		    IC_3<=>IC_2<=>I_F<=>IM_1<=>IM_2			 |
//|		 	||     ||    ||  \\						 |
//|		 	C_3<=> C_2<=>C_1<=>O_p					 |
//|--------------------------------------------------|
//@2017-11-29 Chen Xing-Ji added//
//= == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
//以下这些初值从哪来？
double C_3;// Na C3 rate time constant
double C_2;  //Na C2 rate time constant
double C_1;// Na C3 rate time constant
double O_p;  //Na Op rate time constanttau
double IC_3;    //Na IC3 rate time constant
double IC_2;  //Na IC2 rate time constant
double I_F;  //Na IForward rate time constant
double IM_1;   //Na IM1 rate time constant
double IM_2;    //Na IM2 rate time constant

double a11;//C_3->C_2;IC_3->IC_2;
double a12;//C_2->C_1;IC_2->I_F;
double a13;//C_1->O_p
double b11;//C_2->C_3;IC_2->IC_3;
double b12;//C_1->C_2;I_F->IC_2;
double b13;//O_p->C_1
double a2;//O_p->I_F
double b2;//I_F->O_p
double a3;//I_F->C_1;IC_3->C_3;IC_2->C_2;
double b3;//C_1->I_F;C_3->IC_3;C_2->IC_2;
double a4;//I_F->IM_1
double b4;//IM_1->I_F
double a5;//IM_1->IM_2//If this is a Markov INa99 model ,here are not have the terms about a5 and b5.
double b5;//IM_2->IM_1

double dC_3;    //Na C3 rate time constant
double dC_2;    //Na C2 rate time constant
double dC_1;    //Na C3 rate time constant
double dO_p;    //Na Op rate time constanttau
double dIC_3;    //Na IC3 rate time constant
double dIC_2;    //Na IC2 rate time constant
double dI_F;   //Na IForward rate time constant
double dIM_1;    //Na IM1 rate time constant
double dIM_2;    //Na IM2 rate time constant
double d2C_3;    //Na C3 rate time constant
double d2C_2;    //Na C2 rate time constant
double d2C_1;    //Na C3 rate time constant
double d2O_p;    //Na Op rate time constanttau
double d2IC_3;    //Na IC3 rate time constant
double d2IC_2;   //Na IC2 rate time constant
double d2I_F;   //Na IForward rate time constant
double d2IM_1;    //Na IM1 rate time constant
double d2IM_2;    //Na IM2 rate time constant
//double P1, P2;//total probability
double P_IM, P_IC, P_C;//zone total probability
/* Numerical variables */
int sti_flag = -1;// turn to 1 if stimulus is applied
const double Voffset = 0.1;//mv, a given membrane potential offset
double D1V = 0, D1V_old = 0, D2V = 0;
double D2V_nai, D2V_ki, D2V_cai, D2V_nsr, D2V_djsr;
double v_old;
double dt_rem; //the left time of the current cycle
int b_n = 0;// beat count
const double dt_univ = 0.01;// universal maximum time-step size (ms) 
const double dt_min = 0.001;//ms
double epsilon = 1e-10;

/* List of variables and parameters (this code uses all global variables) */

/* Creation of Data File */
FILE *ap;
FILE *ap94;
FILE *ap99;
FILE *ap02;
FILE *fmaxs;
FILE *fpara;
FILE *CPU_time_CCL;

/* Cell Geometry */
const double l = 0.01; // Length of the cell (cm) 
const double a = 0.0011; // Radius of the cell (cm) 
const double pi = 3.141592; // Pi 
double vcell; // Cell volume (uL) 
double ageo; // Geometric membrane area (cm^2) 
double acap; // Capacitive membrane area (cm^2) 
double vmyo; // Myoplasm volume (uL) 
double vmito; // Mitochondria volume (uL) 
double vsr; // SR volume (uL) 
double vnsr; // NSR volume (uL) 
double vjsr; // JSR volume (uL) 
double vcleft; // Cleft volume (uL) 

/* Voltage */
double v; // Membrane voltage (mV) 
double vnew; // New Voltage (mV) 
double dvdt; // Change in Voltage / Change in Time (mV/ms) 
double dvdtnew; // New dv/dt (mV/ms) 
double flag; // Flag condition to test for dvdtmax

/* Time Step */
double dt; // Time step (ms) 
double t; // Time (ms) 
int steps; // Number of Steps 
int increment; // Loop Control Variable

/* Action Potential Duration and Max. Info */
double vmax[beats]; // Max. Voltage (mV) 
double dvdtmax[beats]; // Max. dv/dt (mV/ms) 
double apd[beats]; // Action Potential Duration 
double toneapd[beats]; // Time of dv/dt Max. 
double ttwoapd[beats]; // Time of 90% Repolarization 
double rmbp[beats]; // Resting Membrane Potential 
double nair[beats]; // Intracellular Na At Rest 
double cair[beats]; // Intracellular Ca At Rest 
double kir[beats]; // Intracellular K At Rest 
double caimax[beats]; // Peak Intracellular Ca 
int i; // Stimulation Counter

/* Total Current and Stimulus */
double st; // Constant Stimulus (uA/cm^2) 
double tstim; // Time Stimulus is Applied (ms) 
double stimtime; // Time period during which stimulus is applied (ms) 
double it; // Total current (uA/cm^2) 
double duration; // duration of the stimulus (ms) 

/* Terms for Solution of Conductance and Reversal Potential */
const double R = 8314; // Universal Gas Constant (J/kmol*K) 
const double frdy = 96485; // Faraday's Constant (C/mol) 
const double temp = 310; // Temperature (K) 

/* Ion Valences */
const double zna = 1; // Na valence 
const double zk = 1; // K valence 
const double zca = 2; // Ca valence 

/* Ion Concentrations */
double nai; // Intracellular Na Concentration (mM) 
double nao; // Extracellular Na Concentration (mM) 
double nabm; // Bulk Medium Na Concentration (mM) 
double dnao; // Change in Cleft Na Concentration (mM) 
double ki; // Intracellular K Concentration (mM) 
double ko; // Extracellular K Concentration (mM) 
double kbm; // Bulk Medium K Concentration (mM) 
double dko; // Change in Cleft K Concentration (mM) 
double cai; // Intracellular Ca Concentration (mM) 
double cao; // Extracellular Ca Concentration (mM) 
double cabm; // Bulk Medium Ca Concentration (mM) 
double dcao; // Change in Cleft Ca Concentration (mM) 
double cmdn; // Calmodulin Buffered Ca Concentration (mM) 
double trpn; // Troponin Buffered Ca Concentration (mM) 
double nsr; // NSR Ca Concentration (mM) 
double jsr; // JSR Ca Concentration (mM) 
double csqn; // Calsequestrin Buffered Ca Concentration (mM) 
const double taudiff = 1000; // Diffusion Constant for Ion Movement from Bulk Medium to Cleft Space

/* Myoplasmic Na Ion Concentration Changes */
double naiont; // Total Na Ion Flow (uA/uF) 
double dnai; // Change in Intracellular Na Concentration (mM) 

/* Myoplasmic K Ion Concentration Changes */
double kiont; // Total K Ion Flow (uA/uF) 
double dki; // Change in Intracellular K Concentration (mM) 

/* NSR Ca Ion Concentration Changes */
double dnsr; // Change in [Ca] in the NSR (mM) 
double iup; // Ca uptake from myo. to NSR (mM/ms) 
double ileak; // Ca leakage from NSR to myo. (mM/ms) 
double kleak; // Rate constant of Ca leakage from NSR (ms^-1) 
const double kmup = 0.00092; // Half-saturation concentration of iup (mM) 
const double iupbar = 0.00875; // Max. current through iup channel (mM/ms) 
const double nsrbar = 15; // Max. [Ca] in NSR (mM)

/* JSR Ca Ion Concentration Changes */
double djsr; // Change in [Ca] in the JSR (mM) 
const double tauon = 0.5; // Time constant of activation of Ca release from JSR (ms) 
const double tauoff = 0.5; // Time constant of deactivation of Ca release from JSR (ms) 
double tcicr; // t=0 at time of CICR (ms) 
double irelcicr; // Ca release from JSR to myo. due to CICR (mM/ms) 
const double csqnth = 8.75; // Threshold for release of Ca from CSQN due to JSR overload (mM) 
const double gmaxrel = 150; // Max. rate constant of Ca release from JSR due to overload (ms^-1) 
double grelbarjsrol; // Rate constant of Ca release from JSR due to overload (ms^-1) 
double greljsrol; // Rate constant of Ca release from JSR due to CICR (ms^-1) 
double tjsrol; // t=0 at time of JSR overload (ms) 
double ireljsrol; // Ca release from JSR to myo. due to JSR overload (mM/ms) 
const double csqnbar = 10; // Max. [Ca] buffered in CSQN (mM) 
const double kmcsqn = 0.8; // Equilibrium constant of buffering for CSQN (mM) 
double bjsr; // b Variable for analytical computation of [Ca] in JSR (mM) 
double cjsr; // c Variable for analytical computation of [Ca] in JSR (mM) 
double on; // Time constant of activation of Ca release from JSR (ms) 
double off; // Time constant of deactivation of Ca release from JSR (ms) 
double magrel; // Magnitude of Ca release 
double dcaiont; // Rate of change of Ca entry 
double dcaiontnew; // New rate of change of Ca entry 
double caiontold; // Old rate of change of Ca entry 

/* Translocation of Ca Ions from NSR to JSR */
double itr; // Translocation current of Ca ions from NSR to JSR (mM/ms) 
const double tautr = 180; // Time constant of Ca transfer from NSR to JSR (ms) 

/* Myoplasmic Ca Ion Concentration Changes */
double caiont; // Total Ca Ion Flow (uA/uF) 
double dcai; // Change in myoplasmic Ca concentration (mM) 
double catotal; // Total myoplasmic Ca concentration (mM) 
double bmyo; // b Variable for analytical computation of [Ca] in myoplasm (mM) 
double cmyo; // c Variable for analytical computation of [Ca] in myoplasm (mM) 
double dmyo; // d Variable for analytical computation of [Ca] in myoplasm (mM) 
double gpig; // Tribute to all the guinea pigs killed for the advancement of knowledge 
const double cmdnbar = 0.050; // Max. [Ca] buffered in CMDN (mM) 
const double trpnbar = 0.070; // Max. [Ca] buffered in TRPN (mM) 
const double kmcmdn = 0.00238; // Equilibrium constant of buffering for CMDN (mM) 
const double kmtrpn = 0.0005; // Equilibrium constant of buffering for TRPN (mM) 

/* Fast Sodium Current (time dependant) */
double ina; // Fast Na Current (uA/uF) 
double gna; // Max. Conductance of the Na Channel (mS/uF) 
double ena; // Reversal Potential of Na (mV) 
double am; // Na alpha-m rate constant (ms^-1) 
double bm; // Na beta-m rate constant (ms^-1) 
double ah; // Na alpha-h rate constant (ms^-1) 
double bh; // Na beta-h rate constant (ms^-1) 
double aj; // Na alpha-j rate constant (ms^-1) 
double bj; // Na beta-j rate constant (ms^-1) 
double mtau; // Na activation 
double htau; // Na inactivation 
double jtau; // Na inactivation 
double mss; // Na activation 
double hss; // Na inactivation 
double jss; // Na inactivation 
double m; // Na activation 
double h; // Na inactivation 
double jj; // Na inactivation 

/* Current through L-type Ca Channel */
double ilca; // Ca current through L-type Ca channel (uA/uF) 
double ilcana; // Na current through L-type Ca channel (uA/uF) 
double ilcak; // K current through L-type Ca channel (uA/uF) 
double ilcatot; // Total current through the L-type Ca channel (uA/uF) 
double ibarca; // Max. Ca current through Ca channel (uA/uF) 
double ibarna; // Max. Na current through Ca channel (uA/uF) 
double ibark; // Max. K current through Ca channel (uA/uF) 
double d; // Voltage dependant activation gate 
double dss; // Steady-state value of activation gate d 
double taud; // Time constant of gate d (ms^-1) 
double f; // Voltage dependant inactivation gate 
double fss; // Steady-state value of inactivation gate f 
double tauf; // Time constant of gate f (ms^-1) 
double fca; // Ca dependant inactivation gate 
const double kmca = 0.0006; // Half-saturation concentration of Ca channel (mM) 
const double pca = 0.00054; // Permeability of membrane to Ca (cm/s) 
const double gacai = 1; // Activity coefficient of Ca 
const double gacao = 0.341; // Activity coefficient of Ca 
const double pna = 0.000000675; // Permeability of membrane to Na (cm/s) 
const double ganai = 0.75; // Activity coefficient of Na 
const double ganao = 0.75; // Activity coefficient of Na 
const double pk = 0.000000193; // Permeability of membrane to K (cm/s) 
const double gaki = 0.75; // Activity coefficient of K 
const double gako = 0.75; // Activity coefficient of K 

/* Current through T-type Ca Channel */
double icat; // Ca current through T-type Ca channel (uA/uF) 
double gcat; // Max. Conductance of the T-type Ca channel (mS/uF) 
double eca; // Reversal Potential of the T-type Ca channel (mV) 
double b; // Voltage dependant activation gate 
double bss; // Steady-state value of activation gate b 
double taub; // Time constant of gate b (ms^-1) 
double g; // Voltage dependant inactivation gate 
double gss; // Steady-state value of inactivation gate g 
double taug; // Time constant of gate g (ms^-1) 

/* Potassium Current (time-dependent) in LRd94*/
//double x;
//double gk;
//double ek;
//double xi;
//double ax;
//double bx;
//double taux;
//double xss;
//double ik;
//const double prnak = 0.01833;

/* Rapidly Activating Potassium Current */
double ikr; // Rapidly Activating K Current (uA/uF) 
double gkr; // Channel Conductance of Rapidly Activating K Current (mS/uF) 
double ekr; // Reversal Potential of Rapidly Activating K Current (mV) 
double xr; // Rapidly Activating K time-dependant activation 
double xrss; // Steady-state value of inactivation gate xr 
double tauxr; // Time constant of gate xr (ms^-1) 
double r; // K time-independent inactivation

/* Slowly Activating Potassium Current */
double iks; // Slowly Activating K Current (uA/uF) 
double gks; // Channel Conductance of Slowly Activating K Current (mS/uF) 
double eks; // Reversal Potential of Slowly Activating K Current (mV) 
double xs1; // Slowly Activating K time-dependant activation 
double xs1ss; // Steady-state value of inactivation gate xs1 
double tauxs1; // Time constant of gate xs1 (ms^-1) 
double xs2; // Slowly Activating K time-dependant activation 
double xs2ss; // Steady-state value of inactivation gate xs2 
double tauxs2; // Time constant of gate xs2 (ms^-1) 
const double prnak = 0.01833; // Na/K Permeability Ratio

/* Potassium Current (time-independent) */
double iki; // Time-independent K current (uA/uF) 
double gki; // Channel Conductance of Time Independant K Current (mS/uF) 
double eki; // Reversal Potential of Time Independant K Current (mV) 
double aki; // K alpha-ki rate constant (ms^-1) 
double bki; // K beta-ki rate constant (ms^-1) 
double kin; // K inactivation 

/* Plateau Potassium Current */
double ikp; // Plateau K current (uA/uF) 
double gkp; // Channel Conductance of Plateau K Current (mS/uF) 
double ekp; // Reversal Potential of Plateau K Current (mV) 
double kp; // K plateau factor 

/* Na-Activated K Channel */
double ikna; // Na activated K channel 
double pona; // Open probability dependant on Nai 
double pov; // Open probability dependant on Voltage 
double ekna; // Reversal potential 
const double gkna = 0.12848; // Maximum conductance (mS/uF) 
const double nkna = 2.8; // Hill coefficient for Na dependance 
const double kdkna = 66; // Dissociation constant for Na dependance(mM)

/* ATP-Sensitive K Channel */
double ikatp; // ATP-sensitive K current (uA/uF) 
double ekatp; // K reversal potential (mV) 
double gkbaratp; // Conductance of the ATP-sensitive K channel (mS/uF) 
double gkatp; // Maximum conductance of the ATP-sensitive K channel (mS/uF) 
double patp; // Percentage availability of open channels 
const double natp = 0.24; // K dependence of ATP-sensitive K current 
const double nicholsarea = 0.00005; // Nichol's area (cm^2) 
const double atpi = 3; // Intracellular ATP concentraion (mM) 
const double hatp = 2; // Hill coefficient 
const double katp = 0.250; // Half-maximal saturation point of ATP-sensitive K current (mM)

/* Ito Transient Outward Current (Dumaine et al. Circ Res 1999;85:803-809) */
double ito; // Transient outward current 
double gitodv; // Maximum conductance of Ito 
double ekdv; // Reversal Potential of Ito 
double rvdv; // Time independent voltage dependence of Ito 
double zdv; // Ito activation 
double azdv; // Ito alpha-z rate constant 
double bzdv; // Ito beta-z rate constant 
double tauzdv; // Time constant of z gate 
double zssdv; // Steady-state value of z gate 
double ydv; // Ito inactivation 
double aydv; // Ito alpha-y rate constant 
double bydv; // Ito beta-y rate constant 
double tauydv; // Time constant of y gate 
double yssdv; // Steady-state value of y gate

/* Sodium-Calcium Exchanger V-S */
double inaca; // NaCa exchanger current (uA/uF) 
const double c1 = 0.00025; // Scaling factor for inaca (uA/uF) 
const double c2 = 0.0001; // Half-saturation concentration of NaCa exhanger (mM) 
const double gammas = 0.15; // Position of energy barrier controlling voltage dependence of inaca 

///* Sodium-Calcium Exchanger V-S in LRd94 be replaced by LRd00*/  
//double inaca; // NaCa exchanger current (uA/uF) 
//const double knaca = 2000; // Scaling factor for inaca (uA/uF) 
//const double ksat = 0.1; // Half-saturation concentration of NaCa exhanger (mM) 
//const double gamma = 0.35; //eta=0.35 in the paper// Position of energy barrier controlling voltage dependence of inaca 
//const double kmna = 87.5; // Half-saturation concentration of Na channel (mM)
//const double kmca2 = 1.38; // Half-saturation concentration of Ca channel (mM)

/* Sodium-Potassium Pump */
double inak; // NaK pump current (uA/uF) 
double fnak; // Voltage-dependence parameter of inak 
double sigma; // [Na] o dependence factor of fnak 
/* In unpublished changes on rudy's website(ref: History of LRd Model Development), ibarnak is changed from 1.5(LRd94) to 2.*/
const double ibarnak = 2.25; // Max. current through Na-K pump (uA/uF) 
const double kmnai = 10; // Half-saturation concentration of NaK pump (mM) 
const double kmko = 1.5; // Half-saturation concentration of NaK pump (mM)

/* Nonspecific Ca-activated Current */
double insna; // Non-specific Na current (uA/uF) 
double insk; // Non-specific K current (uA/uF) 
double ibarnsna; // Max. Na current through NSCa channel (uA/uF) 
double ibarnsk; // Max. K current through NSCa channel (uA/uF) 
const double pnsca = 0.000000175; // Permeability of channel to Na and K (cm/s) 
const double kmnsca = 0.0012; // Half-saturation concentration of NSCa channel (mM) 

/* Sarcolemmal Ca Pump */
double ipca; // Sarcolemmal Ca pump current (uA/uF) 
const double ibarpca = 1.15; // Max. Ca current through sarcolemmal Ca pump (uA/uF) 
const double kmpca = 0.0005; // Half-saturation concentration of sarcolemmal Ca pump (mM)

/* Ca Background Current */
double icab; // Ca background current (uA/uF) 
double gcab; // Max. conductance of Ca background (mS/uF) 
double ecan; // Nernst potential for Ca (mV) 

/* Na Background Current */
double inab; // Na background current (uA/uF) 
double gnab; // Max. conductance of Na background (mS/uF) 
double enan; // Nernst potential for Na (mV) 

/* Ion Current Functions */
void comp_ina(); // Calculates Fast Na Current 
void comp_ical(); // Calculates Currents through L-Type Ca Channel 
void comp_icat(); // Calculates Currents through T-Type Ca Channel 
//comp_ik();
void comp_ikr(); // Calculates Rapidly Activating K Current 
void comp_iks(); // Calculates Slowly Activating K Current 
void comp_iki(); // Calculates Time-Independent K Current 
void comp_ikp(); // Calculates Plateau K Current 
void comp_ikna(); // Calculates Na-activated K Current 
void comp_ikatp(); // Calculates ATP-Sensitive K Current 
void comp_ito(); // Calculates Transient Outward Current 
void comp_inaca(); // Calculates Na-Ca Exchanger Current 
void comp_inak(); // Calculates Na-K Pump Current 
void comp_insca(); // Calculates Non-Specific ca-Activated Current 
void comp_ipca(); // Calculates Sarcolemmal Ca Pump Current 
void comp_icab(); // Calculates Ca Background Current 
void comp_inab(); // Calculates Na Background Current 
void comp_it(); // Calculates Total Current 

/* Ion Concentration Functions */
void conc_nai(); // Calculates new myoplasmic Na ion concentration 
void conc_ki(); // Calculates new myoplasmic K ion concentration 
void conc_nsr(); // Calculates new NSR Ca ion concentration 
void conc_jsr(); // Calculates new JSR Ca ion concentration 
void calc_itr(); // Calculates Translocation of Ca from NSR to JSR 
void conc_cai(); // Calculates new myoplasmic Ca ion concentration 
void conc_cleft(); // Calculates new cleft ion concentrations 

/* Additive funtions */
void comp_revpots();//Calculating new Reversal Potential
void comp_currents();//a package of current functions for updating new currents
void comp_concentrations();//a package for calaculating ion concentrations
void comp_CCL(); // Calculating new time-step size
void Rush_Larsen();//Calculating new Gatings
void comp_voltage();//Calculating new voltage
void prttofile();

void Init_State99();
void Markov_ina99();
void Init_State02();
void Markov_ina02();
int main()
{
	/* Opening of Datafiles */
	ap = fopen("ap", "w");
	fpara = fopen("fpara", "w");
	fmaxs = fopen("fmaxs", "w");
	CPU_time_CCL = fopen("CPU_time_CCL.dat", "w");
	//ap94 = fopen("ap94", "w");
	//Itoap94 = fopen("Itoap94", "w");
	ap99 = fopen("ap99", "w");
	//ap02 = fopen("ap02", "w");
	//Itoap02 = fopen("Itoap02", "w");
	/* Cell Geometry */
	vcell = 1000 * pi*a*a*l; // 3.801e-5 uL 
	ageo = 2 * pi*a*a + 2 * pi*a*l; // 7.671e-5 cm^2 
	acap = ageo * 2; // 1.534e-4 cm^2 
	vmyo = vcell*0.68;
	vmito = vcell*0.26;
	vsr = vcell*0.06;
	vnsr = vcell*0.0552;
	vjsr = vcell*0.0048;
	vcleft = vcell*0.12 / 0.88;

	/* Time Loop Conditions */
	t = 0.0; // Time (ms) 
	dt = 0.01; // Time step (ms) 
	steps = int((bcl*beats) / dt); // Number of ms 
	st = -80.0; // Stimulus 
	tstim = 10.0; // Time to begin stimulus 
	stimtime = 10.0; // Initial Condition for Stimulus 
	duration = 0.5;// duration of the stimulus (ms) 
	v = -88.654973; // Initial Voltage (mv) 

	/* Beginning Ion Concentrations */
	nai = 12.236437; // Initial Intracellular Na (mM) 
	nao = 140; // Initial Extracellular Na (mM) 
	nabm = 140; // Initial Bulk Medium Na (mM) 
	ki = 136.89149; // Initial Intracellular K (mM) 
	ko = 4.5; // Initial Extracellular K (mM) 
	kbm = 4.5; // Initial Bulk Medium K (mM) 
	cai = 0.000079; // Initial Intracellular Ca (mM) 
	cao = 1.8; // Initial Extracellular Ca (mM) 
	cabm = 1.8; // Initial Bulk Medium Ca (mM) 

	/* Initial Gate Conditions */
	m = 0.000838;
	h = 0.993336;
	jj = 0.995484;
	d = 0.000003;
	f = 0.999745;
	xs1 = 0.004503;
	xs2 = 0.004503;
	xr = 0.000129;
	b = 0.000994;
	g = 0.994041;
	zdv = 0.0120892;
	ydv = 0.999978;

	/* Initial Conditions */
	grelbarjsrol = 0;
	tjsrol = 1000;
	tcicr = 1000;
	jsr = 1.179991;
	nsr = 1.179991;
	trpn = 0.0143923;
	cmdn = 0.00257849;
	csqn = 6.97978;
	flag = 0;
	dcaiont = 0;
	i = -1;

	struct timeb start, end;
	int diff;
	ftime(&start);

	comp_revpots();
	comp_currents();
	comp_concentrations();//Calculating new concentrations
	stimtime = stimtime + dt;
	Rush_Larsen();//Calculating new Gatings
	D1V = -it;
	comp_voltage();//Calculating new voltage
	t = t + dt;
	/* Initial State Conditions */
	Init_State99();
	//Init_State02();
	/* Beginning of Time Loop */
	for (increment = 1; t <= bcl*beats; increment++)
	{
		v_old = v;
		comp_revpots();
		comp_currents();
		comp_CCL(); // Calculating new time-step size
		comp_concentrations();//Calculating new concentrations
		stimtime = stimtime + dt;
		Rush_Larsen();//Calculating new Gatings
		comp_voltage();//Calculating new voltage

		vnew = v;
		//vnew = v - it*dt;
		dvdtnew = (vnew - v_old) / dt;
		//v = vnew;
		if (i >= 0){
			if (vnew>vmax[i])
				vmax[i] = vnew;
			if (cai>caimax[i])
				caimax[i] = cai;
			if (dvdtnew>dvdtmax[i]){
				dvdtmax[i] = dvdtnew;
				toneapd[i] = t;
			}
			if (vnew >= (vmax[i] - 0.9*(vmax[i] - rmbp[i])))
				ttwoapd[i] = t;
		}

		if (csqn >= csqnth && tjsrol>50){
			grelbarjsrol = 4;
			tjsrol = 0;
			printf("Spontaneous Release occurred at time: %g\n", t);
			//cout << "Spontaneous Release occurred at time " << t << endl;
		}

		dvdt = dvdtnew;
		caiontold = caiont;
		dcaiont = dcaiontnew;
		t = t + dt;
		prttofile();
	}
	ftime(&end);
	diff = (int)(1000.0*(end.time - start.time) + (end.millitm - start.millitm));

	fprintf(fpara, "%.3f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t"
		"%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n",
		t, v, nai, ki, cai, jsr, nsr, nao, ko, cao, m, h, jj, d, f, xs1, xs2, xr, b, g, tcicr, flag);
	for (i = 0; i<beats; i++)
	{
		apd[i] = ttwoapd[i] - toneapd[i];
		fprintf(fmaxs, "%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", i, vmax[i], dvdtmax[i], apd[i],
			toneapd[i], ttwoapd[i], nair[i], kir[i], cair[i], caimax[i], rmbp[i]);
	}

	fprintf(CPU_time_CCL, "%g\n", diff);
	//cout << "\a\a";
	fclose(ap);
	fclose(ap99);
	//fclose(ap02);
	//fclose(Itoap02);
	return (1);
}


/********************************************************/
//calaculate Reversal Potential
void comp_revpots(){
	ena = ((R*temp) / frdy)*log(nao / nai);
	eca = (R*temp / (2 * frdy))*log(cao / cai);
	ekr = ((R*temp) / frdy)*log(ko / ki);
	eks = ((R*temp) / frdy)*log((ko + prnak*nao) / (ki + prnak*nai));
	eki = ((R*temp) / frdy)*log(ko / ki);
	ekna = ((R*temp) / frdy)*log(ko / ki);
	ekatp = ((R*temp) / frdy)*log(ko / ki);
	//ekdv = ((R*temp) / frdy)*log((ko + prnak*nao) / (ki + prnak*nai)); //LRd94
	ekdv = ((R*temp) / frdy)*log(ko / ki); //LRd99
	ecan = ((R*temp) / (2 * frdy))*log(cao / cai);
}

/*  calaculate new Gatings */
void Rush_Larsen(){
	m = mss - (mss - m)*exp(-dt / mtau);
	h = hss - (hss - h)*exp(-dt / htau);
	jj = jss - (jss - jj)*exp(-dt / jtau);

	d = dss - (dss - d)*exp(-dt / taud);
	f = fss - (fss - f)*exp(-dt / tauf);

	b = bss - (bss - b)*exp(-dt / taub);
	g = gss - (gss - g)*exp(-dt / taug);

	xr = xrss - (xrss - xr)*exp(-dt / tauxr);

	xs1 = xs1ss - (xs1ss - xs1)*exp(-dt / tauxs1);
	xs2 = xs2ss - (xs2ss - xs2)*exp(-dt / tauxs2);
	zdv = zssdv - (zssdv - zdv)*exp(-dt / tauzdv);
	ydv = yssdv - (yssdv - ydv)*exp(-dt / tauydv);
}

//Call This Subroutine After Stimulus
//the new dt should not be larger than the left time of the current cycle
void Comp_Time_Left(){
	if (sti_flag < 0)
		dt_rem = tstim + b_n*bcl - t;
	else
		dt_rem = tstim + (1 + b_n)*bcl - t;
}

void comp_voltage(){
	v = v + dt*D1V + dt*dt * D2V / 2;
	//!WRITE(6, *) "[D1V]", D1V
	//	v = v + dt*D1V + dt**2 * D2V / 2d0
	//	!WRITE(6, *)"{voltage_WU}", "[voltage]", v, "[D1V]", D1V, "[D2V]", D2V

	//	!!nai = nai + dt*nai_dot + dt**2 * D2V_nai / 2d0
	//	!WRITE(6, *)"{nai}", nai

	//	!!ki = ki + dt*ki_dot + dt**2 * D2V_ki / 2d0
	//	!WRITE(6, *)"{ki}", ki

	//	!!!djsr = djsr + dt*djsr_dot + dt**2 * D2V_djsr / 2d0  !!***!!
	//	!WRITE(6, *)"{jsr}", jsr

	//	!!nsr = nsr + dt*nsr_dot + dt**2 * D2V_nsr / 2d0
	//	!WRITE(6, *)"{nsr}", nsr

	//	!!!cai = cai + dt*cai_dot + dt**2 * D2V_cai / 2d0
	//	!WRITE(6, *)"{cai}", cai
}

// Calculating new time-step size
void comp_CCL(){
	double dt_range;
	if (dt_univ > dt * 2){//Min(2*dt,dt_univ)
		dt_range = dt * 2;
	}
	else{
		dt_range = dt_univ;
	}

	//the new dt should not be larger than the left time of the current cycle
	//Comp_Time_Left();
	//if (dt_range > dt_rem){
	//	dt_range = dt_rem;
	//}

	D1V_old = D1V;
	D1V = -(it);

	D2V = (D1V - D1V_old) / dt;
	double DiscriminantP = 0, DiscriminantN = 0, dtz = 0;
	if ((sti_flag >0) && (t + dt_min > bcl*(b_n)-epsilon)){
		//WRITE(6, *) n, t, dt_loc, start + (n)*CL, start + (n)*CL - t
		dt = tstim + bcl*b_n - t;/// not larger than the left time of the cycle
		printf("error !!!!");
		exit(1);
	}
	else{
		if (D1V >= 0){
			DiscriminantP = D1V*D1V + 2 * D2V*Voffset;
			if (D2V>0){
				dt = (-D1V + sqrt(DiscriminantP)) / D2V;
			}
			else if (D2V<0){
				dtz = -D1V / D2V;
				if (DiscriminantP >= 0){
					dt = (-D1V + sqrt(DiscriminantP)) / D2V;
				}
				else{
					dt = dtz;
				}
			}
		}
		else{
			DiscriminantN = D1V*D1V - 2 * D2V*Voffset;
			if (D2V>0){
				dtz = -D1V / D2V;
				if (DiscriminantN >= 0){
					dt = (-D1V - sqrt(DiscriminantN)) / D2V;
				}
				else{
					dt = dtz;
				}
			}
			else if (D2V<0){
				dt = (-D1V - sqrt(DiscriminantN)) / D2V;
			}
		}
		if (dt>dt_range){
			dt = dt_range;
		}
		if (dt<dt_min){
			dt = dt_min;
		}
	}
}

//a package of current functions for updating new currents
void comp_currents(){
	/* List of functions called for each timestep, currents commented
	out are only used when modeling pathological conditions */
	//comp_ina();
	Markov_ina99();
	//Markov_ina02();
	comp_ical();
	comp_icat();
	//comp_ik();
	comp_ikr();
	comp_iks();
	comp_iki();
	comp_ikp();
	//comp_ikna (); 
	//comp_ikatp (); 
	//comp_ito (); 
	comp_inaca();
	comp_inak();
	//comp_insca (); 
	comp_ipca();
	comp_icab();
	comp_inab();
	comp_it();
}

//a package for calaculating ion concentrations
void comp_concentrations(){
	conc_nai();
	conc_ki();
	calc_itr();
	conc_jsr();
	conc_nsr();
	conc_cai();
	//conc_cleft (); 
	/* Cleft Space disabled, if you want to use cleft space, make sure the initial conditions
	of ion concentrations in the bulk medium are the same as the extracellular concentrations */
}

/* Functions that describe the currents begin here */

/* Calculates fast sodium current */
void comp_ina()
{
	gna = 16;
	//ena = ((Rgas*temp) / frdy)*log(nao / nai);

	am = 0.32*(v + 47.13) / (1 - exp(-0.1*(v + 47.13)));
	bm = 0.08*exp(-v / 11);
	if (v < -40)
	{
		ah = 0.135*exp((80 + v) / -6.8);
		bh = 3.56*exp(0.079*v) + 310000 * exp(0.35*v);
		aj = (-127140 * exp(0.2444*v) - 0.00003474*exp(-0.04391*v))*((v + 37.78) / (1 + exp(0.311*(v + 79.23))));
		bj = (0.1212*exp(-0.01052*v)) / (1 + exp(-0.1378*(v + 40.14)));
	}

	else
	{
		ah = 0;
		bh = 1 / (0.13*(1 + exp((v + 10.66) / -11.1)));
		aj = 0;
		bj = (0.3*exp(-0.0000002535*v)) / (1 + exp(-0.1*(v + 32)));
	}
	mtau = 1 / (am + bm);
	htau = 1 / (ah + bh);
	jtau = 1 / (aj + bj);

	mss = am*mtau;
	hss = ah*htau;
	jss = aj*jtau;

	//m = mss - (mss - m)*exp(-dt / mtau);
	//h = hss - (hss - h)*exp(-dt / htau);
	//jj = jss - (jss - jj)*exp(-dt / jtau);

	ina = gna*m*m*m*h*jj*(v - ena);
}

/* Calculates Currents through L-Type Ca Channel */
void comp_ical()
{
	dss = 1 / (1 + exp(-(v + 10) / 6.24));
	taud = dss*(1 - exp(-(v + 10) / 6.24)) / (0.035*(v + 10));

	fss = (1 / (1 + exp((v + 32) / 8))) + (0.6 / (1 + exp((50 - v) / 20)));
	tauf = 1 / (0.0197*exp(-pow(0.0337*(v + 10), 2)) + 0.02);

	//d = dss - (dss - d)*exp(-dt / taud);
	//f = fss - (fss - f)*exp(-dt / tauf);
	ibarca = pca*zca*zca*((v*frdy*frdy) / (R*temp))
		*((gacai*cai*exp((zca*v*frdy) / (R*temp)) - gacao*cao) / (exp((zca*v*frdy) / (R*temp)) - 1));
	ibarna = pna*zna*zna*((v*frdy*frdy) / (R*temp))
		*((ganai*nai*exp((zna*v*frdy) / (R*temp)) - ganao*nao) / (exp((zna*v*frdy) / (R*temp)) - 1));
	ibark = pk*zk*zk*((v*frdy*frdy) / (R*temp))
		*((gaki*ki*exp((zk*v*frdy) / (R*temp)) - gako*ko) / (exp((zk*v*frdy) / (R*temp)) - 1));
	fca = 1 / (1 + cai / kmca);//Hill coefficient=1 in LRd95, Hill coefficient=2 in LRd94
	ilca = d*f*fca*ibarca;
	ilcana = d*f*fca*ibarna;
	ilcak = d*f*fca*ibark;
	ilcatot = ilca + ilcana + ilcak;
}

/* Calculates Currents through T-Type Ca Channel in LRd95 */
void comp_icat()
{
	bss = 1 / (1 + exp(-(v + 14.0) / 10.8));
	taub = 3.7 + 6.1 / (1 + exp((v + 25.0) / 4.5));

	gss = 1 / (1 + exp((v + 60.0) / 5.6));
	if (v <= 0)
		taug = -0.875*v + 12.0;
	else
		taug = 12.0;

	//b = bss - (bss - b)*exp(-dt / taub);
	//g = gss - (gss - g)*exp(-dt / taug);
	gcat = 0.05;
	//eca = (Rgas*temp / (2 * frdy))*log(cao / cai);
	icat = gcat*b*b*g*(v - eca);
}

/* Calculates Potassium Current (time-dependent) in LRd94*/
//void comp_ik()
//{
//	gk = 0.282*sqrt(ko / 5.4);
//	xi = 1 / (1 + exp((v - 56.26) / 32.1));
//	ax = 0.0000719*(v + 30) / (1 - exp(-0.148*(v + 30)));
//	bx = 0.000131*(v + 30) / (-1 + exp(0.0687*(v + 30)));
//	taux = 1 / (ax + bx);
//	xss = ax*taux;
//	x0 = xss - (xss - x)*exp(-dt / taux);
//	ik = gk*xi*x0*x0*(v - ek);
////}

/* Calculates Rapidly Activating K Current in LRd95 */
void comp_ikr()
{
	gkr = 0.02614*sqrt(ko / 5.4);
	//ekr = ((Rgas*temp) / frdy)*log(ko / ki);

	xrss = 1 / (1 + exp(-(v + 21.5) / 7.5));
	tauxr = 1 / (0.00138*(v + 14.2) / (1 - exp(-0.123*(v + 14.2))) + 0.00061*(v + 38.9) / (exp(0.145*(v + 38.9)) - 1));
	//xr = xrss - (xrss - xr)*exp(-dt / tauxr);
	r = 1 / (1 + exp((v + 9) / 22.4));
	ikr = gkr*xr*r*(v - ekr);
}

/* Calculates Slowly Activating K Current in LRd99 Viswanathan */
void comp_iks()
{
	gks = 0.433*(1 + 0.6 / (1 + pow((0.000038 / cai), 1.4)));
	//eks = ((Rgas*temp) / frdy)*log((ko + prnak*nao) / (ki + prnak*nai));

	xs1ss = 1 / (1 + exp(-(v - 1.5) / 16.7));
	xs2ss = xs1ss;
	tauxs1 = 1 / (0.0000719*(v + 30) / (1 - exp(-0.148*(v + 30))) + 0.000131*(v + 30) / (exp(0.0687*(v + 30)) - 1));
	tauxs2 = 4 * tauxs1;
	//xs1 = xs1ss - (xs1ss - xs1)*exp(-dt / tauxs1);
	//xs2 = xs2ss - (xs2ss - xs2)*exp(-dt / tauxs2);
	iks = gks*xs1*xs2*(v - eks);
}

/* Calculates Potassium Current (time-independent) in LRd94*/
void comp_iki()
{
	gki = 0.75*(sqrt(ko / 5.4));
	//eki = ((Rgas*temp) / frdy)*log(ko / ki);

	aki = 1.02 / (1 + exp(0.2385*(v - eki - 59.215)));
	bki = (0.49124*exp(0.08032*(v - eki + 5.476)) + exp(0.06175*(v - eki - 594.31))) / (1 + exp(-0.5143*(v - eki + 4.753)));
	kin = aki / (aki + bki);
	iki = gki*kin*(v - eki);
}

/* Calculates Plateau K Current */
void comp_ikp()
{
	gkp = 0.00552; //the data in LRd94 paper is 0.0183 ??????????????
	//gkp = 0.0183;
	ekp = eki;

	kp = 1 / (1 + exp((7.488 - v) / 5.98));

	ikp = gkp*kp*(v - ekp);
}

/* Calculates Na-activated K Current in LRd00*/
void comp_ikna()
{
	//ekna = ((Rgas*temp) / frdy)*log(ko / ki);
	pona = 0.85 / (1 + pow((kdkna / nai), 2.8));
	pov = 0.8 - 0.65 / (1 + exp((v + 125) / 15));
	ikna = gkna*pona*pov*(v - ekna);
}

/* Calculates ATP-Sensitive K Current in LRd97 */
void comp_ikatp()
{
	/* Note: If you wish to use this current in your simulations, there are additional */
	/* changes which must be made to the code as detailed in Cardiovasc Res 1997;35:256-272 */

	//ekatp = ((Rgas*temp) / frdy)*log(ko / ki);
	gkatp = 0.000195 / nicholsarea;
	patp = 1 / (1 + (pow((atpi / katp), hatp)));
	gkbaratp = gkatp*patp*(pow((ko / 4), natp));

	ikatp = gkbaratp*(v - ekatp);
}

/* Calculates Transient Outward Current in LRd99(Dumaine99) */
void comp_ito()
{
	gitodv = 0.5;
	//ekdv = ((Rgas*temp) / frdy)*log((ko) / (ki));
	rvdv = exp(v / 100);
	azdv = (10 * exp((v - 40) / 25)) / (1 + exp((v - 40) / 25));
	bzdv = (10 * exp(-(v + 90) / 25)) / (1 + exp(-(v + 90) / 25));
	tauzdv = 1 / (azdv + bzdv);
	zssdv = azdv / (azdv + bzdv);
	//zdv = zssdv - (zssdv - zdv)*exp(-dt / tauzdv);

	aydv = 0.015 / (1 + exp((v + 60) / 5));
	bydv = (0.1*exp((v + 25) / 5)) / (1 + exp((v + 25) / 5));
	tauydv = 1 / (aydv + bydv);
	yssdv = aydv / (aydv + bydv);
	//ydv = yssdv - (yssdv - ydv)*exp(-dt / tauydv);
	ito = gitodv*zdv*zdv*zdv*ydv*rvdv*(v - ekdv);
}

/* Calculates Na-Ca Exchanger Current in LRd00*/
void comp_inaca()
{
	inaca = c1*exp((gammas - 1)*v*frdy / (R*temp))
		*((exp(v*frdy / (R*temp))*nai*nai*nai*cao - nao*nao*nao*cai)
		/ (1 + c2*exp((gammas - 1)*v*frdy / (R*temp))*(exp(v*frdy / (R*temp))*nai*nai*nai*cao + nao*nao*nao*cai)));

	//inaca = knaca / (kmna*kmna*kmna + nao*nao*nao) / (kmca2 + cao) / (1 
	//+ ksat*exp((gamma - 1)*v*Frdy / ( R*Temp)))*(exp(gamma*v*Frdy / ( R*Temp))*nai*nai*nai*cao 
	//- exp((gamma - 1)*v*Frdy / ( Rg*Temp))*nao*nao*nao*cai);
}

/* Calculates Sodium-Potassium Pump Current */
void comp_inak()
{
	sigma = (exp(nao / 67.3) - 1) / 7;

	fnak = 1 / (1 + 0.1245*exp((-0.1*v*frdy) / (R*temp)) + 0.0365*sigma*exp((-v*frdy) / (R*temp)));

	inak = ibarnak*fnak*(1 / (1 + pow(kmnai / nai, 2)))*(ko / (ko + kmko));//in LRd94 paper, index is 1.5, not 2 !!!
}

//LRd94, here insk and insna are not included
/* Calculates Non-Specific ca-Activated Current in LRd94*/
void comp_insca()
{
	//pnsca Permeability of channel to Na and K (cm/s) (equally permeable)
	ibarnsna = pnsca*zna*zna*((v*frdy*frdy) / (R*temp))
		*((ganai*nai*exp((zna*v*frdy) / (R*temp)) - ganao*nao) / (exp((zna*v*frdy) / (R*temp)) - 1));
	ibarnsk = pnsca*zk*zk*((v*frdy*frdy) / (R*temp))
		*((gaki*ki*exp((zk*v*frdy) / (R*temp)) - gako*ko) / (exp((zk*v*frdy) / (R*temp)) - 1));

	insna = ibarnsna / (1 + pow(kmnsca / cai, 3));
	insk = ibarnsk / (1 + pow(kmnsca / cai, 3));
}

/* Calculates Sarcolemmal Ca Pump Current */
void comp_ipca()
{
	ipca = (ibarpca*cai) / (kmpca + cai);
}

/* Calculates Ca Background Current */
void comp_icab()
{
	gcab = 0.003016;
	ecan = ((R*temp) / (2 * frdy))*log(cao / cai);
	icab = gcab*(v - ecan);
}

/* Calculates Na Background Current */
void comp_inab()
{
	gnab = 0.004;//Increase of gNa,b max from 0.00141 in LRd94 to 0.004 in LRd99.
	enan = ena;
	inab = gnab*(v - enan);
}

/* Total sum of currents is calculated here, if the time is between
stimtime = 0 and stimtime = duration, a stimulus is applied */
void comp_it()
{
	naiont = ina + inab + ilcana + 3 * inak + 3 * inaca;//insna is removed(LRd94)
	//kiont below: ik is replaced by ikr + iks(LRd95), ikna(LRd00)+ikatp(LRd97)+insk(LRd94) are removed
	kiont = ikr + iks + iki + ikp + ilcak - 2 * inak + ito;
	caiont = ilca + icab + ipca - 2 * inaca + icat;
	if (dvdtnew > 10 && tcicr > 10 && flag == 1){
		flag = 0;
	}

	if (t >= tstim && t<(tstim + dt)){
		stimtime = 0;
		i = i + 1;
		tstim = tstim + bcl;
		rmbp[i] = v;
		nair[i] = nai;
		kir[i] = ki;
		cair[i] = cai;
		printf("t:%g\n", t);
		//cout << t << endl;

		//for CCL
		b_n = b_n + 1;

	}

	if (stimtime >= 0 && stimtime <= duration){
		it = st + naiont + kiont + caiont;

		sti_flag = 1;//for CCL
	}
	else{
		it = naiont + kiont + caiont;
		sti_flag = -1;//for CCL
	}
}

/* Functions that calculate intracellular ion concentrations begins here */

/* Calculates new myoplasmic Na ion concentration */
void conc_nai()
{
	// The units of dnai is in mM. Note that naiont should be multiplied by the 
	// cell capacitance to get the correct units. Since cell capacitance = 1 uF/cm^2, 
	// it doesn't explicitly appear in the equation below. 
	// This holds true for the calculation of dki and dcai. */ 

	dnai = -dt*(naiont*acap) / (vmyo*zna*frdy);
	nai = dnai + nai;
}

/* Calculates new myoplasmic K ion concentration */
void conc_ki()
{
	if (stimtime >= 0 && stimtime <= duration){
		dki = -dt*((kiont + st)*acap) / (vmyo*zk*frdy);
	}
	else{
		dki = -dt*(kiont*acap) / (vmyo*zk*frdy);
	}

	ki = dki + ki;
}

/* Translocation of Ca2+ ions from NSR to JSR: */
void calc_itr()
{
	//Itr, Translocation current of Ca ions from NSR to JSR (mM/ms) 
	itr = (nsr - jsr) / tautr;
}

/* Calculates new JSR Ca ion concentration */
void conc_jsr()
{
	// there is always calcium release triggered by calcium entry, 
	// and the release becomes significant as calcium entry increases.
	//Ca2+ uptake and leakage of NSR: Iup and Ileak
	kleak = iupbar / nsrbar;
	ileak = kleak*nsr;

	iup = iupbar*cai / (cai + kmup);

	dcaiontnew = (caiont - caiontold) / dt;	// Total Ca Ion Flow (uA/uF) 

	//dcaiont; // Rate of change of Ca entry 
	// half-saturation potential(Vh),-35mV.
	if (v>-35 && dcaiontnew>dcaiont && flag == 0){
		flag = 1;// Flag condition to test for dvdtmax
		tcicr = 0;// t=0 at time of CICR (ms).initial value is 1000, make exponential term approached 0
	}

	// Time constant of activation of Ca release from JSR (ms)
	// tcicr, t=0 at time of CICR (ms)
	// ref to Biophys J 78:2392-404, 2000
	on = 1 / (1 + exp((-tcicr + 4) / tauon));
	off = (1 - 1 / (1 + exp((-tcicr + 4) / tauoff)));
	magrel = 1 / (1 + exp(((ilca + icab + ipca - 2 * inaca + icat) + 5) / 0.9)); // Magnitude of Ca release
	irelcicr = gmaxrel*on*off*magrel*(jsr - cai);// Ca release from JSR to myo. due to CICR (mM/ms)

	tcicr = tcicr + dt;

	/********* irelcicr *************/
	////dcai2 not consider Ca released from SR
	//dcai2 = -dt*(((caiont*acap) / (vmyo*zca*Frdy)) + ((iup - ileak)*vnsr / vmyo));
	////vnew = v - it*dt;
	////dvdtnew = (vnew - v) / dt;

	///*Once it is triggered(at 2 milliseconds), the JSR releases almost all of its free Ca2+ contents.
	//In the model, we compute the increase in [Ca2+]	2 milliseconds (rather than 10 milliseconds)
	//from the time of Vmax and denote it as d[Ca2+]i,2 ......or from the time of stimulus..*/
	//if (vnew > v && vnew > vmax_jsr){//vmax_jsr's initial value is -90
	//	vmax_jsr = vnew;
	//	tvmax = 0;
	//	dcaisum = 0;
	//	tcicr = 0;
	//}
	//else if (firstvmax == 1){// find Vmax!!
	//	flagvmax = 1;
	//	firstvmax = 0;
	//}

	//if (flagvmax == 1){
	//	tvmax = tvmax + dt;  //cumulative amount of time, not larger than 2 msec
	//	dcaisum = dcaisum + dcai2; //cumulative amount of Ca2+
	//	if ((tvmax <= 2 && dcaisum > dcaith && firstcicr == 1)){//firstcicr ensure only one CICR in one action potential duration
	//		tcicr = 0;
	//		flagcicr = 1;
	//		firstcicr = 0;
	//	}
	//}

	//if (flagcicr == 1){
	//	grelcicr = grelcicrbar*(dcaisum - dcaith) / (kmrel + dcaisum - dcaith)*(1 - exp(-tcicr / tauon))*exp(-tcicr / tauoff);
	//}
	//else{
	//	grelcicr = 0;
	//}

	//irelcicr = grelcicr*(jsr - cai);
	//tcicr = tcicr + dt;
	/*********** irelcicr ***********/

	/*********** ireljsrol ***********/
	// Rate constant of Ca release from JSR due to overload (ms^-1) 
	// The threshold of overload is set, can be founded by searching "grelbarjsrol"
	greljsrol = grelbarjsrol*(1 - exp(-tjsrol / tauon))*exp(-tjsrol / tauoff);
	ireljsrol = greljsrol*(jsr - cai);
	tjsrol = tjsrol + dt;
	/*********** ireljsrol ***********/

	/*********** jsr ***********/
	//Method 1. Analytical Computation  ref to Circ Res 1995;77:140-152
	csqn = csqnbar*(jsr / (jsr + kmcsqn));//[CSQN],Calsequestrin Buffered Ca Concentration
	djsr = dt*(itr - irelcicr - ireljsrol);//dCa2+_JSR
	bjsr = csqnbar - csqn - djsr - jsr + kmcsqn;
	cjsr = kmcsqn*(csqn + djsr + jsr);

	jsr = (sqrt(bjsr*bjsr + 4 * cjsr) - bjsr) / 2;// Ca buffering in JSR, new
	/*********** jsr ***********/

	/*********** jsr ***********/
	//Method 2. Steffensen's iterative method
	/*djsr = dt*(itr-irelcicr-ireljsrol);
	jsr1 = jsr+djsr*0.05;

	ii=1;
	while (ii <= 20) // Number of iteration
	{
	csqn1 = csqnbar/(1+kmcsqn/jsr1);
	jsr2 = jsr+djsr-(csqn1-csqn);
	csqn2 = csqnbar/(1+kmcsqn/jsr2);
	jsr3 = jsr+djsr-(csqn2-csqn);
	jsrk1 = jsr1-pow((jsr1-jsr2), 2)/(jsr1- 2*jsr2+jsr3);

	if (abs((jsrk1-jsr1)/jsr1)<0.01){
	jsr = jsrk1;
	csqn = csqnbar/(1+kmcsqn/jsr);
	break;
	} else if (((jsr1-2*jsr2+jsr3)/jsr1)<=0.0001){
	jsr = jsr3;
	csqn = csqnbar/(1+kmcsqn/jsr);
	break;
	} else{
	jsr1 = jsrk1;
	}
	ii++;
	}*/
	/*********** jsr ***********/
}

/* Calculates new NSR Ca ion concentration */
void conc_nsr()
{
	//itr;  Translocation current of Ca ions from NSR to JSR (mM/ms) 
	//vjsr; JSR volume, vnsr; NSR volume 
	dnsr = dt*(iup - ileak - itr*vjsr / vnsr);
	nsr = nsr + dnsr;
}

//take care the total Ca2+ concentration,catotal
//Steffensen Iterative Method can also be used here
/* Calculates total myoplasmic Ca concentration (mM) */
void conc_cai()
{
	//Method 1. Steffensen Iterative Method .....

	/*
	//Method 2.  Forward Euler Method
	B1 = -(((caiont*acap) / (vmyo*zca*Frdy)) + ((iup - ileak)*vnsr / vmyo) - (irelcicr*vjsr / vmyo) - (ireljsrol*vjsr / vmyo));
	B2 = 1 + (trpnbar*kmtrpn / pow((cai + kmtrpn), 2)) + (cmdnbar*kmcmdn / pow((cai + kmcmdn), 2));
	dcai = (B1 / B2)*dt; // Change in myoplasmic Ca concentration (mM)
	cai = cai + dcai;
	*/

	//Myoplasmic Ca Ion Concentration Changes
	dcai = -dt*(((caiont*acap) / (vmyo*zca*frdy)) + ((iup - ileak)*vnsr / vmyo)
		- (irelcicr*vjsr / vmyo) - (ireljsrol*vjsr / vmyo));
	trpn = trpnbar*(cai / (cai + kmtrpn));
	cmdn = cmdnbar*(cai / (cai + kmcmdn));

	//Method 3. Analytical Computation  ref to Circ Res 1995;77:140-152
	//Total myoplasmic Ca concentration (mM)
	catotal = trpn + cmdn + dcai + cai;
	bmyo = cmdnbar + trpnbar - catotal + kmtrpn + kmcmdn;
	cmyo = (kmcmdn*kmtrpn) - (catotal*(kmtrpn + kmcmdn)) + (trpnbar*kmcmdn) + (cmdnbar*kmtrpn);
	dmyo = -kmtrpn*kmcmdn*catotal;
	gpig = sqrt(bmyo*bmyo - 3 * cmyo);

	cai = (2 * gpig / 3)*cos(acos((9 * bmyo*cmyo - 2 * bmyo*bmyo*bmyo - 27 * dmyo)
		/ (2 * pow((bmyo*bmyo - 3 * cmyo), 1.5))) / 3) - (bmyo / 3);
}

/* Cleft Space disabled, if you want to use cleft space, make sure the initial conditions
of ion concentrations in the bulk medium are the same as the extracellular concentrations */
/* Calculates new cleft ion concentrations  */
void conc_cleft()
{
	dnao = dt*((nabm - nao) / taudiff + naiont*acap / (vcleft*frdy));
	nao = dnao + nao;

	if (stimtime >= 0 && stimtime <= duration){
		dko = dt*((kbm - ko) / taudiff + (kiont + st)*acap / (vcleft*frdy));
	}
	else{
		dko = dt*((kbm - ko) / taudiff + kiont*acap / (vcleft*frdy));
	}
	ko = dko + ko;

	dcao = dt*((cabm - cao) / taudiff + caiont*acap / (vcleft*frdy * 2));
	cao = dcao + cao;
}

/* Values are printed to a file called ap. The voltage and
currents can be plotted versus time using graphing software. */
void prttofile()
{
	if (t > (0) && t < (bcl*beats)){
		fprintf(ap, "%.3f\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
			t, v, nai, ki, cai, jsr, nsr, ina, ikr, iks, iki, ilca, icat, inab, icab, ipca, inaca, inak);
		//fprintf(Itoap94, "%.3f\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
		//t, v, nai, ki, cai, jsr, nsr, ina, ikr, iks, iki, ilca, icat, inab, icab, ipca, inaca, inak);
		fprintf(ap99, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", t, v, ina, O_p1, I_S, I_F, P_C);
		//fprintf(ap02, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", t, v, ina, O_p1, P_IM, P_IC, P_C, I_F);
		//fprintf(Itoap02, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", t, v, ina, O_p1, P_IM, P_IC, P_C, I_F);
	}
}

void Init_State99(){
	//-----------init value need to rebuild-----2017.12.1---
	loctmp = 0.0;

	d2C_3 = 0.0;
	d2C_2 = 0.0;
	d2C_1 = 0.0;
	d2O_p = 0.0;
	d2I_F = 0.0;
	d2I_S = 0.0;

	dC_3 = 0.0;
	dC_2 = 0.0;
	dC_1 = 0.0;
	dO_p = 0.0;
	dI_F = 0.0;
	dI_S = 0.0;

	C_3 = 0.98627;// 0.82973 + 0.15654;// Na C3 rate time constant
	C_2 = 0.01363;//0.01147 + 0.00216;  //Na C2 rate time constant
	C_1 = 8.478e-5;//6.011e-5+2.467e-5; // Na C3 rate time constant
	O_p = 5e-8;  //Na Op rate time constanttau
	//IC_3 = 0.15654;    //Na IC3 rate time constant
	//IC_2 = 0.00216;  //Na IC2 rate time constant
	I_F = 1.138e-5;  //Na IForward rate time constant
	I_S = 0;//2.467e-5;   //Na IM1 rate time constant
	//P1 = O_p + C_3 + C_2 + C_1 + I_F + I_S ;
}
void Markov_ina99(){
	dtold = dt;
	gna = 16;
	ena = ((R*temp) / frdy)*log(nao / nai);

	a11 = 3.802 / (0.1027 *exp(-v / 17.0) + 0.20 *exp(-v / 150));
	a12 = 3.802 / (0.1027 *exp(-v / 15.0) + 0.23 *exp(-v / 150));
	a13 = 3.802 / (0.1027 *exp(-v / 12.0) + 0.25 *exp(-v / 150));
	b12 = 0.20 *exp(-(v - 5) / 20.3);
	b11 = 0.1917 *exp(-v / 20.3);
	b13 = 0.22 *exp(-(v - 10) / 20.3);
	a3 = 3.7933e-10 * exp(-v / 5.2);
	b3 = (0.0084 + 0.00002 *v);
	a2 = (9.178 *exp(v / 29.68));
	b2 = (a13*a2*a3) / (b13*b3);
	a4 = a2*0.01;
	b4 = a3;

	loctmp = dC_3;
	dC_3 = (C_2*b11 - a11*C_3);
	d2C_3 = (dC_3 - loctmp) / dtold;

	loctmp = dC_2;
	dC_2 = (C_3*a11 + C_1*b12 - (b11 + a12)*C_2);//错误的原因之前少加一项C_1*b12
	d2C_2 = (dC_2 - loctmp) / dtold;

	loctmp = dC_1;
	dC_1 = (C_2*a12 + I_F*a3 + O_p*b13 - (b12 + b3 + a13)*C_1);
	d2C_1 = (dC_1 - loctmp) / dtold;

	loctmp = dO_p;
	dO_p = (C_1*a13 + I_F*b2 - (b13 + a2)*O_p);
	d2O_p = (dO_p - loctmp) / dtold;

	loctmp = dI_F;
	dI_F = (I_S*b4 + O_p*a2 + C_1*b3 - (a4 + b2 + a3)*I_F);
	d2I_F = (dI_F - loctmp) / dtold;

	loctmp = dI_S;
	dI_S = (I_F*a4 - b4*I_S);
	d2I_S = (dI_S - loctmp) / dtold;

	C_3 = C_3 + dC_3*dt + d2C_3*dt*dt / 2;//改为这个方法试一下m = mss - (mss - m)*exp(-dt / mtau);
	C_2 = C_2 + dC_2*dt + d2C_2*dt*dt / 2;
	C_1 = C_1 + dC_1*dt + d2C_1*dt*dt / 2;
	O_p = O_p + dO_p*dt + d2O_p*dt*dt / 2;
	I_F = I_F + dI_F*dt + d2I_F*dt*dt / 2;
	I_S = I_S + dI_S*dt + d2I_S*dt*dt / 2;

	P_I = I_S + I_F;
	P_C = C_3 + C_2 + C_1;
	//P2 = O_p + C_3 + C_2 + C_1  + I_F + I_S ;
	O_p1 = (1.0 - C_3 - C_2 - C_1 - I_F - I_S);
	Ona = O_p1;
	ina = gna * Ona * (v - ena);

	//return 0;
}
void Init_State02(){
	//loctmp = 0.0;

	d2C_3 = 0.0;
	d2C_2 = 0.0;
	d2C_1 = 0.0;
	d2O_p = 0.0;
	d2IC_3 = 0.0;
	d2IC_2 = 0.0;
	d2I_F = 0.0;
	d2IM_1 = 0.0;
	d2IM_2 = 0.0;

	dC_3 = 0.0;
	dC_2 = 0.0;
	dC_1 = 0.0;
	dO_p = 0.0;
	dIC_3 = 0.0;
	dIC_2 = 0.0;
	dI_F = 0.0;
	dIM_1 = 0.0;
	dIM_2 = 0.0;

	C_3 = 0.82973;// Na C3 rate time constant
	C_2 = 0.01147;  //Na C2 rate time constant
	C_1 = 6.011e-5;// Na C3 rate time constant
	O_p = 5e-8;  //Na Op rate time constanttau
	IC_3 = 0.15654;    //Na IC3 rate time constant
	IC_2 = 0.00216;  //Na IC2 rate time constant
	I_F = 1.138e-5;  //Na IForward rate time constant
	IM_1 = 2.467e-5;   //Na IM1 rate time constant
	IM_2 = 0;    //Na IM2 rate time constant
	//P1 = O_p + C_3 + C_2 + C_1 + IC_3 + IC_2 + I_F + IM_1 + IM_2;
}
void Markov_ina02(){
	dtold = dt;
	gna = 16;
	//ena = ((R*temp) / frdy)*log(nao / nai);

	a11 = 3.802 / (0.1027 *exp(-v / 17.0) + 0.20 *exp(-v / 150));
	a12 = 3.802 / (0.1027 *exp(-v / 15.0) + 0.23 *exp(-v / 150));//与Markov INa99/02有出入0.23。
	a13 = 3.802 / (0.1027 *exp(-v / 12.0) + 0.25 *exp(-v / 150));//与Markov INa99/02有出入0.25。
	b12 = 0.20 *exp(-(v - 5) / 20.3);
	b11 = 0.1917 *exp(-v / 20.3);
	b13 = 0.22 *exp(-(v - 10) / 20.3);
	a3 = 3.7933e-7 * exp(-v / 7.7);//与Markov INa99有出入,与Markov INa02一致。
	b3 = (0.0084 + 0.00002 *v);
	a2 = (9.178 *exp(v / 29.68));
	b2 = (a13*a2*a3) / (b13*b3);
	a4 = a2*0.01;
	b4 = a3;
	a5 = a2 / 95000;// 与Markov INa02有出入a5 = a2 / 95000。
	b5 = a3 / 50;  //与Markov INa02一致。

	loctmp = dC_3;
	dC_3 = (IC_3*a3 + C_2*b11 - (b3 + a11)*C_3);
	d2C_3 = (dC_3 - loctmp) / dtold;

	loctmp = dC_2;
	dC_2 = (C_3*a11 + IC_2*a3 + C_1*b12 - (b11 + b3 + a12)*C_2);
	d2C_2 = (dC_2 - loctmp) / dtold;

	loctmp = dC_1;
	dC_1 = (C_2*a12 + I_F*a3 + O_p*b13 - (b12 + b3 + a13)*C_1);
	d2C_1 = (dC_1 - loctmp) / dtold;

	loctmp = dO_p;
	dO_p = (C_1*a13 + I_F*b2 - (b13 + a2)*O_p);
	d2O_p = (dO_p - loctmp) / dtold;

	loctmp = dIC_3;
	dIC_3 = (C_3*b3 + IC_2*b11 - (a3 + a11)*IC_3);
	d2IC_3 = (dIC_3 - loctmp) / dtold;

	loctmp = dIC_2;
	dIC_2 = (IC_3*a11 + I_F*b12 + C_2*b3 - (b11 + a3 + a12)*IC_2);
	d2IC_2 = (dIC_2 - loctmp) / dtold;

	loctmp = dI_F;
	dI_F = (IC_2*a12 + IM_1*b4 + O_p*a2 + C_1*b3 - (b12 + a4 + b2 + a3)*I_F);
	d2I_F = (dI_F - loctmp) / dtold;

	loctmp = dIM_1;
	dIM_1 = (I_F*a4 + IM_2*b5 - (b4 + a5)*IM_1);
	d2IM_1 = (dIM_1 - loctmp) / dtold;

	loctmp = dIM_2;
	dIM_2 = (IM_1*a5 - IM_2*b5);
	d2IM_2 = (dIM_2 - loctmp) / dtold;

	C_3 = C_3 + dC_3*dt + d2C_3*dt*dt / 2;
	C_2 = C_2 + dC_2*dt + d2C_2*dt*dt / 2;
	C_1 = C_1 + dC_1*dt + d2C_1*dt*dt / 2;
	O_p = O_p + dO_p*dt + d2O_p*dt*dt / 2;
	IC_3 = IC_3 + dIC_3*dt + d2IC_3*dt*dt / 2;
	IC_2 = IC_2 + dIC_2*dt + d2IC_2*dt*dt / 2;
	I_F = I_F + dI_F*dt + d2I_F*dt*dt / 2;
	IM_1 = IM_1 + dIM_1*dt + d2IM_1*dt*dt / 2;
	IM_2 = IM_2 + dIM_2*dt + d2IM_2*dt*dt / 2;

	P_IM = IM_1 + IM_2;
	P_IC = IC_3 + IC_2;
	P_C = C_3 + C_2 + C_1;
	//P2 = O_p + C_3 + C_2 + C_1 + IC_3 + IC_2 + I_F + IM_1 + IM_2;
	O_p1 = (1.0 - C_3 - C_2 - C_1 - IC_3 - IC_2 - I_F - IM_1 - IM_2);
	Ona = O_p1;
	ina = gna * Ona * (v - ena);

}