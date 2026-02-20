#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>

// using <mscorlib.dll>

// using namespace std;


// ---------------------------------------------------------------------------------------
// General parameters and settings
// ---------------------------------------------------------------------------------------

int StartYear = 1985; //do not change, ever!
int CurrYear;
int BehavCycleCount;
int STDcycleCount;
int ProjectionTerm =  96; //max 136// 18 for 2000 to 2002 cohort // 36 for up to 2020
const int InitPop = 20000; //do not change, ever!
const int MaxPop = 110000; //60000; change to 100000; if run to 2100;
const int MaxCSWs = 500; //200; change to 500; if go to 2100;

int CofactorType = 0; // Nature of STD cofactors in HIV transmission	: 0 = no cofactors,
					  // 1 = multiplicative cofactors; 2 = saturation cofactors
int HPVCofactor = 0; 

int CalculatePAFs = 0; // 1 if you want to calculate PAFs of STDs for HIV transmission,
					   // 0 otherwise
int SyphilisImmunity = 1; // 1 if you want to allow for resistance to syphilis reinfection
						  // when recovered but still seropositive, 0 otherwise
int KeepCurrYearFixed = 0; // 1 if you want to keep the current year fixed at the start year,
						   // 0 otherwise
int NoViralTransm22 = 0; // 1 if you want to assume no transmission of viral STDs in long-
						 // term mutually monogamous relationships, 0 otherwise
int NoViralTransm12 = 0; // 1 if you want to assume no transmission of viral STDs in long-
						 // term relationships if infected partner is monogamous, 0 otherwise
int NoBacterialTransm22 = 1; // 1 if you want to assume no transmission of non-viral STDs in
							 // long-term mutually monogamous relationships, 0 otherwise
int NoBacterialTransm12 = 0; // 1 if you want to assume no transmission of non-viral STDs in
							 // long-term relationships if infected partner is monogamous

int TransitionCalc = 0; // 0 if calculating dependent probs using tradnal actuarial methods;
						// 1 if calculating dependent probs using exact method
int AllowBalancing = 0; // 0 if no explicit balancing of rate of sexual partner acquisition
						// 1 if there is balancing analagous to TSHISA model
int AllowPartnerRateAdj = 0; // 1 if there is individual variation in rates of starting and
							 // ending partnerships
							 // 0 if there is no indiv variation, analogous to TSHISA model
int AllowHIVsuscepAdj = 0; // 1 if there is individual variation in susceptibility to HIV,
						   // independent of the other sources of variation already allowed for
						   // 0 if there is no indiv variation, analogous to TSHISA model
int ConstantInitialHIV = 1; // 1 if the initial # HIV infections is to be the same across 
							// all simulations, 0 if initial HIV prevalence can vary
int SetInitPrev1990 = 1; // 1 if initial HIV prevalence is assigned in 1990

// Change the following indicators to specify which STDs you want to model (1 = include,
// 0 = exclude):
int HIVind = 1;
int HSVind = 0; //Herpes simplex virus
int TPind = 0; //Syphilis
int HDind = 0; 
int NGind = 0; //Gonorrhoeae
int CTind = 0; //Chlamydia trachomatis
int TVind = 0; //Trichomonas vaginalis
int BVind = 0;
int VCind = 0;
int HPVind = 1;

// Change the following indicators to specify for which STDs you want to generate 
// calibration outputs (1 = include, 0 = exclude)
int HIVcalib = 0;
int HSVcalib = 0;
int TPcalib = 0;
int HDcalib = 0;
int NGcalib = 0;
int CTcalib = 0;
int TVcalib = 0;
int BVcalib = 0;
int VCcalib = 0;
int HPVcalib =0;

int CycleS = 48; // Number of sexual behaviour cycles per year
int CycleD = 48; // Number of STD cycles per year (NB: must be integer multiple of CycleS)

// ---------------------------------------------------------------------------------------
// Arrays for random numbers and doing multiple runs, likelihood calculations
// ---------------------------------------------------------------------------------------

double r[InitPop], rprisk[InitPop], rpID[InitPop], rSTI[MaxPop][100];
double r2[MaxPop], revent[MaxPop], rpAge[MaxPop], rpID2[MaxPop], hiv1618[MaxPop],  hivoth[MaxPop], wane[MaxPop], RandAcceptTxV[MaxPop], Rloss[MaxPop], EfficacyD1Rand[MaxPop], EfficacyD2Rand[MaxPop];
const int ParamCombs = 1; // number of input parameter combinations
const int IterationsPerPC = 1; // number of iterations per parameter combination
const int samplesize =ParamCombs*IterationsPerPC; // number of simulations (must = ParamCombs * IterationsPerPC)
int SeedRecord[ParamCombs][2]; // seeds used when FixedUncertainty = 1
int GetSDfromData = 0; // Set to 1 if you want to calculate the standard deviation in the
// likelihood function for the national prevalence data based on the
// model fit to the data rather than the 95% CIs around the survey
// estimates (similar to the approach in Sevcikova et al (2006)).
int VaryParameters = 1; // 1 if parameters vary across simulations
int FixedUncertainty =1; // 0 if parameters vary stochastically, 1 if previously-
// generated parameters are read in from input files, 0 if UseMedians=1

int CurrSim; // Counter for the current simulation
int FixedPeriod = 0; // # years after start of projection in which we use the same seeds
// (across all simulations), to limit variability in HIV simulations
double TotalLogL; // The log of the likelihood for the current simulation (based on HIV and
// STD prevalence data)
double ANCbiasVar = 0.0; // Variance of the ANC bias term (previously set to 0.0225)
double HSRCbiasVar = 0.0; // Variance of the HSRC bias term (previously set to 0.0225)
double FPCweights[16]; // Weights are the proportions of women using modern contraception

// ----------------------------------------------------------------------------------------
// Sexual behaviour parameters and arrays
// ----------------------------------------------------------------------------------------

double HighPropnM, HighPropnF; // % of males and females with propensity for >1 partner
double HighPropn15to24[2]; // % of sexually experienced 15-24 year olds with propensity for
						   // >1 partner, by sex
double AssortativeM, AssortativeF; // Degree of assortative mixing in males and females
double GenderEquality; // Gender equality factor
double AnnNumberClients; // Average annual number of sex acts FSW has with clients
double SexualDebut[16][2]; // Continuous rate at which individuals in high risk group start
						   // sexual activity (by age and sex)
double DebutAdjLow[2]; // Factor by which the rate of sexual debut in the high risk group is
					   // multiplied in order to get the rate of debut in the low risk group
double PartnershipFormation[2][2]; // Average annual number of new partners, by risk (first 
								   // index) & sex (2nd index)
double BasePartnerAcqH[2]; // Rate at which a single 15-19 yr old in high risk group acquires
						   // next sexual partner, by sex (base for other calcs)
double AgeEffectPartners[16][2]; // Factor by which average rate of partnership formation
								 // is multiplied to get age-specific rates (by age & sex)
double GammaMeanST[2]; // Mean of gamma density used to determine AgeEffectPartners
double GammaStdDevST[2]; // Std deviation of gamma density used to determine AgeEffectPartners
double PartnerEffectNew[2][2]; // Factor by which average rate of partnership formation is
							   // multiplied to get prob of acquiring ADDITIONAL partner
							   // (1st index is type of current rel, 2nd index is sex)
double HIVeffectPartners[5]; // Factor by which average rate of partnership formation is
							 // multiplied to get rates by HIV stage
double MarriageIncidence[16][2]; // Continuous rate at which marriages are formed, by age &
								 // sex (note: this is a departure from the Excel model)
double MarriageRate[2][4]; // Rate at which ST relationships become marriages, by risk group
						   // of partner and sex/risk of individual (MH, ML, FH, FL)
double AgeEffectMarriage[16][4]; // Factor by which average rate of marriage is multiplied
								 // to get age-specific rates (by age & sex/risk)
double FSWcontactConstant; // Constant term used in determining rate at which men in high 
						   // risk group have contact with FSWs
double MeanFSWcontacts; // Ave annual # sex acts with FSWs, per sexually experienced male
double AgeEffectFSWcontact[16]; // Factor by which average rate of FSW contact is multiplied
								// to get age-specific rates
double GammaMeanFSW; // Mean of gamma density used to determine AgeEffectFSWcontact 
double GammaStdDevFSW; // Std deviation of gamma density used to determine AgeEffectFSWcontact
double PartnerEffectFSWcontact[5]; // Factor by which average rate of FSW contact is 
								   // multiplied to take account of current partners (0=>
								   // no partner, 1=>1ST, 2=>1LT, 3=>2ST, 4=>1ST & 1LT)
double InitFSWageDbn[16]; // Initial proportion of sex workers at each age
double FSWentry[16]; // Relative rates of entry into FSW group, by age
double FSWexit[16]; // Rates of exit from the FSW group, by age
double HIVeffectFSWentry[5]; // Effect of HIV stage on rate of entry into FSW group
double HIVeffectFSWexit[5]; // Effect of HIV stage on rate of exit from FSW group
double MeanDurSTrel[2][2]; // Mean duration (in years) of short-term relationships, by male
						   // risk group (1st index) & female risk group (2nd index)
double LTseparation[16][2]; // Annual rate (cts) of separation/divorce for LT unions, by
							// age and sex
double AgePrefF[16][16]; // Proportion of male partners in each age group, for women of each
						 // age (1st index: female age, 2nd index: male age)
double AgePrefM[16][16]; // Proportion of female partners in each age group, for men of each
						 // age (1st index: male age, 2nd index: female age)
double FreqSexST[16][2]; // Ave # sex acts per ST relationship, per sexual behaviour cycle
						 // (by age and sex)
double FreqSexLT[16][2]; // Ave # sex acts per LT relationship, per sexual behaviour cycle
						 // (by age and sex)
double BaselineCondomUse; // % of sex acts protected among 15-19 females in ST rels in 1998
double BaselineCondomSvy; // As above, but before applying bias adjustment
double RelEffectCondom[3]; // Effect of partnership type on odds of condom use at baseline
double AgeEffectCondom[3]; // Effect of age on odds of condom use, by partnership type 
double RatioInitialTo1998[3]; // Ratio of inital odds of condom use to odds in 1998
double RatioUltTo1998[3]; // Ratio of ultimate odds of condom use to odds in 1998
double MedianToBehavChange[3]; // Median time (in years since 1985) to condom behav change
double MedianToBehavChange2[3]; 
double ShapeBehavChange[3]; // Weibull shape parameter determining speed of behaviour change
double CondomUseST[16][2]; // Propn of sex acts that are protected in ST rels, by age and sex
double CondomUseLT[16][2]; // Propn of sex acts that are protected in LT rels, by age and sex
double CondomUseFSW; // Propn of sex acts that are protected in FSW-client relationships
double SDpartnerRateAdj; // Std deviation of rates of partnership formation and dissolution, if
						 // AllowPartnerRateAdj = 1
double BehavBiasVar[3][2]; // Sample variance of bias estimates (on logit scale), by type of
						   // behav (1st index) & sex (2nd index)
double CondomScaling; // Parameter to allow for bias in reporting of condoms (1 = no bias)
int MaxNewPartnerInd; // Indicates if the maximum number of new partners has been formed (1 = yes)
double CumAgePrefM[16][16];
double CumAgePrefF[16][16]; 

// ----------------------------------------------------------------------------------------
// Arrays for balancing male and female sexual activity
// ----------------------------------------------------------------------------------------

double DesiredSTpartners[2][2]; // Desired number of new partners, by risk group (1st index)
								// and sex (2nd index)
double DesiredPartnerRiskM[2][2]; // Desired proportion of partners in high and low risk
								  // groups, by male risk group
double DesiredPartnerRiskF[2][2]; // Desired proportion of partners in high and low risk
								  // groups, by female risk group
double AdjSTrateM[2][2]; // Adjustment to rate at which males form partnerships with females
double AdjSTrateF[2][2]; // Adjustment to rate at which females form partnerships with males
double DesiredMarriagesM[2][2]; // Number of new marriages desired by males
double DesiredMarriagesF[2][2]; // Number of new marriages desired by females
double AdjLTrateM[2][2]; // Adjustment to rate at which males marry females
double AdjLTrateF[2][2]; // Adjustment to rate at which females marry males
double ActualPropnLTH[2][2]; // Proportion of long-term partners who are in the high-risk
							 // group, by risk group (1st index) and sex (2nd index)
double ActualPropnSTH[2][2]; // Proportion of short-term partners who are in the high-risk
							 // group, by risk group (1st index) and sex (2nd index)
double DesiredFSWcontacts; // Total numbers of contacts with sex workers (per annum) by men
						   // in the high risk group
double RequiredNewFSW; // Number of women becoming sex workers in current behaviour cycle,
					   // in order to meet excess male demand
int TotCurrFSW; // Total female sex workers at current time

// ----------------------------------------------------------------------------------------
// Demographic parameters and arrays
// ----------------------------------------------------------------------------------------

double HIVnegFert[7]; // Fertility rates in HIV-negative women
double SexuallyExpFert[7]; // Fertility rates in sexually experienced women
double FertilityTable[7][136]; // Fertility rates in HIV-negative women, by age and year
double InfantMort1st6mM[136]; // Prob of death in 1st 6 months of life (males), by year
double ChildMortM[15][136]; // Non-AIDS mortality rates in male children by age and year
double InfantMort1st6mF[136]; // Prob of death in 1st 6 months of life (females), by year
double ChildMortF[15][136]; // Non-AIDS mortality rates in female children by age and year
double NonAIDSmortM[16][136]; // Non-AIDS mortality rates in males by age and year
double NonAIDSmortF[16][136]; // Non-AIDS mortality rates in females by age and year
double StartPop[91][2]; // Numbers of males and females at each individual age (0, ..., 90)
						// as at the start of the projection
double MaleBirthPropn; // Propn of births that are male
int LTP[320][2]; //starting values for numbers of lifetime partners
// ----------------------------------------------------------------------------------------
// STD parameters not defined in the classes below
// ----------------------------------------------------------------------------------------

double HSVsheddingIncrease[5]; // % increase in HSV-2 shedding by HIV stage
double HSVrecurrenceIncrease[5]; // % increase in recurrence rate by HIV stage
double HSVsymptomInfecIncrease; // Multiple by which HSV-2 infectiousness increased
								// when symptomatic
double InfecIncreaseSyndrome[3][2]; // % by which HIV infectiousness increases when
									// experiencing syndrome
double SuscepIncreaseSyndrome[3][2]; // % by which HIV susceptibility increases when
									// experiencing syndrome
double RatioAsympToAveM; // Ratio of asymptomatic HIV transmission prob to average HIV
						 // transmission prob in males
double RatioAsympToAveF; // Ratio of asymptomatic HIV transmission prob to average HIV
						 // transmission prob in females
double MaxHIVprev; // Maximum HIV prevalence in any cohort in current STD cycle
double RelHIVfertility[6]; // Factor by which fertility rate is multiplied in each HIV stage
double PropnInfectedAtBirth; // Propn of children born to HIV+ mothers who are infected at
							 // or before birth
double PropnInfectedAfterBirth; // Propn of children born to HIV+ mothers who are infected 
								// after birth (through breastfeeding)
double MaleRxRate; // Rate at which adults males seek STD treatment
double MaleTeenRxRate; // Rate at which males aged <20 seek STD treatment
double FemRxRate; // Rate at which adults females seek STD treatment
double FemTeenRxRate; // Rate at which females aged <20 seek STD treatment
double FSWRxRate; // Rate at which female sex workers seek STD treatment
double InitMaleRxRate;
double InitMaleTeenRxRate;
double InitFemRxRate;
double InitFemTeenRxRate;
double InitRecurrenceRateM;
double InitRecurrenceRateF;
double PropnTreatedPublicM; // Propn of male STD cases treated in public health sector
double PropnTreatedPublicF; // Propn of female STD cases treated in public health sector
double PropnTreatedPrivateM; // Propn of male STD cases treated in private health sector
double PropnTreatedPrivateF; // Propn of female STD cases treated in private health sector
double PropnPublicUsingSM[136]; // Propn of providers in public sector using syndromic mngt
double PropnPrivateUsingSM[136]; // Propn of providers in private sector using syndromic mngt
double DrugShortage[136]; // % redn in public sector treatment effectiveness due to drug
						 // shortages, by year
double HAARTaccess[136]; // % of people progressing to AIDS who start HAART, by year
double RateARTstart[136][3][2]; // Rate of starting ART in adults, by year, HIV stage and sex
double ARTinterruption[136]; // Annual rate at which adults interrupt ART
double ARTresumption; // Annual rate at which adults resume ART if previously on ART
double PMTCTaccess[136]; // % of pregnant women who are offered PMTCT for HIV, by year
double PropnCiproResistant[136]; // Propn of NG cases that are ciprofloxacin-resistant, by year
double PropnCiproTreated[136]; // Propn of treated NG cases that are given ciprofloxacin
double RxPhaseIn[136]; // Unscaled rates of phase-in for future STD interentions, by year
double InitDrugEffNG; // Initial effectiveness of drugs for treating NG
double InitCorrectRxHSV;
double InitCorrectRxTVM;
double AcceptScreening; // % of pregnant women offered HIV screening who accept
double AcceptNVP; // % of women testing positive who agree to receive nevirapine
double RednNVP; // % reduction in HIV transmission at/before birth if woman receives NVP
double RednFF; // % reduction in HIV transmission after birth if woman receives formula
			   // feeding OR exclusive breasfeeding
double SecondaryRxMult; // Factor by which the rate of treatment seeking is multiplied
						// when experiencing symptoms of secondary syphilis
double SecondaryCureMult; // Factor by which the probability of cure is multiplied
						  // when treated for secondary syphilis
double FSWasympRxRate; // Rate at which female sex workers seek STD treatment when 
					   // asymptomatically infected with an STD
double FSWasympCure; // Prob that treatment for symptomatic STD in FSW cures other 
					 // asymptomatic STDs (expressed as multiple of the prob of cure if the
					 // STD was symptomatic)
double InitFSWasympRxRate;
double InitFSWasympCure; 
double InitANCpropnScreened;
double InitANCpropnTreated;
double InitHIVprevHigh; // % of high risk group initially HIV-positive (assumed to be asymp)
double HlabisaRatio[7][2]; // Ratio of age- and sex-specific initial prev to 15-49 fem prev
double InitHIVtransm[3][2]; // HIV transm probs per sex act at start of epidemic, by
							// nature of rel (1st index) and sex (2nd index)
double RatioUltToInitHIVtransm; // Ratio of ultimate HIV transmission prob to initial HIV
								// transm prob, in the 'no STD cofactor' scenario
double SDsuscepHIVadj; // Coefficient of variation in HIV susceptibility (after controlling 
					   // for other factors in the model)
double CurrPerinatal; // Propn of HIV-exposed kids infected perinatally in current year
double CurrPostnatal; // Propn of HIV-exposed kids infected postnatally in current year

double PropinCIN3[15][60]; // prop of women per age group who have been in CIN3 for 0 to 60 years
double PropVaccinated[136]; //proportion of girls vaccinated in the national programme. Assume constant for now
int CampaignYear[136]; //Five yearly campaign starting in 2025
double PropVaccinatedWHO; //proportion of girls vaccinated in the national programme. Assume constant for now
double PropVaccinatedHIV; 
double VaccEfficacy[13]; //HPV vaccine efficacy for each of 13 types (there is evidence of cross-protection)
double VaccEfficacyNONA[13]; //HPV vaccine efficacy for each of 13 types (nonavalent vaccine)
double TxVEfficacyD1[13]; // HPV therapeutic vaccine efficacy for each of 13 types HPV infection and CIN 1
double TxVEfficacyD2[13];
double TxVEfficacyD1CIN[13]; // HPV therapeutic vaccine efficacy for each of 13 types for CIN 2/3
double TxVEfficacyD2CIN[13];
// ----------------------------------------------------------------------------------------
// HPV/CC parameters not defined in the classes below
// ----------------------------------------------------------------------------------------
//int CCcalib=1; //if 0, all parameters are drawn from priors, if 1, medians of CIN parameters are read
/*int UpdateStart; //if 1, get new starting conditions for HPV stages. Use medians from posteriors. 
					 //FixedUncertainty should be 0. HIVind should be 0. RoutineScreening and HPVvacc should be 0.
					 //CondomUsage should be se to zero, except for FSW!
int CreateCohort; //if 1, extract all the information to do statistical analyses on HIV/HPV interaction. 
                      //PCs from posterior distributions are read in. FixedUncertainty should be 1, unless UseMedians=1  
int UseMedians;
int OneType =0;  //if 1, only simulate one HPV type specified below, using posterior distributions, ONLY IF FixedUncertainty==1
int WhichType = 0; //0=16 ; 1=18 ; 2=31 ; 3=33 ; 4=35 ; 5=39 ; 6=45 ; 7=51 ; 8=52 ; 9=56 ; 10=58 ; 11=59 ; 12=68
int targets = 0; //if one, need to specify WhichType
int GetMacD = 0;

//Strategy indicators
int ImplementYR; //=2010;
int HPVvacc; //=1; //1 if National HPV vaccination program is active since 2014
int BOYSvacc; //=0; //1 if boys are also included in the National vaccination programme
int RoutineScreening; //=1; //switch off screening =0

int PerfectSchedule=0; //If SA's screening schedule is followed to a T
int HPVDNA=0;
int HPVDNAThermal=0;
double PropThermal=1.0;
int HPVGenotyping=0;
int PapTRIAGE=0;
int Portal3=0; //If 0:yearly follow-up ; If 1:test-and-treat ; NOT IMPLEMENTED YET If 2: PapTRIAGE
//WHO Project indicators
int WHOscenario = 0;  //in WHO scenario, vacc starts at age 10 (as opposed to 9 in SA) 
int WHOvacc =0; //1 if Nonavalent vacc is rolled out 
int WHOScreening=0; 
int S3S4=0;
int S5S6=0;
int S7S11=0;
int CatchUpVaccHIV=0;
int VaccineWane=0;
int VaccDur=20;*/
int CCcalib;   // if 0, all parameters are drawn from priors, if 1, medians of CIN parameters are read
int UpdateStart;// = 0; //if 1, get new starting conditions for HPV stages. Use medians from posteriors. 
					 //FixedUncertainty should be 0. HIVind should be 0. RoutineScreening and HPVvacc should be 0.
					 //CondomUsage should be se to zero, except for FSW!
int CreateCohort;// =0; //if 1, extract all the information to do statistical analyses on HIV/HPV interaction. 
                      //PCs from posterior distributions are read in. FixedUncertainty should be 1, unless UseMedians=1  
int UseMedians;// =0;
int OneType;// =0;  //if 1, only simulate one HPV type specified below, using posterior distributions, ONLY IF FixedUncertainty==1
int WhichType;// = 0; //0=16 ; 1=18 ; 2=31 ; 3=33 ; 4=35 ; 5=39 ; 6=45 ; 7=51 ; 8=52 ; 9=56 ; 10=58 ; 11=59 ; 12=68
int targets;// = 0; //if one, need to specify WhichType
int GetMacD;// = 0;

//Strategy indicators
int ImplementYR;//=2010;
int HPVvacc;//=1; //1 if National HPV vaccination program is active since 2014
int BOYSvacc;//=0; //1 if boys are also included in the National vaccination programme
int RoutineScreening;//=1; //switch off screening =0

int PerfectSchedule;//=0; //If SA's screening schedule is followed to a T
int HPVDNA;//=1;
int MnDCampaign;//=1; //HPVDNA should also be 1 if this is 1
int HPVDNAThermal;//=1; //HPVDNA should also be 1 if this is 1
double PropThermal;//=1.0;
int HPVGenotyping;//=0;
int PapTRIAGE;//=0;
int Portal3;//=0; //If 0:yearly follow-up ; If 1:test-and-treat ; NOT IMPLEMENTED YET If 2: PapTRIAGE
//WHO Project indicators
int WHOscenario;// = 0;  //in WHO scenario, vacc starts at age 10 (as opposed to 9 in SA) 
int WHOvacc;// =0; //1 if Nonavalent vacc is rolled out 
int WHOScreening;//=0; 
int S3S4;//=0;
int S5S6;//=0;
int S7S11;//=0;
int CatchUpVaccHIV;
int CatchUpVacc;//=0;
int CatchUpAgeMIN;//=0;
int CatchUpAgeMAX;//=0;
double CatchUpCoverage;//=0;
int VaccineWane;//=0;
int VaccDur;//=20;
int AdministerMassTxV;
int MassTxVAgeMIN;
int MassTxVAgeMAX;
int MassTxVDur;
int MassTxVWane;
int AdministerMassTxVtoART;
int MassTxVtoARTAgeMIN;
int MassTxVtoARTAgeMAX;
double ReductionFactorHIV;
double ReductionFactorHIVART;
int TxVviaScreeningAlgorithm;
int TxVtoHPVPos;
int TxVtoHPVPosLTFU;
//double ScreenReason[8];
double ScreenReason[8][136]; //Reason for screen by age + HIV status, over time
//int TotScreens[12];
double ScreenProb[8][136]; //Probability of entering screening  by age (4 categories) + HIV/ART status (neg/noART, ART), over time
double AttendColposcopy[3][136];
double PapAdequacy[136];
double WHOcoverage[8][136];
//double ModelCoverage[12][11];

// ---------------------------------------------------------------------------------------
// Classes for STD prevalence data and likelihoods
// ---------------------------------------------------------------------------------------

class PrevalenceData
{
public:
	PrevalenceData();

	double LogL;
	double LogL0;
	double LogL1;
	double LogL2;
	double LogL3;
	int Observations;
	int SexInd; // 0 = males, 1 = females
};

class NationalData : public PrevalenceData
{
public:
	NationalData();

	// Note that in defining the arrays below we are assuming that there would not be more
	// than 100 data points to calibrate to. If for some reason there are (i.e. Observations
	// >100), then one should change the dimensions of the arrays below.

	int StudyYear[100];
	double StudyPrev[100];
	double PrevSE[100];
	int AgeStart[100];
	double ExpSe[100];
	double ExpSp[100];
	double ModelPrev[100];
	double BiasMult[16];
	double ModelVarEst;

	void CalcModelVar();
	void CalcLogL();
};

class AntenatalN : public NationalData
{
public:
	AntenatalN();
};

class HouseholdN : public NationalData
{
public:
	HouseholdN();
	int AgeStart[28];
	int AgeEnd[28];
	int ExclVirgins[28];
	int HIVprevInd[28]; // 1 if HIV prevalence was measured in the study, 0 otherwise
	double HIVprev[28];
	int Cyt_cat[28];
};

class SentinelData : public PrevalenceData
{
public:
	SentinelData();

	// Note that in defining the arrays below we are assuming that there would not be more
	// than 25 data points to calibrate to. If for some reason there are (i.e. Observations
	// >25), then one should change the dimensions of the arrays below.

	int StudyYear[300];
	int StudyN[300];
	double StudyPrev[300];
	double StudyPos[300];
	double ExpSe[300];
	double ExpSp[300];
	double VarSe[300];
	double VarSp[300];
	int HIVprevInd[300]; // 1 if HIV prevalence was measured in the study, 0 otherwise
	double HIVprev[300];
	int Cyt_cat[300];
	int Perfect[300];
	int ARTind[300];
	double ModelPrev[300];
	double ModelNum[300];
	double ModelDenom[300];
	double CurrentModelPrev;
	double BiasMult;
	double VarStudyEffect; // Sigma(b) squared
	double VarStudyEffectAB;
	double VarStudyEffectH;
	void CalcLogL();
};

class Household : public SentinelData
{
public:
	Household();

	int AgeStart[300];
	int AgeEnd[300];
	int ExclVirgins[300];
	
};

class NonHousehold : public SentinelData
{
public:
	NonHousehold();
};

class ANC : public NonHousehold
{
public:
	ANC();
};

class FPC : public NonHousehold
{
public:
	FPC();
};

class GUD : public NonHousehold
{
public:
	GUD();
};

class CSW : public NonHousehold
{
public:
	CSW();
};

class ART : public NonHousehold //HPV - Cari//
{
public:
	ART();

	int AgeStart[300];
	int AgeEnd[300];
};

// ---------------------------------------------------------------------------------------
// Classes for outputs
// ---------------------------------------------------------------------------------------

class PostOutputArray
{
	// We use this to record outputs for different years, from multiple simulations. 

public:
	PostOutputArray(int n);

	int columns;
	double out[samplesize][136]; // None of the arrays require > 136 columns. 

	void RecordSample(const char* filout);
};

class PostOutputArray2
{
	// Same as PostOutputArray, but with dimension ParamCombs instead of samplesize. 

public:
	PostOutputArray2(int n);

	int columns;
	double out[samplesize][136]; // None of the arrays require > 136 columns. 

	void RecordSample(const char* filout);
};

class PostOutputArray3
{
	// We use this to record outputs for different years, for different ages (combining numbers from diff simulations). 

public:
	PostOutputArray3();

	double out[54][136]; // None of the arrays require > 136 columns. 

	void RecordSample(const char* filout, int type);
};
class PostOutputArray4
{
	// We use this to record outputs for different years, for different ages (combining numbers from diff simulations). 

public:
	PostOutputArray4();

	double out[samplesize][136]; // None of the arrays require > 136 columns. 

	void RecordSample(const char* filout, int type);
};


// ---------------------------------------------------------------------------------------
// Classes for STD parameters and probabilities of transition between STD disease states
// ---------------------------------------------------------------------------------------

class STDtransition
{
public:
	STDtransition();

	// NB: The CondomEff and SuscepIncrease members actually belong to the
	// opposite sex. This is potentially confusing, but we have done things this way 
	// to be consistent with the deterministic model.

	int nStates;
	int SexInd; // 0 = male, 1 = female
	double AveDuration[7]; // Average number of weeks spent in state if untreated
	double CondomEff; // Probability that condom prevents transmission
	double SuscepIncrease[16]; // Multiple by which susceptibility to STD increases, by age
	double HIVinfecIncrease[6]; // % by which HIV infectiousness increases, by HIV/STD stage
	double ARTinfectiousness[136]; // decrease of HIV infectiousness of people on ART, by year
	double RelTransmCSW; // Ratio of M->F transm prob per sex act in CSW-client relationships 
						 // to that in non-spousal relationships
	double RelTransmLT; // Ratio of transm prob in spousal relationships to that in
						// non-spousal relationships

	// Objects used for calibration purposes
	ANC ANClogL;
	FPC FPClogL;
	GUD GUDlogL;
	CSW CSWlogL;
	ART ONARTlogL;
	ART NOARTlogL;
	Household HouseholdLogL;
	AntenatalN AntenatalNlogL;
	HouseholdN HouseholdNlogL;
	double CSWprevUnsmoothed[30];

	void ReadPrevData(const char* input);
	void GetCSWprev();
	void SetVarStudyEffect(double Variance);
};

class HIVtransition: public STDtransition
{
public:
	HIVtransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW,
		int ObsHH, int ObsANCN, int ObsHHN, int ObsNOART, int ObsONART);

	double TransmProb[7]; // Transmission probability per sex act (by relationship type)
	double From1to2; // Prob of progressing from acute HIV to asymp per STD cycle
	double From2to3; // Prob of progressing from asymp to pre-AIDS symptoms, per STD cycle
	double From2to5; // Prob of progressing from asymp to ART, per STD cycle
	double From3to4; // Prob of progressing from pre-AIDS to untreated AIDS per STD cycle
	double From3to5; // Prob of progressing from pre-AIDS to ART per STD cycle
	double From4to5; // Prob of starting ART from untreated AIDS per STD cycle
	double From5to6; // Prob of interrupting ART per STD cycle
	double From6to5; // Prob of resuming ART after an interruption, per STD cycle
	double From4toDead; // Prob of dying from AIDS if not receiving HAART
	double From5toDead; // Prob of dying from AIDS while receiving HAART
	double From6toDead; // Prob of dying from AIDS while receiving HAART
	double From5toDeadEarly; // Prob of dying from AIDS while receiving HAART
	double From5toDeadLateShort; // Prob of dying from AIDS while receiving HAART
	double From5toDeadLateLong; // Prob of dying from AIDS while receiving HAART

	void CalcTransitionProbs();
	void GetPrev();
	double GetTransmProb(int ID);
	void GetNewStage(int ID, double p);
};

class NonHIVtransition: public STDtransition
{
public:
	NonHIVtransition();

	double CorrectRxPreSM; // % of cases correctly treated before syndromic management 
	double CorrectRxWithSM; // % of cases correctly treated under syndromic management
	double DrugEff; // Prob of complete cure if effective drugs are prescribed
	double TradnalEff; // Prob of cure if treated by a traditional healer
	double ProbCompleteCure; // Prob that individual seeking treatment is cured
	double HIVsuscepIncrease[6]; // % by which HIV susceptibility increases, by STD stage
	double PropnByStage[320][8]; // 20 behav states (max) x 16 age groups = 320 rows
								 // Max of 8 disease states. This matrix determines initial
								 // % of individuals in different STD stages.

	void CalcProbCure();
	int GetSTDstage(int offset, double r);
};

class TradnalSTDtransition: public NonHIVtransition
{
public:
	TradnalSTDtransition();

	int ImmuneState; // The index of the state representing individuals immune to reinfection
					 // (value is 0 if there is no immune state)
	double TransmProb; // Probability of transmission per act of sex
	double TransmProbSW; // Probability of transmission per act of sex (note that this is
						 // only used for client-to-sex worker transmission)
};

class SyphilisTransition: public TradnalSTDtransition
{
public:
	SyphilisTransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH,
		int ObsANCN, int ObsHHN, int ObsNOART, int ObsONART);

	double ANCpropnScreened; // % of women attending ANCs who are screened for syphilis
	double ANCpropnTreated; // % of women testing positive who receive treatment
	double ANCpropnCured; // % of women attending ANCs who get cured
	double ProbANCcured[16]; // Prob that a woman with syphilis attends an ANC and is cured,
							// per STD cycle (by age: 10-14,15-19, ...)
	double PropnSuscepAfterRx; // Proportion of individuals who are susceptible to reinfection
							   // after cure of primary syphilis

	double From1to2; // Prob that incubating syphilis becomes primary syphilis per STD cycle
	double From2to0; // Prob that primary syphilis is cured and RPR- per STD cycle
	double From2to0T; // Prob that primary syphilis is cured and RPR- per STD cycle, age <20
	double From2to0C; // Prob that primary syphilis is cured and RPR- per STD cycle, in CSWs
	double From2to3; // Prob that primary syphilis becomes secondary syphilis per STD cycle
	double From2to3T; // Prob that primary syphilis becomes secondary syphilis, age <20
	double From2to3C; // Prob that primary syphilis becomes secondary syphilis, in CSWs
	double From2to5; // Prob that primary syphilis is cured but still RPR+ per STD cycle
	double From2to5T; // Prob that primary syphilis is cured but still RPR+, age <20
	double From2to5C; // Prob that primary syphilis is cured but still RPR+, in CSWs
	double From3to4; // Prob that secondary syphilis becomes latent per STD cycle
	double From3to4T; // Prob that secondary syphilis becomes latent per STD cycle, age <20
	double From3to4C; // Prob that secondary syphilis becomes latent per STD cycle, in CSWs
	double From3to5; // Prob that secondary syphilis is cured per STD cycle
	double From3to5T; // Prob that secondary syphilis is cured per STD cycle, age <20
	double From3to5C; // Prob that secondary syphilis is cured per STD cycle, in CSWs
	double From4to6; // Prob that latent syphilis resolves per STD cycle
	double From4to6C; // Prob that latent syphilis resolves per STD cycle, in CSWs
	double From5to0; // Prob of seroreversion if cured in early syphilis, per STD cycle
	double From6to0; // Prob of seroreversion if resolved in latent syphilis, per STD cycle

	void CalcTransitionProbs();
	double GetTransmProb(int ID);
	void GetNewStage(int ID, double p);
};

class HerpesTransition: public TradnalSTDtransition
{
public:
	HerpesTransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH,
		int ObsANCN, int ObsHHN, int ObsNOART, int ObsONART);

	double RecurrenceRate; // Rate at which symptomatic recurrences occur
	double SymptomaticPropn; // % of people who develop primary ulcer
	double From1to2; // Prob that primary ulcer resolves per STD cycle
	double From1to2T; // Prob that primary ulcer resolves per STD cycle, age <20
	double From1to2C; // Prob that primary ulcer resolves per STD cycle, in CSWs
	double From2to3[6]; // Prob of symptomatic recurrence per STD cycle, by HIV stage
	double From3to2; // Prob of resolution of recurrent ulcer per STD cycle
	double From3to2T; // Prob of resolution of recurrent ulcer per STD cycle, age <20
	double From3to2C; // Prob of resolution of recurrent ulcer per STD cycle, in CSWs
	double From2to4; // Prob that transiently asymp indiv becomes permanently asymp

	void CalcTransitionProbs();
	double GetTransmProb(int ID);
	void GetNewStage(int ID, double p);
};

class HPVTransition : public TradnalSTDtransition
{
public:
	HPVTransition();
	//HPVTransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH,
	//int ObsANCN, int ObsHHN, int ObsART);
	void HPVObs(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH, int ObsANCN, int ObsHHN, int ObsNOART, int ObsONART);

	//fixed parameters
	double PropnScreenedCIN[136]; //Proportion screened using Pap, by year
	double PropnAdequateSmear; //Proportion of smears containing endo- and ectocervical cells CHANGE TO BY YEAR
	double PropnDxCIN; //for now, prop of CIN2/3 picked up by PAP. Should be matrix with proportions diagnosed correctly, by cervical disease stage (to determine proportions that need Tx)
	double PropnTxColpos; //Proportion diagnosed who present at colposcopy clinic
	double PropnTxSuccess; //Proportion with successfull treatment
	double TxHPVclear; //Proportion that will clear HPV after treatment
	
	//varying parameters
	double prop1pro; //Proportion in HPV+ who will progress to CIN1
	double prog_fromHPV; //Progression rate from  HPV to CIN1
	double prog_fromCIN1; //Progression rate from  CIN1 to CIN2
	double prog_fromCIN2; //Progression rate from  CIN2 to CIN3
	double reg_fromHPV; //Regression rate from  HPV 
	double reg_fromCIN1; //Regression rate from  CIN1 to Normal
	double reg_fromCIN2; //Regression rate from  CIN2 to CIN1
	double CIN2reg30;
	double CIN2prog30;
	double CIN2prog50;
	double reg_fromCIN3; //Regression rate from  CIN2 to CIN1
	double prop_reg_cl; //proportion that clears HPV infection when regressing from CIN1
	double propLatent; //Proportion who will regress from HPV+ to latent state
	double pImmune; //parameter for partial immunity
	double rr;
	
	//multipliers for reactivation
	double latentHIVreact;
	double lateHIVreact;
	//multipliers for clearance
	double AcuteLateHIVclear;
	double LatentHIVclear; 
	//multipliers for HPV duration in CC calib 
	double durMult;
	//multipliers for progression of reactivated infection 
	double ReactPro;
	//multipliers for HIV
	double CIN2_HIV;
	double CIN1_HIV;
	double prog_ART;
	double CIN2_ART;
	double CIN1_ART;
	double reg_HIV;
	double reg_ART;
	double CIN3HIV;
	double CIN3_50;
	double CIN3HIV_ART;
	double CIN3shape; 
	double CCprogART;
	double CCprogHIV;
	//stages of CC
	double prog_stageI;
	double prog_stageII;
	double prog_stageIII;
	double diag_stageI;
	double diag_stageII;
	double diag_stageIII;
	double diag_stageIV;

	double From1to2ind; // Independent prob of going from HPV+ to CIN1 per STD cycle
	double From1to2dep; // Dependent prob of going from HPV+ to CIN1 per STD cycle
	double From1to6ind; // Independent prob of going from HPV+ to latent per STD cycle
	double From1to6dep; // Dependent prob of going from HPV+ to latent per STD cycle
	double From1to7ind; // Independent prob of going from HPV+ to immune per STD cycle
	double From1to7dep; // Dependent prob of going from HPV+ to immune per STD cycle
	double From1to2indAcute; // Independent prob of going from HPV+ to CIN1 per STD cycle
	double From1to2depAcute; // Dependent prob of going from HPV+ to CIN1 per STD cycle
	double From1to6indAcute; // Independent prob of going from HPV+ to latent per STD cycle
	double From1to6depAcute; // Dependent prob of going from HPV+ to latent per STD cycle
	double From1to7indAcute; // Independent prob of going from HPV+ to immune per STD cycle
	double From1to7depAcute; // Dependent prob of going from HPV+ to immune per STD cycle
	double From1to2indLatent; // Independent prob of going from HPV+ to CIN1 per STD cycle
	double From1to2depLatent; // Dependent prob of going from HPV+ to CIN1 per STD cycle
	double From1to6indLatent; // Independent prob of going from HPV+ to latent per STD cycle
	double From1to6depLatent; // Dependent prob of going from HPV+ to latent per STD cycle
	double From1to7indLatent; // Independent prob of going from HPV+ to immune per STD cycle
	double From1to7depLatent; // Dependent prob of going from HPV+ to immune per STD cycle

	//different rates if reactivated latent infection
	double From1to2depR; // Dependent prob of going from HPV+ to CIN1 per STD cycle
	double From1to6depR; // Dependent prob of going from HPV+ to latent per STD cycle
	double From1to7depR; // Dependent prob of going from HPV+ to immune per STD cycle
	double From1to2depAcuteR; // Dependent prob of going from HPV+ to CIN1 per STD cycle
	double From1to6depAcuteR; // Dependent prob of going from HPV+ to latent per STD cycle
	double From1to7depAcuteR; // Dependent prob of going from HPV+ to immune per STD cycle
	double From1to2depLatentR; // Dependent prob of going from HPV+ to CIN1 per STD cycle
	double From1to6depLatentR; // Dependent prob of going from HPV+ to latent per STD cycle
	double From1to7depLatentR; // Dependent prob of going from HPV+ to immune per STD cycle
	
	//Progression/regression from CIN1 if aged<30
	double From2to1ind; // Prob of regressing from CIN1 to HPV+ per STD cycle
	double From2to6ind; // Prob of regressing from CIN1 to latent per STD cycle
	double From2to7ind; // Prob of regressing from CIN1 to immune per STD cycle
	double From2to3ind; // Prob of progressing from CIN1 to CIN2 per STD cycle
	double From2to1dep; // Prob of regressing from CIN1 to HPV+ per STD cycle
	double From2to6dep; // Prob of regressing from CIN1 to latent per STD cycle
	double From2to7dep; // Prob of regressing from CIN1 to immune per STD cycle
	double From2to3dep; // Prob of progressing from CIN1 to CIN2 per STD cycle
	
	double From2to1depLate; // Prob of regressing from CIN1 to HPV+ per STD cycle
	double From2to3depLate; // Prob of progressing from CIN1 to CIN2 per STD cycle
	double From2to6depLate; // Prob of regressing from CIN1 to latent per STD cycle
	double From2to7depLate; // Prob of regressing from CIN1 to immune per STD cycle
	
	double From2to1depLatent; // Prob of regressing from CIN1 to HPV+ per STD cycle
	double From2to3depLatent; // Prob of progressing from CIN1 to CIN2 per STD cycle
	double From2to6depLatent; // Prob of regressing from CIN1 to latent per STD cycle
	double From2to7depLatent; // Prob of regressing from CIN1 to immune per STD cycle
	
	//Progression/regression from CIN1 if aged>=30
	double From2to1dep30; // Prob of regressing from CIN1 to HPV+ per STD cycle
	double From2to6dep30; // Prob of regressing from CIN1 to latent per STD cycle
	double From2to7dep30; // Prob of regressing from CIN1 to immune per STD cycle
	double From2to3dep30; // Prob of progressing from CIN1 to CIN2 per STD cycle
	
	double From2to1depLate30; // Prob of regressing from CIN1 to HPV+ per STD cycle
	double From2to3depLate30; // Prob of progressing from CIN1 to CIN2 per STD cycle
	double From2to6depLate30; // Prob of regressing from CIN1 to latent per STD cycle
	double From2to7depLate30; // Prob of regressing from CIN1 to immune per STD cycle
	
	double From2to1depLatent30; // Prob of regressing from CIN1 to HPV+ per STD cycle
	double From2to3depLatent30; // Prob of progressing from CIN1 to CIN2 per STD cycle
	double From2to6depLatent30; // Prob of regressing from CIN1 to latent per STD cycle
	double From2to7depLatent30; // Prob of regressing from CIN1 to immune per STD cycle
	//Progression/regression from CIN2 if aged<30
	double From3to2ind; // Prob of regressing from CIN2 to CIN1 per STD cycle
	double From3to4ind; // Prob of progressing from CIN2 to CIN3 per STD cycle
	double From3to2dep; // Prob of regressing from CIN2 to CIN1 per STD cycle
	double From3to4dep; // Prob of progressing from CIN2 to CIN3 per STD cycle
	double From3to2depLate; // Prob of regressing from CIN2 to CIN1 per STD cycle
	double From3to4depLate; // Prob of progressing from CIN2 to CIN3 per STD cycle
	double From3to2depLatent; // Prob of regressing from CIN2 to CIN1 per STD cycle
	double From3to4depLatent; // Prob of progressing from CIN2 to CIN3 per STD cycle
	//Progression/regression from CIN2 if aged>30
	double From3to2dep30; // Prob of regressing from CIN2 to CIN1 per STD cycle
	double From3to4dep30; // Prob of progressing from CIN2 to CIN3 per STD cycle
	double From3to2depLate30; // Prob of regressing from CIN2 to CIN1 per STD cycle
	double From3to4depLate30; // Prob of progressing from CIN2 to CIN3 per STD cycle
	double From3to2depLatent30; // Prob of regressing from CIN2 to CIN1 per STD cycle
	double From3to4depLatent30; // Prob of progressing from CIN2 to CIN3 per STD cycle
	//Progression/regression from CIN2 if aged>50
	double From3to2dep50; // Prob of regressing from CIN2 to CIN1 per STD cycle
	double From3to4dep50; // Prob of progressing from CIN2 to CIN3 per STD cycle
	double From3to2depLate50; // Prob of regressing from CIN2 to CIN1 per STD cycle
	double From3to4depLate50; // Prob of progressing from CIN2 to CIN3 per STD cycle
	double From3to2depLatent50; // Prob of regressing from CIN2 to CIN1 per STD cycle
	double From3to4depLatent50; // Prob of progressing from CIN2 to CIN3 per STD cycle
	
	double From4to3; 
	double From4to5; // Prob of progressing from CIN3 to CC per STD cycle
	double From4to5HIV;
	double From5toDead;
	
	double From5to8ind;
	double From5to11ind;
	double From5to8dep;
	double From5to11dep;
	double From5to8depHIV;
	double From5to11depHIV;
	double From5to8depART;
	double From5to11depART;

	double From8to9ind;
	double From8to12ind;
	double From8to9dep;
	double From8to12dep;
	double From8to9depHIV;
	double From8to12depHIV;
	double From8to9depART;
	double From8to12depART;
	
	double From9to10ind;
	double From9to13ind;
	double From9to10dep;
	double From9to13dep;
	double From9to10depHIV;
	double From9to13depHIV;
	double From9to10depART;
	double From9to13depART;

	double From10toDead;
	double From10to14;
	
	double From6to1; // Prob of going from latent to HPV+ per STD cycle during 
	double From6to1Latent; // Prob of going from latent to HPV+ per STD cycle
	double From6to1Late; // Prob of going from latent to HPV+ per STD cycle
	
	double From7to0; // Prob of going from immune to HPV- per STD cycle
	
	void CalcTransitionProbsM();
	void CalcTransitionProbsF();
	double GetTransmProb(int ID, int type);
	void GetNewStageM(int ID, double p, int type);
	void GetNewStageF(int ID, double p, int type);
};

class HPVtransitionCC : public TradnalSTDtransition
{
public:
	
	HPVtransitionCC(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH,
	int ObsANCN, int ObsHHN, int ObsNOART, int ObsONART);
	void HPVObs(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH, int ObsANCN, int ObsHHN, int ObsNOART, int ObsONART);

	int TotPop30to65neg[136]; 
	int TotPop30to65pos[136];
	int TotPop30to65art[136];
	int  TotPop18to60noart[136];
	int  TotPop18to65art[136];
	
    //Other
    int TotAB30to65neg[136]; 
	int TotAB30to65pos[136]; 
	int TotAB30to65art[136]; 
	int TotAB18to60noart[136];
	int TotAB18to65art[136];
	int TotHSIL30to65neg[136];
	int TotHSIL30to65pos[136];
	int TotHSIL30to65art[136];
	int TotHSIL18to60noart[136];
	int TotHSIL18to65art[136];
	int TotHSILNEG[136];
	int TotHSILPOS[136];
	int TotHSILART[136];
	int TotCCNEG[136];
	int TotCCPOS[136];
	int TotCCART[136];

	//WHO
	int TotPop15upnegF[136]; 
	int TotPop15upposF[136];
	int TotPop15upartF[136];

	int TotPop15to65negF[136]; 
	int TotPop15to65posF[136];
	int TotPop15to65artF[136];
	int TotPop15to65negM[136]; 
	int TotPop15to65posM[136];
	int TotPop15to65artM[136];

	int TotHPV15to65negF[136]; 
	int TotHPV15to65posF[136];
	int TotHPV15to65artF[136];
	int TotHPV15to65negM[136]; 
	int TotHPV15to65posM[136];
	int TotHPV15to65artM[136];

	void GetCurrHPVstage();

};

class OtherSTDtransition: public TradnalSTDtransition
{
public:
	OtherSTDtransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH,
		int ObsANCN, int ObsHHN, int ObsNOART, int ObsONART);

	double SymptomaticPropn; // % of individuals who develop symptoms
	double PropnImmuneAfterRx; // % of individuals who are immune to reinfection after 
							   // successful treatment
	double PropnImmuneAfterSR; // % of individuals who are immune to reinfection after 
							   // spontaneous resolution

	double From1to0; // Prob of symptomatic indiv reverting to susceptible per STD cycle
	double From1to0T; // Prob of symptomatic indiv reverting to susceptible per STD cycle, 
					  // age <20
	double From1to0C; // Prob of symptomatic indiv reverting to susceptible per STD cycle, 
					  // in CSWs
	double From1to3; // Prob of symptomatic indiv becoming immune per STD cycle
	double From1to3T; 
	double From1to3C; 
	double From2to0; // Prob of asymptomatic indiv reverting to susceptible per STD cycle
	double From2to0C; // Prob of asymptomatic indiv reverting to susceptible per STD cycle, 
					  // in CSWs
	double From2to3; // Prob of asymptomatic indiv becoming immune per STD cycle
	double From2to3C; 
	double From3to0; // Prob of immune individual reverting to susceptible per STD cycle

	void CalcTransitionProbs();
	double GetTransmProb(int ID, int STD); // STD 1 for NG, 2 for CT, 3 for TV and 4 for HD 
	void GetNewStage(int ID, double p, int STD); // Same STD codes as before
};

class NonSTDtransition: public NonHIVtransition
{
public:
	NonSTDtransition();

	double DrugPartialEff; // Prob of partial cure if effective drugs are prescribed
	double ProbPartialCure; // Prob that individual seeking treatment is partially cured

	void CalcProbPartialCure();
};

class BVtransition: public NonSTDtransition
{
public:
	BVtransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH, int ObsANCN,
		int ObsHHN, int ObsNOART, int ObsONART);

	double SymptomaticPropn; // Propn of women developing BV who experience symptoms
	double Incidence1; // Weekly incidence of BV in women with intermediate flora, with
					   // one current partner
	double IncidenceMultTwoPartners; // Multiple by which incidence is increased in women
									 // with more than one partner
	double IncidenceMultNoPartners; // Multiple by which incidence is decreased in women
									// with no partners
	double CtsTransition[4][4]; // Continuous rates of transition between states
	double From1to2; // Prob of going from normal vaginal flora to intermediate per STD cycle
	double From2to1ind; // Independent prob of reverting to normal flora per STD cycle
	double From3to1; // Prob of returning to normal flora if BV is currently symptomatic
	double From3to1T; // Prob of returning to normal flora, age <20
	double From3to1C; // Prob of returning to normal flora, in CSWs
	double From3to2; // Prob of returning to intermediate flora if BV is currently symptomatic
	double From3to2T; // Prob of returning to intermediate flora, age <20
	double From3to2C; // Prob of returning to intermediate flora, in CSWs
	double From4to1; // Prob of returning to normal flora if BV is currently asymp
	double From4to1C; // Prob of returning to normal flora if BV is currently asymp, in CSWs
	double From4to2; // Prob of returning to intermediate flora if BV is currently asymp
	double From4to2C; // Prob of returning to intermediate flora if BV is asymp, in CSWs
	double From2to3ind[3]; // Independent prob of developing symptomatic BV per STD cycle, if 
						   // currently intermediate, by # current partners
	double From2to4ind[3]; // Independent prob of developing asymp BV per STD cycle, if 
						   // currently intermediate, by # current partners
	double From2to1dep[3]; // Dependent prob of reverting to normal flora per STD cycle, if 
						   // currently intermediate, by # current partners
	double From2to3dep[3]; // Dependent prob of developing symptomatic BV per STD cycle, if 
						   // currently intermediate, by # current partners
	double From2to4dep[3]; // Dependent prob of developing asymp BV per STD cycle, if 
						   // currently intermediate, by # current partners

	void CalcTransitionProbs();
	void GetNewStage(int ID, double p);
};

class VCtransition: public NonSTDtransition
{
public:
	VCtransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH, int ObsANCN,
		int ObsHHN, int ObsNOART, int ObsONART);

	double RecurrenceRate; // Rate at which symptoms develop in women with yeast colonization
	double Incidence; // Rate at which asymptomatic yeast colonization occurs
	double IncidenceIncrease[5]; // % increase in incidence of yeasts by HIV stage
	double From1to2; // Prob that asymp colonization becomes symptomatic per STD cycle
	double From1to2C; // Prob that asymp colonization becomes symptomatic, in CSWs
	double From1to0; // Prob that asymp colonization resolves per STD cycle
	double From1to0C; // Prob that asymp colonization resolves per STD cycle, in CSWs
	double From2to1; // Prob that symptomatic infection becomes asymp per STD cycle
	double From2to1T; // Prob that symptomatic infection becomes asymp per STD cycle, age <20
	double From2to1C; // Prob that symptomatic infection becomes asymp per STD cycle, in CSWs
	double From2to0; // Prob that symptomatic infection is cured per STD cycle
	double From2to0T; // Prob that symptomatic infection is cured per STD cycle, age <20
	double From2to0C; // Prob that symptomatic infection is cured per STD cycle, in CSWs
	double From0to1[7][6]; // Prob that uninfected woman gets asymptomatically colonized per
						   // STD cycle (by age and HIV stage)
	
	void CalcTransitionProbs();
	void GetNewStage(int ID, double p);
};

class PaedHIV
{
public:
	PaedHIV();

	double PreAIDSmedian;
	double PreAIDSshape;
	double MeanAIDSsurvival;
	double AIDSprob1st6m; // Independent prob of progressing to AIDS in 1st 6 months of life
	double AveYrsOnART; // Average survival after ART initiation (ignoring non-HIV mortality)
	double NonAIDSmort[15];
	double MortProb1st6m; // Independent prob of non-AIDS mortality in 1st 6 months of life

	void GetNewStage(int ID, double p);
};

class Child
{
public:
	Child(int Sex);

	int SexInd; // 0 = male, 1 = female
	double PropnBirths;

	PaedHIV Perinatal;
	PaedHIV Breastmilk;
	double NonAIDSmort[15];
	double MortProb1st6m; // Prob of non-AIDS mortality in 1st 6 months of life
	
	void UpdateMort();
};

// --------------------------------------------------------------------------------------
// Classes for individuals
// --------------------------------------------------------------------------------------

class Indiv
{
public:
	Indiv();

	int AliveInd; // 0 = dead, 1 = alive
	int SexInd; // 0 = male, 1 = female (same as in TSHISA)
	int RiskGroup; // 1 = high risk, 2 = low risk (same as in TSHISA)
	double DOB; // Date of birth (calendar year, with decimal indicating timing of birth)
	int AgeGroup; // 5-year age group (note that in TSHISA, 0 corresponds to 10-14,
				  // but in our model, 0 corresponds to 0-4)
	double AgeExact;
	double Age50;
	int MarriedInd; // 0 = unmarried, 1 = married (same as in TSHISA)
	int FSWind; // 1 => Female sex worker (same as in TSHISA)
	int VirginInd; // 0 = sexually experienced, 1 = virgin (same as in TSHISA)
	int IDprimary; // ID of primary partner (0 => no primary partner)
	int ID2ndary; // ID of secondary partner (0 => no secondary partner)
	// Note that the individual's ID is always 1 + their row number in the array (so that the
	// first individual simulated has ID number 1, not 0).
	int CurrPartners; // Number of current sexual partners (max 2)
	int LifetimePartners; // Cumulative number of sexual partners (ever)
	double NonHIVmortProb; // Prob of dying from causes unrelated to HIV, per behaviour cycle
	double NonHIVfertRate; // Rate of giving birth in HIV-negative women, per behaviour cycle
	double CumSelectionProb; // Prob that this indiv (or someone with lower ID) is selected
	int HIVstage; // 0 for HIV-neg, 1 = acute, 2 = latent, 3 = WHO stage III, 4 = AIDS, 5 = ART, 6=interrupted
	int HIVstageE; // HIV stage at end of STD cycle
	int ARTweeks; //duration of ART
	int ARTstage; //if ART was initiated early (stage 1/2) then 0, if late, then 1
	int CTstage; // 0 for susceptible, 1 = symptomatic, 2 = asymptomatic infected, 3 = immune
	int CTstageE; // CT stage at end of STD cycle
	int HDstage; // 0 for susceptible, 1 = symptomatic, 2 = asymptomatic infected, 3 = immune
	int HDstageE; // HD stage at end of STD cycle
	int NGstage; // 0 for susceptible, 1 = symptomatic, 2 = asymptomatic infected, 3 = immune
	int NGstageE; // NG stage at end of STD cycle
	int TVstage; // 0 for susceptible, 1 = symptomatic, 2 = asymptomatic infected, 3 = immune
	int TVstageE; // TV stage at end of STD cycle
	int TPstage; // 0 for susceptible, 1 = incubation, 2 = primary, 3 = 2ndary, 4 = latent, 5/6 = immune
	int TPstageE; // TP stage at end of STD cycle
	int HSVstage; // 0 for susceptible, 1 = primary ulcer, 2 = transiently asymptomatic,
				  // 3 = recurrent ulcer, 4 = recurrent ulcer
	int HSVstageE; // HSV stage at end of STD cycle
	int BVstage; // 0 = normal, 1 = intermediate, 2 = symptomatic BV, 3 = asymptomatic BV
	int BVstageE; // BV stage at end of cycle
	int VCstage; // 0 = uninfected, 1 = asymptomatic VC, 2 = symptomatic VC
	int VCstageE; // VC stage at end of cycle
	
	//HPV types: 16	18	31	33	35	39	45	51	52	56	58	59	68
	//vector<int> HPVstage;
	int HPVstage[13]; // 0 = uninfected, 1 = infected, 2 = CIN1, 3 = CIN2, 4 = CIN3, 5 = CC-I, 6=latent, 7 = immune, 8 = CC-II, 9 = CC-III, 10 = CC-IV, 
					  // 11 = CC-I symptomatic,  12 = CC-II symptomatic, 13 = CC-III symptomatic, 14 = CC-IV symptomatic, 15 = recovered, 16 = Dead of cancer
	int HPVstageE[13]; // HPV stage at end of cycle
	static const vector<int> allhpv;
	static const vector<int> hpv1618;
	static const vector<int> hpv161845;
	
	static const vector<int> lsil;
	static const vector<int> hsil;
	static const vector<int> cc_un;
	static const vector<int> cc_diag;
	static const vector<int> recover;

	int VaccinationStatus[13];
	int WasLatent[13];
	int GotVacc;
	int GotVaccOffer;
	int TimeVacc;
	int ExpVacc;
	int TimeinCIN3[13]; //number of weeks spent in CIN3
	int WeibullCIN3[13]; //Time to cancer simulated when entering CIN3
	int StageIdeath;
	int StageIIdeath;
	int StageIIIdeath;
	int StageIVdeath;
	int StageIrecover;
	int StageIIrecover;
	int StageIIIrecover;
	int StageIVrecover;
	int TimeinStageI;
	int TimeinStageII;
	int TimeinStageIII;
	int TimeinStageIV;

	//Screening
	int InScreen; //0, change to 1 when person enters screening algorithm
	int ScreenCount;
	int InWHOScreen; //0, change to 1 when HIV+ person had first WHO screen
	int timePassed; //time since last screen 
	int timetoScreen; //time of next screen
	int timetoCol; //time of next screen
	int ScreenResult; //0=Normal; 1=LSIL; 2=HSIL; 3=CC
	int ColResult; //0=Normal; 1=CIN1; 2=CIN2/3; 3=CC
	int TrueStage; //0=Normal; 1=LSIL; 2=HSIL; 3=CC; 4=DiagnosedCC; 5=Recovered/Dead
	int HPVstatus; //0=HPV negative (stages 0, 6, 7); 1=HPV positive
	int repeat; //0=no; 1=yes
	int HPVrepeat;
	double ThermalORPap;
	int reason; //0=routine screening; 1=diagnostic
	int DiagnosedCC; //0=no; 1=yes 
	int ARTnextScreen; //if 1, screen next round
	int GotTxV;
	//int ExpTxV;
	//WHO Screening
	int Scr30;
	int Scr40;
	int Scr50;
	int Scr35;
	int Scr45;
	int Scr16;
	int Scr19;
	int Scr22;
	int Scr25;
	int Scr28;
	int Scr31;
	int Scr34;
	int Scr37;
	int Scr43;
	int Scr46;
	int Scr49;
	int Scr52;
	int Scr55;
	int Scr58;
	int Scr61;
	int Scr64;

	double SuscepHIVadj; // Adjustment to individual's HIV susceptibility
	double DateInfect; // Date at which HIV was acquired
	double PartnerRateAdj; // adjustment to individual's rate of partnership formation and
						   // dissolution (allowing for variation between individuals)
	double DesiredNewPartners; // rate at which new partners are acquired
	int NewStatus; // 0 if indiv has not yet had status updated in current behav cycle
	int UVIprimary; // # unprotected acts of vaginal intercourse with primary partner
	int PVIprimary; // # protected acts of vaginal intercourse with primary partner
	int UVI2ndary; // # unprotected acts of vaginal intercourse with 2ndary partner
	int PVI2ndary; // # protected acts of vaginal intercourse with 2ndary partner
	int UVICSW; // # unprotected acts of vaginal intercourse with CSWs
	int PVICSW; // # protected acts of vaginal intercourse with CSWs
	int IDofCSW; // ID of CSW (if man visits CSW)

	int totUVIprimary; //Total # unprotected acts of vaginal intercourse with primary partners
	int totPVIprimary; //Total # protected acts of vaginal intercourse with primary partners
	int totUVI2ndary; //Total # unprotected acts of vaginal intercourse with secondary partners
	int totPVI2ndary; //Total # protected acts of vaginal intercourse with secondary partners
	int totUVICSW; // Total # unprotected acts of vaginal intercourse with CSWs
	int totPVICSW; // Total # protected acts of vaginal intercourse with CSWs
	int risk_for_start;
	void GetRiskGroup();
	int SexDebutAge;

	int SelectEvent(double rnd); // Output is 0 if no new event, 1 if acquire S1, 2 if acquire S2,
								 // 3 if marry primary, 4 if marry 2ndary, 5 if divorce, 6 if end  
								 // primary S, 7 if end 2ndary S
	double ConvertToDependent1(double rate1);
	double ConvertToDependent2(double rate1, double rate2, int ord);
	double ConvertToDependent3(double rate1, double rate2, double rate3, int ord);
	double ConvertToDependent4(double rate1, double rate2, double rate3, double rate4, int ord);
	void SimulateSexActs(int ID);
	int RandomSexGenerator(double p, double lambda);
	void GetNewHIVstate(int ID, double p, double p2, double p3, double p4, double RandAcceptTxV, double Rloss, double EfficacyD1Rand, double EfficacyD2Rand);
	void GetNewHSVstate(int ID, double p);
	void GetNewTPstate(int ID, double p);
	void GetNewHDstate(int ID, double p);
	void GetNewNGstate(int ID, double p);
	void GetNewCTstate(int ID, double p);
	void GetNewTVstate(int ID, double p);
	void GetNewHPVstate(int ID, double p, int type);
	void GetScreened(int ID, double rea,  double scr, double ade, double tts, double res, double ttC, double CCd, double SI, double SII, double SIII, double SIV, 
							double SId, double SIId, double SIIId, double SIVd, double acceptRand, double efficacyD1Rand, double efficacyD2Rand, double D2LossRand);	
	void AdministerTherapeuticVaccine(int ID, double acceptRand, double efficacyD1Rand, double efficacyD2Rand, double D2LossRand);
	void ScreenAlgorithm(int ID, double rea, double ade, double tts, double res, double ttC, double CCd, double SI, double SII, double SIII, double SIV, 
							double SId, double SIId, double SIIId, double SIVd, double acceptRand, double efficacyD1Rand, double efficacyD2Rand, double D2LossRand);
	void MnDScreenAlgorithm(int ID, double rea, double ade, double tts, double res, double ttC, double CCd, double SI, double SII, double SIII, double SIV, 
							double SId, double SIId, double SIIId, double SIVd);
	void HPVScreenAlgorithm(int ID, double rea, double ade, double tts, double res, double ttC, double CCd, double SI, double SII, double SIII, double SIV, 
							double SId, double SIId, double SIIId, double SIVd);
	void HPV_ThermalScreenAlgorithm(int ID, double rea, double ade, double tts, double res, double ttC, double CCd, double SI, double SII, double SIII, double SIV, 
							double SId, double SIId, double SIIId, double SIVd);
	void WHOGetScreened(int ID, double rea, double scr, double ade, double tts, double res, double ttC, double CCd, double SI, double SII, double SIII, double SIV, double clr, 
							double SId, double SIId, double SIIId, double SIVd, double acceptRand, double efficacyD1Rand, double efficacyD2Rand, double D2LossRand);
	void WHOScreenAlgorithm(int ID,  double tts, double ttC,double clr, double SI, double SII, double SIII, double SIV, 
							double SId, double SIId, double SIIId, double SIVd);
	void GetTreated(int ID, double res, double trt, double clr, double regr, double tts, double SI, double SII, double SIII, double SIV, 
							double SId, double SIId, double SIIId, double SIVd);
	void PerfectGetScreened(int ID, double rea,  double scr, double ade, double tts, double res, double ttC, double CCd, double SI, double SII, double SIII, double SIV, 
							double SId, double SIId, double SIIId, double SIVd, double acceptRand, double efficacyD1Rand, double efficacyD2Rand, double D2LossRand);
	
	void AssignTimeinCIN3(int age_group, double p, int type );
	static bool AnyHPV(const int* XXX, const vector<int> & type_subset, const vector<int> & stage_subset);
	
};

class Pop
{
public:
	Pop(int xxx);
	int PopPyramid61[61][136];
	int PopVaxx61[61][136];
	int PopPyramid[54][136];
	int PopPyramid9[6][136];
	int PopPyramidMale[54][136];
	int HPVprevVT[54][136];
	int HPVprevAll[54][136];
	int CIN2prev[54][136];
	
	int PopPyramidAll[18][136];
	int PopPyramidAllART[136];
	int BehavPyramid[18][6]; // virgin M, experienced unmarried M, married M (then F)
	int PrefMatrix[16][16]; // 1st index is age of male, 2nd index is age of female
	int NumberPartners[16][10]; // unmarried 0, unmarried 1, unmarried 2, married 1,
								// married 2 (males then females)
	int TotSex[16][8]; // ST UVI, ST PVI, LT UVI, LT PVI (males then females)
	int AdultHIVstageTrend[136][12]; // 2nd index is HIV stage in M (0-5) and F (6-11)
	int AdultHPVstageTrend[136][208]; // 2nd index is HPV stage in M (0-7) and F (8-15)
	//int AdultHPVstageTrendAge[18][208];// 2nd index is HPV stage in M (0-7) and F (8-15) (16) for each of the 13 types: 16*13=208
	int AdultHPVstageTrendAge[320][208];// 1st index is 16 age groups*20 risk groups; 2nd index is HPV stage in M (0-7) and F (8-15) (16) for each of the 13 types: 16*13=208
	double MacDprev[500][14];
	int NewScreen[54][136];
	int NewHPVScreen[54][136];
	int GetReferred[54][136];
	int NewColposcopy[54][136];
	int NewLLETZ[54][136];
	int NewUnnecessary[54][136];
	int NewVAT[54][136];
	int NewThermal[54][136];
	int NewVACC[108][136];
	int NewTxV[54][136];
	int TotScreens[8];
	int NewCancer[18][136];
	int StageIdiag[54][136];
	int StageIIdiag[54][136];
	int StageIIIdiag[54][136];
	int StageIVdiag[54][136];
	int WeeksInStageI[54][101];
	int WeeksInStageII[54][101];
	int WeeksInStageIII[54][101];
	int WeeksInStageIV[54][101];
	int NewCancerDeath[18][136];
	int NewCancerDeathHIV[54][136];
	int NewCancerART[136];
	int NewDiagCancer[54][136];
	int NewDiagCancer1618[18][136];
	int NewDiagCancerDeath[18][136];
	int HIVDeath[18][136];
	int HIVDeathM[18][136];
	int NewDiagCancerART[136];
	int StageDiag[4][136];
	int StageDiagNeg[4][136];
	int StageDiagNoART[4][136];
	int StageDiagOnART[4][136];
	//double ScreenProb[8][136]; //Probability of entering screening  by age (4 categories) + HIV/ART status (neg/noART, ART), over time
	double ModelCoverage[54][136];
	double ModelTxVCoverage[54][136];
	double ModelHPVCoverage[54][136];
	double ModelVaccCoverage[108][136];
	double ModelColpCoverage[54][136];
	double ModelLLETZCoverage[54][136];
	double ModelUnnecessaryCoverage[54][136];
	double ModelVATCoverage[54][136];
	double ModelThermalCoverage[54][136];
	double ModelGetReferred[54][136];
	//int TimeInCIN3mat[10000][2*samplesize];
	int StartingConditionsUpdate[288][208]; //first index is 16 age categories*20 risk groups
	int InitHPVStage[18][8];
	int CSWregister[MaxCSWs]; // Records IDs of women who are currently sex workers
	double EligibleByAgeM[16][2], EligibleByAgeF[16][2]; // Sum of DesiredNewPartners
		// in individuals eligible to acquire new partners, by age and risk group

	// Functions to assign initial profile
	void AssignAgeSex();
	void GetAllPartnerRates();
	void AssignBehav();
	void ChooseLTpartner(int ID, double rnd1, double rnd2);
	void ChooseSTpartner(int ID, double rnd1, double rnd2);
	void AssignHIV();
	void AssignHIV1990();
	void AssignVacc2024();
	void AssignSTIs();

	// Output functions
	void GetPopPyramid();
	void GetBehavPyramid();
	void GetPrefMatrix(int type); // type = 1 for married, 2 for unmarried
	void GetNumberPartners();
	void GetTotSex();
	void GetNumbersByHIVstage();
	void GetStartingConditions();
	
	void GetInitHIVprevH();
	void GetHHprevByAge();
	void GetANCprevByAge();
	void GetInitSTIprev();
	void GetCurrSTIprev();
	double GetQstatistic(int Cmatrix[3][3], int MatDim);
	void GetSTIconcordance();
	void SavePopPyramid(const char* filout);
	void SavePopPyramid9(const char* filout);
	void SaveCancerCases(const char* filout);
	void SaveCancerDeaths(const char* filout);
	void SaveBehavPyramid(char* filout);
	void SavePrefMatrix(char* filout);
	void SaveNumberPartners(char* filout);
	void SaveLifetimePartners(char* filout);
	void SaveLifetimeCIN3(const char* filout);
	void SaveTotSex(char* filout);
	void SaveMarriageRates(char* filout);
	void SaveAdultHIVstage(char* filout);
	void SaveStartingConditions(char* filout);
	// Functions updated every year
	void OneYear();
	void ResetFlow();
	void UpdateAgeGroup();
	void CalcNonHIVmort();
	void CalcNonHIVfert();
	void CalcMTCTrates();
	void UpdateProbCure();
	void UpdateSTDtransitionProbs();
	void UpdateCondomUse();
	void GetANCprev(STDtransition* a, int STDind); // 1 = HSV, 2 = TP, 3 = HD, 4 = NG, 5 = CT, 6 = TV, 7 = BV, 8 = VC, 9=HPV
	void GetFPCprev(STDtransition* a, int STDind);
	void GetCSWprev(STDtransition* a, int STDind);
	
	//HPV
	void SaveDurationofCIN3(char* filout);
	void SaveInitHPVstage(char* filout, int type);
	void GetHHprev(STDtransition* a, int STDind);
	void GetHHNprev(STDtransition* a, int STDind);
	void GetNOARTprev(STDtransition* a, int STDind); 
	void GetONARTprev(STDtransition* a, int STDind);
	void GetCurrHPVprev(int WhichType);
	void GetCurrHPVstage();
	void GetNumbersByHPVstage();
	void GetNumbersByHPVstageAge();
	void SaveAdultHPVstage(char* filout);
	void GetMacDprev();
	void SaveMacDprev(const char* filout);
	void SaveNewScreen(const char* filout);
	void SaveWeeksInStage(const char* filout);
	void SaveStagediag(const char* filout);
	void SaveAdultHPVstageAge(const char* filout);
	void GetDurationofCIN3();
	//void GetDurationofCIN3age();
	void GetInitHPVstage(int type);
	void HitTargets(int type);
    void CalcModelCoverage();
	void CalcLSE();
	void CalcStageatDiag();

	// Functions updated every sexual behaviour cycle
	void OneBehavCycle();
	void UpdateBirths();
	void UpdateNonHIVmort();
	void UpdateFSWprofile();
	void UpdateMarriageIncidence();
	void BalanceSexualPartners();
	void CalcPartnerTransitions();
	void GetNewPartner(int ID, double rnd, double rnd2, int rsk);
	void SetNewStatusTo1(int ID);
	void SetToDead(int ID);
	void NewBirth(int ID);

	// Functions updated every STD cycle
	void OneSTDcycle();
	void GetSexActs();
	void GetHIVtransitions();
	void GetSTDtransitions();
};

class Partner
{
public:
	Partner();

	int ID;
	double DesiredNewPartners; // rate at which new partners are acquired
};

class PartnerCohort
{
	// This class defines a collection of individuals of the same age group, risk group 
	// and sex, who are eligible to form new partnerships.

public:
	PartnerCohort();
	
	std::vector<Partner> Pool;
	double TotalDesire;

	void Reset(); // Note that this function should be redundant in the new model.
	void AddMember(int PID);
	int SamplePool(double rand1); // Returns the ID of a randomly chosen partner
};

// --------------------------------------------------------------------------------------
// General functions
// --------------------------------------------------------------------------------------

// Functions for reading input parameters
void ReadCCStrategies();
void ReadSexAssumps(const char* input);
void ReadSTDepi(const char* input);
void ReadRatesByYear();
void ReadMortTables();
void ReadFertTables();
void ReadOneStartProfileM(std::ifstream* file, int group);
void ReadOneStartProfileF(std::ifstream* file, int group);
void ReadStartProfile(const char* input);
void ReadInitHIV();
void ReadStartPop();
void ReadSTDprev();
void ReadSTDparameters();
void ReadTimeinCIN3();
void ReadAllInputFiles();
void ReadLHSamples();
void ReadHPVparams();
void ReadScreenData();
void ReadLifeTimePartners();

// Mathematical functions

//double GetQstatistic(int Cmatrix[3][3], int MatDim);

// Calibration/likelihood functions

void SetCalibParameters();
void CalcTotalLogL();

// Functions for running simulations and uncertainty analysis
void InitialiseHPV();
void OneSimulation();
void RunSimulations();
void StoreOutputs();
void AggregateSims();
void SimulateParameters();
void SimulateHIVparams();
void SimulateTPparameters();
void SimulateHSVparameters();
void SimulateNGparameters();
void SimulateCTparameters();
void SimulateTVparameters();
void SimulateBVparameters();
void SimulateVCparameters();

void SimulateHPVparamsType();
void SimulateHPVparams16();
void SimulateHPVparams18();
void SimulateHPVparamsOTHER(int type);
//void SimulateScreeningProbs();
void ProspectiveCohort(int cycle);
// --------------------------------------------------------------------------------------
// Objects created from classes defined above
// --------------------------------------------------------------------------------------

PostOutputArray BVprev15to49F(30);
PostOutputArray CTprev15to49F(30);
PostOutputArray CTprev15to49M(30);
PostOutputArray CTprevCSW(30);
PostOutputArray HDprev15to49F(21);
PostOutputArray HDprev15to49M(21);
PostOutputArray HIVprev15to49F(136);
PostOutputArray HIVprev15to49all(136);
PostOutputArray HIVprev15to49M(136);
PostOutputArray ARTcov15to49M(136);
PostOutputArray ARTcov15to49F(136);
PostOutputArray HPVprev15to49F(136);
PostOutputArray HPVprev15to49M(136);

PostOutputArray HSVprev15to49F(40);
PostOutputArray HSVprev15to49M(40);
PostOutputArray HSVprevCSW(30);
PostOutputArray HSVprevANC(30);
PostOutputArray NGprev15to49F(30);
PostOutputArray NGprev15to49M(30);
PostOutputArray NGprevCSW(30);
PostOutputArray TPprev15to49F(40);
PostOutputArray TPprev15to49M(40);
PostOutputArray TPprevCSW(30);
PostOutputArray TPprevANC(30);
PostOutputArray TVprev15to49F(30);
PostOutputArray TVprev15to49M(30);
PostOutputArray TVprevCSW(30);
PostOutputArray VCprev15to49F(30);
PostOutputArray NewCT(30);
//PostOutputArray NewHIV(30);
PostOutputArray NewHSV(30);
PostOutputArray NewNG(30);
PostOutputArray NewTP(30);
PostOutputArray NewTV(30);

PostOutputArray PrevPreg15(23);
PostOutputArray PrevPreg20(23);
PostOutputArray PrevPreg25(23);
PostOutputArray PrevPreg30(23);
PostOutputArray PrevPreg35(23);
PostOutputArray PrevPreg40(23);
PostOutputArray PrevPregTotal(23);
PostOutputArray PrevHH2005(18);
PostOutputArray PrevHH2008(18);
PostOutputArray PrevHH2012(18);
PostOutputArray PrevHH2017(18);

/*PostOutputArray CTconcordance(21);
PostOutputArray HDconcordance(21);
PostOutputArray HIVconcordance(21);
PostOutputArray HSVconcordance(21);
PostOutputArray NGconcordance(21);
PostOutputArray TPconcordance(21);
PostOutputArray TVconcordance(21);*/

PostOutputArray2 HIVparamsLogL(9);
PostOutputArray2 TPparamsLogL(11);
PostOutputArray2 HSVparamsLogL(12);
PostOutputArray2 NGparamsLogL(11);
PostOutputArray2 CTparamsLogL(11);
PostOutputArray2 TVparamsLogL(13);
PostOutputArray2 BVparamsLogL(9);
PostOutputArray2 VCparamsLogL(7);

PostOutputArray2 RandomUniformHIV(8);
PostOutputArray2 RandomUniformTP(10);
PostOutputArray2 RandomUniformHSV(11);
PostOutputArray2 RandomUniformNG(10);
PostOutputArray2 RandomUniformCT(10);
PostOutputArray2 RandomUniformTV(12);
PostOutputArray2 RandomUniformBV(8);
PostOutputArray2 RandomUniformVC(6);
PostOutputArray2 OutANCbias(1);
PostOutputArray2 OutModelVarANC(1);
PostOutputArray2 OutModelVarHH(2);

HIVtransition HIVtransitionM(0, 0, 0, 0, 0, 0, 0, 27, 0, 0);
HIVtransition HIVtransitionF(1, 0, 0, 0, 0, 0, 80, 27, 0, 0);
SyphilisTransition TPtransitionM(0, 0, 0, 0, 0, 6, 0, 0, 0, 0);
SyphilisTransition TPtransitionF(1, 22, 6, 0, 7, 5, 0, 0, 0, 0);
HerpesTransition HSVtransitionM(0, 0, 0, 0, 0, 10, 0, 0, 0, 0);
HerpesTransition HSVtransitionF(1, 2, 1, 0, 1, 10, 0, 0, 0, 0);
OtherSTDtransition NGtransitionM(0, 0, 0, 0, 0, 9, 0, 0, 0, 0);
OtherSTDtransition NGtransitionF(1, 11, 7, 0, 7, 9, 0, 0, 0, 0);
OtherSTDtransition CTtransitionM(0, 0, 0, 0, 0, 10, 0, 0, 0, 0);
OtherSTDtransition CTtransitionF(1, 11, 7, 0, 7, 10, 0, 0, 0, 0);
OtherSTDtransition HDtransitionM(0, 0, 0, 11, 0, 0, 0, 0, 0, 0);
OtherSTDtransition HDtransitionF(1, 0, 0, 3, 0, 0, 0, 0, 0, 0);
OtherSTDtransition TVtransitionM(0, 0, 0, 0, 0, 2, 0, 0, 0, 0);
OtherSTDtransition TVtransitionF(1, 11, 7, 0, 2, 3, 0, 0, 0, 0);
BVtransition BVtransitionF(1, 5, 4, 0, 0, 0, 0, 0, 0, 0);
VCtransition VCtransitionF(1, 5, 3, 0, 2, 0, 0, 0, 0, 0);
//HPVtransitionCC HPVtransitionCC(1, 0, 13, 0, 0, 222, 0, 0, 55, 44);
HPVtransitionCC HPVtransitionCC(1, 0, 13, 0, 0, 138, 0, 0, 28, 15);


//(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH, int ObsANCN, int ObsHHN, int ObsNOART, , int ObsONART)
Child MaleChild(0);
Child FemChild(1);

std::vector<Indiv> Register(InitPop);
Pop RSApop(0);
PartnerCohort MalePartners[16][2]; // Cohorts of men eligible to enter relationships, by age and sex
PartnerCohort FemPartners[16][2]; // Cohorts of women eligible to enter relationships, by age and sex

std::vector<HPVTransition> HPVTransM(13);
std::vector<HPVTransition> HPVTransF(13);

std::vector<PostOutputArray3> NewHPV(13);
std::vector<PostOutputArray3> NewCC(13);

//PostOutputArray3 TimeinCIN3age(30);
//PostOutputArray3 EverinCIN3age(30);
//PostOutputArray3 FemPop(136);
//PostOutputArray3 NewHIV;
//PostOutputArray2 NewScreen(136);

//std::vector<PostOutputArray4> HPVparamsLogL(13);
PostOutputArray2 HPVparamsLogL(39);
//std::vector<PostOutputArray4> RandomUniformHPV(13);
PostOutputArray2 RandomUniformHPV(34);
PostOutputArray2 CC_ASR(136);
PostOutputArray2 CC_diag_ASR(136);
PostOutputArray2 CC_diag_IR(136);
PostOutputArray2 CC_diag_ASR1618(136);
PostOutputArray2 CC_diag_ASR_death(136);
PostOutputArray2 CC_ASR2(136);
PostOutputArray2 CC_diag_ASR2(136);
PostOutputArray2 CC_diag_ASR16182(136);
PostOutputArray2 CC_diag_ASR_death2(136);
PostOutputArray2 CC_ASR_death2(136);
PostOutputArray2 CC_ASR3(136);
PostOutputArray2 CC_diag_ASR3(136);
PostOutputArray2 CC_diag_ASR16183(136);
PostOutputArray2 CC_diag_ASR_death3(136);
PostOutputArray2 CC_ASR_death3(136);
PostOutputArray2 NewCIN2neg(136);
PostOutputArray2 NewCIN2pos(136);
PostOutputArray2 NewCIN2art(136);

PostOutputArray2 CC_diag_ASR_ART(136);
PostOutputArray2 CC_20(136);
PostOutputArray2 CC_25(136);
PostOutputArray2 CC_30(136);
PostOutputArray2 CC_35(136);
PostOutputArray2 CC_40(136);
PostOutputArray2 CC_45(136);
PostOutputArray2 CC_50(136);
PostOutputArray2 CC_55(136);
PostOutputArray2 CC_60(136);
PostOutputArray2 CC_65(136);
PostOutputArray2 CC_70(136);
PostOutputArray2 CC_20diag(136);
PostOutputArray2 CC_25diag(136);
PostOutputArray2 CC_30diag(136);
PostOutputArray2 CC_35diag(136);
PostOutputArray2 CC_40diag(136);
PostOutputArray2 CC_45diag(136);
PostOutputArray2 CC_50diag(136);
PostOutputArray2 CC_55diag(136);
PostOutputArray2 CC_60diag(136);
PostOutputArray2 CC_65diag(136);
PostOutputArray2 CC_70diag(136);
PostOutputArray2 CC_75diag(136);
PostOutputArray2 StageI(136);
PostOutputArray2 StageII(136);
PostOutputArray2 StageIII(136);
PostOutputArray2 StageIV(136);

/*PostOutputArray4 NewHPVM;
PostOutputArray4 NewHPVF;
PostOutputArray4 HPVprev15to64HIVposM;
PostOutputArray4 HPVprev15to64HIVnegM;
PostOutputArray4 HPVprev15to64HIVposF;
PostOutputArray4 HPVprev15to64HIVnegF;
PostOutputArray4 HPVprev15to24HIVposF;
PostOutputArray4 HPVprev15to24HIVnegF;
PostOutputArray4 HPVprev15to24HIVposM;
PostOutputArray4 HPVprev15to24HIVnegM;
PostOutputArray4 HPVprev15to64F;
PostOutputArray4 HPVprev15to64M;
PostOutputArray4 HPVprev15to64ARTF;
PostOutputArray4 HPVprev15to64ARTM;
PostOutputArray4 HPVprevFHIVneg;
PostOutputArray4 HPVprevMHIVneg;
PostOutputArray4 HPVprevFHIVpos;
PostOutputArray4 HPVprevMHIVpos;
PostOutputArray4 HPVprevF;
PostOutputArray4 HPVprevM;

std::vector<PostOutputArray4> Targets(20);*/

//Other:
PostOutputArray ABprev30to65neg(136);
PostOutputArray ABprev30to65pos(136);
PostOutputArray ABprev30to65art(136);
PostOutputArray HSIL_ABprev30to65neg(136);
PostOutputArray HSIL_ABprev30to65pos(136);
PostOutputArray HSIL_ABprev30to65art(136);
PostOutputArray ABprev18to60noart(136);
PostOutputArray HSIL_ABprev18to60noart(136);
PostOutputArray ABprev18to65art(136);
PostOutputArray HSIL_ABprev18to65art(136);
PostOutputArray HSILprev18to65art(136);

PostOutputArray HPVprev15to64HIVposF(136);
PostOutputArray HPVprev15to64HIVnegF(136);
PostOutputArray HPVprev15to24HIVposF(136);
PostOutputArray HPVprev15to24HIVnegF(136);
PostOutputArray HPVprev15to64F(136);
PostOutputArray HPVprev15to64ARTF(136);
PostOutputArray HPVprevFHIVneg(136);
PostOutputArray HPVprevFHIVpos(136);
PostOutputArray HPVprevF(136);

PostOutputArray DiagCCPost2000(2);
PostOutputArray TrueStageByYear(136);
//PostOutputArray GetReferred(136);
PostOutputArray GetTreatment(136);
//PostOutputArray Targets(160);

//WHO outputs
PostOutputArray HPVprev15to65allF(136);
PostOutputArray HPVprev15to65allM(136);
PostOutputArray HPVprev15to65negF(136);
PostOutputArray HPVprev15to65negM(136);
PostOutputArray HPVprev15to65posF(136);
PostOutputArray HPVprev15to65posM(136);
PostOutputArray HPVprev15to65noartF(136);
PostOutputArray HPVprev15to65noartM(136);
PostOutputArray HPVprev15to65artF(136);
PostOutputArray HPVprev15to65artM(136);
PostOutputArray HSILprevALL(136);
PostOutputArray HSILprevNEG(136);
PostOutputArray HSILprevPOS(136);
PostOutputArray HSILprevNOART(136);
PostOutputArray HSILprevART(136);
PostOutputArray CCprevALL(136);
PostOutputArray CCprevNEG(136);
PostOutputArray CCprevPOS(136);
PostOutputArray CCprevNOART(136);
PostOutputArray CCprevART(136);
