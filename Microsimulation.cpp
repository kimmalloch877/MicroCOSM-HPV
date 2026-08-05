// This is the main project file for VC++ application project
// generated using an Application Wizard.

//#include "stdafx.h"
 using namespace std;

// using <mscorlib.dll>

// using namespace System;
#include "Microsimulation.h"
#include "StatFunctions.h"
#include "randomc.h"
#include <time.h>
#include <cstring>
#include <sstream>
#include <string>
#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;	

CRandomMersenne rg(111);
int process_num;

int main(int argc, char *argv[])
{
	clock_t start, finish;
	double elapsed_time;

	start = clock();

	if (argc > 1) { // There is at least one argument
		process_num = atoi(argv[1]); // Convert the first argument to an integer
		rg = CRandomMersenne(process_num); // Reinitialize the RNG
	}
	else {
		process_num = 0; // No arguments, so set the process_num to zero.
	}
	
	RunSimulations();


	finish = clock();
	elapsed_time = (finish - start);
	std::cout << "Time taken: " << elapsed_time << endl;
	//system("PAUSE");
	return 0;
	
}
PrevalenceData::PrevalenceData(){
	LogL = 0.0;
	LogL0 = 0.0;
	LogL1 = 0.0;
	LogL2 = 0.0;
	LogL3 = 0.0;
}

NationalData::NationalData()
{
	int ia;

	for (ia = 0; ia<16; ia++){
		BiasMult[ia] = 1.0;
	}
	ModelVarEst = 0.0;
}

void NationalData::CalcModelVar()
{
	// This function should ONLY get called if GetSDfromData = 1

	int ia, ic;
	double AdjPrev, StochasticVar;
	double VarLogitPrev[100]; // Variance of the logit-transformed prevalence estimates 
	// (sampling variation only). Change dimension if n>100.
	double SpSum, BiasSum, BiasLevel;

	BiasSum = 0.0;
	for (ia = 0; ia<16; ia++){
		BiasSum += log(BiasMult[ia]);
	}
	SpSum = 0.0;
	for (ic = 0; ic<Observations; ic++){
		SpSum += 1.0 - ExpSp[ic]; //Specificity
	}

	BiasLevel = 0.0;
	if (BiasSum > 0.0 || SpSum > 0.0){ // ANC data
		for (ic = 0; ic<Observations; ic++){
			BiasLevel += log(StudyPrev[ic] / (1.0 - StudyPrev[ic])) - log(ModelPrev[ic] / (1.0 -
				ModelPrev[ic]));
		}
		BiasLevel = BiasLevel / Observations;
		if (BiasLevel < 0.0){ BiasLevel = 0.0; }
	}

	ModelVarEst = 0.0;
	for (ic = 0; ic<Observations; ic++){
		AdjPrev = 1.0 / (1.0 + (1.0 / ModelPrev[ic] - 1.0) * exp(-BiasLevel));
		StochasticVar = 0.0;
		VarLogitPrev[ic] = pow(PrevSE[ic] / (StudyPrev[ic] * (1.0 - StudyPrev[ic])), 2.0) + StochasticVar;
		if (BiasSum > 0.0 || SpSum > 0.0){ VarLogitPrev[ic] += ANCbiasVar; }
		else{ VarLogitPrev[ic] += HSRCbiasVar; }
		ModelVarEst += pow(log(StudyPrev[ic] / (1.0 - StudyPrev[ic])) - log(AdjPrev / (1.0 -
			AdjPrev)), 2.0) - VarLogitPrev[ic];
	}
	ModelVarEst = ModelVarEst / Observations;
}

void NationalData::CalcLogL()
{
	int ia, ic, SimCount2;
	double AdjPrev, StochasticVar;
	double VarLogitPrev[100]; // Variance of the logit-transformed prevalence estimates 
	// (sampling variation only). Change dimension if n>100.
	double SpSum, BiasSum, BiasLevel;

	SimCount2 = (CurrSim - 1) / IterationsPerPC;
	BiasSum = 0.0;
	for (ia = 0; ia<16; ia++){
		BiasSum += log(BiasMult[ia]);
	}
	SpSum = 0.0;
	for (ic = 0; ic<Observations; ic++){
		SpSum += 1.0 - ExpSp[ic];
	}

	BiasLevel = 0.0;
	if (BiasSum > 0.0 || SpSum > 0.0){ // ANC data
		for (ic = 0; ic<Observations; ic++){
			BiasLevel += log(StudyPrev[ic] / (1.0 - StudyPrev[ic])) - log(ModelPrev[ic] / (1.0 -
				ModelPrev[ic]));
		}
		BiasLevel = BiasLevel / Observations;
		if (BiasLevel < 0.0){ BiasLevel = 0.0; }
		if (FixedUncertainty == 1){ OutANCbias.out[SimCount2][0] = exp(BiasLevel); }
	}

	LogL = 0.0;
	for (ic = 0; ic<Observations; ic++){
		AdjPrev = 1.0 / (1.0 + (1.0 / ModelPrev[ic] - 1.0) * exp(-BiasLevel));
		if (GetSDfromData == 0){
			LogL += -0.5 * log(2.0 * 3.141592654) - log(PrevSE[ic])
				- 0.5 * pow((StudyPrev[ic] - AdjPrev) / PrevSE[ic], 2.0);
			//cout << ic << " " << 	ModelPrev[ic] << " " << StudyPrev[ic] << " " << BiasLevel << " " << AdjPrev << " " << -0.5 * log(2.0 * 3.141592654) - log(PrevSE[ic])
				//- 0.5 * pow((StudyPrev[ic] - AdjPrev) / PrevSE[ic], 2.0) << endl;
		}
		else{
			StochasticVar = 0.0;
			//StochasticVar = exp(2.0 * (111.271 - 0.0561 * StudyYear[ic] - 4.299 * ModelPrev[ic] +
			//	3.342 * pow(ModelPrev[ic], 2.0))) / IterationsPerPC;
			VarLogitPrev[ic] = pow(PrevSE[ic] / (StudyPrev[ic] * (1.0 - StudyPrev[ic])), 2.0) + StochasticVar;
			if (BiasSum > 0.0 || SpSum > 0.0){ VarLogitPrev[ic] += ANCbiasVar; }
			else{ VarLogitPrev[ic] += HSRCbiasVar; }
		}
		ModelPrev[ic] = AdjPrev;
	}

	if (GetSDfromData == 1){
		if (FixedUncertainty == 1){
			if (BiasSum > 0.0 || SpSum > 0.0){
				OutModelVarANC.out[SimCount2][0] = ModelVarEst;
			}
			else if (SexInd == 1){
				OutModelVarHH.out[SimCount2][0] = ModelVarEst;
			}
			else{
				OutModelVarHH.out[SimCount2][1] = ModelVarEst;
			}
		}
		for (ic = 0; ic<Observations; ic++){
			LogL += -0.5 * (log(2.0 * 3.141592654 * (VarLogitPrev[ic] + ModelVarEst)) +
				pow(log(StudyPrev[ic] / (1.0 - StudyPrev[ic])) - log(ModelPrev[ic] / (1.0 -
				ModelPrev[ic])), 2.0) / (VarLogitPrev[ic] + ModelVarEst));
			
		}
	}
}

AntenatalN::AntenatalN(){}

HouseholdN::HouseholdN(){}

SentinelData::SentinelData()
{
	BiasMult = 1.0;
}

void SentinelData::CalcLogL()
{
	int ic,xx;
	double Mean, Var; // The mean and variance of theta(i), the modelled prevalence
	// of the STD in study i
	double MeanFb, VarFb; // The mean and variance of f(bi), where bi is the 'random
	// effect' for study i
	double alpha, beta; // The alpha and beta parameters for the beta prior on theta(i)
	double a, b, c, d; // Arguments for the gamma functions
	double LogLikelihood;

	
	//ofstream file("test.txt", std::ios::app);
	
	LogL = 0.0;
	LogL0 = 0.0;
	LogL1 = 0.0;
	LogL2 = 0.0;
	LogL3 = 0.0;
	for (ic = 0; ic<Observations; ic++){
		
		ModelPrev[ic] = 1.0*ModelNum[ic]/ModelDenom[ic];

		if (ModelPrev[ic]<0.0001){ ModelPrev[ic] = 0.0001; } // To be consistent with GetPrev function for HIV
		if (ModelPrev[ic]>0.9999){ ModelPrev[ic] = 0.9999; }

		if(Perfect[ic]==1){
			if(Cyt_cat[ic]==1||Cyt_cat[ic]==2) {VarStudyEffect = VarStudyEffectH;}
			else if(Cyt_cat[ic]==0){VarStudyEffect = VarStudyEffectAB;}
			else{
				for (int xx = 0; xx < 13; xx++){ 
					if(Cyt_cat[ic]==xx+3) VarStudyEffect = HPVTransM[xx].HouseholdLogL.VarStudyEffect; 
				}
			}
			MeanFb = ModelPrev[ic] + VarStudyEffect * ModelPrev[ic] * (1.0 - ModelPrev[ic]) *
			(0.5 - ModelPrev[ic]); //Third order Taylor approximation Eq 4.9 in PhD
			VarFb = pow(ModelPrev[ic] * (1.0 - ModelPrev[ic]), 2.0) * (VarStudyEffect +
			(1.5 - 8.0 * ModelPrev[ic] + 8.0 * pow(ModelPrev[ic], 2.0)) * pow(VarStudyEffect,
			2.0) + pow(pow(ModelPrev[ic], 2.0) - ModelPrev[ic] + 1.0 / 6.0, 2.0) * 15.0 *
			pow(VarStudyEffect, 3.0)); //Third order Taylor approximation Eq 4.11 in PhD
			
			Mean = MeanFb; 
			Var = VarFb; //Eq 4.13+4.14 in PhD
		}
		
		StudyPos[ic] = StudyN[ic] * StudyPrev[ic];
		
		if (Var>0){
			alpha = Mean * (Mean * (1.0 - Mean) / Var - 1.0); //Eq 4.15 in PhD
			beta = (1.0 - Mean) * (Mean * (1.0 - Mean) / Var - 1.0); //Eq 4.15 in PhD

			a = alpha + beta;
			b = alpha + StudyPos[ic];
			c = beta + StudyN[ic] - StudyPos[ic];
			d = alpha + beta + StudyN[ic];

			LogLikelihood = gamma_log(&a) + gamma_log(&b) + gamma_log(&c) - gamma_log(&d) -
				gamma_log(&alpha) - gamma_log(&beta); //Eq 4.19 in PhD
			//cout << Mean <<" "<<Var<<" "<<gamma_log(&a) << " " << gamma_log(&b) <<" "<<gamma_log(&c)<<" "<<gamma_log(&d)<<" "<<gamma_log(&alpha)<<" "<<gamma_log(&beta)<<endl;
		}
		else{
			// In this case, the mean is fixed, so the likelihood is just the binomial.
			a = StudyN[ic] + 1.0;
			b = StudyPos[ic] + 1.0;
			c = StudyN[ic] - StudyPos[ic] + 1.0;
			LogLikelihood = gamma_log(&a) - gamma_log(&b) - gamma_log(&c) + StudyPos[ic] *
				log(Mean) + (StudyN[ic] - StudyPos[ic]) * log(1.0 - Mean);
			//cout << gamma_log(&a) << " " << gamma_log(&b) <<" "<<gamma_log(&c)<<" "<<StudyPos[ic] *	log(Mean)<<" "<<(StudyN[ic] - StudyPos[ic]) * log(1.0 - Mean)<<endl;	
		}
		LogL += LogLikelihood;
		if(Cyt_cat[ic]==0) { LogL0 += LogLikelihood;}
		if(Cyt_cat[ic]==1) { LogL1 += LogLikelihood;}
		if(Cyt_cat[ic]==2) { LogL2 += LogLikelihood;}
		if(Cyt_cat[ic]>2) { LogL3 += LogLikelihood;}
		//file << ic << " " << Cyt_cat[ic] << " " << ModelPrev[ic] << " " << ModelDenom[ic] << " " <<  StudyN[ic] << " " <<  StudyPrev[ic] << " " <<  Mean << " " << LogLikelihood << endl;
		//cout << ic << " " << ModelPrev[ic] << " " <<ModelDenom[ic] << " " <<  StudyN[ic] << " " <<  StudyPrev[ic] << " " <<  Mean << " " << LogLikelihood << endl;
		ModelNum[ic]=0;
		ModelDenom[ic]=0;
	}
	//file.close();
}

Household::Household(){}

NonHousehold::NonHousehold(){}

ANC::ANC(){}

FPC::FPC(){}

GUD::GUD(){}

CSW::CSW(){}

ART::ART(){}

PostOutputArray::PostOutputArray(int n)
{
	columns = n;
}

void PostOutputArray::RecordSample(const char *filout)
{
	int i, c;
	ostringstream s;

	if (process_num > 0){
		s << process_num << filout;
	}
	else{
		s << filout;
	}

	//ofstream file(filout);
	string path = "./output/" + s.str();
	ofstream file(path.c_str()); // Converts s to a C string

	for (i = 0; i<samplesize; i++){
		file << setw(6) << right << i << "	";
		for (c = 0; c<columns; c++){
			file << "	" << setw(10) << right << out[i][c];
		}
		file << endl;
	}
	file.close();
}

PostOutputArray2::PostOutputArray2(int n)
{
	columns = n;
}

void PostOutputArray2::RecordSample(const char *filout)
{
	int i, c;
	ostringstream s;

	if (process_num > 0){
		s <<  process_num << filout;
	}
	else{
		s << filout;
	}
	
	//ofstream file(filout);
	string path = "./output/" + s.str();
	ofstream file(path.c_str()); // Converts s to a C string

	for (i = 0; i<ParamCombs; i++){
		file << setw(6) << right << i << "	";
		for (c = 0; c<columns; c++){
			file << "	" << setw(10) << right << setprecision(9) << out[i][c];
		}
		file << endl;
	}
	file.close();
}

STDtransition::STDtransition(){}

void STDtransition::ReadPrevData(const char *input)
{
	int ic;
	ifstream file;

	file.open(input);
	if (file.fail()) {
		cerr << "Could not open input file.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	if (ANClogL.Observations>0){
		for (ic = 0; ic<ANClogL.Observations; ic++){
			file >> ANClogL.StudyYear[ic] >> ANClogL.StudyN[ic] >> ANClogL.StudyPrev[ic] >>
				ANClogL.ExpSe[ic] >> ANClogL.ExpSp[ic] >> ANClogL.VarSe[ic] >> ANClogL.VarSp[ic] >>
				ANClogL.HIVprevInd[ic] >> ANClogL.HIVprev[ic];
		}
		file.ignore(255, '\n');
	}
	file.ignore(255, '\n');
	/*if (FPClogL.Observations>0){
		for (ic = 0; ic<FPClogL.Observations; ic++){
			file >> FPClogL.StudyYear[ic] >> FPClogL.StudyN[ic] >> FPClogL.StudyPrev[ic] >>
				FPClogL.ExpSe[ic] >> FPClogL.ExpSp[ic] >> FPClogL.VarSe[ic] >> FPClogL.VarSp[ic] >>
				FPClogL.HIVprevInd[ic] >> FPClogL.HIVprev[ic] >> FPClogL.Cyt_cat[ic];
		}
		file.ignore(255, '\n');
	}*/
	if (FPClogL.Observations>0){
		for (ic = 0; ic<FPClogL.Observations; ic++){
			file >> FPClogL.StudyYear[ic] >> FPClogL.StudyN[ic] >> FPClogL.StudyPrev[ic] >>
				FPClogL.HIVprevInd[ic] >> FPClogL.HIVprev[ic] >> FPClogL.Cyt_cat[ic] >>
				FPClogL.Perfect[ic] >> FPClogL.ARTind[ic];
				//cout << ic << " " << FPClogL.StudyPrev[ic] << endl;
			}
		file.ignore(255, '\n');
	}
	file.ignore(255, '\n');
	if (GUDlogL.Observations>0){
		for (ic = 0; ic<GUDlogL.Observations; ic++){
			file >> GUDlogL.StudyYear[ic] >> GUDlogL.StudyN[ic] >> GUDlogL.StudyPrev[ic] >>
				GUDlogL.ExpSe[ic] >> GUDlogL.ExpSp[ic] >> GUDlogL.VarSe[ic] >> GUDlogL.VarSp[ic] >>
				GUDlogL.HIVprevInd[ic] >> GUDlogL.HIVprev[ic];
		}
		file.ignore(255, '\n');
	}
	file.ignore(255, '\n');
	if (CSWlogL.Observations>0){
		for (ic = 0; ic<CSWlogL.Observations; ic++){
			file >> CSWlogL.StudyYear[ic] >> CSWlogL.StudyN[ic] >> CSWlogL.StudyPrev[ic] >>
				CSWlogL.ExpSe[ic] >> CSWlogL.ExpSp[ic] >> CSWlogL.VarSe[ic] >> CSWlogL.VarSp[ic] >>
				CSWlogL.HIVprevInd[ic] >> CSWlogL.HIVprev[ic];
		}
		file.ignore(255, '\n');
	}
	file.ignore(255, '\n');
	/*if (HouseholdLogL.Observations>0){
		for (ic = 0; ic<HouseholdLogL.Observations; ic++){
			file >> HouseholdLogL.StudyYear[ic] >> HouseholdLogL.StudyN[ic] >>
				HouseholdLogL.StudyPrev[ic] >> HouseholdLogL.ExpSe[ic] >>
				HouseholdLogL.ExpSp[ic] >> HouseholdLogL.VarSe[ic] >> HouseholdLogL.VarSp[ic] >>
				HouseholdLogL.HIVprevInd[ic] >> HouseholdLogL.HIVprev[ic] >>
				HouseholdLogL.AgeStart[ic] >> HouseholdLogL.AgeEnd[ic] >>
				HouseholdLogL.ExclVirgins[ic] >> HouseholdLogL.Cyt_cat[ic];
		}
		file.ignore(255, '\n');
	}*/
	if (HouseholdLogL.Observations>0){
		for (ic = 0; ic<HouseholdLogL.Observations; ic++){
			file >> HouseholdLogL.StudyYear[ic] >> HouseholdLogL.StudyN[ic] >>
				HouseholdLogL.StudyPrev[ic] >> 
				HouseholdLogL.HIVprevInd[ic] >> HouseholdLogL.HIVprev[ic] >>
				HouseholdLogL.AgeStart[ic] >> HouseholdLogL.AgeEnd[ic] >>
				HouseholdLogL.ExclVirgins[ic] >> HouseholdLogL.Cyt_cat[ic] >> 
				HouseholdLogL.Perfect[ic] >> HouseholdLogL.ARTind[ic];
			}
		file.ignore(255, '\n');
	}
	file.ignore(255, '\n');
	if (AntenatalNlogL.Observations>0){
		for (ic = 0; ic<AntenatalNlogL.Observations; ic++){
			file >> AntenatalNlogL.StudyYear[ic] >> AntenatalNlogL.AgeStart[ic] >>
				AntenatalNlogL.StudyPrev[ic] >> AntenatalNlogL.PrevSE[ic] >>
				AntenatalNlogL.ExpSe[ic] >> AntenatalNlogL.ExpSp[ic];
		}
		file.ignore(255, '\n');
	}
	file.ignore(255, '\n');
	if (HouseholdNlogL.Observations>0){
		for (ic = 0; ic<HouseholdNlogL.Observations; ic++){
			file >> HouseholdNlogL.StudyYear[ic] >> //HouseholdNlogL.StudyN[ic] >>
				HouseholdNlogL.StudyPrev[ic] >> HouseholdNlogL.PrevSE[ic] >> HouseholdNlogL.ExpSe[ic] >>
				HouseholdNlogL.ExpSp[ic] >> //HouseholdNlogL.VarSe[ic] >> HouseholdNlogL.VarSp[ic] >>
				HouseholdNlogL.HIVprevInd[ic] >> HouseholdNlogL.HIVprev[ic] >>
				HouseholdNlogL.AgeStart[ic] >> HouseholdNlogL.AgeEnd[ic] >>
				HouseholdNlogL.ExclVirgins[ic] >> HouseholdNlogL.Cyt_cat[ic];

				/*	file >> HouseholdNlogL.StudyYear[ic] >> HouseholdNlogL.AgeStart[ic] >>
				HouseholdNlogL.StudyPrev[ic] >> HouseholdNlogL.PrevSE[ic] >>
				HouseholdNlogL.ExpSe[ic] >> HouseholdNlogL.ExpSp[ic];*/
		}
		file.ignore(255, '\n');
	}
	file.ignore(255, '\n');
	/*if (NOARTlogL.Observations>0){
		for (ic = 0; ic<NOARTlogL.Observations; ic++){
			file >> NOARTlogL.StudyYear[ic] >> NOARTlogL.StudyN[ic] >>
				NOARTlogL.StudyPrev[ic] >> NOARTlogL.ExpSe[ic] >>
				NOARTlogL.ExpSp[ic] >> NOARTlogL.VarSe[ic] >> NOARTlogL.VarSp[ic] >>
				NOARTlogL.AgeStart[ic] >> NOARTlogL.AgeEnd[ic]>> NOARTlogL.Cyt_cat[ic];
			
			}
		file.ignore(255, '\n');
	}*/
	if (NOARTlogL.Observations>0){
		for (ic = 0; ic<NOARTlogL.Observations; ic++){
			file >> NOARTlogL.StudyYear[ic] >> NOARTlogL.StudyN[ic] >>
				NOARTlogL.StudyPrev[ic] >> NOARTlogL.HIVprevInd[ic] >> NOARTlogL.HIVprev[ic] >>
				NOARTlogL.AgeStart[ic] >> NOARTlogL.AgeEnd[ic]>> NOARTlogL.Cyt_cat[ic] >> 
				NOARTlogL.Perfect[ic] >> NOARTlogL.ARTind[ic];
			}
		file.ignore(255, '\n');
	}
	file.ignore(255, '\n');
	/*if (ONARTlogL.Observations>0){
		for (ic = 0; ic<ONARTlogL.Observations; ic++){
			file >> ONARTlogL.StudyYear[ic] >> ONARTlogL.StudyN[ic] >>
				ONARTlogL.StudyPrev[ic] >> ONARTlogL.ExpSe[ic] >>
				ONARTlogL.ExpSp[ic] >> ONARTlogL.VarSe[ic] >> ONARTlogL.VarSp[ic] >>
				ONARTlogL.AgeStart[ic] >> ONARTlogL.AgeEnd[ic] >> ONARTlogL.Cyt_cat[ic];
			
			}
		file.ignore(255, '\n');
	}*/
	if (ONARTlogL.Observations>0){
		for (ic = 0; ic<ONARTlogL.Observations; ic++){
			file >> ONARTlogL.StudyYear[ic] >> ONARTlogL.StudyN[ic] >>
				ONARTlogL.StudyPrev[ic] >> ONARTlogL.HIVprevInd[ic] >> ONARTlogL.HIVprev[ic] >>
				ONARTlogL.AgeStart[ic] >> ONARTlogL.AgeEnd[ic] >> ONARTlogL.Cyt_cat[ic]>> 
				ONARTlogL.Perfect[ic]>> ONARTlogL.ARTind[ic];
			}
		file.ignore(255, '\n');
	}
	file.close();
}

void STDtransition::GetCSWprev()
{
	int ic, iy;

	if (CSWlogL.Observations>0){
		for (ic = 0; ic<CSWlogL.Observations; ic++){
			iy = CSWlogL.StudyYear[ic] - StartYear;
			CSWlogL.ModelPrev[ic] = 0.05 * CSWprevUnsmoothed[iy - 3];
			CSWlogL.ModelPrev[ic] += 0.12 * CSWprevUnsmoothed[iy - 2];
			CSWlogL.ModelPrev[ic] += 0.20 * CSWprevUnsmoothed[iy - 1];
			CSWlogL.ModelPrev[ic] += 0.26 * CSWprevUnsmoothed[iy];
			CSWlogL.ModelPrev[ic] += 0.20 * CSWprevUnsmoothed[iy + 1];
			CSWlogL.ModelPrev[ic] += 0.12 * CSWprevUnsmoothed[iy + 2];
			CSWlogL.ModelPrev[ic] += 0.05 * CSWprevUnsmoothed[iy + 3];
		}
	}
}

void STDtransition::SetVarStudyEffect(double Variance)
{
	ANClogL.VarStudyEffect = Variance;
	FPClogL.VarStudyEffect = Variance;
	CSWlogL.VarStudyEffect = Variance;
	GUDlogL.VarStudyEffect = Variance;
	HouseholdLogL.VarStudyEffect = Variance;
	ONARTlogL.VarStudyEffect = Variance;
	NOARTlogL.VarStudyEffect = Variance;
}

NonHIVtransition::NonHIVtransition(){}

HIVtransition::HIVtransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD,
	int ObsCSW, int ObsHH, int ObsANCN, int ObsHHN, int ObsNOART, int ObsONART)
{
	SexInd = Sex;
	nStates = 6;

	ANClogL.Observations = ObsANC;
	FPClogL.Observations = ObsFPC;
	GUDlogL.Observations = ObsGUD;
	CSWlogL.Observations = ObsCSW;
	HouseholdLogL.Observations = ObsHH;
	AntenatalNlogL.Observations = ObsANCN;
	HouseholdNlogL.Observations = ObsHHN;
	NOARTlogL.Observations = ObsNOART;
	ONARTlogL.Observations = ObsONART;
	HouseholdNlogL.SexInd = Sex;
}

void HIVtransition::CalcTransitionProbs()
{
	int iy;

	iy = CurrYear - StartYear;
	From1to2 = 1.0 - exp(-(1.0/AveDuration[0]) * 52.0/CycleD);
	From2to3 = (1.0 - exp(-(1.0/AveDuration[1]) * 52.0/CycleD)) *
		exp(-RateARTstart[iy][0][SexInd] / CycleD);
	From2to5 = 1.0 - exp(-RateARTstart[iy][0][SexInd] / CycleD);
	From3to4 = (1.0 - exp(-(1.0 / AveDuration[2]) * 52.0 / CycleD)) * 
		exp(-RateARTstart[iy][1][SexInd] / CycleD);
	From3to5 = 1.0 - exp(-RateARTstart[iy][1][SexInd] / CycleD);
	From4to5 = (1.0 - exp(-RateARTstart[iy][2][SexInd] / CycleD)) * 
		exp(-(1.0 / AveDuration[3]) * 52.0 / CycleD);
	From5to6 = (1.0 - exp(-ARTinterruption[iy] / CycleD)) *
		exp(-(1.0 / AveDuration[4]) * 52.0 / CycleD);
	From6to5 = (1.0 - exp(-0.9 / CycleD)) *
		exp(-(1.0 / AveDuration[4]) * 52.0 / CycleD);
	From4toDead = (1.0 - exp(-(1.0/AveDuration[3]) * 52.0/CycleD));
	//From5toDead = (1.0 - exp(-(1.0/AveDuration[4]) * 52.0/CycleD));
	From5toDeadEarly =  (1.0 - exp(-0.011/ CycleD)); 
	From5toDeadLateShort =  (1.0 - exp(-0.05/ CycleD)); 
	From5toDeadLateLong =  (1.0 - exp(-0.033/ CycleD)); 
	From6toDead = (1.0 - exp(-(1.0/AveDuration[4]) * 52.0/CycleD));
	
}

void HIVtransition::GetPrev()
{
	// Aggregates the results for the relevant parameter combination

	int ii, ia, iy;

	if (SexInd == 1){
		for (ii = 0; ii < AntenatalNlogL.Observations; ii++){
			AntenatalNlogL.ModelPrev[ii] = 0.0;
		}
	}
	for (ii = 0; ii < HouseholdNlogL.Observations; ii++){
		HouseholdNlogL.ModelPrev[ii] = 0.0;
	}

	for (ii = 0; ii < IterationsPerPC; ii++){
		// Antenatal results
		if (SexInd == 1){
			for (iy = 0; iy < 16; iy++){
				AntenatalNlogL.ModelPrev[iy * 5] += PrevPreg15.out[CurrSim - 1 - ii][iy + 7];
				AntenatalNlogL.ModelPrev[iy * 5 + 1] += PrevPreg20.out[CurrSim - 1 - ii][iy + 7];
				AntenatalNlogL.ModelPrev[iy * 5 + 2] += PrevPreg25.out[CurrSim - 1 - ii][iy + 7];
				AntenatalNlogL.ModelPrev[iy * 5 + 3] += PrevPreg30.out[CurrSim - 1 - ii][iy + 7];
				AntenatalNlogL.ModelPrev[iy * 5 + 4] += PrevPreg35.out[CurrSim - 1 - ii][iy + 7];
			}
		}
		// Household survey results
		for (ia = 0; ia < 9; ia++){
			HouseholdNlogL.ModelPrev[ia] += PrevHH2005.out[CurrSim - 1 - ii][ia + 9 * SexInd];
			HouseholdNlogL.ModelPrev[ia + 9] += PrevHH2008.out[CurrSim - 1 - ii][ia + 9 * SexInd];
			HouseholdNlogL.ModelPrev[ia + 18] += PrevHH2012.out[CurrSim - 1 - ii][ia + 9 * SexInd];
		}
	}

	// Take averages
	if (SexInd == 1){
		for (ii = 0; ii < AntenatalNlogL.Observations; ii++){
			AntenatalNlogL.ModelPrev[ii] = AntenatalNlogL.ModelPrev[ii] / IterationsPerPC;
			if (AntenatalNlogL.ModelPrev[ii] == 0.0){ AntenatalNlogL.ModelPrev[ii] = 0.0001; }
		}
	}
	for (ii = 0; ii < HouseholdNlogL.Observations; ii++){
		HouseholdNlogL.ModelPrev[ii] = HouseholdNlogL.ModelPrev[ii] / IterationsPerPC;
		if (HouseholdNlogL.ModelPrev[ii] == 0.0){ HouseholdNlogL.ModelPrev[ii] = 0.0001; }
	}
}

double HIVtransition::GetTransmProb(int ID)
{
	int PID1, PID2,  is, ig, ia, ic, iy;
	int IRisk, PRisk;
	double NoTransmProb, SingleActProb;

	IRisk = Register[ID - 1].RiskGroup;
	NoTransmProb = 1.0;
	iy = CurrYear-StartYear;

	// Infection from primary partner
	if (Register[ID - 1].CurrPartners>0){
		PID1 = Register[ID - 1].IDprimary;
		
		if (Register[PID1 - 1].HIVstage>0){
			// Get base probability depending on partnership type, risk group
			PRisk = Register[PID1 - 1].RiskGroup;
			if (Register[ID - 1].MarriedInd == 1){
				if (PRisk == 1 && IRisk == 1){
					SingleActProb = TransmProb[4];
				}
				else if (PRisk == 2 && IRisk == 2){
					SingleActProb = TransmProb[6];
				}
				else{
					SingleActProb = TransmProb[5];
				}
				if (PRisk == 2 && NoViralTransm12 == 1){
					SingleActProb = 0.0;
				}
				if (PRisk == 2 && IRisk == 2 && NoViralTransm22 == 1){
					SingleActProb = 0.0;
				}
			}
			else{
				if (PRisk == 1 && IRisk == 1){
					SingleActProb = TransmProb[1];
				}
				else if (PRisk == 2 && IRisk == 2){
					SingleActProb = TransmProb[3];
				}
				else{
					SingleActProb = TransmProb[2];
				}
			}
			// Make adjustment for HIV stage of partner
			is = Register[PID1 - 1].HIVstage;
			if(is == 5 ) {SingleActProb *= (1.0 + ARTinfectiousness[iy]);}
			else {SingleActProb *= (1.0 + HIVinfecIncrease[is - 1]);}

			ig = Register[PID1 - 1].SexInd;
			if (ig == 0){
				SingleActProb *= RatioAsympToAveM;
			}
			else{
				SingleActProb *= RatioAsympToAveF;
			}
			// Make adjustment for age of individual
			ia = Register[ID - 1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			// Make adjustment for STI cofactors
			if (CofactorType>0){
				// Still need to add code here
				// Note that we will need to also change the previous code to refer to
				// InitHIVtransm when allowing for STI cofactors.
			}
			// Make adjustment for other sources of heterogeneity
			if (AllowHIVsuscepAdj == 1){
				SingleActProb *= Register[ID - 1].SuscepHIVadj;
			}
			if (SingleActProb>1.0){
				SingleActProb = 1.0;
			}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVIprimary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PVIprimary);
		}
	}

	// Infection from 2ndary partner
	if (Register[ID - 1].CurrPartners == 2){
		PID2 = Register[ID - 1].ID2ndary;
		if (Register[PID2 - 1].HIVstage>0){
			// Get base probability depending on risk group
			PRisk = Register[PID2 - 1].RiskGroup;
			if (PRisk == 1 && IRisk == 1){
				SingleActProb = TransmProb[1];
			}
			else if (PRisk == 2 && IRisk == 2){
				SingleActProb = TransmProb[3];
			}
			else{
				SingleActProb = TransmProb[2];
			}
			// Make adjustment for HIV stage of partner
			is = Register[PID2 - 1].HIVstage;
			if(is == 5 ) {SingleActProb *= (1.0 + ARTinfectiousness[iy]);}
			else {SingleActProb *= (1.0 + HIVinfecIncrease[is - 1]);}
			ig = Register[PID2 - 1].SexInd;
			if (ig == 0){
				SingleActProb *= RatioAsympToAveM;
			}
			else{
				SingleActProb *= RatioAsympToAveF;
			}
			// Make adjustment for age of individual
			ia = Register[ID - 1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			// Make adjustment for STI cofactors
			if (CofactorType>0){
				// Still need to add code here
			}
			// Make adjustment for other sources of heterogeneity
			if (AllowHIVsuscepAdj == 1){
				SingleActProb *= Register[ID - 1].SuscepHIVadj;
			}
			if (SingleActProb>1.0){
				SingleActProb = 1.0;
			}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVI2ndary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PVI2ndary);
		}
	}

	// Infection from CSW (relevant only to high-risk men)
	if (Register[ID - 1].UVICSW + Register[ID - 1].PVICSW > 0){
		PID2 = Register[ID - 1].IDofCSW;
		if (Register[PID2 - 1].HIVstage>0){
			SingleActProb = TransmProb[0];
			is = Register[PID2 - 1].HIVstage;
			if(is == 5 ) {SingleActProb *= (1.0 + ARTinfectiousness[iy]);}
			else {SingleActProb *= (1.0 + HIVinfecIncrease[is - 1]);}
			SingleActProb *= RatioAsympToAveF;
			// Make adjustment for age of individual
			ia = Register[ID - 1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			// Make adjustment for STI cofactors
			if (CofactorType>0){
				// Still need to add code here
			}
			// Make adjustment for other sources of heterogeneity
			if (AllowHIVsuscepAdj == 1){
				SingleActProb *= Register[ID - 1].SuscepHIVadj;
			}
			if (SingleActProb>1.0){
				SingleActProb = 1.0;
			}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVICSW);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PVICSW);
		}
	}
	int tpp = Register.size();
	// Infection from client (relevant only to CSWs)
	if (Register[ID - 1].FSWind == 1){
		ia = Register[ID - 1].AgeGroup - 2;
		
			for (ic = 0; ic<tpp; ic++){
			if (Register[ic].IDofCSW == ID && (Register[ic].UVICSW + Register[ic].PVICSW) > 0){
				if (Register[ic].HIVstage>0){
					SingleActProb = TransmProb[0];
					is = Register[ic].HIVstage;
					if(is == 5 ) {SingleActProb *= (1.0 + ARTinfectiousness[iy]);}
					else {SingleActProb *= (1.0 + HIVinfecIncrease[is - 1]);}
					SingleActProb *= RatioAsympToAveM * SuscepIncrease[ia];
					if (AllowHIVsuscepAdj == 1){
						SingleActProb *= Register[ID - 1].SuscepHIVadj;
					}
					// Make adjustment for STI cofactors
					if (CofactorType>0){
						// Still need to add code here
					}
					if (SingleActProb>1.0){
						SingleActProb = 1.0;
					}
					NoTransmProb *= pow(1.0 - SingleActProb, Register[ic].UVICSW);
					NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
						Register[ic].PVICSW);
				}
			}
		}
	}

	return 1.0 - NoTransmProb;
}

void HIVtransition::GetNewStage(int ID, double p)
{
	int is;

	is = Register[ID-1].HIVstage;
	if(is==1){
		if(p<From1to2){Register[ID-1].HIVstageE = 2;}
		else{Register[ID-1].HIVstageE = 1;}
	}
	if(is==2){
		if(p<From2to3){Register[ID-1].HIVstageE = 3;}
		else if (p<(From2to3+From2to5)){ Register[ID-1].HIVstageE = 5; }
		else{Register[ID-1].HIVstageE = 2;}
	}
	if(is==3){
		if(p<From3to4){
			Register[ID-1].HIVstageE = 4;
		}
		else if(p<From3to4+From3to5){Register[ID-1].HIVstageE = 5;}
		else{Register[ID-1].HIVstageE = 3;}
	}
	if(is==4){
		if(p<From4toDead){
			RSApop.SetToDead(ID);
			if(Register[ID-1].SexInd==1){RSApop.HIVDeath[Register[ID-1].AgeGroup][CurrYear-StartYear] += 1;}
			if(Register[ID-1].SexInd==0){RSApop.HIVDeathM[Register[ID-1].AgeGroup][CurrYear-StartYear] += 1;}
		}
		else if (p<From4toDead + From4to5){ Register[ID - 1].HIVstageE = 5; }
		else{Register[ID-1].HIVstageE = 4;}
	}
	if(is==5){
		if(Register[ID-1].ARTstage==0){
			if(p<From5toDeadEarly){
				RSApop.SetToDead(ID);
				if(Register[ID-1].SexInd==1){RSApop.HIVDeath[Register[ID-1].AgeGroup][CurrYear-StartYear] += 1;}
				if(Register[ID-1].SexInd==0){RSApop.HIVDeathM[Register[ID-1].AgeGroup][CurrYear-StartYear] += 1;}
			}
			else if (p<From5toDeadEarly + From5to6){ Register[ID - 1].HIVstageE = 6; }
			else{Register[ID-1].HIVstageE = 5;}
		}
		if(Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks<104){
			if(p<From5toDeadLateShort){
				RSApop.SetToDead(ID);
				if(Register[ID-1].SexInd==1){RSApop.HIVDeath[Register[ID-1].AgeGroup][CurrYear-StartYear] += 1;}
				if(Register[ID-1].SexInd==0){RSApop.HIVDeathM[Register[ID-1].AgeGroup][CurrYear-StartYear] += 1;}
			}
			else if (p<From5toDeadLateShort + From5to6){ Register[ID - 1].HIVstageE = 6; }
			else{Register[ID-1].HIVstageE = 5;}
		}
		if(Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks>=104){
			if(p<From5toDeadLateLong){
				RSApop.SetToDead(ID);
				if(Register[ID-1].SexInd==1){RSApop.HIVDeath[Register[ID-1].AgeGroup][CurrYear-StartYear] += 1;}
				if(Register[ID-1].SexInd==0){RSApop.HIVDeathM[Register[ID-1].AgeGroup][CurrYear-StartYear] += 1;}
			}
			else if (p<From5toDeadLateLong + From5to6){ Register[ID - 1].HIVstageE = 6; }
			else{Register[ID-1].HIVstageE = 5;}
		}
	}
	if (is == 6){
		if (p<From6toDead){ 
			RSApop.SetToDead(ID); 
			if(Register[ID-1].SexInd==1){RSApop.HIVDeath[Register[ID-1].AgeGroup][CurrYear-StartYear] += 1;}
			if(Register[ID-1].SexInd==0){RSApop.HIVDeathM[Register[ID-1].AgeGroup][CurrYear-StartYear] += 1;}
		}
		else if (p<From6toDead + From6to5){ Register[ID - 1].HIVstageE = 5; }
		else{ Register[ID - 1].HIVstageE = 6; }
	}
}

void NonHIVtransition::CalcProbCure()
{
	double PropnTreatedPublic, PropnTreatedPrivate;
	int iy;

	if (SexInd == 0){
		PropnTreatedPublic = PropnTreatedPublicM;
		PropnTreatedPrivate = PropnTreatedPrivateM;
	}
	else{
		PropnTreatedPublic = PropnTreatedPublicF;
		PropnTreatedPrivate = PropnTreatedPrivateF;
	}
	iy = CurrYear - StartYear;

	ProbCompleteCure = (PropnTreatedPublic * (PropnPublicUsingSM[iy] * CorrectRxWithSM +
		(1.0 - PropnPublicUsingSM[iy]) * CorrectRxPreSM) * (1.0 - DrugShortage[iy]) +
		PropnTreatedPrivate * (PropnPrivateUsingSM[iy] * CorrectRxWithSM +
		(1.0 - PropnPrivateUsingSM[iy]) * CorrectRxPreSM)) * DrugEff + TradnalEff *
		(1.0 - PropnTreatedPublic - PropnTreatedPrivate);
}

int NonHIVtransition::GetSTDstage(int offset, double r)
{
	int is, stage;
	double CumProb;

	CumProb = 0.0;
	for (is = 0; is<nStates; is++){
		CumProb += PropnByStage[offset][is];
		if (r<CumProb){
			stage = is;
			break;
		}
	}

	return stage;
}

TradnalSTDtransition::TradnalSTDtransition(){}

SyphilisTransition::SyphilisTransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD,
	int ObsCSW, int ObsHH, int ObsANCN, int ObsHHN, int ObsNOART, int ObsONART)
{
	SexInd = Sex;
	nStates = 7;
	ImmuneState = 1;

	ANClogL.Observations = ObsANC;
	FPClogL.Observations = ObsFPC;
	GUDlogL.Observations = ObsGUD;
	CSWlogL.Observations = ObsCSW;
	HouseholdLogL.Observations = ObsHH;
	AntenatalNlogL.Observations = ObsANCN;
	HouseholdNlogL.Observations = ObsHHN;
	NOARTlogL.Observations = ObsNOART;
	ONARTlogL.Observations = ObsONART;
}

void SyphilisTransition::CalcTransitionProbs()
{
	int ia;
	double RxRate, TeenRxRate, RednAsympDurFSW;
	double Adj2ndary; // Adjustment to prob of cure (for primary syphilis) in individuals 
	// with 2ndary syphilis

	if (SexInd == 0){
		RxRate = MaleRxRate;
		TeenRxRate = MaleTeenRxRate;
	}
	else{
		RxRate = FemRxRate;
		TeenRxRate = FemTeenRxRate;
		RednAsympDurFSW = 1.0 / (1.0 + 1.0 / (FSWasympRxRate * ProbCompleteCure * FSWasympCure *
			AveDuration[3]));
	}

	ANCpropnCured = ANCpropnScreened * ANCpropnTreated * DrugEff;
	ProbANCcured[0] = 0;
	if (SexInd == 1){
		for (ia = 1; ia<8; ia++){
			ProbANCcured[ia] = (1.0 - pow(1.0 - SexuallyExpFert[ia - 1], 1.0 / CycleD))*
				ANCpropnCured;
		}
	}

	From1to2 = 1.0 - exp(-(1.0 / AveDuration[0]) * 52.0 / CycleD);
	From2to3 = (1.0 - exp(-(1.0 / AveDuration[1]) * 52.0 / CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-RxRate * ProbCompleteCure * 52.0 / CycleD)));
	From2to3T = (1.0 - exp(-(1.0 / AveDuration[1]) * 52.0 / CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-TeenRxRate * ProbCompleteCure * 52.0 / CycleD)));
	From2to3C = (1.0 - exp(-(1.0 / AveDuration[1]) * 52.0 / CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-FSWRxRate * ProbCompleteCure * 52.0 / CycleD)));
	From2to5 = (1.0 - exp(-RxRate * ProbCompleteCure * 52.0 / CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-(1.0 / AveDuration[1]) * 52.0 / CycleD)));
	From2to0 = From2to5 * PropnSuscepAfterRx;
	From2to5 *= (1.0 - PropnSuscepAfterRx);
	From2to5T = (1.0 - exp(-TeenRxRate * ProbCompleteCure * 52.0 / CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-(1.0 / AveDuration[1]) * 52.0 / CycleD)));
	From2to0T = From2to5T * PropnSuscepAfterRx;
	From2to5T *= (1.0 - PropnSuscepAfterRx);
	From2to5C = (1.0 - exp(-FSWRxRate * ProbCompleteCure * 52.0 / CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-(1.0 / AveDuration[1]) * 52.0 / CycleD)));
	From2to0C = From2to5C * PropnSuscepAfterRx;
	From2to5C *= (1.0 - PropnSuscepAfterRx);
	Adj2ndary = (1.0 - SecondaryRxMult) * (1.0 - SecondaryCureMult);
	From3to4 = (1.0 - exp(-(1.0 / AveDuration[2]) * 52.0 / CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-RxRate * ProbCompleteCure * Adj2ndary * 52.0 / CycleD)));
	From3to4T = (1.0 - exp(-(1.0 / AveDuration[2]) * 52.0 / CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-TeenRxRate * ProbCompleteCure * Adj2ndary * 52.0 / CycleD)));
	From3to4C = (1.0 - exp(-(1.0 / AveDuration[2]) * 52.0 / CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-FSWRxRate * ProbCompleteCure * Adj2ndary * 52.0 / CycleD)));
	From3to5 = (1.0 - exp(-RxRate * ProbCompleteCure * Adj2ndary * 52.0 / CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-(1.0 / AveDuration[2]) * 52.0 / CycleD)));
	From3to5T = (1.0 - exp(-TeenRxRate * ProbCompleteCure * Adj2ndary * 52.0 / CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-(1.0 / AveDuration[2]) * 52.0 / CycleD)));
	From3to5C = (1.0 - exp(-FSWRxRate * ProbCompleteCure * Adj2ndary * 52.0 / CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-(1.0 / AveDuration[2]) * 52.0 / CycleD)));
	From4to6 = 1.0 - exp(-(1.0 / AveDuration[3]) * 52.0 / CycleD);
	From4to6C = 1.0 - exp(-(1.0 / (AveDuration[3] * (1 - RednAsympDurFSW))) * 52.0 / CycleD);
	From5to0 = 1.0 - exp(-(1.0 / AveDuration[4]) * 52.0 / CycleD);
	From6to0 = 1.0 - exp(-(1.0 / AveDuration[5]) * 52.0 / CycleD);
}

double SyphilisTransition::GetTransmProb(int ID)
{
	// ID is the ID of the suscpetible partner, but the SyphilisTransition object is for
	// the opposite sex (the sex of the infected partner).

	int PID1, PID2,   ia, ic;
	int IRisk, PRisk;
	double NoTransmProb, SingleActProb;

	IRisk = Register[ID - 1].RiskGroup;
	NoTransmProb = 1.0;

	// Infection from primary partner
	if (Register[ID - 1].CurrPartners>0){
		PID1 = Register[ID - 1].IDprimary;
		if (Register[PID1 - 1].TPstage == 2 || Register[PID1 - 1].TPstage == 3){
			// Get base probability depending on partnership type
			SingleActProb = TransmProb;
			PRisk = Register[PID1 - 1].RiskGroup;
			if (Register[ID - 1].MarriedInd == 1){
				SingleActProb *= RelTransmLT;
				if (PRisk == 2 && NoBacterialTransm12 == 1){
					SingleActProb = 0.0;
				}
				if (PRisk == 2 && IRisk == 2 && NoBacterialTransm22 == 1){
					SingleActProb = 0.0;
				}
			}
			// Make adjustment for age of individual
			ia = Register[ID - 1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if (SingleActProb>1.0){
				SingleActProb = 1.0;
			}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVIprimary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PVIprimary);
		}
	}

	// Infection from 2ndary partner
	if (Register[ID - 1].CurrPartners == 2){
		PID2 = Register[ID - 1].ID2ndary;
		if (Register[PID2 - 1].TPstage == 2 || Register[PID2 - 1].TPstage == 3){
			SingleActProb = TransmProb;
			// Make adjustment for age of individual
			ia = Register[ID - 1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if (SingleActProb>1.0){
				SingleActProb = 1.0;
			}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVI2ndary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PVI2ndary);
		}
	}

	// Infection from CSW (relevant only to high-risk men)
	if (Register[ID - 1].UVICSW + Register[ID - 1].PVICSW > 0){
		PID2 = Register[ID - 1].IDofCSW;
		if (Register[PID2 - 1].TPstage == 2 || Register[PID2 - 1].TPstage == 3){
			SingleActProb = TransmProb;
			// F-to-M transmission prob is assumed to be the same in commercial sex
			// Make adjustment for age of individual
			ia = Register[ID - 1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if (SingleActProb>1.0){
				SingleActProb = 1.0;
			}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVICSW);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PVICSW);
		}
	}
	int tpp = Register.size();
	// Infection from client (relevant only to CSWs)
	if (Register[ID - 1].FSWind == 1){
		ia = Register[ID - 1].AgeGroup - 2;
	
		for (ic = 0; ic<tpp; ic++){
			if (Register[ic].IDofCSW == ID && (Register[ic].UVICSW + Register[ic].PVICSW) > 0){
				if (Register[ic].TPstage == 2 || Register[ic].TPstage == 3){
					SingleActProb = TransmProbSW;
					SingleActProb *= SuscepIncrease[ia];
					if (SingleActProb>1.0){
						SingleActProb = 1.0;
					}
					NoTransmProb *= pow(1.0 - SingleActProb, Register[ic].UVICSW);
					NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
						Register[ic].PVICSW);
				}
			}
		}
	}

	return 1.0 - NoTransmProb;
}

void SyphilisTransition::GetNewStage(int ID, double p)
{
	int is;

	is = Register[ID - 1].TPstage;
	if (is == 1){
		if (p<From1to2){ Register[ID - 1].TPstageE = 2; }
		else{ Register[ID - 1].TPstageE = 1; }
	}
	if (is == 2){
		if (Register[ID - 1].FSWind == 1){
			if (p<From2to0C){ Register[ID - 1].TPstageE = 0; }
			else if (p<From2to0C + From2to3C){ Register[ID - 1].TPstageE = 3; }
			else if (p<From2to0C + From2to3C + From2to5C){
				Register[ID - 1].TPstageE = 5;
			}
			else{ Register[ID - 1].TPstageE = 2; }
		}
		else if (Register[ID - 1].AgeGroup<4){
			if (p<From2to0T){ Register[ID - 1].TPstageE = 0; }
			else if (p<From2to0T + From2to3T){ Register[ID - 1].TPstageE = 3; }
			else if (p<From2to0T + From2to3T + From2to5T){
				Register[ID - 1].TPstageE = 5;
			}
			else{ Register[ID - 1].TPstageE = 2; }
		}
		else{
			if (p<From2to0){ Register[ID - 1].TPstageE = 0; }
			else if (p<From2to0 + From2to3){ Register[ID - 1].TPstageE = 3; }
			else if (p<From2to0 + From2to3 + From2to5){
				Register[ID - 1].TPstageE = 5;
			}
			else{ Register[ID - 1].TPstageE = 2; }
		}
	}
	if (is == 3){
		if (Register[ID - 1].FSWind == 1){
			if (p<From3to4C){ Register[ID - 1].TPstageE = 4; }
			else if (p<From3to4C + From3to5C){ Register[ID - 1].TPstageE = 5; }
			else{ Register[ID - 1].TPstageE = 3; }
		}
		else if (Register[ID - 1].AgeGroup<4){
			if (p<From3to4T){ Register[ID - 1].TPstageE = 4; }
			else if (p<From3to4T + From3to5T){ Register[ID - 1].TPstageE = 5; }
			else{ Register[ID - 1].TPstageE = 3; }
		}
		/*else{
			if (p<From3to4){ Register[ID - 1].TPstageE = 4; }
			else if (p<From3to4 + From3to5){ Register[ID - 1].TPstageE = 5; }
			else{ Register[ID - 1].TPstageE = 3; }
		}*/
	}
	if (is == 4){
		if (Register[ID - 1].FSWind == 1){
			if (p<From4to6C){ Register[ID - 1].TPstageE = 6; }
			else{ Register[ID - 1].TPstageE = 4; }
		}
		else{
			if (p<From4to6){ Register[ID - 1].TPstageE = 6; }
			else{ Register[ID - 1].TPstageE = 4; }
		}
	}
	if (is == 5){
		if (p<From5to0){ Register[ID - 1].TPstageE = 0; }
		else{ Register[ID - 1].TPstageE = 5; }
	}
	if (is == 6){
		if (p<From6to0){ Register[ID - 1].TPstageE = 0; }
		else{ Register[ID - 1].TPstageE = 6; }
	}
}

HerpesTransition::HerpesTransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD,
	int ObsCSW, int ObsHH, int ObsANCN, int ObsHHN, int ObsNOART, int ObsONART)
{
	SexInd = Sex;
	nStates = 5;

	ANClogL.Observations = ObsANC;
	FPClogL.Observations = ObsFPC;
	GUDlogL.Observations = ObsGUD;
	CSWlogL.Observations = ObsCSW;
	HouseholdLogL.Observations = ObsHH;
	AntenatalNlogL.Observations = ObsANCN;
	HouseholdNlogL.Observations = ObsHHN;
	NOARTlogL.Observations = ObsNOART;
	ONARTlogL.Observations = ObsONART;
}

void HerpesTransition::CalcTransitionProbs()
{
	int is;
	double RxRate, TeenRxRate;

	if (SexInd == 0){
		RxRate = MaleRxRate;
		TeenRxRate = MaleTeenRxRate;
	}
	else{
		RxRate = FemRxRate;
		TeenRxRate = FemTeenRxRate;
	}

	From1to2 = 1.0 - exp(-((1.0 / AveDuration[0]) + RxRate * ProbCompleteCure)*52.0 / CycleD);
	From1to2T = 1.0 - exp(-((1.0 / AveDuration[0]) + TeenRxRate * ProbCompleteCure)*
		52.0 / CycleD);
	From2to3[0] = (1.0 - exp(-RecurrenceRate*52.0 / CycleD)) *
		(1.0 - 0.5*(1.0 - exp(-(1.0 / AveDuration[1])*52.0 / CycleD)));
	for (is = 1; is<6; is++){
		From2to3[is] = 1.0 - exp(-RecurrenceRate * (1.0 + HSVrecurrenceIncrease[is - 1])*
			52.0 / CycleD);
	}
	From3to2 = 1.0 - exp(-((1.0 / AveDuration[2]) + RxRate * ProbCompleteCure)*52.0 / CycleD);
	From3to2T = 1.0 - exp(-((1.0 / AveDuration[2]) + TeenRxRate * ProbCompleteCure)*
		52.0 / CycleD);
	From2to4 = (1 - exp(-(1.0 / AveDuration[1])*52.0 / CycleD)) *
		(1.0 - 0.5*(1.0 - exp(-RecurrenceRate*52.0 / CycleD)));
	if (SexInd == 1){
		From1to2C = 1.0 - exp(-((1.0 / AveDuration[0]) + FSWRxRate * ProbCompleteCure)*
			52.0 / CycleD);
		From3to2C = 1.0 - exp(-((1.0 / AveDuration[2]) + FSWRxRate * ProbCompleteCure)*
			52.0 / CycleD);
	}
}

double HerpesTransition::GetTransmProb(int ID)
{
	// ID is the ID of the suscpetible partner, but the HerpesTransition object is for
	// the opposite sex (the sex of the infected partner).

	int PID1, PID2,  is,  ia, ic;
	int IRisk, PRisk;
	double NoTransmProb, SingleActProb;

	IRisk = Register[ID - 1].RiskGroup;
	NoTransmProb = 1.0;

	// Infection from primary partner
	if (Register[ID - 1].CurrPartners>0){
		PID1 = Register[ID - 1].IDprimary;
		if (Register[PID1 - 1].HSVstage>0){
			// Get base probability depending on partnership type
			SingleActProb = TransmProb;
			PRisk = Register[PID1 - 1].RiskGroup;
			if (Register[ID - 1].MarriedInd == 1){
				SingleActProb *= RelTransmLT;
				if (PRisk == 2 && NoViralTransm12 == 1){
					SingleActProb = 0.0;
				}
				if (PRisk == 2 && IRisk == 2 && NoViralTransm22 == 1){
					SingleActProb = 0.0;
				}
			}
			// Make adjustment for higher HSV-2 shedding when symptomatic
			if (Register[PID1 - 1].HSVstage == 1 || Register[PID1 - 1].HSVstage == 3){
				SingleActProb *= HSVsymptomInfecIncrease;
			}
			// Make adjustment for higher HSV-2 shedding if HIV-positive
			if (Register[PID1 - 1].HIVstage>0){
				is = Register[PID1 - 1].HIVstage;
				SingleActProb *= (1.0 + HSVsheddingIncrease[is - 1]);
			}
			// Make adjustment for age of individual
			ia = Register[ID - 1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if (SingleActProb>1.0){
				SingleActProb = 1.0;
			}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVIprimary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PVIprimary);
		}
	}

	// Infection from 2ndary partner
	if (Register[ID - 1].CurrPartners == 2){
		PID2 = Register[ID - 1].ID2ndary;
		if (Register[PID2 - 1].HSVstage>0){
			SingleActProb = TransmProb;
			// Make adjustment for higher HSV-2 shedding when symptomatic
			if (Register[PID2 - 1].HSVstage == 1 || Register[PID2 - 1].HSVstage == 3){
				SingleActProb *= HSVsymptomInfecIncrease;
			}
			// Make adjustment for higher HSV-2 shedding if HIV-positive
			if (Register[PID2 - 1].HIVstage>0){
				is = Register[PID2 - 1].HIVstage;
				SingleActProb *= (1.0 + HSVsheddingIncrease[is - 1]);
			}
			// Make adjustment for age of individual
			ia = Register[ID - 1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if (SingleActProb>1.0){
				SingleActProb = 1.0;
			}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVI2ndary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PVI2ndary);
		}
	}

	// Infection from CSW (relevant only to high-risk men)
	if (Register[ID - 1].UVICSW + Register[ID - 1].PVICSW > 0){
		PID2 = Register[ID - 1].IDofCSW;
		if (Register[PID2 - 1].HSVstage>0){
			SingleActProb = TransmProb;
			// F-to-M transmission prob is assumed to be the same in commercial sex
			// Make adjustment for higher HSV-2 shedding when symptomatic
			if (Register[PID2 - 1].HSVstage == 1 || Register[PID2 - 1].HSVstage == 3){
				SingleActProb *= HSVsymptomInfecIncrease;
			}
			// Make adjustment for higher HSV-2 shedding if HIV-positive
			if (Register[PID2 - 1].HIVstage>0){
				is = Register[PID2 - 1].HIVstage;
				SingleActProb *= (1.0 + HSVsheddingIncrease[is - 1]);
			}
			// Make adjustment for age of individual
			ia = Register[ID - 1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if (SingleActProb>1.0){
				SingleActProb = 1.0;
			}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVICSW);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PVICSW);
		}
	}
	int tpp = Register.size();
	// Infection from client (relevant only to CSWs)
	if (Register[ID - 1].FSWind == 1){
		ia = Register[ID - 1].AgeGroup - 2;
		
		for (ic = 0; ic<tpp; ic++){
			if (Register[ic].IDofCSW == ID && (Register[ic].UVICSW + Register[ic].PVICSW) > 0){
				if (Register[ic].HSVstage>0){
					SingleActProb = TransmProbSW;
					// Make adjustment for higher HSV-2 shedding when symptomatic
					if (Register[ic].HSVstage == 1 || Register[ic].HSVstage == 3){
						SingleActProb *= HSVsymptomInfecIncrease;
					}
					// Make adjustment for higher HSV-2 shedding if HIV-positive
					if (Register[ic].HIVstage>0){
						is = Register[ic].HIVstage;
						SingleActProb *= (1.0 + HSVsheddingIncrease[is - 1]);
					}
					// Make adjustment for age of individual
					SingleActProb *= SuscepIncrease[ia];
					if (SingleActProb>1.0){
						SingleActProb = 1.0;
					}
					NoTransmProb *= pow(1.0 - SingleActProb, Register[ic].UVICSW);
					NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
						Register[ic].PVICSW);
				}
			}
		}
	}

	return 1.0 - NoTransmProb;
}

void HerpesTransition::GetNewStage(int ID, double p)
{
	int is, is2;

	is = Register[ID - 1].HSVstage;
	if (is == 1){
		if (Register[ID - 1].FSWind == 1){
			if (p<From1to2C){ Register[ID - 1].HSVstageE = 2; }
			else{ Register[ID - 1].HSVstageE = 1; }
		}
		else if (Register[ID - 1].AgeGroup<4){
			if (p<From1to2T){ Register[ID - 1].HSVstageE = 2; }
			else{ Register[ID - 1].HSVstageE = 1; }
		}
		else{
			if (p<From1to2){ Register[ID - 1].HSVstageE = 2; }
			else{ Register[ID - 1].HSVstageE = 1; }
		}
	}
	if (is == 2){
		is2 = Register[ID - 1].HIVstage;
		if (p<From2to3[is2]){ Register[ID - 1].HSVstageE = 3; }
		else if (is2 == 0 && (p<From2to3[0] + From2to4)){
			Register[ID - 1].HSVstageE = 4;
		}
		else{
			Register[ID - 1].HSVstageE = 2;
		}
	}
	if (is == 3){
		if (Register[ID - 1].FSWind == 1){
			if (p<From3to2C){ Register[ID - 1].HSVstageE = 2; }
			else{ Register[ID - 1].HSVstageE = 3; }
		}
		else if (Register[ID - 1].AgeGroup<4){
			if (p<From3to2T){ Register[ID - 1].HSVstageE = 2; }
			else{ Register[ID - 1].HSVstageE = 3; }
		}
		else{
			if (p<From3to2){ Register[ID - 1].HSVstageE = 2; }
			else{ Register[ID - 1].HSVstageE = 3; }
		}
	}
	if (is == 4){
		Register[ID - 1].HSVstageE = 4;
	}
}

OtherSTDtransition::OtherSTDtransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD,
	int ObsCSW, int ObsHH, int ObsANCN, int ObsHHN, int ObsNOART, int ObsONART)
{
	SexInd = Sex;
	nStates = 4;
	ImmuneState = 3;

	ANClogL.Observations = ObsANC;
	FPClogL.Observations = ObsFPC;
	GUDlogL.Observations = ObsGUD;
	CSWlogL.Observations = ObsCSW;
	HouseholdLogL.Observations = ObsHH;
	AntenatalNlogL.Observations = ObsANCN;
	HouseholdNlogL.Observations = ObsHHN;
	NOARTlogL.Observations = ObsNOART;
	ONARTlogL.Observations = ObsONART;
}

void OtherSTDtransition::CalcTransitionProbs()
{
	double RxRate, TeenRxRate; // , RednAsympDurFSW;

	if (SexInd == 0){
		RxRate = MaleRxRate;
		TeenRxRate = MaleTeenRxRate;
	}
	else{
		RxRate = FemRxRate;
		TeenRxRate = FemTeenRxRate;
		//RednAsympDurFSW = 1.0/(1.0 + 1.0/(FSWasympRxRate * ProbCompleteCure * FSWasympCure *
		//	AveDuration[1]));
	}

	From1to0 = 1.0 - exp(-(((1.0 - PropnImmuneAfterSR) / AveDuration[0]) + RxRate *
		ProbCompleteCure * (1.0 - PropnImmuneAfterRx))*52.0 / CycleD);
	From1to0T = 1.0 - exp(-(((1.0 - PropnImmuneAfterSR) / AveDuration[0]) + TeenRxRate *
		ProbCompleteCure * (1.0 - PropnImmuneAfterRx))*52.0 / CycleD);
	From2to0 = 1.0 - exp(-((1.0 - PropnImmuneAfterSR) / AveDuration[1])*52.0 / CycleD);
	From1to3 = 1.0 - exp(-((PropnImmuneAfterSR / AveDuration[0]) + RxRate *
		ProbCompleteCure * PropnImmuneAfterRx)*52.0 / CycleD);
	From1to3T = 1.0 - exp(-((PropnImmuneAfterSR / AveDuration[0]) + TeenRxRate *
		ProbCompleteCure * PropnImmuneAfterRx)*52.0 / CycleD);
	From2to3 = 1.0 - exp(-(PropnImmuneAfterSR / AveDuration[1])*52.0 / CycleD);
	From3to0 = 1.0 - exp(-(1.0 / AveDuration[2])*52.0 / CycleD);
	// Transition probabilities for sex workers
	if (SexInd == 1){
		From1to0C = 1.0 - exp(-(((1.0 - PropnImmuneAfterSR) / AveDuration[0]) + FSWRxRate *
			ProbCompleteCure * (1.0 - PropnImmuneAfterRx))*52.0 / CycleD);
		From2to0C = 1.0 - exp(-(((1.0 - PropnImmuneAfterSR) / AveDuration[1]) + FSWasympRxRate *
			ProbCompleteCure * FSWasympCure * (1.0 - PropnImmuneAfterRx))*52.0 / CycleD);
		From1to3C = 1.0 - exp(-((PropnImmuneAfterSR / AveDuration[0]) + FSWRxRate *
			ProbCompleteCure * PropnImmuneAfterRx)*52.0 / CycleD);
		From2to3C = 1.0 - exp(-((PropnImmuneAfterSR / AveDuration[1]) + FSWasympRxRate *
			ProbCompleteCure * FSWasympCure * PropnImmuneAfterRx)*52.0 / CycleD);
	}
}

double OtherSTDtransition::GetTransmProb(int ID, int STD)
{
	// ID is the ID of the suscpetible partner, but the OtherSTDtransition object is for
	// the opposite sex (the sex of the infected partner).
	// STD indicator is 1 for NG, 2 for CT, 3 for TV and 4 for HD

	int PID1, PID2,  ia, ic;
	int IRisk, PRisk;
	double NoTransmProb, SingleActProb;

	IRisk = Register[ID - 1].RiskGroup;
	NoTransmProb = 1.0;

	// Infection from primary partner
	if (Register[ID - 1].CurrPartners>0){
		PID1 = Register[ID - 1].IDprimary;
		if ((STD == 1 && Register[PID1 - 1].NGstage>0 && Register[PID1 - 1].NGstage<3) ||
			(STD == 2 && Register[PID1 - 1].CTstage>0 && Register[PID1 - 1].CTstage<3) ||
			(STD == 3 && Register[PID1 - 1].TVstage>0 && Register[PID1 - 1].TVstage<3) ||
			(STD == 4 && Register[PID1 - 1].HDstage>0 && Register[PID1 - 1].HDstage<3)){
			// Get base probability depending on partnership type
			SingleActProb = TransmProb;
			PRisk = Register[PID1 - 1].RiskGroup;
			if (Register[ID - 1].MarriedInd == 1){
				SingleActProb *= RelTransmLT;
				if (PRisk == 2 && NoBacterialTransm12 == 1){
					SingleActProb = 0.0;
				}
				if (PRisk == 2 && IRisk == 2 && NoBacterialTransm22 == 1){
					SingleActProb = 0.0;
				}
			}
			// Make adjustment for age of individual
			ia = Register[ID - 1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if (SingleActProb>1.0){
				SingleActProb = 1.0;
			}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVIprimary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PVIprimary);
		}
	}

	// Infection from 2ndary partner
	if (Register[ID - 1].CurrPartners == 2){
		PID2 = Register[ID - 1].ID2ndary;
		if ((STD == 1 && Register[PID2 - 1].NGstage>0 && Register[PID2 - 1].NGstage<3) ||
			(STD == 2 && Register[PID2 - 1].CTstage>0 && Register[PID2 - 1].CTstage<3) ||
			(STD == 3 && Register[PID2 - 1].TVstage>0 && Register[PID2 - 1].TVstage<3) ||
			(STD == 4 && Register[PID2 - 1].HDstage>0 && Register[PID2 - 1].HDstage<3)){
			SingleActProb = TransmProb;
			// Make adjustment for age of individual
			ia = Register[ID - 1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if (SingleActProb>1.0){
				SingleActProb = 1.0;
			}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVI2ndary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PVI2ndary);
		}
	}

	// Infection from CSW (relevant only to high-risk men)
	if (Register[ID - 1].UVICSW + Register[ID - 1].PVICSW > 0){
		PID2 = Register[ID - 1].IDofCSW;
		if ((STD == 1 && Register[PID2 - 1].NGstage>0 && Register[PID2 - 1].NGstage<3) ||
			(STD == 2 && Register[PID2 - 1].CTstage>0 && Register[PID2 - 1].CTstage<3) ||
			(STD == 3 && Register[PID2 - 1].TVstage>0 && Register[PID2 - 1].TVstage<3) ||
			(STD == 4 && Register[PID2 - 1].HDstage>0 && Register[PID2 - 1].HDstage<3)){
			SingleActProb = TransmProb;
			// F-to-M transmission prob is assumed to be the same in commercial sex
			// Make adjustment for age of individual
			ia = Register[ID - 1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if (SingleActProb>1.0){
				SingleActProb = 1.0;
			}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVICSW);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PVICSW);
		}
	}

	int  tpp = Register.size();
	// Infection from client (relevant only to CSWs)
	if (Register[ID - 1].FSWind == 1){
		ia = Register[ID - 1].AgeGroup - 2;
		
		for (ic = 0; ic<tpp; ic++){
			if (Register[ic].IDofCSW == ID && (Register[ic].UVICSW + Register[ic].PVICSW) > 0){
				if ((STD == 1 && Register[ic].NGstage>0 && Register[ic].NGstage<3) ||
					(STD == 2 && Register[ic].CTstage>0 && Register[ic].CTstage<3) ||
					(STD == 3 && Register[ic].TVstage>0 && Register[ic].TVstage<3) ||
					(STD == 4 && Register[ic].HDstage>0 && Register[ic].HDstage<3)){
					SingleActProb = TransmProbSW;
					SingleActProb *= SuscepIncrease[ia];
					if (SingleActProb>1.0){
						SingleActProb = 1.0;
					}
					NoTransmProb *= pow(1.0 - SingleActProb, Register[ic].UVICSW);
					NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
						Register[ic].PVICSW);
				}
			}
		}
	}

	return 1.0 - NoTransmProb;
}

void OtherSTDtransition::GetNewStage(int ID, double p, int STD)
{
	int StartStage, EndStage;

	if (STD == 1){ StartStage = Register[ID - 1].NGstage; }
	if (STD == 2){ StartStage = Register[ID - 1].CTstage; }
	if (STD == 3){ StartStage = Register[ID - 1].TVstage; }
	if (STD == 4){ StartStage = Register[ID - 1].HDstage; }

	if (StartStage == 1){
		if (Register[ID - 1].FSWind == 1){
			if (p<From1to0C){ EndStage = 0; }
			else if (p<From1to0C + From1to3C){ EndStage = 3; }
			else{ EndStage = 1; }
		}
		else if (Register[ID - 1].AgeGroup<4){
			if (p<From1to0T){ EndStage = 0; }
			else if (p<From1to0T + From1to3T){ EndStage = 3; }
			else{ EndStage = 1; }
		}
		else{
			if (p<From1to0){ EndStage = 0; }
			else if (p<From1to0 + From1to3){ EndStage = 3; }
			else{ EndStage = 1; }
		}
	}
	if (StartStage == 2){
		if (Register[ID - 1].FSWind == 1){
			if (p<From2to0C){ EndStage = 0; }
			else if (p<From2to0C + From2to3C){ EndStage = 3; }
			else{ EndStage = 2; }
		}
		else{
			if (p<From2to0){ EndStage = 0; }
			else if (p<From2to0 + From2to3){ EndStage = 3; }
			else{ EndStage = 2; }
		}
	}
	if (StartStage == 3){
		if (p<From3to0){ EndStage = 0; }
		else{ EndStage = 3; }
	}

	if (STD == 1){ Register[ID - 1].NGstageE = EndStage; }
	if (STD == 2){ Register[ID - 1].CTstageE = EndStage; }
	if (STD == 3){ Register[ID - 1].TVstageE = EndStage; }
	if (STD == 4){ Register[ID - 1].HDstageE = EndStage; }
}

NonSTDtransition::NonSTDtransition(){}

void NonSTDtransition::CalcProbPartialCure()
{
	int iy;

	iy = CurrYear - StartYear;

	ProbPartialCure = (PropnTreatedPublicF * (PropnPublicUsingSM[iy] * CorrectRxWithSM +
		(1.0 - PropnPublicUsingSM[iy]) * CorrectRxPreSM) * (1.0 - DrugShortage[iy]) +
		PropnTreatedPrivateF * (PropnPrivateUsingSM[iy] * CorrectRxWithSM +
		(1.0 - PropnPrivateUsingSM[iy]) * CorrectRxPreSM)) * DrugPartialEff + TradnalEff *
		(1.0 - PropnTreatedPublicF - PropnTreatedPrivateF);
}

BVtransition::BVtransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD,
	int ObsCSW, int ObsHH, int ObsANCN, int ObsHHN, int ObsNOART, int ObsONART)
{
	SexInd = Sex;
	nStates = 4;

	ANClogL.Observations = ObsANC;
	FPClogL.Observations = ObsFPC;
	GUDlogL.Observations = ObsGUD;
	CSWlogL.Observations = ObsCSW;
	HouseholdLogL.Observations = ObsHH;
	AntenatalNlogL.Observations = ObsANCN;
	HouseholdNlogL.Observations = ObsHHN;
	NOARTlogL.Observations = ObsNOART;
	ONARTlogL.Observations = ObsONART;
}

void BVtransition::CalcTransitionProbs()
{
	int ij;
	double RednAsympDurFSW;

	RednAsympDurFSW = 1.0 / (1.0 + 1.0 / (FSWasympRxRate * (ProbCompleteCure + ProbPartialCure) *
		FSWasympCure * AveDuration[3]));
	From1to2 = 1.0 - exp(-CtsTransition[0][1] * 52.0 / CycleD);
	From2to1ind = 1.0 - exp(-CtsTransition[1][0] * 52.0 / CycleD);
	From3to1 = (1.0 - exp(-(CtsTransition[2][0] + FemRxRate * ProbCompleteCure) * 52.0 /
		CycleD)) * (1.0 - 0.5*(1.0 - exp(-(CtsTransition[2][1] + FemRxRate *
		ProbPartialCure) * 52.0 / CycleD)));
	From3to1T = (1.0 - exp(-(CtsTransition[2][0] + FemTeenRxRate * ProbCompleteCure) * 52.0 /
		CycleD)) * (1.0 - 0.5*(1.0 - exp(-(CtsTransition[2][1] + FemTeenRxRate *
		ProbPartialCure) * 52.0 / CycleD)));
	From3to1C = (1.0 - exp(-(CtsTransition[2][0] + FSWRxRate * ProbCompleteCure) * 52.0 /
		CycleD)) * (1.0 - 0.5*(1.0 - exp(-(CtsTransition[2][1] + FSWRxRate *
		ProbPartialCure) * 52.0 / CycleD)));
	From3to2 = (1.0 - exp(-(CtsTransition[2][1] + FemRxRate * ProbPartialCure) * 52.0 /
		CycleD)) * (1.0 - 0.5*(1.0 - exp(-(CtsTransition[2][0] + FemRxRate *
		ProbCompleteCure) * 52.0 / CycleD)));
	From3to2T = (1.0 - exp(-(CtsTransition[2][1] + FemTeenRxRate * ProbPartialCure) * 52.0 /
		CycleD)) * (1.0 - 0.5*(1.0 - exp(-(CtsTransition[2][0] + FemTeenRxRate *
		ProbCompleteCure) * 52.0 / CycleD)));
	From3to2C = (1.0 - exp(-(CtsTransition[2][1] + FSWRxRate * ProbPartialCure) * 52.0 /
		CycleD)) * (1.0 - 0.5*(1.0 - exp(-(CtsTransition[2][0] + FSWRxRate *
		ProbCompleteCure) * 52.0 / CycleD)));
	From4to1 = (1.0 - exp(-CtsTransition[3][0] * 52.0 / CycleD)) *
		(1.0 - 0.5 * (1.0 - exp(-CtsTransition[3][1] * 52.0 / CycleD)));
	From4to1C = (1.0 - exp(-(CtsTransition[3][0] / (1.0 - RednAsympDurFSW))*52.0 / CycleD)) *
		(1.0 - 0.5 * (1.0 - exp(-(CtsTransition[3][1] / (1.0 - RednAsympDurFSW))*52.0 / CycleD)));
	From4to2 = (1.0 - exp(-CtsTransition[3][1] * 52.0 / CycleD)) *
		(1.0 - 0.5 * (1.0 - exp(-CtsTransition[3][0] * 52.0 / CycleD)));
	From4to2C = (1.0 - exp(-(CtsTransition[3][1] / (1.0 - RednAsympDurFSW))*52.0 / CycleD)) *
		(1.0 - 0.5 * (1.0 - exp(-(CtsTransition[3][0] / (1.0 - RednAsympDurFSW))*52.0 / CycleD)));

	From2to3ind[0] = 1.0 - exp(-CtsTransition[1][2] * (1.0 - IncidenceMultNoPartners)*52.0 /
		CycleD);
	From2to3ind[1] = 1.0 - exp(-CtsTransition[1][2] * 52.0 / CycleD);
	From2to3ind[2] = 1.0 - exp(-CtsTransition[1][2] * (1.0 + IncidenceMultTwoPartners)*52.0 /
		CycleD);
	From2to4ind[0] = 1.0 - exp(-CtsTransition[1][3] * (1.0 - IncidenceMultNoPartners)*52.0 /
		CycleD);
	From2to4ind[1] = 1.0 - exp(-CtsTransition[1][3] * 52.0 / CycleD);
	From2to4ind[2] = 1.0 - exp(-CtsTransition[1][3] * (1.0 + IncidenceMultTwoPartners)*52.0 /
		CycleD);
	for (ij = 0; ij<3; ij++){
		From2to1dep[ij] = From2to1ind * (1.0 - 0.5 * (From2to3ind[ij] + From2to4ind[ij]) +
			From2to3ind[ij] * From2to4ind[ij] / 3.0);
		From2to3dep[ij] = From2to3ind[ij] * (1.0 - 0.5 * (From2to1ind + From2to4ind[ij]) +
			From2to1ind * From2to4ind[ij] / 3.0);
		From2to4dep[ij] = From2to4ind[ij] * (1.0 - 0.5 * (From2to1ind + From2to3ind[ij]) +
			From2to1ind * From2to3ind[ij] / 3.0);
	}
}

void BVtransition::GetNewStage(int ID, double p)
{
	int partners;

	if (Register[ID - 1].BVstage == 1){
		if (p<From1to2){ Register[ID - 1].BVstageE = 2; }
		else{ Register[ID - 1].BVstageE = 1; }
	}
	if (Register[ID - 1].BVstage == 2){
		partners = Register[ID - 1].CurrPartners;
		if (p<From2to1dep[partners]){ Register[ID - 1].BVstageE = 1; }
		else if (p<From2to1dep[partners] + From2to3dep[partners]){
			Register[ID - 1].BVstageE = 3;
		}
		else if (p<From2to1dep[partners] + From2to3dep[partners] + From2to4dep[partners]){
			Register[ID - 1].BVstageE = 4;
		}
		else{ Register[ID - 1].BVstageE = 2; }
	}
	if (Register[ID - 1].BVstage == 3){
		if (Register[ID - 1].FSWind == 1){
			if (p<From3to1C){ Register[ID - 1].BVstageE = 1; }
			else if (p<From3to1C + From3to2C){ Register[ID - 1].BVstageE = 2; }
			else{ Register[ID - 1].BVstageE = 3; }
		}
		else if (Register[ID - 1].AgeGroup<4){
			if (p<From3to1T){ Register[ID - 1].BVstageE = 1; }
			else if (p<From3to1T + From3to2T){ Register[ID - 1].BVstageE = 2; }
			else{ Register[ID - 1].BVstageE = 3; }
		}
		else{
			if (p<From3to1){ Register[ID - 1].BVstageE = 1; }
			else if (p<From3to1 + From3to2){ Register[ID - 1].BVstageE = 2; }
			else{ Register[ID - 1].BVstageE = 3; }
		}
	}
	if (Register[ID - 1].BVstage == 4){
		if (Register[ID - 1].FSWind == 1){
			if (p<From4to1C){ Register[ID - 1].BVstageE = 1; }
			else if (p<From4to1C + From4to2C){ Register[ID - 1].BVstageE = 2; }
			else{ Register[ID - 1].BVstageE = 4; }
		}
		else{
			if (p<From4to1){ Register[ID - 1].BVstageE = 1; }
			else if (p<From4to1 + From4to2){ Register[ID - 1].BVstageE = 2; }
			else{ Register[ID - 1].BVstageE = 4; }
		}
	}
}

VCtransition::VCtransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD,
	int ObsCSW, int ObsHH, int ObsANCN, int ObsHHN, int ObsNOART, int ObsONART)
{
	SexInd = Sex;
	nStates = 3;

	ANClogL.Observations = ObsANC;
	FPClogL.Observations = ObsFPC;
	GUDlogL.Observations = ObsGUD;
	CSWlogL.Observations = ObsCSW;
	HouseholdLogL.Observations = ObsHH;
	AntenatalNlogL.Observations = ObsANCN;
	HouseholdNlogL.Observations = ObsHHN;
	NOARTlogL.Observations = ObsNOART;
	ONARTlogL.Observations = ObsONART;
}

void VCtransition::CalcTransitionProbs()
{
	int ia, is;
	double RednAsympDurFSW; // Although it's not strictly necessary to include this (since
	// we set it to 0), we may want to change this in future, so I
	// have included it.

	RednAsympDurFSW = 0.0;
	From1to2 = (1.0 - exp(-RecurrenceRate * 52.0 / CycleD)) *
		(1.0 - 0.5 * (1.0 - exp(-(1.0 / AveDuration[0]) * 52.0 / CycleD)));
	From1to2C = (1.0 - exp(-RecurrenceRate * 52.0 / CycleD)) * (1.0 - 0.5 *
		(1.0 - exp(-(1.0 / (AveDuration[0] * (1.0 - RednAsympDurFSW))) * 52.0 / CycleD)));
	From1to0 = (1.0 - exp(-(1.0 / AveDuration[0]) * 52.0 / CycleD)) *
		(1.0 - 0.5 * (1.0 - exp(-RecurrenceRate * 52.0 / CycleD)));
	From1to0C = (1.0 - exp(-(1.0 / (AveDuration[0] * (1.0 - RednAsympDurFSW))) * 52.0 / CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-RecurrenceRate * 52.0 / CycleD)));
	From2to1 = (1.0 - exp(-((1.0 / AveDuration[1]) + FemRxRate * ProbPartialCure)*52.0 /
		CycleD))*(1.0 - 0.5 * (1.0 - exp(-FemRxRate * ProbCompleteCure * 52.0 / CycleD)));
	From2to1T = (1.0 - exp(-((1.0 / AveDuration[1]) + FemTeenRxRate * ProbPartialCure)*52.0 /
		CycleD))*(1.0 - 0.5 * (1.0 - exp(-FemTeenRxRate * ProbCompleteCure * 52.0 / CycleD)));
	From2to1C = (1.0 - exp(-((1.0 / AveDuration[1]) + FSWRxRate * ProbPartialCure)*52.0 /
		CycleD))*(1.0 - 0.5 * (1.0 - exp(-FSWRxRate * ProbCompleteCure * 52.0 / CycleD)));
	From2to0 = (1.0 - exp(-FemRxRate * ProbCompleteCure * 52.0 / CycleD)) * (1.0 - 0.5 *
		(1.0 - exp(-((1.0 / AveDuration[1]) + FemRxRate * ProbPartialCure)*52.0 / CycleD)));
	From2to0T = (1.0 - exp(-FemTeenRxRate * ProbCompleteCure * 52.0 / CycleD)) * (1.0 - 0.5 *
		(1.0 - exp(-((1.0 / AveDuration[1]) + FemTeenRxRate * ProbPartialCure)*52.0 / CycleD)));
	From2to0C = (1.0 - exp(-FSWRxRate * ProbCompleteCure * 52.0 / CycleD)) * (1.0 - 0.5 *
		(1.0 - exp(-((1.0 / AveDuration[1]) + FSWRxRate * ProbPartialCure)*52.0 / CycleD)));

	for (ia = 0; ia<7; ia++){
		From0to1[ia][0] = 1.0 - exp(-Incidence * (FertilityTable[ia][0] /
			FertilityTable[0][0])*52.0 / CycleD);
	}
	if (HIVind == 1){
		for (ia = 0; ia<7; ia++){
			for (is = 1; is<6; is++){
				From0to1[ia][is] = 1.0 - exp(-Incidence * (FertilityTable[ia][0] /
					FertilityTable[0][0]) * (1.0 + IncidenceIncrease[is - 1]) * 52.0 / CycleD);
			}
		}
	}
}

void VCtransition::GetNewStage(int ID, double p)
{
	int ia, is;

	if (Register[ID - 1].VCstage == 0){
		ia = Register[ID - 1].AgeGroup - 3;
		is = Register[ID - 1].HIVstage;
		if (ia >= 0 && ia<7){
			if (p<From0to1[ia][is]){ Register[ID - 1].VCstageE = 1; }
			else{ Register[ID - 1].VCstageE = 0; }
		}
		else{ Register[ID - 1].VCstageE = 0; }
	}
	if (Register[ID - 1].VCstage == 1){
		if (Register[ID - 1].FSWind == 1){
			if (p<From1to0C){ Register[ID - 1].VCstageE = 0; }
			else if (p<From1to0C + From1to2C){ Register[ID - 1].VCstageE = 2; }
			else{ Register[ID - 1].VCstageE = 1; }
		}
		else{
			if (p<From1to0){ Register[ID - 1].VCstageE = 0; }
			else if (p<From1to0 + From1to2){ Register[ID - 1].VCstageE = 2; }
			else{ Register[ID - 1].VCstageE = 1; }
		}
	}
	if (Register[ID - 1].VCstage == 2){
		if (Register[ID - 1].FSWind == 1){
			if (p<From2to0C){ Register[ID - 1].VCstageE = 0; }
			else if (p<From2to0C + From2to1C){ Register[ID - 1].VCstageE = 1; }
			else{ Register[ID - 1].VCstageE = 2; }
		}
		else if (Register[ID - 1].AgeGroup<4){
			if (p<From2to0T){ Register[ID - 1].VCstageE = 0; }
			else if (p<From2to0T + From2to1T){ Register[ID - 1].VCstageE = 1; }
			else{ Register[ID - 1].VCstageE = 2; }
		}
		else{
			if (p<From2to0){ Register[ID - 1].VCstageE = 0; }
			else if (p<From2to0 + From2to1){ Register[ID - 1].VCstageE = 1; }
			else{ Register[ID - 1].VCstageE = 2; }
		}
	}
}

PaedHIV::PaedHIV(){}

void PaedHIV::GetNewStage(int ID, double p)
{
	double t1, t2;
	double AIDSprog, From2to4, From4toDead, From5toDead;

	if (Register[ID - 1].HIVstage == 2){
		t2 = CurrYear + 0.5 + 1.0 * (STDcycleCount / CycleD) +
			1.0 * (BehavCycleCount - 1.0) / CycleS - Register[ID - 1].DOB;
		t1 = t2 - 1.0 / CycleD;
		if (t1<0.0){ t1 = 0.0; }
		AIDSprog = pow(0.5, pow(t2 / PreAIDSmedian, PreAIDSshape) -
			pow(t1 / PreAIDSmedian, PreAIDSshape));
		From2to4 = AIDSprog * (1.0 - HAARTaccess[CurrYear - StartYear]);
		if (p < From2to4){
			Register[ID - 1].HIVstageE = 4;
		}
		else if (p < AIDSprog){
			Register[ID - 1].HIVstageE = 5;
		}
		else{
			Register[ID - 1].HIVstageE = 2;
		}
	}
	if (Register[ID - 1].HIVstage == 4){
		From4toDead = 1.0 - exp(-1.0 / (MeanAIDSsurvival * CycleD));
		if (p < From4toDead){ RSApop.SetToDead(ID); }
		else{ Register[ID - 1].HIVstageE = 4; }
	}
	if (Register[ID - 1].HIVstage == 5){
		From5toDead = 1.0 - exp(-1.0 / (AveYrsOnART * CycleD));
		if (p < From5toDead){ RSApop.SetToDead(ID); }
		else{ Register[ID - 1].HIVstageE = 5; }
	}
}

Child::Child(int Sex)
{
	SexInd = Sex;
}

void Child::UpdateMort()
{
	int ia, yr;

	yr = CurrYear - StartYear;

	if (SexInd == 0){
		MortProb1st6m = InfantMort1st6mM[yr];
	}
	else{
		MortProb1st6m = InfantMort1st6mF[yr];
	}
	Perinatal.MortProb1st6m = MortProb1st6m;
	Breastmilk.MortProb1st6m = MortProb1st6m;
	for (ia = 0; ia<15; ia++){
		if (SexInd == 0){
			NonAIDSmort[ia] = ChildMortM[ia][yr];
		}
		else{
			NonAIDSmort[ia] = ChildMortF[ia][yr];
		}
		Perinatal.NonAIDSmort[ia] = NonAIDSmort[ia];
		Breastmilk.NonAIDSmort[ia] = NonAIDSmort[ia];
	}
}

Indiv::Indiv()
{
	HIVstage = 0;
	DateInfect = 0.0;
	DesiredNewPartners = 0.0;
	NewStatus = 0;
	UVIprimary = 0;
	PVIprimary = 0;
	UVI2ndary = 0;
	PVI2ndary = 0;
	UVICSW = 0;
	PVICSW = 0;
	IDofCSW = 0;
	SuscepHIVadj = 1.0;

	totUVIprimary = 0; //Total # unprotected acts of vaginal intercourse with primary partners
	totPVIprimary = 0; //Total # protected acts of vaginal intercourse with primary partners
	totUVI2ndary = 0; //Total # unprotected acts of vaginal intercourse with secondary partners
	totPVI2ndary = 0; //Total # protected acts of vaginal intercourse with secondary partners
	totUVICSW = 0; // Total # unprotected acts of vaginal intercourse with CSWs
	totPVICSW = 0; // Total # protected acts of vaginal intercourse with CSWs
	SexDebutAge = 0;

	reason=0;
	repeat=0;
	HPVrepeat=0;
	DiagnosedCC=0;
	
}

const vector<int> Indiv::allhpv = {0,1,2,3,4,5,6,7,8,9,10,11,12};
const vector<int> Indiv::hpv1618 = {0,1};
const vector<int> Indiv::hpv161845 = {0,1,6};

const vector<int> Indiv::lsil = {2};
const vector<int> Indiv::hsil = {3, 4};
const vector<int> Indiv::cc_un12 = {5, 8};
const vector<int> Indiv::cc_un34 = {9, 10};
const vector<int> Indiv::cc_diag = {11, 12, 13, 14};
const vector<int> Indiv::recover = {15};

int Indiv::SelectEvent(double rnd)
{
	int ia, ii, EventType;
	double AcquireS1, AcquireS2, Marry1, Marry2, Divorce, LoseS1, LoseS2;
	// Dependent versions of above rates:
	double DAcquireS1, DAcquireS2, DMarry1, DMarry2, DDivorce, DLoseS1, DLoseS2;
	int prisk, temp;
	double CumEventProb[8];

	AcquireS1 = 0.0;
	AcquireS2 = 0.0;
	Marry1 = 0.0;
	Marry2 = 0.0;
	Divorce = 0.0;
	LoseS1 = 0.0;
	LoseS2 = 0.0;
	DAcquireS1 = 0.0;
	DAcquireS2 = 0.0;
	DMarry1 = 0.0;
	DMarry2 = 0.0;
	DDivorce = 0.0;
	DLoseS1 = 0.0;
	DLoseS2 = 0.0;

	// Calculate rates at which new ST partners are acquired
	if (DesiredNewPartners>0.0 && VirginInd == 0){
		if (SexInd == 0){
			AcquireS1 = DesiredNewPartners * DesiredPartnerRiskM[RiskGroup - 1][0] *
				AdjSTrateM[RiskGroup - 1][0] / CycleS;
			AcquireS2 = DesiredNewPartners * DesiredPartnerRiskM[RiskGroup - 1][1] *
				AdjSTrateM[RiskGroup - 1][1] / CycleS;
		}
		else{
			AcquireS1 = DesiredNewPartners * DesiredPartnerRiskF[RiskGroup - 1][0] *
				AdjSTrateF[RiskGroup - 1][0] / CycleS;
			AcquireS2 = DesiredNewPartners * DesiredPartnerRiskF[RiskGroup - 1][1] *
				AdjSTrateF[RiskGroup - 1][1] / CycleS;
		}
	}
	else if (DesiredNewPartners>0.0 && VirginInd == 1){ // No adjustment factor applied
		if (SexInd == 0){
			AcquireS1 = DesiredNewPartners * DesiredPartnerRiskM[RiskGroup - 1][0] / CycleS;
			AcquireS2 = DesiredNewPartners * DesiredPartnerRiskM[RiskGroup - 1][1] / CycleS;
		}
		else{
			AcquireS1 = DesiredNewPartners * DesiredPartnerRiskF[RiskGroup - 1][0] / CycleS;
			AcquireS2 = DesiredNewPartners * DesiredPartnerRiskF[RiskGroup - 1][1] / CycleS;
		}
	}

	// Calculate rates at which marriage occurs
	if (CurrPartners>0.0 && MarriedInd == 0){
		ia = AgeGroup - 2;
		temp = SexInd * 2 + RiskGroup - 1;
		if (Register[IDprimary - 1].MarriedInd == 0 &&
			Register[IDprimary - 1].NewStatus == 0){
			prisk = Register[IDprimary - 1].RiskGroup - 1;
			Marry1 = MarriageRate[prisk][temp] * AgeEffectMarriage[ia][temp] / CycleS;
			if (SexInd == 0){
				Marry1 *= AdjLTrateM[RiskGroup - 1][prisk];
			}
			else{
				Marry1 *= AdjLTrateF[RiskGroup - 1][prisk];
			}
		}
		if (CurrPartners == 2.0){
			if (Register[ID2ndary - 1].MarriedInd == 0 &&
				Register[ID2ndary - 1].NewStatus == 0){
				prisk = Register[ID2ndary - 1].RiskGroup - 1;
				Marry2 = MarriageRate[prisk][temp] * AgeEffectMarriage[ia][temp] / CycleS;
				if (SexInd == 0){
					Marry2 *= AdjLTrateM[RiskGroup - 1][prisk];
				}
				else{
					Marry2 *= AdjLTrateF[RiskGroup - 1][prisk];
				}
			}
		}
	}

	// Calculate rate at which divorce occurs
	if (MarriedInd == 1){
		if (Register[IDprimary - 1].NewStatus == 0){ // Messy
			Divorce = LTseparation[AgeGroup - 2][SexInd] / CycleS;
		}
	}

	// Calculate rate at which primary partnership is terminated if it is ST
	// (Note that the code differs from the corresponding function in TSHISA, where we weren't
	// worried about the order in which partnerships occur).
	if (CurrPartners>0 && MarriedInd == 0){
		prisk = Register[IDprimary - 1].RiskGroup - 1; // Messy
		if (Register[IDprimary - 1].NewStatus == 1){
			LoseS1 = 0.0;
		}
		else if (SexInd == 0){
			LoseS1 = PartnerRateAdj / MeanDurSTrel[RiskGroup - 1][prisk];
		}
		else{
			LoseS1 = PartnerRateAdj / MeanDurSTrel[prisk][RiskGroup - 1];
		}
		LoseS1 = LoseS1 / CycleS;
	}

	// Calculate rate at which 2ndary partnership is terminated
	// (Note that the code differs from the corresponding function in TSHISA, where we weren't
	// worried about the order in which partnerships occur).
	if (CurrPartners == 2){
		prisk = Register[ID2ndary - 1].RiskGroup - 1; // Messy
		if (Register[ID2ndary - 1].NewStatus == 1){
			LoseS2 = 0.0;
		}
		else if (SexInd == 0){
			LoseS2 = PartnerRateAdj / MeanDurSTrel[RiskGroup - 1][prisk];
		}
		else{
			LoseS2 = PartnerRateAdj / MeanDurSTrel[prisk][RiskGroup - 1];
		}
		LoseS2 = LoseS2 / CycleS;
	}

	// Calculate dependent probabilities of different events
	if (CurrPartners == 0){
		DAcquireS1 = ConvertToDependent2(AcquireS1, AcquireS2, 1);
		DAcquireS2 = ConvertToDependent2(AcquireS1, AcquireS2, 2);
	}
	if (CurrPartners == 1 && MarriedInd == 0){
		if (RiskGroup == 1){
			DAcquireS1 = ConvertToDependent4(AcquireS1, AcquireS2, LoseS1, Marry1, 1);
			DAcquireS2 = ConvertToDependent4(AcquireS1, AcquireS2, LoseS1, Marry1, 2);
			DLoseS1 = ConvertToDependent4(AcquireS1, AcquireS2, LoseS1, Marry1, 3);
			DMarry1 = ConvertToDependent4(AcquireS1, AcquireS2, LoseS1, Marry1, 4);
		}
		else{
			DLoseS1 = ConvertToDependent2(LoseS1, Marry1, 1);
			DMarry1 = ConvertToDependent2(LoseS1, Marry1, 2);
		}
	}
	if (CurrPartners == 1 && MarriedInd == 1){
		if (RiskGroup == 1){
			DAcquireS1 = ConvertToDependent3(AcquireS1, AcquireS2, Divorce, 1);
			DAcquireS2 = ConvertToDependent3(AcquireS1, AcquireS2, Divorce, 2);
			DDivorce = ConvertToDependent3(AcquireS1, AcquireS2, Divorce, 3);
		}
		else{
			DDivorce = ConvertToDependent1(Divorce);
		}
	}
	if (CurrPartners == 2 && MarriedInd == 0){
		DLoseS1 = ConvertToDependent4(LoseS1, LoseS2, Marry1, Marry2, 1);
		DLoseS2 = ConvertToDependent4(LoseS1, LoseS2, Marry1, Marry2, 2);
		DMarry1 = ConvertToDependent4(LoseS1, LoseS2, Marry1, Marry2, 3);
		DMarry2 = ConvertToDependent4(LoseS1, LoseS2, Marry1, Marry2, 4);
	}
	if (CurrPartners == 2 && MarriedInd == 1){
		DDivorce = ConvertToDependent2(Divorce, LoseS2, 1);
		DLoseS2 = ConvertToDependent2(Divorce, LoseS2, 2);
	}

	// Determine which event occurs
	CumEventProb[7] = 1.0;
	CumEventProb[6] = CumEventProb[7] - DLoseS2;
	CumEventProb[5] = CumEventProb[6] - DLoseS1;
	CumEventProb[4] = CumEventProb[5] - DDivorce;
	CumEventProb[3] = CumEventProb[4] - DMarry2;
	CumEventProb[2] = CumEventProb[3] - DMarry1;
	CumEventProb[1] = CumEventProb[2] - DAcquireS2;
	CumEventProb[0] = CumEventProb[1] - DAcquireS1;

	for (ii = 0; ii<10; ii++){
		if (rnd<CumEventProb[ii]){
			EventType = ii;
			break;
		}
	}

	return EventType;
}

double Indiv::ConvertToDependent1(double rate1)
{
	double prob;
	prob = 1.0 - exp(-rate1);
	return prob;
}

double Indiv::ConvertToDependent2(double rate1, double rate2, int ord)
{
	double prob1, prob2, prob;

	if (TransitionCalc == 0){
		prob1 = (1.0 - exp(-rate1)) * (1.0 - 0.5 * (1.0 - exp(-rate2)));
		prob2 = (1.0 - exp(-rate2)) * (1.0 - 0.5 * (1.0 - exp(-rate1)));
	}
	else{
		prob1 = (1.0 - exp(-rate1 - rate2)) * rate1 / (rate1 + rate2);
		prob2 = (1.0 - exp(-rate1 - rate2)) * rate2 / (rate1 + rate2);
	}

	if (ord == 1){
		prob = prob1;
	}
	if (ord == 2){
		prob = prob2;
	}

	return prob;
}

double Indiv::ConvertToDependent3(double rate1, double rate2, double rate3, int ord)
{
	double prob, totalrate;

	if (TransitionCalc == 0){
		if (ord == 1){
			prob = (1.0 - exp(-rate1)) * (1.0 - 0.5 * (2.0 - exp(-rate2) - exp(-rate3)) +
				(1.0 - exp(-rate2)) * (1.0 - exp(-rate3)) / 3.0);
		}
		if (ord == 2){
			prob = (1.0 - exp(-rate2)) * (1.0 - 0.5 * (2.0 - exp(-rate1) - exp(-rate3)) +
				(1.0 - exp(-rate1)) * (1.0 - exp(-rate3)) / 3.0);
		}
		if (ord == 3){
			prob = (1.0 - exp(-rate3)) * (1.0 - 0.5 * (2.0 - exp(-rate1) - exp(-rate2)) +
				(1.0 - exp(-rate1)) * (1.0 - exp(-rate2)) / 3.0);
		}
	}
	else{
		totalrate = rate1 + rate2 + rate3;
		if (ord == 1){
			prob = (1.0 - exp(-totalrate)) * rate1 / totalrate;
		}
		if (ord == 2){
			prob = (1.0 - exp(-totalrate)) * rate2 / totalrate;
		}
		if (ord == 3){
			prob = (1.0 - exp(-totalrate)) * rate3 / totalrate;
		}
	}

	return prob;
}

double Indiv::ConvertToDependent4(double rate1, double rate2, double rate3, double rate4, int ord)
{
	double prob, totalrate;
	double iprob1, iprob2, iprob3, iprob4;

	if (TransitionCalc == 0){
		iprob1 = 1.0 - exp(-rate1);
		iprob2 = 1.0 - exp(-rate2);
		iprob3 = 1.0 - exp(-rate3);
		iprob4 = 1.0 - exp(-rate4);
		if (ord == 1){
			prob = iprob1 * (1.0 - 0.5 * (iprob2 + iprob3 + iprob4) + (iprob2 * iprob3 +
				iprob2 * iprob4 + iprob3 * iprob4) / 3.0 - 0.25 * iprob2 * iprob3 * iprob4);
		}
		if (ord == 2){
			prob = iprob2 * (1.0 - 0.5 * (iprob1 + iprob3 + iprob4) + (iprob1 * iprob3 +
				iprob1 * iprob4 + iprob3 * iprob4) / 3.0 - 0.25 * iprob1 * iprob3 * iprob4);
		}
		if (ord == 3){
			prob = iprob3 * (1.0 - 0.5 * (iprob1 + iprob2 + iprob4) + (iprob1 * iprob2 +
				iprob1 * iprob4 + iprob2 * iprob4) / 3.0 - 0.25 * iprob1 * iprob2 * iprob4);
		}
		if (ord == 4){
			prob = iprob4 * (1.0 - 0.5 * (iprob1 + iprob2 + iprob3) + (iprob1 * iprob2 +
				iprob1 * iprob3 + iprob2 * iprob3) / 3.0 - 0.25 * iprob1 * iprob2 * iprob3);
		}
	}
	else{
		totalrate = rate1 + rate2 + rate3 + rate4;
		if (ord == 1){
			prob = (1.0 - exp(-totalrate)) * rate1 / totalrate;
		}
		if (ord == 2){
			prob = (1.0 - exp(-totalrate)) * rate2 / totalrate;
		}
		if (ord == 3){
			prob = (1.0 - exp(-totalrate)) * rate3 / totalrate;
		}
		if (ord == 4){
			prob = (1.0 - exp(-totalrate)) * rate4 / totalrate;
		}
	}

	return prob;
}

void Indiv::SimulateSexActs(int ID)
{
	// Note that this function only gets called if the individual is alive and is male.
	// Also note that in this function we're assuming that the specified frequencies of sex in
	// the SexAssumps sheet are monthly.
	// From Raikov's theorem we know that if the sum of 2 independent RVs is Poisson-distributed
	// then so is each of the RVs. So we can simulate the # unprotected and the # protected sex
	// acts as separate Poisson processes.

	double r[10]; // random numbers
	double ExpActs, ExpCondomUse;
	int ii, IndAge, PartnerAge, x;
	double xlam;

	// Note that I've arbitrarily fixed the seed for now.
	//int seed = 1272 + ID * 8 + (2025 - CurrYear) * 27 + BehavCycleCount * 31 + STDcycleCount * 95; 
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for (ii = 0; ii<10; ii++){
		r[ii] = rg.Random();
	}

	IndAge = AgeGroup - 2;

	// Sex acts with primary partner
	if (CurrPartners>0){
		PartnerAge = Register[IDprimary - 1].AgeGroup - 2;
		if (MarriedInd == 1){
			ExpActs = (GenderEquality * FreqSexLT[PartnerAge][1] + (1.0 - GenderEquality) *
				FreqSexLT[IndAge][0]) * 12.0 / CycleD;
			//if(AllowPartnerRateAdj==1){
			//	ExpActs = (GenderEquality * FreqSexLT[PartnerAge][1] * Register[IDprimary-1].PartnerRateAdj + 
			//		(1.0 - GenderEquality) * FreqSexLT[IndAge][0] * PartnerRateAdj) * 12.0/CycleD;
			//}
			ExpCondomUse = GenderEquality * CondomUseLT[PartnerAge][1] + (1.0 - GenderEquality) *
				CondomUseLT[IndAge][0];
		}
		else{
			ExpActs = (GenderEquality * FreqSexST[PartnerAge][1] + (1.0 - GenderEquality) *
				FreqSexST[IndAge][0]) * 12.0 / CycleD;
			//if(AllowPartnerRateAdj==1){
			//	ExpActs = (GenderEquality * FreqSexST[PartnerAge][1] * Register[IDprimary-1].PartnerRateAdj + 
			//		(1.0 - GenderEquality) * FreqSexST[IndAge][0] * PartnerRateAdj) * 12.0/CycleD;
			//}
			ExpCondomUse = GenderEquality * CondomUseST[PartnerAge][1] + (1.0 - GenderEquality) *
				CondomUseST[IndAge][0];
		}

		// Simulate # unprotected sex acts
		xlam = ExpActs * (1.0 - ExpCondomUse);
		x = RandomSexGenerator(r[0], xlam); //this function returns the number of sex acts 
		UVIprimary = x;
		totUVIprimary += x;
		if (Register[IDprimary - 1].IDprimary == ID){
			Register[IDprimary - 1].UVIprimary = x;
			Register[IDprimary - 1].totUVIprimary += x;
		}
		else{
			Register[IDprimary - 1].UVI2ndary = x;
			Register[IDprimary - 1].totUVI2ndary += x;
		}

		// Simulate # protected sex acts
		xlam = ExpActs * ExpCondomUse;
		x = RandomSexGenerator(r[1], xlam);
		PVIprimary = x;
		totPVIprimary += x;
		if (Register[IDprimary - 1].IDprimary == ID){
			Register[IDprimary - 1].PVIprimary = x;
			Register[IDprimary - 1].totPVIprimary += x;
		}
		else{
			Register[IDprimary - 1].PVI2ndary = x;
			Register[IDprimary - 1].totPVI2ndary += x;
		}
	}
	else{
		UVIprimary = 0;
		PVIprimary = 0;
	}

	// Sex acts with 2ndary partner
	if (CurrPartners == 2){
		PartnerAge = Register[ID2ndary - 1].AgeGroup - 2;
		ExpActs = (GenderEquality * FreqSexST[PartnerAge][1] + (1.0 - GenderEquality) *
			FreqSexST[IndAge][0]) * 12.0 / CycleD;
		//if(AllowPartnerRateAdj==1){
		//	ExpActs = (GenderEquality * FreqSexST[PartnerAge][1] * Register[ID2ndary-1].PartnerRateAdj + 
		//		(1.0 - GenderEquality) * FreqSexST[IndAge][0] * PartnerRateAdj) * 12.0/CycleD;
		//}
		ExpCondomUse = GenderEquality * CondomUseST[PartnerAge][1] + (1.0 - GenderEquality) *
			CondomUseST[IndAge][0];

		// Simulate # unprotected sex acts
		xlam = ExpActs * (1.0 - ExpCondomUse);
		x = RandomSexGenerator(r[2], xlam);
		UVI2ndary = x;
		totUVI2ndary += x;
		if (Register[ID2ndary - 1].IDprimary == ID){
			Register[ID2ndary - 1].UVIprimary = x;
			Register[ID2ndary - 1].totUVIprimary += x;
		}
		else{
			Register[ID2ndary - 1].UVI2ndary = x;
			Register[ID2ndary - 1].totUVI2ndary += x;
		}

		// Simulate # protected sex acts
		xlam = ExpActs * ExpCondomUse;
		x = RandomSexGenerator(r[3], xlam);
		PVI2ndary = x;
		totPVI2ndary += x;
		if (Register[ID2ndary - 1].IDprimary == ID){
			Register[ID2ndary - 1].PVIprimary = x;
			Register[ID2ndary - 1].totPVIprimary += x;
		}
		else{
			Register[ID2ndary - 1].PVI2ndary = x;
			Register[ID2ndary - 1].totPVI2ndary += x;
		}
	}
	else{
		UVI2ndary = 0;
		PVI2ndary = 0;
	}

	// Sex acts with sex workers
	if (RiskGroup == 1 && TotCurrFSW>0){
		ExpActs = FSWcontactConstant * AgeEffectFSWcontact[IndAge] / CycleD;
		if (CurrPartners == 0){
			ExpActs *= PartnerEffectFSWcontact[0];
		}
		if (CurrPartners == 1 && MarriedInd == 0){
			ExpActs *= PartnerEffectFSWcontact[1];
		}
		if (CurrPartners == 1 && MarriedInd == 1){
			ExpActs *= PartnerEffectFSWcontact[2];
		}
		if (CurrPartners == 2 && MarriedInd == 0){
			ExpActs *= PartnerEffectFSWcontact[3];
		}
		if (CurrPartners == 2 && MarriedInd == 1){
			ExpActs *= PartnerEffectFSWcontact[4];
		}
		ExpCondomUse = CondomUseFSW;

		// Simulate # unprotected sex acts
		xlam = ExpActs * (1.0 - ExpCondomUse);
		x = RandomSexGenerator(r[4], xlam);
		UVICSW = x;
		totUVICSW += x;

		// Simulate # protected sex acts
		xlam = ExpActs * ExpCondomUse;
		x = RandomSexGenerator(r[5], xlam);
		PVICSW = x;
		totPVICSW += x;

		// Determine which sex worker the man has sex with.
		// (For simplicity we're assuming no man has contacts
		// with >1 CSW in a single STD cycle.)
		IDofCSW = static_cast<int>(r[6] * TotCurrFSW);
	}
	else{
		UVICSW = 0;
		PVICSW = 0;
	}
}

int Indiv::RandomSexGenerator(double p, double lambda)
{
	// Due to rounding there is a very small risk that the algorithm for returning
	// the Poisson-distributed RV might get into an infinite loop. To avoid this, I've
	// limited the # sex acts that can occur to at most 720 per year (which is a lot).

	int ii, SexActs, MaxSex;
	double CumProb, NextProb;

	MaxSex = 720 / CycleD;
	NextProb = exp(-lambda);
	if (p<NextProb){
		SexActs = 0;
	}
	else{
		CumProb = NextProb;
		SexActs = MaxSex;
		for (ii = 1; ii<MaxSex; ii++){
			NextProb *= lambda / ii;
			CumProb += NextProb;
			if (p<CumProb){
				SexActs = ii;
				break;
			}
		}
	}

	return SexActs;
}

void Indiv::GetNewHIVstate(int ID, double p, double p2, double p3, double p4, double RandAcceptTxV, double Rloss, double EfficacyD1Rand, double EfficacyD2Rand)
{
	double Prob1;
	int iy, xx, zz;
	double adjustedTxVEfficacyD1[13];
	double adjustedTxVEfficacyD2[13];
	if (HIVstage==0) {zz=0;}
		else if(HIVstage==5) {zz=1;} 
		else  {zz=2;}	
	if (HIVstage == 0)
	{// Note that we refer to opposite sex when getting transm prob
		if (SexInd == 0)	{Prob1 = HIVtransitionF.GetTransmProb(ID);}
		else
		{Prob1 = HIVtransitionM.GetTransmProb(ID);}
		if (p < Prob1)
		{HIVstageE = 1;
			DateInfect = CurrYear + 0.5 + (STDcycleCount - 1.0) / CycleD +
						 (BehavCycleCount - 1.0) / CycleS + (p / Prob1) / CycleD;
			// if (FixedUncertainty == 0){ NewHIV.out[CurrSim - 1][CurrYear - StartYear] += 1; }
			iy = AgeGroup;
			// NewHIV.out[iy][CurrYear - StartYear] += 1;
		}else{HIVstageE = 0;}
	}else
	{
		if (SexInd == 0){HIVtransitionM.GetNewStage(ID, p);}
			else{HIVtransitionF.GetNewStage(ID, p);}}
	if ((HIVstage == 1 || HIVstage == 2 || HIVstage == 6) && HIVstageE == 5){ARTstage = 0;}
	else if ((HIVstage == 3 || HIVstage == 4) && HIVstageE == 5){ARTstage = 1;}
	if (HIVstage == 5 && HIVstageE == 5){ARTweeks += 1;}
	if (HIVstage == 5 && HIVstageE == 6){ARTweeks = 0;}

	if (HIVstage != 5 && HIVstage != 6 && HIVstageE == 5 && timetoScreen > 0 && ScreenResult == 0 &&
		((WHOScreening == 0 && PerfectSchedule == 0) || CurrYear < ImplementYR))
	{
		timetoScreen *= 0.2;
	}
	for (xx = 0; xx < 13; xx++)
	{
		if (HIVstage == 0 && HIVstageE == 1 && HPVstage[xx] == 4 && WeibullCIN3[xx] > 0)
		{
			WeibullCIN3[xx] *= HPVTransF[xx].CIN3HIV;
		}
	}
	//  deliver TxV to those who are starting ART
	// if(CatchUpVaccHIV==1 && HIVstage == 0 && HIVstageE == 1 && CurrYear<=(ImplementYR+16) && CurrYear>=ImplementYR &&
/*	if (CatchUpVaccHIV == 1 && HIVstage != 5 && HIVstage != 6 && HIVstageE == 5 && CurrYear >= ImplementYR && // && CurrYear<=(ImplementYR+16)
		(AgeExact >= CatchUpAgeMIN && AgeExact < CatchUpAgeMAX) && AliveInd == 1 && SexInd == 1 && GotVaccOffer == 0)
	{
		GotVaccOffer = 1;
		if (p2 < CatchUpCoverage)
		{GotVacc = 1;
			RSApop.NewVACC[18 * SexInd + AgeGroup][CurrYear - StartYear] += 1;
			if (WHOvacc == 1)
			{
				VaccinationStatus[0] = 1;
				VaccinationStatus[1] = 1;
				VaccinationStatus[2] = 1;
				VaccinationStatus[3] = 1;
				VaccinationStatus[6] = 1;
				VaccinationStatus[8] = 1;
				VaccinationStatus[10] = 1;
			}
			else
			{
				VaccinationStatus[0] = 1;
				VaccinationStatus[1] = 1;
				if (p3 < 0.5)
				{
					VaccinationStatus[2] = 1;
					VaccinationStatus[3] = 1;
					VaccinationStatus[6] = 1;
				}
			}
			if (VaccineWane == 1)
			{
				TimeVacc = 0;
				ExpVacc = 48 * 20 + (-48 * VaccDur * log(p4));
				if (ExpVacc == 0)
				{
					ExpVacc = 1;
				}
			}
		}
	}
*/
	if (AdministerMassTxVtoART == 1 && HIVstage != 5 && HIVstage != 6 && HIVstageE == 5 && CurrYear >= 2030 && // && CurrYear<=(ImplementYR+16)
		(AgeExact >= MassTxVtoARTAgeMIN && AgeExact < MassTxVtoARTAgeMAX) && AliveInd == 1 && SexInd == 1)
	{	
		{
			if ((RandAcceptTxV < 0.9) && ((CatchUpVaccHIV == 1 && GotTxV == 0) || (CatchUpVaccHIV == 0)))
			{ 
				if (CatchUpVaccHIV ==1 ){
					RSApop.NewVACC[18 * SexInd + AgeGroup][CurrYear - StartYear] += 1;
					GotVacc = 1;
					}
				RSApop.NewTxV[zz*18 + AgeGroup][CurrYear - StartYear] += 1;
				GotTxV = 1;
				for (int yy = 0; yy < 13; yy++)
				{
					if (HPVstage[yy] > 2 && HPVstage[yy] <= 4)
					{
						adjustedTxVEfficacyD1[yy] = TxVEfficacyD1CIN[yy];
						adjustedTxVEfficacyD2[yy] = TxVEfficacyD2CIN[yy];
					}
					else
					{
						adjustedTxVEfficacyD1[yy] = TxVEfficacyD1[yy];
						adjustedTxVEfficacyD2[yy] = TxVEfficacyD2[yy];
					}
					if (!(HIVstage == 0 || HIVstage == 5)) // || Register[ic].HIVstage == 5 && artweeks 
					{
						adjustedTxVEfficacyD1[yy] *= ReductionFactorHIV;
						adjustedTxVEfficacyD2[yy] *= ReductionFactorHIV;
					}
					if (HIVstage == 5)
					{
						adjustedTxVEfficacyD1[yy] *= ReductionFactorHIVART;
						adjustedTxVEfficacyD2[yy] *= ReductionFactorHIVART;
					}
					if (EfficacyD1Rand < adjustedTxVEfficacyD1[yy])
					{
						if ((HPVstage[yy] >= 1 && HPVstage[yy] <= 4) || (HPVstage[yy] == 6))
						{ 
							HPVstage[yy] = 0;
							HPVstageE[yy] = 0;
							if (CatchUpVaccHIV ==1 ){
								VaccinationStatus[yy] = 1;
								}
						}
					}
				} 
				if (Rloss < 0.7)
				{RSApop.NewTxV[zz*18 + AgeGroup][CurrYear - StartYear] += 1;
					GotTxV = 2;
					if (CatchUpVaccHIV ==1 ){
						RSApop.NewVACC[18 * SexInd + AgeGroup][CurrYear - StartYear] += 1;
						GotVacc = 1;
						}
					for (int yy = 0; yy < 13; yy++)
					{
						if (EfficacyD2Rand < adjustedTxVEfficacyD2[yy])
						{
							if ((HPVstage[yy] >= 1 && HPVstage[yy] <= 4) || (HPVstage[yy] == 6))
							{
								HPVstage[yy] = 0;
								HPVstageE[yy] = 0;
								if (CatchUpVaccHIV ==1 ){
									VaccinationStatus[yy] = 1;
									}
							}
						}
					}
				}
			}
		}
	}
}

void Indiv::GetNewHSVstate(int ID, double p)
{
	double Prob1, Symptomatic;

	if (HSVstage == 0){
		// Note that we refer to opposite sex when getting transm prob
		if (SexInd == 0){
			Prob1 = HSVtransitionF.GetTransmProb(ID);
		}
		else{
			Prob1 = HSVtransitionM.GetTransmProb(ID);
		}
		if (p<Prob1){
			if (SexInd == 0){ Symptomatic = HSVtransitionM.SymptomaticPropn; }
			else{ Symptomatic = HSVtransitionF.SymptomaticPropn; }
			if (p / Prob1 < Symptomatic){ HSVstageE = 1; }
			else{ HSVstageE = 4; }
			if (FixedUncertainty == 1){ NewHSV.out[CurrSim - 1][CurrYear - StartYear] += 1; }
		}
		else{
			HSVstageE = 0;
		}
	}
	else{
		if (SexInd == 0){
			HSVtransitionM.GetNewStage(ID, p);
		}
		else{
			HSVtransitionF.GetNewStage(ID, p);
		}
	}
}

void Indiv::GetNewTPstate(int ID, double p)
{
	double Prob1;

	if (TPstage == 0){
		// Note that we refer to opposite sex when getting transm prob
		if (SexInd == 0){
			Prob1 = TPtransitionF.GetTransmProb(ID);
		}
		else{
			Prob1 = TPtransitionM.GetTransmProb(ID);
		}
		if (p<Prob1){
			TPstageE = 1;
			if (CurrYear > 2009) { cout << CurrYear << " success" << endl; }
			if (FixedUncertainty == 1){ NewTP.out[CurrSim - 1][CurrYear - StartYear] += 1; }
		}
		else{ TPstageE = 0; }
	}
	else{
		if (SexInd == 0){
			TPtransitionM.GetNewStage(ID, p);
		}
		else{
			TPtransitionF.GetNewStage(ID, p);
		}
	}
}

void Indiv::GetNewHDstate(int ID, double p)
{
	double Prob1, Symptomatic;

	if (HDstage == 0){
		// Note that we refer to opposite sex when getting transm prob
		if (SexInd == 0){
			Prob1 = HDtransitionF.GetTransmProb(ID, 4);
		}
		else{
			Prob1 = HDtransitionM.GetTransmProb(ID, 4);
		}
		if (p<Prob1){
			if (SexInd == 0){ Symptomatic = HDtransitionM.SymptomaticPropn; }
			else{ Symptomatic = HDtransitionF.SymptomaticPropn; }
			if (p / Prob1 < Symptomatic){ HDstageE = 1; }
			else{ HDstageE = 2; }
		}
		else{
			HDstageE = 0;
		}
	}
	else{
		if (SexInd == 0){
			HDtransitionM.GetNewStage(ID, p, 4);
		}
		else{
			HDtransitionF.GetNewStage(ID, p, 4);
		}
	}
}

void Indiv::GetNewNGstate(int ID, double p)
{
	double Prob1, Symptomatic;

	if (NGstage == 0){
		// Note that we refer to opposite sex when getting transm prob
		if (SexInd == 0){
			Prob1 = NGtransitionF.GetTransmProb(ID, 1);
		}
		else{
			Prob1 = NGtransitionM.GetTransmProb(ID, 1);
		}
		if (p<Prob1){
			if (SexInd == 0){ Symptomatic = NGtransitionM.SymptomaticPropn; }
			else{ Symptomatic = NGtransitionF.SymptomaticPropn; }
			if (p / Prob1 < Symptomatic){ NGstageE = 1; }
			else{ NGstageE = 2; }
			if (FixedUncertainty == 1){ NewNG.out[CurrSim - 1][CurrYear - StartYear] += 1; }
		}
		else{
			NGstageE = 0;
		}
	}
	else{
		if (SexInd == 0){
			NGtransitionM.GetNewStage(ID, p, 1);
		}
		else{
			NGtransitionF.GetNewStage(ID, p, 1);
		}
	}
}

void Indiv::GetNewCTstate(int ID, double p)
{
	double Prob1, Symptomatic;

	if (CTstage == 0){
		// Note that we refer to opposite sex when getting transm prob
		if (SexInd == 0){
			Prob1 = CTtransitionF.GetTransmProb(ID, 2);
		}
		else{
			Prob1 = CTtransitionM.GetTransmProb(ID, 2);
		}
		if (p<Prob1){
			if (SexInd == 0){ Symptomatic = CTtransitionM.SymptomaticPropn; }
			else{ Symptomatic = CTtransitionF.SymptomaticPropn; }
			if (p / Prob1 < Symptomatic){ CTstageE = 1; }
			else{ CTstageE = 2; }
			if (FixedUncertainty == 1){ NewCT.out[CurrSim - 1][CurrYear - StartYear] += 1; }
		}
		else{
			CTstageE = 0;
		}
	}
	else{
		if (SexInd == 0){
			CTtransitionM.GetNewStage(ID, p, 2);
		}
		else{
			CTtransitionF.GetNewStage(ID, p, 2);
		}
	}
}

void Indiv::GetNewTVstate(int ID, double p)
{
	double Prob1, Symptomatic;

	if (TVstage == 0){
		// Note that we refer to opposite sex when getting transm prob
		if (SexInd == 0){
			Prob1 = TVtransitionF.GetTransmProb(ID, 3);
		}
		else{
			Prob1 = TVtransitionM.GetTransmProb(ID, 3);
		}
		if (p<Prob1){
			if (SexInd == 0){ Symptomatic = TVtransitionM.SymptomaticPropn; }
			else{ Symptomatic = TVtransitionF.SymptomaticPropn; }
			if (p / Prob1 < Symptomatic){ TVstageE = 1; }
			else{ TVstageE = 2; }
			if (FixedUncertainty == 1){ NewTV.out[CurrSim - 1][CurrYear - StartYear] += 1; }
		}
		else{
			TVstageE = 0;
		}
	}
	else{
		if (SexInd == 0){
			TVtransitionM.GetNewStage(ID, p, 3);
		}
		else{
			TVtransitionF.GetNewStage(ID, p, 3);
		}
	}
}

Pop::Pop(int xxx){}

void Pop::AssignAgeSex()
{
	int ic, ir, ia, male, female;
	male = 0;
	female = 0;
	double TotalPop;
	double AgeExact;
	vector<Indiv>::const_iterator ic2;

	TotalPop = 0.0;
	for (ia = 0; ia<91; ia++){
		TotalPop += StartPop[ia][0] + StartPop[ia][1];
	}
	StartPop[0][0] = StartPop[0][0] / TotalPop;
	for (ia = 1; ia<91; ia++){
		StartPop[ia][0] = StartPop[ia - 1][0] + StartPop[ia][0] / TotalPop;
	}
	StartPop[0][1] = StartPop[90][0] + StartPop[0][1] / TotalPop;
	for (ia = 1; ia<91; ia++){
		StartPop[ia][1] = StartPop[ia - 1][1] + StartPop[ia][1] / TotalPop;
	}

	//int seed = 5481; // Note that I've arbitrarily fixed the seed for now.
	//if(FixedPeriod==0){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for (ir = 0; ir<InitPop; ir++){
		r[ir] = rg.Random();
	}

	int tpp = Register.size();

	// Simulate sex and date of birth for each individual
	for (ic = 0; ic<tpp; ic++){
		Register[ic].AliveInd = 1;
		if (r[ic] <= StartPop[90][0]){
			Register[ic].SexInd = 0;
			male += 1;
			if (r[ic] <= StartPop[0][0]){
				Register[ic].DOB = 0.5 + StartYear - r[ic] / StartPop[0][0];
			}
			else for (ir = 1; ir<91; ir++){
				if (r[ic]>StartPop[ir - 1][0] && r[ic] <= StartPop[ir][0]){
					Register[ic].DOB = 0.5 + StartYear - ir - (r[ic] -
						StartPop[ir - 1][0]) / (StartPop[ir][0] - StartPop[ir - 1][0]);
					break;
				}
			}
		}
		else{
			Register[ic].SexInd = 1;
			female += 1;
			if (r[ic] <= StartPop[0][1]){
				Register[ic].DOB = 0.5 + StartYear - (r[ic] - StartPop[90][0]) /
					(StartPop[0][1] - StartPop[90][0]);
			}
			else for (ir = 1; ir<91; ir++){
				if (r[ic]>StartPop[ir - 1][1] && r[ic] <= StartPop[ir][1]){
					Register[ic].DOB = 0.5 + StartYear - ir - (r[ic] -
						StartPop[ir - 1][1]) / (StartPop[ir][1] - StartPop[ir - 1][1]);
					break;
				}
			}
		}
	}

	// Calculate age group for each individual
	for (ic = 0; ic<InitPop; ic++){
		AgeExact = StartYear + 0.5 - Register[ic].DOB;
		if (AgeExact < 90.0){
			Register[ic].AgeGroup = static_cast<int> (AgeExact / 5);
		}
		else{
			Register[ic].AgeGroup = 17;
		}
		if(AgeExact < 50.0){
			Register[ic].Age50 = 0;
		}
		else{
			Register[ic].Age50 = 1;
		}
	}
}

void Pop::GetAllPartnerRates()
{
	// The same as the GetAllPartnerRates function in TSHISA model

	int ia, ig;
	double lambda, alpha;
	double GammaScaling[2];

	// Calculate the age adjustment factors and the high risk partner acquisition rates
	for (ig = 0; ig<2; ig++){
		lambda = (GammaMeanST[ig] - 10.0) / pow(GammaStdDevST[ig], 2.0);
		alpha = (GammaMeanST[ig] - 10.0) * lambda;
		for (ia = 0; ia<16; ia++){
			AgeEffectPartners[ia][ig] = pow(lambda, alpha) * pow(ia*5.0 + 2.5, alpha - 1.0) *
				exp(-lambda * (ia*5.0 + 2.5));
		}
		GammaScaling[ig] = BasePartnerAcqH[ig] / AgeEffectPartners[1][ig];
		for (ia = 0; ia<16; ia++){
			AgeEffectPartners[ia][ig] *= GammaScaling[ig];
		}
		PartnershipFormation[0][ig] = 1.0;
	}
}

void Pop::AssignBehav()
{
	int ic, ia, ib, ig, ii, ij;
	int ExactAge, found; // , index1;
	double MarriageProb;
	int UnpartneredCount; // # married individuals not yet assigned a partner
	double TotalSelectionProb, partnerrisk, partnerID, temp1;
	int STpartners, ind;
	double x,  a, b, p, q; //y,
	vector<int> indices(InitPop);

	// In the arrays defined below, the I- prefix indicates individual ages, to be
	// consistent with the TSHISA model. Age starts at 10 (not 0), since 10 is the
	// assumed youngest age of sexual activity.
	double IMarriedPropn[81][2], IVirginPropnH[81][2], IVirginPropnL[81][2];
	//double IVirginsH[81][2]; // IVirginsL[81][2];
	//double ISexuallyExpH[81][2] , ISexuallyExpL[81][2];
	//double IMarriedH[81][2], IMarriedL[81][2];
	//double ISexuallyExpUnmarriedH[81][2], ISexuallyExpUnmarriedL[81][2];
	double SpouseMortM, SpouseMortF;

	double cHD, cLD; // Rates of partnership formation x ave duration in high & low risk
	double GUnmarriedH0[16][2], GUnmarriedH1[16][2], GUnmarriedH2[16][2];
	double GUnmarriedL0[16][2], GUnmarriedL1[16][2];
	double GMarriedH1[16][2], GMarriedH2[16][2];

	double alpha, lambda; // Parameters for the gamma dbn used to determine relative
	// frequencies of FSW contact at different ages
	double FSWcontactBase[16], EligibleForCSW[16]; //FSWcontactDemand[16], 
	double BaseFSWdemand, MalePop15to49;

	// First calculate IMarriedPropn and IVirginPropn
	IMarriedPropn[0][0] = 0.0;
	IMarriedPropn[0][1] = 0.0;
	IVirginPropnH[0][0] = 1.0;
	IVirginPropnH[0][1] = 1.0;
	IVirginPropnL[0][0] = 1.0;
	IVirginPropnL[0][1] = 1.0;
	for (ia = 1; ia<81; ia++){
		SpouseMortM = 0;
		SpouseMortF = 0;
		for (ib = 0; ib<16; ib++){
			SpouseMortM += AgePrefM[(ia - 1) / 5][ib] * NonAIDSmortF[ib][0];
			SpouseMortF += AgePrefF[(ia - 1) / 5][ib] * NonAIDSmortM[ib][0];
		}
		IMarriedPropn[ia][0] = (1.0 - IMarriedPropn[ia - 1][0]) * (1.0 - exp(
			-MarriageIncidence[(ia - 1) / 5][0])) + IMarriedPropn[ia - 1][0] * exp(
			-LTseparation[(ia - 1) / 5][0]) * (1.0 - SpouseMortM);
		IMarriedPropn[ia][1] = (1.0 - IMarriedPropn[ia - 1][1]) * (1.0 - exp(
			-MarriageIncidence[(ia - 1) / 5][1])) + IMarriedPropn[ia - 1][1] * exp(
			-LTseparation[(ia - 1) / 5][1]) * (1.0 - SpouseMortF);
		for (ig = 0; ig<2; ig++){
			if (ia<20){
				IVirginPropnH[ia][ig] = IVirginPropnH[ia - 1][ig] *
					exp(-SexualDebut[(ia - 1) / 5][ig]);
				IVirginPropnL[ia][ig] = IVirginPropnL[ia - 1][ig] *
					exp(-SexualDebut[(ia - 1) / 5][ig] * DebutAdjLow[ig]);
			}
			else{
				IVirginPropnH[ia][ig] = 0;
				IVirginPropnL[ia][ig] = 0;
			}
		}
	}

	// Now assign risk group to each individual

	//int seed = 2072; // Note that I've arbitrarily fixed the seed for now.
	//if(FixedPeriod==0){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		r[ic] = rg.Random();
	}

	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].SexInd == 0){
			if (r[ic] < HighPropnM){
				Register[ic].RiskGroup = 1;
			}
			else{
				Register[ic].RiskGroup = 2;
			}
		}
		else{
			if (r[ic] < HighPropnF){
				Register[ic].RiskGroup = 1;
			}
			else{
				Register[ic].RiskGroup = 2;
			}
		}
	}

	// Next assign virgin status to each individual

	//int seed2 = 9290; // Note that I've arbitrarily fixed the seed for now.
	//if(FixedPeriod==0){
	//	seed2 += CurrSim;}
	//CRandomMersenne rg2(seed2);
	for (ic = 0; ic<tpp; ic++){
		r[ic] = rg.Random();
	}

	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AgeGroup < 2){
			Register[ic].VirginInd = 1;
		}
		else if (Register[ic].AgeGroup >= 6){
			Register[ic].VirginInd = 0;
		}
		else{
			ExactAge = static_cast<int> (StartYear + 0.5 - Register[ic].DOB - 10);
			ig = Register[ic].SexInd;
			if (Register[ic].RiskGroup == 1){
				if (r[ic] < IVirginPropnH[ExactAge][ig]){
					Register[ic].VirginInd = 1;
				}
				else{
					Register[ic].VirginInd = 0;
				}
			}
			else{
				if (r[ic] < IVirginPropnL[ExactAge][ig]){
					Register[ic].VirginInd = 1;
				}
				else{
					Register[ic].VirginInd = 0;
				}
			}
		}
	}

	// Next assign marital status to each individual

	//int seed3 = 59; // Note that I've arbitrarily fixed the seed for now.
	//if(FixedPeriod==0){
	//	seed3 += CurrSim;}
	//CRandomMersenne rg3(seed3);
	for (ic = 0; ic<tpp; ic++){
		r[ic] = rg.Random();
	}

	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].VirginInd == 1){
			Register[ic].MarriedInd = 0;
		}
		else{
			ExactAge = static_cast<int> (StartYear + 0.5 - Register[ic].DOB - 10);
			ig = Register[ic].SexInd;
			MarriageProb = IMarriedPropn[ExactAge][ig];
			if (ExactAge<20){
				if (Register[ic].RiskGroup == 1){
					MarriageProb = MarriageProb / (1.0 - IVirginPropnH[ExactAge][ig]);
				}
				else{
					MarriageProb = MarriageProb / (1.0 - IVirginPropnL[ExactAge][ig]);
				}
			}
			if (r[ic] < MarriageProb){
				Register[ic].MarriedInd = 1;
			}
			else{
				Register[ic].MarriedInd = 0;
			}
		}
	}

	// Next calculate sexual mixing in married individuals (copied from code in TSHISA)

	DesiredMarriagesM[0][0] = ((1.0 - AssortativeM) + AssortativeM * HighPropnF) * HighPropnM;
	DesiredMarriagesM[0][1] = AssortativeM * (1.0 - HighPropnF) * HighPropnM;
	DesiredMarriagesM[1][1] = ((1.0 - AssortativeM) + AssortativeM * (1.0 - HighPropnF)) *
		(1.0 - HighPropnM);
	DesiredMarriagesM[1][0] = AssortativeM * HighPropnF * (1 - HighPropnM);
	DesiredMarriagesF[0][0] = ((1.0 - AssortativeF) + AssortativeF * HighPropnM) * HighPropnF;
	DesiredMarriagesF[0][1] = AssortativeF * (1.0 - HighPropnM) * HighPropnF;
	DesiredMarriagesF[1][1] = ((1.0 - AssortativeF) + AssortativeF * (1.0 - HighPropnM)) *
		(1.0 - HighPropnF);
	DesiredMarriagesF[1][0] = AssortativeF * HighPropnM * (1 - HighPropnF);

	for (ii = 0; ii<2; ii++){
		for (ij = 0; ij<2; ij++){
			if (AllowBalancing == 1){
				AdjLTrateM[ii][ij] = (GenderEquality * DesiredMarriagesF[ij][ii] +
					(1.0 - GenderEquality) * DesiredMarriagesM[ii][ij]);
				AdjLTrateF[ii][ij] = (GenderEquality * DesiredMarriagesF[ii][ij] +
					(1.0 - GenderEquality) * DesiredMarriagesM[ij][ii]);
			}
			else{
				AdjLTrateM[ii][ij] = DesiredMarriagesM[ii][ij];
				AdjLTrateF[ii][ij] = DesiredMarriagesF[ii][ij];
			}
		}
		ActualPropnLTH[ii][0] = AdjLTrateM[ii][0] / (AdjLTrateM[ii][0] + AdjLTrateM[ii][1]);
		ActualPropnLTH[ii][1] = AdjLTrateF[ii][0] / (AdjLTrateF[ii][0] + AdjLTrateF[ii][1]);
	}

	// Assign partner IDs to each married individual

	//int seed4 = 480; // Note that I've arbitrarily fixed the seeds for now.
	//if(FixedPeriod==0){
	//	seed4 += CurrSim;}
	//CRandomMersenne rg4(seed4);
	for (ic = 0; ic<tpp; ic++){
		r[ic] = rg.Random();
	}
	//int seed5 = 2771; 
	//if(FixedPeriod==0){
	//	seed5 += CurrSim;}
	//CRandomMersenne rg5(seed5);
	for (ic = 0; ic<tpp; ic++){
		rprisk[ic] = rg.Random();
	}
	//int seed6 = 11052; 
	//if(FixedPeriod==0){
	//	seed6 += CurrSim;}
	//CRandomMersenne rg6(seed6);
	for (ic = 0; ic<tpp; ic++){
		rpID[ic] = rg.Random();
	}

	clock_t start2, finish2, finish3;
	double elapsed_time2, elapsed_time3;
	elapsed_time2 = 0.0;
	elapsed_time3 = 0.0;

	MaxNewPartnerInd = 0;
	//for(ii=0; ii<InitPop; ii++){
	//	indices[ii] = ii;}
	ii = 0;
	while (MaxNewPartnerInd == 0 && ii<tpp){
		start2 = clock();
		/*index1 = indices[r[ii] * indices.size()];
		ic = indices[index1];
		if(Register[ic].MarriedInd==1 && Register[ic].IDprimary==0){
		partnerrisk = rprisk[ii];
		partnerID = rpID[ii];
		ChooseLTpartner(ic+1, partnerrisk, partnerID);
		}
		indices[index1] = indices[indices.size()-1];
		indices.pop_back();
		ii += 1;*/
		UnpartneredCount = 0;
		for (ic = 0; ic<tpp; ic++){
			if (Register[ic].MarriedInd == 1 && Register[ic].IDprimary == 0){
				UnpartneredCount += 1;
			}
		}
		TotalSelectionProb = 0.0;
		for (ic = 0; ic<tpp; ic++){
			if (Register[ic].MarriedInd == 1 && Register[ic].IDprimary == 0){
				TotalSelectionProb += 1.0 / UnpartneredCount;
				if (r[ii] <= TotalSelectionProb && r[ii] >(TotalSelectionProb -
					1.0 / UnpartneredCount)){
					finish2 = clock();
					partnerrisk = rprisk[ii];
					partnerID = rpID[ii];
					ChooseLTpartner(ic + 1, partnerrisk, partnerID);
				}
			}
		}
		finish3 = clock();
		elapsed_time2 += (finish2 - start2);
		elapsed_time3 += (finish3 - finish2);
		ii += 1;
	}
	//cout<<"Time to select married individuals with no spouse ID: "<<elapsed_time2<<endl;
	//cout<<"Time to select spouse ID: "<<elapsed_time3<<endl;
	// For married people who haven't been assigned a partner ID, change status to unmarried.
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].MarriedInd == 1 && Register[ic].IDprimary == 0){
			Register[ic].MarriedInd = 0;
		}
	}

	// Assign number of short-term partners to each individual
	for (ig = 0; ig<2; ig++){ // Modified from GetStartProfile function in TSHISA
		for (ia = 0; ia<16; ia++){
			cHD = PartnershipFormation[0][ig] * AgeEffectPartners[ia][ig] * MeanDurSTrel[0][0];
			cLD = PartnershipFormation[1][ig] * AgeEffectPartners[ia][ig] * MeanDurSTrel[0][0];
			GUnmarriedH0[ia][ig] = 1.0 / (1.0 + cHD + 0.5 * PartnerEffectNew[0][ig] * pow(cHD, 2.0));
			GUnmarriedH1[ia][ig] = cHD / (1.0 + cHD + 0.5 * PartnerEffectNew[0][ig] * pow(cHD, 2.0));
			GUnmarriedH2[ia][ig] = 1.0 - GUnmarriedH0[ia][ig] - GUnmarriedH1[ia][ig];
			GUnmarriedL0[ia][ig] = 1.0 / (1.0 + cLD);
			GUnmarriedL1[ia][ig] = cLD / (1.0 + cLD);
			GMarriedH1[ia][ig] = 1.0 / (1.0 + cHD * PartnerEffectNew[1][ig]);
			GMarriedH2[ia][ig] = 1.0 - GMarriedH1[ia][ig];
		}
	}

	//int seed7 = 9927; // Note that I've arbitrarily fixed the seeds for now.
	//if(FixedPeriod==0){
	//	seed7 += CurrSim;}
	//CRandomMersenne rg7(seed7);
	for (ic = 0; ic<tpp; ic++){
		r[ic] = rg.Random();
	}

	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].VirginInd == 1){
			Register[ic].CurrPartners = 0;
		}
		else{
			ia = Register[ic].AgeGroup - 2;
			ig = Register[ic].SexInd;
			if (Register[ic].RiskGroup == 1){
				if (Register[ic].MarriedInd == 1){
					if (r[ic]<GMarriedH1[ia][ig]){
						Register[ic].CurrPartners = 1;
					}
					else{
						Register[ic].CurrPartners = 2;
					}
				}
				else{
					if (r[ic]<GUnmarriedH0[ia][ig]){
						Register[ic].CurrPartners = 0;
					}
					else if (r[ic]<1.0 - GUnmarriedH2[ia][ig]){
						Register[ic].CurrPartners = 1;
					}
					else{
						Register[ic].CurrPartners = 2;
					}
				}
			}
			else{
				if (Register[ic].MarriedInd == 1){
					Register[ic].CurrPartners = 1;
				}
				else{
					if (r[ic]<GUnmarriedL0[ia][ig]){
						Register[ic].CurrPartners = 0;
					}
					else{
						Register[ic].CurrPartners = 1;
					}
				}
			}
		}
	}

	// Next calculate sexual mixing in ST partnerships (adapted from code in TSHISA)
	DesiredSTpartners[0][0] = 0;
	DesiredSTpartners[0][1] = 0;
	DesiredSTpartners[1][0] = 0;
	DesiredSTpartners[1][1] = 0;
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].CurrPartners>0){
			STpartners = Register[ic].CurrPartners - Register[ic].MarriedInd;
			ib = Register[ic].RiskGroup - 1;
			ig = Register[ic].SexInd;
			DesiredSTpartners[ib][ig] += STpartners;
		}
	}

	DesiredPartnerRiskM[0][0] = (1.0 - AssortativeM) + AssortativeM *
		DesiredSTpartners[0][1] / (DesiredSTpartners[0][1] + DesiredSTpartners[1][1]);
	DesiredPartnerRiskM[0][1] = 1.0 - DesiredPartnerRiskM[0][0];
	DesiredPartnerRiskM[1][1] = (1.0 - AssortativeM) + AssortativeM *
		DesiredSTpartners[1][1] / (DesiredSTpartners[0][1] + DesiredSTpartners[1][1]);
	DesiredPartnerRiskM[1][0] = 1.0 - DesiredPartnerRiskM[1][1];
	DesiredPartnerRiskF[0][0] = (1.0 - AssortativeF) + AssortativeF *
		DesiredSTpartners[0][0] / (DesiredSTpartners[0][0] + DesiredSTpartners[1][0]);
	DesiredPartnerRiskF[0][1] = 1.0 - DesiredPartnerRiskF[0][0];
	DesiredPartnerRiskF[1][1] = (1.0 - AssortativeF) + AssortativeF *
		DesiredSTpartners[1][0] / (DesiredSTpartners[0][0] + DesiredSTpartners[1][0]);
	DesiredPartnerRiskF[1][0] = 1.0 - DesiredPartnerRiskF[1][1];

	for (ii = 0; ii<2; ii++){
		for (ij = 0; ij<2; ij++){
			if (AllowBalancing == 1){
				AdjSTrateM[ii][ij] = (GenderEquality * DesiredSTpartners[ij][1] *
					DesiredPartnerRiskF[ij][ii] + (1.0 - GenderEquality) *
					DesiredSTpartners[ii][0] * DesiredPartnerRiskM[ii][ij]);
				AdjSTrateF[ii][ij] = (GenderEquality * DesiredSTpartners[ii][1] *
					DesiredPartnerRiskF[ii][ij] + (1.0 - GenderEquality) *
					DesiredSTpartners[ij][0] * DesiredPartnerRiskM[ij][ii]);
			}
			else{
				AdjSTrateM[ii][ij] = DesiredSTpartners[ii][0] * DesiredPartnerRiskM[ii][ij];
				AdjSTrateF[ii][ij] = DesiredSTpartners[ii][1] * DesiredPartnerRiskF[ii][ij];
			}
		}
		ActualPropnSTH[ii][0] = AdjSTrateM[ii][0] / (AdjSTrateM[ii][0] + AdjSTrateM[ii][1]);
		ActualPropnSTH[ii][1] = AdjSTrateF[ii][0] / (AdjSTrateF[ii][0] + AdjSTrateF[ii][1]);
	}

	// Assign partner ID(s) to each individual in a ST relationship

	//int seed8 = 2717; // Note that I've arbitrarily fixed the seeds for now.
	//if(FixedPeriod==0){
	//	seed8 += CurrSim;}
	//CRandomMersenne rg8(seed8);
	for (ic = 0; ic<tpp; ic++){ // It's theoretically possible that # ST partnerships may
		// be > InitPop, but very unlikely.
		r[ic] = rg.Random();
	}
	//int seed9 = 2115; 
	//if(FixedPeriod==0){
	//	seed9 += CurrSim;}
	//CRandomMersenne rg9(seed9);
	for (ic = 0; ic<tpp; ic++){
		rprisk[ic] = rg.Random();
	}
	//int seed10 = 7192; 
	//if(FixedPeriod==0){
	//	seed10 += CurrSim;}
	//CRandomMersenne rg10(seed10);
	for (ic = 0; ic<tpp; ic++){
		rpID[ic] = rg.Random();
	}

	MaxNewPartnerInd = 0;
	ii = 0;
	while (MaxNewPartnerInd == 0 && ii<tpp){
		UnpartneredCount = 0;
		for (ic = 0; ic<tpp; ic++){
			if (Register[ic].CurrPartners == 1 && Register[ic].IDprimary == 0){
				UnpartneredCount += 1;
			}
			if (Register[ic].CurrPartners == 2 && Register[ic].IDprimary == 0){
				UnpartneredCount += 2;
			}
			if (Register[ic].CurrPartners == 2 && Register[ic].IDprimary>0
				&& Register[ic].ID2ndary == 0){
				UnpartneredCount += 1;
			}
		}
		TotalSelectionProb = 0.0;
		found = 0;
		for (ic = 0; ic<tpp && !found; ic++){
			if (Register[ic].CurrPartners == 1 && Register[ic].IDprimary == 0){
				ij = 1;
			}
			else if (Register[ic].CurrPartners == 2 && Register[ic].IDprimary == 0){
				ij = 2;
			}
			else if (Register[ic].CurrPartners == 2 && Register[ic].ID2ndary == 0){
				ij = 1;
			}
			else{
				ij = 0;
			}
			if (ij > 0){
				temp1 = TotalSelectionProb;
				TotalSelectionProb += 1.0 * ij / UnpartneredCount;
				if (r[ii] <= TotalSelectionProb && r[ii]>temp1){
					partnerrisk = rprisk[ii];
					partnerID = rpID[ii];
					ChooseSTpartner(ic + 1, partnerrisk, partnerID);
					found = 1;
				}
			}
		}
		ii += 1;
	}
	// For people in ST relationships who haven't been assigned a partner ID, reduce # partners.
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].CurrPartners>0 && Register[ic].IDprimary == 0){
			Register[ic].CurrPartners = 0;
		}
		else if (Register[ic].CurrPartners == 2 && Register[ic].ID2ndary == 0){
			Register[ic].CurrPartners = 1;
		}
	}

	// Determine rate at which men visit CSWs (copied from code in TSHISA)
	lambda = GammaMeanFSW / pow(GammaStdDevFSW, 2.0);
	alpha = GammaMeanFSW * lambda;
	for (ia = 0; ia<16; ia++){
		AgeEffectFSWcontact[ia] = pow(lambda, alpha) * pow(ia*5.0 + 2.5, alpha - 1.0) *
			exp(-lambda * (ia*5.0 + 2.5));
		FSWcontactBase[ia] = 0.0;
	}
	MalePop15to49 = 0.0;
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].SexInd == 0 && Register[ic].RiskGroup == 1 && Register[ic].VirginInd == 0){
			ia = Register[ic].AgeGroup - 2;
			if (Register[ic].CurrPartners == 0){
				FSWcontactBase[ia] += AgeEffectFSWcontact[ia] * PartnerEffectFSWcontact[0];
			}
			if (Register[ic].CurrPartners == 1 && Register[ic].MarriedInd == 0){
				FSWcontactBase[ia] += AgeEffectFSWcontact[ia] * PartnerEffectFSWcontact[1];
			}
			if (Register[ic].CurrPartners == 1 && Register[ic].MarriedInd == 1){
				FSWcontactBase[ia] += AgeEffectFSWcontact[ia] * PartnerEffectFSWcontact[2];
			}
			if (Register[ic].CurrPartners == 2 && Register[ic].MarriedInd == 0){
				FSWcontactBase[ia] += AgeEffectFSWcontact[ia] * PartnerEffectFSWcontact[3];
			}
			if (Register[ic].CurrPartners == 2 && Register[ic].MarriedInd == 1){
				FSWcontactBase[ia] += AgeEffectFSWcontact[ia] * PartnerEffectFSWcontact[4];
			}
		}
		if (Register[ic].SexInd == 0 && Register[ic].AgeGroup >= 3 && Register[ic].AgeGroup<10){
			MalePop15to49 += 1.0;
		}
	}
	BaseFSWdemand = 0.0;
	/*for (ia = 0; ia<16; ia++){
	BaseFSWdemand += FSWcontactBase[ia];}
	FSWcontactConstant = MeanFSWcontacts / (pow(lambda, alpha) * pow(11.0, alpha - 1.0) * exp(-lambda * 11.0));*/
	for (ia = 1; ia<8; ia++){
		BaseFSWdemand += FSWcontactBase[ia];
	}
	FSWcontactConstant = MeanFSWcontacts * MalePop15to49 / BaseFSWdemand;
	BaseFSWdemand += FSWcontactBase[0];
	for (ia = 8; ia<16; ia++){
		BaseFSWdemand += FSWcontactBase[ia];
	}
	DesiredFSWcontacts = BaseFSWdemand * FSWcontactConstant;
	RequiredNewFSW = DesiredFSWcontacts / AnnNumberClients;

	// Assign women to commercial sex work
	for (ia = 0; ia<16; ia++){
		EligibleForCSW[ia] = 0;
	}
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].SexInd == 1 && Register[ic].RiskGroup == 1 &&
			Register[ic].VirginInd == 0 && Register[ic].CurrPartners == 0){
			ia = Register[ic].AgeGroup - 2;
			EligibleForCSW[ia] += 1;
		}
	}

	//int seed11 = 3099; // Note that I've arbitrarily fixed the seeds for now.
	//if(FixedPeriod==0){
	//	seed11 += CurrSim;}
	//CRandomMersenne rg11(seed11);
	for (ic = 0; ic<tpp; ic++){
		r[ic] = rg.Random();
	}
	TotCurrFSW = 0;
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].SexInd == 1 && Register[ic].RiskGroup == 1 &&
			Register[ic].VirginInd == 0 && Register[ic].CurrPartners == 0){
			ia = Register[ic].AgeGroup - 2;
			if (r[ic] < RequiredNewFSW * InitFSWageDbn[ia] / EligibleForCSW[ia]){
				Register[ic].FSWind = 1;
				CSWregister[TotCurrFSW] = ic + 1;
				TotCurrFSW += 1;
			}
		}
	}

	// I haven't copied the code from TSHISA for calculating the FSWentry array values,
	// because at baseline these would be very sensitive to the relative numbers of unpartnered
	// women in the high risk group at each age, and it would make more sense to calculate
	// the FSWentry values at the start of each sexual behaviour cycle (i.e. so that the results
	// are less sensitive to the baseline numbers of unpartnered women by age).

	// Assign PartnerRateAdj

	if (AllowPartnerRateAdj == 1){
		//int seed12 = 1015; // Note that I've arbitrarily fixed the seeds for now.
		//if(FixedPeriod==0){
		//	seed12 += CurrSim;}
		//CRandomMersenne rg12(seed12);
		for (ic = 0; ic<tpp; ic++){
			r[ic] = rg.Random();
		}
		ind = 2;
		// Note that the following formulas for a and b apply only when the gamma mean is 1.
		a = 1.0 / pow(SDpartnerRateAdj, 2.0);
		b = a;
		for (ic = 0; ic<tpp; ic++){
			p = r[ic];
			q = 1 - r[ic];
			cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
			Register[ic].PartnerRateAdj = x;
		}
	}
	else{
		for (ic = 0; ic<tpp; ic++){
			Register[ic].PartnerRateAdj = 1.0;
		}
	}
}

void Pop::ChooseLTpartner(int ID, double rnd1, double rnd2)
{
	// This function is used ONLY at baseline, to match married individuals.

	int ic, ia;
	int iage, page; // individual's age group - 2 (so that 0 ==> age 10-14), partner's age group - 2
	int igender, pgender; // individual's gender and partner's gender
	int irisk, prisk; // individual's risk group -1 (so that 0 ==> high risk), partner's risk
	int TotalRisk[2]; // total potential LT partners by risk group
	int EligibleByAge[16];
	double WeightsByAge[16];
	double Normalizer, TotalSelectionProb2;

	iage = Register[ID - 1].AgeGroup - 2;
	igender = Register[ID - 1].SexInd;
	pgender = 1 - igender;
	irisk = Register[ID - 1].RiskGroup - 1;

	int tpp = Register.size();

	// (a) First choose risk group of LT partner
	// =========================================

	TotalRisk[0] = 0;
	TotalRisk[1] = 0;
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].MarriedInd == 1 && Register[ic].IDprimary == 0 && Register[ic].SexInd == pgender){
			if (Register[ic].RiskGroup == 1){
				TotalRisk[0] += 1;
			}
			else{
				TotalRisk[1] += 1;
			}
		}
	}
	if (rnd1 < ActualPropnLTH[irisk][igender]){
		if (TotalRisk[0] > 0){ prisk = 1; }
		else if (TotalRisk[1] > 0){ prisk = 2; }
		else{ MaxNewPartnerInd = 1; }
	}
	else{
		if (TotalRisk[1] > 0){ prisk = 2; }
		else if (TotalRisk[0] > 0){ prisk = 1; }
		else{ MaxNewPartnerInd = 1; }
	}

	// (b) Next choose individual partner
	// ==================================

	// (b1) Calculate age weights
	if (MaxNewPartnerInd == 0){
		for (ia = 0; ia<16; ia++){
			EligibleByAge[ia] = 0;
		}
		for (ic = 0; ic<tpp; ic++){
			if (Register[ic].MarriedInd == 1 && Register[ic].IDprimary == 0 &&
				Register[ic].SexInd == pgender && Register[ic].RiskGroup == prisk){
				ia = Register[ic].AgeGroup - 2; // -2 because age groups start at 10
				EligibleByAge[ia] += 1;
			}
		}
		Normalizer = 0.0;
		for (ia = 0; ia<16; ia++){
			if (igender == 0 && EligibleByAge[ia] > 0){
				WeightsByAge[ia] = AgePrefM[iage][ia] / EligibleByAge[ia];
				Normalizer += AgePrefM[iage][ia];
			}
			else if (igender == 1 && EligibleByAge[ia] > 0){
				WeightsByAge[ia] = AgePrefF[iage][ia] / EligibleByAge[ia];
				Normalizer += AgePrefF[iage][ia];
			}
			else{
				WeightsByAge[ia] = 0.0;
			}
		}
		for (ia = 0; ia<16; ia++){
			WeightsByAge[ia] = WeightsByAge[ia] / Normalizer;
		}
	}

	// (b2) Calculate individual weights
	if (MaxNewPartnerInd == 0 && Normalizer>0.0){
		TotalSelectionProb2 = 0.0;
		for (ic = 0; ic<tpp; ic++){
			if (Register[ic].MarriedInd == 1 && Register[ic].IDprimary == 0 &&
				Register[ic].SexInd == pgender && Register[ic].RiskGroup == prisk){
				page = Register[ic].AgeGroup - 2;
				TotalSelectionProb2 += WeightsByAge[page];
				Register[ic].CumSelectionProb = TotalSelectionProb2;
			}
		}
		if (rnd2>TotalSelectionProb2){ // Due to rounding there is a very small chance
			// that rnd2 may be > TotalSelectionProb2
			rnd2 = TotalSelectionProb2;
		}
	}

	// (b3) Choose partner
	if (MaxNewPartnerInd == 0){
		if (Normalizer == 0){ // The only available partners are in ineligble age groups
			Register[ID - 1].MarriedInd = 0;
		}
		else{
			for (ic = 0; ic<tpp; ic++){
				if (Register[ic].MarriedInd == 1 && Register[ic].IDprimary == 0 &&
					Register[ic].SexInd == pgender && Register[ic].RiskGroup == prisk
					&& rnd2 <= Register[ic].CumSelectionProb){
					Register[ID - 1].IDprimary = ic + 1;
					Register[ic].IDprimary = ID;
					break;
				}
			}
		}
	}
	else{
		Register[ID - 1].MarriedInd = 0;
	}
}

void Pop::ChooseSTpartner(int ID, double rnd1, double rnd2)
{
	// This function is used ONLY at baseline, to match individuals in short-term relationships.

	int ic, ia;
	int iage, page; // individual's age group - 2 (so that 0 ==> age 10-14), partner's age group - 2
	int igender, pgender; // individual's gender and partner's gender
	int irisk, prisk; // individual's risk group -1 (so that 0 ==> high risk), partner's risk
	int TotalRisk[2]; // total potential ST partners by risk group
	int EligibleByAge[16];
	double WeightsByAge[16];
	double Normalizer, TotalSelectionProb2;

	iage = Register[ID - 1].AgeGroup - 2;
	igender = Register[ID - 1].SexInd;
	pgender = 1 - igender;
	irisk = Register[ID - 1].RiskGroup - 1;

	// (a) First choose risk group of ST partner
	// =========================================

	TotalRisk[0] = 0;
	TotalRisk[1] = 0;
	for (ic = 0; ic<InitPop; ic++){
		if (Register[ic].CurrPartners>0 && Register[ic].IDprimary == 0 && Register[ic].SexInd == pgender){
			if (Register[ic].RiskGroup == 1){
				TotalRisk[0] += 1;
			}
			else{
				TotalRisk[1] += 1;
			}
		}
		if (Register[ic].CurrPartners == 2 && Register[ic].ID2ndary == 0 && Register[ic].SexInd == pgender){
			TotalRisk[0] += 1;
		}
	}
	if (rnd1 < ActualPropnSTH[irisk][igender]){
		if (TotalRisk[0] > 0){ prisk = 1; }
		else if (TotalRisk[1] > 0){ prisk = 2; }
		else{ MaxNewPartnerInd = 1; }
	}
	else{
		if (TotalRisk[1] > 0){ prisk = 2; }
		else if (TotalRisk[0] > 0){ prisk = 1; }
		else{ MaxNewPartnerInd = 1; }
	}

	// (b) Next choose individual partner
	// ==================================

	// (b1) Calculate age weights
	if (MaxNewPartnerInd == 0){
		for (ia = 0; ia<16; ia++){
			EligibleByAge[ia] = 0;
		}
		for (ic = 0; ic<InitPop; ic++){
			if (Register[ic].CurrPartners>0 && Register[ic].IDprimary == 0 &&
				Register[ic].SexInd == pgender && Register[ic].RiskGroup == prisk){
				ia = Register[ic].AgeGroup - 2; // -2 because age groups start at 10
				EligibleByAge[ia] += 1;
			}
			if (Register[ic].CurrPartners == 2 && Register[ic].ID2ndary == 0 &&
				Register[ic].SexInd == pgender && Register[ic].RiskGroup == prisk){
				ia = Register[ic].AgeGroup - 2;
				EligibleByAge[ia] += 1;
			}
		}
		Normalizer = 0.0;
		for (ia = 0; ia<16; ia++){
			if (igender == 0 && EligibleByAge[ia] > 0){
				WeightsByAge[ia] = AgePrefM[iage][ia] / EligibleByAge[ia];
				Normalizer += AgePrefM[iage][ia];
			}
			else if (igender == 1 && EligibleByAge[ia] > 0){
				WeightsByAge[ia] = AgePrefF[iage][ia] / EligibleByAge[ia];
				Normalizer += AgePrefF[iage][ia];
			}
			else{
				WeightsByAge[ia] = 0.0;
			}
		}
		for (ia = 0; ia<16; ia++){
			WeightsByAge[ia] = WeightsByAge[ia] / Normalizer;
		}
	}

	// (b2) Calculate individual weights
	if (MaxNewPartnerInd == 0 && Normalizer>0.0){
		TotalSelectionProb2 = 0.0;
		for (ic = 0; ic<InitPop; ic++){
			if (Register[ic].CurrPartners>0 && Register[ic].SexInd == pgender &&
				Register[ic].RiskGroup == prisk){
				page = Register[ic].AgeGroup - 2;
				if (Register[ic].IDprimary == 0){
					TotalSelectionProb2 += WeightsByAge[page];
				}
				if (Register[ic].CurrPartners == 2 && Register[ic].ID2ndary == 0){
					TotalSelectionProb2 += WeightsByAge[page];
				}
				Register[ic].CumSelectionProb = TotalSelectionProb2;
			}
		}
		if (rnd2>TotalSelectionProb2){ // Due to rounding there is a very small chance
			// that rnd2 may be > TotalSelectionProb2
			rnd2 = TotalSelectionProb2;
		}
	}

	// (b3) Choose partner
	if (MaxNewPartnerInd == 0){
		if (Normalizer == 0.0){ // The only available partners are in ineligble age groups
			if (Register[ID - 1].CurrPartners>0 && Register[ID - 1].IDprimary == 0){
				Register[ID - 1].CurrPartners = 0;
			}
			else if (Register[ID - 1].CurrPartners == 2 && Register[ID - 1].ID2ndary == 0){
				Register[ID - 1].CurrPartners = 1;
			}
			else{
				cout << "Error 1:indiv shouldn't be assigned a new partner" << endl;
			}
		}
		else{
			for (ic = 0; ic<InitPop; ic++){
				if (Register[ic].CurrPartners>0 && Register[ic].SexInd == pgender &&
					Register[ic].RiskGroup == prisk && rnd2 <= Register[ic].CumSelectionProb){
					if (Register[ID - 1].IDprimary == 0){
						Register[ID - 1].IDprimary = ic + 1;
					}
					else{
						Register[ID - 1].ID2ndary = ic + 1;
					}
					if (Register[ic].IDprimary == 0){
						Register[ic].IDprimary = ID;
					}
					else{
						Register[ic].ID2ndary = ID;
					}
					break;
				}
			}
		}
	}
	else{ // No partners available in either risk group
		if (Register[ID - 1].CurrPartners>0 && Register[ID - 1].IDprimary == 0){
			Register[ID - 1].CurrPartners = 0;
		}
		else if (Register[ID - 1].CurrPartners == 2 && Register[ID - 1].ID2ndary == 0){
			Register[ID - 1].CurrPartners = 1;
		}
		else{
			cout << "Error 2:indiv shouldn't be assigned a new partner" << endl;
		}
	}
}

void Pop::AssignHIV()
{
	int ic, InitHigh15to49, InitHighMarried, InitHighUnmarried;
	int CumHighRisk1, CumHighRisk2, index1, nextID;

	//int seed = 6810; // Note that I've arbitrarily fixed the seeds for now.
	//if(FixedPeriod==0){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for (ic = 0; ic<InitPop; ic++){
		r[ic] = rg.Random();
	}
	//HPV
	for (ic = 0; ic<InitPop; ic++){
		Register[ic].ARTweeks = 0;
		Register[ic].ARTstage = 0;
	}

	for (ic = 0; ic<InitPop; ic++){
		if (Register[ic].RiskGroup == 1 && Register[ic].VirginInd == 0 && Register[ic].AgeGroup >= 3 &&
			Register[ic].AgeGroup<10 && r[ic]<InitHIVprevHigh && ConstantInitialHIV == 0){
			Register[ic].HIVstage = 2; // To be consistent with TSHISA
			Register[ic].DateInfect = StartYear; // Assuming they acquired HIV 6 months prior
			// to the start of the simulation 
		}
		else{
			Register[ic].HIVstage = 0;
		}
	}
	if (ConstantInitialHIV == 1){
		InitHigh15to49 = 0;
		InitHighMarried = 0;
		for (ic = 0; ic<InitPop; ic++){
			if (Register[ic].RiskGroup == 1 && Register[ic].VirginInd == 0 && Register[ic].AgeGroup >= 3 &&
				Register[ic].AgeGroup<10){
				InitHigh15to49 += 1;
				if (Register[ic].MarriedInd == 1){
					InitHighMarried += 1;
				}
			}
		}
		InitHighUnmarried = InitHigh15to49 - InitHighMarried;
		vector<int> indices1(InitHighUnmarried);
		vector<int> indices2(InitHighMarried);
		CumHighRisk1 = 0;
		CumHighRisk2 = 0;
		for (ic = 0; ic<InitPop; ic++){
			if (Register[ic].RiskGroup == 1 && Register[ic].VirginInd == 0 && Register[ic].AgeGroup >= 3 &&
				Register[ic].AgeGroup<10){
				// Define indices vectors to contain IDs of people in the high risk group
				if (Register[ic].MarriedInd == 0){
					indices1[CumHighRisk1] = ic + 1;
					CumHighRisk1 += 1;
				}
				else{
					indices2[CumHighRisk2] = ic + 1;
					CumHighRisk2 += 1;
				}
			}
		}
		ic = 0;
		while (ic < (InitHIVprevHigh * InitHighUnmarried)){
			index1 = static_cast<int> (r[ic] * indices1.size() );
			nextID = indices1[index1];
			Register[nextID - 1].HIVstage = 2; // To be consistent with TSHISA
			Register[nextID - 1].DateInfect = StartYear;
			indices1[index1] = indices1[indices1.size() - 1];
			indices1.pop_back();
			ic += 1;
		}
		while (ic < (InitHIVprevHigh * InitHigh15to49)){
			index1 = static_cast<int> (r[ic] * indices2.size());
			nextID = indices2[index1];
			Register[nextID - 1].HIVstage = 2; // To be consistent with TSHISA
			Register[nextID - 1].DateInfect = StartYear;
			indices2[index1] = indices2[indices2.size() - 1];
			indices2.pop_back();
			ic += 1;
		}
	}

	// Assign SuscepHIVadj
	int ind;
	double x,  a, b, p, q; //,y;

	//int seed2 = 3814; // Note that I've arbitrarily fixed the seeds for now.
	//if(FixedPeriod==0){
	//	seed2 += CurrSim;}
	//CRandomMersenne rg2(seed2);
	for (ic = 0; ic<InitPop; ic++){
		r[ic] = rg.Random();
	}
	for (ic = 0; ic<InitPop; ic++){
		if (AllowHIVsuscepAdj == 1){
			ind = 2;
			// Note that the following formulas for a and b apply only when the gamma mean is 1.
			a = 1.0 / pow(SDsuscepHIVadj, 2.0);
			b = a;
			p = r[ic];
			q = 1 - r[ic];
			cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
			Register[ic].SuscepHIVadj = x;
		}
		else{
			Register[ic].SuscepHIVadj = 1.0;
		}
	}
}

void Pop::AssignHIV1990()
{
	int ia, ig, ic, CumNewHIV, NewID;
	double ExpectedHIV, Temp, HighPrevByAge[7][2];

	//int seed = 6810; // Note that I've arbitrarily fixed the seeds for now.
	//if (FixedPeriod <= 5){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	int tpp = Register.size();
	for (ic = 0; ic < tpp; ic++){
		r2[ic] = rg.Random();
	}
	for (ia = 0; ia < 7; ia++){
		for (ig = 0; ig < 2; ig++){
			HighPrevByAge[ia][ig] = HlabisaRatio[ia][ig] * InitHIVprevHigh;
		}
	}
	//HPV
	for (ic = 0; ic<InitPop; ic++){
		Register[ic].ARTweeks = 0;
		Register[ic].ARTstage = 0;
	}
	
	if (ConstantInitialHIV == 0){ // Allowing stochastic variation in # infections
		for (ic = 0; ic < tpp; ic++){
			if (Register[ic].RiskGroup == 1 && Register[ic].VirginInd == 0 &&
				Register[ic].AgeGroup >= 3 && Register[ic].AgeGroup < 10 && Register[ic].AliveInd == 1){
				ia = Register[ic].AgeGroup - 3;
				ig = Register[ic].SexInd;
				if (r2[ic] < HighPrevByAge[ia][ig]){
					Register[ic].HIVstage = 2; // To be consistent with TSHISA
					Register[ic].DateInfect = 1990.0; // Assuming they acquired HIV 
					// 6 months prior to mid-1990 
				}
				else{
					Register[ic].HIVstage = 0;
				}
			}
			else{
				Register[ic].HIVstage = 0;
			}
		}
	}
	
	if (ConstantInitialHIV == 1){ // Initial # HIV infections is fixed
		ExpectedHIV = 0.0;
		for (ic = 0; ic<tpp; ic++){
			if (Register[ic].RiskGroup == 1 && Register[ic].VirginInd == 0 && Register[ic].AgeGroup >= 3 &&
				Register[ic].AgeGroup<10 && Register[ic].AliveInd == 1){
				ia = Register[ic].AgeGroup - 3;
				ig = Register[ic].SexInd;
				ExpectedHIV += HighPrevByAge[ia][ig];
				Register[ic].CumSelectionProb = ExpectedHIV;
			}
		}
		CumNewHIV = 0;
		while (CumNewHIV < ExpectedHIV){
			for (ic = 0; ic < tpp; ic++){
				Temp = Register[ic].CumSelectionProb / ExpectedHIV;
				if (r2[CumNewHIV]<Temp && Register[ic].RiskGroup == 1 &&
					Register[ic].VirginInd == 0 && Register[ic].AgeGroup >= 3 &&
					Register[ic].AgeGroup < 10 && Register[ic].AliveInd == 1){
					NewID = ic;
					break;
				}
			}
			CumNewHIV += 1;
			Register[NewID].HIVstage = 2; // To be consistent with TSHISA
			Register[NewID].DateInfect = 1990.0; // Assuming they acquired HIV 6 months prior to mid-1990
		}
	}

	// Assign SuscepHIVadj
	int ind;
	double x,  a, b, p, q; //y;

	//int seed2 = 3814; // Note that I've arbitrarily fixed the seeds for now.
	//if (FixedPeriod <= 5){
	//	seed2 += CurrSim;}
	//CRandomMersenne rg2(seed2);
	
	for (ic = 0; ic<tpp; ic++){
		r2[ic] = rg.Random();
	}
	for (ic = 0; ic<tpp; ic++){
		if (AllowHIVsuscepAdj == 1){
			ind = 2;
			// Note that the following formulas for a and b apply only when the gamma mean is 1.
			a = 1.0 / pow(SDsuscepHIVadj, 2.0);
			b = a;
			p = r2[ic];
			q = 1 - r2[ic];
			cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
			Register[ic].SuscepHIVadj = x;
		}
		else{
			Register[ic].SuscepHIVadj = 1.0;
		}
	}
}

void Pop::AssignSTIs()
{
	int ic, id, group, offset, IDprimary, ID2ndary, xx;
	int RiskPrimary, Risk2ndary; // risk groups of primary and 2ndary partners
	double yy;
	
	//ofstream file("testCIN3.txt", std::ios::app);
	
	//int seed = 781; // Note that I've arbitrarily fixed the seeds for now.
	//if(FixedPeriod==0){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for (ic = 0; ic<InitPop; ic++){
		for (id = 0; id<50; id++){ // changed back to 59...
			rSTI[ic][id] = rg.Random();
		}
	}
	for (ic = 0; ic<InitPop; ic++){
		// First determine the behaviour group
		if (Register[ic].VirginInd == 1){
			if (Register[ic].RiskGroup == 1){ group = 0; }
			else{ group = 1; }
		}
		else{
			if (Register[ic].RiskGroup == 1){ group = 2; }
			else{ group = 14; }
			if (Register[ic].CurrPartners == 1){
				group += 1;
				IDprimary = Register[ic].IDprimary;
				RiskPrimary = Register[IDprimary - 1].RiskGroup;
				if (Register[ic].MarriedInd == 1){
					group += 2;
				}
				group += RiskPrimary - 1;
			}
			if (Register[ic].CurrPartners == 2){
				group += 5;
				IDprimary = Register[ic].IDprimary;
				ID2ndary = Register[ic].ID2ndary;
				RiskPrimary = Register[IDprimary - 1].RiskGroup;
				Risk2ndary = Register[ID2ndary - 1].RiskGroup;
				if (Register[ic].MarriedInd == 0){
					group += RiskPrimary + Risk2ndary - 2;
				}
				else{
					group += (RiskPrimary - 1) * 2 + Risk2ndary + 2;
				}
			}
			if (Register[ic].FSWind == 1){ group = 19; }
		}
		// Determine offset
		offset = group * 16 + Register[ic].AgeGroup - 2;
		//Assign number of lifetime partners
		if(Register[ic].AgeGroup>1){
			if(Register[ic].SexInd==0)	{Register[ic].LifetimePartners = LTP[offset][0];}
			else {Register[ic].LifetimePartners = LTP[offset][1];} 
		}
		else {Register[ic].LifetimePartners=0;}
		// Assign status for individual STIs
		if (HSVind == 1){
			if (Register[ic].AgeGroup<2){
				Register[ic].HSVstage = 0;
			}
			else if (Register[ic].SexInd == 0){
				Register[ic].HSVstage = HSVtransitionM.GetSTDstage(offset, rSTI[ic][0]);
			}
			else{
				Register[ic].HSVstage = HSVtransitionF.GetSTDstage(offset, rSTI[ic][0]);
			}
		}
		if (TPind == 1){
			if (Register[ic].AgeGroup<2){
				Register[ic].TPstage = 0;
			}
			else if (Register[ic].SexInd == 0){
				Register[ic].TPstage = TPtransitionM.GetSTDstage(offset, rSTI[ic][1]);
			}
			else{
				Register[ic].TPstage = TPtransitionF.GetSTDstage(offset, rSTI[ic][1]);
			}
		}
		if (HDind == 1){
			if (Register[ic].AgeGroup<2){
				Register[ic].HDstage = 0;
			}
			else if (Register[ic].SexInd == 0){
				Register[ic].HDstage = HDtransitionM.GetSTDstage(offset, rSTI[ic][2]);
			}
			else{
				Register[ic].HDstage = HDtransitionF.GetSTDstage(offset, rSTI[ic][2]);
			}
		}
		if (NGind == 1){
			if (Register[ic].AgeGroup<2){
				Register[ic].NGstage = 0;
			}
			else if (Register[ic].SexInd == 0){
				Register[ic].NGstage = NGtransitionM.GetSTDstage(offset, rSTI[ic][3]);
			}
			else{
				Register[ic].NGstage = NGtransitionF.GetSTDstage(offset, rSTI[ic][3]);
			}
		}
		if (CTind == 1){
			if (Register[ic].AgeGroup<2){
				Register[ic].CTstage = 0;
			}
			else if (Register[ic].SexInd == 0){
				Register[ic].CTstage = CTtransitionM.GetSTDstage(offset, rSTI[ic][4]);
			}
			else{
				Register[ic].CTstage = CTtransitionF.GetSTDstage(offset, rSTI[ic][4]);
			}
		}
		if (TVind == 1){
			if (Register[ic].AgeGroup<2){
				Register[ic].TVstage = 0;
			}
			else if (Register[ic].SexInd == 0){
				Register[ic].TVstage = TVtransitionM.GetSTDstage(offset, rSTI[ic][5]);
			}
			else{
				Register[ic].TVstage = TVtransitionF.GetSTDstage(offset, rSTI[ic][5]);
			}
		}
		if (BVind == 1){
			if (Register[ic].AgeGroup<2 || Register[ic].SexInd == 0){
				Register[ic].BVstage = 1;
			}
			else{
				Register[ic].BVstage = BVtransitionF.GetSTDstage(offset, rSTI[ic][6]) + 1;
			}
		}
		if (VCind == 1){
			if (Register[ic].AgeGroup<2 || Register[ic].SexInd == 0){
				Register[ic].VCstage = 0;
			}
			else{
				Register[ic].VCstage = VCtransitionF.GetSTDstage(offset, rSTI[ic][7]);
			}
		}
		if (HPVind == 1){
			Register[ic].GotVacc=0;
			Register[ic].GotVaccOffer=0;
			Register[ic].InScreen=0;
			Register[ic].InWHOScreen=0;
			Register[ic].ScreenCount=0;
			Register[ic].Scr30=0;
			Register[ic].Scr50=0;
			Register[ic].Scr35=0;
			Register[ic].Scr45=0;
			Register[ic].Scr16=0;
			Register[ic].Scr19=0;
			Register[ic].Scr22=0;
			Register[ic].Scr25=0;
			Register[ic].Scr28=0;
			Register[ic].Scr31=0;
			Register[ic].Scr34=0;
			Register[ic].Scr37=0;
			Register[ic].Scr40=0;
			Register[ic].Scr43=0;
			Register[ic].Scr46=0;
			Register[ic].Scr49=0;
			Register[ic].Scr52=0;
			Register[ic].Scr55=0;
			Register[ic].Scr58=0;
			Register[ic].Scr61=0;
			Register[ic].Scr64=0;
			Register[ic].GotTxV=0;

			Register[ic].ScreenResult=0;
			Register[ic].timePassed=0;

			Register[ic].StageIdeath=0;
			Register[ic].StageIIdeath=0;
			Register[ic].StageIIIdeath=0;
			Register[ic].StageIVdeath=0;
			Register[ic].StageIrecover=0;
			Register[ic].StageIIrecover=0;
			Register[ic].StageIIIrecover=0;
			Register[ic].StageIVrecover=0;
			Register[ic].TimeinStageI=0;
			Register[ic].TimeinStageII=0;
			Register[ic].TimeinStageIII=0;
			Register[ic].TimeinStageIV=0;
			Register[ic].ThermalORPap=rSTI[ic][41];

			for (xx = 0; xx < 13; xx++){
				Register[ic].VaccinationStatus[xx]=0;
				if (Register[ic].AgeGroup < 2){
					Register[ic].HPVstage[xx] = 0;
				}
				else if (Register[ic].SexInd == 0){
					Register[ic].HPVstage[xx] = HPVTransM[xx].GetSTDstage(offset, rSTI[ic][8 + xx]);
				}
				else{
					Register[ic].HPVstage[xx] = HPVTransF[xx].GetSTDstage(offset, rSTI[ic][8 + xx]);
				}

				if (Register[ic].HPVstage[xx] == 4 && Register[ic].SexInd == 1){
					Register[ic].AssignTimeinCIN3(Register[ic].AgeGroup, rSTI[ic][21 + xx], xx);
					if(Register[ic].Age50==0){ yy = exp(-pow((Register[ic].TimeinCIN3[xx] / HPVTransF[xx].AveDuration[3]), HPVTransF[xx].CIN3shape));}
					else {yy = exp(-pow((Register[ic].TimeinCIN3[xx] /(HPVTransF[xx].CIN3_50 * HPVTransF[xx].AveDuration[3])), HPVTransF[xx].CIN3shape));}
					//file << yy << " " ;
					yy = yy - yy*rSTI[ic][40];	
					if(Register[ic].Age50==0){ 
						Register[ic].WeibullCIN3[xx] = static_cast<int> (exp((1.0 / HPVTransF[xx].CIN3shape)*log(-log(yy)) + log(HPVTransF[xx].AveDuration[3])));
					}
					else {
						Register[ic].WeibullCIN3[xx] = static_cast<int> (exp((1.0 / HPVTransF[xx].CIN3shape)*log(-log(yy)) + 
						log(HPVTransF[xx].CIN3_50 * HPVTransF[xx].AveDuration[3])));
					}
					
					//cout << Register[ic].WeibullCIN3[xx] << endl;
					//file << CurrSim << "	" << ic << "	" << Register[ic].AgeGroup << "	" << xx << " " << rSTI[ic][21 + xx] << " " ;
					//file <<  Register[ic].TimeinCIN3[xx] / 52 << "	" << yy << " " << Register[ic].WeibullCIN3[xx] / 52 << endl;
				}
			}
				
		}
	}
	//file.close();	
}

void Pop::GetPopPyramid()
{
	int ia, ig, ic, ih, ix;

	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 ){
			ia = Register[ic].AgeGroup;
			ix = int (Register[ic].AgeExact);
			if(Register[ic].HIVstage==0) {ih=0;}
			else if(Register[ic].HIVstage==5 ){ih=1;} 
			else {ih=2;}

			//ig = Register[ic].SexInd;
			ig = CurrYear-StartYear;
			if(Register[ic].SexInd==1){
				PopPyramid[18*ih + ia][ig] += 1;
				PopPyramidAll[ia][ig] += 1;
				if(Register[ic].HIVstage==5){PopPyramidAllART[ig] += 1;}
				if(Register[ic].AgeExact>=CatchUpAgeMIN && Register[ic].AgeExact<CatchUpAgeMAX) {PopPyramid9[ih][ig] += 1;}
				if(Indiv::AnyHPV(Register[ic].HPVstage,  {0,1}, {1,2,3,4})){
					HPVprevVT[18*ih + ia][ig] += 1;
				}
				if(Register[ic].HPVstatus==1){
					HPVprevAll[18*ih + ia][ig] += 1;
				}
				if(Register[ic].TrueStage>1 ){
					CIN2prev[18*ih + ia][ig] += 1;
				}
				if(Register[ic].AgeExact<60){
					PopPyramid61[ix][ig] += 1;
					if(Register[ic].GotVacc==1){PopVaxx61[ix][ig] += 1;}
				} else {
					PopPyramid61[60][ig] += 1;
					if(Register[ic].GotVacc==1){PopVaxx61[60][ig] += 1;}
				}
			}
			else {
				PopPyramidMale[18*ih + ia][ig] += 1;
			}
			if(Register[ic].SexInd==0 && Register[ic].AgeExact>=CatchUpAgeMIN && Register[ic].AgeExact<CatchUpAgeMAX) {PopPyramid9[3 + ih][ig] += 1;}

			//std::cout << PopPyramid[ia][ig] << std::endl;
		}
	}
	
}

void Pop::GetBehavPyramid()
{
	int ia, ig, is, ic;

	for (ia = 0; ia<18; ia++){
		for (ig = 0; ig<2; ig++){
			for (is = 0; is<3; is++){
				BehavPyramid[ia][ig * 3 + is] = 0;
			}
		}
	}
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1){
			ia = Register[ic].AgeGroup;
			ig = Register[ic].SexInd;
			if (Register[ic].VirginInd == 1){
				BehavPyramid[ia][ig * 3] += 1;
			}
			else if (Register[ic].MarriedInd == 1){
				BehavPyramid[ia][ig * 3 + 2] += 1;
			}
			else{
				BehavPyramid[ia][ig * 3 + 1] += 1;
			}
		}
	}
}

void Pop::GetPrefMatrix(int type)
{
	int ic;
	int iage, page; // indiv age group (-2) and partner age group (-2)
	int partnerID;

	for (iage = 0; iage<16; iage++){
		for (page = 0; page<16; page++){
			PrefMatrix[iage][page] = 0;
		}
	}
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 1){
			iage = Register[ic].AgeGroup - 2;
			if (type == 1 && Register[ic].MarriedInd == 1){
				partnerID = Register[ic].IDprimary;
				page = Register[partnerID - 1].AgeGroup - 2;
				PrefMatrix[iage][page] += 1;
			}
			if (type == 2 && Register[ic].MarriedInd == 0 &&
				Register[ic].IDprimary>0){
				partnerID = Register[ic].IDprimary;
				page = Register[partnerID - 1].AgeGroup - 2;
				PrefMatrix[iage][page] += 1;
			}
			if (type == 2 && Register[ic].ID2ndary>0){
				partnerID = Register[ic].ID2ndary;
				page = Register[partnerID - 1].AgeGroup - 2;
				PrefMatrix[iage][page] += 1;
			}
		}
	}
}

void Pop::GetNumberPartners()
{
	int ia, ig, is, ic;

	for (ia = 0; ia<16; ia++){
		for (ig = 0; ig<2; ig++){
			for (is = 0; is<5; is++){
				NumberPartners[ia][ig * 5 + is] = 0;
			}
		}
	}
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup >= 2 &&
			Register[ic].VirginInd == 0){
			ia = Register[ic].AgeGroup - 2;
			ig = Register[ic].SexInd;
			is = Register[ic].CurrPartners;
			if (Register[ic].MarriedInd == 1){
				NumberPartners[ia][ig * 5 + 2 + is] += 1;
			}
			else{
				NumberPartners[ia][ig * 5 + is] += 1;
			}
		}
	}
}

void Pop::GetTotSex()
{
	int ia, ig, ii, il, ic;
	int CSWID, CSWsexActsTot[MaxCSWs][2];

	for (ia = 0; ia<16; ia++){
		for (ig = 0; ig<2; ig++){
			for (il = 0; il<4; il++){
				TotSex[ia][ig * 4 + il] = 0;
			}
		}
	}
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup >= 2 &&
			Register[ic].VirginInd == 0){
			ia = Register[ic].AgeGroup - 2;
			ig = Register[ic].SexInd;
			if (Register[ic].MarriedInd == 1){
				TotSex[ia][ig * 4 + 2] += Register[ic].UVIprimary;
				TotSex[ia][ig * 4 + 3] += Register[ic].PVIprimary;
			}
			else{
				TotSex[ia][ig * 4] += Register[ic].UVIprimary;
				TotSex[ia][ig * 4 + 1] += Register[ic].PVIprimary;
			}
			TotSex[ia][ig * 4] += Register[ic].UVI2ndary;
			TotSex[ia][ig * 4 + 1] += Register[ic].PVI2ndary;
			if (ig == 0 && Register[ic].RiskGroup == 1){
				if (Register[ic].UVICSW + Register[ic].PVICSW >0){
					CSWID = Register[ic].IDofCSW;
					for (ii = 0; ii<TotCurrFSW; ii++){
						if (CSWregister[ii] == CSWID){
							CSWsexActsTot[ii][0] += Register[ic].UVICSW;
							CSWsexActsTot[ii][1] += Register[ic].PVICSW;
							break;
						}
					}
				}
			}
		}
	}

	// Commercial sex acts
	/*for(ic=0; ic<TotCurrFSW; ic++){
	cout<<"Unprotected sex acts with CSW "<<ic+1<<": "<<CSWsexActsTot[ic][0]<<endl;
	cout<<"Protected sex acts with CSW "<<ic+1<<": "<<CSWsexActsTot[ic][1]<<endl;
	}*/
}

void Pop::GetNumbersByHIVstage()
{
	int ic, iy, ig, is;

	iy = CurrYear - StartYear;

	for (ig = 0; ig < 2; ig++){
		for (is = 0; is < 6; is++){
			AdultHIVstageTrend[iy][ig * 6 + is] = 0;
		}
	}
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup >= 3){
			// Note that this gives total pop at ages 15+
			ig = Register[ic].SexInd;
			is = Register[ic].HIVstage;
			AdultHIVstageTrend[iy][ig * 6 + is] += 1;
		}
	}
}

void Pop::GetInitHIVprevH()
{
	int Positive[2], Negative[2];
	int ic, ig;

	Positive[0] = 0;
	Negative[0] = 0;
	Positive[1] = 0;
	Negative[1] = 0;
	for (ic = 0; ic<InitPop; ic++){
		if (Register[ic].RiskGroup == 1 && Register[ic].AgeGroup>2 &&
			Register[ic].AgeGroup<10 && Register[ic].VirginInd == 0){
			ig = Register[ic].SexInd;
			if (Register[ic].HIVstage == 0){
				Negative[ig] += 1;
			}
			else{
				Positive[ig] += 1;
			}
		}
	}
	cout << "HIV-positive males in high risk group, 15-49: " << Positive[0] << endl;
	cout << "HIV-negative males in high risk group, 15-49: " << Negative[0] << endl;
	cout << "HIV-positive females in high risk group, 15-49: " << Positive[1] << endl;
	cout << "HIV-negative females in high risk group, 15-49: " << Negative[1] << endl;
}

void Pop::GetHHprevByAge()
{
	int ic, ia, ig;
	int Total[9][2], Positive[9][2];

	for (ia = 0; ia < 9; ia++){
		for (ig = 0; ig < 2; ig++){
			Total[ia][ig] = 0;
			Positive[ia][ig] = 0;
		}
	}
	int tpp = Register.size();
	for (ic = 0; ic < tpp; ic++){
		if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 12 && Register[ic].AliveInd == 1){
			ia = Register[ic].AgeGroup - 3;
			ig = Register[ic].SexInd;
			Total[ia][ig] += 1;
			if (Register[ic].HIVstage>0){
				Positive[ia][ig] += 1;
			}
		}
	}
	for (ia = 0; ia < 9; ia++){
		for (ig = 0; ig < 2; ig++){
			if (CurrYear == 2005){
				PrevHH2005.out[CurrSim - 1][ia + 9 * ig] = 1.0 * Positive[ia][ig] / Total[ia][ig];
			}
			if (CurrYear == 2008){
				PrevHH2008.out[CurrSim - 1][ia + 9 * ig] = 1.0 * Positive[ia][ig] / Total[ia][ig];
			}
			if (CurrYear == 2012){
				PrevHH2012.out[CurrSim - 1][ia + 9 * ig] = 1.0 * Positive[ia][ig] / Total[ia][ig];
			}
			if (CurrYear == 2017){
				PrevHH2017.out[CurrSim - 1][ia + 9 * ig] = 1.0 * Positive[ia][ig] / Total[ia][ig];
			}
		}
	}
}

void Pop::GetANCprevByAge()
{
	int ia, ic, is;
	double TempFert, Total[7], Positive[7], temp1, temp2;

	for (ia = 0; ia < 7; ia++){
		Total[ia] = 0.0;
		Positive[ia] = 0.0;
	}
	int tpp = Register.size();
	for (ic = 0; ic < tpp; ic++){
		if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 10 &&
			Register[ic].AliveInd == 1 && Register[ic].SexInd == 1){
			ia = Register[ic].AgeGroup - 3;
			is = Register[ic].HIVstage;
			TempFert = Register[ic].NonHIVfertRate;
			if (is>0){
				TempFert *= RelHIVfertility[is - 1];
			}
			Total[ia] += TempFert;
			if (Register[ic].HIVstage>0){
				Positive[ia] += TempFert;
			}
		}
	}
	PrevPreg15.out[CurrSim - 1][CurrYear - 1990] = Positive[0] / Total[0];
	PrevPreg20.out[CurrSim - 1][CurrYear - 1990] = Positive[1] / Total[1];
	PrevPreg25.out[CurrSim - 1][CurrYear - 1990] = Positive[2] / Total[2];
	PrevPreg30.out[CurrSim - 1][CurrYear - 1990] = Positive[3] / Total[3];
	PrevPreg35.out[CurrSim - 1][CurrYear - 1990] = Positive[4] / Total[4];
	PrevPreg40.out[CurrSim - 1][CurrYear - 1990] = (Positive[5] + Positive[6]) / (Total[5] + Total[6]);
	temp1 = 0.0;
	temp2 = 0.0;
	for (ia = 0; ia < 7; ia++){
		temp1 += Total[ia];
		temp2 += Positive[ia];
	}
	PrevPregTotal.out[CurrSim - 1][CurrYear - 1990] = temp2 / temp1;
}

void Pop::GetInitSTIprev() //don't bother adding HPV - I write out initial HPV stages by age
{
	// Returns prevalence in 15-49 year old males and females

	int TotPop[2], TotHSV[2], TotTP[2], TotHD[2], TotNG[2];
	int TotCT[2], TotTV[2], TotBV, TotVC;
	int ic, ig;

	for (ig = 0; ig<2; ig++){
		TotPop[ig] = 0;
		TotHSV[ig] = 0;
		TotTP[ig] = 0;
		TotHD[ig] = 0;
		TotNG[ig] = 0;
		TotCT[ig] = 0;
		TotTV[ig] = 0;
		TotBV = 0;
		TotVC = 0;
	}
	for (ic = 0; ic<InitPop; ic++){
		if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup<10){
			ig = Register[ic].SexInd;
			TotPop[ig] += 1;
			if (Register[ic].HSVstage>0){ TotHSV[ig] += 1; }
			if (Register[ic].TPstage>0 && Register[ic].TPstage<5){ TotTP[ig] += 1; }
			if (Register[ic].HDstage>0 && Register[ic].HDstage<3){ TotHD[ig] += 1; }
			if (Register[ic].NGstage>0 && Register[ic].NGstage<3){ TotNG[ig] += 1; }
			if (Register[ic].CTstage>0 && Register[ic].CTstage<3){ TotCT[ig] += 1; }
			if (Register[ic].TVstage>0 && Register[ic].TVstage<3){ TotTV[ig] += 1; }
			if (Register[ic].BVstage>2 && ig == 1){ TotBV += 1; }
			if (Register[ic].VCstage>0 && ig == 1){ TotVC += 1; }
			
			}
	}
	cout << "Male HSV-2 prevalence, 15-49: " << TotHSV[0] << endl;
	cout << "Female HSV-2 prevalence, 15-49: " << TotHSV[1] << endl;
	cout << "Male TP prevalence, 15-49: " << TotTP[0] << endl;
	cout << "Female TP prevalence, 15-49: " << TotTP[1] << endl;
	cout << "Male HD prevalence, 15-49: " << TotHD[0] << endl;
	cout << "Female HD prevalence, 15-49: " << TotHD[1] << endl;
	cout << "Male NG prevalence, 15-49: " << TotNG[0] << endl;
	cout << "Female NG prevalence, 15-49: " << TotNG[1] << endl;
	cout << "Male CT prevalence, 15-49: " << TotCT[0] << endl;
	cout << "Female CT prevalence, 15-49: " << TotCT[1] << endl;
	cout << "Male TV prevalence, 15-49: " << TotTV[0] << endl;
	cout << "Female TV prevalence, 15-49: " << TotTV[1] << endl;
	cout << "Female BV prevalence, 15-49: " << TotBV << endl;
	cout << "Female VC prevalence, 15-49: " << TotVC << endl;
	cout << "Male pop 15-49: " << TotPop[0] << endl;
	cout << "Female pop 15-49: " << TotPop[1] << endl;
}

void Pop::GetCurrSTIprev()
{
	// Calculates prevalence in 15-49 year old males and females
	// Very similar to GetInitSTIprev function.

	int TotPop[2], TotHSV[2], TotTP[2], TotHD[2], TotNG[2];
	int TotCT[2], TotTV[2], TotBV, TotVC, TotHIV[2], TotHPV[2], TotART[2];
	int ic, ig, xx;

	for (ig = 0; ig<2; ig++){
		TotPop[ig] = 0;
		TotHIV[ig] = 0;
		/*TotHSV[ig] = 0;
		TotTP[ig] = 0;
		TotHD[ig] = 0;
		TotNG[ig] = 0;
		TotCT[ig] = 0;
		TotTV[ig] = 0;
		TotBV = 0;
		TotVC = 0;*/
		TotHPV[ig] = 0;
		TotART[ig] = 0;
	}
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AgeGroup>2  && Register[ic].AliveInd == 1){ //
			ig = Register[ic].SexInd;
			TotPop[ig] += 1;
			if (Register[ic].HIVstage>0){ TotHIV[ig] += 1; }
			if (Register[ic].HIVstage==5){ TotART[ig] += 1; }
		}
	}
	HIVprev15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * (TotHIV[0]) / (TotPop[0]);
	HIVprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHIV[1] / TotPop[1];
	/*for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AgeGroup>2  && Register[ic].AliveInd == 1&& Register[ic].AgeGroup<10){ //
			ig = Register[ic].SexInd;
			TotPop[ig] += 1;
			if (Register[ic].HIVstage>0){ TotHIV[ig] += 1; }
			if (Register[ic].HIVstage==5){ TotART[ig] += 1; }
		}
	}
	HIVprev15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * (TotHIV[0]) / (TotPop[0]);
	HIVprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHIV[1] / TotPop[1];
	HIVprev15to49all.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * (TotHIV[1]+TotHIV[0]) / (TotPop[1]+TotPop[0]);
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AgeGroup>=10  && Register[ic].AliveInd == 1){ //
			ig = Register[ic].SexInd;
			if (Register[ic].HIVstage>0){ TotHIV[ig] += 1; }
			if (Register[ic].HIVstage==5){ TotART[ig] += 1; }
		}
	}*/
	ARTcov15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotART[0] / TotHIV[0];
	ARTcov15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotART[1] / TotHIV[1];
	/*HSVprev15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHSV[0] / TotPop[0];
	HSVprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHSV[1] / TotPop[1];
	TPprev15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotTP[0] / TotPop[0];
	TPprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotTP[1] / TotPop[1];
	HDprev15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHD[0] / TotPop[0];
	HDprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHD[1] / TotPop[1];
	NGprev15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotNG[0] / TotPop[0];
	NGprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotNG[1] / TotPop[1];
	CTprev15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotCT[0] / TotPop[0];
	CTprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotCT[1] / TotPop[1];
	TVprev15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotTV[0] / TotPop[0];
	TVprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotTV[1] / TotPop[1];
	BVprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotBV / TotPop[1];
	VCprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotVC / TotPop[1];
	//HPVprev15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHPV[0] / TotPop[0];
	//HPVprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHPV[1] / TotPop[1];

	// Calculate the STI prevalence required for likelihood calculation purposes
	// (I haven't added code for HD, BV or VC as yet.)
	if (HSVcalib == 1){
		GetANCprev(&HSVtransitionF, 1);
		GetHHprev(&HSVtransitionF, 1);
		GetHHprev(&HSVtransitionM, 1);
		GetFPCprev(&HSVtransitionF, 1);
		GetCSWprev(&HSVtransitionF, 1);
	}
	if (TPcalib == 1){
		GetANCprev(&TPtransitionF, 2);
		GetHHprev(&TPtransitionF, 2);
		GetHHprev(&TPtransitionM, 2);
		GetFPCprev(&TPtransitionF, 2);
		GetCSWprev(&TPtransitionF, 2);
	}
	if (NGcalib == 1){
		GetANCprev(&NGtransitionF, 4);
		GetHHprev(&NGtransitionF, 4);
		GetHHprev(&NGtransitionM, 4);
		GetFPCprev(&NGtransitionF, 4);
		GetCSWprev(&NGtransitionF, 4);
	}
	if (CTcalib == 1){
		GetANCprev(&CTtransitionF, 5);
		GetHHprev(&CTtransitionF, 5);
		GetHHprev(&CTtransitionM, 5);
		GetFPCprev(&CTtransitionF, 5);
		GetCSWprev(&CTtransitionF, 5);
	}
	if (TVcalib == 1){
		GetANCprev(&TVtransitionF, 6);
		GetHHprev(&TVtransitionF, 6);
		GetHHprev(&TVtransitionM, 6);
		GetFPCprev(&TVtransitionF, 6);
		GetCSWprev(&TVtransitionF, 6);
	}*/
	if (HPVcalib == 1){
				
		GetHHprev(&HPVtransitionCC, 9);
		GetHHNprev(&HPVtransitionCC, 9);
		GetFPCprev(&HPVtransitionCC, 9);
		GetNOARTprev(&HPVtransitionCC, 9);
		GetONARTprev(&HPVtransitionCC, 9);
	}
}

double Pop::GetQstatistic(int Cmatrix[3][3], int MatDim)
{
	int ir, ic;
	double RowTotal, trace, Qstatistic;

	trace = 0.0;
	for (ir = 0; ir<MatDim; ir++){
		RowTotal = 0.0;
		for (ic = 0; ic<MatDim; ic++){
			RowTotal += Cmatrix[ir][ic];
		}
		trace += Cmatrix[ir][ir] / RowTotal;
	}

	Qstatistic = (trace - 1.0) / (MatDim - 1.0);

	return Qstatistic;
}

void Pop::GetSTIconcordance()
{
	int ic, pID, pID2, iSTI, pSTI;
	int HIVmatrix[3][3], HSVmatrix[3][3], TPmatrix[3][3], HDmatrix[3][3];
	int NGmatrix[3][3], CTmatrix[3][3], TVmatrix[3][3];

	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AgeGroup>2 && Register[ic].AliveInd == 1 &&
			Register[ic].CurrPartners>0 && Register[ic].SexInd == 1){
			pID = Register[ic].IDprimary;
			pID2 = Register[ic].ID2ndary;
			if (HIVind == 1){
				if (Register[ic].HIVstage>0){ iSTI = 1; }
				else{ iSTI = 0; }
				if (Register[pID - 1].HIVstage>0){ pSTI = 1; }
				else{ pSTI = 0; }
				HIVmatrix[iSTI][pSTI] += 1;
				if (Register[ic].CurrPartners>1){
					if (Register[pID2 - 1].HIVstage>0){ pSTI = 1; }
					else{ pSTI = 0; }
					HIVmatrix[iSTI][pSTI] += 1;
				}
			}
			if (HSVind == 1){
				if (Register[ic].HSVstage>0){ iSTI = 1; }
				else{ iSTI = 0; }
				if (Register[pID - 1].HSVstage>0){ pSTI = 1; }
				else{ pSTI = 0; }
				HSVmatrix[iSTI][pSTI] += 1;
				if (Register[ic].CurrPartners>1){
					if (Register[pID2 - 1].HSVstage>0){ pSTI = 1; }
					else{ pSTI = 0; }
					HSVmatrix[iSTI][pSTI] += 1;
				}
			}
			if (TPind == 1){
				if (Register[ic].TPstage>0){ iSTI = 1; }
				else{ iSTI = 0; }
				if (Register[ic].TPstage>4){ iSTI = 0; } // immune
				if (Register[pID - 1].TPstage>0){ pSTI = 1; }
				else{ pSTI = 0; }
				if (Register[pID - 1].TPstage>4){ pSTI = 0; }
				TPmatrix[iSTI][pSTI] += 1;
				if (Register[ic].CurrPartners>1){
					if (Register[pID2 - 1].TPstage>0){ pSTI = 1; }
					else{ pSTI = 0; }
					if (Register[pID2 - 1].TPstage>4){ pSTI = 0; }
					TPmatrix[iSTI][pSTI] += 1;
				}
			}
			if (HDind == 1){
				if (Register[ic].HDstage>0){ iSTI = 1; }
				else{ iSTI = 0; }
				if (Register[ic].HDstage>2){ iSTI = 0; } // immune
				if (Register[pID - 1].HDstage>0){ pSTI = 1; }
				else{ pSTI = 0; }
				if (Register[pID - 1].HDstage>2){ pSTI = 0; }
				HDmatrix[iSTI][pSTI] += 1;
				if (Register[ic].CurrPartners>1){
					if (Register[pID2 - 1].HDstage>0){ pSTI = 1; }
					else{ pSTI = 0; }
					if (Register[pID2 - 1].HDstage>2){ pSTI = 0; }
					HDmatrix[iSTI][pSTI] += 1;
				}
			}
			if (NGind == 1){
				if (Register[ic].NGstage>0){ iSTI = 1; }
				else{ iSTI = 0; }
				if (Register[ic].NGstage>2){ iSTI = 0; } // immune
				if (Register[pID - 1].NGstage>0){ pSTI = 1; }
				else{ pSTI = 0; }
				if (Register[pID - 1].NGstage>2){ pSTI = 0; }
				NGmatrix[iSTI][pSTI] += 1;
				if (Register[ic].CurrPartners>1){
					if (Register[pID2 - 1].NGstage>0){ pSTI = 1; }
					else{ pSTI = 0; }
					if (Register[pID2 - 1].NGstage>2){ pSTI = 0; }
					NGmatrix[iSTI][pSTI] += 1;
				}
			}
			if (CTind == 1){
				if (Register[ic].CTstage>0){ iSTI = 1; }
				else{ iSTI = 0; }
				if (Register[ic].CTstage>2){ iSTI = 0; } // immune
				if (Register[pID - 1].CTstage>0){ pSTI = 1; }
				else{ pSTI = 0; }
				if (Register[pID - 1].CTstage>2){ pSTI = 0; }
				CTmatrix[iSTI][pSTI] += 1;
				if (Register[ic].CurrPartners>1){
					if (Register[pID2 - 1].CTstage>0){ pSTI = 1; }
					else{ pSTI = 0; }
					if (Register[pID2 - 1].CTstage>2){ pSTI = 0; }
					CTmatrix[iSTI][pSTI] += 1;
				}
			}
			if (TVind == 1){
				if (Register[ic].TVstage>0){ iSTI = 1; }
				else{ iSTI = 0; }
				if (Register[ic].TVstage>2){ iSTI = 0; } // immune
				if (Register[pID - 1].TVstage>0){ pSTI = 1; }
				else{ pSTI = 0; }
				if (Register[pID - 1].TVstage>2){ pSTI = 0; }
				TVmatrix[iSTI][pSTI] += 1;
				if (Register[ic].CurrPartners>1){
					if (Register[pID2 - 1].TVstage>0){ pSTI = 1; }
					else{ pSTI = 0; }
					if (Register[pID2 - 1].TVstage>2){ pSTI = 0; }
					TVmatrix[iSTI][pSTI] += 1;
				}
			}
		}
	}

	/*if(HIVind==1){
	HIVconcordance.out[CurrSim-1][CurrYear-StartYear] = GetQstatistic(HIVmatrix, 2);}
	if(HSVind==1){
	HSVconcordance.out[CurrSim-1][CurrYear-StartYear] = GetQstatistic(HSVmatrix, 2);}
	if(TPind==1){
	TPconcordance.out[CurrSim-1][CurrYear-StartYear] = GetQstatistic(TPmatrix, 2);}
	if(HDind==1){
	HDconcordance.out[CurrSim-1][CurrYear-StartYear] = GetQstatistic(HDmatrix, 2);}
	if(NGind==1){
	NGconcordance.out[CurrSim-1][CurrYear-StartYear] = GetQstatistic(NGmatrix, 2);}
	if(CTind==1){
	CTconcordance.out[CurrSim-1][CurrYear-StartYear] = GetQstatistic(CTmatrix, 2);}
	if(TVind==1){
	TVconcordance.out[CurrSim-1][CurrYear-StartYear] = GetQstatistic(TVmatrix, 2);}*/
}

void Pop::SavePopPyramid(const char* filout)
{
	int ia, is;
	ostringstream s;

		if (process_num >0){
			s << process_num << "_"  << filout;
		}
		else{
			s <<  filout;
		}
				
		string path = "./output/" + s.str();
	ofstream file(path.c_str()); // Converts s to a C string

	/*for (ia = 0; ia<18; ia++){
		file << right << PopPyramid[ia][0] << "	" << PopPyramid[ia][1] << endl;
	}*/
	for (ia = 0; ia < 54; ia++){
		for (is = 0; is<136; is++){
			file << right << PopPyramid[ia][is] << "	";
		}
		file << endl;
	}
	/*for (ia = 0; ia < 54; ia++){
		for (is = 0; is<136; is++){
			file << right << PopPyramidMale[ia][is] << "	";
		}
		file << endl;
	}*/
	for (ia = 0; ia < 54; ia++){
		for (is = 0; is<136; is++){
			file << right << HPVprevVT[ia][is] << "	";
		}
		file << endl;
	}
	for (ia = 0; ia < 54; ia++){
		for (is = 0; is<136; is++){
			file << right << HPVprevAll[ia][is] << "	";
		}
		file << endl;
	}
	for (ia = 0; ia < 54; ia++){
		for (is = 0; is<136; is++){
			file << right << CIN2prev[ia][is] << "	";
		}
		file << endl;
	}
	/*for (ia = 0; ia < 61; ia++){
		for (is = 0; is<136; is++){
			file << right << PopPyramid61[ia][is] << "	";
		}
		file << endl;
	}
	for (ia = 0; ia < 61; ia++){
		for (is = 0; is<136; is++){
			file << right << PopVaxx61[ia][is] << "	";
		}
		file << endl;
	}*/
	
	file.close();
	
}

void Pop::SaveBehavPyramid(char* filout)
{
	int ia, is;
	ofstream file(filout);

	for (ia = 0; ia<18; ia++){
		for (is = 0; is<6; is++){
			file << right << BehavPyramid[ia][is] << "	";
		}
		file << endl;
	}
	file.close();
}

void Pop::SavePrefMatrix(char* filout)
{
	int ia, is;
	ofstream file(filout);

	for (ia = 0; ia<16; ia++){
		for (is = 0; is<16; is++){
			file << right << PrefMatrix[ia][is] << "	";
		}
		file << endl;
	}
	file.close();
}

void Pop::SaveNumberPartners(char* filout)
{
	int ia, is;
	ofstream file(filout);

	for (ia = 0; ia<16; ia++){
		for (is = 0; is<10; is++){
			file << right << NumberPartners[ia][is] << "	";
		}
		file << endl;
	}
	file.close();
}

void Pop::SaveLifetimePartners(char* filout)
{
	int ic;
	double AgeNow;
	ofstream file(filout);
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 0 && Register[ic].CurrPartners>0 &&
			Register[ic].AgeGroup >= 3 && Register[ic].AgeGroup<10){
			AgeNow = 0.5 + CurrYear - Register[ic].DOB;
			file << right << ic << "	" << AgeNow << "	" << Register[ic].LifetimePartners << endl;
		}
	}
	file.close();
}

void Pop::SaveTotSex(char* filout)
{
	int ia, is;
	ofstream file(filout);

	for (ia = 0; ia<16; ia++){
		for (is = 0; is<8; is++){
			file << right << TotSex[ia][is] << "	";
		}
		file << endl;
	}
	file.close();
}

void Pop::SaveMarriageRates(char* filout)
{
	int ia, is;
	ofstream file(filout);

	for (ia = 0; ia<16; ia++){
		for (is = 0; is<4; is++){
			file << right << AgeEffectMarriage[ia][is] << "	";
		}
		file << endl;
	}
	file.close();
}

void Pop::SaveAdultHIVstage(char* filout)
{
	int iy, is;
	ofstream file(filout);

	for (iy = 0; iy <= 40; iy++){
		for (is = 0; is<12; is++){
			file << right << AdultHIVstageTrend[iy][is] << "	";
		}
		file << endl;
	}
	file.close();
}

void Pop::OneYear()
{
	int ii;

	// Update age group and reset flow variables
	//if (FixedUncertainty == 1){ResetFlow();  }
	ResetFlow();

	if (CurrYear>StartYear){
		UpdateAgeGroup();
	}
	
	// Calculate all stock variable outputs at the start of the year
	//GetNumbersByHIVstage();
	//GetNumbersByHPVstage();
	
	if (UpdateStart == 1 && CurrYear==2030){ GetNumbersByHPVstageAge();}
	//if(WHOscenario==1 && WHOvacc==1 && CurrYear==ImplementYR){
	if((AdministerMassTxV == 1  || AdministerMassTxVtoART == 1) && (CampaignYear[CurrYear-StartYear] ==1)){//CampaignYear[CurrYear-StartYear]==1){//#CurrYear==ImplementYR){
		RSApop.AssignVacc2024();
	}
	GetPopPyramid();
	
	//CalcModelCoverage();

	//GetSTIconcordance();
	/*if ((CurrYear == 2005 || CurrYear == 2008 || CurrYear == 2012 || CurrYear == 2017)){ // && HIVcalib == 1){
		GetHHprevByAge();
	}*/

	CalcNonHIVmort();

	CalcNonHIVfert();
	
	GetCurrSTIprev();
	if (FixedUncertainty == 1 && CCcalib==0) {
		GetCurrHPVstage();
	}

	if (CurrYear >= 1990 && OneType == 1){
		GetCurrHPVprev(WhichType); 
	}
	if (CurrYear >= 2000 && targets == 1){
		HitTargets(WhichType);
	}
	
	if (CurrYear <= 2012 && CurrYear >= 1990 && HIVcalib == 1){
		GetANCprevByAge();
	}
	if(CurrYear == 2000 && GetMacD==1){
		GetMacDprev();
	}


	CalcMTCTrates();
	UpdateProbCure();
	UpdateSTDtransitionProbs();
	UpdateCondomUse();
		
	// Project changes in behaviour, disease & demography over year
	BehavCycleCount = 0;
	for (ii = 0; ii<CycleS; ii++){

		OneBehavCycle();
	    //if(ii==24) { GetPopPyramid();}
		//if (CreateCohort == 1 && CurrYear > 1999 && (ii == 0 || ii == 12 || ii == 24 || ii == 36)){
		//if (CreateCohort == 1 && CurrYear >= 2010 && (ii == 0) ){

		if (CreateCohort == 1 && (CurrYear >= 2010) && (CurrYear <= 2030) && (ii == 24)  ){ //
			ProspectiveCohort(ii);
		}	
	}

	//GetPopPyramid();
	CalcModelCoverage();
	//{
	//	double total_pap = 0, total_hpv = 0, total_thermal = 0, total_colp = 0;
	//	for(int ig=0; ig<54; ig++){
	//		total_pap += RSApop.ModelCoverage[ig][CurrYear - StartYear];
	//		total_hpv += RSApop.ModelHPVCoverage[ig][CurrYear - StartYear];
	///		total_thermal += RSApop.ModelThermalCoverage[ig][CurrYear - StartYear];
	//		total_colp += RSApop.ModelColpCoverage[ig][CurrYear - StartYear];
	//	}
	//	cout << "Year " << CurrYear << ": Pap smears = " << total_pap
	//	     << ", HPV DNA tests = " << total_hpv
	//	     << ", Thermal ablations = " << total_thermal
	//	     << ", Colposcopies = " << total_colp << endl;
	//}
	CurrYear += 1;

}

void Pop::ResetFlow()
{
	int iy;
	//NewHIV.out[CurrSim - 1][CurrYear - StartYear] = 0;
	NewHSV.out[CurrSim - 1][CurrYear - StartYear] = 0;
	NewTP.out[CurrSim - 1][CurrYear - StartYear] = 0;
	NewNG.out[CurrSim - 1][CurrYear - StartYear] = 0;
	NewCT.out[CurrSim - 1][CurrYear - StartYear] = 0;
	NewTV.out[CurrSim - 1][CurrYear - StartYear] = 0;
	for (iy=0; iy<54; iy++){
		RSApop.NewScreen[iy][CurrYear - StartYear]=0;
		RSApop.NewReflex[iy][CurrYear - StartYear]=0;
		RSApop.NewTxV[iy][CurrYear - StartYear]=0;
		RSApop.NewHPVScreen[iy][CurrYear - StartYear]=0;
		RSApop.NewColposcopy[iy][CurrYear - StartYear]=0;
		RSApop.NewLLETZ[iy][CurrYear - StartYear]=0;
		RSApop.NewUnnecessary[iy][CurrYear - StartYear]=0;
		RSApop.NewVAT[iy][CurrYear - StartYear]=0;
		RSApop.NewThermal[iy][CurrYear - StartYear]=0;
		RSApop.GetReferred[iy][CurrYear - StartYear]=0;
	}
	for (iy=0; iy<36; iy++){
		RSApop.NewVACC[iy][CurrYear - StartYear]=0;
	}	
}

void Pop::UpdateAgeGroup()
{
	int ic;
	
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1){
			Register[ic].AgeExact = 1.0 * CurrYear + 0.5 - Register[ic].DOB;
			if (Register[ic].AgeExact < 90.0){
				Register[ic].AgeGroup = static_cast<int>(Register[ic].AgeExact / 5);
			}
			else{
				Register[ic].AgeGroup = 17;
			}
			if(Register[ic].AgeExact>=50 && Register[ic].Age50==0 && Register[ic].SexInd==1){
				for(int xx = 0; xx<13; xx++){
					if(Register[ic].WeibullCIN3[xx]>0  && (Register[ic].HIVstage==0)){ //|| 
						//(Register[ic].HIVstage==5 && Register[ic].ARTstage==0)||
						//(Register[ic].HIVstage==5  && Register[ic].ARTstage==1 && Register[ic].ARTweeks>=104))) {
							Register[ic].WeibullCIN3[xx] *= HPVTransF[xx].CIN3_50;
						}
				}
				Register[ic].Age50 = 1;
			}
		}
	}
}

void Pop::CalcNonHIVmort()
{
	int ic, ia;
	//double AnnMortProb;
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1){
			// Get annual rate
			if (Register[ic].AgeGroup >= 2){
				ia = Register[ic].AgeGroup - 2;
				if (Register[ic].SexInd == 0){
					Register[ic].NonHIVmortProb = NonAIDSmortM[ia][CurrYear - StartYear];
				}
				else{
					Register[ic].NonHIVmortProb = NonAIDSmortF[ia][CurrYear - StartYear];
				}
			}
			else{
				ia = static_cast<int>(CurrYear + 0.5 - Register[ic].DOB);
				if (Register[ic].SexInd == 0){
					Register[ic].NonHIVmortProb = ChildMortM[ia][CurrYear - StartYear];
				}
				else{
					Register[ic].NonHIVmortProb = ChildMortF[ia][CurrYear - StartYear];
				}
			}
			// Convert annual rate into rate per behav cycle
			Register[ic].NonHIVmortProb = 1.0 - pow(1.0 - Register[ic].NonHIVmortProb,
				1.0 / CycleS);
		}
	}
}

void Pop::CalcNonHIVfert()
{
	int ia, ic;
	double VirginSum[7], TotalFemSum[7];

	// First calculate SexuallyExpFert
	for (ia = 0; ia<7; ia++){
		VirginSum[ia] = 0.0;
		TotalFemSum[ia] = 0.0;
		HIVnegFert[ia] = FertilityTable[ia][CurrYear - StartYear];
	}
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 1 &&
			Register[ic].AgeGroup>2 && Register[ic].AgeGroup<10){
			ia = Register[ic].AgeGroup - 3;
			TotalFemSum[ia] += 1.0;
			if (Register[ic].VirginInd == 1){
				VirginSum[ia] += 1.0;
			}
		}
	}
	for (ia = 0; ia<7; ia++){
		SexuallyExpFert[ia] = HIVnegFert[ia] / (1.0 - VirginSum[ia] / TotalFemSum[ia]);
	}

	// Then calculate NonHIVfertRate
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 1 &&
			Register[ic].AgeGroup>2 && Register[ic].AgeGroup<10 &&
			Register[ic].VirginInd == 0){
			ia = Register[ic].AgeGroup - 3;
			Register[ic].NonHIVfertRate = SexuallyExpFert[ia] / CycleS;
		}
		else{
			Register[ic].NonHIVfertRate = 0.0;
		}
	}
}

void Pop::CalcMTCTrates()
{
	int yr;

	yr = CurrYear - StartYear;
	CurrPerinatal = PropnInfectedAtBirth * (1.0 - PMTCTaccess[yr] *
		AcceptScreening * AcceptNVP * RednNVP);
	CurrPostnatal = (1.0 - CurrPerinatal) * PropnInfectedAfterBirth *
		(1.0 - PMTCTaccess[yr] * AcceptScreening * RednFF);
}

void Pop::UpdateProbCure()
{
	/*MaleRxRate = InitMaleRxRate * (1.0 + RxPhaseIn[CurrYear-StartYear]);
	MaleTeenRxRate = InitMaleTeenRxRate * (1.0 + RxPhaseIn[CurrYear-StartYear]);
	FemRxRate = InitFemRxRate * (1.0 + RxPhaseIn[CurrYear-StartYear]);
	FemTeenRxRate = InitFemTeenRxRate * (1.0 + RxPhaseIn[CurrYear-StartYear]);*/
	/*FSWasympRxRate = InitFSWasympRxRate + RxPhaseIn[CurrYear-StartYear] * 0.5 *
	(0.25 - InitFSWasympRxRate);
	FSWasympCure = (InitFSWasympCure * (1.0 - RxPhaseIn[CurrYear-StartYear] * 0.5) *
	InitFSWasympRxRate + 1.0 * RxPhaseIn[CurrYear-StartYear] * 0.5 * 0.25)/
	((1.0 - RxPhaseIn[CurrYear-StartYear] * 0.5) * InitFSWasympRxRate +
	RxPhaseIn[CurrYear-StartYear] * 0.5 * 0.25);*/
	if (HSVind == 1){
		/*HSVtransitionM.CorrectRxWithSM = InitCorrectRxHSV + (0.9 - InitCorrectRxHSV) *
		RxPhaseIn[CurrYear-StartYear];
		HSVtransitionF.CorrectRxWithSM = HSVtransitionM.CorrectRxWithSM;*/
		HSVtransitionM.CalcProbCure();
		HSVtransitionF.CalcProbCure();
	}
	if (TPind == 1){
		/*TPtransitionF.ANCpropnScreened = InitANCpropnScreened + RxPhaseIn[CurrYear-StartYear] *
		(1.0 - InitANCpropnScreened);
		TPtransitionF.ANCpropnTreated = InitANCpropnTreated + RxPhaseIn[CurrYear-StartYear] *
		(1.0 - InitANCpropnTreated);*/
		TPtransitionM.CalcProbCure();
		TPtransitionF.CalcProbCure();
	}
	if (HDind == 1){
		HDtransitionM.CalcProbCure();
		HDtransitionF.CalcProbCure();
	}
	if (NGind == 1){
		/*NGtransitionM.DrugEff = (1.0 - PropnCiproResistant[CurrYear-StartYear] *
		PropnCiproTreated[CurrYear-StartYear]) * InitDrugEffNG;
		NGtransitionF.DrugEff = NGtransitionM.DrugEff;*/
		NGtransitionM.CalcProbCure();
		NGtransitionF.CalcProbCure();
	}
	if (CTind == 1){
		CTtransitionM.CalcProbCure();
		CTtransitionF.CalcProbCure();
	}
	if (TVind == 1){
		//TVtransitionM.CorrectRxWithSM = InitCorrectRxTVM + (0.9 - InitCorrectRxTVM) * 
		//	RxPhaseIn[CurrYear-StartYear];
		TVtransitionM.CalcProbCure();
		TVtransitionF.CalcProbCure();
	}
	if (BVind == 1){
		BVtransitionF.CalcProbCure();
		BVtransitionF.CalcProbPartialCure();
	}
	if (VCind == 1){
		VCtransitionF.CalcProbCure();
		VCtransitionF.CalcProbPartialCure();
	}
}

void Pop::UpdateSTDtransitionProbs()
{
	int  xx; 
	if (HIVind == 1){
		HIVtransitionM.CalcTransitionProbs();
		HIVtransitionF.CalcTransitionProbs();
	}
	if (HSVind == 1){
		/*HSVtransitionM.RecurrenceRate = InitRecurrenceRateM * (1.0 - 0.8 *
		RxPhaseIn[CurrYear-StartYear]);
		HSVtransitionF.RecurrenceRate = InitRecurrenceRateF * (1.0 - 0.8 *
		RxPhaseIn[CurrYear-StartYear]);*/
		HSVtransitionM.CalcTransitionProbs();
		HSVtransitionF.CalcTransitionProbs();
	}
	if (TPind == 1){
		TPtransitionM.CalcTransitionProbs();
		TPtransitionF.CalcTransitionProbs();
	}
	if (HDind == 1){
		HDtransitionM.CalcTransitionProbs();
		HDtransitionF.CalcTransitionProbs();
	}
	if (NGind == 1){
		NGtransitionM.CalcTransitionProbs();
		NGtransitionF.CalcTransitionProbs();
	}
	if (CTind == 1){
		CTtransitionM.CalcTransitionProbs();
		CTtransitionF.CalcTransitionProbs();
	}
	if (TVind == 1){
		TVtransitionM.CalcTransitionProbs();
		TVtransitionF.CalcTransitionProbs();
	}
	if (BVind == 1){
		BVtransitionF.CalcTransitionProbs();
	}
	if (VCind == 1){
		VCtransitionF.CalcTransitionProbs();
	}
	if (HPVind == 1){
		//if (OneType == 0 && targets == 0){
			for (xx = 0; xx < 13; xx++){
				HPVTransM[xx].CalcTransitionProbsM();
				HPVTransF[xx].CalcTransitionProbsF();
				
			}
		//}
		//else{
		//	HPVTransM[WhichType].CalcTransitionProbsM();
		//	HPVTransF[WhichType].CalcTransitionProbsF();
		//}
	}
}

void Pop::UpdateCondomUse()
{
	int ia, ib;
	double Rate15to19[3], x;

	BaselineCondomUse = 0.104; //0.08 + CondomScaling * (0.2 - 0.08);
	//cout << "BaselineCondomUse " << BaselineCondomUse << endl;
	//RatioUltTo1998[0] = 3.0 + CondomScaling * (15.0 - 3.0);
	RatioUltTo1998[0] = 4.6; // 2.0 + CondomScaling * (15.0 - 2.0);
	//RatioUltTo1998[1] = 1.5 + CondomScaling * (7.0 - 1.5);
	RatioUltTo1998[1] = 2.16; //1.2 + CondomScaling * (6.0 - 1.2);
	RatioUltTo1998[2] = 3.8;
	ShapeBehavChange[0] = 4.6; //3.8 + CondomScaling * (2.8 - 3.8);
	ShapeBehavChange[1] = 2.16; //3.6 + CondomScaling * (1.8 - 3.6);
	ShapeBehavChange[2] = 3.18; 
	MedianToBehavChange[0] = 13.0 * exp(-log(log(1.0 - log(RatioInitialTo1998[0]) /
		log(RatioUltTo1998[0])) / log(2.0)) / ShapeBehavChange[0]);
	MedianToBehavChange[1] = 13.0 * exp(-log(log(1.0 - log(RatioInitialTo1998[1]) /
		log(RatioUltTo1998[1])) / log(2.0)) / ShapeBehavChange[1]);

	for (ib = 0; ib<3; ib++){
		x = (BaselineCondomUse / (1.0 - BaselineCondomUse)) * RelEffectCondom[ib] *
			RatioInitialTo1998[ib] * 
			pow(RatioUltTo1998[ib] / RatioInitialTo1998[ib], 
			1.0 - pow(0.5, pow((CurrYear - 1985) / MedianToBehavChange[ib], ShapeBehavChange[ib])))/
			pow(RatioUltTo1998[ib] / exp(0.606) , 
			1.0 - pow(0.5, pow((CurrYear - 1985) / MedianToBehavChange2[ib], 2.0*ShapeBehavChange[ib])));
		Rate15to19[ib] = x / (1.0 + x); //0;//if UpdateStart==1, set to 0
	}
	// Condom usage for females in ST relationships
	for (ia = 0; ia < 16; ia++){
		x = (Rate15to19[0] / (1.0 - Rate15to19[0])) * exp(5.0 * (ia - 1) * AgeEffectCondom[0]);
		CondomUseST[ia][1] = x / (1.0 + x);
		//cout << CurrYear << " " <<	ia << " " << CondomUseST[ia][1] << endl;
	}
	// Condom usage for females in LT relationships
	for (ia = 0; ia < 16; ia++){
		x = (Rate15to19[1] / (1.0 - Rate15to19[1])) * exp(5.0 * (ia - 1) * AgeEffectCondom[1]);
		CondomUseLT[ia][1] = x / (1.0 + x);
	}
	// Condom use in males
	for (ia = 0; ia<16; ia++){
		CondomUseST[ia][0] = 0;
		CondomUseLT[ia][0] = 0;
		for (ib = 0; ib<16; ib++){
			CondomUseST[ia][0] += AgePrefM[ia][ib] * CondomUseST[ib][1];
			CondomUseLT[ia][0] += AgePrefM[ia][ib] * CondomUseLT[ib][1];
		}
	}
	// Condom use in FSW-client relationships
	CondomUseFSW =Rate15to19[2]; //0.1;// if UpdateStart==1, set to 0.1
}

void Pop::GetANCprev(STDtransition* a, int STDind)
{
	int ic, ii, is;
	double numerator, denominator, TempFert;

	if (a->ANClogL.Observations>0){
		numerator = 0.0;
		denominator = 0.0;
		int tpp = Register.size();
		for (ii = 0; ii < tpp; ii++){
			if (Register[ii].AgeGroup>2 && Register[ii].AgeGroup < 10 &&
				Register[ii].AliveInd == 1 && Register[ii].SexInd == 1){
				is = Register[ii].HIVstage;
				TempFert = Register[ii].NonHIVfertRate;
				if (is>0){
					TempFert *= RelHIVfertility[is - 1];
				}
				denominator += TempFert;
				if (STDind == 1 && Register[ii].HSVstage>0){ numerator += TempFert; }
				// Note that for calibration purposes we include the women who are TP-immune but exclude incubation
				if (STDind == 2 && Register[ii].TPstage>1){ numerator += TempFert; }
				if (STDind == 3 && Register[ii].HDstage>0 && Register[ii].HDstage<3){ numerator += TempFert; }
				if (STDind == 4 && Register[ii].NGstage>0 && Register[ii].NGstage<3){ numerator += TempFert; }
				if (STDind == 5 && Register[ii].CTstage>0 && Register[ii].CTstage<3){ numerator += TempFert; }
				if (STDind == 6 && Register[ii].TVstage>0 && Register[ii].TVstage<3){ numerator += TempFert; }
				if (STDind == 7 && Register[ii].BVstage>1){ numerator += TempFert; }
				if (STDind == 8 && Register[ii].VCstage>0){ numerator += TempFert; }
			}
		}
		a->ANClogL.CurrentModelPrev = numerator / denominator;
		// Add code for storing prevalence in current year?
		if (STDind == 1 && FixedUncertainty == 1 && HSVcalib == 1){
			HSVprevANC.out[CurrSim - 1][CurrYear - StartYear] = numerator / denominator;
		}
		if (STDind == 2 && FixedUncertainty == 1 && TPcalib == 1){
			TPprevANC.out[CurrSim - 1][CurrYear - StartYear] = numerator / denominator;
		}
		for (ic = 0; ic<a->ANClogL.Observations; ic++){
			if (a->ANClogL.StudyYear[ic] == CurrYear){
				a->ANClogL.ModelPrev[ic] = a->ANClogL.CurrentModelPrev;
			}
		}
	}
}

void Pop::GetCSWprev(STDtransition* a, int STDind)
{
	int ii;
	double numerator, denominator;

	if (a->CSWlogL.Observations>0){
		numerator = 0.0;
		denominator = 0.0;
		int tpp = Register.size();
		for (ii = 0; ii < tpp; ii++){
			if (Register[ii].FSWind == 1 && Register[ii].AliveInd == 1 &&
				Register[ii].SexInd == 1){
				denominator += 1.0;
				if (STDind == 1 && Register[ii].HSVstage>0){ numerator += 1.0; }
				// Note that for calibration purposes we include the women who are TP-immune but exclude incubation
				if (STDind == 2 && Register[ii].TPstage>1){ numerator += 1.0; }
				if (STDind == 3 && Register[ii].HDstage>0 && Register[ii].HDstage<3){ numerator += 1.0; }
				if (STDind == 4 && Register[ii].NGstage>0 && Register[ii].NGstage<3){ numerator += 1.0; }
				if (STDind == 5 && Register[ii].CTstage>0 && Register[ii].CTstage<3){ numerator += 1.0; }
				if (STDind == 6 && Register[ii].TVstage>0 && Register[ii].TVstage<3){ numerator += 1.0; }
				if (STDind == 7 && Register[ii].BVstage>1){ numerator += 1.0; }
				if (STDind == 8 && Register[ii].VCstage>0){ numerator += 1.0; }
			}
		}
		a->CSWprevUnsmoothed[CurrYear - StartYear] = numerator / denominator;
	}
}

void Pop::OneBehavCycle()
{
	int ic, STDcyclesPerBehavCycle;
	int tpp = Register.size();

	BehavCycleCount += 1;
	
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup >= 2){
			Register[ic].NewStatus = 0;
		}
	}
	
	UpdateBirths();
	
	// Call the OneSTDcycle the appropriate number of times
	STDcyclesPerBehavCycle = CycleD / CycleS;
	STDcycleCount = 0;
	for (ic = 0; ic<STDcyclesPerBehavCycle; ic++){
		OneSTDcycle();
	}
	
	// Calculate non-HIV mortality
	UpdateNonHIVmort();
	
	// Calculate movements between sexual behaviour classes
	UpdateFSWprofile();
	UpdateMarriageIncidence();
	BalanceSexualPartners();
	CalcPartnerTransitions();
}

void Pop::UpdateBirths()
{
	int ic,  is;
	double r[MaxPop];
	double TempFert;

	//int seed = 8821 + CurrYear * 33 + BehavCycleCount * 71;
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for (ic = 0; ic<MaxPop; ic++){
		r[ic] = rg.Random();
	}
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 1){
			is = Register[ic].HIVstage;
			TempFert = Register[ic].NonHIVfertRate;
			if (is>0){
				TempFert *= RelHIVfertility[is - 1];
			}
			if (r[ic]<TempFert){
				NewBirth(ic + 1);
			}
		}
	}
}

void Pop::UpdateNonHIVmort()
{
	int ic;
	double rmort[MaxPop], tt;

	//int seed = 5020 + CurrYear * 27 + BehavCycleCount * 85;
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for (ic = 0; ic<MaxPop; ic++){
		rmort[ic] = rg.Random();
	}
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1){
			if (Register[ic].AgeGroup == 0){
				if (Register[ic].DOB >(CurrYear + 0.5 + (BehavCycleCount - 1.0) / CycleS)){
					// Time from start of behav cycle to DOB (in years)
					tt = Register[ic].DOB - (CurrYear + 0.5 + (BehavCycleCount - 1.0) / CycleS);
					// Time from DOB to end of behav cycle (in behav cycles)
					tt = 1.0 - (tt * CycleS);
					Register[ic].NonHIVmortProb = 1.0 - pow(1.0 -
						Register[ic].NonHIVmortProb, tt);
				}
			}
			if (rmort[ic]<Register[ic].NonHIVmortProb){
				SetToDead(ic + 1);
			}
		}
	}
}

void Pop::UpdateFSWprofile()
{
	int ia, ic,  ij, is;
	double rexit[MaxCSWs], rentry[MaxCSWs];
	double ExitProb, EntryProb, temp, normalizer, DoubleIntDif, IntPart;
	double SingleHRfem[16];
	int CSWID, offset, ActualNewCSW, tempID[MaxCSWs], TotFSWbyAge[16], EntryInd;

	// Firstly randomly determine which sex workers retire
	//int seed1 = 3773 + CurrYear * 34 + BehavCycleCount * 92;
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed1 += CurrSim;}
	//CRandomMersenne rg1(seed1);
	for (ic = 0; ic<TotCurrFSW; ic++){
		rexit[ic] = rg.Random();
	}
	for (ic = 0; ic<MaxCSWs; ic++){
		tempID[ic] = CSWregister[ic];
	}
	offset = 0;
	for (ic = 0; ic<MaxCSWs; ic++){
		CSWID = CSWregister[ic];
		if (CSWID == 0){
			break;
		}
		ia = Register[CSWID - 1].AgeGroup - 2;
		is = Register[CSWID - 1].HIVstage;
		ExitProb = exp(-FSWexit[ia] / CycleS);
		if (is>0){
			ExitProb = pow(ExitProb, HIVeffectFSWexit[is - 1]);
		}
		ExitProb = 1.0 - ExitProb;
		if (rexit[ic]<ExitProb){
			TotCurrFSW = TotCurrFSW - 1;
			Register[CSWID - 1].FSWind = 0;
			for (ij = 0; ij<MaxCSWs - ic - 1 - offset; ij++){
				tempID[ic + ij - offset] = tempID[ic + ij + 1 - offset];
			}
			offset += 1;
		}
		Register[CSWID - 1].NewStatus = 1;
	}
	for (ic = 0; ic<MaxCSWs; ic++){
		CSWregister[ic] = tempID[ic];
	}

	// Secondly get rates of entry into CSW and age-specific weights
	TotCurrFSW = 0;
	DesiredFSWcontacts = 0.0;
	for (ia = 0; ia<16; ia++){
		FSWentry[ia] = 0.0;
		TotFSWbyAge[ia] = 0;
		SingleHRfem[ia] = 0.0;
	}
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup >= 2 && Register[ic].RiskGroup == 1){
			if (Register[ic].SexInd == 1 && Register[ic].CurrPartners == 0 &&
				Register[ic].VirginInd == 0){
				ia = Register[ic].AgeGroup - 2;
				if (Register[ic].FSWind == 0){
					is = Register[ic].HIVstage;
					if (is == 0){
						SingleHRfem[ia] += 1.0;
					}
					else{
						SingleHRfem[ia] += HIVeffectFSWentry[is - 1];
					}
				}
				else{
					TotFSWbyAge[ia] += 1;
					TotCurrFSW += 1;
				}
			}
			if (Register[ic].SexInd == 0 && Register[ic].VirginInd == 0){
				ia = Register[ic].AgeGroup - 2;
				temp = FSWcontactConstant * AgeEffectFSWcontact[ia];
				if (Register[ic].CurrPartners == 0){
					DesiredFSWcontacts += temp * PartnerEffectFSWcontact[0];
				}
				if (Register[ic].CurrPartners == 1 && Register[ic].MarriedInd == 0){
					DesiredFSWcontacts += temp * PartnerEffectFSWcontact[1];
				}
				if (Register[ic].CurrPartners == 1 && Register[ic].MarriedInd == 1){
					DesiredFSWcontacts += temp * PartnerEffectFSWcontact[2];
				}
				if (Register[ic].CurrPartners == 2 && Register[ic].MarriedInd == 0){
					DesiredFSWcontacts += temp * PartnerEffectFSWcontact[3];
				}
				if (Register[ic].CurrPartners == 2 && Register[ic].MarriedInd == 1){
					DesiredFSWcontacts += temp * PartnerEffectFSWcontact[4];
				}
			}
		}
	}
	if (TotCurrFSW<DesiredFSWcontacts / AnnNumberClients){
		normalizer = 0.0;
		for (ia = 0; ia<16; ia++){
			if (SingleHRfem[ia]>0.0){
				FSWentry[ia] = ((DesiredFSWcontacts / AnnNumberClients) *
					InitFSWageDbn[ia] - TotFSWbyAge[ia]) / SingleHRfem[ia];
			}
			else{
				FSWentry[ia] = 0.0;
			}
			if (FSWentry[ia]<0.0){
				FSWentry[ia] = 0.0;
			}
			else{
				normalizer += FSWentry[ia] * SingleHRfem[ia];
			}
		}
	}
	else{
		for (ia = 0; ia<16; ia++){
			FSWentry[ia] = 0.0;
		}
	}

	// Lastly determine which women become sex workers
	if (TotCurrFSW<DesiredFSWcontacts / AnnNumberClients){
		//int seed2 = 6100  + CurrYear * 32 + BehavCycleCount * 73; 
		//if(CurrYear>=StartYear+FixedPeriod){
		//	seed2 += CurrSim;}
		//CRandomMersenne rg2(seed2);
		for (ic = 0; ic<MaxCSWs; ic++){
			rentry[ic] = rg.Random();
		}

		// (a) Determine the number of new sex workers
		RequiredNewFSW = (DesiredFSWcontacts / AnnNumberClients) - TotCurrFSW;
		DoubleIntDif = modf(RequiredNewFSW, &IntPart);
		ActualNewCSW = static_cast<int>(RequiredNewFSW - DoubleIntDif);
		if (rentry[0]<DoubleIntDif){
			ActualNewCSW += 1;
		}

		// (b) Determine which women become sex workers
		temp = 0.0;
		int tpp = Register.size();
		for (ic = 0; ic<tpp; ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup >= 2 && Register[ic].RiskGroup == 1 &&
				Register[ic].SexInd == 1 && Register[ic].VirginInd == 0 &&
				Register[ic].CurrPartners == 0 && Register[ic].FSWind == 0){
				ia = Register[ic].AgeGroup - 2;
				is = Register[ic].HIVstage;
				EntryInd = 0;
				EntryProb = FSWentry[ia] / normalizer;
				if (is>0){
					EntryProb *= HIVeffectFSWentry[is - 1];
				}
				for (ij = 1; ij <= ActualNewCSW; ij++){
					if (rentry[ij]>temp && rentry[ij] <= (temp + EntryProb)){
						EntryInd = 1;
					}
				}
				if (EntryInd == 1){
					Register[ic].NewStatus = 1;
					// Not nec to call SetNewStatusTo1 since EligibleByAge only gets 
					// calculated later (in the CalcPartnerTransition function).
					Register[ic].FSWind = 1;
					CSWregister[TotCurrFSW] = ic + 1;
					TotCurrFSW += 1;
				}
				temp += EntryProb;
			}
		}
	}
}

void Pop::UpdateMarriageIncidence()
{
	int ia, ic, ig, ii, ij;
	int ID1, ID2; // IDs of primary and secondary partners
	int IR, PR1, PR2; // Risk groups of indiv, primary partner and 2ndary partners
	int UnmarriedPyr[16][4]; // # unmarried individuals, by age & risk group (MH, ML, FH, FL)
	int STpartnersPyr[16][8]; // # ST partnerships that can lead to marriage, by age & risk group 
	// (MH-FH, MH-FL, ML-FH, ML-FL, FH-MH, FH-ML, FL-MH, FL-ML)
	double CurrentSTpartners[2][2]; // Current ST partnerships that can lead to marriage, 
	// by male risk group (1st index) and female risk group (second index)
	double temp;

	// First determine the denominators
	for (ia = 0; ia<16; ia++){
		for (ii = 0; ii<8; ii++){
			STpartnersPyr[ia][ii] = 0;
		}
		for (ii = 0; ii<4; ii++){
			UnmarriedPyr[ia][ii] = 0;
		}
	}
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1){
			if (Register[ic].AgeGroup >= 2 && Register[ic].MarriedInd == 0){
				// We include the condition MarriedInd==0 because we are assuming people
				// cannot get married to a 2nd person if they are already married
				ia = Register[ic].AgeGroup - 2;
				ig = Register[ic].SexInd;
				IR = Register[ic].RiskGroup - 1;
				UnmarriedPyr[ia][ig * 2 + IR] += 1;
				if (Register[ic].CurrPartners>0){
					ID1 = Register[ic].IDprimary;
					if (ID1 == 0){
						cout << "Error: indiv " << ic + 1 << " has partner with 0 ID" << endl;
					}
					PR1 = Register[ID1 - 1].RiskGroup - 1;
					if (Register[ID1 - 1].MarriedInd == 0){
						// Note that this condition wasn't included in the deterministic
						// model because we weren't categorizing individuals according to
						// whether their ST partners are married or not.
						STpartnersPyr[ia][ig * 4 + IR * 2 + PR1] += 1;
					}
				}
				if (Register[ic].CurrPartners == 2){
					ID2 = Register[ic].ID2ndary;
					PR2 = Register[ID2 - 1].RiskGroup - 1;
					if (Register[ID2 - 1].MarriedInd == 0){
						STpartnersPyr[ia][ig * 4 + IR * 2 + PR2] += 1;
					}
				}
			}
		}
	}

	// Then determine the MarriageRate adjustment factors
	for (ii = 0; ii<2; ii++){
		for (ij = 0; ij<2; ij++){
			CurrentSTpartners[ii][ij] = 0;
		}
	}
	for (ia = 0; ia<16; ia++){
		CurrentSTpartners[0][0] += STpartnersPyr[ia][0];
		CurrentSTpartners[0][1] += STpartnersPyr[ia][1];
		CurrentSTpartners[1][0] += STpartnersPyr[ia][2];
		CurrentSTpartners[1][1] += STpartnersPyr[ia][3];
	}
	ActualPropnSTH[0][0] = CurrentSTpartners[0][0] / (CurrentSTpartners[0][0] +
		CurrentSTpartners[0][1]);
	ActualPropnSTH[1][0] = CurrentSTpartners[1][0] / (CurrentSTpartners[1][0] +
		CurrentSTpartners[1][1]);
	ActualPropnSTH[0][1] = CurrentSTpartners[0][0] / (CurrentSTpartners[0][0] +
		CurrentSTpartners[1][0]);
	ActualPropnSTH[1][1] = CurrentSTpartners[0][1] / (CurrentSTpartners[0][1] +
		CurrentSTpartners[1][1]);
	MarriageRate[0][0] = ((1.0 - ActualPropnSTH[0][0]) / ActualPropnSTH[0][0]) /
		((1.0 - ActualPropnLTH[0][0]) / ActualPropnLTH[0][0]);
	MarriageRate[0][1] = ((1.0 - ActualPropnSTH[1][0]) / ActualPropnSTH[1][0]) /
		((1.0 - ActualPropnLTH[1][0]) / ActualPropnLTH[1][0]);
	MarriageRate[0][2] = ((1.0 - ActualPropnSTH[0][1]) / ActualPropnSTH[0][1]) /
		((1.0 - ActualPropnLTH[0][1]) / ActualPropnLTH[0][1]);
	MarriageRate[0][3] = ((1.0 - ActualPropnSTH[1][1]) / ActualPropnSTH[1][1]) /
		((1.0 - ActualPropnLTH[1][1]) / ActualPropnLTH[1][1]);
	for (ij = 0; ij<4; ij++){
		MarriageRate[1][ij] = 1.0;
	}

	// Then calculate AgeEffectMarriage
	for (ia = 0; ia<16; ia++){
		temp = STpartnersPyr[ia][0] * MarriageRate[0][0] + STpartnersPyr[ia][1] * MarriageRate[1][0];
		if (temp>0){ AgeEffectMarriage[ia][0] = MarriageIncidence[ia][0] * UnmarriedPyr[ia][0] / temp; }
		else{ AgeEffectMarriage[ia][0] = 0.0; }
		temp = STpartnersPyr[ia][2] * MarriageRate[0][1] + STpartnersPyr[ia][3] * MarriageRate[1][1];
		if (temp>0){ AgeEffectMarriage[ia][1] = MarriageIncidence[ia][0] * UnmarriedPyr[ia][1] / temp; }
		else{ AgeEffectMarriage[ia][1] = 0.0; }
		temp = STpartnersPyr[ia][4] * MarriageRate[0][2] + STpartnersPyr[ia][5] * MarriageRate[1][2];
		if (temp>0){ AgeEffectMarriage[ia][2] = MarriageIncidence[ia][1] * UnmarriedPyr[ia][2] / temp; }
		else{ AgeEffectMarriage[ia][2] = 0.0; }
		temp = STpartnersPyr[ia][6] * MarriageRate[0][3] + STpartnersPyr[ia][7] * MarriageRate[1][3];
		if (temp>0){ AgeEffectMarriage[ia][3] = MarriageIncidence[ia][1] * UnmarriedPyr[ia][3] / temp; }
		else{ AgeEffectMarriage[ia][3] = 0.0; }
	}
}

void Pop::BalanceSexualPartners()
{
	int ia, ic, ig, ii, ij, ir, is;
	int ID1, ID2, PR1, PR2;
	//int TotFSWbyAge[16];
	//double temp;

	for (ir = 0; ir<2; ir++){
		for (ig = 0; ig<2; ig++){
			DesiredSTpartners[ir][ig] = 0.0;
		}
	}

	// First calculate the desired rates at which new partnerships are formed
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup >= 2){
			if (Register[ic].CurrPartners == 0 || (Register[ic].CurrPartners == 1 && Register[ic].RiskGroup == 1)){
				ia = Register[ic].AgeGroup - 2;
				ig = Register[ic].SexInd;
				ir = Register[ic].RiskGroup - 1;
				is = Register[ic].HIVstage;
				if (Register[ic].VirginInd == 0){
					Register[ic].DesiredNewPartners = PartnershipFormation[ir][ig] *
						AgeEffectPartners[ia][ig] * Register[ic].PartnerRateAdj;
					if (Register[ic].CurrPartners == 1 && Register[ic].MarriedInd == 0){
						Register[ic].DesiredNewPartners *= PartnerEffectNew[0][ig];
					}
					if (Register[ic].CurrPartners == 1 && Register[ic].MarriedInd == 1){
						Register[ic].DesiredNewPartners *= PartnerEffectNew[1][ig];
					}
					if (is>0){
						Register[ic].DesiredNewPartners *= HIVeffectPartners[is - 1];
					}
					if (Register[ic].FSWind == 1){
						Register[ic].DesiredNewPartners = 0.0;
					}
					DesiredSTpartners[ir][ig] += Register[ic].DesiredNewPartners; // Excludes virgins
				}
				else{
					if (ia<4){
						Register[ic].DesiredNewPartners = SexualDebut[ia][ig];
					}
					else{
						Register[ic].DesiredNewPartners = SexualDebut[3][ig];
					}
					if (ir == 1){
						Register[ic].DesiredNewPartners *= DebutAdjLow[ig];
					}
				}
			}
			else{
				Register[ic].DesiredNewPartners = 0.0;
			}
		}
	}

	// Then calculate the balancing factors for the rates of ST partner acquisition
	DesiredPartnerRiskM[0][0] = (1.0 - AssortativeM) + AssortativeM *
		DesiredSTpartners[0][1] / (DesiredSTpartners[0][1] + DesiredSTpartners[1][1]);
	DesiredPartnerRiskM[0][1] = 1.0 - DesiredPartnerRiskM[0][0];
	DesiredPartnerRiskM[1][1] = (1.0 - AssortativeM) + AssortativeM *
		DesiredSTpartners[1][1] / (DesiredSTpartners[0][1] + DesiredSTpartners[1][1]);
	DesiredPartnerRiskM[1][0] = 1.0 - DesiredPartnerRiskM[1][1];

	DesiredPartnerRiskF[0][0] = (1.0 - AssortativeF) + AssortativeF *
		DesiredSTpartners[0][0] / (DesiredSTpartners[0][0] + DesiredSTpartners[1][0]);
	DesiredPartnerRiskF[0][1] = 1.0 - DesiredPartnerRiskF[0][0];
	DesiredPartnerRiskF[1][1] = (1.0 - AssortativeF) + AssortativeF *
		DesiredSTpartners[1][0] / (DesiredSTpartners[0][0] + DesiredSTpartners[1][0]);
	DesiredPartnerRiskF[1][0] = 1.0 - DesiredPartnerRiskF[1][1];

	if (AllowBalancing == 1){
		for (ii = 0; ii<2; ii++){
			for (ij = 0; ij<2; ij++){
				AdjSTrateM[ii][ij] = (GenderEquality * DesiredSTpartners[ij][1] *
					DesiredPartnerRiskF[ij][ii] + (1.0 - GenderEquality) *
					DesiredSTpartners[ii][0] * DesiredPartnerRiskM[ii][ij]) /
					(DesiredSTpartners[ii][0] * DesiredPartnerRiskM[ii][ij]);
				AdjSTrateF[ii][ij] = (GenderEquality * DesiredSTpartners[ii][1] *
					DesiredPartnerRiskF[ii][ij] + (1.0 - GenderEquality) *
					DesiredSTpartners[ij][0] * DesiredPartnerRiskM[ij][ii]) /
					(DesiredSTpartners[ii][1] * DesiredPartnerRiskF[ii][ij]);
			}
		}
	}
	else{
		for (ii = 0; ii<2; ii++){
			for (ij = 0; ij<2; ij++){
				AdjSTrateM[ii][ij] = 1.0;
				AdjSTrateF[ii][ij] = 1.0;
			}
		}
	}

	// Secondly calculate the desired rates at which marriages are formed
	for (ii = 0; ii<2; ii++){
		for (ij = 0; ij<2; ij++){
			DesiredMarriagesM[ii][ij] = 0.0;
			DesiredMarriagesF[ii][ij] = 0.0;
		}
	}
	
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup >= 2){
			if (Register[ic].CurrPartners>0 && Register[ic].MarriedInd == 0){
				ia = Register[ic].AgeGroup - 2;
				ig = Register[ic].SexInd;
				ir = Register[ic].RiskGroup - 1;
				ID1 = Register[ic].IDprimary;
				if (Register[ID1 - 1].MarriedInd == 0){ // This condition isn't included in THISA
					PR1 = Register[ID1 - 1].RiskGroup - 1;
					if (ig == 0){
						DesiredMarriagesM[ir][PR1] += MarriageRate[PR1][ir] *
							AgeEffectMarriage[ia][ir];
					}
					else{
						DesiredMarriagesF[ir][PR1] += MarriageRate[PR1][2 + ir] *
							AgeEffectMarriage[ia][2 + ir];
					}
				}
				if (Register[ic].CurrPartners == 2){
					ID2 = Register[ic].ID2ndary;
					if (Register[ID2 - 1].MarriedInd == 0){
						PR2 = Register[ID2 - 1].RiskGroup - 1;
						if (ig == 0){
							DesiredMarriagesM[ir][PR2] += MarriageRate[PR2][ir] *
								AgeEffectMarriage[ia][ir];
						}
						else{
							DesiredMarriagesF[ir][PR2] += MarriageRate[PR2][2 + ir] *
								AgeEffectMarriage[ia][2 + ir];
						}
					}
				}
			}
		}
	}

	// Then calculate the balancing factors for the rates of marriage
	if (AllowBalancing == 1){
		for (ii = 0; ii<2; ii++){
			for (ij = 0; ij<2; ij++){
				AdjLTrateM[ii][ij] = (GenderEquality * DesiredMarriagesF[ij][ii] +
					(1.0 - GenderEquality) * DesiredMarriagesM[ii][ij]) / DesiredMarriagesM[ii][ij];
				AdjLTrateF[ii][ij] = (GenderEquality * DesiredMarriagesF[ii][ij] +
					(1.0 - GenderEquality) * DesiredMarriagesM[ij][ii]) / DesiredMarriagesF[ii][ij];
			}
		}
	}
	else{
		for (ii = 0; ii<2; ii++){
			for (ij = 0; ij<2; ij++){
				AdjLTrateM[ii][ij] = 1.0;
				AdjLTrateF[ii][ij] = 1.0;
			}
		}
	}
}

void Pop::CalcPartnerTransitions()
{
	int ia, ic, ii, ir;
	int TotUnassigned, nextID, EventType, pID, pID2, prisk, index1;
	double  rr; //TotalSelectionProb,
	vector<int> indices(Register.size());

	//int seed1 = 4727 + CurrYear * 29 + BehavCycleCount * 51;
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed1 += CurrSim;}
	//CRandomMersenne rg1(seed1);
	for (ic = 0; ic<MaxPop; ic++){
		r2[ic] = rg.Random();
	}
	//int seed2 = 8309  + CurrYear * 57 + BehavCycleCount * 31; 
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed2 += CurrSim;}
	//CRandomMersenne rg2(seed2);
	for (ic = 0; ic<MaxPop; ic++){
		revent[ic] = rg.Random();
	}
	// New addition for sampling partner age;
	//int seed3 = 927  + CurrYear * 37 + BehavCycleCount * 92;
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed3 += CurrSim;}
	//CRandomMersenne rg3(seed3);
	for (ic = 0; ic<MaxPop; ic++){
		rpAge[ic] = rg.Random();
	}
	//int seed4 = 1598  + CurrYear * 17 + BehavCycleCount * 42;
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed4 += CurrSim;}
	//CRandomMersenne rg4(seed4);
	for (ic = 0; ic<MaxPop; ic++){
		rpID2[ic] = rg.Random();
	}

	TotUnassigned = 0; // Number who haven't yet been assigned a status
	for (ia = 0; ia<16; ia++){
		//EligibleByAgeM[ia][0] = 0.0;
		//EligibleByAgeF[ia][0] = 0.0;
		//EligibleByAgeM[ia][1] = 0.0;
		//EligibleByAgeF[ia][1] = 0.0;
		MalePartners[ia][0].Reset();
		MalePartners[ia][1].Reset();
		FemPartners[ia][0].Reset();
		FemPartners[ia][1].Reset();
	}
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup >= 2 && Register[ic].NewStatus == 0){
			TotUnassigned += 1;
			if (Register[ic].DesiredNewPartners>0 && Register[ic].VirginInd == 0){
				ia = Register[ic].AgeGroup - 2; // -2 because age groups start at 10
				ir = Register[ic].RiskGroup - 1;
				if (Register[ic].SexInd == 0){
					MalePartners[ia][ir].AddMember(ic + 1);
				}
				else{
					FemPartners[ia][ir].AddMember(ic + 1);
				}
			}
		}
	}
	int ipp = indices.size();
	for (ii = 0; ii<ipp; ii++){
		indices[ii] = ii;
	}
	ii = 0;
	while (TotUnassigned>0){
		nextID = 0;
		while (nextID == 0){
			index1 = static_cast<int>(r2[ii] * indices.size());
			nextID = indices[index1] + 1;
			if (Register[nextID - 1].AliveInd == 0 || Register[nextID - 1].AgeGroup<2 ||
				Register[nextID - 1].NewStatus == 1){
				nextID = 0;
				ii += 1;
			}
			indices[index1] = indices[indices.size() - 1];
			indices.pop_back();
		}
		rr = revent[ii];
		EventType = Register[nextID - 1].SelectEvent(rr);
		SetNewStatusTo1(nextID);
		TotUnassigned = TotUnassigned - 1;
		if (EventType == 0){ // No change in behav status
			if (Register[nextID - 1].CurrPartners>0){
				pID = Register[nextID - 1].IDprimary;
				prisk = Register[pID - 1].RiskGroup;
				if (prisk == 2 && Register[pID - 1].NewStatus == 0){
					SetNewStatusTo1(pID);
					TotUnassigned = TotUnassigned - 1;
				}
			}
			if (Register[nextID - 1].CurrPartners == 2){
				pID = Register[nextID - 1].ID2ndary;
				prisk = Register[pID - 1].RiskGroup;
				if (prisk == 2 && Register[pID - 1].NewStatus == 0){
					SetNewStatusTo1(pID);
					TotUnassigned = TotUnassigned - 1;
				}
			}
		}
		if (EventType == 1){ // Acquire ST partner in high risk group
			GetNewPartner(nextID, rpAge[ii], rpID2[ii], 1);
			if (MaxNewPartnerInd == 1){
				GetNewPartner(nextID, rpAge[ii], rpID2[ii], 2);
			}
			if (MaxNewPartnerInd != 1){ // Indiv has been assigned a partner
				TotUnassigned = TotUnassigned - 1;
			}
			if (MaxNewPartnerInd == 1){ // Will only happen if there are no 
				// potential partners in either high or low risk group. In this
				// case, assume there is no change in behavioural status.
				if (Register[nextID - 1].CurrPartners>0){
					pID = Register[nextID - 1].IDprimary;
					prisk = Register[pID - 1].RiskGroup;
					if (prisk == 2 && Register[pID - 1].NewStatus == 0){
						SetNewStatusTo1(pID);
						TotUnassigned = TotUnassigned - 1;
					}
				}
			}
		}
		if (EventType == 2){ // Acquire ST partner in low risk group
			GetNewPartner(nextID, rpAge[ii], rpID2[ii], 2);
			if (MaxNewPartnerInd == 1){
				GetNewPartner(nextID, rpAge[ii], rpID2[ii], 1);
			}
			if (MaxNewPartnerInd != 1){ // Indiv has been assigned a partner
				TotUnassigned = TotUnassigned - 1;
			}
			if (MaxNewPartnerInd == 1){ // Will only happen if there are no 
				// potential partners in either high or low risk group. In this
				// case, assume there is no change in behavioural status.
				if (Register[nextID - 1].CurrPartners>0){
					pID = Register[nextID - 1].IDprimary;
					prisk = Register[pID - 1].RiskGroup;
					if (prisk == 2 && Register[pID - 1].NewStatus == 0){
						SetNewStatusTo1(pID);
						TotUnassigned = TotUnassigned - 1;
					}
				}
			}
		}
		if (EventType == 3){ // Marry primary partner
			pID = Register[nextID - 1].IDprimary;
			Register[nextID - 1].MarriedInd = 1;
			Register[pID - 1].MarriedInd = 1;
			SetNewStatusTo1(pID);
			TotUnassigned = TotUnassigned - 1;
			if (Register[pID - 1].ID2ndary == nextID){
				Register[pID - 1].ID2ndary = Register[pID - 1].IDprimary;
				Register[pID - 1].IDprimary = nextID;
			}
			if (Register[pID - 1].CurrPartners == 2){
				pID2 = Register[pID - 1].ID2ndary;
				prisk = Register[pID2 - 1].RiskGroup;
				if (prisk == 2 && Register[pID2 - 1].NewStatus == 0){
					SetNewStatusTo1(pID2);
					TotUnassigned = TotUnassigned - 1;
				}
			}
			if (Register[nextID - 1].CurrPartners == 2){
				pID = Register[nextID - 1].ID2ndary;
				prisk = Register[pID - 1].RiskGroup;
				if (prisk == 2 && Register[pID - 1].NewStatus == 0){
					SetNewStatusTo1(pID);
					TotUnassigned = TotUnassigned - 1;
				}
			}
		}
		if (EventType == 4){ // Marry 2ndary partner
			pID = Register[nextID - 1].ID2ndary;
			Register[nextID - 1].MarriedInd = 1;
			Register[pID - 1].MarriedInd = 1;
			SetNewStatusTo1(pID);
			TotUnassigned = TotUnassigned - 1;
			Register[nextID - 1].ID2ndary = Register[nextID - 1].IDprimary;
			Register[nextID - 1].IDprimary = pID;
			if (Register[pID - 1].ID2ndary == nextID){
				Register[pID - 1].ID2ndary = Register[pID - 1].IDprimary;
				Register[pID - 1].IDprimary = nextID;
			}
			if (Register[pID - 1].CurrPartners == 2){
				pID2 = Register[pID - 1].ID2ndary;
				prisk = Register[pID2 - 1].RiskGroup;
				if (prisk == 2 && Register[pID2 - 1].NewStatus == 0){
					SetNewStatusTo1(pID2);
					TotUnassigned = TotUnassigned - 1;
				}
			}
			pID = Register[nextID - 1].ID2ndary; // Not the same as before
			prisk = Register[pID - 1].RiskGroup;
			if (prisk == 2 && Register[pID - 1].NewStatus == 0){
				SetNewStatusTo1(pID);
				TotUnassigned = TotUnassigned - 1;
			}
		}
		if (EventType == 5){ // Divorce
			pID = Register[nextID - 1].IDprimary;
			Register[nextID - 1].MarriedInd = 0;
			Register[pID - 1].MarriedInd = 0;
			SetNewStatusTo1(pID);
			TotUnassigned = TotUnassigned - 1;
			if (Register[nextID - 1].CurrPartners == 1){
				Register[nextID - 1].CurrPartners = 0;
				Register[nextID - 1].IDprimary = 0; 
			}
			else{
				Register[nextID - 1].CurrPartners = 1;
				Register[nextID - 1].IDprimary = Register[nextID - 1].ID2ndary;
				Register[nextID - 1].ID2ndary = 0;
			}
			if (Register[pID - 1].CurrPartners == 1){
				Register[pID - 1].CurrPartners = 0;
				Register[pID - 1].IDprimary = 0; 
			}
			else{
				Register[pID - 1].CurrPartners = 1;
				Register[pID - 1].IDprimary = Register[pID - 1].ID2ndary;
				Register[pID - 1].ID2ndary = 0;
			}
		}
		if (EventType == 6){ // End primary ST relationship
			pID = Register[nextID - 1].IDprimary;
			if (Register[nextID - 1].CurrPartners == 1){
				Register[nextID - 1].IDprimary = 0; 
				Register[nextID - 1].CurrPartners = 0;
			}
			else{
				Register[nextID - 1].IDprimary = Register[nextID - 1].ID2ndary;
				Register[nextID - 1].ID2ndary = 0;
				Register[nextID - 1].CurrPartners = 1;
			}
			if (Register[pID - 1].CurrPartners == 1){
				Register[pID - 1].IDprimary = 0; 
				Register[pID - 1].CurrPartners = 0;
			}
			else{
				if (Register[pID - 1].IDprimary == nextID){
					Register[pID - 1].IDprimary = Register[pID - 1].ID2ndary;
				}
				Register[pID - 1].ID2ndary = 0;
				Register[pID - 1].CurrPartners = 1;
			}
			SetNewStatusTo1(pID);
			TotUnassigned = TotUnassigned - 1;
		}
		if (EventType == 7){ // End secondary ST relationship
			pID = Register[nextID - 1].ID2ndary;
			Register[nextID - 1].ID2ndary = 0;
			Register[nextID - 1].CurrPartners = 1;
			if (Register[pID - 1].CurrPartners == 1){
				Register[pID - 1].IDprimary = 0; 
				Register[pID - 1].CurrPartners = 0;
			}
			else{
				if (Register[pID - 1].IDprimary == nextID){
					Register[pID - 1].IDprimary = Register[pID - 1].ID2ndary;
				}
				Register[pID - 1].ID2ndary = 0;
				Register[pID - 1].CurrPartners = 1;
			}
			SetNewStatusTo1(pID);
			TotUnassigned = TotUnassigned - 1;
		}
		/*TotUnassigned = 0;
		for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1 && Register[ic].AgeGroup>=2 && Register[ic].NewStatus==0){
		TotUnassigned += 1;}
		}*/
		ii += 1;
	}
}

void Pop::GetNewPartner(int ID, double rnd, double rnd2, int rsk)
{
	// Similar to the ChooseSTpartner function.

	int ic, ia, ib;
	int iage, page; // individual's age group - 2 (so that 0 ==> age 10-14), partner's age group - 2
	int igender, pgender; // individual's gender and partner's gender
	int irisk; // individual's risk group -1 (so that 0 ==> high risk)
	//double EligibleByAge[16];
	double rnd3; // , WeightsByAge[16];
	//double Normalizer, TotalSelectionProb2;

	iage = Register[ID - 1].AgeGroup - 2;
	igender = Register[ID - 1].SexInd;
	pgender = 1 - igender;
	irisk = Register[ID - 1].RiskGroup - 1;

	MaxNewPartnerInd = 0;

	// (1) Select partner age
	if (igender == 0){
		for (ia = 0; ia<16; ia++){
			//if(rnd<CumAgePrefM[iage][ia] && FemPartners[ia][rsk-1].TotalDesire>0.0){
			if (rnd<CumAgePrefM[iage][ia] && FemPartners[ia][rsk - 1].Pool.size()>0){
				page = ia;
				break;
			}
			//if(rnd<CumAgePrefM[iage][ia] && FemPartners[ia][rsk-1].TotalDesire==0.0){
			if (rnd<CumAgePrefM[iage][ia] && FemPartners[ia][rsk - 1].Pool.size() == 0){
				// Sample another age group
				if (ia == 0){
					rnd3 = rnd / AgePrefM[iage][0];
				}
				else{
					rnd3 = (rnd - CumAgePrefM[iage][ia - 1]) / AgePrefM[iage][ia];
				}
				for (ib = 0; ib<16; ib++){
					if (rnd3<CumAgePrefM[iage][ib]){
						//if(FemPartners[ib][rsk-1].TotalDesire>0.0){
						if (FemPartners[ib][rsk - 1].Pool.size()>0){
							page = ib;
						}
						else{
							// No luck after 2 tries, so give up.
							MaxNewPartnerInd = 1;
						}
						break;
					}
				}
				break;
			}
		}
	}
	else{
		for (ia = 0; ia<16; ia++){
			//if(rnd<CumAgePrefF[iage][ia] && MalePartners[ia][rsk-1].TotalDesire>0.0){
			if (rnd<CumAgePrefF[iage][ia] && MalePartners[ia][rsk - 1].Pool.size()>0){
				page = ia;
				break;
			}
			//if(rnd<CumAgePrefF[iage][ia] && MalePartners[ia][rsk-1].TotalDesire==0.0){
			if (rnd<CumAgePrefF[iage][ia] && MalePartners[ia][rsk - 1].Pool.size() == 0){
				// Sample another age group
				if (ia == 0){
					rnd3 = rnd / AgePrefF[iage][0];
				}
				else{
					rnd3 = (rnd - CumAgePrefF[iage][ia - 1]) / AgePrefF[iage][ia];
				}
				for (ib = 0; ib<16; ib++){
					if (rnd3<CumAgePrefF[iage][ib]){
						//if(MalePartners[ib][rsk-1].TotalDesire>0.0){
						if (MalePartners[ib][rsk - 1].Pool.size()>0){
							page = ib;
						}
						else{
							// No luck after 2 tries, so give up.
							MaxNewPartnerInd = 1;
						}
						break;
					}
				}
				break;
			}
		}
	}

	// (2) Select partner ID
	if (MaxNewPartnerInd == 0){
		if (igender == 0){
			ic = FemPartners[page][rsk - 1].SamplePool(rnd2);
		}
		else{
			ic = MalePartners[page][rsk - 1].SamplePool(rnd2);
		}
		// Because of rounding errors, there's a small chance that an individual
		// may be selected even if they're not eligible - so check eligibility.
		if (Register[ic - 1].NewStatus == 1){
			MaxNewPartnerInd = 1;
		}
		else{
			if (Register[ID - 1].IDprimary == 0){
				Register[ID - 1].IDprimary = ic;
			}
			else{
				Register[ID - 1].ID2ndary = ic;
			}
			if (Register[ic - 1].IDprimary == 0){
				Register[ic - 1].IDprimary = ID;
			}
			else{
				Register[ic - 1].ID2ndary = ID;
			}
			Register[ID - 1].CurrPartners += 1;
			Register[ic - 1].CurrPartners += 1;
			Register[ID - 1].LifetimePartners += 1;
			Register[ic - 1].LifetimePartners += 1;
			if (Register[ID - 1].VirginInd == 1){
				Register[ID - 1].VirginInd = 0;
				Register[ID - 1].SexDebutAge = Register[ID - 1].AgeGroup;
			}
			if (Register[ic - 1].VirginInd == 1){
				Register[ic - 1].VirginInd = 0;
				Register[ic - 1].SexDebutAge = Register[ic - 1].AgeGroup;
			}
			SetNewStatusTo1(ic);
		}
	}
}

void Pop::SetNewStatusTo1(int ID)
{
	int ia, ir;
	double Decrement;

	Register[ID - 1].NewStatus = 1;
	Decrement = Register[ID - 1].DesiredNewPartners;
	if (Decrement>0.0 && Register[ID - 1].VirginInd == 0){
		ia = Register[ID - 1].AgeGroup - 2;
		ir = Register[ID - 1].RiskGroup - 1;
		if (Register[ID - 1].SexInd == 0){
			MalePartners[ia][ir].TotalDesire =
				MalePartners[ia][ir].TotalDesire - Decrement;
		}
		else{
			FemPartners[ia][ir].TotalDesire =
				FemPartners[ia][ir].TotalDesire - Decrement;
		}
	}
}

void Pop::SetToDead(int ID)
{
	// Differs from SetNewStatusTo 1 in that (a) we have to adjust partner profile(s)
	// but (b) we don't have to correct EligibleByAge (since it hasn't yet been calculated).

	int ic, ij, pID1, pID2;

	Register[ID - 1].AliveInd = 0;
	Register[ID - 1].NewStatus = 1;
	if (Register[ID - 1].IDprimary>0){
		pID1 = Register[ID - 1].IDprimary;
		if (Register[pID1 - 1].CurrPartners == 1){
			Register[pID1 - 1].IDprimary = 0; 
			Register[pID1 - 1].CurrPartners = 0;
		}
		else{
			if (Register[pID1 - 1].IDprimary == ID){
				Register[pID1 - 1].IDprimary = Register[pID1 - 1].ID2ndary;
			}
			Register[pID1 - 1].ID2ndary = 0;
			Register[pID1 - 1].CurrPartners = 1;
		}
		Register[pID1 - 1].NewStatus = 1;
		if (Register[ID - 1].MarriedInd == 1){
			Register[pID1 - 1].MarriedInd = 0;
		}
	}
	if (Register[ID - 1].ID2ndary>0){
		pID2 = Register[ID - 1].ID2ndary;
		if (Register[pID2 - 1].CurrPartners == 1){
			Register[pID2 - 1].IDprimary = 0; 
			Register[pID2 - 1].CurrPartners = 0;
		}
		else{
			if (Register[pID2 - 1].IDprimary == ID){
				Register[pID2 - 1].IDprimary = Register[pID2 - 1].ID2ndary;
			}
			Register[pID2 - 1].ID2ndary = 0;
			Register[pID2 - 1].CurrPartners = 1;
		}
		Register[pID2 - 1].NewStatus = 1;
	}
	if (Register[ID - 1].FSWind == 1){
		TotCurrFSW = TotCurrFSW - 1;
		// Update CSWregister
		for (ic = 0; ic<MaxCSWs; ic++){
			if (CSWregister[ic] == ID){
				for (ij = ic; ij<MaxCSWs - 1; ij++){
					CSWregister[ij] = CSWregister[ij + 1];
				}
				CSWregister[MaxCSWs - 1] = 0;
			}
		}
	}
	Register[ID - 1].CurrPartners = 0;
	Register[ID - 1].IDprimary = 0;
	Register[ID - 1].ID2ndary = 0;
}

void Pop::NewBirth(int ID)
{
	int ic, ind, MatAgeGrp; //popsize, 
	Indiv KidA;
	double r[30];
	double x, a, b, p, q; // y,

	memset(&KidA, 0, sizeof(KidA));
	//int seed = 4475 + CurrYear * 62 + BehavCycleCount * 98 + ID * 23;
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for (ic = 0; ic<30; ic++){
		r[ic] = rg.Random();
	}

	KidA.AliveInd = 1;
	KidA.VirginInd = 1;
	KidA.MarriedInd = 0;
	KidA.CurrPartners = 0;
	KidA.FSWind = 0;
	KidA.IDprimary = 0;
	KidA.ID2ndary = 0;
	KidA.NonHIVfertRate = 0.0;
	//HPV
	KidA.ARTweeks = 0;
	KidA.ARTstage = 0;
	KidA.InScreen = 0;
	KidA.ScreenCount = 0;
	KidA.InWHOScreen = 0;
	KidA.Scr30 = 0;
	KidA.Scr50 = 0;
	KidA.Scr35 = 0;
	KidA.Scr45 = 0;
	KidA.Scr16 = 0;
	KidA.Scr19 = 0;
	KidA.Scr22 = 0;
	KidA.Scr25 = 0;
	KidA.Scr28 = 0;
	KidA.Scr31 = 0;
	KidA.Scr34 = 0;
	KidA.Scr37 = 0;
	KidA.Scr40 = 0;
	KidA.Scr43 = 0;
	KidA.Scr46 = 0;
	KidA.Scr49 = 0;
	KidA.Scr52=0;
	KidA.Scr55=0;
	KidA.Scr58=0;
	KidA.Scr61=0;
	KidA.Scr64=0;
	KidA.ScreenResult = 0;
	KidA.timePassed =0;
	KidA.reason=0;
	KidA.repeat=0;
	KidA.HPVrepeat=0;
	KidA.DiagnosedCC=0;
	KidA.TrueStage=0;
	KidA.HPVstatus=0;
	KidA.Age50=0;
	KidA.GotTxV=0;

	// Assign sex, risk group, NonHIVmortProb and initial STD states
	if (r[0]<MaleBirthPropn){
		KidA.SexInd = 0;
		if (r[1]<HighPropnM){ KidA.RiskGroup = 1; }
		else{ KidA.RiskGroup = 2; }
		KidA.NonHIVmortProb = 1.0 - pow(1.0 -
			InfantMort1st6mM[CurrYear - StartYear], 2.0 / CycleS);
		KidA.BVstage = 0;
	}
	else{
		KidA.SexInd = 1;
		if (r[1]<HighPropnF){ KidA.RiskGroup = 1; }
		else{ KidA.RiskGroup = 2; }
		KidA.NonHIVmortProb = 1.0 - pow(1.0 -
			InfantMort1st6mF[CurrYear - StartYear], 2.0 / CycleS);
		KidA.BVstage = 1;
	}
	KidA.CTstage = 0;
	KidA.HDstage = 0;
	KidA.HSVstage = 0;
	KidA.NGstage = 0;
	KidA.TPstage = 0;
	KidA.TVstage = 0;
	KidA.VCstage = 0;
	for (int xx = 0; xx < 13; xx++) { 
		KidA.HPVstage[xx] = 0;
		KidA.HPVstageE[xx] = 0; 
		KidA.TimeinCIN3[xx] = 0; 
		KidA.WeibullCIN3[xx] = 0; 
		KidA.VaccinationStatus[xx]=0;
	}
	KidA.TimeinStageI=0;
	KidA.TimeinStageII=0;
	KidA.TimeinStageIII=0;
	KidA.TimeinStageIV=0;
	KidA.StageIdeath=0;
	KidA.StageIIdeath=0;
	KidA.StageIIIdeath=0;
	KidA.StageIVdeath=0;
	KidA.StageIrecover=0;
	KidA.StageIIrecover=0;
	KidA.StageIIIrecover=0;
	KidA.StageIVrecover=0;
	KidA.ThermalORPap=r[22];
	
	// Assign exact date of birth and age group
	KidA.DOB = 0.5 + CurrYear + (BehavCycleCount - r[2]) / CycleS;
	KidA.AgeGroup = 0;

	// Handle situation in which mother is HIV-positive
	MatAgeGrp = Register[ID - 1].AgeGroup - 3;
	if (Register[ID - 1].HIVstage>0){
		if (r[3]<CurrPerinatal){
			KidA.HIVstage = 2;
			KidA.DateInfect = KidA.DOB;
		}
		else if (r[3]<CurrPerinatal + CurrPostnatal){
			KidA.HIVstage = 2;
			// We arbitrarily assume age of postnatal infection is 6 months.
			// At this stage this is only used for determining if the transmission
			// was perinatal or postnatal, so the exact value is not important.
			KidA.DateInfect = KidA.DOB + 0.5;
		}
		else{
			KidA.HIVstage = 0;
		}
	}

	// Assign PartnerRateAdj
	if (AllowPartnerRateAdj == 1){
		ind = 2;
		// Note that the following formulas for a and b apply only when the gamma mean is 1.
		a = 1.0 / pow(SDpartnerRateAdj, 2.0);
		b = a;
		p = r[4];
		q = 1 - r[4];
		cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
		KidA.PartnerRateAdj = x;
	}
	else{
		KidA.PartnerRateAdj = 1.0;
	}

	// Assign SuscepHIVadj
	if (AllowHIVsuscepAdj == 1){
		ind = 2;
		// Note that the following formulas for a and b apply only when the gamma mean is 1.
		a = 1.0 / pow(SDsuscepHIVadj, 2.0);
		b = a;
		p = r[5];
		q = 1 - r[5];
		cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
		KidA.SuscepHIVadj = x;
	}
	else{
		KidA.SuscepHIVadj = 1.0;
	}
	//Decide if KidA will be vaccinated against HPV (girls born after 2004 will be 9 starting in 2014 and eligible for vaccination)
	if(HPVvacc==1 && WHOvacc==0 && CurrYear>2004){
		//KidA.GotVaccOffer=1;
		if( r[6]<PropVaccinated[CurrYear-StartYear] && KidA.SexInd==1){
			KidA.GotVacc=1;
			for(ic=7 ; ic<20; ic++){
				if(r[ic]<VaccEfficacy[ic-7]){KidA.VaccinationStatus[ic-7]=1;}
			}	
		}
		else if (r[6]<PropVaccinated[CurrYear-StartYear] && KidA.SexInd==0 && BOYSvacc==1 && CurrYear>=(ImplementYR-10)){
			KidA.GotVacc=1;
			for(ic=7 ; ic<20; ic++){
				if(r[ic]<VaccEfficacy[ic-7]){KidA.VaccinationStatus[ic-7]=1;}	
			}	
		}
	}
	if(CurrYear>2004 && HPVvacc==1 && WHOvacc==1){
		//KidA.GotVaccOffer=1;
		if(WHOscenario==0) {
			if( CurrYear<(ImplementYR-9) &&  r[6]<PropVaccinated[CurrYear-StartYear] && KidA.SexInd==1){
				KidA.GotVacc=1;
				for(ic=7 ; ic<20; ic++){
					if(r[ic]<VaccEfficacy[ic-7]){KidA.VaccinationStatus[ic-7]=1;}
				}
			}
			if(CurrYear>=(ImplementYR-9)  &&  r[6]<PropVaccinated[CurrYear-StartYear] && KidA.SexInd==1){
				KidA.GotVacc=1;
				for(ic=7 ; ic<20; ic++){
					if(r[ic]<VaccEfficacyNONA[ic-7]){KidA.VaccinationStatus[ic-7]=1;}
				}	
			}
			if(CurrYear>=(ImplementYR-9)  &&  r[6]<PropVaccinated[CurrYear-StartYear] && KidA.SexInd==0 && BOYSvacc==1){
				KidA.GotVacc=1;
				for(ic=7 ; ic<20; ic++){
					if(r[ic]<VaccEfficacyNONA[ic-7]){KidA.VaccinationStatus[ic-7]=1;}
				}	
			}	
		}
		else if(CurrYear>=(ImplementYR-9)  &&  r[6]<PropVaccinatedWHO	 && KidA.SexInd==1){
			KidA.GotVacc=1;
			for(ic=7 ; ic<20; ic++){
				if(r[ic]<VaccEfficacyNONA[ic-7]){KidA.VaccinationStatus[ic-7]=1;}
			}	
		}
	}
	if(KidA.GotVacc==1 && VaccineWane==1){
		KidA.TimeVacc = 0;
		KidA.ExpVacc = 48*20 + (-48 * VaccDur * log(r[22]));
		if(KidA.ExpVacc==0){KidA.ExpVacc=1;}
	} 
	

	Register.push_back(KidA);
}

void Pop::OneSTDcycle()
{
	int ic, xx;

	STDcycleCount += 1;

	GetSexActs();
		
	if (HIVind == 1){ GetHIVtransitions(); }
	if (HSVind == 1 || TPind == 1 || HDind == 1 || NGind == 1 || CTind == 1 ||
		TVind == 1 || BVind == 1 || VCind == 1 || HPVind == 1){
		GetSTDtransitions();
	}
	
	// Update disease states
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		
		if (Register[ic].AliveInd == 1){
			if (HIVind == 1){
				Register[ic].HIVstage = Register[ic].HIVstageE;
			}
			if (HSVind == 1){
				Register[ic].HSVstage = Register[ic].HSVstageE;
			}
			if (TPind == 1){
				Register[ic].TPstage = Register[ic].TPstageE;
			}
			if (HDind == 1){
				Register[ic].HDstage = Register[ic].HDstageE;
			}
			if (NGind == 1){
				Register[ic].NGstage = Register[ic].NGstageE;
			}
			if (CTind == 1){
				Register[ic].CTstage = Register[ic].CTstageE;
			}
			if (TVind == 1){
				Register[ic].TVstage = Register[ic].TVstageE;
			}
			if (Register[ic].SexInd == 1){
				if (BVind == 1 && Register[ic].AgeGroup >= 2){
					Register[ic].BVstage = Register[ic].BVstageE;
				}
				if (VCind == 1){
					Register[ic].VCstage = Register[ic].VCstageE;
				}
			}
			if (HPVind == 1){
				//if (OneType == 0 && targets == 0){
					for (xx = 0; xx < 13; xx++){
						Register[ic].HPVstage[xx] = Register[ic].HPVstageE[xx];
					}
				//}
				//else{ Register[ic].HPVstage[WhichType] = Register[ic].HPVstageE[WhichType]; }
			}
		}
	}	
}

void Pop::GetSexActs()
{
	int ic, CSWID;
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].VirginInd == 0){
			if (Register[ic].SexInd == 0){
				Register[ic].SimulateSexActs(ic + 1);
				if (Register[ic].UVICSW + Register[ic].PVICSW > 0){
					CSWID = Register[ic].IDofCSW;
					Register[ic].IDofCSW = CSWregister[CSWID];
				}
			}
			else{
				if (Register[ic].CurrPartners == 0){
					Register[ic].UVIprimary = 0;
					Register[ic].PVIprimary = 0;
				}
				if (Register[ic].CurrPartners<2){
					Register[ic].UVI2ndary = 0;
					Register[ic].PVI2ndary = 0;
				}
			}
		}
	}
}

void Pop::GetHIVtransitions()
{
	int ic;
	int  seedy;
	seedy = CurrSim * 92 + process_num * 7928 + CurrYear + 25;
	CRandomMersenne rg1(seedy);
	//int seed = 7819 + CurrYear * 35 + BehavCycleCount * 73 + STDcycleCount * 22;
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		r2[ic] = rg.Random();
		if(CurrYear>=ImplementYR){ 
			hiv1618[ic] = rg1.Random();
			hivoth[ic] = rg1.Random();
			wane[ic] = rg1.Random();
			RandAcceptTxV[ic] = rg1.Random();
			Rloss[ic] = rg1.Random();
			EfficacyD1Rand[ic] = rg1.Random();
			EfficacyD2Rand[ic] = rg1.Random();
		}
	}

	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].HIVstage>0 &&
			Register[ic].AgeGroup<3 && (Register[ic].DateInfect -
			Register[ic].DOB) < 3.0){
			// Vertically-infected children 
			if (Register[ic].SexInd == 0 &&
				Register[ic].DOB == Register[ic].DateInfect){
				MaleChild.Perinatal.GetNewStage(ic + 1, r2[ic]);
			}
			else if (Register[ic].SexInd == 0){
				MaleChild.Breastmilk.GetNewStage(ic + 1, r2[ic]);
			}
			if (Register[ic].SexInd == 1 &&
				Register[ic].DOB == Register[ic].DateInfect){
				FemChild.Perinatal.GetNewStage(ic + 1, r2[ic]);
			}
			else if (Register[ic].SexInd == 1){
				FemChild.Breastmilk.GetNewStage(ic + 1, r2[ic]);
			}
		}
		else if (Register[ic].AliveInd == 1 && (Register[ic].HIVstage>0 ||
			Register[ic].VirginInd == 0)){
			Register[ic].GetNewHIVstate(ic + 1, r2[ic], hiv1618[ic], hivoth[ic], wane[ic], RandAcceptTxV[ic], Rloss[ic], EfficacyD1Rand[ic], EfficacyD2Rand[ic]);
		}
	}
}

void Pop::GetSTDtransitions()
{
	int ic, id, xx;
	int seed, SimCount2;
	//int seed = 3705 + CurrYear * 37 + BehavCycleCount * 81 + STDcycleCount * 53;
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){
		for (id = 0; id<59; id++){
			rSTI[ic][id] = rg.Random();
		}
	}

	for (ic = 0; ic<tpp; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup >= 2){
			if(VaccineWane==1 && Register[ic].GotVacc==1){ 
				if(Register[ic].ExpVacc>0 && Register[ic].TimeVacc<Register[ic].ExpVacc){
					Register[ic].TimeVacc += 1 ;
				}
				if(Register[ic].TimeVacc>=Register[ic].ExpVacc){
					Register[ic].GotVacc = 0;
					for (xx = 0; xx < 13; xx++){
							Register[ic].VaccinationStatus[xx] = 0;
					}
				}	
			}
			if (Register[ic].VirginInd == 0){
				//if (HSVind == 1){ Register[ic].GetNewHSVstate(ic + 1, rSTI[ic][0]); }
				//if (TPind == 1){ Register[ic].GetNewTPstate(ic + 1, rSTI[ic][1]); }
				//if (HDind == 1){ Register[ic].GetNewHDstate(ic + 1, rSTI[ic][2]); }
				//if (NGind == 1){ Register[ic].GetNewNGstate(ic + 1, rSTI[ic][3]); }
				//if (CTind == 1){ Register[ic].GetNewCTstate(ic + 1, rSTI[ic][4]); }
				//if (TVind == 1){ Register[ic].GetNewTVstate(ic + 1, rSTI[ic][5]); }
				if (HPVind == 1){ 
					for (xx = 0; xx < 13; xx++){
						Register[ic].GetNewHPVstate(ic + 1, rSTI[ic][8 + xx], xx);
					}
					if(Register[ic].SexInd==1){
						Register[ic].TrueStage=0;
						if ((Indiv::AnyHPV(Register[ic].HPVstage,  Register[ic].allhpv, Register[ic].lsil)) && Register[ic].AliveInd==1) { 
							Register[ic].TrueStage=1;}
						if ((Indiv::AnyHPV(Register[ic].HPVstage,  Register[ic].allhpv, Register[ic].hsil)) && Register[ic].AliveInd==1) { 
							Register[ic].TrueStage=2;}
						if ((Indiv::AnyHPV(Register[ic].HPVstage,  Register[ic].allhpv, Register[ic].cc_un12)) && Register[ic].AliveInd==1) { 
							Register[ic].TrueStage=3;}
						if ((Indiv::AnyHPV(Register[ic].HPVstage,  Register[ic].allhpv, Register[ic].cc_diag)) && Register[ic].AliveInd==1) { 
							Register[ic].TrueStage=4;}
						if ((Indiv::AnyHPV(Register[ic].HPVstage,  Register[ic].allhpv, {15})) && Register[ic].AliveInd==1) { 
							Register[ic].TrueStage=5;}
						if ((Indiv::AnyHPV(Register[ic].HPVstage,  Register[ic].allhpv, Register[ic].cc_un34)) && Register[ic].AliveInd==1) { 
							Register[ic].TrueStage=6;}
	

						Register[ic].HPVstatus=0;
						if(Register[ic].TrueStage>0||Indiv::AnyHPV(Register[ic].HPVstage,  Register[ic].allhpv, {1})) {
							Register[ic].HPVstatus=1;}
					}				
					//else{ Register[ic].GetNewHPVstate(ic + 1, rSTI[ic][8 + WhichType], WhichType); }
				}
			}
			//if (Register[ic].SexInd == 1){
			//	if (BVind == 1){ BVtransitionF.GetNewStage(ic + 1, rSTI[ic][6]); }
			//	if (VCind == 1){ VCtransitionF.GetNewStage(ic + 1, rSTI[ic][7]); }
			//}
		}
		if (RoutineScreening==1 && CurrYear>1999 && HPVind == 1 && Register[ic].AliveInd == 1 && Register[ic].AgeGroup >= 3 && Register[ic].AgeGroup < 12 && 
			Register[ic].SexInd == 1 && Register[ic].DiagnosedCC == 0){
			if (WHOScreening==0 && PerfectSchedule==0 ){
				Register[ic].GetScreened(ic + 1, rSTI[ic][25], rSTI[ic][26], rSTI[ic][27], rSTI[ic][28], rSTI[ic][29], rSTI[ic][30], 
												rSTI[ic][31], rSTI[ic][32], rSTI[ic][33], rSTI[ic][34] , rSTI[ic][35],
												rSTI[ic][45], rSTI[ic][46], rSTI[ic][47], rSTI[ic][48], rSTI[ic][55], rSTI[ic][56], rSTI[ic][57], rSTI[ic][58]);
				if(Register[ic].timetoCol > 0) {
					Register[ic].GetTreated(ic + 1, rSTI[ic][36], rSTI[ic][37], rSTI[ic][38], rSTI[ic][39], rSTI[ic][40],
											rSTI[ic][41], rSTI[ic][42], rSTI[ic][43], rSTI[ic][44],
											rSTI[ic][49], rSTI[ic][50], rSTI[ic][51], rSTI[ic][52]);
				}
			}
			if (WHOScreening==1 && PerfectSchedule==0){
				if(CurrYear>1999 && CurrYear<ImplementYR){
					Register[ic].GetScreened(ic + 1, rSTI[ic][25], rSTI[ic][26], rSTI[ic][27], rSTI[ic][28], rSTI[ic][29], rSTI[ic][30], 
													rSTI[ic][31], rSTI[ic][32], rSTI[ic][33], rSTI[ic][34] , rSTI[ic][35],
												rSTI[ic][45], rSTI[ic][46], rSTI[ic][47], rSTI[ic][48], rSTI[ic][55], rSTI[ic][56],rSTI[ic][57],rSTI[ic][58]);
					if(Register[ic].timetoCol > 0) {
						Register[ic].GetTreated(ic + 1, rSTI[ic][36], rSTI[ic][37], rSTI[ic][38], rSTI[ic][39], rSTI[ic][40],
											rSTI[ic][41], rSTI[ic][42], rSTI[ic][43], rSTI[ic][44],
											rSTI[ic][49], rSTI[ic][50], rSTI[ic][51], rSTI[ic][52]);
					}
				}
				if(CurrYear>=ImplementYR ){
					Register[ic].WHOGetScreened(ic + 1, rSTI[ic][25], rSTI[ic][26], rSTI[ic][27], rSTI[ic][28], rSTI[ic][29], rSTI[ic][30], 
													rSTI[ic][31], rSTI[ic][32],  rSTI[ic][33], rSTI[ic][34] , rSTI[ic][35], rSTI[ic][36],
												rSTI[ic][46], rSTI[ic][47], rSTI[ic][48], rSTI[ic][49], rSTI[ic][55], rSTI[ic][56], rSTI[ic][57],rSTI[ic][58]);	
					if(Register[ic].timetoCol > 0) {
						Register[ic].GetTreated(ic + 1,  rSTI[ic][37], rSTI[ic][38], rSTI[ic][39], rSTI[ic][40],
												rSTI[ic][41], rSTI[ic][42], rSTI[ic][43], rSTI[ic][44], rSTI[ic][45],
											rSTI[ic][50], rSTI[ic][51], rSTI[ic][52], rSTI[ic][53]);
					}
				}
			}
			if (WHOScreening==0 && PerfectSchedule==1){
				if(CurrYear>1999 && CurrYear<ImplementYR){
					Register[ic].GetScreened(ic + 1, rSTI[ic][25], rSTI[ic][26], rSTI[ic][27], rSTI[ic][28], rSTI[ic][29], rSTI[ic][30], 
													rSTI[ic][31], rSTI[ic][32], rSTI[ic][33], rSTI[ic][34] , rSTI[ic][35],
												rSTI[ic][45], rSTI[ic][46], rSTI[ic][47], rSTI[ic][48], rSTI[ic][55], rSTI[ic][56], rSTI[ic][57],rSTI[ic][58]);
					if(Register[ic].timetoCol > 0) {
						Register[ic].GetTreated(ic + 1, rSTI[ic][36], rSTI[ic][37], rSTI[ic][38], rSTI[ic][39], rSTI[ic][40],
						rSTI[ic][41], rSTI[ic][42], rSTI[ic][43], rSTI[ic][44],
											rSTI[ic][49], rSTI[ic][50], rSTI[ic][51], rSTI[ic][52]);
					}
				}
				if(CurrYear>=ImplementYR){
					Register[ic].PerfectGetScreened(ic + 1, rSTI[ic][25], rSTI[ic][26], rSTI[ic][27], rSTI[ic][28], rSTI[ic][29], rSTI[ic][30], 
												rSTI[ic][31], rSTI[ic][32], rSTI[ic][33], rSTI[ic][34], rSTI[ic][35],
												rSTI[ic][45], rSTI[ic][46], rSTI[ic][47], rSTI[ic][48], rSTI[ic][55], rSTI[ic][56], rSTI[ic][57],rSTI[ic][58]);
					if(Register[ic].timetoCol > 0) {
						Register[ic].GetTreated(ic + 1, rSTI[ic][36], rSTI[ic][37], rSTI[ic][38], rSTI[ic][39], rSTI[ic][40],
											rSTI[ic][41], rSTI[ic][42], rSTI[ic][43], rSTI[ic][44],
											rSTI[ic][49], rSTI[ic][50], rSTI[ic][51], rSTI[ic][52]);
					}
				}
			}	
		}			
	}
}

Partner::Partner(){}

PartnerCohort::PartnerCohort()
{
	TotalDesire = 0.0;
}

void PartnerCohort::Reset()
{
	TotalDesire = 0.0;
	Pool.clear();
}

void PartnerCohort::AddMember(int PID)
{
	Partner A;
	A.ID = PID;
	A.DesiredNewPartners = Register[PID - 1].DesiredNewPartners;

	Pool.push_back(A);
	TotalDesire += A.DesiredNewPartners;
}

int PartnerCohort::SamplePool(double rand1)
{
	int ReachedEnd; // Indicator of whether search algorithm has finished
	int ic, PartnerID;
	double CumProb;

	rand1 *= TotalDesire;
	ic = 0;
	CumProb = 0.0;
	ReachedEnd = 0;
	while (ReachedEnd == 0){
		PartnerID = Pool[ic].ID;
		if (Register[PartnerID - 1].NewStatus == 0){
			CumProb += Pool[ic].DesiredNewPartners;
			if (rand1 < CumProb){
				ReachedEnd = 1;
			}
		}
		if (ic == Pool.size() - 1 && ReachedEnd == 0){
			// Prevent infinite loops that may arise due to rounding errors
			ReachedEnd = 1;
		}
		else{
			ic += 1;
		}
	}

	return PartnerID;
}

void ReadSexAssumps(const char *input)
{
	int ia, ib,  is;
	ifstream file;

	stringstream s;
	s << input;
	string path = "./input/" + s.str();
	
	file.open(path.c_str());
	if (file.fail()) {
		cerr << "Could not open input file.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	file >> HighPropnM >> HighPropnF;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> AssortativeM >> AssortativeF;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> GenderEquality;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> AnnNumberClients;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<4; ia++){
		file >> SexualDebut[ia][0] >> SexualDebut[ia][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> DebutAdjLow[0] >> DebutAdjLow[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PartnershipFormation[0][0] >> PartnershipFormation[1][0] >>
		PartnershipFormation[0][1] >> PartnershipFormation[1][1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BasePartnerAcqH[0] >> BasePartnerAcqH[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<16; ia++){
		file >> AgeEffectPartners[ia][0];
	}
	for (ia = 0; ia<16; ia++){
		file >> AgeEffectPartners[ia][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> GammaMeanST[0] >> GammaMeanST[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> GammaStdDevST[0] >> GammaStdDevST[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PartnerEffectNew[0][0] >> PartnerEffectNew[0][1] >> PartnerEffectNew[1][0] >>
		PartnerEffectNew[1][1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is<5; is++){
		file >> HIVeffectPartners[is];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<16; ia++){
		for (ib = 0; ib<2; ib++){
			file >> MarriageIncidence[ia][ib];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MeanFSWcontacts >> GammaMeanFSW >> GammaStdDevFSW;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ib = 0; ib<5; ib++){
		file >> PartnerEffectFSWcontact[ib];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<16; ia++){
		file >> InitFSWageDbn[ia] >> FSWexit[ia];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is<5; is++){
		file >> HIVeffectFSWentry[is];
	}
	for (is = 0; is<5; is++){
		file >> HIVeffectFSWexit[is];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MeanDurSTrel[0][0] >> MeanDurSTrel[0][1] >> MeanDurSTrel[1][0] >> MeanDurSTrel[1][1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<16; ia++){
		file >> LTseparation[ia][0] >> LTseparation[ia][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<16; ia++){
		for (ib = 0; ib<16; ib++){
			file >> AgePrefF[ia][ib];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<16; ia++){
		for (ib = 0; ib<16; ib++){
			file >> AgePrefM[ia][ib];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<16; ia++){
		file >> FreqSexST[ia][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<16; ia++){
		file >> FreqSexLT[ia][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BaselineCondomSvy;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ib = 0; ib<3; ib++){
		file >> RelEffectCondom[ib];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ib = 0; ib<3; ib++){
		file >> AgeEffectCondom[ib];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ib = 0; ib<3; ib++){
		file >> RatioInitialTo1998[ib];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ib = 0; ib<3; ib++){
		file >> RatioUltTo1998[ib];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ib = 0; ib<3; ib++){
		file >> MedianToBehavChange[ib];
	}
	file.ignore(255, '\n');
	for (ib = 0; ib<3; ib++){
		file >> MedianToBehavChange2[ib];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> SDpartnerRateAdj;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CondomScaling;
	file.close();

	// Calculate frequency of sex in men
	double sumxy;
	for (ia = 0; ia<16; ia++){
		sumxy = 0;
		for (ib = 0; ib<16; ib++){
			sumxy += AgePrefM[ia][ib] * FreqSexST[ib][1];
		}
		FreqSexST[ia][0] = sumxy;
	}
	for (ia = 0; ia<16; ia++){
		sumxy = 0;
		for (ib = 0; ib<16; ib++){
			sumxy += AgePrefM[ia][ib] * FreqSexLT[ib][1];
		}
		FreqSexLT[ia][0] = sumxy;
	}

	// Set BaselineCondomUse
	BaselineCondomUse = BaselineCondomSvy;

	// Calculate Weibull shape parameters for pace of behaviour change
	for (ib = 0; ib<3; ib++){
		ShapeBehavChange[ib] = log(log(1.0 - log(RatioInitialTo1998[ib]) /
			log(RatioUltTo1998[ib])) / log(2.0)) / log(13.0 / MedianToBehavChange[ib]);
	}

	// Calculate CumAgePrefM and CumAgePrefF
	for (ia = 0; ia<16; ia++){
		CumAgePrefM[ia][0] = AgePrefM[ia][0];
		CumAgePrefF[ia][0] = AgePrefF[ia][0];
		for (ib = 1; ib<15; ib++){
			CumAgePrefM[ia][ib] = CumAgePrefM[ia][ib - 1] + AgePrefM[ia][ib];
			CumAgePrefF[ia][ib] = CumAgePrefF[ia][ib - 1] + AgePrefF[ia][ib];
		}
		CumAgePrefM[ia][15] = 1.0;
		CumAgePrefF[ia][15] = 1.0;
	}
}

void ReadSTDepi(const char *input)
{
	int ia, is, iz, xx;
	ifstream file;
	stringstream s;
	s << input;
	string path = "./input/" + s.str();
	
	file.open(path.c_str());
	if (file.fail()) {
		cerr << "Could not open input file.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	file >> HSVtransitionM.AveDuration[0] >> HSVtransitionM.AveDuration[1] >>
		HSVtransitionM.AveDuration[2] >> HSVtransitionF.AveDuration[0] >>
		HSVtransitionF.AveDuration[1] >> HSVtransitionF.AveDuration[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iz = 0; iz<6; iz++){
		file >> TPtransitionM.AveDuration[iz];
	}
	for (iz = 0; iz<6; iz++){
		file >> TPtransitionF.AveDuration[iz];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HDtransitionM.AveDuration[0] >> HDtransitionM.AveDuration[1] >> HDtransitionM.AveDuration[2] >>
		HDtransitionF.AveDuration[0] >> HDtransitionF.AveDuration[1] >> HDtransitionF.AveDuration[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> NGtransitionM.AveDuration[0] >> NGtransitionM.AveDuration[1] >> NGtransitionM.AveDuration[2] >>
		NGtransitionF.AveDuration[0] >> NGtransitionF.AveDuration[1] >> NGtransitionF.AveDuration[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CTtransitionM.AveDuration[0] >> CTtransitionM.AveDuration[1] >> CTtransitionM.AveDuration[2] >>
		CTtransitionF.AveDuration[0] >> CTtransitionF.AveDuration[1] >> CTtransitionF.AveDuration[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> TVtransitionM.AveDuration[0] >> TVtransitionM.AveDuration[1] >> TVtransitionM.AveDuration[2] >>
		TVtransitionF.AveDuration[0] >> TVtransitionF.AveDuration[1] >> TVtransitionF.AveDuration[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is<5; is++){
		file >> HIVtransitionM.AveDuration[is];
	}
	for (is = 0; is<5; is++){
		file >> HIVtransitionF.AveDuration[is];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> VCtransitionF.AveDuration[0] >> VCtransitionF.AveDuration[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iz = 1; iz<4; iz++){
		file >> BVtransitionF.CtsTransition[iz][0];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BVtransitionF.CtsTransition[0][1] >> BVtransitionF.CtsTransition[2][1] >>
		BVtransitionF.CtsTransition[3][1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HSVtransitionM.RecurrenceRate >> HSVtransitionF.RecurrenceRate >>
		VCtransitionF.RecurrenceRate;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HDtransitionM.PropnImmuneAfterRx >> NGtransitionM.PropnImmuneAfterRx >>
		CTtransitionM.PropnImmuneAfterRx >> TVtransitionM.PropnImmuneAfterRx >>
		HDtransitionF.PropnImmuneAfterRx >> NGtransitionF.PropnImmuneAfterRx >>
		CTtransitionF.PropnImmuneAfterRx >> TVtransitionF.PropnImmuneAfterRx;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HDtransitionM.PropnImmuneAfterSR >> NGtransitionM.PropnImmuneAfterSR >>
		CTtransitionM.PropnImmuneAfterSR >> TVtransitionM.PropnImmuneAfterSR >>
		HDtransitionF.PropnImmuneAfterSR >> NGtransitionF.PropnImmuneAfterSR >>
		CTtransitionF.PropnImmuneAfterSR >> TVtransitionF.PropnImmuneAfterSR;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HSVtransitionM.SymptomaticPropn >> HDtransitionM.SymptomaticPropn >>
		NGtransitionM.SymptomaticPropn >> CTtransitionM.SymptomaticPropn >>
		TVtransitionM.SymptomaticPropn >> HSVtransitionF.SymptomaticPropn >>
		HDtransitionF.SymptomaticPropn >> NGtransitionF.SymptomaticPropn >>
		CTtransitionF.SymptomaticPropn >> TVtransitionF.SymptomaticPropn >>
		BVtransitionF.SymptomaticPropn;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HSVtransitionM.TransmProb >> TPtransitionM.TransmProb >> HDtransitionM.TransmProb >>
		NGtransitionM.TransmProb >> CTtransitionM.TransmProb >> TVtransitionM.TransmProb >>
		HSVtransitionF.TransmProb >> TPtransitionF.TransmProb >> HDtransitionF.TransmProb >>
		NGtransitionF.TransmProb >> CTtransitionF.TransmProb >> TVtransitionF.TransmProb >>
		HSVtransitionM.TransmProbSW >> TPtransitionM.TransmProbSW >> HDtransitionM.TransmProbSW >>
		NGtransitionM.TransmProbSW >> CTtransitionM.TransmProbSW >> TVtransitionM.TransmProbSW;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HSVtransitionM.RelTransmLT >> TPtransitionM.RelTransmLT >> HDtransitionM.RelTransmLT >>
		NGtransitionM.RelTransmLT >> CTtransitionM.RelTransmLT >> TVtransitionM.RelTransmLT >>
		HIVtransitionM.RelTransmLT >> HSVtransitionF.RelTransmLT >> TPtransitionF.RelTransmLT >>
		HDtransitionF.RelTransmLT >> NGtransitionF.RelTransmLT >> CTtransitionF.RelTransmLT >>
		TVtransitionF.RelTransmLT >> HIVtransitionF.RelTransmLT;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is<3; is++){
		file >> InitHIVtransm[is][0];
	}
	for (is = 0; is<3; is++){
		file >> InitHIVtransm[is][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	// Note that in the next line we are reading the male parameters into the female
	// arrays and the female parameters into the male arrays. This is deliberate; it makes
	// things a lot simpler when calculating the InfectProb arrays from the TransProb arrays
	// (in the STDtransition class).
	file >> HSVtransitionF.CondomEff >> TPtransitionF.CondomEff >> HDtransitionF.CondomEff >>
		NGtransitionF.CondomEff >> CTtransitionF.CondomEff >> TVtransitionF.CondomEff >>
		HIVtransitionF.CondomEff >>  HSVtransitionM.CondomEff >> TPtransitionM.CondomEff >>
		HDtransitionM.CondomEff >> NGtransitionM.CondomEff >> CTtransitionM.CondomEff >>
		TVtransitionM.CondomEff >> HIVtransitionM.CondomEff ;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (xx = 0; xx < 13; xx++){
		file >> HPVTransF[xx].CondomEff;
	}
	for (xx = 0; xx < 13; xx++){
		file >> HPVTransM[xx].CondomEff;
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iz = 0; iz<4; iz++){
		file >> HSVtransitionM.HIVinfecIncrease[iz];
	}
	for (iz = 0; iz<4; iz++){
		file >> HSVtransitionF.HIVinfecIncrease[iz];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iz = 0; iz<6; iz++){
		file >> TPtransitionM.HIVinfecIncrease[iz];
	}
	for (iz = 0; iz<6; iz++){
		file >> TPtransitionF.HIVinfecIncrease[iz];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HDtransitionM.HIVinfecIncrease[0] >> HDtransitionM.HIVinfecIncrease[1] >>
		HDtransitionF.HIVinfecIncrease[0] >> HDtransitionF.HIVinfecIncrease[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> NGtransitionM.HIVinfecIncrease[0] >> NGtransitionM.HIVinfecIncrease[1] >>
		NGtransitionF.HIVinfecIncrease[0] >> NGtransitionF.HIVinfecIncrease[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CTtransitionM.HIVinfecIncrease[0] >> CTtransitionM.HIVinfecIncrease[1] >>
		CTtransitionF.HIVinfecIncrease[0] >> CTtransitionF.HIVinfecIncrease[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> TVtransitionM.HIVinfecIncrease[0] >> TVtransitionM.HIVinfecIncrease[1] >>
		TVtransitionF.HIVinfecIncrease[0] >> TVtransitionF.HIVinfecIncrease[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is<6; is++){
		file >> HIVtransitionM.HIVinfecIncrease[is];
	}
	for (is = 0; is<6; is++){
		file >> HIVtransitionF.HIVinfecIncrease[is];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> VCtransitionF.HIVinfecIncrease[0] >> VCtransitionF.HIVinfecIncrease[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iz = 0; iz<3; iz++){
		file >> BVtransitionF.HIVinfecIncrease[iz];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iz = 0; iz<4; iz++){
		file >> HSVtransitionM.HIVsuscepIncrease[iz];
	}
	for (iz = 0; iz<4; iz++){
		file >> HSVtransitionF.HIVsuscepIncrease[iz];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iz = 0; iz<6; iz++){
		file >> TPtransitionM.HIVsuscepIncrease[iz];
	}
	for (iz = 0; iz<6; iz++){
		file >> TPtransitionF.HIVsuscepIncrease[iz];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HDtransitionM.HIVsuscepIncrease[0] >> HDtransitionM.HIVsuscepIncrease[1] >>
		HDtransitionF.HIVsuscepIncrease[0] >> HDtransitionF.HIVsuscepIncrease[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> NGtransitionM.HIVsuscepIncrease[0] >> NGtransitionM.HIVsuscepIncrease[1] >>
		NGtransitionF.HIVsuscepIncrease[0] >> NGtransitionF.HIVsuscepIncrease[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CTtransitionM.HIVsuscepIncrease[0] >> CTtransitionM.HIVsuscepIncrease[1] >>
		CTtransitionF.HIVsuscepIncrease[0] >> CTtransitionF.HIVsuscepIncrease[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> TVtransitionM.HIVsuscepIncrease[0] >> TVtransitionM.HIVsuscepIncrease[1] >>
		TVtransitionF.HIVsuscepIncrease[0] >> TVtransitionF.HIVsuscepIncrease[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> VCtransitionF.HIVsuscepIncrease[0] >> VCtransitionF.HIVsuscepIncrease[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iz = 0; iz<3; iz++){
		file >> BVtransitionF.HIVsuscepIncrease[iz];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	// Note that in the next few lines we are reading the male parameters into the female
	// arrays and the female parameters into the male arrays. This is deliberate; it makes
	// things a lot simpler when calculating the InfectProb arrays from the TransProb arrays
	// (in the STDtransition class).
	for (ia = 0; ia<16; ia++){
		file >> HSVtransitionF.SuscepIncrease[ia] >> TPtransitionF.SuscepIncrease[ia] >>
			HDtransitionF.SuscepIncrease[ia] >> NGtransitionF.SuscepIncrease[ia] >>
			CTtransitionF.SuscepIncrease[ia] >> TVtransitionF.SuscepIncrease[ia] >>
			HIVtransitionF.SuscepIncrease[ia];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<16; ia++){
		file >> HSVtransitionM.SuscepIncrease[ia] >> TPtransitionM.SuscepIncrease[ia] >>
			HDtransitionM.SuscepIncrease[ia] >> NGtransitionM.SuscepIncrease[ia] >>
			CTtransitionM.SuscepIncrease[ia] >> TVtransitionM.SuscepIncrease[ia] >>
			HIVtransitionM.SuscepIncrease[ia];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is<5; is++){
		file >> HSVsheddingIncrease[is];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is<5; is++){
		file >> HSVrecurrenceIncrease[is];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is<5; is++){
		file >> VCtransitionF.IncidenceIncrease[is];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> VCtransitionF.Incidence >> BVtransitionF.Incidence1;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BVtransitionF.IncidenceMultTwoPartners >> BVtransitionF.IncidenceMultNoPartners;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HSVsymptomInfecIncrease;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iz = 0; iz<3; iz++){
		file >> InfecIncreaseSyndrome[iz][0];
	}
	for (iz = 0; iz<3; iz++){
		file >> InfecIncreaseSyndrome[iz][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iz = 0; iz<3; iz++){
		file >> SuscepIncreaseSyndrome[iz][0];
	}
	for (iz = 0; iz<3; iz++){
		file >> SuscepIncreaseSyndrome[iz][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is<6; is++){
		file >> RelHIVfertility[is];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PropnInfectedAtBirth >> PropnInfectedAfterBirth;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MaleChild.Perinatal.PreAIDSmedian >> FemChild.Perinatal.PreAIDSmedian >>
		MaleChild.Breastmilk.PreAIDSmedian >> FemChild.Breastmilk.PreAIDSmedian;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MaleChild.Perinatal.PreAIDSshape >> FemChild.Perinatal.PreAIDSshape >>
		MaleChild.Breastmilk.PreAIDSshape >> FemChild.Breastmilk.PreAIDSshape;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MaleChild.Perinatal.MeanAIDSsurvival >> FemChild.Perinatal.MeanAIDSsurvival >>
		MaleChild.Breastmilk.MeanAIDSsurvival >> FemChild.Breastmilk.MeanAIDSsurvival;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MaleTeenRxRate >> MaleRxRate >> FemTeenRxRate >> FemRxRate >> FSWRxRate;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PropnTreatedPublicM >> PropnTreatedPublicF >> PropnTreatedPrivateM >>
		PropnTreatedPrivateF;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HSVtransitionF.CorrectRxPreSM >> TPtransitionF.CorrectRxPreSM >>
		HDtransitionF.CorrectRxPreSM >> NGtransitionF.CorrectRxPreSM >>
		CTtransitionF.CorrectRxPreSM >> TVtransitionF.CorrectRxPreSM >>
		BVtransitionF.CorrectRxPreSM >> VCtransitionF.CorrectRxPreSM;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HSVtransitionM.CorrectRxWithSM >> TPtransitionM.CorrectRxWithSM >>
		HDtransitionM.CorrectRxWithSM >> NGtransitionM.CorrectRxWithSM >>
		CTtransitionM.CorrectRxWithSM >> TVtransitionM.CorrectRxWithSM >>
		HSVtransitionF.CorrectRxWithSM >> TPtransitionF.CorrectRxWithSM >>
		HDtransitionF.CorrectRxWithSM >> NGtransitionF.CorrectRxWithSM >>
		CTtransitionF.CorrectRxWithSM >> TVtransitionF.CorrectRxWithSM >>
		BVtransitionF.CorrectRxWithSM;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HSVtransitionF.DrugEff >> TPtransitionF.DrugEff >> HDtransitionF.DrugEff >>
		NGtransitionF.DrugEff >> CTtransitionF.DrugEff >> TVtransitionF.DrugEff >>
		BVtransitionF.DrugEff >> VCtransitionF.DrugEff;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BVtransitionF.DrugPartialEff >> VCtransitionF.DrugPartialEff;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HSVtransitionF.TradnalEff >> TPtransitionF.TradnalEff >> HDtransitionF.TradnalEff >>
		NGtransitionF.TradnalEff >> CTtransitionF.TradnalEff >> TVtransitionF.TradnalEff >>
		BVtransitionF.TradnalEff >> VCtransitionF.TradnalEff;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> TPtransitionF.ANCpropnScreened >> TPtransitionF.ANCpropnTreated;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> AcceptScreening >> AcceptNVP >> RednNVP >> RednFF;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> SecondaryRxMult >> SecondaryCureMult >> TPtransitionM.PropnSuscepAfterRx >>
		TPtransitionF.PropnSuscepAfterRx;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> FSWasympRxRate >> FSWasympCure;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> InitHIVprevHigh;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RatioUltToInitHIVtransm;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> SDsuscepHIVadj;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iz = 0; iz<13; iz++){
		file >> VaccEfficacy[iz];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PropVaccinatedWHO;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PropVaccinatedHIV;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iz = 0; iz<13; iz++){
		file >> VaccEfficacyNONA[iz];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (int iz = 0; iz < 13; iz++) {
		file >> TxVEfficacyD1[iz];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (int iz = 0; iz < 13; iz++) {
		file >> TxVEfficacyD1CIN[iz]; 
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (int iz = 0; iz < 13; iz++) {
		file >> TxVEfficacyD2[iz];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (int iz = 0; iz < 13; iz++) {
		file >> TxVEfficacyD2CIN[iz]; 
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file.close();
	
	if (VaryParameters == 1){
		SimulateParameters();
	}
	// Set the effectiveness of SM for VC
	VCtransitionF.CorrectRxWithSM = VCtransitionF.CorrectRxPreSM;

	// Specify parameter values that apply to both sexes.
	HSVtransitionM.CorrectRxPreSM = HSVtransitionF.CorrectRxPreSM;
	TPtransitionM.CorrectRxPreSM = TPtransitionF.CorrectRxPreSM;
	HDtransitionM.CorrectRxPreSM = HDtransitionF.CorrectRxPreSM;
	NGtransitionM.CorrectRxPreSM = NGtransitionF.CorrectRxPreSM;
	CTtransitionM.CorrectRxPreSM = CTtransitionF.CorrectRxPreSM;
	TVtransitionM.CorrectRxPreSM = TVtransitionF.CorrectRxPreSM;

	HSVtransitionM.DrugEff = HSVtransitionF.DrugEff;
	TPtransitionM.DrugEff = TPtransitionF.DrugEff;
	HDtransitionM.DrugEff = HDtransitionF.DrugEff;
	NGtransitionM.DrugEff = NGtransitionF.DrugEff;
	CTtransitionM.DrugEff = CTtransitionF.DrugEff;
	TVtransitionM.DrugEff = TVtransitionF.DrugEff;

	HSVtransitionM.TradnalEff = HSVtransitionF.TradnalEff;
	TPtransitionM.TradnalEff = TPtransitionF.TradnalEff;
	HDtransitionM.TradnalEff = HDtransitionF.TradnalEff;
	NGtransitionM.TradnalEff = NGtransitionF.TradnalEff;
	CTtransitionM.TradnalEff = CTtransitionF.TradnalEff;
	TVtransitionM.TradnalEff = TVtransitionF.TradnalEff;

	// Convert annualized recurrence and incidence rates into weekly rates
	//HSVtransitionM.RecurrenceRate = HSVtransitionM.RecurrenceRate/52.0;
	//HSVtransitionF.RecurrenceRate = HSVtransitionF.RecurrenceRate/52.0;
	VCtransitionF.RecurrenceRate = VCtransitionF.RecurrenceRate / 52.0;
	VCtransitionF.Incidence = VCtransitionF.Incidence / 52.0;

	// Calculate remaining elements of the CtsTransition matrix and AveDuration for BV
	BVtransitionF.CtsTransition[1][2] = BVtransitionF.Incidence1 *
		BVtransitionF.SymptomaticPropn;
	BVtransitionF.CtsTransition[1][3] = BVtransitionF.Incidence1 *
		(1.0 - BVtransitionF.SymptomaticPropn);
	BVtransitionF.AveDuration[0] = 1.0 / BVtransitionF.CtsTransition[0][1];
	BVtransitionF.AveDuration[1] = 1.0 / (BVtransitionF.CtsTransition[1][0] +
		BVtransitionF.CtsTransition[1][2] + BVtransitionF.CtsTransition[1][3]);
	BVtransitionF.AveDuration[2] = 1.0 / (BVtransitionF.CtsTransition[2][0] +
		BVtransitionF.CtsTransition[2][1]);
	BVtransitionF.AveDuration[3] = 1.0 / (BVtransitionF.CtsTransition[3][0] +
		BVtransitionF.CtsTransition[3][1]);

	// Calculate RatioAsympToAveM and RatioAsympToAveF
	double sumx, sumxy;
	sumx = 0.0;
	sumxy = 0.0;
	for (is = 0; is<4; is++){
		sumx += HIVtransitionM.AveDuration[is];
		sumxy += HIVtransitionM.AveDuration[is] * HIVtransitionM.HIVinfecIncrease[is];
	}
	RatioAsympToAveM = sumx / (sumx + sumxy);
	sumx = 0.0;
	sumxy = 0.0;
	for (is = 0; is<4; is++){
		sumx += HIVtransitionF.AveDuration[is];
		sumxy += HIVtransitionF.AveDuration[is] * HIVtransitionF.HIVinfecIncrease[is];
	}
	RatioAsympToAveF = sumx / (sumx + sumxy);

	// Calculate average ART survival in kids
	MaleChild.Perinatal.AveYrsOnART = HIVtransitionM.AveDuration[4] / 52.0;
	MaleChild.Breastmilk.AveYrsOnART = HIVtransitionM.AveDuration[4] / 52.0;
	FemChild.Perinatal.AveYrsOnART = HIVtransitionF.AveDuration[4] / 52.0;
	FemChild.Breastmilk.AveYrsOnART = HIVtransitionF.AveDuration[4] / 52.0;

	// Set the HIV transmission probs in current year to their initial values
	HIVtransitionM.TransmProb[0] = InitHIVtransm[0][0];
	HIVtransitionF.TransmProb[0] = InitHIVtransm[0][1];
	for (is = 0; is<3; is++){
		HIVtransitionM.TransmProb[is + 1] = InitHIVtransm[1][0];
		HIVtransitionM.TransmProb[is + 4] = InitHIVtransm[2][0];
		HIVtransitionF.TransmProb[is + 1] = InitHIVtransm[1][1];
		HIVtransitionF.TransmProb[is + 4] = InitHIVtransm[2][1];
	}
	if (CofactorType == 0){
		HIVtransitionM.TransmProb[1] *= 2.0;
		HIVtransitionM.TransmProb[3] *= 0.5;
		HIVtransitionM.TransmProb[4] *= 2.0;
		HIVtransitionM.TransmProb[6] *= 0.5;
		HIVtransitionF.TransmProb[1] *= 2.0;
		HIVtransitionF.TransmProb[3] *= 0.5;
		HIVtransitionF.TransmProb[4] *= 2.0;
		HIVtransitionF.TransmProb[6] *= 0.5;
	}

	// Set the RelTransmCSW values
	HSVtransitionM.RelTransmCSW = HSVtransitionM.TransmProbSW / HSVtransitionM.TransmProb;
	TPtransitionM.RelTransmCSW = TPtransitionM.TransmProbSW / TPtransitionM.TransmProb;
	HDtransitionM.RelTransmCSW = HDtransitionM.TransmProbSW / HDtransitionM.TransmProb;
	NGtransitionM.RelTransmCSW = NGtransitionM.TransmProbSW / NGtransitionM.TransmProb;
	CTtransitionM.RelTransmCSW = CTtransitionM.TransmProbSW / CTtransitionM.TransmProb;
	TVtransitionM.RelTransmCSW = TVtransitionM.TransmProbSW / TVtransitionM.TransmProb;

	// Set the initial STD treatment parameters
	InitDrugEffNG = NGtransitionM.DrugEff;
	InitMaleRxRate = MaleRxRate;
	InitMaleTeenRxRate = MaleTeenRxRate;
	InitFemRxRate = FemRxRate;
	InitFemTeenRxRate = FemTeenRxRate;
	InitFSWasympRxRate = FSWasympRxRate;
	InitFSWasympCure = FSWasympCure;
	InitCorrectRxHSV = HSVtransitionM.CorrectRxWithSM;
	InitCorrectRxTVM = TVtransitionM.CorrectRxWithSM;
	InitANCpropnScreened = TPtransitionF.ANCpropnScreened;
	InitANCpropnTreated = TPtransitionF.ANCpropnTreated;
	//InitRecurrenceRateM = HSVtransitionM.RecurrenceRate;
	//InitRecurrenceRateF = HSVtransitionF.RecurrenceRate;
}

void ReadRatesByYear()
{
	int iy, ic;
	ifstream file;
	stringstream s;
	s << "RatesByYear.txt";
	string path = "./input/" + s.str();
	
	file.open(path.c_str());
	if (file.fail()) {
		cerr << "Could not open RatesByYear.txt\n";
		exit(1);
	}
	file.ignore(255,'\n');
	for(iy=0; iy<136; iy++){
		file>>PropnPrivateUsingSM[iy];}
	for(iy=0; iy<136; iy++){
		file>>PropnPublicUsingSM[iy];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iy=0; iy<136; iy++){
		file>>DrugShortage[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ic = 0; ic < 3; ic++){
		for (iy = 0; iy < 136; iy++){
			file >> RateARTstart[iy][ic][0];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ic = 0; ic < 3; ic++){
		for (iy = 0; iy < 136; iy++){
			file >> RateARTstart[iy][ic][1];}
	}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iy=0; iy<136; iy++){
		file>>HAARTaccess[iy];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iy=0; iy<136; iy++){
		file>>PMTCTaccess[iy];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iy=0; iy<136; iy++){
		file>>PMTCTaccess[iy];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iy=0; iy<136; iy++){
		file>>PropnCiproResistant[iy];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iy=0; iy<136; iy++){
		file>>PropnCiproTreated[iy];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iy=0; iy<136; iy++){
		file>>RxPhaseIn[iy];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iy=0; iy<136; iy++){
		file>>HIVtransitionM.ARTinfectiousness[iy];}
	for(iy=0; iy<136; iy++){
		file>>HIVtransitionF.ARTinfectiousness[iy];}	
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iy=0; iy<136; iy++){
		file>>ARTinterruption[iy];
	}	
	file.close();
	
}

void ReadMortTables()
{
	int ia, iy;
	ifstream file;
	stringstream s;
	s << "MortTables.txt";
	string path = "./input/" + s.str();
	
	file.open(path.c_str());
	if (file.fail()) {
		cerr << "Could not open MortTables.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	for (ia = 0; ia<16; ia++){
		for (iy = 0; iy<136; iy++){
			file >> NonAIDSmortM[ia][iy];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<16; ia++){
		for (iy = 0; iy<136; iy++){
			file >> NonAIDSmortF[ia][iy];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<136; iy++){
		file >> InfantMort1st6mM[iy];
	}
	for (ia = 0; ia<15; ia++){
		for (iy = 0; iy<136; iy++){
			file >> ChildMortM[ia][iy];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<136; iy++){
		file >> InfantMort1st6mF[iy];
	}
	for (ia = 0; ia<15; ia++){
		for (iy = 0; iy<136; iy++){
			file >> ChildMortF[ia][iy];
		}
	}
	file.close();
}

void ReadFertTables()
{
	int ia, iy;
	ifstream file;
	stringstream s;
	s << "FertTables.txt";
	string path = "./input/" + s.str();
	
	file.open(path.c_str());
	if (file.fail()) {
		cerr << "Could not open FertTables.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	for (ia = 0; ia<7; ia++){
		for (iy = 0; iy<136; iy++){
			file >> FertilityTable[ia][iy];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MaleBirthPropn;
	file.close();
}

void ReadOneStartProfileM(ifstream* file, int group)
{
	int ia, iz, offset, xx;
	double dummy;

	offset = 16 * group;
	for (ia = 0; ia<16; ia++){
		*file >> dummy;
		for (iz = 0; iz<5; iz++){
			*file >> HSVtransitionM.PropnByStage[ia + offset][iz];
		}
		for (iz = 0; iz<7; iz++){
			*file >> TPtransitionM.PropnByStage[ia + offset][iz];
		}
		for (iz = 0; iz<4; iz++){
			*file >> HDtransitionM.PropnByStage[ia + offset][iz];
		}
		for (iz = 0; iz<4; iz++){
			*file >> NGtransitionM.PropnByStage[ia + offset][iz];
		}
		for (iz = 0; iz<4; iz++){
			*file >> CTtransitionM.PropnByStage[ia + offset][iz];
		}
		for (iz = 0; iz<4; iz++){
			*file >> TVtransitionM.PropnByStage[ia + offset][iz];
		}
		for (xx = 0; xx < 13; xx++){
			for (iz = 0; iz < 8; iz++){
				*file >> HPVTransM[xx].PropnByStage[ia + offset][iz];
			}
		}
	}
}

void ReadOneStartProfileF(ifstream* file, int group)
{
	int ia, iz, offset, xx;
	double dummy;

	offset = 16 * group;
	for (ia = 0; ia<16; ia++){
		*file >> dummy;
		for (iz = 0; iz<5; iz++){
			*file >> HSVtransitionF.PropnByStage[ia + offset][iz];
		}
		for (iz = 0; iz<7; iz++){
			*file >> TPtransitionF.PropnByStage[ia + offset][iz];
		}
		for (iz = 0; iz<4; iz++){
			*file >> HDtransitionF.PropnByStage[ia + offset][iz];
		}
		for (iz = 0; iz<4; iz++){
			*file >> NGtransitionF.PropnByStage[ia + offset][iz];
		}
		for (iz = 0; iz<4; iz++){
			*file >> CTtransitionF.PropnByStage[ia + offset][iz];
		}
		for (iz = 0; iz<4; iz++){
			*file >> TVtransitionF.PropnByStage[ia + offset][iz];
		}
		for (iz = 0; iz<3; iz++){
			*file >> VCtransitionF.PropnByStage[ia + offset][iz];
		}
		for (iz = 0; iz<4; iz++){
			*file >> BVtransitionF.PropnByStage[ia + offset][iz];
		}
		
		for (xx = 0; xx < 13; xx++){
			for (iz = 0; iz < 8; iz++){
				*file >> HPVTransF[xx].PropnByStage[ia + offset][iz];
			}
		}
		
	}
}

void ReadStartProfile(const char *input)
{
	int ia, iz, xx;
	double dummy;
	ifstream file;
	stringstream s;
	s << input;
	string path = "./input/" + s.str();
	
	file.open(path.c_str());
	if (file.fail()) {
		cerr << "File open error\n";
		exit(1);
	}
	file.ignore(255, '\n');
	for (ia = 0; ia<4; ia++){
		file >> dummy;
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<4; ia++){
		file >> dummy;
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<4; ia++){
		file >> dummy;
		for (iz = 0; iz<3; iz++){
			file >> VCtransitionF.PropnByStage[ia][iz];
		}
		for (iz = 0; iz<4; iz++){
			file >> BVtransitionF.PropnByStage[ia][iz];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<4; ia++){
		file >> dummy;
		for (iz = 0; iz<3; iz++){
			file >> VCtransitionF.PropnByStage[16 + ia][iz];
		}
		for (iz = 0; iz<4; iz++){
			file >> BVtransitionF.PropnByStage[16 + ia][iz];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 2); // MH, 0
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 3); // MH, 1, S
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 4); // MH, 2, S
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 5); // MH, 1, L
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 6); // MH, 2, L
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 7); // MH, 11, S
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 8); // MH, 12, S
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 9); // MH, 22, S
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 10); // MH, 11, L
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 11); // MH, 12, L
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 12); // MH, 21, L
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 13); // MH, 22, L
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 14); // ML, 0
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 15); // ML, 1, S
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 16); // ML, 2, S
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 17); // ML, 1, L
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileM(&file, 18); // ML, 2, L
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 19); // FSW
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 2); // FH, 0
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 3); // FH, 1, S
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 4); // FH, 2, S
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 5); // FH, 1, L
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 6); // FH, 2, L
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 7); // FH, 11, S
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 8); // FH, 12, S
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 9); // FH, 22, S
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 10); // FH, 11, L
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 11); // FH, 12, L
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 12); // FH, 21, L
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 13); // FH, 22, L
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 14); // FL, 0
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 15); // FL, 1, S
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 16); // FL, 2, S
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 17); // FL, 1, L
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	ReadOneStartProfileF(&file, 18); // FL, 2, L
	// Note that we aren't reading the rest of the file, as the initial child numbers
	// are now in StartPop and the male birth propn is in FertTables.
	file.close();
	
	// Set the proportions of uninfected virgins to 100%.
	for (ia = 0; ia<4; ia++){
		if (HSVind == 1){
			HSVtransitionM.PropnByStage[ia][0] = 1.0;
			HSVtransitionM.PropnByStage[ia + 16][0] = 1.0;
			HSVtransitionF.PropnByStage[ia][0] = 1.0;
			HSVtransitionF.PropnByStage[ia + 16][0] = 1.0;
		}
		if (TPind == 1){
			TPtransitionM.PropnByStage[ia][0] = 1.0;
			TPtransitionM.PropnByStage[ia + 16][0] = 1.0;
			TPtransitionF.PropnByStage[ia][0] = 1.0;
			TPtransitionF.PropnByStage[ia + 16][0] = 1.0;
		}
		if (HDind == 1){
			HDtransitionM.PropnByStage[ia][0] = 1.0;
			HDtransitionM.PropnByStage[ia + 16][0] = 1.0;
			HDtransitionF.PropnByStage[ia][0] = 1.0;
			HDtransitionF.PropnByStage[ia + 16][0] = 1.0;
		}
		if (NGind == 1){
			NGtransitionM.PropnByStage[ia][0] = 1.0;
			NGtransitionM.PropnByStage[ia + 16][0] = 1.0;
			NGtransitionF.PropnByStage[ia][0] = 1.0;
			NGtransitionF.PropnByStage[ia + 16][0] = 1.0;
		}
		if (CTind == 1){
			CTtransitionM.PropnByStage[ia][0] = 1.0;
			CTtransitionM.PropnByStage[ia + 16][0] = 1.0;
			CTtransitionF.PropnByStage[ia][0] = 1.0;
			CTtransitionF.PropnByStage[ia + 16][0] = 1.0;
		}
		if (TVind == 1){
			TVtransitionM.PropnByStage[ia][0] = 1.0;
			TVtransitionM.PropnByStage[ia + 16][0] = 1.0;
			TVtransitionF.PropnByStage[ia][0] = 1.0;
			TVtransitionF.PropnByStage[ia + 16][0] = 1.0;
		}
		if (HPVind == 1){
			for (xx = 0; xx < 13; xx++){
				HPVTransM[xx].PropnByStage[ia][0] = 1.0;
				HPVTransM[xx].PropnByStage[ia + 16][0] = 1.0;
				HPVTransF[xx].PropnByStage[ia][0] = 1.0;
				HPVTransF[xx].PropnByStage[ia + 16][0] = 1.0;
			}
		}
	}
}

void ReadInitHIV()
{
	int ia;
	ifstream file;
	stringstream s;
	s << "InitHIV.txt";
	string path = "./input/" + s.str();
	
	file.open(path.c_str());
	if (file.fail()) {
		cerr << "Could not open InitHIV.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	for (ia = 0; ia<7; ia++){
		file >> HlabisaRatio[ia][0] >> HlabisaRatio[ia][1];
	}
	file.close();
}

void ReadStartPop()
{
	int ia;
	ifstream file;
	stringstream s;
	s << "StartPop.txt";
	string path = "./input/" + s.str();
	
	file.open(path.c_str());
	if (file.fail()) {
		cerr << "Could not open StartPop.txt\n";
		exit(1);
	}
	for (ia = 0; ia<91; ia++){
		file >> StartPop[ia][0] >> StartPop[ia][1];
	}
	file.close();
}

void ReadSTDprev() 
{
	if (HIVcalib == 1){
		HIVtransitionM.ReadPrevData("HIVdataM.txt");
		HIVtransitionF.ReadPrevData("HIVdataF.txt");
	}
	if (HSVcalib == 1){
		HSVtransitionM.ReadPrevData("HSVdataM.txt");
		HSVtransitionF.ReadPrevData("HSVdataF.txt");
	}
	if (TPcalib == 1){
		TPtransitionM.ReadPrevData("TPdataM.txt");
		TPtransitionF.ReadPrevData("TPdataF.txt");
	}
	if (HDcalib == 1){
		HDtransitionM.ReadPrevData("HDdataM.txt");
		HDtransitionF.ReadPrevData("HDdataF.txt");
	}
	if (NGcalib == 1){
		NGtransitionM.ReadPrevData("NGdataM.txt");
		NGtransitionF.ReadPrevData("NGdataF.txt");
	}
	if (CTcalib == 1){
		CTtransitionM.ReadPrevData("CTdataM.txt");
		CTtransitionF.ReadPrevData("CTdataF.txt");
	}
	if (TVcalib == 1){
		TVtransitionM.ReadPrevData("TVdataM.txt");
		TVtransitionF.ReadPrevData("TVdataF.txt");
	}
	if (BVcalib == 1){ BVtransitionF.ReadPrevData("BVdataF.txt"); }
	if (VCcalib == 1){ VCtransitionF.ReadPrevData("VCdataF.txt"); }
	
	if (HPVcalib == 1){
		HPVtransitionCC.ReadPrevData("CCStageData.txt");
	}
}

void ReadSTDparameters()
{
	int ic, iz; // , dummy;
	ifstream file, file2, file3, file4, file5, file6, file7, file8, file9;

	if (HIVcalib == 1){
		file.open("RandomUniformHIV.txt");
		if (file.fail()) {
			cerr << "Could not open RandomUniformHIV.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformHIV.columns; iz++){
				file >> RandomUniformHIV.out[ic][iz];
			}
		}
		file.close();
	}

	if (HSVcalib == 1){
		file2.open("RandomUniformHSV.txt");
		if (file2.fail()) {
			cerr << "Could not open RandomUniformHSV.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file2 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformHSV.columns; iz++){
				file2 >> RandomUniformHSV.out[ic][iz];
			}
		}
		file2.close();
	}

	if (NGcalib == 1){
		file4.open("RandomUniformNG.txt");
		if (file4.fail()) {
			cerr << "Could not open RandomUniformNG.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file4 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformNG.columns; iz++){
				file4 >> RandomUniformNG.out[ic][iz];
			}
		}
		file4.close();
	}

	if (CTcalib == 1){
		file5.open("RandomUniformCT.txt");
		if (file5.fail()) {
			cerr << "Could not open RandomUniformCT.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file5 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformCT.columns; iz++){
				file5 >> RandomUniformCT.out[ic][iz];
			}
		}
		file5.close();
	}

	if (TVcalib == 1){
		file6.open("RandomUniformTV.txt");
		if (file6.fail()) {
			cerr << "Could not open RandomUniformTV.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file6 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformTV.columns; iz++){
				file6 >> RandomUniformTV.out[ic][iz];
			}
		}
		file6.close();
	}

	if (TPcalib == 1){
		file3.open("RandomUniformTP.txt");
		if (file3.fail()) {
			cerr << "Could not open RandomUniformTP.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file3 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformTP.columns; iz++){
				file3 >> RandomUniformTP.out[ic][iz];
			}
		}
		file3.close();
	}

	if (BVcalib == 1){
		file7.open("RandomUniformBV.txt");
		if (file7.fail()) {
			cerr << "Could not open RandomUniformBV.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file7 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformBV.columns; iz++){
				file7 >> RandomUniformBV.out[ic][iz];
			}
		}
		file7.close();
	}

	if (VCcalib == 1){
		file8.open("RandomUniformVC.txt");
		if (file8.fail()) {
			cerr << "Could not open RandomUniformVC.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file8 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformVC.columns; iz++){
				file8 >> RandomUniformVC.out[ic][iz];
			}
		}
		file8.close();
	}
	if (HPVind == 1){
			stringstream s;
			ifstream file9;
			s << process_num << "RandomUniformHPV.txt";
			
			string path = "./randomHPV/" + s.str();
				
			file9.open(path.c_str()); // Converts s to a C string
			
			if (file9.fail()) {
				cerr << "Could not open RandomUniformHPV.txt\n";
				exit(1);
			}
			for (ic = 0; ic < ParamCombs; ic++){
				file9 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
				for (iz = 0; iz < 34; iz++){
					file9 >> RandomUniformHPV.out[ic][iz];
				}
			}
			file9.close();
	}
	/*if (HPVcalib == 1){
		file9.open("RandomUniformHPV.txt");
		if (file9.fail()) {
			cerr << "Could not open RandomUniformHPV.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file9 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < 38; iz++){
				file9 >> RandomUniformHPV.out[ic][iz];
			}
		}
		file9.close();
	}*/	
	/*
	if (HPVcalib == 1){
		if (OneType == 1 || targets == 1){
			//stringstream s;
			ifstream file9;
			//s << WhichType << "RandomUniformHPV.txt";
			
			//file9.open(s.str().c_str()); // Converts s to a C string
			
			if(WhichType==0) {file9.open("0RandomUniformHPV.txt");}
            else if(WhichType==1) {file9.open("1RandomUniformHPV.txt");}
            else if(WhichType==2) {file9.open("2RandomUniformHPV.txt");}
            else if(WhichType==3) {file9.open("3RandomUniformHPV.txt");}
            else if(WhichType==4) {file9.open("4RandomUniformHPV.txt");}
            else if(WhichType==5) {file9.open("5RandomUniformHPV.txt");}
            else if(WhichType==6) {file9.open("6RandomUniformHPV.txt");}
            else if(WhichType==7) {file9.open("7RandomUniformHPV.txt");}
            else if(WhichType==8) {file9.open("8RandomUniformHPV.txt");}
            else if(WhichType==9) {file9.open("9RandomUniformHPV.txt");}
            else if(WhichType==10) {file9.open("10RandomUniformHPV.txt");}
            else if(WhichType==11) {file9.open("11RandomUniformHPV.txt");}
            else if(WhichType==12) {file9.open("12RandomUniformHPV.txt");}

			if (file9.fail()) {
				cerr << "Could not open RandomUniformHPV.txt\n";
				exit(1);
			}
			for (ic = 0; ic < ParamCombs; ic++){
				file9 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
				for (iz = 0; iz < 16; iz++){
					file9 >> RandomUniformHPV[WhichType].out[ic][iz];
				}
			}
			file9.close();
		}
		else{
			ifstream file0, file1, file2, file3, file4, file5, file6, file7, file8, file9, file10, file11, file12;
			file0.open("0RandomUniformHPV.txt");
			for (ic = 0; ic < ParamCombs; ic++){
				file0 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
				for (iz = 0; iz < 16; iz++){
					file0 >> RandomUniformHPV[0].out[ic][iz];
				}
			}
			file0.close();
			file1.open("1RandomUniformHPV.txt");
			for (ic = 0; ic < ParamCombs; ic++){
				file1 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
				for (iz = 0; iz < 16; iz++){
					file1 >> RandomUniformHPV[1].out[ic][iz];
				}
			}
			file1.close();
			file2.open("2RandomUniformHPV.txt");
			for (ic = 0; ic < ParamCombs; ic++){
				file2 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
				for (iz = 0; iz < 16; iz++){
					file2 >> RandomUniformHPV[2].out[ic][iz];
				}
			}
			file2.close();
			file3.open("3RandomUniformHPV.txt");
			for (ic = 0; ic < ParamCombs; ic++){
				file3 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
				for (iz = 0; iz < 16; iz++){
					file3 >> RandomUniformHPV[3].out[ic][iz];
				}
			}
			file3.close();
			file4.open("4RandomUniformHPV.txt");
			for (ic = 0; ic < ParamCombs; ic++){
				file4 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
				for (iz = 0; iz < 16; iz++){
					file4 >> RandomUniformHPV[4].out[ic][iz];
				}
			}
			file4.close();
			file5.open("5RandomUniformHPV.txt");
			for (ic = 0; ic < ParamCombs; ic++){
				file5 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
				for (iz = 0; iz < 16; iz++){
					file5 >> RandomUniformHPV[5].out[ic][iz];
				}
			}
			file5.close();
			file6.open("6RandomUniformHPV.txt");
			for (ic = 0; ic < ParamCombs; ic++){
				file6 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
				for (iz = 0; iz < 16; iz++){
					file6 >> RandomUniformHPV[6].out[ic][iz];
				}
			}
			file6.close();
			file7.open("7RandomUniformHPV.txt");
			for (ic = 0; ic < ParamCombs; ic++){
				file7 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
				for (iz = 0; iz < 16; iz++){
					file7 >> RandomUniformHPV[7].out[ic][iz];
				}
			}
			file7.close();
			file8.open("8RandomUniformHPV.txt");
			for (ic = 0; ic < ParamCombs; ic++){
				file8 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
				for (iz = 0; iz < 16; iz++){
					file8 >> RandomUniformHPV[8].out[ic][iz];
				}
			}
			file8.close();
			file9.open("9RandomUniformHPV.txt");
			for (ic = 0; ic < ParamCombs; ic++){
				file9 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
				for (iz = 0; iz < 16; iz++){
					file9 >> RandomUniformHPV[9].out[ic][iz];
				}
			}
			file9.close();
			file10.open("10RandomUniformHPV.txt"); 
			for (ic = 0; ic < ParamCombs; ic++){
				file10 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
				for (iz = 0; iz < 16; iz++){
					file10 >> RandomUniformHPV[10].out[ic][iz];
				}
			}
			file10.close();
			file11.open("11RandomUniformHPV.txt");
			for (ic = 0; ic < ParamCombs; ic++){
				file11 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
				for (iz = 0; iz < 16; iz++){
					file11 >> RandomUniformHPV[11].out[ic][iz];
				}
			}
			file11.close();
			file12.open("12RandomUniformHPV.txt");
			for (ic = 0; ic < ParamCombs; ic++){
				file12 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
				for (iz = 0; iz < 16; iz++){
					file12 >> RandomUniformHPV[12].out[ic][iz];
				}
			}
			file12.close();

		}
		
	}*/
}

void ReadAllInputFiles()
{
	ReadCCStrategies();
	ReadSexAssumps("SexAssumps.txt");
	ReadSTDepi("STDepidemiology.txt");
	ReadRatesByYear();
	ReadMortTables();
	ReadFertTables();
	ReadStartProfile("StartProfile.txt");
	ReadInitHIV();
	ReadStartPop();
	ReadSTDprev();
	ReadTimeinCIN3();
	ReadScreenData();
	ReadLifeTimePartners();
	
}

void SetCalibParameters()
{
	int ia;
	double AveOR;

	AveOR = 0.68; // Posterior mean from the ASSA2002 uncertainty analysis

	// Values in the 5 lines below are from the simulateparameters function in the C++
	// version of the ASSA2002 model (see derivation in section 3.3.1.3 of the uncertainty
	// analysis report)
	/*HIVtransitionF.AntenatalNlogL.BiasMult[1] = 0.885450 / (AveOR*0.044021 + 0.841431);
	HIVtransitionF.AntenatalNlogL.BiasMult[2] = 0.902648 / (AveOR*0.080730 + 0.810868);
	HIVtransitionF.AntenatalNlogL.BiasMult[3] = 0.870821 / (AveOR*0.108702 + 0.707487);
	HIVtransitionF.AntenatalNlogL.BiasMult[4] = 0.921259 / (AveOR*0.142296 + 0.698347);
	HIVtransitionF.AntenatalNlogL.BiasMult[5] = 0.937727 / (AveOR*0.117837 + 0.784181);*/

	// FPC weights are the proportions of women currently using any modern contraceptive
	// method, as estimated from the 1998 DHS. By assuming the same proportions for single
	// women and women currently in partnerships, we may be under-estimating STD prevalence,
	// but by setting the proportions to 0 in the <15 and 50+ age categories, we may be
	// over-estimating STD prevalence. So the biases cancel out to some extent.
				FPCweights[0] = 0.0;
				FPCweights[1] = 0.285;
				FPCweights[2] = 0.565;
                FPCweights[3] = 0.578;
                FPCweights[4] = 0.586;
                FPCweights[5] = 0.558;
                FPCweights[6] = 0.502;
                FPCweights[7] = 0.380;
                for (ia = 8; ia < 16; ia++){
                    FPCweights[ia] = 0.0;
                }


	// Setting the default values for the variance of the study effects (only relevant to
	// the sentinel surveillance data)
	if (HIVcalib == 1 && HIVtransitionF.ANClogL.VarStudyEffect == 0){
		HIVtransitionM.SetVarStudyEffect(0.09);
		HIVtransitionF.SetVarStudyEffect(0.09);
	}
	if (HSVcalib == 1 && HSVtransitionF.ANClogL.VarStudyEffect == 0){
		HSVtransitionM.SetVarStudyEffect(0.48);
		HSVtransitionF.SetVarStudyEffect(0.48);
	}
	if (TPcalib == 1 && TPtransitionF.ANClogL.VarStudyEffect == 0){
		TPtransitionM.SetVarStudyEffect(0.14);
		TPtransitionF.SetVarStudyEffect(0.14);
	}
	if (HDcalib == 1 && HDtransitionF.ANClogL.VarStudyEffect == 0){
		HDtransitionM.SetVarStudyEffect(0.30);
		HDtransitionF.SetVarStudyEffect(0.30);
	}
	if (NGcalib == 1 && NGtransitionF.ANClogL.VarStudyEffect == 0){
		NGtransitionM.SetVarStudyEffect(0.09);
		NGtransitionF.SetVarStudyEffect(0.09);
	}
	if (CTcalib == 1 && CTtransitionF.ANClogL.VarStudyEffect == 0){
		CTtransitionM.SetVarStudyEffect(0.11);
		CTtransitionF.SetVarStudyEffect(0.11);
	}
	if (TVcalib == 1 && TVtransitionF.ANClogL.VarStudyEffect == 0){
		TVtransitionM.SetVarStudyEffect(0.28);
		TVtransitionF.SetVarStudyEffect(0.28);
	}
	if (BVcalib == 1 && BVtransitionF.ANClogL.VarStudyEffect == 0){
		BVtransitionF.SetVarStudyEffect(0.0529);
	}
	if (VCcalib == 1 && VCtransitionF.ANClogL.VarStudyEffect == 0){
		VCtransitionF.SetVarStudyEffect(0.0529);
	}
}

void CalcTotalLogL()
{
	double ModelErrorVar;
	int xx;
	TotalLogL = 0.0;
	if (HIVcalib == 1){
		//HIVtransitionF.CSWlogL.CalcLogL();
		if (GetSDfromData == 1){
			HIVtransitionF.AntenatalNlogL.CalcModelVar();
			HIVtransitionF.HouseholdNlogL.CalcModelVar();
			HIVtransitionM.HouseholdNlogL.CalcModelVar();
			ModelErrorVar = (HIVtransitionF.HouseholdNlogL.ModelVarEst * HIVtransitionF.HouseholdNlogL.Observations +
				HIVtransitionM.HouseholdNlogL.ModelVarEst * HIVtransitionM.HouseholdNlogL.Observations +
				HIVtransitionF.AntenatalNlogL.ModelVarEst * HIVtransitionF.AntenatalNlogL.Observations) /
				(HIVtransitionF.HouseholdNlogL.Observations + HIVtransitionM.HouseholdNlogL.Observations +
				HIVtransitionF.AntenatalNlogL.Observations);
			if (ModelErrorVar < 0.0){ ModelErrorVar = 0.0; }
			HIVtransitionF.HouseholdNlogL.ModelVarEst = ModelErrorVar;
			HIVtransitionM.HouseholdNlogL.ModelVarEst = ModelErrorVar;
			HIVtransitionF.AntenatalNlogL.ModelVarEst = ModelErrorVar;
		}
		HIVtransitionF.AntenatalNlogL.CalcLogL();
		HIVtransitionF.HouseholdNlogL.CalcLogL();
		HIVtransitionM.HouseholdNlogL.CalcLogL();
		TotalLogL += HIVtransitionF.CSWlogL.LogL + HIVtransitionF.AntenatalNlogL.LogL +
			HIVtransitionF.HouseholdNlogL.LogL + HIVtransitionM.HouseholdNlogL.LogL;
	}
	if (HSVcalib == 1){
		HSVtransitionF.ANClogL.CalcLogL();
		HSVtransitionF.FPClogL.CalcLogL();
		HSVtransitionF.CSWlogL.CalcLogL();
		HSVtransitionF.HouseholdLogL.CalcLogL();
		HSVtransitionM.HouseholdLogL.CalcLogL();
		HSVparamsLogL.out[CurrSim - 1][11] = HSVtransitionF.ANClogL.LogL + HSVtransitionF.FPClogL.LogL +
			HSVtransitionF.CSWlogL.LogL + HSVtransitionF.HouseholdLogL.LogL +
			HSVtransitionM.HouseholdLogL.LogL;
	}
	if (TPcalib == 1){
		TPtransitionF.AntenatalNlogL.CalcLogL();
		TPtransitionF.ANClogL.CalcLogL();
		TPtransitionF.FPClogL.CalcLogL();
		TPtransitionF.CSWlogL.CalcLogL();
		TPtransitionF.HouseholdLogL.CalcLogL();
		TPtransitionM.HouseholdLogL.CalcLogL();
		TPparamsLogL.out[CurrSim - 1][10] = TPtransitionF.AntenatalNlogL.LogL + TPtransitionF.ANClogL.LogL +
			TPtransitionF.FPClogL.LogL + TPtransitionF.CSWlogL.LogL +
			TPtransitionF.HouseholdLogL.LogL + TPtransitionM.HouseholdLogL.LogL;
	}
	/*if (HDcalib == 1){
	HDtransitionF.GUDlogL.CalcLogL();
	HDtransitionM.GUDlogL.CalcLogL();
	TotalLogL += HDtransitionF.GUDlogL.LogL + HDtransitionM.GUDlogL.LogL;
	}*/
	if (NGcalib == 1){
		NGtransitionF.ANClogL.CalcLogL();
		NGtransitionF.FPClogL.CalcLogL();
		NGtransitionF.CSWlogL.CalcLogL();
		NGtransitionF.HouseholdLogL.CalcLogL();
		NGtransitionM.HouseholdLogL.CalcLogL();
		NGparamsLogL.out[CurrSim - 1][10] = NGtransitionF.ANClogL.LogL + NGtransitionF.FPClogL.LogL +
			NGtransitionF.CSWlogL.LogL + NGtransitionF.HouseholdLogL.LogL +
			NGtransitionM.HouseholdLogL.LogL;
	}
	if (CTcalib == 1){
		CTtransitionF.ANClogL.CalcLogL();
		CTtransitionF.FPClogL.CalcLogL();
		CTtransitionF.CSWlogL.CalcLogL();
		CTtransitionF.HouseholdLogL.CalcLogL();
		CTtransitionM.HouseholdLogL.CalcLogL();
		CTparamsLogL.out[CurrSim - 1][10] = CTtransitionF.ANClogL.LogL + CTtransitionF.FPClogL.LogL +
			CTtransitionF.CSWlogL.LogL + CTtransitionF.HouseholdLogL.LogL +
			CTtransitionM.HouseholdLogL.LogL;
	}
	if (TVcalib == 1){
		TVtransitionF.ANClogL.CalcLogL();
		TVtransitionF.FPClogL.CalcLogL();
		TVtransitionF.CSWlogL.CalcLogL();
		TVtransitionF.HouseholdLogL.CalcLogL();
		TVtransitionM.HouseholdLogL.CalcLogL();
		TVparamsLogL.out[CurrSim - 1][12] = TVtransitionF.ANClogL.LogL + TVtransitionF.FPClogL.LogL +
			TVtransitionF.CSWlogL.LogL + TVtransitionF.HouseholdLogL.LogL +
			TVtransitionM.HouseholdLogL.LogL;
	}
	int SimCount2 = (CurrSim) / IterationsPerPC;
		
	if (HPVcalib == 1){ 
		if (CurrSim == (SimCount2 * IterationsPerPC)){HPVparamsLogL.out[SimCount2-1][0] = 0.0;}
		HPVtransitionCC.HouseholdLogL.CalcLogL();
		HPVtransitionCC.NOARTlogL.CalcLogL();
		HPVtransitionCC.ONARTlogL.CalcLogL();
		HPVtransitionCC.FPClogL.CalcLogL();
		//HPVtransitionCC.HouseholdNlogL.CalcLogL();
		if (CurrSim == (SimCount2 * IterationsPerPC)){
			HPVparamsLogL.out[SimCount2-1][0] += HPVtransitionCC.FPClogL.LogL + HPVtransitionCC.HouseholdLogL.LogL +
			 HPVtransitionCC.NOARTlogL.LogL + HPVtransitionCC.ONARTlogL.LogL ;
			 HPVparamsLogL.out[SimCount2-1][1] += HPVtransitionCC.FPClogL.LogL0 + HPVtransitionCC.HouseholdLogL.LogL0 +
			 HPVtransitionCC.NOARTlogL.LogL0 + HPVtransitionCC.ONARTlogL.LogL0 ;
			 HPVparamsLogL.out[SimCount2-1][2] += HPVtransitionCC.FPClogL.LogL1 + HPVtransitionCC.HouseholdLogL.LogL1 +
			 HPVtransitionCC.NOARTlogL.LogL1 + HPVtransitionCC.ONARTlogL.LogL1 ;
			 HPVparamsLogL.out[SimCount2-1][3] += HPVtransitionCC.FPClogL.LogL2 + HPVtransitionCC.HouseholdLogL.LogL2 +
			 HPVtransitionCC.NOARTlogL.LogL2 + HPVtransitionCC.ONARTlogL.LogL2 ;
			 HPVparamsLogL.out[SimCount2-1][4] += HPVtransitionCC.FPClogL.LogL3 + HPVtransitionCC.HouseholdLogL.LogL3 +
			 HPVtransitionCC.NOARTlogL.LogL3 + HPVtransitionCC.ONARTlogL.LogL3 ; 
		}
	}

	if (BVcalib == 1){
		BVtransitionF.ANClogL.CalcLogL();
		BVtransitionF.FPClogL.CalcLogL();
		BVtransitionF.CSWlogL.CalcLogL();
		TotalLogL += BVtransitionF.ANClogL.LogL + BVtransitionF.FPClogL.LogL +
		BVtransitionF.CSWlogL.LogL;
	}
	if (VCcalib == 1){
		VCtransitionF.ANClogL.CalcLogL();
		VCtransitionF.FPClogL.CalcLogL();
		VCtransitionF.CSWlogL.CalcLogL();
		TotalLogL += VCtransitionF.ANClogL.LogL + VCtransitionF.FPClogL.LogL +
		VCtransitionF.CSWlogL.LogL;
	}
}

void OneSimulation()
{
	int seed, SimCount2;
	int tpp = Register.size();
	for (int i = 0; i < tpp; ++i){
		memset(&Register[i], 0, sizeof(Register[0]));
	}
	ReadAllInputFiles();
	
	// Reset the seed
	if (FixedUncertainty == 0){
		SimCount2 = (CurrSim - 1) / IterationsPerPC;
		seed = SimCount2 * 91 + process_num * 7927 + (CurrSim - SimCount2 * IterationsPerPC);
	}
	else {
		if(CurrSim==1){
			SimCount2 = (CurrSim - 1) / IterationsPerPC;
			seed = SeedRecord[SimCount2][1] * 91 + SeedRecord[SimCount2][0] * 7927 + (CurrSim - SimCount2 * IterationsPerPC);
		}
		else{
			seed = SeedRecord[0][1] * 91 + SeedRecord[0][0] * 7927  + CurrSim;
		}
	}
	
	//seed=729;
	std::cout << seed << std::endl;
	rg.RandomInit(seed);
	
	RSApop.AssignAgeSex();
	
	RSApop.GetAllPartnerRates();
	
	RSApop.AssignBehav();
	
	SetCalibParameters();
	if (SetInitPrev1990 == 0){
		RSApop.AssignHIV();
	}
	
	if (HSVind == 1 || TPind == 1 || HDind == 1 || NGind == 1 || CTind == 1 ||
		TVind == 1 || BVind == 1 || VCind == 1 || HPVind == 1){
		RSApop.AssignSTIs();
	}
	/*if (CurrSim == samplesize){
		for (xx = 0; xx < 13; xx++){
			RSApop.GetInitHPVstage(xx);
			RSApop.SaveInitHPVstage("InitHPVstage.txt", xx);
		}
	}*/
	CurrYear = StartYear;
	
	for (int ii = 0; ii<ProjectionTerm; ii++){
		if (CurrYear == 1990 && SetInitPrev1990 == 1){
			RSApop.AssignHIV1990();
		}
		//cout << CurrYear << " " << Register.size() << endl;
		RSApop.OneYear();		
	}
	
	RSApop.UpdateAgeGroup();
	
}

void RunSimulations()
{
	int i, ic, SimCount2, SimCount3;
	
	InitialiseHPV();
	if (VaryParameters == 1 && FixedUncertainty == 1){
		ReadSTDparameters();
	}

	CurrSim = 0;
	for (i = 0; i<samplesize; i++){
		CurrSim += 1;
		if (i>0){
			Register.clear();
			Register.resize(InitPop);
			Register.reserve(MaxPop);
			TotCurrFSW = 0;
			for (ic = 0; ic < MaxCSWs; ic++){
				RSApop.CSWregister[ic] = 0;
			}
		}

		
		OneSimulation();
		cout << "Completed simulation " << CurrSim << endl;
		SimCount2 = CurrSim / IterationsPerPC;
		if (CurrSim == (SimCount2 * IterationsPerPC)){
			AggregateSims();
		}
		
		if(UpdateStart==1) {RSApop.SaveLifetimeCIN3("CIN3test.txt");}

		SimCount3 = CurrSim / 500;
		if (CurrSim == (SimCount3 * 500)){
			StoreOutputs();
		}
		
	}

	if (CurrSim != (SimCount3 * 500)){
		StoreOutputs();
	}
}

void StoreOutputs()
{
	int xx;
	if (HPVind == 1){ 
		
				
		if (UpdateStart == 1) {  
			RSApop.SaveAdultHPVstageAge("HPVstagesAgeRisk.txt"); 
			HSILprevNEG.RecordSample("HSILprevNEG.txt");
		}
		if (FixedUncertainty==0 && UpdateStart == 0){ 
			//RandomUniformHPV.RecordSample("RandomUniformHPV.txt");
			//HPVparamsLogL.RecordSample("HPVparamsLogL.txt");
				if(ParamCombs==1) {
					RSApop.SaveNewScreen("ModelCoverage.txt");
					RSApop.SaveCancerCases("CancerCases.txt");
				}
			if(CCcalib==1){
				/*CC_diag_IR.RecordSample("CC_diag_IR.txt");
				CC_ASR.RecordSample("CC_ASR.txt");
				CC_diag_ASR.RecordSample("CC_diag_ASR.txt");
				CC_diag_ASR1618.RecordSample("CC_diag_ASR1618.txt");
				CC_diag_ASR_death.RecordSample("CC_diag_ASR_death.txt");*/
				
				CC_ASR2.RecordSample("CC_ASR2.txt");
				/*CC_diag_ASR2.RecordSample("CC_diag_ASR2.txt");
				CC_diag_ASR16182.RecordSample("CC_diag_ASR16182.txt");
				CC_diag_ASR_death2.RecordSample("CC_diag_ASR_death2.txt");*/
				
				/*CC_ASR3.RecordSample("CC_ASR3.txt");
				CC_diag_ASR3.RecordSample("CC_diag_ASR3.txt");
				CC_diag_ASR16183.RecordSample("CC_diag_ASR16183.txt");
				CC_diag_ASR_death3.RecordSample("CC_diag_ASR_death3.txt");
				
				CC_diag_ASR_ART.RecordSample("CC_diag_ASR_ART.txt");*/
				/*CC_20.RecordSample("CC_20.txt");
				CC_25.RecordSample("CC_25.txt");
				CC_30.RecordSample("CC_30.txt");
				CC_35.RecordSample("CC_35.txt");
				CC_40.RecordSample("CC_40.txt");
				CC_45.RecordSample("CC_45.txt");
				CC_50.RecordSample("CC_50.txt");
				CC_55.RecordSample("CC_55.txt");
				CC_60.RecordSample("CC_60.txt");
				CC_65.RecordSample("CC_65.txt");
				CC_70.RecordSample("CC_70.txt");
				CC_20diag.RecordSample("CC_20diag.txt");
				CC_25diag.RecordSample("CC_25diag.txt");
				CC_30diag.RecordSample("CC_30diag.txt");
				CC_35diag.RecordSample("CC_35diag.txt");
				CC_40diag.RecordSample("CC_40diag.txt");
				CC_45diag.RecordSample("CC_45diag.txt");
				CC_50diag.RecordSample("CC_50diag.txt");
				CC_55diag.RecordSample("CC_55diag.txt");
				CC_60diag.RecordSample("CC_60diag.txt");
				CC_65diag.RecordSample("CC_65diag.txt");
				CC_70diag.RecordSample("CC_70diag.txt");
				CC_75diag.RecordSample("CC_75diag.txt");*/
				
				/*StageI.RecordSample("StageI.txt");
				StageII.RecordSample("StageII.txt");
				StageIII.RecordSample("StageIII.txt");
				StageIV.RecordSample("StageIV.txt");
				DiagCCPost2000.RecordSample("DiagCCPost2000.txt");*/

				
			}		
		
		}
		if (FixedUncertainty == 1 && GetMacD==1){RSApop.SaveMacDprev("MacDprev.txt");}
		/*	if (FixedUncertainty == 1 && OneType==1){
			NewHPVM.RecordSample("NewHPVM.txt", WhichType);
			NewHPVF.RecordSample("NewHPVF.txt", WhichType);
		}*/
		if (FixedUncertainty == 1){
			if(CCcalib==1){
				
				//CC_diag_IR.RecordSample("CC_diag_IR.txt");
				/*CC_ASR.RecordSample("CC_ASR.txt");
				CC_diag_ASR.RecordSample("CC_diag_ASR.txt");
				CC_diag_ASR_death.RecordSample("CC_diag_ASR_death.txt");
				CC_diag_ASR1618.RecordSample("CC_diag_ASR1618.txt");
				CC_diag_ASR_ART.RecordSample("CC_diag_ASR_ART.txt");*/
				/*CC_ASR2.RecordSample("CC_ASR2.txt");
				CC_diag_ASR2.RecordSample("CC_diag_ASR2.txt");
				CC_diag_ASR_death2.RecordSample("CC_diag_ASR_death2.txt");
				CC_ASR_death2.RecordSample("CC_ASR_death2.txt");
				
				CC_ASR3.RecordSample("CC_ASR3.txt");
				CC_diag_ASR3.RecordSample("CC_diag_ASR3.txt");
				CC_diag_ASR_death3.RecordSample("CC_diag_ASR_death3.txt");
				CC_ASR_death3.RecordSample("CC_ASR_death3.txt");*/
				
				/*CC_20.RecordSample("CC_20.txt");
				CC_25.RecordSample("CC_25.txt");
				CC_30.RecordSample("CC_30.txt");
				CC_35.RecordSample("CC_35.txt");
				CC_40.RecordSample("CC_40.txt");
				CC_45.RecordSample("CC_45.txt");
				CC_50.RecordSample("CC_50.txt");
				CC_55.RecordSample("CC_55.txt");
				CC_60.RecordSample("CC_60.txt");
				CC_65.RecordSample("CC_65.txt");
				CC_70.RecordSample("CC_70.txt");
				CC_20diag.RecordSample("CC_20diag.txt");
				CC_25diag.RecordSample("CC_25diag.txt");
				CC_30diag.RecordSample("CC_30diag.txt");
				CC_35diag.RecordSample("CC_35diag.txt");
				CC_40diag.RecordSample("CC_40diag.txt");
				CC_45diag.RecordSample("CC_45diag.txt");
				CC_50diag.RecordSample("CC_50diag.txt");
				CC_55diag.RecordSample("CC_55diag.txt");
				CC_60diag.RecordSample("CC_60diag.txt");
				CC_65diag.RecordSample("CC_65diag.txt");
				CC_70diag.RecordSample("CC_70diag.txt");
				CC_75diag.RecordSample("CC_75diag.txt");*/
				
				//StageI.RecordSample("StageI.txt");
				//StageII.RecordSample("StageII.txt");
				//StageIII.RecordSample("StageIII.txt");
				//StageIV.RecordSample("StageIV.txt");
				/*
			
				DiagCCPost2000.RecordSample("DiagCCPost2000.txt");*/

				/*CCprevALL.RecordSample("CCprevALL.txt");
				CCprevNEG.RecordSample("CCprevNEG.txt");
				CCprevPOS.RecordSample("CCprevPOS.txt");
				CCprevNOART.RecordSample("CCprevNOART.txt");
				CCprevART.RecordSample("CCprevART.txt");*/
				if(ParamCombs==1) {
					RSApop.SavePopPyramid("PopulationPyramid.txt");	
					RSApop.SavePopPyramid9("PopulationPyramid9.txt");	
					RSApop.SaveNewScreen("ModelCoverage.txt");
					RSApop.SaveWeeksInStage("StageWeeks.txt");
					RSApop.SaveStagediag("Stagediag.txt");
					
					RSApop.SaveCancerCases("CancerCases.txt");
					RSApop.SaveCancerDeaths("CancerDeaths.txt");
					
					//GetReferred.RecordSample("GetReferred.txt");
					//GetTreatment.RecordSample("GetTreatment.txt");
					for(xx=0; xx<13; xx++){
						NewCC[xx].RecordSample("NewCC.txt",xx);
						/*NewHPV[xx].RecordSample("NewHPV.txt",xx);
						*/
					}
				}
			}
			
			if(CCcalib==0){
				/*ABprev18to65art.RecordSample("ABprev18to65art.txt");
				HSIL_ABprev18to65art.RecordSample("HSIL_ABprev18to65art.txt");
				HSILprev18to65art.RecordSample("HSILprev18to65art.txt");
				ABprev30to65neg.RecordSample("ABprev30to65neg.txt");
				ABprev30to65pos.RecordSample("ABprev30to65pos.txt");
				ABprev30to65art.RecordSample("ABprev30to65art.txt");
				HSIL_ABprev30to65neg.RecordSample("HSIL_ABprev30to65neg.txt");
				HSIL_ABprev30to65pos.RecordSample("HSIL_ABprev30to65pos.txt");
				HSIL_ABprev30to65art.RecordSample("HSIL_ABprev30to65art.txt");*/
				/*HPVprev15to65allF.RecordSample("HPVprev15to65allF.txt");
				HPVprev15to65allM.RecordSample("HPVprev15to65allM.txt");
				HPVprev15to65posF.RecordSample("HPVprev15to65posF.txt");
				HPVprev15to65posM.RecordSample("HPVprev15to65posM.txt");
				HPVprev15to65negF.RecordSample("HPVprev15to65negF.txt");
				HPVprev15to65negM.RecordSample("HPVprev15to65negM.txt");
				HPVprev15to65noartF.RecordSample("HPVprev15to65noartF.txt");
				HPVprev15to65noartM.RecordSample("HPVprev15to65noartM.txt");
				HPVprev15to65artF.RecordSample("HPVprev15to65artF.txt");
				HPVprev15to65artM.RecordSample("HPVprev15to65artM.txt");
				HSILprevALL.RecordSample("HSILprevALL.txt");
				HSILprevNEG.RecordSample("HSILprevNEG.txt");
				HSILprevPOS.RecordSample("HSILprevPOS.txt");
				HSILprevNOART.RecordSample("HSILprevNOART.txt");
				HSILprevART.RecordSample("HSILprevART.txt");
				NewCIN2neg.RecordSample("NewCIN2neg.txt");
				NewCIN2pos.RecordSample("NewCIN2pos.txt");
				NewCIN2art.RecordSample("NewCIN2art.txt");*/
				if(ParamCombs==1) {
					RSApop.SavePopPyramid("PopulationPyramid.txt");	
					RSApop.SaveNewScreen("ModelCoverage.txt");
				}
			
			}
		}
	}
	if(HIVind==1){ // && FixedUncertainty==1){
	//	HIVprev15to49F.RecordSample("HIVprev15plus.txt");
		/*HIVprev15to49M.RecordSample("HIVprev15to49M.txt");
		HIVprev15to49all.RecordSample("HIVprev15to49all.txt");*/
		//ARTcov15to49F.RecordSample("ARTcov15plusF.txt");
		//ARTcov15to49M.RecordSample("ARTcov15plusM.txt");
		/*HIVprev15to49F.RecordSample("HIVprev15to49F.txt");
		RSApop.SaveCancerCases("CancerCases.txt");
		RSApop.SavePopPyramid("PopulationPyramid.txt");	
		PrevHH2005.RecordSample("PrevHH2005.txt");
		PrevHH2008.RecordSample("PrevHH2008.txt");
		PrevHH2012.RecordSample("PrevHH2012.txt");
		PrevHH2017.RecordSample("PrevHH2017.txt");*/
	}	
}

void AggregateSims()
{
	int SimCount2, ix,  iy, NewDiag[136], TotPop[136];

	SimCount2 = (CurrSim - 1) / IterationsPerPC;
	//WHO 2001
	double WORLD[18] = { 8857.0,8687.0,8597.0,8467.0,8217.0,7927.0,7607.0,7148.0,6588.0,6038.0,5368.0,4548.0,3719.0,2959.0,2209.0,1519.0,910.0,635.0 };
	//Canfell
	double WORLD2[18] = { 8895.0,8508.0,8082.0,7850.0,7974.0,8191.0,7444.0,6756.0,6565.0,6198.0,5510.0,4701.0,4115.0,3092.0,2249.0,1763.0,1154.0,954.0 };
	//SEGI 1966
	double WORLD3[18] = { 12000.0,10000.0,9000.0,9000.0,8000.0,8000.0,6000.0,6000.0,6000.0,6000.0,5000.0,4000.0,4000.0,3000.0,2000.0,1000.0,500.0,500.0 };
	
	if(HPVcalib==1 && FixedUncertainty == 0){
		CalcTotalLogL();
	}

	/*if(HPVind==1){
		for(iy=0; iy<136; iy++){
				CC_diag_IR.out[SimCount2][iy] = 0.0;
				CC_diag_ASR.out[SimCount2][iy] = 0.0;
				CC_diag_ASR1618.out[SimCount2][iy] = 0.0;
				CC_diag_ASR_ART.out[SimCount2][iy] = 0.0;
				CC_diag_ASR_death.out[SimCount2][iy] = 0.0;
				CC_ASR.out[SimCount2][iy] = 0.0;
				CC_diag_ASR2.out[SimCount2][iy] = 0.0;
				CC_diag_ASR16182.out[SimCount2][iy] = 0.0;
				CC_diag_ASR_death2.out[SimCount2][iy] = 0.0;
				CC_ASR_death2.out[SimCount2][iy] = 0.0;
				CC_ASR2.out[SimCount2][iy] = 0.0;
				CC_diag_ASR3.out[SimCount2][iy] = 0.0;
				CC_diag_ASR16183.out[SimCount2][iy] = 0.0;
				CC_diag_ASR_death3.out[SimCount2][iy] = 0.0;
				CC_ASR_death3.out[SimCount2][iy] = 0.0;
				CC_ASR3.out[SimCount2][iy] = 0.0;
				CC_20.out[SimCount2][iy] = 0.0;
				CC_25.out[SimCount2][iy] = 0.0;
				CC_30.out[SimCount2][iy] = 0.0;
				CC_35.out[SimCount2][iy] = 0.0;
				CC_40.out[SimCount2][iy] = 0.0;
				CC_45.out[SimCount2][iy] = 0.0;
				CC_50.out[SimCount2][iy] = 0.0;
				CC_55.out[SimCount2][iy] = 0.0;
				CC_60.out[SimCount2][iy] = 0.0;
				CC_65.out[SimCount2][iy] = 0.0;
				CC_70.out[SimCount2][iy] = 0.0;
				CC_20diag.out[SimCount2][iy] = 0.0;
				CC_25diag.out[SimCount2][iy] = 0.0;
				CC_30diag.out[SimCount2][iy] = 0.0;
				CC_35diag.out[SimCount2][iy] = 0.0;
				CC_40diag.out[SimCount2][iy] = 0.0;
				CC_45diag.out[SimCount2][iy] = 0.0;
				CC_50diag.out[SimCount2][iy] = 0.0;
				CC_55diag.out[SimCount2][iy] = 0.0;
				CC_60diag.out[SimCount2][iy] = 0.0;
				CC_65diag.out[SimCount2][iy] = 0.0;
				CC_70diag.out[SimCount2][iy] = 0.0;
				CC_75diag.out[SimCount2][iy] = 0.0;
				StageI.out[SimCount2][iy] = 0.0;
				StageII.out[SimCount2][iy] = 0.0;
				StageIII.out[SimCount2][iy] = 0.0;
				StageIV.out[SimCount2][iy] = 0.0;
				NewDiag[iy] = 0;
				TotPop[iy] = 0;				
			}
		
	
		for(ix=3; ix<18; ix++){
			for(iy=0; iy<ProjectionTerm; iy++){
				CC_ASR.out[SimCount2][iy] += RSApop.NewCancer[ix][iy]*WORLD[ix]/RSApop.PopPyramidAll[ix][iy];
				CC_diag_ASR.out[SimCount2][iy] += (RSApop.NewDiagCancer[ix][iy]+RSApop.NewDiagCancer[ix+18][iy]+RSApop.NewDiagCancer[ix+36][iy])*WORLD[ix]/RSApop.PopPyramidAll[ix][iy];
				CC_diag_ASR1618.out[SimCount2][iy] += RSApop.NewDiagCancer1618[ix][iy]*WORLD[ix]/RSApop.PopPyramidAll[ix][iy];
				CC_diag_ASR_death.out[SimCount2][iy] += RSApop.NewDiagCancerDeath[ix][iy]*WORLD[ix]/RSApop.PopPyramidAll[ix][iy];

				CC_ASR2.out[SimCount2][iy] += RSApop.NewCancer[ix][iy]*WORLD2[ix]/RSApop.PopPyramidAll[ix][iy];
				CC_diag_ASR2.out[SimCount2][iy] += (RSApop.NewDiagCancer[ix][iy]+RSApop.NewDiagCancer[ix+18][iy]+RSApop.NewDiagCancer[ix+36][iy])*WORLD2[ix]/RSApop.PopPyramidAll[ix][iy];
				CC_diag_ASR16182.out[SimCount2][iy] += RSApop.NewDiagCancer1618[ix][iy]*WORLD2[ix]/RSApop.PopPyramidAll[ix][iy];
				CC_diag_ASR_death2.out[SimCount2][iy] += RSApop.NewDiagCancerDeath[ix][iy]*WORLD2[ix]/RSApop.PopPyramidAll[ix][iy];
				CC_ASR_death2.out[SimCount2][iy] += RSApop.NewCancerDeath[ix][iy]*WORLD2[ix]/RSApop.PopPyramidAll[ix][iy];

				CC_ASR3.out[SimCount2][iy] += RSApop.NewCancer[ix][iy]*WORLD3[ix]/RSApop.PopPyramidAll[ix][iy];
				CC_diag_ASR3.out[SimCount2][iy] += (RSApop.NewDiagCancer[ix][iy]+RSApop.NewDiagCancer[ix+18][iy]+RSApop.NewDiagCancer[ix+36][iy])*WORLD3[ix]/RSApop.PopPyramidAll[ix][iy];
				CC_diag_ASR16183.out[SimCount2][iy] += RSApop.NewDiagCancer1618[ix][iy]*WORLD3[ix]/RSApop.PopPyramidAll[ix][iy];
				CC_diag_ASR_death3.out[SimCount2][iy] += RSApop.NewDiagCancerDeath[ix][iy]*WORLD3[ix]/RSApop.PopPyramidAll[ix][iy];
				CC_ASR_death3.out[SimCount2][iy] += RSApop.NewCancerDeath[ix][iy]*WORLD3[ix]/RSApop.PopPyramidAll[ix][iy];

			}
		}
		for(iy=2005-StartYear; iy<ProjectionTerm; iy++){
			CC_diag_ASR_ART.out[SimCount2][iy] += RSApop.NewDiagCancerART[iy]*100000.0/RSApop.PopPyramidAllART[iy];			
		}
		for(ix=3; ix<18; ix++){
			for(iy=0; iy<ProjectionTerm; iy++){
				NewDiag[iy] += (RSApop.NewDiagCancer[ix][iy]+RSApop.NewDiagCancer[ix+18][iy]+RSApop.NewDiagCancer[ix+36][iy]);
				TotPop[iy] += RSApop.PopPyramidAll[ix][iy];
			}
		}
		for(iy=0; iy<ProjectionTerm; iy++){
			CC_diag_IR.out[SimCount2][iy] = NewDiag[iy]*100000.0/TotPop[iy];
			
			StageI.out[SimCount2][iy] = RSApop.StageDiag[0][iy];
			StageII.out[SimCount2][iy] = RSApop.StageDiag[1][iy];
			StageIII.out[SimCount2][iy] = RSApop.StageDiag[2][iy];
			StageIV.out[SimCount2][iy] = RSApop.StageDiag[3][iy];
			
			CC_20.out[SimCount2][iy] = RSApop.NewCancer[4][iy]*100000.0/RSApop.PopPyramidAll[4][iy];
			CC_25.out[SimCount2][iy] = RSApop.NewCancer[5][iy]*100000.0/RSApop.PopPyramidAll[5][iy];
			CC_30.out[SimCount2][iy] = RSApop.NewCancer[6][iy]*100000.0/RSApop.PopPyramidAll[6][iy];
			CC_35.out[SimCount2][iy] = RSApop.NewCancer[7][iy]*100000.0/RSApop.PopPyramidAll[7][iy];
			CC_40.out[SimCount2][iy] = RSApop.NewCancer[8][iy]*100000.0/RSApop.PopPyramidAll[8][iy];
			CC_45.out[SimCount2][iy] = RSApop.NewCancer[9][iy]*100000.0/RSApop.PopPyramidAll[9][iy];
			CC_50.out[SimCount2][iy] = RSApop.NewCancer[10][iy]*100000.0/RSApop.PopPyramidAll[10][iy];
			CC_55.out[SimCount2][iy] = RSApop.NewCancer[11][iy]*100000.0/RSApop.PopPyramidAll[11][iy];
			CC_60.out[SimCount2][iy] = RSApop.NewCancer[12][iy]*100000.0/RSApop.PopPyramidAll[12][iy];
			CC_65.out[SimCount2][iy] = RSApop.NewCancer[13][iy]*100000.0/RSApop.PopPyramidAll[13][iy];
			CC_70.out[SimCount2][iy] = RSApop.NewCancer[14][iy]*100000.0/RSApop.PopPyramidAll[14][iy];
			CC_20diag.out[SimCount2][iy] = (RSApop.NewDiagCancer[4][iy]+RSApop.NewDiagCancer[4+18][iy]+RSApop.NewDiagCancer[4+36][iy])*100000.0/RSApop.PopPyramidAll[4][iy];
			CC_25diag.out[SimCount2][iy] = (RSApop.NewDiagCancer[5][iy]+RSApop.NewDiagCancer[5+18][iy]+RSApop.NewDiagCancer[5+36][iy])*100000.0/RSApop.PopPyramidAll[5][iy];
			CC_30diag.out[SimCount2][iy] = (RSApop.NewDiagCancer[6][iy]+RSApop.NewDiagCancer[6+18][iy]+RSApop.NewDiagCancer[6+36][iy])*100000.0/RSApop.PopPyramidAll[6][iy];
			CC_35diag.out[SimCount2][iy] = (RSApop.NewDiagCancer[7][iy]+RSApop.NewDiagCancer[7+18][iy]+RSApop.NewDiagCancer[7+36][iy])*100000.0/RSApop.PopPyramidAll[7][iy];
			CC_40diag.out[SimCount2][iy] = (RSApop.NewDiagCancer[8][iy]+RSApop.NewDiagCancer[8+18][iy]+RSApop.NewDiagCancer[8+36][iy])*100000.0/RSApop.PopPyramidAll[8][iy];
			CC_45diag.out[SimCount2][iy] = (RSApop.NewDiagCancer[9][iy]+RSApop.NewDiagCancer[9+18][iy]+RSApop.NewDiagCancer[9+36][iy])*100000.0/RSApop.PopPyramidAll[9][iy];
			CC_50diag.out[SimCount2][iy] = (RSApop.NewDiagCancer[10][iy]+RSApop.NewDiagCancer[10+18][iy]+RSApop.NewDiagCancer[10+36][iy])*100000.0/RSApop.PopPyramidAll[10][iy];
			CC_55diag.out[SimCount2][iy] = (RSApop.NewDiagCancer[11][iy]+RSApop.NewDiagCancer[11+18][iy]+RSApop.NewDiagCancer[11+36][iy])*100000.0/RSApop.PopPyramidAll[11][iy];
			CC_60diag.out[SimCount2][iy] = (RSApop.NewDiagCancer[12][iy]+RSApop.NewDiagCancer[12+18][iy]+RSApop.NewDiagCancer[12+36][iy])*100000.0/RSApop.PopPyramidAll[12][iy];
			CC_65diag.out[SimCount2][iy] = (RSApop.NewDiagCancer[13][iy]+RSApop.NewDiagCancer[13+18][iy]+RSApop.NewDiagCancer[13+36][iy])*100000.0/RSApop.PopPyramidAll[13][iy];
			CC_70diag.out[SimCount2][iy] = (RSApop.NewDiagCancer[14][iy]+RSApop.NewDiagCancer[14+18][iy]+RSApop.NewDiagCancer[14+36][iy])*100000.0/RSApop.PopPyramidAll[14][iy];
			CC_75diag.out[SimCount2][iy] = (RSApop.NewDiagCancer[15][iy]+RSApop.NewDiagCancer[15+18][iy]+RSApop.NewDiagCancer[15+36][iy]+
											RSApop.NewDiagCancer[16][iy]+RSApop.NewDiagCancer[16+18][iy]+RSApop.NewDiagCancer[16+36][iy]+
											RSApop.NewDiagCancer[17][iy]+RSApop.NewDiagCancer[17+18][iy]+RSApop.NewDiagCancer[17+36][iy])*100000.0/(RSApop.PopPyramidAll[15][iy]+RSApop.PopPyramidAll[16][iy]+RSApop.PopPyramidAll[17][iy]);
		}
		
		if(ParamCombs > 1){

			for(ix=3; ix<18; ix++){
				for(iy=0; iy<136; iy++){
					RSApop.NewCancer[ix][iy]=0;
					RSApop.NewDiagCancer1618[ix][iy]=0;
					RSApop.NewDiagCancerDeath[ix][iy]=0;
					RSApop.NewCancerDeath[ix][iy]=0;
					RSApop.PopPyramidAll[ix][iy]=0;
				}
			}
			for(ix=3; ix<54; ix++){
				for(iy=0; iy<136; iy++){
					RSApop.NewDiagCancer[ix][iy]=0;
				}
			}
			for(iy=0; iy<136; iy++){
				RSApop.PopPyramidAllART[iy]=0;
				RSApop.NewDiagCancerART[iy]=0;
				RSApop.StageDiag[0][iy]=0;
				RSApop.StageDiag[1][iy]=0;
				RSApop.StageDiag[2][iy]=0;
				RSApop.StageDiag[3][iy]=0;
				
			}
		}		
	}
	*/
	/*if(HPVind==1){ 
		for(iy=0; iy<136; iy++){
			
			ABprev18to65art.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotAB18to65art[iy]/HPVtransitionCC.TotPop18to65art[iy];
			HSIL_ABprev18to65art.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotHSIL18to65art[iy]/HPVtransitionCC.TotAB18to65art[iy];
			HSILprev18to65art.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotHSIL18to65art[iy]/HPVtransitionCC.TotPop18to65art[iy];
			HSILprevALL.out[SimCount2][iy] = 1.0 *(HPVtransitionCC.TotHSILNEG[iy]+HPVtransitionCC.TotHSILPOS[iy])/
											(HPVtransitionCC.TotPop15to65negF[iy]+HPVtransitionCC.TotPop15to65posF[iy]);
			HSILprevNEG.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotHSILNEG[iy]/HPVtransitionCC.TotPop15to65negF[iy];
			HSILprevPOS.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotHSILPOS[iy]/HPVtransitionCC.TotPop15to65posF[iy];
			HSILprevNOART.out[SimCount2][iy] = 1.0 *(HPVtransitionCC.TotHSILPOS[iy]-HPVtransitionCC.TotHSILART[iy])/
												(HPVtransitionCC.TotPop15to65posF[iy]-HPVtransitionCC.TotPop15to65artF[iy]);
			HSILprevART.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotHSILART[iy]/HPVtransitionCC.TotPop15to65artF[iy];
			CCprevALL.out[SimCount2][iy] = 1.0 *(HPVtransitionCC.TotCCNEG[iy]+HPVtransitionCC.TotCCPOS[iy])/
												(HPVtransitionCC.TotPop15upnegF[iy]+HPVtransitionCC.TotPop15upposF[iy]);
			CCprevNEG.out[SimCount2][iy] = 1.0 *(HPVtransitionCC.TotCCNEG[iy])/(HPVtransitionCC.TotPop15upnegF[iy]);
			CCprevPOS.out[SimCount2][iy] = 1.0 *(HPVtransitionCC.TotCCPOS[iy])/(HPVtransitionCC.TotPop15upposF[iy]);
			CCprevNOART.out[SimCount2][iy] = 1.0 *(HPVtransitionCC.TotCCPOS[iy]-HPVtransitionCC.TotCCART[iy])/
												(HPVtransitionCC.TotPop15upposF[iy]-HPVtransitionCC.TotPop15upartF[iy]);
			CCprevART.out[SimCount2][iy] = 1.0 *(HPVtransitionCC.TotCCART[iy])/(HPVtransitionCC.TotPop15upartF[iy]);
		
			
			
			//Other:
			ABprev30to65neg.out[SimCount2][iy] = 1.0 * HPVtransitionCC.TotAB30to65neg[iy]/HPVtransitionCC.TotPop30to65neg[iy];
			ABprev30to65pos.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotAB30to65pos[iy]/HPVtransitionCC.TotPop30to65pos[iy];
			ABprev30to65art.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotAB30to65art[iy]/HPVtransitionCC.TotPop30to65art[iy];
			HSIL_ABprev30to65neg.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotHSIL30to65neg[iy]/ HPVtransitionCC.TotAB30to65neg[iy];
			HSIL_ABprev30to65pos.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotHSIL30to65pos[iy]/HPVtransitionCC.TotAB30to65pos[iy];
			HSIL_ABprev30to65art.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotHSIL30to65art[iy]/HPVtransitionCC.TotAB30to65art[iy];
			ABprev18to60noart.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotAB18to60noart[iy]/HPVtransitionCC.TotPop18to60noart[iy];
			HSIL_ABprev18to60noart.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotHSIL18to60noart[iy]/HPVtransitionCC.TotAB18to60noart[iy];
			
			HPVprev15to65allF.out[SimCount2][iy] = 1.0 *(HPVtransitionCC.TotHPV15to65negF[iy]+HPVtransitionCC.TotHPV15to65posF[iy])/
												(HPVtransitionCC.TotPop15to65negF[iy]+HPVtransitionCC.TotPop15to65posF[iy]);
			HPVprev15to65allM.out[SimCount2][iy] = 1.0 *(HPVtransitionCC.TotHPV15to65negM[iy]+HPVtransitionCC.TotHPV15to65posM[iy])/
												(HPVtransitionCC.TotPop15to65negM[iy]+HPVtransitionCC.TotPop15to65posM[iy]);
			
			HPVprev15to65negF.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotHPV15to65negF[iy]/HPVtransitionCC.TotPop15to65negF[iy];
			HPVprev15to65negM.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotHPV15to65negM[iy]/HPVtransitionCC.TotPop15to65negM[iy];
			
			HPVprev15to65posF.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotHPV15to65posF[iy]/HPVtransitionCC.TotPop15to65posF[iy];
			HPVprev15to65posM.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotHPV15to65posM[iy]/HPVtransitionCC.TotPop15to65posM[iy];
			
			HPVprev15to65noartF.out[SimCount2][iy] = 1.0 *(HPVtransitionCC.TotHPV15to65posF[iy]-HPVtransitionCC.TotHPV15to65artF[iy])/
													(HPVtransitionCC.TotPop15to65posF[iy]-HPVtransitionCC.TotPop15to65artF[iy]);
			HPVprev15to65noartM.out[SimCount2][iy] = 1.0 *(HPVtransitionCC.TotHPV15to65posM[iy]-HPVtransitionCC.TotHPV15to65artM[iy])/
													(HPVtransitionCC.TotPop15to65posM[iy]-HPVtransitionCC.TotPop15to65artM[iy]);
			
			HPVprev15to65artF.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotHPV15to65artF[iy]/HPVtransitionCC.TotPop15to65artF[iy];
			HPVprev15to65artM.out[SimCount2][iy] = 1.0 *HPVtransitionCC.TotHPV15to65artM[iy]/HPVtransitionCC.TotPop15to65artM[iy];
			

			HPVtransitionCC.TotPop30to65neg[iy]=0; HPVtransitionCC.TotPop30to65pos[iy]=0; 
			HPVtransitionCC.TotPop30to65art[iy]=0; HPVtransitionCC.TotPop18to60noart[iy]=0 ;
			HPVtransitionCC.TotPop18to65art[iy]=0 ; 
			
			HPVtransitionCC.TotPop15to65negF[iy]=0; HPVtransitionCC.TotPop15to65negM[iy]=0 ;
			HPVtransitionCC.TotPop15to65posF[iy]=0; HPVtransitionCC.TotPop15to65posM[iy]=0 ;
			HPVtransitionCC.TotPop15to65artF[iy]=0; HPVtransitionCC.TotPop15to65artM[iy]=0 ;
			HPVtransitionCC.TotHPV15to65negF[iy]=0; HPVtransitionCC.TotHPV15to65negM[iy]=0 ;
			HPVtransitionCC.TotHPV15to65posF[iy]=0; HPVtransitionCC.TotHPV15to65posM[iy]=0 ;
			HPVtransitionCC.TotHPV15to65artF[iy]=0; HPVtransitionCC.TotHPV15to65artM[iy]=0 ;
			
			//Other
			HPVtransitionCC.TotAB30to65neg[iy]=0; HPVtransitionCC.TotAB30to65pos[iy]=0; 
			HPVtransitionCC.TotAB30to65art[iy]=0; HPVtransitionCC.TotAB18to60noart[iy]=0;
			HPVtransitionCC.TotHSIL30to65neg[iy]=0; HPVtransitionCC.TotHSIL30to65pos[iy]=0; 
			HPVtransitionCC.TotHSIL30to65art[iy]=0; HPVtransitionCC.TotHSIL18to60noart[iy]=0;
			HPVtransitionCC.TotAB18to65art[iy]=0 ;HPVtransitionCC.TotHSIL18to65art[iy]=0;	
			HPVtransitionCC.TotHSILNEG[iy]=0;	HPVtransitionCC.TotHSILPOS[iy]=0;	
			HPVtransitionCC.TotHSILART[iy]=0;			

			HPVtransitionCC.TotCCNEG[iy]=0;	HPVtransitionCC.TotCCPOS[iy]=0;
			HPVtransitionCC.TotCCART[iy]=0; HPVtransitionCC.TotPop15upnegF[iy]=0; 
			HPVtransitionCC.TotPop15upposF[iy]=0; HPVtransitionCC.TotPop15upartF[iy]=0;

		}
	}*/
	/*if (HIVcalib == 1){
		HIVtransitionM.GetPrev();
		HIVtransitionF.GetPrev();
	}
	if (HSVcalib == 1){ HSVtransitionF.GetCSWprev(); }
	if (TPcalib == 1){ TPtransitionF.GetCSWprev(); }
	if (NGcalib == 1){ NGtransitionF.GetCSWprev(); }
	if (CTcalib == 1){ CTtransitionF.GetCSWprev(); }
	if (TVcalib == 1){ TVtransitionF.GetCSWprev(); }
	if (FixedUncertainty == 1){
		for (iy = 3; iy < 27; iy++){
			if (HSVcalib == 1){
				//HSVprevCSW.out[CurrSim - 1][iy] = HSVtransitionF.CSWprevUnsmoothed[iy];
				HSVprevCSW.out[CurrSim - 1][iy] = 0.05 * HSVtransitionF.CSWprevUnsmoothed[iy - 3] +
					0.12 * HSVtransitionF.CSWprevUnsmoothed[iy - 2] + 0.20 * HSVtransitionF.CSWprevUnsmoothed[iy - 1] +
					0.26 * HSVtransitionF.CSWprevUnsmoothed[iy] + 0.20 * HSVtransitionF.CSWprevUnsmoothed[iy + 1] +
					0.12 * HSVtransitionF.CSWprevUnsmoothed[iy + 2] + 0.05 * HSVtransitionF.CSWprevUnsmoothed[iy + 3];
			}
			if (TPcalib == 1){
				TPprevCSW.out[CurrSim - 1][iy] = 0.05 * TPtransitionF.CSWprevUnsmoothed[iy - 3] +
					0.12 * TPtransitionF.CSWprevUnsmoothed[iy - 2] + 0.20 * TPtransitionF.CSWprevUnsmoothed[iy - 1] +
					0.26 * TPtransitionF.CSWprevUnsmoothed[iy] + 0.20 * TPtransitionF.CSWprevUnsmoothed[iy + 1] +
					0.12 * TPtransitionF.CSWprevUnsmoothed[iy + 2] + 0.05 * TPtransitionF.CSWprevUnsmoothed[iy + 3];
			}
			if (NGcalib == 1){
				//NGprevCSW.out[CurrSim - 1][iy] = NGtransitionF.CSWprevUnsmoothed[iy];
				NGprevCSW.out[CurrSim - 1][iy] = 0.05 * NGtransitionF.CSWprevUnsmoothed[iy - 3] +
					0.12 * NGtransitionF.CSWprevUnsmoothed[iy - 2] + 0.20 * NGtransitionF.CSWprevUnsmoothed[iy - 1] +
					0.26 * NGtransitionF.CSWprevUnsmoothed[iy] + 0.20 * NGtransitionF.CSWprevUnsmoothed[iy + 1] +
					0.12 * NGtransitionF.CSWprevUnsmoothed[iy + 2] + 0.05 * NGtransitionF.CSWprevUnsmoothed[iy + 3];
			}
			if (CTcalib == 1){
				//CTprevCSW.out[CurrSim - 1][iy] = CTtransitionF.CSWprevUnsmoothed[iy];
				CTprevCSW.out[CurrSim - 1][iy] = 0.05 * CTtransitionF.CSWprevUnsmoothed[iy - 3] +
					0.12 * CTtransitionF.CSWprevUnsmoothed[iy - 2] + 0.20 * CTtransitionF.CSWprevUnsmoothed[iy - 1] +
					0.26 * CTtransitionF.CSWprevUnsmoothed[iy] + 0.20 * CTtransitionF.CSWprevUnsmoothed[iy + 1] +
					0.12 * CTtransitionF.CSWprevUnsmoothed[iy + 2] + 0.05 * CTtransitionF.CSWprevUnsmoothed[iy + 3];
			}
			if (TVcalib == 1){
				//TVprevCSW.out[CurrSim - 1][iy] = TVtransitionF.CSWprevUnsmoothed[iy];
				TVprevCSW.out[CurrSim - 1][iy] = 0.05 * TVtransitionF.CSWprevUnsmoothed[iy - 3] +
					0.12 * TVtransitionF.CSWprevUnsmoothed[iy - 2] + 0.20 * TVtransitionF.CSWprevUnsmoothed[iy - 1] +
					0.26 * TVtransitionF.CSWprevUnsmoothed[iy] + 0.20 * TVtransitionF.CSWprevUnsmoothed[iy + 1] +
					0.12 * TVtransitionF.CSWprevUnsmoothed[iy + 2] + 0.05 * TVtransitionF.CSWprevUnsmoothed[iy + 3];
			}
		}
	}*/
	/*if (HIVcalib == 1){
		HIVparamsLogL.out[SimCount2][8] = TotalLogL;
	}*/
}

void SimulateParameters()
{
	/*if (HIVcalib == 1){ SimulateHIVparams(); }
	if (TPcalib == 1){ SimulateTPparameters(); }
	if (HSVcalib == 1){ SimulateHSVparameters(); }
	if (NGcalib == 1){ SimulateNGparameters(); }
	if (CTcalib == 1){ SimulateCTparameters(); }
	if (TVcalib == 1){ SimulateTVparameters(); }
	if (BVcalib == 1){ SimulateBVparameters(); }
	if (VCcalib == 1){ SimulateVCparameters(); }
	if (HPVcalib == 1){ SimulateHPVparamsType(); }*/
	if (HIVind == 1){ SimulateHIVparams(); }
	if (TPind == 1){ SimulateTPparameters(); }
	if (HSVind == 1){ SimulateHSVparameters(); }
	if (NGind == 1){ SimulateNGparameters(); }
	if (CTind == 1){ SimulateCTparameters(); }
	if (TVind == 1){ SimulateTVparameters(); }
	if (BVind == 1){ SimulateBVparameters(); }
	if (VCind == 1){ SimulateVCparameters(); }
	if (HPVind == 1){ SimulateHPVparamsType();}
	
	
}

void SimulateTPparameters()
{
	// Simulate M->F transmission prob
	TPtransitionM.TransmProb = 0.228;
	TPtransitionM.RelTransmCSW = 1.0;
	TPtransitionM.TransmProbSW = 0.228;
	// Simulate F->M transmission prob
	TPtransitionF.TransmProb = 0.155;
	// Simulate the average duration of primary syphilis (mean 6.6, std dev 2)
	TPtransitionM.AveDuration[1] = 6.59;
	TPtransitionF.AveDuration[1] = 6.59;
	// Simulate the average duration of 2ndary syphilis (mean 15.6, std dev 4)
	TPtransitionM.AveDuration[2] = 15.1;
	TPtransitionF.AveDuration[2] = 15.1;
	// Simulate the average duration of latent syphilis (mean 520, std dev 150)
	TPtransitionM.AveDuration[3] = 595;
	TPtransitionF.AveDuration[3] = 595;
	// Simulate the average duration of early immunity (mean 26, std dev 8)
	TPtransitionM.AveDuration[4] = 23.0;
	TPtransitionF.AveDuration[4] = 23.0;
	// Simulate the average duration of late immunity (mean 52, std dev 16)
	TPtransitionM.AveDuration[5] = 53.8;
	TPtransitionF.AveDuration[5] = 53.8;
	// Simulate the propn of cases correctly treated prior to introduction of SM
	TPtransitionF.CorrectRxPreSM = 0.711;
	TPtransitionM.CorrectRxPreSM = 0.711;
	// Simulate the propn of primary syphilis cases that are RPR-neg after Rx
	TPtransitionF.PropnSuscepAfterRx = 0.398;
	TPtransitionM.PropnSuscepAfterRx = 0.398;
	SecondaryRxMult = 1 - 0.682;
	
}

void SimulateHSVparameters()
{
	// Simulate M->F transmission prob in ST relationships (mean 0.0095, SD 0.0038)
	HSVtransitionM.TransmProb = 0.016;
	// Simulate F->M transmission prob in ST relationships (mean 0.0065, SD 0.0026)
	HSVtransitionF.TransmProb = 0.0035;
	// Simulate M->F transmission prob in LT relationships (mean 0.0009, SD 0.00036)
	HSVtransitionM.RelTransmLT = 0.0009 / HSVtransitionM.TransmProb;
	// Simulate F->M transmission prob in LT relationships (mean 0.00015, SD 0.00006)
	HSVtransitionF.RelTransmLT = 0.0001 / HSVtransitionF.TransmProb;
	// Simulate M->F transmission prob in CSW-client relationships (mean 0.002, SD 0.0005)
	HSVtransitionM.TransmProbSW = 0.016;
	HSVtransitionM.RelTransmCSW = 0.016 / HSVtransitionM.TransmProb;
	// Simulate the symptomatic proportion (mean 0.15, SD 0.05)
	HSVtransitionM.SymptomaticPropn = 0.14;
	HSVtransitionF.SymptomaticPropn = 0.14;
	// Simulate the average male frequency of reactivation pa (mean 6, std dev 1)
	HSVtransitionM.RecurrenceRate = 5.95 / 52.0;
	// Simulate the average female frequency of reactivation pa (mean 3, std dev 0.5)
	HSVtransitionF.RecurrenceRate = 2.88 / 52.0;
	// Simulate the increase in infectiousness during symptomatic episodes (mean 15, SD 5)
	HSVsymptomInfecIncrease = 13.9;
	// Simulate the rate of transition out of early HSV stage (mean 0.1, std dev 0.02)
	HSVtransitionM.AveDuration[1] = 52.0 / 0.098;
	HSVtransitionF.AveDuration[1] = 52.0 / 0.098;
	
}

void SimulateNGparameters()
{
	// Simulate M->F transmission prob
	NGtransitionM.TransmProb = 0.459;
	NGtransitionM.RelTransmCSW = 1.0;
	NGtransitionM.TransmProbSW = 0.459;
	// Simulate F->M transmission prob
	NGtransitionF.TransmProb = 0.237;
	// Simulate % of male NG cases that become symptomatic
	NGtransitionM.SymptomaticPropn = 0.866;
	// Simulate % of female NG cases that become symptomatic
	NGtransitionF.SymptomaticPropn = 0.299;
	// Simulate the average duration of untreated NG in males
	NGtransitionM.AveDuration[0] = 34.0;
	NGtransitionM.AveDuration[1] = 34.0;
	// Simulate the average duration of untreated NG in females
	NGtransitionF.AveDuration[0] = 33.6;
	NGtransitionF.AveDuration[1] = 33.6;
	// Simulate the proportion of treated NG cases that are immune to reinfection
	NGtransitionM.PropnImmuneAfterRx = 0.401;
	NGtransitionF.PropnImmuneAfterRx = 0.401;
	// Simulate the average duration of immunity
	NGtransitionM.AveDuration[2] = 48.8;
	NGtransitionF.AveDuration[2] = 48.8;
	// Simulate the proportion of cases correctly treated prior to syndromic mngt
	NGtransitionM.CorrectRxPreSM = 0.704;
	NGtransitionF.CorrectRxPreSM = 0.704;

}

void SimulateCTparameters()
{
	// Simulate M->F transmission prob
	CTtransitionM.TransmProb = 0.162;
	CTtransitionM.RelTransmCSW = 1.0;
	CTtransitionM.TransmProbSW = 0.162;
	// Simulate F->M transmission prob
	CTtransitionF.TransmProb = 0.0975;
	// Simulate % of male CT cases that become symptomatic
	CTtransitionM.SymptomaticPropn = 0.367;
	// Simulate % of female CT cases that become symptomatic
	CTtransitionF.SymptomaticPropn = 0.116;
	// Simulate the average duration of untreated symptomatic CT
	CTtransitionM.AveDuration[0] = 15.0;
	CTtransitionF.AveDuration[0] = 15.0;
	// Simulate the average duration of untreated asymptomatic CT
	CTtransitionM.AveDuration[1] = 106.6;
	CTtransitionF.AveDuration[1] = 106.6;
	// Simulate the proportion of treated CT cases that are immune to reinfection
	CTtransitionM.PropnImmuneAfterRx = 0.732;
	CTtransitionF.PropnImmuneAfterRx = 0.732;
	// Simulate the average duration of immunity
	CTtransitionM.AveDuration[2] = 295;
	CTtransitionF.AveDuration[2] = 295;
	// Simulate the proportion of cases correctly treated prior to syndromic mngt
	CTtransitionM.CorrectRxPreSM = 0.711;
	CTtransitionF.CorrectRxPreSM = 0.711;

}

void SimulateTVparameters()
{
	// Simulate M->F transmission prob
	TVtransitionM.TransmProb = 0.192;
	TVtransitionM.RelTransmCSW = 1.0;
	TVtransitionM.TransmProbSW = 0.192;
	// Simulate F->M transmission prob
	TVtransitionF.TransmProb = 0.0386;
	// Simulate % of male TV cases that become symptomatic
	TVtransitionM.SymptomaticPropn = 0.391;
	// Simulate % of female TV cases that become symptomatic
	TVtransitionF.SymptomaticPropn = 0.275;
	// Simulate the average duration of untreated symptomatic TV in men
	TVtransitionM.AveDuration[0] = 2.08;
	// Simulate the average duration of untreated symptomatic TV in women
	TVtransitionF.AveDuration[0] = 14.6;
	// Simulate the average duration of untreated asymptomatic TV in men
	TVtransitionM.AveDuration[1] = 20.1;
	// Simulate the average duration of untreated asymptomatic TV in women
	TVtransitionF.AveDuration[1] = 249;
	// Simulate the average duration of immunity
	TVtransitionM.AveDuration[2] = 30.0;
	TVtransitionF.AveDuration[2] = 30.0;
	// Simulate propn of cases correctly treated pre-SM
	TVtransitionM.CorrectRxPreSM = 0.447;
	TVtransitionF.CorrectRxPreSM = 0.447;
	// Simulate the proportion of treated TV cases that are immune to reinfection
	TVtransitionM.PropnImmuneAfterRx = 0.741;
	TVtransitionF.PropnImmuneAfterRx = 0.741; 
	

}

void SimulateBVparameters()
{
	int ind, i;
	double x, y, a, b, p, q;
	double r[8]; // Random variables from U(0, 1)

	ind = 2;

	if (FixedUncertainty == 0){
		int seed = (CurrSim - 1) * 91 + process_num * 7927;
		CRandomMersenne rg(seed);
		for (i = 0; i<8; i++){
			r[i] = rg.Random();
			RandomUniformBV.out[CurrSim - 1][i] = r[i];
		}
	}
	else{
		for (i = 0; i<8; i++){
			r[i] = RandomUniformBV.out[CurrSim - 1][i];
		}
	}

	// Simulate % of female BV cases that become symptomatic
	a = 4.4375;
	b = 13.313;
	p = r[0];
	q = 1 - r[0];
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
	BVtransitionF.SymptomaticPropn = x;
	BVparamsLogL.out[CurrSim - 1][0] = x;

	// Simulate the weekly incidence of BV in women with intermediate flora (mean 0.1, SD 0.03)
	a = 11.111;
	b = 111.11;
	p = r[1];
	q = 1 - r[1];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	BVtransitionF.Incidence1 = x;
	BVparamsLogL.out[CurrSim - 1][1] = x;

	// Simulate the transition rate from BV to normal flora (mean 0.008, SD 0.003)
	a = 7.1111;
	b = 888.89;
	p = r[2];
	q = 1 - r[2];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	BVtransitionF.CtsTransition[2][0] = x;
	BVtransitionF.CtsTransition[3][0] = x;
	BVparamsLogL.out[CurrSim - 1][2] = x;

	// Simulate the transition rate from BV to intermediate flora (mean 0.051, SD 0.015)
	a = 11.56;
	b = 226.67;
	p = r[3];
	q = 1 - r[3];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	BVtransitionF.CtsTransition[2][1] = x;
	BVtransitionF.CtsTransition[3][1] = x;
	BVparamsLogL.out[CurrSim - 1][3] = x;

	// Simulate the transition rate from normal to intermediate flora (mean 0.03, SD 0.01)
	a = 9.0;
	b = 300.0;
	p = r[4];
	q = 1 - r[4];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	BVtransitionF.CtsTransition[0][1] = x;
	BVparamsLogL.out[CurrSim - 1][4] = x;

	// Simulate the transition rate from intermediate to normal flora (mean 0.069, SD 0.02)
	a = 11.903;
	b = 172.5;
	p = r[5];
	q = 1 - r[5];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	BVtransitionF.CtsTransition[1][0] = x;
	BVparamsLogL.out[CurrSim - 1][5] = x;

	// Simulate the std deviation of the study effect
	/*a = 4.0;
	b = 13.333;
	p = r[6];
	q = 1 - r[6];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	BVtransitionF.HouseholdLogL.VarStudyEffect = pow(x, 2.0);
	BVtransitionF.ANClogL.VarStudyEffect = pow(x, 2.0);
	BVtransitionF.FPClogL.VarStudyEffect = pow(x, 2.0);
	BVtransitionF.CSWlogL.VarStudyEffect = pow(x, 2.0);
	BVparamsLogL.out[CurrSim-1][6] = pow(x, 2.0);*/

	// Simulate propn of cases correctly treated pre-SM
	a = 3.8667;
	b = 5.8;
	p = r[7];
	q = 1 - r[7];
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
	BVtransitionF.CorrectRxPreSM = x;
	BVparamsLogL.out[CurrSim - 1][7] = x;

	// Calculate remaining elements of the CtsTransition matrix and AveDuration for BV
	/*BVtransitionF.CtsTransition[1][2] = BVtransitionF.Incidence1 *
	BVtransitionF.SymptomaticPropn;
	BVtransitionF.CtsTransition[1][3] = BVtransitionF.Incidence1 *
	(1.0 - BVtransitionF.SymptomaticPropn);
	BVtransitionF.AveDuration[0] = 1.0/BVtransitionF.CtsTransition[0][1];
	BVtransitionF.AveDuration[1] = 1.0/(BVtransitionF.CtsTransition[1][0] +
	BVtransitionF.CtsTransition[1][2] + BVtransitionF.CtsTransition[1][3]);
	BVtransitionF.AveDuration[2] = 1.0/(BVtransitionF.CtsTransition[2][0] +
	BVtransitionF.CtsTransition[2][1]);
	BVtransitionF.AveDuration[3] = 1.0/(BVtransitionF.CtsTransition[3][0] +
	BVtransitionF.CtsTransition[3][1]);*/
}

void SimulateVCparameters()
{
	int ind, i;
	double x, y, a, b, p, q;
	double r[6]; // Random variables from U(0, 1)

	ind = 2;

	if (FixedUncertainty == 0){
		int seed = (CurrSim - 1) * 91 + process_num * 7927;
		CRandomMersenne rg(seed);
		for (i = 0; i<6; i++){
			r[i] = rg.Random();
			RandomUniformVC.out[CurrSim - 1][i] = r[i];
		}
	}
	else{
		for (i = 0; i<6; i++){
			r[i] = RandomUniformVC.out[CurrSim - 1][i];
		}
	}

	// Simulate the ave dur of symptomatic VC (mean 12, SD 3)
	a = 16.0;
	b = 1.3333;
	p = r[0];
	q = 1 - r[0];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	VCtransitionF.AveDuration[1] = x;
	VCparamsLogL.out[CurrSim - 1][0] = x;

	// Simulate the ave dur of asymptomatic VC (mean 26, SD 6)
	a = 18.778;
	b = 0.7222;
	p = r[1];
	q = 1 - r[1];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	VCtransitionF.AveDuration[0] = x;
	VCparamsLogL.out[CurrSim - 1][1] = x;

	// Simulate the incidence of Candida colonization (mean 0.8, SD 0.2)
	a = 16.0;
	b = 20.0;
	p = r[2];
	q = 1 - r[2];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	VCtransitionF.Incidence = x;
	VCparamsLogL.out[CurrSim - 1][2] = x;

	// Simulate the incidence of symptomatic VC in all asymp women (mean 0.15, SD 0.05)
	a = 9.0;
	b = 60.0;
	p = r[3];
	q = 1 - r[3];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	VCtransitionF.RecurrenceRate = x * (VCtransitionF.Incidence + 1.0 /
		VCtransitionF.AveDuration[0]) / VCtransitionF.Incidence;
	VCparamsLogL.out[CurrSim - 1][3] = x;

	// Simulate the std deviation of the study effect
	/*a = 4.0;
	b = 13.333;
	p = r[4];
	q = 1 - r[4];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	VCtransitionF.HouseholdLogL.VarStudyEffect = pow(x, 2.0);
	VCtransitionF.ANClogL.VarStudyEffect = pow(x, 2.0);
	VCtransitionF.FPClogL.VarStudyEffect = pow(x, 2.0);
	VCtransitionF.CSWlogL.VarStudyEffect = pow(x, 2.0);
	VCparamsLogL.out[CurrSim-1][4] = pow(x, 2.0);*/

	// Simulate the proportion of VC cases correctly treated (mean 0.5, SD 0.2)
	a = 2.625;
	b = 2.625;
	p = r[5];
	q = 1 - r[5];
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
	VCtransitionF.CorrectRxPreSM = x;
	VCtransitionF.CorrectRxWithSM = x;
	VCparamsLogL.out[CurrSim - 1][5] = x;
}

void SimulateHIVparams()
{
	// Simulate relative infectiousness during acute HIV
	HIVtransitionM.HIVinfecIncrease[0] = 18.3; // x - 1.0;
	HIVtransitionF.HIVinfecIncrease[0] = 18.3; //19.3; // x - 1.0;
    // Simulate relative infectiousness during advanced HIV
	HIVtransitionM.HIVinfecIncrease[2] = pow(6.90 , 0.5) - 1.0; //pow(x, 0.5) - 1.0;
	HIVtransitionF.HIVinfecIncrease[2] = pow(6.90 , 0.5) - 1.0;  //pow(x, 0.5) - 1.0;
	HIVtransitionM.HIVinfecIncrease[3] = 5.90; // x - 1.0;
	HIVtransitionF.HIVinfecIncrease[3] = 5.90; // x - 1.0;
    HIVtransitionM.HIVinfecIncrease[5] = pow(6.90 , 0.5) - 1.0; //pow(x, 0.5) - 1.0;
	HIVtransitionF.HIVinfecIncrease[5] = pow(6.90 , 0.5) - 1.0;  //pow(x, 0.5) - 1.0;
	// Simulate M->F transmission prob per sex act in short-term partnerships
	InitHIVtransm[1][0] = 0.0087;  //x; //0.0081
	// Simulate F->M transmission prob per sex act in short-term partnerships
	InitHIVtransm[1][1] = 0.0035; // x; //0.0036;
	InitHIVtransm[0][1] = 0.0035; // x; //0.0036;
	// Simulate M->F transmission prob per sex act in spousal partnerships (low risk spouse)
	InitHIVtransm[2][0] = 0.0022; // x;//0.0019
	// Simulate F->M transmission prob per sex act in spousal partnerships (low risk spouse)
	InitHIVtransm[2][1] = 0.0017; // x;
	// Simulate bias in reporting of condom use
	CondomScaling = 0.3617665;//r[6]; 	
	// Simulate initial HIV prevalence
	InitHIVprevHigh = 0.0231; //0.6568825; // (r[7] * 0.02) + 0.01;

}


///HPV - Cari///
HPVTransition::HPVTransition(){}
void HPVTransition::HPVObs(int Sex, int ObsANC, int ObsFPC, int ObsGUD,
	int ObsCSW, int ObsHH, int ObsANCN, int ObsHHN, int ObsNOART, int ObsONART) //fix
{
	SexInd = Sex;
	nStates = 8; //Leave this as 8 - will not simulate diagnosed states from 1985

	ANClogL.Observations = ObsANC;
	FPClogL.Observations = ObsFPC;
	GUDlogL.Observations = ObsGUD;
	CSWlogL.Observations = ObsCSW;
	HouseholdLogL.Observations = ObsHH;
	AntenatalNlogL.Observations = ObsANCN;
	HouseholdNlogL.Observations = ObsHHN;
	NOARTlogL.Observations = ObsNOART;
	ONARTlogL.Observations = ObsONART;
}
HPVtransitionCC::HPVtransitionCC(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH,
	int ObsANCN, int ObsHHN, int ObsNOART, int ObsONART)
{
	SexInd = Sex;
	nStates = 8; //Leave this as 8 - will not simulate diagnosed states from 1985

	ANClogL.Observations = ObsANC;
	FPClogL.Observations = ObsFPC;
	GUDlogL.Observations = ObsGUD;
	CSWlogL.Observations = ObsCSW;
	HouseholdLogL.Observations = ObsHH;
	AntenatalNlogL.Observations = ObsANCN;
	HouseholdNlogL.Observations = ObsHHN;
	NOARTlogL.Observations = ObsNOART;
	ONARTlogL.Observations = ObsONART;
}
void HPVTransition::CalcTransitionProbsF()
{
	//Progression/Regression from HPV+ if not reactivated infection before
	From1to2ind = (1 - exp(-(prog_fromHPV)* 52.0 / CycleD));
	From1to6ind = (1 - exp(-(reg_fromHPV)* 52.0 / CycleD))*propLatent;
	From1to7ind = (1 - exp(-(reg_fromHPV)* 52.0 / CycleD))*(1 - propLatent);
	
	From1to6dep = (1.0 - exp(-From1to6ind)) * (1.0 - 0.5 * (2.0 - exp(-From1to7ind) - exp(-From1to2ind)) + (1.0 - exp(-From1to7ind)) * (1.0 - exp(-From1to2ind)) / 3.0); 
	From1to7dep = (1.0 - exp(-From1to7ind)) * (1.0 - 0.5 * (2.0 - exp(-From1to6ind) - exp(-From1to2ind)) + (1.0 - exp(-From1to6ind)) * (1.0 - exp(-From1to2ind)) / 3.0); 
	From1to2dep = (1.0 - exp(-From1to2ind)) * (1.0 - 0.5 * (2.0 - exp(-From1to6ind) - exp(-From1to7ind)) + (1.0 - exp(-From1to6ind)) * (1.0 - exp(-From1to7ind)) / 3.0); 
	
	From1to2indAcute = (1 - exp(-(CIN1_HIV*prog_fromHPV)* 52.0 / CycleD));
	From1to6indAcute = (1 - exp(-(reg_HIV*reg_fromHPV)* 52.0 / CycleD))*propLatent;
	From1to7indAcute = (1 - exp(-(reg_HIV*reg_fromHPV)* 52.0 / CycleD))*(1 - propLatent);
	
	From1to6depAcute = (1.0 - exp(-From1to6indAcute)) * (1.0 - 0.5 * (2.0 - exp(-From1to7indAcute) - exp(-From1to2indAcute)) + (1.0 - exp(-From1to7indAcute)) * (1.0 - exp(-From1to2indAcute)) / 3.0); 
	From1to7depAcute = (1.0 - exp(-From1to7indAcute)) * (1.0 - 0.5 * (2.0 - exp(-From1to6indAcute) - exp(-From1to2indAcute)) + (1.0 - exp(-From1to6indAcute)) * (1.0 - exp(-From1to2indAcute)) / 3.0); 
	From1to2depAcute = (1.0 - exp(-From1to2indAcute)) * (1.0 - 0.5 * (2.0 - exp(-From1to6indAcute) - exp(-From1to7indAcute)) + (1.0 - exp(-From1to6indAcute)) * (1.0 - exp(-From1to7indAcute)) / 3.0); 
	
	From1to2indLatent = (1 - exp(-(CIN1_ART*prog_fromHPV)* 52.0 / CycleD));
	From1to6indLatent = (1 - exp(-(reg_ART*reg_fromHPV)* 52.0 / CycleD))*propLatent;
	From1to7indLatent = (1 - exp(-(reg_ART*reg_fromHPV)* 52.0 / CycleD))*(1 - propLatent);
		
	From1to6depLatent = (1.0 - exp(-From1to6indLatent)) * (1.0 - 0.5 * (2.0 - exp(-From1to7indLatent) - exp(-From1to2indLatent)) + (1.0 - exp(-From1to7indLatent)) * (1.0 - exp(-From1to2indLatent)) / 3.0); 
	From1to7depLatent = (1.0 - exp(-From1to7indLatent)) * (1.0 - 0.5 * (2.0 - exp(-From1to6indLatent) - exp(-From1to2indLatent)) + (1.0 - exp(-From1to6indLatent)) * (1.0 - exp(-From1to2indLatent)) / 3.0); 
	From1to2depLatent = (1.0 - exp(-From1to2indLatent)) * (1.0 - 0.5 * (2.0 - exp(-From1to6indLatent) - exp(-From1to7indLatent)) + (1.0 - exp(-From1to6indLatent)) * (1.0 - exp(-From1to7indLatent)) / 3.0); 
	
	//latent to infected
	From6to1 = (1.0 - exp(-(1 / AveDuration[5])* 52.0 / CycleD));
	From6to1Latent = (1.0 - exp(-latentHIVreact*(1 / AveDuration[5])* 52.0 / CycleD));
	From6to1Late = (1.0 - exp(-lateHIVreact*(1 / AveDuration[5])* 52.0 / CycleD));
	//immune to susceptible
	From7to0 = 1 - exp(-(1 / AveDuration[6])* 52.0 / CycleD);

	//Progression/regression from CIN1 <30
	From2to1ind = (1.0 - exp(-(reg_fromCIN1/CycleD)))*(1-prop_reg_cl); 
	From2to6ind = (1.0 - exp(-(reg_fromCIN1/CycleD)))*prop_reg_cl*propLatent; 
	From2to7ind = (1.0 - exp(-(reg_fromCIN1/CycleD)))*prop_reg_cl*(1-propLatent);
	From2to3ind = (1.0 - exp(-(prog_fromCIN1/CycleD)));
	
	From2to1dep = From2to1ind * (1.0 - 0.5 * (From2to6ind + From2to7ind + From2to3ind) + (From2to6ind * From2to7ind +
				From2to6ind * From2to3ind + From2to7ind * From2to3ind) / 3.0 - 0.25 * From2to6ind * From2to7ind * From2to3ind);
	From2to6dep = From2to6ind * (1.0 - 0.5 * (From2to1ind + From2to7ind + From2to3ind) + (From2to1ind * From2to7ind +
				From2to1ind * From2to3ind + From2to7ind * From2to3ind) / 3.0 - 0.25 * From2to1ind * From2to7ind * From2to3ind);
	From2to7dep = From2to7ind * (1.0 - 0.5 * (From2to1ind + From2to6ind + From2to3ind) + (From2to1ind * From2to6ind +
				From2to1ind * From2to3ind + From2to6ind * From2to3ind) / 3.0 - 0.25 * From2to1ind * From2to6ind * From2to3ind);
	From2to3dep = From2to3ind * (1.0 - 0.5 * (From2to1ind + From2to6ind + From2to7ind) + (From2to1ind * From2to6ind +
				From2to1ind * From2to7ind + From2to6ind * From2to7ind) / 3.0 - 0.25 * From2to1ind * From2to6ind * From2to7ind);
	
	From2to1ind = (1.0 - exp(-(reg_fromCIN1*reg_HIV/CycleD)))*(1-prop_reg_cl); 
	From2to6ind = (1.0 - exp(-(reg_fromCIN1*reg_HIV/CycleD)))*prop_reg_cl*propLatent; 
	From2to7ind = (1.0 - exp(-(reg_fromCIN1*reg_HIV/CycleD)))*prop_reg_cl*(1-propLatent);
	From2to3ind = (1.0 - exp(-(prog_fromCIN1*CIN2_HIV/CycleD)));
	
	From2to1depLate = From2to1ind * (1.0 - 0.5 * (From2to6ind + From2to7ind + From2to3ind) + (From2to6ind * From2to7ind +
				From2to6ind * From2to3ind + From2to7ind * From2to3ind) / 3.0 - 0.25 * From2to6ind * From2to7ind * From2to3ind);
	From2to6depLate = From2to6ind * (1.0 - 0.5 * (From2to1ind + From2to7ind + From2to3ind) + (From2to1ind * From2to7ind +
				From2to1ind * From2to3ind + From2to7ind * From2to3ind) / 3.0 - 0.25 * From2to1ind * From2to7ind * From2to3ind);
	From2to7depLate = From2to7ind * (1.0 - 0.5 * (From2to1ind + From2to6ind + From2to3ind) + (From2to1ind * From2to6ind +
				From2to1ind * From2to3ind + From2to6ind * From2to3ind) / 3.0 - 0.25 * From2to1ind * From2to6ind * From2to3ind);
	From2to3depLate = From2to3ind * (1.0 - 0.5 * (From2to1ind + From2to6ind + From2to7ind) + (From2to1ind * From2to6ind +
				From2to1ind * From2to7ind + From2to6ind * From2to7ind) / 3.0 - 0.25 * From2to1ind * From2to6ind * From2to7ind);
	
	From2to1ind = (1.0 - exp(-(reg_fromCIN1*reg_ART/CycleD)))*(1-prop_reg_cl); 
	From2to6ind = (1.0 - exp(-(reg_fromCIN1*reg_ART/CycleD)))*prop_reg_cl*propLatent; 
	From2to7ind = (1.0 - exp(-(reg_fromCIN1*reg_ART/CycleD)))*prop_reg_cl*(1-propLatent);
	From2to3ind = (1.0 - exp(-(prog_fromCIN1*CIN2_ART/CycleD)));
	
	From2to1depLatent = From2to1ind * (1.0 - 0.5 * (From2to6ind + From2to7ind + From2to3ind) + (From2to6ind * From2to7ind +
				From2to6ind * From2to3ind + From2to7ind * From2to3ind) / 3.0 - 0.25 * From2to6ind * From2to7ind * From2to3ind);
	From2to6depLatent = From2to6ind * (1.0 - 0.5 * (From2to1ind + From2to7ind + From2to3ind) + (From2to1ind * From2to7ind +
				From2to1ind * From2to3ind + From2to7ind * From2to3ind) / 3.0 - 0.25 * From2to1ind * From2to7ind * From2to3ind);
	From2to7depLatent = From2to7ind * (1.0 - 0.5 * (From2to1ind + From2to6ind + From2to3ind) + (From2to1ind * From2to6ind +
				From2to1ind * From2to3ind + From2to6ind * From2to3ind) / 3.0 - 0.25 * From2to1ind * From2to6ind * From2to3ind);
	From2to3depLatent = From2to3ind * (1.0 - 0.5 * (From2to1ind + From2to6ind + From2to7ind) + (From2to1ind * From2to6ind +
				From2to1ind * From2to7ind + From2to6ind * From2to7ind) / 3.0 - 0.25 * From2to1ind * From2to6ind * From2to7ind);
	
	//Progression/regression from CIN2 if aged<30 	 
	From3to2ind = (1.0 - exp(-(reg_fromCIN2/CycleD)));
	From3to4ind = (1.0 - exp(-(prog_fromCIN2/CycleD)));
	From3to2dep = From3to2ind * (1.0 - 0.5*(From3to4ind));
	From3to4dep = From3to4ind * (1.0 - 0.5*(From3to2ind));
	//cout << "From3to2dep " << From3to2dep << " From3to4dep " <<From3to4dep<< endl;
	
	From3to2ind = (1.0 - exp(-(reg_fromCIN2*reg_HIV/CycleD)));
	From3to4ind = (1.0 - exp(-(prog_fromCIN2*CIN2_HIV/CycleD)));
	From3to2depLate = From3to2ind * (1.0 - 0.5*(From3to4ind));
	From3to4depLate = From3to4ind * (1.0 - 0.5*(From3to2ind));
	//cout << "From3to2depLate " << From3to2depLate << " From3to4depLate " <<From3to4depLate<< endl;
	
	From3to2ind = (1.0 - exp(-(reg_fromCIN2*reg_ART/CycleD)));
	From3to4ind = (1.0 - exp(-(prog_fromCIN2*CIN2_ART/CycleD)));
	From3to2depLatent = From3to2ind * (1.0 - 0.5*(From3to4ind));
	From3to4depLatent = From3to4ind * (1.0 - 0.5*(From3to2ind));
	//cout << "From3to2depLatent " << From3to2depLatent << " From3to4depLatent " <<From3to4depLatent<< endl;
	
	//Progression/regression from CIN2 if aged>=30
	From3to2ind = (1.0 - exp(-(CIN2reg30*reg_fromCIN2/CycleD)));
	From3to4ind = (1.0 - exp(-(CIN2prog30*prog_fromCIN2/CycleD)));
	From3to2dep30 = From3to2ind * (1.0 - 0.5*(From3to4ind));
	From3to4dep30 = From3to4ind * (1.0 - 0.5*(From3to2ind));
	//cout << "From3to2dep30 " << From3to2dep30 << " From3to4dep30 " <<From3to4dep30<< endl;
	
	From3to2ind = (1.0 - exp(-(CIN2reg30*reg_fromCIN2*reg_HIV/CycleD)));
	From3to4ind = (1.0 - exp(-(CIN2prog30*prog_fromCIN2*CIN2_HIV/CycleD)));
	From3to2depLate30 = From3to2ind * (1.0 - 0.5*(From3to4ind));
	From3to4depLate30 = From3to4ind * (1.0 - 0.5*(From3to2ind));
	//cout << "From3to2depLate30 " << From3to2depLate30 << " From3to4depLate30 " <<From3to4depLate30<< endl;
	
	From3to2ind = (1.0 - exp(-(CIN2reg30*reg_fromCIN2*reg_ART/CycleD)));
	From3to4ind = (1.0 - exp(-(CIN2prog30*prog_fromCIN2*CIN2_ART/CycleD)));
	From3to2depLatent30 = From3to2ind * (1.0 - 0.5*(From3to4ind));
	From3to4depLatent30 = From3to4ind * (1.0 - 0.5*(From3to2ind));
	//cout << "From3to2depLatent30 " << From3to2depLatent30 << " From3to4depLatent30 " <<From3to4depLatent30<< endl;
	
	//Progression/regression from CIN2 if aged>=50
    From3to2ind = (1.0 - exp(-(CIN2reg30*reg_fromCIN2/CycleD)));
    From3to4ind = (1.0 - exp(-(CIN2prog50*prog_fromCIN2/CycleD)));
    From3to2dep50 = From3to2ind * (1.0 - 0.5*(From3to4ind));
    From3to4dep50 = From3to4ind * (1.0 - 0.5*(From3to2ind));
	//cout << "From3to2dep50 " << From3to2dep50 << " From3to4dep50 " <<From3to4dep50<< endl;
	
    
    From3to2ind = (1.0 - exp(-(CIN2reg30*reg_fromCIN2*reg_HIV/CycleD)));
    From3to4ind = (1.0 - exp(-(CIN2prog50*prog_fromCIN2*CIN2_HIV/CycleD)));
    From3to2depLate50 = From3to2ind * (1.0 - 0.5*(From3to4ind));
    From3to4depLate50 = From3to4ind * (1.0 - 0.5*(From3to2ind));
	//cout << "From3to2depLate50 " << From3to2depLate50 << " From3to4depLate50 " <<From3to4depLate50<< endl;
	
    From3to2ind = (1.0 - exp(-(CIN2reg30*reg_fromCIN2*reg_ART/CycleD)));
    From3to4ind = (1.0 - exp(-(CIN2prog50*prog_fromCIN2*CIN2_ART/CycleD)));
    From3to2depLatent50 = From3to2ind * (1.0 - 0.5*(From3to4ind));
    From3to4depLatent50 = From3to4ind * (1.0 - 0.5*(From3to2ind));
	//cout << "From3to2depLatent50 " << From3to2depLatent50 << " From3to4depLatent50 " <<From3to4depLatent50<< endl;
	
	 
	//if(CIN3distribution==0){
	//	From4to5 = 1 - exp(-(1 / AveDuration[3])* 52.0 / CycleD);
	//	From4to5HIV = 1 - exp(-(1 /CIN3HIV*AveDuration[3])* 52.0 / CycleD);
	//}
	
	//Progression through cancer stages
	//From5toDead = 1.0 - exp(-(1 / AveDuration[4])* 52.0 / CycleD);
	From5to8ind = (1.0 - exp(-(prog_stageI/CycleD)));
	From5to11ind = 1.0 - exp(-(-log(1.0-diag_stageI)/CycleD)); // ann prob worked to ann rate to weekly prob
	From5to8dep = From5to8ind * (1.0 - 0.5*(From5to11ind));
	From5to11dep = From5to11ind * (1.0 - 0.5*(From5to8ind));

	From5to8ind = (1.0 - exp(-(CIN2_HIV*prog_stageI/CycleD)));
	From5to8depHIV = From5to8ind * (1.0 - 0.5*(From5to11ind));
	From5to11depHIV = From5to11ind * (1.0 - 0.5*(From5to8ind));

	From5to8ind = (1.0 - exp(-(CIN2_ART*prog_stageI/CycleD)));
	From5to8depART = From5to8ind * (1.0 - 0.5*(From5to11ind));
	From5to11depART = From5to11ind * (1.0 - 0.5*(From5to8ind));

	From8to9ind = (1.0 - exp(-(prog_stageII/CycleD)));
	From8to12ind = 1.0 - exp(-(-log(1.0-diag_stageII)/CycleD)); // ann prob worked to ann rate to weekly prob
	From8to9dep = From8to9ind * (1.0 - 0.5*(From8to12ind));
	From8to12dep = From8to12ind * (1.0 - 0.5*(From8to9ind));
	
	From8to9ind = (1.0 - exp(-(CIN2_HIV*prog_stageII/CycleD)));
	From8to9depHIV = From8to9ind * (1.0 - 0.5*(From8to12ind));
	From8to12depHIV = From8to12ind * (1.0 - 0.5*(From8to9ind));

	From8to9ind = (1.0 - exp(-(CIN2_ART*prog_stageII/CycleD)));
	From8to9depART = From8to9ind * (1.0 - 0.5*(From8to12ind));
	From8to12depART = From8to12ind * (1.0 - 0.5*(From8to9ind));

	From9to10ind = (1.0 - exp(-(prog_stageIII/CycleD)));
	From9to13ind = 1.0 - exp(-(-log(1.0-diag_stageIII)/CycleD)); // ann prob worked to ann rate to weekly prob
	From9to10dep = From9to10ind * (1.0 - 0.5*(From9to13ind));
	From9to13dep = From9to13ind * (1.0 - 0.5*(From9to10ind));
	
	From9to10ind = (1.0 - exp(-(CIN2_HIV*prog_stageIII/CycleD)));
	From9to10depHIV = From9to10ind * (1.0 - 0.5*(From9to13ind));
	From9to13depHIV = From9to13ind * (1.0 - 0.5*(From9to10ind));
	
	From9to10ind = (1.0 - exp(-(CIN2_ART*prog_stageIII/CycleD)));
	From9to10depART = From9to10ind * (1.0 - 0.5*(From9to13ind));
	From9to13depART = From9to13ind * (1.0 - 0.5*(From9to10ind));
	
	From10to14 = (1.0 - exp(-(-log(1.0-diag_stageIV)/CycleD)))*(1 - 0.5*(1.0 - exp(-(1.0 / 26.0)*52.0 / CycleD)));
	From10toDead = (1.0 - exp(-(1.0 / 26.0)*52.0 / CycleD))*(1 - 0.5*(1.0 - exp(-(-log(1.0-diag_stageIV)/CycleD)))); 
		
}
void HPVTransition::CalcTransitionProbsM()
{
	From1to6ind = (1 - exp(-(1 / AveDuration[0])* 52.0 / CycleD))*propLatent;
	From1to7ind = (1 - exp(-(1 / AveDuration[0])* 52.0 / CycleD))*(1 - propLatent);
	From1to6dep = From1to6ind * (1.0 - 0.5*(From1to7ind));
	From1to7dep = From1to7ind * (1.0 - 0.5*(From1to6ind));

	From1to6indAcute = (1 - exp(-(1 / (AcuteLateHIVclear*AveDuration[0]))* 52.0 / CycleD))*propLatent;
	From1to7indAcute = (1 - exp(-(1 / (AcuteLateHIVclear*AveDuration[0]))* 52.0 / CycleD))*(1 - propLatent);
	From1to6depAcute = From1to6indAcute * (1.0 - 0.5*(From1to7indAcute));
	From1to7depAcute = From1to7indAcute * (1.0 - 0.5*(From1to6indAcute));

	From1to6indLatent = (1 - exp(-(1 / (LatentHIVclear*AveDuration[0]))* 52.0 / CycleD))*propLatent;
	From1to7indLatent = (1 - exp(-(1 / (LatentHIVclear*AveDuration[0]))* 52.0 / CycleD))*(1 - propLatent);
	From1to6depLatent = From1to6indLatent * (1.0 - 0.5*(From1to7indLatent));
	From1to7depLatent = From1to7indLatent * (1.0 - 0.5*(From1to6indLatent));

	From6to1 = (1.0 - exp(-(1 / AveDuration[5])* 52.0 / CycleD));
	From6to1Latent = (1.0 - exp(-latentHIVreact*(1 / AveDuration[5])* 52.0 / CycleD));
	From6to1Late = (1.0 - exp(-lateHIVreact*(1 / AveDuration[5])* 52.0 / CycleD));
		
	From7to0 = 1 - exp(-(1 / AveDuration[6])* 52.0 / CycleD);
}

double HPVTransition::GetTransmProb(int ID, int type)
{
	// ID is the ID of the suscpetible partner, but the HPVTransition object is for
	// the opposite sex (the sex of the infected partner).
	
	int PID1, PID2,  ic;
	int IRisk; // , PRisk;
	double NoTransmProb, SingleActProb;

	IRisk = Register[ID - 1].RiskGroup;
	NoTransmProb = 1.0;
	
	if(Register[ID - 1].VaccinationStatus[type]==1){ NoTransmProb = 1.0; }
	
	else {
		// Infection from primary partner	
		if (Register[ID - 1].CurrPartners > 0 ){
			PID1 = Register[ID - 1].IDprimary;
			if (Register[PID1 - 1].HPVstage[type] == 1 || Register[PID1 - 1].HPVstage[type] == 2 || Register[PID1 - 1].HPVstage[type] == 3 || Register[PID1 - 1].HPVstage[type] == 4){    
				// Get base probability depending on partnership type
				SingleActProb = TransmProb;
				NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVIprimary);		
				NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
					Register[ID - 1].PVIprimary);				
			}
		}

		// Infection from 2ndary partner
		if (Register[ID - 1].CurrPartners == 2){
			PID2 = Register[ID - 1].ID2ndary;
			if (Register[PID2 - 1].HPVstage[type] == 1 || Register[PID2 - 1].HPVstage[type] == 2 || Register[PID2 - 1].HPVstage[type] == 3 || Register[PID2 - 1].HPVstage[type] == 4){
				SingleActProb = TransmProb;
				
				NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVI2ndary);
				NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
					Register[ID - 1].PVI2ndary);
			
			}
		}

		// Infection from CSW (relevant only to high-risk men)
		if (Register[ID - 1].UVICSW + Register[ID - 1].PVICSW > 0){
			PID2 = Register[ID - 1].IDofCSW;
		
			if (Register[PID2 - 1].HPVstage[type] == 1 || Register[PID2 - 1].HPVstage[type] == 2 || Register[PID2 - 1].HPVstage[type] == 3 || Register[PID2 - 1].HPVstage[type] == 4){
				SingleActProb = TransmProb;
				
				// F-to-M transmission prob is assumed to be the same in commercial sex
					NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVICSW);
					NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
					Register[ID - 1].PVICSW);
			}
		}	

		// Infection from client (relevant only to CSWs)
		int tpp = Register.size();
		if (Register[ID - 1].FSWind == 1){
			for (ic = 0; ic<tpp; ic++){
				if (Register[ic].IDofCSW == ID && (Register[ic].UVICSW + Register[ic].PVICSW) > 0){
				
					if (Register[ic].HPVstage[type] == 1 || Register[ic].HPVstage[type] == 2 || Register[ic].HPVstage[type] == 3 || Register[ic].HPVstage[type] == 4){
						SingleActProb = TransmProb;
						NoTransmProb *= pow(1.0 - SingleActProb, Register[ic].UVICSW);
						NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
							Register[ic].PVICSW);
					}
				}
			}
		}
	}
	
	return 1.0 - NoTransmProb;
}

void HPVTransition::GetNewStageF(int ID, double p, int type)
{
	int is, xx, ih, SimCount2; 

	SimCount2 = (CurrSim - 1) / IterationsPerPC;

	xx = type;

	is = Register[ID - 1].HPVstage[xx];
	ih = Register[ID - 1].HIVstage;
	if (is == 1 ){
			if (ih == 0){
				if (p < From1to7dep){ Register[ID - 1].HPVstageE[xx] = 7; }
				else if (p < From1to7dep + From1to6dep){ Register[ID - 1].HPVstageE[xx] = 6; }
				else if (p < From1to7dep + From1to6dep + From1to2dep){ Register[ID - 1].HPVstageE[xx] = 2; }
				else{ Register[ID - 1].HPVstageE[xx] = 1; }
			}
			if (ih == 2 ||ih == 1 || ih==3 || ih==4 || (ih==5  && Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks<104)||ih==6){ //
				if (p < From1to7depAcute){ Register[ID - 1].HPVstageE[xx] = 7; }
				else if (p < From1to7depAcute + From1to6depAcute){ Register[ID - 1].HPVstageE[xx] = 6; }
				else if (p < From1to7depAcute + From1to6depAcute + From1to2depAcute){ Register[ID - 1].HPVstageE[xx] = 2; }
				else{ Register[ID - 1].HPVstageE[xx] = 1; }
			}
			if ( (ih==5  && Register[ID-1].ARTstage==0)||
				(ih==5  && Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks>=104)){ 
				if (p < From1to7depLatent){ Register[ID - 1].HPVstageE[xx] = 7; }
				else if (p < From1to7depLatent + From1to6depLatent){ Register[ID - 1].HPVstageE[xx] = 6; }
				else if (p < From1to7depLatent + From1to6depLatent + From1to2depLatent){ Register[ID - 1].HPVstageE[xx] = 2; }
				else{ Register[ID - 1].HPVstageE[xx] = 1; }
			}
		}

	if (is == 2 ){
			if(ih == 0 ) { //|| (ih == 5 && Register[ID - 1].ARTweeks>=104)
				if (p < From2to1dep){ Register[ID - 1].HPVstageE[xx] = 1; }
				else if (p < From2to1dep + From2to3dep){ Register[ID - 1].HPVstageE[xx] = 3; }
				else if (p < From2to1dep + From2to3dep + From2to6dep){ Register[ID - 1].HPVstageE[xx] = 6; }
				else if (p < From2to1dep + From2to3dep + From2to6dep + From2to7dep){ Register[ID - 1].HPVstageE[xx] = 7; }
				else{ Register[ID - 1].HPVstageE[xx] = 2; }
			}
			if (ih == 2 ||ih == 1 || ih==3 || ih==4 ||  (ih==5  && Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks<104)||ih==6){ //&& Register[ID-1].ARTweeks<104
				if (p < From2to1depLate){ Register[ID - 1].HPVstageE[xx] = 1; }
				else if (p < From2to1depLate + From2to3depLate){ Register[ID - 1].HPVstageE[xx] = 3; }
				else if (p < From2to1depLate + From2to3depLate + From2to6depLate){ Register[ID - 1].HPVstageE[xx] = 6; }
				else if (p < From2to1depLate + From2to3depLate + From2to6depLate + From2to7depLate){ Register[ID - 1].HPVstageE[xx] = 7; }
				else{ Register[ID - 1].HPVstageE[xx] = 2; }
			}
			if ( (ih==5  && Register[ID-1].ARTstage==0)||
				(ih==5  && Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks>=104)){ //&& Register[ID-1].ARTweeks<104
				if (p < From2to1depLatent){ Register[ID - 1].HPVstageE[xx] = 1; }
				else if (p < From2to1depLatent + From2to3depLatent){ Register[ID - 1].HPVstageE[xx] = 3; }
				else if (p < From2to1depLatent + From2to3depLatent + From2to6depLatent){ Register[ID - 1].HPVstageE[xx] = 6; }
				else if (p < From2to1depLatent + From2to3depLatent + From2to6depLatent + From2to7depLatent){ Register[ID - 1].HPVstageE[xx] = 7; }
				else{ Register[ID - 1].HPVstageE[xx] = 2; }
			}
		}

	if (is == 3 && Register[ID - 1].AgeExact<30.0){
		if(ih == 0 ) { //|| (ih == 5 && Register[ID - 1].ARTweeks>=104)
			if (p < From3to2dep){ Register[ID - 1].HPVstageE[xx] = 2; }
			else if (p < From3to4dep + From3to2dep){ Register[ID - 1].HPVstageE[xx] = 4; }
			else{ Register[ID - 1].HPVstageE[xx] = 3; }
		}
		if (ih == 2 || ih == 1 || ih==3 || ih==4 ||  (ih==5  && Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks<104)||ih==6){ //&& Register[ID-1].ARTweeks<104
			if (p < From3to2depLate){ Register[ID - 1].HPVstageE[xx] = 2; }
			else if (p < From3to4depLate + From3to2depLate){ Register[ID - 1].HPVstageE[xx] = 4; }
			else{ Register[ID - 1].HPVstageE[xx] = 3; }
		}
		if ((ih==5  && Register[ID-1].ARTstage==0)||
			(ih==5  && Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks>=104)){ //&& Register[ID-1].ARTweeks<104
			if (p < From3to2depLatent){ Register[ID - 1].HPVstageE[xx] = 2; }
			else if (p < From3to4depLatent + From3to2depLatent){ Register[ID - 1].HPVstageE[xx] = 4; }
			else{ Register[ID - 1].HPVstageE[xx] = 3; }
		}
	}
	if (is == 3 && Register[ID - 1].AgeExact>=30.0 && Register[ID - 1].AgeExact<50.0){
            if(ih == 0 ) { //|| (ih == 5 && Register[ID - 1].ARTweeks>=104)
                if (p < From3to2dep30){ Register[ID - 1].HPVstageE[xx] = 2; }
                else if (p < From3to4dep30 + From3to2dep30){ Register[ID - 1].HPVstageE[xx] = 4; }
                else{ Register[ID - 1].HPVstageE[xx] = 3; }
            }
            if (ih == 2 ||ih == 1 || ih==3 || ih==4 ||  (ih==5  && Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks<104)||ih==6){ //&& Register[ID-1].ARTweeks<104
                if (p < From3to2depLate30){ Register[ID - 1].HPVstageE[xx] = 2; }
                else if (p < From3to4depLate30 + From3to2depLate30){ Register[ID - 1].HPVstageE[xx] = 4; }
                else{ Register[ID - 1].HPVstageE[xx] = 3; }
            }
            if ((ih==5  && Register[ID-1].ARTstage==0)||
                (ih==5  && Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks>=104)){ //&& Register[ID-1].ARTweeks<104
                if (p < From3to2depLatent30){ Register[ID - 1].HPVstageE[xx] = 2; }
                else if (p < From3to4depLatent30 + From3to2depLatent30){ Register[ID - 1].HPVstageE[xx] = 4; }
                else{ Register[ID - 1].HPVstageE[xx] = 3; }
            }
        }
        if (is == 3 && Register[ID - 1].AgeExact>=50.0 ){
            if(ih == 0 ) { 
                if (p < From3to2dep50){ Register[ID - 1].HPVstageE[xx] = 2; }
                else if (p < From3to4dep50 + From3to2dep50){ Register[ID - 1].HPVstageE[xx] = 4; }
                else{ Register[ID - 1].HPVstageE[xx] = 3; }
            }
            if (ih == 2 ||ih == 1 || ih==3 || ih==4 ||  (ih==5  && Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks<104)||ih==6){ //&& Register[ID-1].ARTweeks<104
                if (p < From3to2depLate50){ Register[ID - 1].HPVstageE[xx] = 2; }
                else if (p < From3to4depLate50 + From3to2depLate50){ Register[ID - 1].HPVstageE[xx] = 4; }
                else{ Register[ID - 1].HPVstageE[xx] = 3; }
            }
            if ((ih==5  && Register[ID-1].ARTstage==0)||
                (ih==5  && Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks>=104)){ //&& Register[ID-1].ARTweeks<104
                if (p < From3to2depLatent50){ Register[ID - 1].HPVstageE[xx] = 2; }
                else if (p < From3to4depLatent50 + From3to2depLatent50){ Register[ID - 1].HPVstageE[xx] = 4; }
                else{ Register[ID - 1].HPVstageE[xx] = 3; }
            }
        }   
    

	
	if (is == 5 && Register[ID - 1].DiagnosedCC==0){
			if(ih == 0 ) {
				if (p < From5to8dep){ Register[ID - 1].HPVstageE[xx] = 8;  }
				else if (p < From5to8dep + From5to11dep){ Register[ID - 1].HPVstageE[xx] = 11; }
				else{ Register[ID - 1].HPVstageE[xx] = 5; }
			}
			if (ih == 2 ||ih == 1 || ih==3 || ih==4 ||  (ih==5  && Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks<104)||ih==6){
				if (p < From5to8depHIV){ Register[ID - 1].HPVstageE[xx] = 8;  }
				else if (p < From5to8depHIV + From5to11depHIV){ Register[ID - 1].HPVstageE[xx] = 11; }
				else{ Register[ID - 1].HPVstageE[xx] = 5; }
			}
			if ((ih==5  && Register[ID-1].ARTstage==0)||
				(ih==5  && Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks>=104)){
				if (p < From5to8depART){ Register[ID - 1].HPVstageE[xx] = 8;  }
				else if (p < From5to8depART + From5to11depART){ Register[ID - 1].HPVstageE[xx] = 11; }
				else{ Register[ID - 1].HPVstageE[xx] = 5; }
			}
	}
	if (is == 8 && Register[ID - 1].DiagnosedCC==0){
			if(ih == 0 ) {
				if (p < From8to9dep){ Register[ID - 1].HPVstageE[xx] = 9;  }
				else if (p < From8to9dep + From8to12dep){ Register[ID - 1].HPVstageE[xx] = 12; }
				else{ Register[ID - 1].HPVstageE[xx] = 8; }
			}
			if (ih == 2 ||ih == 1 || ih==3 || ih==4 ||  (ih==5  && Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks<104)||ih==6){
				if (p < From8to9depHIV){ Register[ID - 1].HPVstageE[xx] = 9;  }
				else if (p < From8to9depHIV + From8to12depHIV){ Register[ID - 1].HPVstageE[xx] = 12; }
				else{ Register[ID - 1].HPVstageE[xx] = 8; }
			}
			if ((ih==5  && Register[ID-1].ARTstage==0)||
				(ih==5  && Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks>=104)){
				if (p < From8to9depART){ Register[ID - 1].HPVstageE[xx] = 9;  }
				else if (p < From8to9depART + From8to12depART){ Register[ID - 1].HPVstageE[xx] = 12; }
				else{ Register[ID - 1].HPVstageE[xx] = 8; }
			}	
		}
	if (is == 9 && Register[ID - 1].DiagnosedCC==0){
			if(ih == 0 ) {
				if (p < From9to10dep){ Register[ID - 1].HPVstageE[xx] = 10;  }
				else if (p < From9to10dep + From9to13dep){ Register[ID - 1].HPVstageE[xx] = 13; }
				else{ Register[ID - 1].HPVstageE[xx] = 9; }
			}
			if (ih == 2 ||ih == 1 || ih==3 || ih==4 ||  (ih==5  && Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks<104)||ih==6){
				if (p < From9to10depHIV){ Register[ID - 1].HPVstageE[xx] = 10;  }
				else if (p < From9to10depHIV + From9to13depHIV){ Register[ID - 1].HPVstageE[xx] = 13; }
				else{ Register[ID - 1].HPVstageE[xx] = 9; }
			}
			if ((ih==5  && Register[ID-1].ARTstage==0)||
				(ih==5  && Register[ID-1].ARTstage==1 && Register[ID-1].ARTweeks>=104)){
				if (p < From9to10depART){ Register[ID - 1].HPVstageE[xx] = 10;  }
				else if (p < From9to10depART + From9to13depART){ Register[ID - 1].HPVstageE[xx] = 13; }
				else{ Register[ID - 1].HPVstageE[xx] = 9; }
			}
				
		}
	if (is == 10 && Register[ID - 1].DiagnosedCC==0){
			if (p < From10to14) {Register[ID - 1].HPVstageE[xx] = 14;}
			else if (p < From10to14 + From10toDead){ 
				RSApop.SetToDead(ID);
				Register[ID - 1].HPVstageE[xx] = 16;
				DiagCCPost2000.out[SimCount2][0] += 1;
				RSApop.NewCancerDeath[Register[ID - 1].AgeGroup][CurrYear-StartYear] += 1;
				if(Register[ID - 1].HIVstage==0){RSApop.NewCancerDeathHIV[Register[ID - 1].AgeGroup][CurrYear-StartYear] += 1;}
				if(Register[ID - 1].HIVstage>0 && Register[ID - 1].HIVstage!=5){RSApop.NewCancerDeathHIV[2*18 + Register[ID - 1].AgeGroup][CurrYear-StartYear] += 1;}
				if(Register[ID - 1].HIVstage==5){RSApop.NewCancerDeathHIV[18 + Register[ID - 1].AgeGroup][CurrYear-StartYear] += 1;}
			}
			else{ Register[ID - 1].HPVstageE[xx] = 10; }
		}
		
	if (is == 6){
			if (ih == 0 || (ih == 5 && Register[ID - 1].ARTweeks >= 104)){
				if (p < From6to1){ 
					Register[ID - 1].HPVstageE[xx] = 1;
				}
				else{ Register[ID - 1].HPVstageE[xx] = 6; }
			}
			if (ih == 1 || ih == 3 || ih == 4 || (ih == 5 && Register[ID - 1].ARTweeks<104)||ih==6){
				if (p < From6to1Late){ 
					Register[ID - 1].HPVstageE[xx] = 1; 
			}
				else{ Register[ID - 1].HPVstageE[xx] = 6; }
			}
			if (ih == 2){
				if (p < From6to1Latent){ 
					Register[ID - 1].HPVstageE[xx] = 1; 
				}
				else{ Register[ID - 1].HPVstageE[xx] = 6; }
			}
		}
	if (is == 7){
			if (p < From7to0){ Register[ID - 1].HPVstageE[xx] = 0; }
			else{ Register[ID - 1].HPVstageE[xx] = 7; }
		}	
}
void HPVTransition::GetNewStageM(int ID, double p, int type)
{
	int is, xx, ih;
	xx = type;

		is = Register[ID - 1].HPVstage[xx];
		ih = Register[ID - 1].HIVstage;
		if (is == 1){
			if (ih == 0 || (ih == 5 && Register[ID - 1].ARTweeks>=104)){
				if (p < From1to7dep){ Register[ID - 1].HPVstageE[xx] = 7; }
				else if (p < From1to7dep + From1to6dep){ Register[ID - 1].HPVstageE[xx] = 6; }
				else{ Register[ID - 1].HPVstageE[xx] = 1; }
			}
			if (ih == 1 || ih==3 || ih==4 || (ih==5 && Register[ID-1].ARTweeks<104)||ih==6){
				if (p < From1to7depAcute){ Register[ID - 1].HPVstageE[xx] = 7; }
				else if (p < From1to7depAcute + From1to6depAcute){ Register[ID - 1].HPVstageE[xx] = 6; }
				else{ Register[ID - 1].HPVstageE[xx] = 1; }
			}
			if (ih == 2){
				if (p < From1to7depLatent){ Register[ID - 1].HPVstageE[xx] = 7; }
				else if (p < From1to7depLatent + From1to6depLatent){ Register[ID - 1].HPVstageE[xx] = 6; }
				else{ Register[ID - 1].HPVstageE[xx] = 1; }
			}
		}
		if (is == 6){
			if (ih == 0 || (ih == 5 && Register[ID - 1].ARTweeks >= 104)){
				if (p < From6to1){ 
					Register[ID - 1].HPVstageE[xx] = 1;
					//if (FixedUncertainty == 1){ NewHPVM.out[CurrSim - 1][CurrYear - StartYear] += 1; } 
				}
				else{ Register[ID - 1].HPVstageE[xx] = 6; }
			}
			if (ih == 1 || ih == 3 || ih == 4 || (ih == 5 && Register[ID - 1].ARTweeks<104)||ih==6){
				if (p < From6to1Late){ 
					Register[ID - 1].HPVstageE[xx] = 1; 
					//if (FixedUncertainty == 1){ NewHPVM.out[CurrSim - 1][CurrYear - StartYear] += 1; }
				}
				else{ Register[ID - 1].HPVstageE[xx] = 6; }
			}
			if (ih == 2){
				if (p < From6to1Latent){ 
					Register[ID - 1].HPVstageE[xx] = 1; 
					//if (FixedUncertainty == 1){ NewHPVM.out[CurrSim - 1][CurrYear - StartYear] += 1; }
				}
				else{ Register[ID - 1].HPVstageE[xx] = 6; }
			}
		}
		
		if (is == 7){
			if (p < From7to0){ Register[ID - 1].HPVstageE[xx] = 0; }
			else{ Register[ID - 1].HPVstageE[xx] = 7; }
		}	
}

void Indiv::GetNewHPVstate(int ID, double p, int type) 
{
	double Prob1, xx, yy; // , rr;
	int iy, SimCount2, zz; //, ic, id;  , is;

	if (HIVstage==0) {zz=0;}
	else if(HIVstage==5) {zz=1;} //||HIVstage==6
	else  {zz=2;}	

	SimCount2 = (CurrSim - 1) / IterationsPerPC;
	double rr[9];

	for (int ii = 0; ii < 9; ii++){
		rr[ii] = rg.Random();		
	}	
	
	if (HPVstage[type] == 0 ){
		// Note that we refer to opposite sex when getting transm prob
		if (SexInd == 0){
			Prob1 = HPVTransF[type].GetTransmProb(ID, type); 
		}
		else{
			Prob1 = HPVTransM[type].GetTransmProb(ID, type);
		}
		if (p<Prob1){
			HPVstageE[type] = 1;
			//if (FixedUncertainty == 1 && SexInd==0){ NewHPVM.out[CurrSim - 1][CurrYear - StartYear] += 1; }
			//if ( SexInd==1){ NewHPVF.out[CurrSim - 1][CurrYear - StartYear] += 1; }
		}
		else{ HPVstageE[type] = 0; }
	}
	//if(CIN3distribution == 1){
		if (HPVstage[type] == 4 && HPVstage[0] != 5 && HPVstage[1] != 5 && HPVstage[2] != 5 && HPVstage[3] != 5 && HPVstage[4] != 5 &&
			HPVstage[5] != 5 && HPVstage[6] != 5 && HPVstage[7] != 5 && HPVstage[8] != 5 && HPVstage[9] != 5 && HPVstage[10] != 5 && 
			HPVstage[11] != 5 && HPVstage[12] != 5 && HPVstage[0] < 8 && HPVstage[1] < 8 && HPVstage[2] < 8 && HPVstage[3] < 8 && HPVstage[4] < 8 &&
			HPVstage[5] < 8 && HPVstage[6] < 8 && HPVstage[7] < 8 && HPVstage[8] < 8 && HPVstage[9] < 8 && HPVstage[10] < 8 && 
			HPVstage[11] < 8 && HPVstage[12] < 8){
			if (TimeinCIN3[type] > WeibullCIN3[type]) { 
				HPVstageE[type] = 5; 
				
			}
			else { HPVstageE[type] = 4; }
		}	
	//}
	else {
		if (SexInd == 0){
				HPVTransM[type].GetNewStageM(ID, p, type);
		}
		else{
				HPVTransF[type].GetNewStageF(ID, p, type); 	
		}
	}
	if(SexInd==1){
		if ((HPVstage[type] == 0 && HPVstageE[type] == 1)  ){  //||(HPVstage[type] == 6 && HPVstageE[type] == 1)
			NewHPV[type].out[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
		}
		if (HPVstage[type] == 4 && HPVstageE[type] == 5){
			DiagCCPost2000.out[SimCount2][1] += 1;
			RSApop.NewCancer[AgeGroup][CurrYear-StartYear] += 1;
			NewCC[type].out[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
		}	
		if ((HPVstage[type] == 5 && HPVstageE[type] == 11)||
			(HPVstage[type] == 8 && HPVstageE[type] == 12)||
			(HPVstage[type] == 9 && HPVstageE[type] == 13)||
			(HPVstage[type] == 10 && HPVstageE[type] == 14)){
				DiagnosedCC=1;
				RSApop.NewDiagCancer[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
				if(type==0||type==1) {RSApop.NewDiagCancer1618[AgeGroup][CurrYear-StartYear] += 1;}
				if(HIVstage==5 ) {RSApop.NewDiagCancerART[CurrYear-StartYear] += 1;}
		}
		//if(CIN3distribution == 1){
			if (HPVstage[type] == 3 && HPVstageE[type] == 4){
				if(HIVstage==0 && Age50==0){
					WeibullCIN3[type] = static_cast<int> (HPVTransF[type].AveDuration[3]*pow(-log(rr[0]), 1/HPVTransF[type].CIN3shape));
				}
				if(HIVstage==0 && Age50==1){
					WeibullCIN3[type] = static_cast<int> (HPVTransF[type].AveDuration[3]* HPVTransF[type].CIN3_50* pow(-log(rr[0]), 1/HPVTransF[type].CIN3shape));
				}
				else if((HIVstage==5 && ARTstage==0)||(HIVstage==5  && ARTstage==1 && ARTweeks>=104)) {
					WeibullCIN3[type] = static_cast<int> (HPVTransF[type].AveDuration[3]* (1.0/HPVTransF[type].CIN2_ART)*pow(-log(rr[0]), 1/HPVTransF[type].CIN3shape));
				}	
				else {
					WeibullCIN3[type] = static_cast<int> (HPVTransF[type].AveDuration[3]* (1.0/HPVTransF[type].CIN2_HIV)* pow(-log(rr[0]), 1/HPVTransF[type].CIN3shape));
				}
			}
			
			//CIN2+incidenceStageI
			if (HPVstage[type] == 2 && HPVstageE[type] == 3){
				if(HIVstage==0){
					NewCIN2neg.out[SimCount2][CurrYear-StartYear] += 1;
				}
				if(HIVstage==5 ) {
					NewCIN2art.out[SimCount2][CurrYear-StartYear] += 1;
				}	
				if(HIVstage>0 ){
					NewCIN2pos.out[SimCount2][CurrYear-StartYear] += 1;
				}
			}
		//}
		if (HPVstage[type] == 5 && HPVstageE[type] == 5 && CurrYear>=ImplementYR){
			RSApop.WeeksInStageI[zz*18 + AgeGroup][CurrYear-ImplementYR] += 1;
		}
		if (HPVstage[type] == 8 && HPVstageE[type] == 8 && CurrYear>=ImplementYR){
			RSApop.WeeksInStageII[zz*18 + AgeGroup][CurrYear-ImplementYR] += 1;
		}
		if (HPVstage[type] == 9 && HPVstageE[type] == 9 && CurrYear>=ImplementYR){
			RSApop.WeeksInStageIII[zz*18 + AgeGroup][CurrYear-ImplementYR] += 1;
		}
		if (HPVstage[type] == 10 && HPVstageE[type] == 10 && CurrYear>=ImplementYR){
			RSApop.WeeksInStageIV[zz*18 + AgeGroup][CurrYear-ImplementYR] += 1;
		}

		if (HPVstage[type] == 5 && HPVstageE[type] == 11){
			RSApop.StageDiag[0][CurrYear-StartYear] += 1;
			RSApop.StageIdiag[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
			DiagnosedCC=1;
			//DiagCCPost2000.out[CurrSim-1][0] += 1;
			if(rr[1]<0.192){StageIdeath = static_cast<int> (48.0 * 3.08 * pow(-log(rr[2]), 1.0/1.23));}
			else{StageIrecover = 8 ;}
			//StageIdeath = static_cast<int> (48.0 * 126.5 * pow(-log(rr[1]), 1.0/0.61));
		}
		if (HPVstage[type] == 8 && HPVstageE[type] == 12){
			RSApop.StageDiag[1][CurrYear-StartYear] += 1;
			RSApop.StageIIdiag[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
			
			DiagnosedCC=1;
			//DiagCCPost2000.out[CurrSim-1][0] += 1;
			if(rr[3]<0.466){StageIIdeath = static_cast<int> (48.0 * 2.39 * pow(-log(rr[4]), 1.0/1.17));}
			else{StageIIrecover = 24 ;}
			//StageIIdeath = static_cast<int> (48.0 * 16.28 * pow(-log(rr[2]), 1.0/0.67));
		}
		if (HPVstage[type] == 9 && HPVstageE[type] == 13){
			RSApop.StageDiag[2][CurrYear-StartYear] += 1;
			RSApop.StageIIIdiag[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
			
			DiagnosedCC=1;
			//DiagCCPost2000.out[CurrSim-1][0] += 1;
			if(rr[5]<0.715){StageIIIdeath = static_cast<int> (48.0 * 1.18 * pow(-log(rr[6]), 1.0/0.91));}
			else{StageIIIrecover = 24 ;}
			//StageIIIdeath = static_cast<int> (48.0 * 3.91 * pow(-log(rr[3]), 1.0/0.56));
		}
		if (HPVstage[type] == 10 && HPVstageE[type] == 14){
			RSApop.StageDiag[3][CurrYear-StartYear] += 1;
			RSApop.StageIVdiag[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
			
			DiagnosedCC=1;
			//DiagCCPost2000.out[CurrSim-1][0] += 1;
			if(rr[7]<0.931){StageIVdeath = static_cast<int> (48.0 * 0.46 * pow(-log(rr[8]), 1.0/0.9));}
			else{StageIVrecover = 24 ;}
			//StageIVdeath = static_cast<int> (48.0 * 0.53 * pow(-log(rr[4]), 1.0/0.78));
		}
		if (HPVstage[type] == 11 && HPVstageE[type] == 11){
			TimeinStageI += 1;
			if(CurrYear>=ImplementYR) {
				RSApop.WeeksInStageI[zz*18 + AgeGroup][CurrYear-ImplementYR] += 1;
			}
			
			if(TimeinStageI > StageIdeath && TimeinStageI < 480 && StageIdeath>0) { //has to die from cancer within 10 years 10*48=480 
				RSApop.SetToDead(ID); 
				RSApop.NewDiagCancerDeath[AgeGroup][CurrYear-StartYear] += 1;
				RSApop.NewCancerDeath[AgeGroup][CurrYear-StartYear] += 1;
				RSApop.NewCancerDeathHIV[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
				HPVstageE[type] = 16;
			}
			if(TimeinStageI > StageIrecover && StageIrecover>0){
				HPVstageE[type] = 15;	
			}
		}
		if (HPVstage[type] == 12 && HPVstageE[type] == 12){
			TimeinStageII += 1;
			if(CurrYear>=ImplementYR) {
				RSApop.WeeksInStageII[zz*18 + AgeGroup][CurrYear-ImplementYR] += 1;
			}
			if(TimeinStageII > StageIIdeath && TimeinStageII < 480 && StageIIdeath>0) { //has to die from cancer within 10 years 10*48=480 
				RSApop.SetToDead(ID);
				RSApop.NewDiagCancerDeath[AgeGroup][CurrYear-StartYear] += 1;
				RSApop.NewCancerDeath[AgeGroup][CurrYear-StartYear] += 1;
				RSApop.NewCancerDeathHIV[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
				HPVstageE[type] = 16;
			}
			if(TimeinStageII > StageIIrecover && StageIIrecover>0){
				HPVstageE[type] = 15;	
			}
		}
		if (HPVstage[type] == 13 && HPVstageE[type] == 13){
			TimeinStageIII += 1;
			if(CurrYear>=ImplementYR) {
				RSApop.WeeksInStageIII[zz*18 + AgeGroup][CurrYear-ImplementYR] += 1;
			}
			if(TimeinStageIII > StageIIIdeath && TimeinStageIII < 480 && StageIIIdeath>0) { //has to die from cancer within 10 years 10*48=480 
				RSApop.SetToDead(ID); 
				RSApop.NewDiagCancerDeath[AgeGroup][CurrYear-StartYear] += 1;
				RSApop.NewCancerDeath[AgeGroup][CurrYear-StartYear] += 1;
				RSApop.NewCancerDeathHIV[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
				HPVstageE[type] = 16;
			}
			if(TimeinStageIII > StageIIIrecover && StageIIIrecover>0){
				HPVstageE[type] = 15;	
			}
		}
		if (HPVstage[type] == 14 && HPVstageE[type] == 14){
			TimeinStageIV += 1;
			if(CurrYear>=ImplementYR) {
				RSApop.WeeksInStageIII[zz*18 + AgeGroup][CurrYear-ImplementYR] += 1;
			}		
			if(TimeinStageIV > StageIVdeath && TimeinStageIV < 480 && StageIVdeath>0) { //has to die from cancer within 10 years 10*48=480 
				RSApop.SetToDead(ID); 
				RSApop.NewDiagCancerDeath[AgeGroup][CurrYear-StartYear] += 1;
				RSApop.NewCancerDeath[AgeGroup][CurrYear-StartYear] += 1;
				RSApop.NewCancerDeathHIV[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
				HPVstageE[type] = 16;
			}
			if(TimeinStageIV > StageIVrecover && StageIVrecover>0){
				HPVstageE[type] = 15;	
			}
		}
			
		if (HPVstage[type] == 4 && HPVstageE[type] == 4 && HPVstage[0] != 5 && HPVstage[1] != 5 && HPVstage[2] != 5 && HPVstage[3] != 5 && HPVstage[4] != 5 &&
			HPVstage[5] != 5 && HPVstage[6] != 5 && HPVstage[7] != 5 && HPVstage[8] != 5 && HPVstage[9] != 5 && HPVstage[10] != 5 && 
			HPVstage[11] != 5 && HPVstage[12] != 5 && HPVstage[0] < 8 && HPVstage[1] < 8 && HPVstage[2] < 8 && HPVstage[3] < 8 && HPVstage[4] < 8 &&
			HPVstage[5] < 8 && HPVstage[6] < 8 && HPVstage[7] < 8 && HPVstage[8] < 8 && HPVstage[9] < 8 && HPVstage[10] < 8 && 
			HPVstage[11] < 8 && HPVstage[12] < 8 ){
			TimeinCIN3[type] += 1;
		}
	}
}

void Pop::GetNumbersByHPVstageAge()
{
	int ic, iy, ig, is, xx, ir;
	int tpp = Register.size();
	for (xx = 0; xx < 13; xx++){
		for (ic = 0; ic < tpp; ic++){
			if (Register[ic].AliveInd == 1   && Register[ic].AgeGroup>1){
				Register[ic].GetRiskGroup();
				iy = Register[ic].AgeGroup-2;
				ig = Register[ic].SexInd;
				is = Register[ic].HPVstage[xx];
				if(is>7){is=5;}
				ir = Register[ic].risk_for_start;
				AdultHPVstageTrendAge[iy + 16 * ir][is + xx * 8 + 8*13 * ig] += 1;
			}
		}
	}
}
void Pop::SaveAdultHPVstageAge(const char *filout)
{
	int iy, is;
	ostringstream s;

	if (process_num >0){
		s << process_num << "_" << filout;
	}
	else{
		s << filout;
	}
	
	ofstream file(s.str().c_str());
	
	for (iy = 0; iy < 304; iy++){
		for (is = 0; is<208; is++){
			file << right << AdultHPVstageTrendAge[iy][is] << "	";
			//file << right << AdultHPVstageTrend[iy][is] << "	";
		}
		file << endl;
	}
	file.close();
}

void Pop::GetNumbersByHPVstage()
{
	int ic, iy, ig, is,  xx;
	
	iy = CurrYear - StartYear;
	int tpp = Register.size();
	/*for (ic = 0; ic < tpp; ic++){
		for (xx = 0; xx < 13; xx++){
			if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup >= 3){
				// Note that this gives total pop at ages 15+
				ig = Register[ic].SexInd;
				is = Register[ic].HPVstage[xx];
				AdultHPVstageTrend[iy][ig * 8 + is + 16*xx] += 1;
			}
		}
	}*/
	for (ic = 0; ic < tpp; ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup >= 3 && Register[ic].SexInd==1){
				// Note that this gives total pop at ages 15+
				is = Register[ic].TrueStage;
				AdultHPVstageTrend[is][iy] += 1;
			}
	}
	//cout << CurrYear << " "<< AdultHPVstageTrend[0][iy] << " " <<AdultHPVstageTrend[1][iy] << " " <<AdultHPVstageTrend[2][iy] << " " <<AdultHPVstageTrend[3][iy] << endl;
}
void Pop::SaveAdultHPVstage(char* filout)
{
	int iy, is;
	ofstream file(filout);

	for (iy = 0; iy < 136; iy++){
		for (is = 0; is<208; is++){
			file << right << AdultHPVstageTrend[iy][is] << "	";
		}
		file << endl;
	}
	file.close();
}

PostOutputArray3::PostOutputArray3(){}

void PostOutputArray3::RecordSample(const char *filout, int type)
{
	int i, c;
	ostringstream s;

	if (process_num >0){
		s << process_num << "_" << type << filout;
	}
	else{
		s << type << filout;
	}
	
	string path = "./output/" + s.str();
	ofstream file(path.c_str()); // Converts s to a C string

	for (i = 0; i<54; i++){
		file << setw(6) << right << i << "	";
		for (c = 0; c<136; c++){
			file << "	" << setw(10) << right << out[i][c];
		}
		file << endl;
	}
	file.close();
	
}

PostOutputArray4::PostOutputArray4(){}

void PostOutputArray4::RecordSample(const char *filout, int type)
{
	int i, c;
	ostringstream s;

	if (process_num >0){
		s << process_num << "_" << type << filout;
	}
	else{
		s << type << filout;
	}
	
	string path = "./output/" + s.str();
	ofstream file(path.c_str()); // Converts s to a C string

	for (i = 0; i<ParamCombs; i++){
		file << setw(6) << right << i << "	";
		for (c = 0; c<136; c++){
			file << "	" << setw(10) << right << out[i][c];
		}
		file << endl;
	}
	file.close();

}

void Pop::SaveLifetimeCIN3(const char* filout)
{
	int ic, xx;
	double AgeNow;
	int tpp = Register.size();
	
	ostringstream s;

	if (process_num >0){
		s << process_num << "_" << filout;
	}
	else{
		s << filout;
	}
	
	ofstream file(s.str().c_str(), std::ios::app);
	//ofstream file;
	//file.open(filout, std::ios::app);

	for (ic = 0; ic<tpp; ic++){
		for (xx = 0; xx < 13; xx++){
			if (Register[ic].SexInd == 1 && (Register[ic].HPVstage[xx] == 4 || Register[ic].HPVstage[xx] == 5) &&
				Register[ic].AgeGroup >= 3 && Register[ic].AliveInd==1){
				AgeNow = 0.5 + CurrYear - Register[ic].DOB;
				file << CurrSim << "	" << ic << "	" << xx << "	" << AgeNow << "	" << Register[ic].HPVstage[xx] << "	" << Register[ic].TimeinCIN3[xx] << "	" << Register[ic].WeibullCIN3[xx] << endl;
			}
		}
	}
	file.close();
}

void Pop::GetInitHPVstage(int type)
{
	int ic, iy, yy;
	for (iy = 0; iy < 18; iy++){
		for (yy = 0; yy < 8; yy++){
			InitHPVStage[iy][yy] = 0;
		}
	}
	int tpp = Register.size();
	for (ic = 0; ic < tpp; ic++){
		iy = Register[ic].AgeGroup;
		yy = Register[ic].HPVstage[type] ;
		if (Register[ic].SexInd == 1 ) {
			InitHPVStage[iy][yy] += 1;
		}
	}
}

void Pop::SaveInitHPVstage(char* filout, int type)
{
	int iy, is;
	ostringstream s;

	s << type << filout;
	
	//ofstream file(filout);
	ofstream file(s.str().c_str());
	
	for (iy = 0; iy < 18; iy++){
		for (is = 0; is<8; is++){
			file << right << InitHPVStage[iy][is] << "	";
		}
		file << endl;
	}
	file.close();
}

void Indiv::AssignTimeinCIN3(int age_group, double p, int type)
{
	int a = age_group;
	int  iy;
	double CumProb;
	if (a < 3){ TimeinCIN3[type] = 0; } //no one younger than 15 spends any time in CIN3 
	else{
		CumProb = 0.0;
		for (iy = 0; iy < 60; iy++){

			CumProb += PropinCIN3[a - 3][iy];
			//cout << a << " " << iy << " " <<PropinCIN3[a - 3][iy] << " " << CumProb << " " << p << endl;
			if (p < CumProb){
				TimeinCIN3[type] = static_cast<int> (iy * 52);
				break;
			}
		}
		//cout << age_group << " " << TimeinCIN3[type] << endl;
	}
}

void ReadTimeinCIN3()
{
	int ia, iy;
	ifstream file;
	stringstream s;
	s << "timeCIN3Prop.txt";
	string path = "./input/" + s.str();
	
	file.open(path.c_str());
	if (file.fail()) {
		cerr << "Could not open timeCIN3Prop.txt\n";
		exit(1);
	}
	for (ia = 0; ia<15; ia++){
		for (iy = 0; iy<60; iy++){
			file >> PropinCIN3[ia][iy];
			// std::cout << "age	" << ia << "	time	" << iy << "	proportion	" << PropinCIN3[ia][iy] << std::endl; 
		}
	}

	file.close();
}

void InitialiseHPV(){

	int xx, ii, jj, iy, is;
	HPVTransM[0].HPVObs(0, 0, 0, 0, 0, 4, 0, 0, 0,0);
	HPVTransF[0].HPVObs(1, 0, 1, 0, 0, 13, 0, 0, 2,0);
	HPVTransM[1].HPVObs(0, 0, 0, 0, 0, 4, 0, 0, 0,0);
	HPVTransF[1].HPVObs(1, 0, 1, 0, 0, 13, 0, 0, 2,0);
	HPVTransM[2].HPVObs(0, 0, 0, 0, 0, 4, 0, 0, 0,0);
	HPVTransF[2].HPVObs(1, 0, 1, 0, 0, 11, 0, 0, 2,0);
	HPVTransM[3].HPVObs(0, 0, 0, 0, 0, 4, 0, 0, 0,0);
	HPVTransF[3].HPVObs(1, 0, 1, 0, 0, 11, 0, 0, 2,0);
	HPVTransM[4].HPVObs(0, 0, 0, 0, 0, 4, 0, 0, 0,0);
	HPVTransF[4].HPVObs(1, 0, 1, 0, 0, 11, 0, 0, 2,0);
	HPVTransM[5].HPVObs(0, 0, 0, 0, 0, 4, 0, 0, 0,0);
	HPVTransF[5].HPVObs(1, 0, 1, 0, 0, 11, 0, 0, 2,0);
	HPVTransM[6].HPVObs(0, 0, 0, 0, 0, 4, 0, 0, 0,0);
	HPVTransF[6].HPVObs(1, 0, 1, 0, 0, 11, 0, 0, 2,0);
	HPVTransM[7].HPVObs(0, 0, 0, 0, 0, 4, 0, 0, 0,0);
	HPVTransF[7].HPVObs(1, 0, 1, 0, 0, 11, 0, 0, 2,0);
	HPVTransM[8].HPVObs(0, 0, 0, 0, 0, 4, 0, 0, 0,0);
	HPVTransF[8].HPVObs(1, 0, 1, 0, 0, 11, 0, 0, 2,0);
	HPVTransM[9].HPVObs(0, 0, 0, 0, 0, 4, 0, 0, 0,0);
	HPVTransF[9].HPVObs(1, 0, 1, 0, 0, 11, 0, 0, 2,0);
	HPVTransM[10].HPVObs(0, 0, 0, 0, 0, 4, 0, 0, 0,0);
	HPVTransF[10].HPVObs(1, 0, 1, 0, 0, 11, 0, 0, 2,0);
	HPVTransM[11].HPVObs(0, 0, 0, 0, 0, 4, 0, 0, 0,0);
	HPVTransF[11].HPVObs(1, 0, 1, 0, 0, 11, 0, 0, 2,0);
	HPVTransM[12].HPVObs(0, 0, 0, 0, 0, 3, 0, 0, 0,0);
	HPVTransF[12].HPVObs(1, 0, 1, 0, 0, 11, 0, 0, 2,0);

	for (xx = 0; xx < 13; xx++){
		for (ii = 0; ii < 54; ii++){
			for (jj = 0; jj < 136; jj++){
				NewCC[xx].out[ii][jj] = 0;
				NewHPV[xx].out[ii][jj] = 0;
			}
		}
	}
	for (ii = 0; ii < samplesize; ii++){
			for (jj = 0; jj < 136; jj++){
				NewCIN2neg.out[ii][jj] = 0; 
				NewCIN2pos.out[ii][jj] = 0; 
				NewCIN2art.out[ii][jj] = 0; 
				RandomUniformHPV.out[ii][jj] = 0; 
				HPVparamsLogL.out[ii][jj] = 0; 
			}
	}
	for (ii = 0; ii < 54; ii++){
		for (jj = 0; jj < 136; jj++){
			//NewHIV.out[ii][jj] = 0;
			RSApop.StageIdiag[ii][jj] = 0;
			RSApop.StageIIdiag[ii][jj] = 0;
			RSApop.StageIIIdiag[ii][jj] = 0;
			RSApop.StageIVdiag[ii][jj] = 0;
		}
	}
	for (iy = 0; iy < 18; iy++){
		for (is = 0; is<208; is++){
			RSApop.AdultHPVstageTrendAge[iy][is] = 0;
		}
	}

	for (iy = 0; iy < 136; iy++){
		for (is = 0; is<208; is++){
			RSApop.AdultHPVstageTrend[iy][is] = 0;
		}
	}
	for (iy = 0; iy < 500; iy++){
		for (is = 0; is<14; is++){
			RSApop.MacDprev[iy][is] =0;
		}
	}
	for (iy = 0; iy < 54; iy++){
		for (is = 0; is<136; is++){
			RSApop.NewScreen[iy][is] =0;
			RSApop.NewReflex[iy][is] =0;
			RSApop.NewTxV[iy][is] =0;	
			RSApop.NewHPVScreen[iy][is] =0;
			RSApop.NewColposcopy[iy][is] =0;
			RSApop.NewLLETZ[iy][is] =0;
			RSApop.NewUnnecessary[iy][is] =0;
			RSApop.NewVAT[iy][is] =0;
			RSApop.NewThermal[iy][is] =0;
			RSApop.GetReferred[iy][is] =0;
		}
	}
	for (iy = 0; iy < 36; iy++){
		for (is = 0; is<136; is++){
			RSApop.NewVACC[iy][is] =0;
			RSApop.ModelVaccCoverage[iy][is]=0;
		}
	}
	for (iy = 0; iy < 8; iy++){
		for (is = 0; is<136; is++){
			ScreenProb[iy][is] =0.0;
		}
	}
	for (iy = 0; iy < 136; iy++){
		for (is = 0; is<54; is++){
			RSApop.ModelCoverage[is][iy]=0;
			RSApop.ModelReflexCoverage[is][iy]=0;
			RSApop.ModelTxVCoverage[is][iy]=0;
			RSApop.ModelHPVCoverage[is][iy]=0;
			RSApop.ModelColpCoverage[is][iy]=0;
			RSApop.ModelLLETZCoverage[is][iy]=0;
			RSApop.ModelUnnecessaryCoverage[is][iy]=0;
			RSApop.ModelVATCoverage[is][iy]=0;
			RSApop.ModelThermalCoverage[is][iy]=0;
		}
	}

	for (iy = 0; iy<54; iy++){
		for (is = 0; is<136; is++){
			RSApop.PopPyramid[iy][is] = 0;
			RSApop.PopPyramidMale[iy][is] = 0;
			RSApop.HPVprevVT[iy][is] = 0;
			RSApop.HPVprevAll[iy][is] = 0;
			RSApop.CIN2prev[iy][is] = 0;
			RSApop.NewCancerDeathHIV[iy][is] = 0;
			RSApop.NewDiagCancer[iy][is] = 0;

		}

	}
	for (iy = 0; iy<6; iy++){
		for (is = 0; is<136; is++){
			RSApop.PopPyramid9[iy][is] = 0;
		}

	}
	for (iy = 0; iy<61; iy++){
		for (is = 0; is<136; is++){
			RSApop.PopPyramid61[iy][is] = 0;
			RSApop.PopVaxx61[iy][is] = 0;
		}
	}
	
	for (iy = 0; iy<18; iy++){
		for (is = 0; is<136; is++){
			RSApop.PopPyramidAll[iy][is] = 0;
			
			RSApop.NewDiagCancer1618[iy][is] = 0;
			RSApop.NewCancer[iy][is] = 0;
			RSApop.NewCancerDeath[iy][is] = 0;
		}
	}

	for (iy = 0; iy<54; iy++){
		for (is = 0; is<101; is++){
			RSApop.WeeksInStageI[iy][is] = 0;
			RSApop.WeeksInStageII[iy][is] = 0;
			RSApop.WeeksInStageIII[iy][is] = 0;
			RSApop.WeeksInStageIV[iy][is] = 0;
		}
	}

	for (is = 0; is<136; is++){
			RSApop.PopPyramidAllART[is] = 0;
			RSApop.NewDiagCancerART[is] = 0;
	}
	for (is = 0; is<136; is++){
		for(iy=0; iy<4; iy++){
			RSApop.StageDiag[iy][is] = 0;
		}	
	}
	for (iy = 0; iy < 18; iy++){
		for (is = 0; is<136; is++){
			RSApop.HIVDeath[iy][is] =0;
		}
	}
	for (iy = 0; iy < 18; iy++){
		for (is = 0; is<136; is++){
			RSApop.HIVDeathM[iy][is] =0;
		}
	}
	for(iy=0; iy<8; iy++){RSApop.TotScreens[iy]=0;}

	for(iy=0; iy<136; iy++){
		HPVtransitionCC.TotPop30to65neg[iy]=0; HPVtransitionCC.TotPop30to65pos[iy]=0; 
		HPVtransitionCC.TotPop30to65art[iy]=0; HPVtransitionCC.TotPop18to60noart[iy]=0 ;
		HPVtransitionCC.TotAB30to65neg[iy]=0; HPVtransitionCC.TotAB30to65pos[iy]=0; 
		 HPVtransitionCC.TotAB30to65art[iy]=0; HPVtransitionCC.TotAB18to60noart[iy]=0;
		HPVtransitionCC.TotHSIL30to65neg[iy]=0; HPVtransitionCC.TotHSIL30to65pos[iy]=0;
		HPVtransitionCC.TotHSIL30to65art[iy]=0; HPVtransitionCC.TotHSIL18to60noart[iy]=0;

		HPVtransitionCC.TotPop15to65negF[iy]=0; HPVtransitionCC.TotPop15to65negM[iy]=0 ;
		HPVtransitionCC.TotPop15to65posF[iy]=0; HPVtransitionCC.TotPop15to65posM[iy]=0 ;
		HPVtransitionCC.TotPop15to65artF[iy]=0; HPVtransitionCC.TotPop15to65artM[iy]=0 ;
		HPVtransitionCC.TotHPV15to65negF[iy]=0; HPVtransitionCC.TotHPV15to65negM[iy]=0 ;
		HPVtransitionCC.TotHPV15to65posF[iy]=0; HPVtransitionCC.TotHPV15to65posM[iy]=0 ;
		HPVtransitionCC.TotHPV15to65artF[iy]=0; HPVtransitionCC.TotHPV15to65artM[iy]=0 ;
		HPVtransitionCC.TotHSILNEG[iy]=0; HPVtransitionCC.TotHSILPOS[iy]=0; HPVtransitionCC.TotHSILART[iy]=0;
		HPVtransitionCC.TotCCNEG[iy]=0;	HPVtransitionCC.TotCCPOS[iy]=0;
		HPVtransitionCC.TotCCART[iy]=0; HPVtransitionCC.TotPop15upnegF[iy]=0; 
		HPVtransitionCC.TotPop15upposF[iy]=0; HPVtransitionCC.TotPop15upartF[iy]=0;
	}
}

void SimulateHPVparamsType()
{
	int SimCount2, seed, ii, ind, xx;
	double x, y, a, b, p, q;
	double r[34];

	ind = 2;
	SimCount2 = (CurrSim-1) / IterationsPerPC;
		
	if (FixedUncertainty == 0){
		seed = SimCount2 * 91 + process_num * 7927 ;
		
		CRandomMersenne rg(seed);
		for (ii = 0; ii < 34; ii++){
			r[ii] = rg.Random();
			if (CurrSim == (SimCount2 * IterationsPerPC + 1 )){
				RandomUniformHPV.out[SimCount2][ii] = r[ii];
			}	
		}
	}
	else{
		for (ii = 0; ii < 34; ii++){
			r[ii] = RandomUniformHPV.out[SimCount2][ii];
		}
	}
	
	//Read the following median values from the prevalence calibration:
	// The transmission prob M 
	// The transmission prob F 
	// The average duration of HPV DNA positive infection for males 
	// The average duration of HPV DNA positive infection for females
	// The average duration of latency for males
	// The average duration of latency for females
	// The average duration of immunity for males
	// The average duration of immunity for females
	// The proportion that will become latently infected after DNA clearance
	// The impact of HIV stage on reactivation of latency
	// The impact of HIV stage on HPV infection duration
	// The std deviation of the study effect
	for(xx = 0; xx < 13; xx++){
		HPVTransF[xx].prog_stageI = 0.225;  
		HPVTransF[xx].prog_stageII = 0.3;  
		HPVTransF[xx].prog_stageIII = 0.45;
	}			
	ReadHPVparams();
	
	if(UseMedians==0){
		
		if(CCcalib==0){
			// Simulate a HPV duration multiplier 
			HPVTransF[0].durMult =  0.2 + 0.8 * r[0]; 
				HPVTransF[0].AveDuration[0]  = HPVTransF[0].AveDuration[0] * HPVTransF[0].durMult ;	
			
			HPVTransF[1].durMult =  0.2 + 0.8 * r[1]; 
				HPVTransF[1].AveDuration[0]  = HPVTransF[1].AveDuration[0] * HPVTransF[1].durMult ;	
			
			for (xx = 2; xx < 13; xx++){
				HPVTransF[xx].durMult =  0.2 + 0.8 * r[2]; 
				HPVTransF[xx].AveDuration[0]  = HPVTransF[xx].AveDuration[0] * HPVTransF[xx].durMult ;	
			}
			
			// Simulate the proportion that will progress from normal to CIN1 (looks like 25% from landon's data, 
			// 23% from Insinga 2011, 32% from Insinga 2007)
			a = 4.742;
			b = 13.498;
			//Type 16
			p = r[3];
			q = 1.0 - r[3];
			cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
			HPVTransF[0].prop1pro = x; 
			
			//Type 18
			HPVTransF[1].prop1pro = x *(1.0 - 0.5 * r[4]); 	
			
			//Other
			for(xx = 2; xx < 13; xx++){
				HPVTransF[xx].prop1pro =  x *(1.0 - 0.5 * r[5]); 
			}
			
			for(xx = 0; xx < 13; xx++){
				HPVTransF[xx].prog_fromHPV = (1.0/HPVTransF[xx].AveDuration[0])*HPVTransF[xx].prop1pro; 
				HPVTransF[xx].reg_fromHPV = (1.0/HPVTransF[xx].AveDuration[0])*(1.0-HPVTransF[xx].prop1pro); 
			}

			// Simulate the progression rate from CIN1 to CIN2  
			a = 3.24;
			b = 36.0;
			//Type 16
			p = r[6];
			q = 1.0 - r[6];
			cdfgam(&ind, &p, &q, &x,  &a, &b, 0, 0);
			HPVTransF[0].prog_fromCIN1 = x ; 
			
			//Type 18
			HPVTransF[1].prog_fromCIN1 = x * r[7]  ; //anything between 0 and 1 times the rate for 16; //
			//Other
			for(xx = 2; xx < 13; xx++){
				HPVTransF[xx].prog_fromCIN1 =x * r[8]  ; 
			}

			// Simulate the rate of regressing from CIN1
			a = 4.6225;
			b = 10.75;
			p = r[12];
			q = 1.0 - r[12];
			cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
			HPVTransF[0].reg_fromCIN1 = x  ; 
			//Type 18
			HPVTransF[1].reg_fromCIN1 =  x * (1.0 + r[13]) ; 
			//Other
			for(xx = 2; xx < 13; xx++){
				HPVTransF[xx].reg_fromCIN1 =  x * (1.0 + r[14]) ;  
			}

			// Simulate the rate of regressing from CIN2
			a = 131.1025;
			b = 286.25;
			p = r[15];
			q = 1.0 - r[15];
			cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
			//x=0.0;
			HPVTransF[0].reg_fromCIN2 = x  ; 
			//Type 18
			HPVTransF[1].reg_fromCIN2 =  x * (1.0 + r[16]) ; 
			//Other
			for(xx = 2; xx < 13; xx++){
				HPVTransF[xx].reg_fromCIN2 =  x * (1.0 + r[17]) ; 
			}

			for(xx = 0; xx < 13; xx++){
				HPVTransF[xx].CIN2reg30 = 0.45 + 0.3*r[31]; //1; //
				HPVTransF[xx].CIN1_HIV = 2 + 3.32*r[19];  
				HPVTransF[xx].prog_ART =  0.55 + 0.35*r[21];  
				HPVTransF[xx].CIN1_ART = max(HPVTransF[xx].CIN1_HIV * HPVTransF[xx].prog_ART, 1.0);
				HPVTransF[xx].reg_HIV = 0.56 + 0.26*r[22]; 
				HPVTransF[xx].reg_ART = min(HPVTransF[xx].reg_HIV * (1.3 + 0.7*r[23]), 1.0);
			}
			
			// Read the std deviation of the study effect (ABnormal prevalence)
			a = 16.0;
			b = 26.6667;
			p = r[26];
			q = 1.0 - r[26];
			cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
			HPVtransitionCC.FPClogL.VarStudyEffectAB = pow(x, 2.0);
			HPVtransitionCC.HouseholdLogL.VarStudyEffectAB = pow(x, 2.0);
			HPVtransitionCC.NOARTlogL.VarStudyEffectAB = pow(x, 2.0);
			HPVtransitionCC.ONARTlogL.VarStudyEffectAB = pow(x, 2.0);

		}
		
		if(CCcalib==1){
			// Simulate the progression rate from CIN2 to CIN3 
			a = 3.738;
			b = 64.444;
			//Type 16
			p = r[9];
			q = 1.0 - r[9];
			cdfgam(&ind, &p, &q, &x,  &a, &b, 0, 0);
			//x=1.0;
			HPVTransF[0].prog_fromCIN2 = x ; 
			//Type 18
			HPVTransF[1].prog_fromCIN2 = x * r[10] ; 
			//Other
			for(xx = 2; xx < 13; xx++){
				HPVTransF[xx].prog_fromCIN2 =x * r[11]  ;  
			}
			
			//The impact of HIV stage on CIN3 duration 
			a = 3.0426;
			b = 1.5674;
			p = r[18];
			q = 1.0 - r[18];
			cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
			
			int yy = r[24]*150;
			int zz = r[25]*10;

			for(xx = 0; xx < 13; xx++){
				HPVTransF[xx].CIN2prog30 = 2.0 + 1.0*r[32];
				HPVTransF[xx].CIN2prog50 = HPVTransF[xx].CIN2prog30*(1.0+r[33]);
				HPVTransF[xx].CIN3_50 = 1.0; //0.5 + 0.5*r[33];
				HPVTransF[xx].CIN2_HIV =  1.1 + 0.9*r[20]; 
				HPVTransF[xx].CIN2_ART = max(HPVTransF[xx].CIN2_HIV * HPVTransF[xx].prog_ART, 1.0);
				HPVTransF[xx].CIN3HIV = 1.0/HPVTransF[xx].CIN2_HIV; //0.66; //x; 
				HPVTransF[xx].CIN3HIV_ART = 1.0/HPVTransF[xx].CIN2_ART; //  //	
				HPVTransF[xx].AveDuration[3] = (5.0 + 15.0*(yy/150.0)) * 52.0; //10.0*52.0; //
				HPVTransF[xx].CIN3shape = 2.0 + 1.0*(zz/10.0); //2.5; //
				HPVTransF[xx].diag_stageI = r[27]*0.05; //0.15;  
				HPVTransF[xx].diag_stageII = 0.05 + 0.15*r[28]; //r[28]*0.225;  
				HPVTransF[xx].diag_stageIII = 0.4 + r[29]*0.4; //0.6;  
				HPVTransF[xx].diag_stageIV = 0.85 + r[30]*0.15;   //0.8;  
			}
		}
	}
	
	if (CurrSim == (SimCount2 * IterationsPerPC +1)){
		HPVparamsLogL.out[SimCount2][5] = HPVTransF[0].durMult;
		HPVparamsLogL.out[SimCount2][6] = HPVTransF[1].durMult;
		HPVparamsLogL.out[SimCount2][7] = HPVTransF[2].durMult;
		HPVparamsLogL.out[SimCount2][8] = HPVTransF[0].prop1pro;
		HPVparamsLogL.out[SimCount2][9] = HPVTransF[1].prop1pro;
		HPVparamsLogL.out[SimCount2][10] = HPVTransF[2].prop1pro;
		HPVparamsLogL.out[SimCount2][11] = HPVTransF[0].prog_fromCIN1;
		HPVparamsLogL.out[SimCount2][12] = HPVTransF[1].prog_fromCIN1;
		HPVparamsLogL.out[SimCount2][13] = HPVTransF[2].prog_fromCIN1;
		HPVparamsLogL.out[SimCount2][14] = HPVTransF[0].prog_fromCIN2;
		HPVparamsLogL.out[SimCount2][15] = HPVTransF[1].prog_fromCIN2;
		HPVparamsLogL.out[SimCount2][16] = HPVTransF[2].prog_fromCIN2;
		HPVparamsLogL.out[SimCount2][17] = HPVTransF[0].reg_fromCIN1;
		HPVparamsLogL.out[SimCount2][18] = HPVTransF[1].reg_fromCIN1;
		HPVparamsLogL.out[SimCount2][19] = HPVTransF[2].reg_fromCIN1;
		HPVparamsLogL.out[SimCount2][20] = HPVTransF[0].reg_fromCIN2;
		HPVparamsLogL.out[SimCount2][21] = HPVTransF[1].reg_fromCIN2;
		HPVparamsLogL.out[SimCount2][22] = HPVTransF[2].reg_fromCIN2;
		HPVparamsLogL.out[SimCount2][23] = HPVTransF[0].CIN3HIV;
		HPVparamsLogL.out[SimCount2][24] = HPVTransF[0].AveDuration[3];
		HPVparamsLogL.out[SimCount2][25] = HPVTransF[0].CIN3shape;
		HPVparamsLogL.out[SimCount2][26] = HPVTransF[0].CIN1_HIV;
		HPVparamsLogL.out[SimCount2][27] = HPVTransF[0].CIN2_HIV;
		HPVparamsLogL.out[SimCount2][28] = HPVTransF[0].prog_ART;
		HPVparamsLogL.out[SimCount2][29] = HPVTransF[0].reg_HIV; 
		HPVparamsLogL.out[SimCount2][30] = HPVTransF[0].reg_ART;
		HPVparamsLogL.out[SimCount2][31] = HPVtransitionCC.HouseholdLogL.VarStudyEffectAB;
		HPVparamsLogL.out[SimCount2][32] = HPVTransF[0].diag_stageI;  
		HPVparamsLogL.out[SimCount2][33] = HPVTransF[0].diag_stageII;  
		HPVparamsLogL.out[SimCount2][34] = HPVTransF[0].diag_stageIII;
		HPVparamsLogL.out[SimCount2][35] = HPVTransF[0].diag_stageIV;
		HPVparamsLogL.out[SimCount2][36] = HPVTransF[0].CIN2reg30;
		HPVparamsLogL.out[SimCount2][37] = HPVTransF[0].CIN2prog30;
		HPVparamsLogL.out[SimCount2][38] = HPVTransF[0].CIN2prog50; 
	}
		// Read the std deviation of the study effect (HSIL prevalence or HSIL given abnormal)
		/*a = 10.92574;
		b = 16.52705;
		p = r[24];
		q = 1.0 - r[24];
		cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
		HPVtransitionCC.FPClogL.VarStudyEffectH = 0.0; //pow(x, 2.0);
		HPVtransitionCC.HouseholdLogL.VarStudyEffectH = 0.0; //pow(x, 2.0);
		HPVtransitionCC.NOARTlogL.VarStudyEffectH = 0.0; //pow(x, 2.0);
		HPVtransitionCC.ONARTlogL.VarStudyEffectH = 0.0; //pow(x, 2.0);

		//Sensitivity among HIV negative (AB)
		a = 20.617; //15.9; //
		b = 12.476; //29.5286; //
		p = r[25];
		q = 1.0 - r[25];
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
		//HPVtransitionCC.ExpSeL = x;
		HPVtransitionCC.ExpSeL = 0.45 + r[25]*(0.84 - 0.45);

		//Sensitivity among HIV negative (HSIL|AB)
		a = 18.966; 
		b = 6.664; 
		p = r[26];
		q = 1.0 - r[26];
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
		//HPVtransitionCC.ExpSeH_L = x;
		HPVtransitionCC.ExpSeH_L = 0.66 + r[26]*(0.9 - 0.66);

		//Specificity among HIV negative  (AB)*
		a = 345.041; //24.101; //
		b = 16.637; //0.926;  //
		p = r[27];
		q = 1.0 - r[27];
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
		//HPVtransitionCC.ExpSpL = x;
		HPVtransitionCC.ExpSpL = 0.82 + r[27]*(1.0 - 0.82);

		//Specificity among HIV negative (HSIL|AB)
		a = 6.895;
		b = 2.668;
		p = r[28];
		q = 1.0 - r[28];
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
		//HPVtransitionCC.ExpSpH_L = x;
		HPVtransitionCC.ExpSpH_L = 0.57 + r[28]*(0.88 - 0.57);

		//Sensitivity among HIV positive (AB)
		a = 11.647;
		b = 4.454;
		p = r[29];
		q = 1.0 - r[29];
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
		//HPVtransitionCC.ExpSeHIVL = x;
		HPVtransitionCC.ExpSeHIVL = 0.45 + r[29]*(0.84 - 0.45);
		
		//Sensitivity among HIV positive (HSIL|AB)
		a = 11.242;
		b = 2.935;
		p = r[30];
		q = 1.0 - r[30];
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
		//HPVtransitionCC.ExpSeHIVH_L = x;
		HPVtransitionCC.ExpSeHIVH_L = 0.66 + r[30]*(0.9 - 0.66);
		
		//Specificity among HIV positive (AB)
		a = 15.525;
		b = 2.744;
		p = r[31];
		q = 1.0 - r[31];
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
		//HPVtransitionCC.ExpSpHIVL = x;
		HPVtransitionCC.ExpSpHIVL = 0.82 + r[31]*(1.0 - 0.82);

		//Specificity among HIV positive (HSIL|AB) 
		a = 18.575;
		b = 6.028;
		p = r[32];
		q = 1.0 - r[32];
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
		//HPVtransitionCC.ExpSpHIVH_L = x;
		HPVtransitionCC.ExpSpHIVH_L = 0.57 + r[32]*(0.88 - 0.57);
		
		//Sensitivity among HIV positive ART (AB)
		a = 11.647;
		b = 4.454;
		p = r[33];
		q = 1.0 - r[33];
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
		//HPVtransitionCC.ExpSeARTL = x;
		HPVtransitionCC.ExpSeARTL = 0.45 + r[33]*(0.84 - 0.45);

		//Sensitivity among HIV positive ART (HSIL|AB)
		a = 11.242;
		b = 2.935;
		p = r[34];
		q = 1.0 - r[34];
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
		//HPVtransitionCC.ExpSeARTH_L = x;
		HPVtransitionCC.ExpSeARTH_L = 0.66 + r[34]*(0.9 - 0.66);

		//Specificity among HIV positive ART (AB)
		a = 15.525;
		b = 2.744;
		p = r[35];
		q = 1.0 - r[35];
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
		//HPVtransitionCC.ExpSpARTL = x;
		HPVtransitionCC.ExpSpARTL = 0.82 + r[35]*(1.0 - 0.82);
		
		//Specificity among HIV positive ART (HSIL|AB) 
		a = 18.575;
		b = 6.028;
		p = r[36];
		q = 1.0 - r[36];
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
		//HPVtransitionCC.ExpSpARTH_L = x;
		HPVtransitionCC.ExpSpARTH_L = 0.57 + r[36]*(0.88 - 0.57);

		//Bias multiplier for WC data
		a = 0.0;
		b = 0.25;
		p = r[37];
		q = 1.0 - r[37];
		cdfnor(&ind, &p, &q, &x, &a, &b,  0, 0);
		HPVtransitionCC.WCBias = x;
		
		for (xx = 0; xx < 13; xx++){ HPVTransF[xx].ReactPro =  1.0;}
		 //r[38];}
		
	//cout << HPVtransitionCC.WCBias << endl;
	if (CurrSim == (SimCount2 * IterationsPerPC +1)){
		HPVparamsLogL.out[SimCount2][1] = HPVTransF[0].durMult;
		HPVparamsLogL.out[SimCount2][2] = HPVTransF[1].durMult;
		HPVparamsLogL.out[SimCount2][3] = HPVTransF[2].durMult;
		HPVparamsLogL.out[SimCount2][4] = HPVTransF[0].prop1pro;
		HPVparamsLogL.out[SimCount2][5] = HPVTransF[1].prop1pro;
		HPVparamsLogL.out[SimCount2][6] = HPVTransF[2].prop1pro;
		HPVparamsLogL.out[SimCount2][7] = HPVTransF[0].prog_fromCIN1;
		HPVparamsLogL.out[SimCount2][8] = HPVTransF[1].prog_fromCIN1;
		HPVparamsLogL.out[SimCount2][9] = HPVTransF[2].prog_fromCIN1;
		HPVparamsLogL.out[SimCount2][10] = HPVTransF[0].prog_fromCIN2;
		HPVparamsLogL.out[SimCount2][11] = HPVTransF[1].prog_fromCIN2;
		HPVparamsLogL.out[SimCount2][12] = HPVTransF[2].prog_fromCIN2;
		HPVparamsLogL.out[SimCount2][13] = HPVTransF[0].reg_fromCIN1;
		HPVparamsLogL.out[SimCount2][14] = HPVTransF[1].reg_fromCIN1;
		HPVparamsLogL.out[SimCount2][15] = HPVTransF[2].reg_fromCIN1;
		HPVparamsLogL.out[SimCount2][16] = HPVTransF[0].reg_fromCIN2;
		HPVparamsLogL.out[SimCount2][17] = HPVTransF[1].reg_fromCIN2;
		HPVparamsLogL.out[SimCount2][18] = HPVTransF[2].reg_fromCIN2;
		HPVparamsLogL.out[SimCount2][19] = HPVTransF[0].prop_reg_cl;
		HPVparamsLogL.out[SimCount2][20] = HPVTransF[0].CIN3HIV;
		HPVparamsLogL.out[SimCount2][21] = HPVTransF[0].AveDuration[3];
		HPVparamsLogL.out[SimCount2][22] = HPVTransF[0].CIN1_HIV;
		HPVparamsLogL.out[SimCount2][23] = HPVTransF[0].CIN2_HIV;
		HPVparamsLogL.out[SimCount2][24] = HPVTransF[0].prog_ART;
		HPVparamsLogL.out[SimCount2][25] = HPVTransF[0].reg_HIV;
		HPVparamsLogL.out[SimCount2][26] = HPVtransitionCC.HouseholdLogL.VarStudyEffectAB;
		HPVparamsLogL.out[SimCount2][27] = HPVtransitionCC.HouseholdLogL.VarStudyEffectH;
		HPVparamsLogL.out[SimCount2][28] = HPVtransitionCC.ExpSeL;
		HPVparamsLogL.out[SimCount2][29] = HPVtransitionCC.ExpSeH_L;
		HPVparamsLogL.out[SimCount2][30] = HPVtransitionCC.ExpSpL;
		HPVparamsLogL.out[SimCount2][31] = HPVtransitionCC.ExpSpH_L;
		HPVparamsLogL.out[SimCount2][32] = HPVtransitionCC.ExpSeHIVL;
		HPVparamsLogL.out[SimCount2][33] = HPVtransitionCC.ExpSeHIVH_L;
		HPVparamsLogL.out[SimCount2][34] = HPVtransitionCC.ExpSpHIVL;
		HPVparamsLogL.out[SimCount2][35] = HPVtransitionCC.ExpSpHIVH_L;
		HPVparamsLogL.out[SimCount2][36] = HPVtransitionCC.ExpSeARTL;
		HPVparamsLogL.out[SimCount2][37] = HPVtransitionCC.ExpSeARTH_L;
		HPVparamsLogL.out[SimCount2][38] = HPVtransitionCC.ExpSpARTL;
		HPVparamsLogL.out[SimCount2][39] = HPVtransitionCC.ExpSpARTH_L;
		HPVparamsLogL.out[SimCount2][40] = HPVtransitionCC.WCBias;
	}
	*/
	//cout << HPVtransitionCC.ExpSpARTH_L << endl;
	
	/*if (OneType == 0 && targets == 0){
		//if (CreateCohort == 1 && UseMedians == 1){
		if (UseMedians == 1){
			ReadHPVparams();
		}
		else{
			int xx;
			SimulateHPVparams16();
			SimulateHPVparams18();

			for (xx = 2; xx < 13; xx++){
				SimulateHPVparamsOTHER(xx);
			}
		}
	}
	else{
		if (WhichType == 0){ SimulateHPVparams16(); }
		else {
			if (WhichType == 1){ SimulateHPVparams18(); }
			else{ SimulateHPVparamsOTHER(WhichType); }
		}
	}*/

}

void Indiv::GetRiskGroup()
{
		
	int RiskPrimary;
	int Risk2ndary;
	int ic, ix;
	
	if (VirginInd == 1){
		if (RiskGroup == 1){ risk_for_start = 0; }
		else{ risk_for_start = 1; }
	}
	else{
		if (RiskGroup == 1){ risk_for_start = 2; }
		else{ risk_for_start = 14; }
		if (CurrPartners == 1){
			risk_for_start += 1;
			ic = IDprimary;
			RiskPrimary = Register[ic - 1].RiskGroup;
			if (MarriedInd == 1){
				risk_for_start += 2;
			}
			risk_for_start += RiskPrimary - 1;
		}
		if (CurrPartners == 2){
			risk_for_start += 5;
			ic = IDprimary;
			ix = ID2ndary;
			RiskPrimary = Register[ic - 1].RiskGroup;
			Risk2ndary = Register[ix - 1].RiskGroup;
			if (MarriedInd == 0){
				risk_for_start += RiskPrimary + Risk2ndary - 2;
			}
			else{
				risk_for_start += (RiskPrimary - 1) * 2 + Risk2ndary + 2;
			}
		}
		if (FSWind == 1){ risk_for_start = 19; }
	}
}

void ProspectiveCohort(int cycle)
{
	
		int tss = Register.size();
		int SimCount2 = (CurrSim - 1)/IterationsPerPC;

		ostringstream s;

		if (process_num >0){
			s << process_num << "_" << SimCount2 << "ProspectiveCohort.txt";
		}
		else{
			s << SimCount2 << "ProspectiveCohort.txt";
		}
				
		ofstream file(s.str().c_str(), std::ios::app); 

		for (int i = 0; i < tss; i++){
			if (Register[i].VirginInd == 0 && Register[i].AliveInd == 1 && Register[i].AgeExact >= 18 &&  Register[i].AgeExact < 25 && Register[i].HIVstage == 0 && Register[i].SexInd==1 ){// && Register[i].VirginInd==0 ){  //&& Register[i].TrueStage>2
				//){   
				file << CurrSim << " " << CurrYear << " " <<  i+1 << " " <<  Register[i].AgeExact << " ";
				//file << Register[i].SexInd << " "; // << Register[i].VirginInd << " " ; 
				for (int xx = 0; xx < 13; xx++){
					file << Register[i].HPVstage[xx] << " " ;
				}
				//file << Register[i].TrueStage << " " << Register[i].HIVstage << " "; // << Register[i].RiskGroup <<endl;//" " ; //
				//file << Register[i].VaccinationStatus[0] << " " << Register[i].VaccinationStatus[1] << " " ; 
			
				file << endl;
					
			}
		}
		
		file.close();	
}

void ReadHPVparams()
{
	ifstream file;
	stringstream s;
	s << "MediansHPV.txt";
	string path = "./input/" + s.str();
	
	file.open(path.c_str());
	//Read transmission prob M
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransM[xx].TransmProb; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	//Read transmission prob F
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].TransmProb; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	//Read the average duration of HPV DNA positive infection for males
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransM[xx].AveDuration[0]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	//Read the average duration of HPV DNA positive infection for females
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].AveDuration[0]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	//Read the average duration of latency for males
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransM[xx].AveDuration[5]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	//Read the average duration of latency for females
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].AveDuration[5]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	//Read the average duration of immunity for males
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransM[xx].AveDuration[6]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	//Read the average duration of immunity for females
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].AveDuration[6]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	//Read the proportion that will become latently infected after DNA clearance
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].propLatent; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');

	//Read the impact of HIV stage on reactivation of latency
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransM[xx].latentHIVreact; }
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransM[xx].lateHIVreact; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	//Read the impact of HIV stage on HPV infection duration
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransM[xx].AcuteLateHIVclear; }
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransM[xx].LatentHIVclear; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	//Read the std deviation of the study effect
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransM[xx].HouseholdLogL.VarStudyEffect; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].durMult;}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].prop1pro; }
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].prog_fromCIN1; }
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].reg_fromCIN1; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].prog_fromCIN2; }
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].reg_fromCIN2; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].prop_reg_cl; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].CIN1_HIV ;}
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].CIN2_HIV ;}
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].prog_ART ;}
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].reg_HIV ;}
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].reg_ART ;}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].CIN3HIV ; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].AveDuration[3]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].CIN3shape ; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].CIN3HIV_ART ;}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HPVtransitionCC.FPClogL.VarStudyEffectAB;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HPVtransitionCC.FPClogL.VarStudyEffectH;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].diag_stageI;}
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].diag_stageII;}
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].diag_stageIII;}
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].diag_stageIV;}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].CIN2reg30;}
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].CIN2prog30;}
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].CIN2prog50;}	
	for (int xx = 0; xx < 13; xx++){ file >> HPVTransF[xx].CIN3_50;}	
	
	file.close();
	//cout << HPVTransF[12].CIN2reg30 << " " <<HPVTransF[12].CIN2prog30 << " " <<HPVTransF[12].CIN2prog50 << endl;
	for (int xx = 0; xx < 13; xx++){ 
		HPVTransM[xx].propLatent = HPVTransF[xx].propLatent; 
		HPVTransF[xx].latentHIVreact = HPVTransM[xx].latentHIVreact;
		HPVTransF[xx].lateHIVreact =  HPVTransM[xx].lateHIVreact;
		HPVTransF[xx].AcuteLateHIVclear = HPVTransM[xx].AcuteLateHIVclear;
		HPVTransF[xx].LatentHIVclear = HPVTransM[xx].LatentHIVclear;
		HPVTransF[xx].HouseholdLogL.VarStudyEffect = HPVTransM[xx].HouseholdLogL.VarStudyEffect;
		HPVTransF[xx].NOARTlogL.VarStudyEffect = HPVTransM[xx].HouseholdLogL.VarStudyEffect;
		HPVTransF[xx].FPClogL.VarStudyEffect = HPVTransM[xx].HouseholdLogL.VarStudyEffect;
		//HPVTransF[xx].AveDuration[3] = HPVTransF[xx].AveDuration[3] * 52.0; 
		//HPVTransF[xx].CIN3shape = 2.5 ;
		HPVTransF[xx].AveDuration[4] = 0.5 * 52.0 ; //3.0 * 52.0; Does not get used!!
		HPVTransF[xx].prog_fromHPV = (1.0/(HPVTransF[xx].durMult*HPVTransF[xx].AveDuration[0]))*HPVTransF[xx].prop1pro; 
		HPVTransF[xx].reg_fromHPV = (1.0/(HPVTransF[xx].durMult*HPVTransF[xx].AveDuration[0]))*(1.0-HPVTransF[xx].prop1pro);
		HPVTransF[xx].CIN1_ART = max(HPVTransF[xx].CIN1_HIV * HPVTransF[xx].prog_ART, 1.0);
		HPVTransF[xx].CIN2_ART = max(HPVTransF[xx].CIN2_HIV * HPVTransF[xx].prog_ART, 1.0);
		HPVTransF[xx].CIN3HIV = 1.0/HPVTransF[xx].CIN2_HIV ;
		
	}
	HPVtransitionCC.HouseholdLogL.VarStudyEffectAB = HPVtransitionCC.FPClogL.VarStudyEffectAB;
	HPVtransitionCC.NOARTlogL.VarStudyEffectAB = HPVtransitionCC.FPClogL.VarStudyEffectAB;
	HPVtransitionCC.ONARTlogL.VarStudyEffectAB = HPVtransitionCC.FPClogL.VarStudyEffectAB;
	HPVtransitionCC.HouseholdLogL.VarStudyEffectH = HPVtransitionCC.FPClogL.VarStudyEffectH;
	HPVtransitionCC.NOARTlogL.VarStudyEffectH = HPVtransitionCC.FPClogL.VarStudyEffectH;
	HPVtransitionCC.ONARTlogL.VarStudyEffectH = HPVtransitionCC.FPClogL.VarStudyEffectH;

}

void Pop::GetCurrHPVprev(int WhichType)
{
	
	int TotPop,  TotPopART, TotPopHIVneg, TotPopHIVpos, TotPopHIVneg1to2, TotPopHIVpos1to2;
	int TotPopHIVneg1to11, TotPopHIVpos1to11, TotPop1to11;
	int	TotHPV,  TotHPVART, TotHPVHIVneg, TotHPVHIVpos, TotHPVHIVneg1to2, TotHPVHIVpos1to2;
	int TotHPVHIVneg1to11, TotHPVHIVpos1to11, TotHPV1to11;
	int ic;

	TotPop = 0; TotHPV = 0;
	TotPopART = 0; TotHPVART = 0;
	TotPopHIVneg = 0; TotPopHIVpos = 0; TotPopHIVneg1to2 = 0; TotPopHIVpos1to2 = 0; 
	TotPopHIVneg1to11=0; TotPopHIVpos1to11=0; TotPop1to11=0;
	TotHPVHIVneg = 0; TotHPVHIVpos = 0; TotHPVHIVneg1to2 = 0; TotHPVHIVpos1to2 = 0;
	TotHPVHIVneg1to11=0; TotHPVHIVpos1to11=0; TotHPV1to11=0;
	
	int tpp = Register.size();

	for (ic = 0; ic<tpp; ic++){
		if(Register[ic].AliveInd == 1 && Register[ic].SexInd == 1){
			if (Register[ic].AgeGroup>2  &&  Register[ic].VirginInd == 0){
				TotPop += 1;
				if (Register[ic].HPVstage[WhichType] == 1||Register[ic].HPVstage[WhichType] == 2||
					Register[ic].HPVstage[WhichType] == 3||Register[ic].HPVstage[WhichType] == 4||Register[ic].HPVstage[WhichType] == 5 )
					{TotHPV += 1; }	
			}
			if (Register[ic].AgeGroup>2  && Register[ic].HIVstage == 0 && Register[ic].VirginInd == 0){
				TotPopHIVneg += 1;
				if (Register[ic].HPVstage[WhichType] == 1||Register[ic].HPVstage[WhichType] == 2||
					Register[ic].HPVstage[WhichType] == 3||Register[ic].HPVstage[WhichType] == 4||Register[ic].HPVstage[WhichType] == 5 )
					{ TotHPVHIVneg += 1; }
			}
			if (Register[ic].AgeGroup>2  && Register[ic].HIVstage > 0 && Register[ic].VirginInd == 0){
				TotPopHIVpos += 1;
				if (Register[ic].HPVstage[WhichType] == 1||Register[ic].HPVstage[WhichType] == 2||
					Register[ic].HPVstage[WhichType] == 3||Register[ic].HPVstage[WhichType] == 4||Register[ic].HPVstage[WhichType] == 5 )
					{ TotHPVHIVpos += 1; }
			}
			if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup<5  && Register[ic].HIVstage == 0 && Register[ic].VirginInd == 0){
				TotPopHIVneg1to2 += 1;
				if (Register[ic].HPVstage[WhichType] == 1||Register[ic].HPVstage[WhichType] == 2||
					Register[ic].HPVstage[WhichType] == 3||Register[ic].HPVstage[WhichType] == 4||Register[ic].HPVstage[WhichType] == 5 )
					{ TotHPVHIVneg1to2 += 1; }
			}
			if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup<5  && Register[ic].HIVstage > 0 && Register[ic].VirginInd == 0){
				TotPopHIVpos1to2 += 1;
				if (Register[ic].HPVstage[WhichType] == 1||Register[ic].HPVstage[WhichType] == 2||
					Register[ic].HPVstage[WhichType] == 3||Register[ic].HPVstage[WhichType] == 4||Register[ic].HPVstage[WhichType] == 5 )
					{ TotHPVHIVpos1to2 += 1; }
			}
			if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup<13  && Register[ic].VirginInd == 0){
				TotPop1to11 += 1;
				if (Register[ic].HPVstage[WhichType] == 1 ||Register[ic].HPVstage[WhichType] == 2||
					Register[ic].HPVstage[WhichType] == 3||Register[ic].HPVstage[WhichType] == 4||Register[ic].HPVstage[WhichType] == 5 )
					{TotHPV1to11 += 1; }	
			}
			if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup<13  && Register[ic].HIVstage == 0 && Register[ic].VirginInd == 0){
				TotPopHIVneg1to11 += 1;
				if (Register[ic].HPVstage[WhichType] == 1||Register[ic].HPVstage[WhichType] == 2||
					Register[ic].HPVstage[WhichType] == 3||Register[ic].HPVstage[WhichType] == 4||Register[ic].HPVstage[WhichType] == 5 )
					{ TotHPVHIVneg1to11 += 1; }
			}
			if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup<13  && Register[ic].HIVstage > 0 && Register[ic].VirginInd == 0){
				TotPopHIVpos1to11 += 1;
				if (Register[ic].HPVstage[WhichType] == 1||Register[ic].HPVstage[WhichType] == 2||
					Register[ic].HPVstage[WhichType] == 3||Register[ic].HPVstage[WhichType] == 4||Register[ic].HPVstage[WhichType] == 5 )
					{ TotHPVHIVpos1to11 += 1; }
			}
		}
	}
		HPVprev15to64HIVposF.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHPVHIVpos1to11/TotPopHIVpos1to11;
		HPVprev15to64HIVnegF.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHPVHIVneg1to11/TotPopHIVneg1to11;
		HPVprev15to24HIVposF.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHPVHIVpos1to2/TotPopHIVpos1to2;
		HPVprev15to24HIVnegF.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHPVHIVneg1to2/TotPopHIVneg1to2;
		HPVprev15to64F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHPV1to11/TotPop1to11;
		HPVprev15to64ARTF.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHPVART/TotPopART;
		HPVprevFHIVpos.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHPVHIVpos/TotPopHIVpos;
		HPVprevFHIVneg.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHPVHIVneg/TotPopHIVneg;
		HPVprevF.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHPV/TotPop;
}

void Pop::GetHHprev(STDtransition* a, int STDind) 
{
	int ic, ii, xx;
	double	ia;
	double numerator, denominator;
	int tpp = Register.size();
	
	if (a->HouseholdLogL.Observations>0){
		for (ic = 0; ic<a->HouseholdLogL.Observations; ic++){
			if (a->HouseholdLogL.StudyYear[ic] == CurrYear){
				numerator = 0.0;
				denominator = 0.0;
				for (ii = 0; ii < tpp; ii++){
					ia = Register[ii].AgeExact;
					
					//virgins were not excluded and HIV not tested
					if (ia >= (a->HouseholdLogL.AgeStart[ic]) && ia < (a->HouseholdLogL.AgeEnd[ic]) &&
						Register[ii].AliveInd == 1 && Register[ii].SexInd == a->SexInd && a->HouseholdLogL.ExclVirgins[ic] == 0 && 
						a->HouseholdLogL.HIVprevInd[ic] == 0 )
					{
						if(a->HouseholdLogL.Cyt_cat[ic] == 0){
							denominator += 1.0;
							if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3 ))
									{ numerator += 1.0; }	
						}
						else if(a->HouseholdLogL.Cyt_cat[ic] == 2){
							denominator += 1.0;
							if (STDind == 9 && (Register[ii].TrueStage==2||Register[ii].TrueStage==3 ))
									{ numerator += 1.0; }	
						}
						else if(a->HouseholdLogL.Cyt_cat[ic] == 1){
							if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3))
								{ denominator += 1.0;
								if(Register[ii].TrueStage==2||Register[ii].TrueStage==3)
								{ numerator += 1.0; }}
						}
						else if(a->HouseholdLogL.Cyt_cat[ic] ==3){
								denominator += 1.0;
								if(Register[ii].HPVstage[0]==1||Register[ii].HPVstage[0]==2||Register[ii].HPVstage[0]==3||Register[ii].HPVstage[0]==4||Register[ii].HPVstage[0]==5)
								{ numerator += 1.0; }
						}
						else if(a->HouseholdLogL.Cyt_cat[ic] ==4){
								denominator += 1.0;
								if(Register[ii].HPVstage[1]==1||Register[ii].HPVstage[1]==2||Register[ii].HPVstage[1]==3||Register[ii].HPVstage[1]==4||Register[ii].HPVstage[1]==5)
								{ numerator += 1.0; }
						}
					}	
					//virgins were excluded, HIV was tested and all HIV negative
					if (ia >=(a->HouseholdLogL.AgeStart[ic]) && ia < (a->HouseholdLogL.AgeEnd[ic]) &&
						Register[ii].AliveInd == 1 && Register[ii].SexInd == a->SexInd &&
						a->HouseholdLogL.ExclVirgins[ic] == 1 && Register[ii].VirginInd == 0 &&
						Register[ii].HIVstage == 0 && a->HouseholdLogL.HIVprevInd[ic] == 1 && a->HouseholdLogL.HIVprev[ic] == 0 )
					{
						if(a->HouseholdLogL.Cyt_cat[ic] == 0){
							denominator += 1.0;
							if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3 ))
									{ numerator += 1.0; }	
						}
						else if(a->HouseholdLogL.Cyt_cat[ic] == 2){
							denominator += 1.0;
							if (STDind == 9 && (Register[ii].TrueStage==2||Register[ii].TrueStage==3 ))
									{ numerator += 1.0; }	
						}
						else if(a->HouseholdLogL.Cyt_cat[ic] == 1){
							if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3))
								{ denominator += 1.0;
								if(Register[ii].TrueStage==2||Register[ii].TrueStage==3)
								{ numerator += 1.0; }}
						}
						else {
							for(xx = 0; xx < 13; xx++){
								if(a->HouseholdLogL.Cyt_cat[ic] == xx+3){
									denominator += 1.0;
									if(Register[ii].HPVstage[xx]==1||Register[ii].HPVstage[xx]==2||Register[ii].HPVstage[xx]==3||Register[ii].HPVstage[xx]==4||Register[ii].HPVstage[xx]==5)
									{ numerator += 1.0; }
								}
							}		
						}
					}
					//virgins were not excluded, HIV was tested and all HIV negative (MacDonald)
					if (ia >=(a->HouseholdLogL.AgeStart[ic]) && ia < (a->HouseholdLogL.AgeEnd[ic]) &&
						Register[ii].AliveInd == 1 && Register[ii].SexInd == a->SexInd && a->HouseholdLogL.ExclVirgins[ic] == 0 &&
						Register[ii].HIVstage == 0 && a->HouseholdLogL.HIVprevInd[ic] == 1 && a->HouseholdLogL.HIVprev[ic] == 0 )
						{
							if(a->HouseholdLogL.Cyt_cat[ic] == 0){
								denominator += 1.0;
								if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3 ))
									{ numerator += 1.0; }	
							}
							else if(a->HouseholdLogL.Cyt_cat[ic] == 1){
								if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3))
									{ denominator += 1.0;  
									 if(Register[ii].TrueStage==2||Register[ii].TrueStage==3)
										{ numerator += 1.0; }
									}
							}
							else if(a->HouseholdLogL.Cyt_cat[ic] == 2){
								denominator += 1.0;
								if (STDind == 9 && (Register[ii].TrueStage==2||Register[ii].TrueStage==3 ))
									{ numerator += 1.0; }
							}
							else {
								for(xx = 0; xx < 13; xx++){
									if(a->HouseholdLogL.Cyt_cat[ic] == xx+3){
										denominator += 1.0;
										if(Register[ii].HPVstage[xx]==1||Register[ii].HPVstage[xx]==2||Register[ii].HPVstage[xx]==3||Register[ii].HPVstage[xx]==4||Register[ii].HPVstage[xx]==5)
										{ numerator += 1.0; }
									}
								}		
							}
						}
					//virgins were excluded, HIV was tested and all HIV positive 
					if (ia >=(a->HouseholdLogL.AgeStart[ic]) && ia < (a->HouseholdLogL.AgeEnd[ic]) &&
						Register[ii].AliveInd == 1 && Register[ii].SexInd == a->SexInd &&
						a->HouseholdLogL.ExclVirgins[ic] == 1 && Register[ii].VirginInd == 0 &&
						Register[ii].HIVstage > 0 && a->HouseholdLogL.HIVprevInd[ic] == 1 && a->HouseholdLogL.HIVprev[ic] == 1 )
						{
							if(a->HouseholdLogL.Cyt_cat[ic] == 0){
							denominator += 1.0;
							if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3 ))
									{ numerator += 1.0; }	
							}
							else if(a->HouseholdLogL.Cyt_cat[ic] == 2){
							denominator += 1.0;
							if (STDind == 9 && (Register[ii].TrueStage==2||Register[ii].TrueStage==3 ))
									{ numerator += 1.0; }	
							}
							else if(a->HouseholdLogL.Cyt_cat[ic] == 1){
								if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3))
									{ denominator += 1.0; 
								if(Register[ii].TrueStage==2||Register[ii].TrueStage==3)
									{ numerator += 1.0; }}
							}
							else {
								for(xx = 0; xx < 13; xx++){
									if(a->HouseholdLogL.Cyt_cat[ic] == xx+3){
										denominator += 1.0;
										if(Register[ii].HPVstage[xx]==1||Register[ii].HPVstage[xx]==2||Register[ii].HPVstage[xx]==3||Register[ii].HPVstage[xx]==4||Register[ii].HPVstage[xx]==5)
										{ numerator += 1.0; }
									}
								}		
							}
						}
					//virgins were not excluded, HIV was tested and all HIV positive (MacDonald)
					if (ia >=(a->HouseholdLogL.AgeStart[ic]) && ia < (a->HouseholdLogL.AgeEnd[ic]) &&
						Register[ii].AliveInd == 1 && Register[ii].SexInd == a->SexInd && a->HouseholdLogL.ExclVirgins[ic] == 0 &&
						Register[ii].HIVstage > 0 && a->HouseholdLogL.HIVprevInd[ic] == 1 && a->HouseholdLogL.HIVprev[ic] == 1 )
						{
							if(a->HouseholdLogL.Cyt_cat[ic] == 0){
								denominator += 1.0;
								if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3 ))
									{ numerator += 1.0; }	
							}
							else if(a->HouseholdLogL.Cyt_cat[ic] == 1){
								if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3))
									{ denominator += 1.0;  
									 if(Register[ii].TrueStage==2||Register[ii].TrueStage==3)
										{ numerator += 1.0; }
									}
							}
							else if(a->HouseholdLogL.Cyt_cat[ic] == 2){
								denominator += 1.0;
								if (STDind == 9 && (Register[ii].TrueStage==2||Register[ii].TrueStage==3))
									{ numerator += 1.0; }
							}
							else {
								for(xx = 0; xx < 13; xx++){
									if(a->HouseholdLogL.Cyt_cat[ic] == xx+3){
										denominator += 1.0;
										if(Register[ii].HPVstage[xx]==1||Register[ii].HPVstage[xx]==2||Register[ii].HPVstage[xx]==3||Register[ii].HPVstage[xx]==4||Register[ii].HPVstage[xx]==5)
										{ numerator += 1.0; }
									}
								}		
							}
						}
				}
				//a->HouseholdLogL.ModelPrev[ic] = 1.0*numerator / denominator;
				a->HouseholdLogL.ModelNum[ic] += numerator; 
				a->HouseholdLogL.ModelDenom[ic] += denominator;
			}	
		}
	}
}
void Pop::GetHHNprev(STDtransition* a, int STDind) 
{
	int ic, ii;
	double	ia;
	double numerator, denominator;
	int tpp = Register.size();
	if (a->HouseholdNlogL.Observations>0){
		for (ic = 0; ic<a->HouseholdNlogL.Observations; ic++){
			if (a->HouseholdNlogL.StudyYear[ic] == CurrYear){
				numerator = 0.0;
				denominator = 0.0;
				for (ii = 0; ii < tpp; ii++){
					ia = Register[ii].AgeExact;
					
					//virgins were not excluded and HIV not tested
					if (ia >= (a->HouseholdNlogL.AgeStart[ic]) && ia < (a->HouseholdNlogL.AgeEnd[ic]) &&
						Register[ii].AliveInd == 1 && Register[ii].SexInd == a->SexInd && a->HouseholdNlogL.ExclVirgins[ic] == 0 && 
						a->HouseholdNlogL.HIVprevInd[ic] == 0 )
					{
						if(a->HouseholdNlogL.Cyt_cat[ic] == 0){
							denominator += 1.0;
							if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3 ))
									{ numerator += 1.0; }	
						}
						else{
							if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3))
								{ denominator += 1.0;
								if(Register[ii].TrueStage==2||Register[ii].TrueStage==3)
								{ numerator += 1.0; }}
						}
					}	
		
				}
				a->HouseholdNlogL.ModelPrev[ic] = 1.0*numerator / denominator;
			}
		}
	}
}
void Pop::GetFPCprev(STDtransition* a, int STDind)
{
	int ic, ii, ia, xx;
	double numerator, denominator, TempFP;
	
	if (a->FPClogL.Observations>0){
		
		for (ic = 0; ic < a->FPClogL.Observations; ic++){
			if (a->FPClogL.StudyYear[ic] == CurrYear){
				numerator = 0.0;
				denominator = 0.0;
				int tpp = Register.size();
				for (ii = 0; ii < tpp; ii++){
				
					if (Register[ii].AgeGroup>1 && Register[ii].AliveInd == 1 &&
						Register[ii].SexInd == 1 && 
						Register[ii].HIVstage == 0 && a->FPClogL.HIVprevInd[ic] == 1 && a->FPClogL.HIVprev[ic] == 0){
							ia = Register[ii].AgeGroup - 2;
							TempFP = FPCweights[ia];
							if(a->FPClogL.Cyt_cat[ic] == 0){
								denominator += TempFP;
								if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3  ))
									{ numerator += TempFP; }		
							}
							else if(a->FPClogL.Cyt_cat[ic] == 1){
								if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3  ))
									{ denominator += TempFP; 
								if(Register[ii].TrueStage==2||Register[ii].TrueStage==3)
									{ numerator += TempFP; }}
							}
						
					}
					if (Register[ii].AgeGroup>1 && Register[ii].AliveInd == 1 &&
						Register[ii].SexInd == 1 && 
						Register[ii].HIVstage >0  && a->FPClogL.HIVprevInd[ic] == 1 && a->FPClogL.HIVprev[ic] == 1){
							ia = Register[ii].AgeGroup - 2;
							TempFP = FPCweights[ia];
							if(a->FPClogL.Cyt_cat[ic] == 0){
								denominator += TempFP;
								if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3  ))
									{ numerator += TempFP; }		
							}
							else if(a->FPClogL.Cyt_cat[ic] == 1){
								if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3  ))
									{ denominator += TempFP; 
								if(Register[ii].TrueStage==2||Register[ii].TrueStage==3)
									{ numerator += TempFP; }}
							}
					}
					if (Register[ii].AgeGroup>1 && Register[ii].AliveInd == 1 &&
						Register[ii].SexInd == 1 && 
						a->FPClogL.HIVprevInd[ic] == 0 && a->FPClogL.HIVprev[ic] == 0){
							if(a->FPClogL.Cyt_cat[ic] > 2){
								for(xx = 0; xx < 13; xx++){
									if(a->FPClogL.Cyt_cat[ic] == xx+3){
										denominator += 1.0;
										if(Register[ii].HPVstage[xx]==1||Register[ii].HPVstage[xx]==2||Register[ii].HPVstage[xx]==3||Register[ii].HPVstage[xx]==4||Register[ii].HPVstage[xx]==5)
										{ numerator += 1.0; }
									}
								}
							}		
					}
				}
				
				a->FPClogL.ModelNum[ic] += numerator; 
				a->FPClogL.ModelDenom[ic] += denominator;
			}
		}
	}
}
void Pop::GetNOARTprev(STDtransition* a, int STDind)
{
	int ic, ii, xx;
	double ia;
	double numerator, denominator;
	int tpp = Register.size();

	if (a->NOARTlogL.Observations>0){
		for (ic = 0; ic < a->NOARTlogL.Observations; ic++){
			

			if (a->NOARTlogL.StudyYear[ic] == CurrYear){
				numerator = 0.0;
				denominator = 0.0;
				for (ii = 0; ii < tpp; ii++){
					ia = Register[ii].AgeExact;
					if (ia >=(a->NOARTlogL.AgeStart[ic] - 1) && ia <= (a->NOARTlogL.AgeEnd[ic] + 1)&&
						Register[ii].AliveInd == 1 && Register[ii].SexInd == a->SexInd && 
						(Register[ii].HIVstage == 3 || Register[ii].HIVstage == 4)){ //|| Register[ii].HIVstage == 6 || (Register[ii].HIVstage == 5 && Register[ii].ARTweeks<104))){ //Moodley + Soweto starting ART//
							if(a->NOARTlogL.Cyt_cat[ic] == 0){
								denominator += 1.0;
								if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3))
									{ numerator += 1.0; }
							}
							else if (a->NOARTlogL.Cyt_cat[ic] == 2){
								denominator += 1.0;							
								 if(Register[ii].TrueStage==2||Register[ii].TrueStage==3)
										{ numerator += 1.0; }
							}	
							else if (a->NOARTlogL.Cyt_cat[ic] == 1){
								if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3))
									{ 
										denominator += 1.0; 
									 if(Register[ii].TrueStage==2||Register[ii].TrueStage==3)
										{ numerator += 1.0; }}
							}
							else if(a->NOARTlogL.Cyt_cat[ic] > 2){
								for(xx = 0; xx < 13; xx++){
									if(a->NOARTlogL.Cyt_cat[ic] == xx+3){
										denominator += 1.0;
										if(Register[ii].HPVstage[xx]==1||Register[ii].HPVstage[xx]==2||Register[ii].HPVstage[xx]==3||Register[ii].HPVstage[xx]==4||Register[ii].HPVstage[xx]==5)
										{ numerator += 1.0; }
									}
								}
							}
					}
				}
				//a->NOARTlogL.ModelPrev[ic] = numerator / denominator;	
				a->NOARTlogL.ModelNum[ic] += numerator; 
				a->NOARTlogL.ModelDenom[ic] += denominator;
			}
		}
	}
}
void Pop::GetONARTprev(STDtransition* a, int STDind)
{
	int ic, ii, xx;
	double ia;
	double numerator, denominator;
	int tpp = Register.size();

	if (a->ONARTlogL.Observations>0){
		for (ic = 0; ic < a->ONARTlogL.Observations; ic++){
			if (a->ONARTlogL.StudyYear[ic] == CurrYear){
				numerator = 0.0;
				denominator = 0.0;
				for (ii = 0; ii < tpp; ii++){
					ia = Register[ii].AgeExact;
					if (ia >=(a->ONARTlogL.AgeStart[ic] - 1) && ia <= (a->ONARTlogL.AgeEnd[ic] + 1)&&
						Register[ii].AliveInd == 1 && Register[ii].SexInd == a->SexInd && 
						(Register[ii].HIVstage == 5 || Register[ii].HIVstage == 6)){ 
							if(a->ONARTlogL.Cyt_cat[ic] == 0){
								denominator += 1.0;
								if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3))
									{ numerator += 1.0; }
							}
							else if(a->ONARTlogL.Cyt_cat[ic] == 2){
								denominator += 1.0;
								
								if (STDind == 9 && (Register[ii].TrueStage==2||Register[ii].TrueStage==3))
									{ numerator += 1.0; }
							}
							else if(a->ONARTlogL.Cyt_cat[ic] == 1){
								if (STDind == 9 && (Register[ii].TrueStage==1||Register[ii].TrueStage==2||Register[ii].TrueStage==3))
									{ denominator += 1.0; 
								 if(Register[ii].TrueStage==2||Register[ii].TrueStage==3)
									{ numerator += 1.0; }}
							}
							else if(a->ONARTlogL.Cyt_cat[ic] > 2){
								for(xx = 0; xx < 13; xx++){
									if(a->ONARTlogL.Cyt_cat[ic] == xx+3){
										denominator += 1.0;
										if(Register[ii].HPVstage[xx]==1||Register[ii].HPVstage[xx]==2||Register[ii].HPVstage[xx]==3||Register[ii].HPVstage[xx]==4||Register[ii].HPVstage[xx]==5)
										{ numerator += 1.0; }
									}
								}
							}
					}
				}
				//a->ONARTlogL.ModelPrev[ic] = 1.0*numerator / denominator;	
				a->ONARTlogL.ModelNum[ic] += numerator; 
				a->ONARTlogL.ModelDenom[ic] += denominator;	
			}
		}
	}
}

void Pop::AssignVacc2024()
{
	int ic, SimCount2, zz;
	double adjustedTxVEfficacyD1[13];
	double adjustedTxVEfficacyD2[13];
	double rcatch[MaxPop];
	double wane1[MaxPop];
	double cross[MaxPop];
	double Rloss1[MaxPop];
	double EfficacyD1Rand1[MaxPop];
	double EfficacyD2Rand1[MaxPop];
	double RandAcceptTxV1[MaxPop];
	int seedy;
	seedy = CurrSim * 92 + process_num * 7928 + CurrYear;
	SimCount2 = (CurrSim - 1) / IterationsPerPC;
	CRandomMersenne rg2(seedy);
	int tpp = Register.size();

	for (ic = 0; ic < tpp; ic++)
	{
		rcatch[ic] = rg2.Random();
		wane1[ic] = rg2.Random();
		cross[ic] = rg2.Random();
		RandAcceptTxV1[ic] = rg2.Random();
		Rloss1[ic] = rg2.Random();
		EfficacyD1Rand1[ic] = rg2.Random();
		EfficacyD2Rand1[ic] = rg2.Random();
	}
	// ofstream file("VaccDur.txt", std::ios::app);

	for (ic = 0; ic < tpp; ic++)
	{
		if (Register[ic].HIVstage == 0) {zz=0;}
			else if (Register[ic].HIVstage == 5){zz=1;}
			else {zz = 2;}

			if (AdministerMassTxV == 1 &&
				(Register[ic].AgeExact >= MassTxVAgeMIN &&
				Register[ic].AgeExact < MassTxVAgeMAX) &&
				Register[ic].AliveInd == 1 &&
				Register[ic].SexInd == 1)
			{
				if (RandAcceptTxV1[ic] < 0.5) //&& ((CatchUpVacc == 1 && Register[ic].GotTxV == 0) || (CatchUpVacc == 0)))
				{
				RSApop.NewTxV[zz*18 + Register[ic].AgeGroup][CurrYear - StartYear] += 1;
				Register[ic].GotTxV = 1;
				if (CatchUpVacc ==1 ){
					RSApop.NewVACC[18 * Register[ic].SexInd + Register[ic].AgeGroup][CurrYear - StartYear] += 1;
					Register[ic].GotVacc = 1;
					 }
					for (int yy = 0; yy < 13; yy++)
					{
						if (Register[ic].HPVstage[yy] > 2 && Register[ic].HPVstage[yy] <= 4)
						{
							adjustedTxVEfficacyD1[yy] = TxVEfficacyD1CIN[yy];
							adjustedTxVEfficacyD2[yy] = TxVEfficacyD2CIN[yy];
						}
						else
						{
							adjustedTxVEfficacyD1[yy] = TxVEfficacyD1[yy];
							adjustedTxVEfficacyD2[yy] = TxVEfficacyD2[yy];
						}
						if (!(Register[ic].HIVstage == 0 || Register[ic].HIVstage == 5))
						{
							adjustedTxVEfficacyD1[yy] *= ReductionFactorHIV;
							adjustedTxVEfficacyD2[yy] *= ReductionFactorHIV;
						}
						if (Register[ic].HIVstage == 5)
						{
							adjustedTxVEfficacyD1[yy] *= ReductionFactorHIVART;
							adjustedTxVEfficacyD2[yy] *= ReductionFactorHIVART;
						}
						
						if (EfficacyD1Rand1[ic] < adjustedTxVEfficacyD1[yy])
						{
							if ((Register[ic].HPVstage[yy] >= 1 && Register[ic].HPVstage[yy] <= 4) || Register[ic].HPVstage[yy] == 6)
							{
								 Register[ic].HPVstage[yy] = 0;
								 Register[ic].HPVstageE[yy] = 0;
								 if (CatchUpVacc == 1 ){
									Register[ic].VaccinationStatus[yy] = 1;
								 }
							}
						}
					}
					if (Rloss1[ic] < 0.7)
					{	RSApop.NewTxV[zz*18 + Register[ic].AgeGroup][CurrYear - StartYear] += 1;
						Register[ic].GotTxV = 2;
						if (CatchUpVacc ==1 ){
							RSApop.NewVACC[18 * Register[ic].SexInd + Register[ic].AgeGroup][CurrYear - StartYear] += 1;
							Register[ic].GotVacc = 1;
							}
						for (int yy = 0; yy < 13; yy++)
						{
							if (EfficacyD2Rand1[ic] < adjustedTxVEfficacyD2[yy])
							{
								if ((Register[ic].HPVstage[yy] >= 1 && Register[ic].HPVstage[yy] <= 4) || Register[ic].HPVstage[yy] == 6)
								{
									Register[ic].HPVstage[yy] = 0;
									Register[ic].HPVstageE[yy] = 0;
									if (CatchUpVacc == 1 ){
										Register[ic].VaccinationStatus[yy] = 1;
									 }
								}
							}
						}
					}
				}
			}
		/*
		if (CatchUpVacc == 1 && (Register[ic].AgeExact >= CatchUpAgeMIN && Register[ic].AgeExact < CatchUpAgeMAX ) &&
			Register[ic].AliveInd == 1 && Register[ic].GotVacc == 0 && Register[ic].SexInd == 1 && rcatch[ic] < CatchUpCoverage)
		{
			RSApop.NewVACC[18 * Register[ic].SexInd + Register[ic].AgeGroup][CurrYear - StartYear] += 1;
			Register[ic].GotVacc = 1;
			if (WHOvacc == 1)
			{
				Register[ic].VaccinationStatus[0] = 1;
				Register[ic].VaccinationStatus[1] = 1;
				Register[ic].VaccinationStatus[2] = 1;
				Register[ic].VaccinationStatus[3] = 1;
				Register[ic].VaccinationStatus[6] = 1;
				Register[ic].VaccinationStatus[8] = 1;
				Register[ic].VaccinationStatus[10] = 1;
			}
			else
			{
				Register[ic].VaccinationStatus[0] = 1;
				Register[ic].VaccinationStatus[1] = 1;
				if (cross[ic] < 0.5)
				{
					Register[ic].VaccinationStatus[2] = 1;
					Register[ic].VaccinationStatus[3] = 1;
					Register[ic].VaccinationStatus[6] = 1;
				}
			}
			if (VaccineWane == 1)
			{
				Register[ic].TimeVacc = 0;
				Register[ic].ExpVacc = 48 * 20 + (-48 * VaccDur * log(wane1[ic]));
				if (Register[ic].ExpVacc == 0)
				{
					Register[ic].ExpVacc = 1;
				}
				// file << ic << " " << Register[ic].ExpVacc << endl;
			}
		}*/ // commented out for TxV Implementation
		
		if (AdministerMassTxVtoART == 1 &&
			Register[ic].AgeExact >= MassTxVtoARTAgeMIN &&
			Register[ic].AgeExact < MassTxVtoARTAgeMAX &&
			Register[ic].AliveInd == 1 &&
			Register[ic].SexInd == 1 &&
			Register[ic].HIVstage == 5)
		{
			if (RandAcceptTxV1[ic] < 0.9) //&& ((CatchUpVaccHIV == 1 && Register[ic].GotTxV == 0) || (CatchUpVaccHIV == 0)))
			{ 
				RSApop.NewTxV[zz*18 + Register[ic].AgeGroup][CurrYear - StartYear] += 1;
				Register[ic].GotTxV = 1;
				if (CatchUpVacc ==1){
					RSApop.NewVACC[18 * Register[ic].SexInd + Register[ic].AgeGroup][CurrYear - StartYear] += 1;
					Register[ic].GotVacc = 1;
					}
				for (int yy = 0; yy < 13; yy++)
				{
					if (Register[ic].HPVstage[yy] > 2 && Register[ic].HPVstage[yy] <= 4)
					{
						adjustedTxVEfficacyD1[yy] = TxVEfficacyD1CIN[yy];
						adjustedTxVEfficacyD2[yy] = TxVEfficacyD2CIN[yy];
					}
					else
					{
						adjustedTxVEfficacyD1[yy] = TxVEfficacyD1[yy];
						adjustedTxVEfficacyD2[yy] = TxVEfficacyD2[yy];
					}
					if (EfficacyD1Rand1[ic] < adjustedTxVEfficacyD1[yy])
					{
						if ((Register[ic].HPVstage[yy] >= 1 && Register[ic].HPVstage[yy] <= 4) || Register[ic].HPVstage[yy] == 6 )//|| Register[ic].HPVstage[yy] == 7)
						{
							Register[ic].HPVstage[yy] = 0;
							Register[ic].HPVstageE[yy] = 0;
							if (CatchUpVacc ==1 ){
								Register[ic].VaccinationStatus[yy] = 1;
							 }
						}
					}
				}
				if (Rloss1[ic] < 0.7)
				{
					RSApop.NewTxV[zz*18 + Register[ic].AgeGroup][CurrYear - StartYear] += 1;
					Register[ic].GotTxV = 2;
					if (CatchUpVacc ==1 ){
									RSApop.NewVACC[18 * Register[ic].SexInd + Register[ic].AgeGroup][CurrYear - StartYear] += 1;
									Register[ic].GotVacc = 1;
					}
					for (int yy = 0; yy < 13; yy++)
					{
						if ((Register[ic].HPVstage[yy] >= 1 && Register[ic].HPVstage[yy] <= 4) || Register[ic].HPVstage[yy] == 6)
						{
							if (EfficacyD2Rand1[ic] < adjustedTxVEfficacyD2[yy])
							{
								Register[ic].HPVstage[yy] = 0;
								Register[ic].HPVstageE[yy] = 0;
								if (CatchUpVacc ==1 ){
									Register[ic].VaccinationStatus[yy] = 1;
								 }
							}
						}
					}
				}
			}
		}
	}
}

			void Pop::HitTargets(int WhichType)
			{

				double TotPop[20], TotHPV[20];
				int ic, ig, xx;

				for (ig = 0; ig < 20; ig++)
				{
					TotPop[ig] = 0;
					TotHPV[ig] = 0;
				}
				int tpp = Register.size();

				for (ic = 0; ic < tpp; ic++)
				{
					if (Register[ic].AliveInd == 1)
					{
						// McDonald 2014, F, HIVneg
						if (Register[ic].AgeExact >= 18 && Register[ic].AgeExact < 65 && Register[ic].SexInd == 1 && Register[ic].HIVstage == 0)
						{
							TotPop[0] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[0] += 1;
							}
						}
						// McDonald 2014, F, HIVpos
						if (Register[ic].AgeExact >= 18 && Register[ic].AgeExact < 65 && Register[ic].SexInd == 1 && Register[ic].HIVstage > 0)
						{
							TotPop[1] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[1] += 1;
							}
						}
						// Giuliano 2012, F, HIVneg
						if (Register[ic].AgeExact >= 15 && Register[ic].AgeExact < 25 && Register[ic].SexInd == 1 && Register[ic].HIVstage == 0 && Register[ic].VirginInd == 0)
						{
							TotPop[2] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[2] += 1;
							}
						}
						// Snyman 2011, F, not tested
						if (Register[ic].AgeExact >= 18 && Register[ic].AgeExact < 65 && Register[ic].SexInd == 1)
						{
							TotPop[3] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[3] += 1;
							}
						}
						// Snyman 2012, F, not tested
						if (Register[ic].AgeExact >= 18 && Register[ic].AgeExact < 65 && Register[ic].SexInd == 1)
						{
							TotPop[4] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[4] += 1;
							}
						}
						// Adler 2013, F, HIV neg
						if (Register[ic].AgeExact >= 17 && Register[ic].AgeExact < 22 && Register[ic].SexInd == 1 && Register[ic].HIVstage == 0 && Register[ic].VirginInd == 0)
						{
							TotPop[5] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[5] += 1;
							}
						}
						// Adler 2013, F, HIV pos
						if (Register[ic].AgeExact >= 17 && Register[ic].AgeExact < 22 && Register[ic].SexInd == 1 && Register[ic].HIVstage > 0 && Register[ic].VirginInd == 0)
						{
							TotPop[6] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[6] += 1;
							}
						}
						// Mubulawa 2015, F, HIV neg
						if (Register[ic].AgeExact >= 16 && Register[ic].AgeExact < 23 && Register[ic].SexInd == 1 && Register[ic].HIVstage == 0 && Register[ic].VirginInd == 0)
						{
							TotPop[7] += 1;
							TotPop[8] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[7] += 1;
								TotHPV[8] += 1;
							}
						}
						// Mbulawa 2006, F, HIVpos
						if (Register[ic].AgeExact >= 18 && Register[ic].AgeExact < 65 && Register[ic].SexInd == 1 && Register[ic].HIVstage > 0 && Register[ic].VirginInd == 0)
						{
							TotPop[9] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[9] += 1;
							}
						}
						// Mbulawa 2006, F, HIVneg
						if (Register[ic].AgeExact >= 18 && Register[ic].AgeExact < 65 && Register[ic].SexInd == 1 && Register[ic].HIVstage == 0 && Register[ic].VirginInd == 0)
						{
							TotPop[10] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[10] += 1;
							}
						}
						// Denny 2002, F, HIVpos
						if (Register[ic].AgeExact >= 18 && Register[ic].AgeExact < 53 && Register[ic].SexInd == 1 && Register[ic].HIVstage > 0 && Register[ic].VirginInd == 0)
						{
							TotPop[11] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[11] += 1;
							}
						}
						// CAPRISA 2007, F, HIV neg
						if (Register[ic].AgeExact >= 18 && Register[ic].AgeExact < 40 && Register[ic].SexInd == 1 && Register[ic].HIVstage == 0 && Register[ic].VirginInd == 0)
						{
							TotPop[12] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[12] += 1;
							}
						}
						// Vardas 2005, M, HIV neg
						if (Register[ic].AgeExact >= 16 && Register[ic].AgeExact < 25 && Register[ic].SexInd == 0 && Register[ic].HIVstage == 0 &&
							Register[ic].VirginInd == 0 && Register[ic].LifetimePartners < 6)
						{
							TotPop[13] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[13] += 1;
							}
						}
						// Mbulawa 2006, M, HIV pos
						if (Register[ic].AgeExact >= 18 && Register[ic].AgeExact < 65 && Register[ic].SexInd == 0 && Register[ic].HIVstage > 0 && Register[ic].VirginInd == 0)
						{
							TotPop[14] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[14] += 1;
							}
						}
						// Mbulawa 2006, M, HIV neg
						if (Register[ic].AgeExact >= 18 && Register[ic].AgeExact < 65 && Register[ic].SexInd == 0 && Register[ic].HIVstage == 0 && Register[ic].VirginInd == 0)
						{
							TotPop[15] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[15] += 1;
							}
						}
						// Chikandiwa2015, M, HIV pos
						if (Register[ic].AgeExact >= 18 && Register[ic].AgeExact < 65 && Register[ic].SexInd == 0 && Register[ic].HIVstage > 0 && Register[ic].VirginInd == 0)
						{
							TotPop[16] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[16] += 1;
							}
						}

						if (Register[ic].AgeExact >= 18 && Register[ic].AgeExact < 65 && Register[ic].SexInd == 1 && Register[ic].HIVstage > 0 &&
							Register[ic].VirginInd == 0 && Register[ic].ARTweeks < 104)
						{
							TotPop[17] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[17] += 1;
							}
						}
						if (Register[ic].AgeExact >= 18 && Register[ic].AgeExact < 65 && Register[ic].SexInd == 1 && Register[ic].HIVstage > 0 &&
							Register[ic].VirginInd == 0 && Register[ic].ARTweeks < 104)
						{
							TotPop[18] += 1;
							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[18] += 1;
							}
						}
						if (Register[ic].AgeExact >= 15 && Register[ic].AgeExact < 50 && Register[ic].SexInd == 1 && Register[ic].VirginInd == 0)
						{
							ig = Register[ic].AgeGroup - 2;
							TotPop[19] += FPCweights[ig];

							if (Register[ic].HPVstage[WhichType] == 1)
							{
								TotHPV[19] += FPCweights[ig];
							}
						}
					}
				}
				// for(ig=0; ig<20; ig++){
				//	Targets[ig].out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHPV[ig] / TotPop[ig];
				// }
			}

			void Pop::GetMacDprev()
			{
				int ic, iy, ig, is, xx, ir;
				int TotPop[14], TotHPV[14];
				for (xx = 0; xx < 14; xx++)
				{
					TotPop[xx] = 0;
					TotHPV[xx] = 0;
				}
				int tpp = Register.size();
				for (ic = 0; ic < tpp; ic++)
				{
					if (Register[ic].AliveInd == 1 && CurrYear == 2000 && Register[ic].AgeExact >= 17 && Register[ic].AgeExact < 50 && Register[ic].SexInd == 1 &&
						Register[ic].VirginInd == 0)
					{
						iy = Register[ic].AgeGroup - 3;
						if (Register[ic].HIVstage == 0)
						{
							TotPop[iy] += 1;
							if (Indiv::AnyHPV(Register[ic].HPVstage, Register[ic].allhpv, {1}))
							{
								TotHPV[iy] += 1;
							}
						}
						else
						{
							TotPop[7 + iy] += 1;
							if (Indiv::AnyHPV(Register[ic].HPVstage, Register[ic].allhpv, {1}))
							{
								TotHPV[7 + iy] += 1;
							}
						}
					}
				}
				for (xx = 0; xx < 14; xx++)
				{
					MacDprev[CurrSim - 1][xx] = 1.0 * TotHPV[xx] / TotPop[xx];
				}
			}
			void Pop::SaveMacDprev(const char *filout)
			{
				int iy, is;
				ofstream file(filout);

				for (iy = 0; iy < 500; iy++)
				{
					for (is = 0; is < 14; is++)
					{
						file << right << MacDprev[iy][is] << "	";
					}
					file << endl;
				}
				file.close();
			}

			void Indiv::GetScreened(
				int ID, double rea, double scr, double ade,
				double tts, double res, double ttC, double CCd,
				double SI, double SII, double SIII, double SIV,
				double SId, double SIId, double SIIId, double SIVd,
				double acceptRand, double efficacyD1Rand,
				double efficacyD2Rand, double D2LossRand)
			{

				int xx, yy, zz;
				// Get age that matches coverage age
				yy = 0;
				if (AgeGroup == 6 || AgeGroup == 7)
				{
					yy = 1;
				}
				else if (AgeGroup == 8 || AgeGroup == 9)
				{
					yy = 2;
				}
				else if (AgeGroup >= 10)
				{
					yy = 3;
				}
				if (HIVstage == 5 || HIVstage == 6)
				{
					zz = 1;
				}
				else
				{
					zz = 0;
				}
				if (InScreen == 1 && timetoCol == 0)
				{
					if (timePassed > timetoScreen && timetoScreen > 0)
					{
						if ((CurrYear == 2020 || CurrYear == 2021) && SIVd < 0.3)
						{	// remember to change rates of
							// entering screening in
							// ScreeningByYear.txt
							if (PerfectSchedule == 0 || CurrYear < ImplementYR)
							{
								if (HIVstage == 5)
								{
									timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
								}
								else
								{
									timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
								}
								if (timetoScreen == 0)
								{
									timetoScreen = 1;
								}
							}
							if (PerfectSchedule == 1 && CurrYear >= ImplementYR)
							{
								timetoScreen = 48;
							}
							repeat = 0;
						}
						if (rea < ScreenReason[zz * 4 + yy][CurrYear - StartYear])
						{
							reason = 0;
						}
						else
						{
							reason = 1;
						}
						// timetoScreen=0;
						if ((InvCaseAlgorithm == 1) && (acceptRand < HPVScaleUpProp[CurrYear - StartYear]))
						{
							HPVScreenAlgorithm_InvCase(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV, SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand); // reusing random numbers to save computational time
						}
						else if (InvCaseAlgorithm == 1)
						{
							ScreenAlgorithm(ID, rea, ade, tts, res, ttC,
											CCd, SI, SII, SIII, SIV, SId,
											SIId, SIIId, SIVd, acceptRand,
											efficacyD1Rand, efficacyD2Rand,
											D2LossRand);
						}
						else if ((HPVDNA == 0 || CurrYear < ImplementYR))
						{
							ScreenAlgorithm(ID, rea, ade, tts, res, ttC,
											CCd, SI, SII, SIII, SIV, SId,
											SIId, SIIId, SIVd, acceptRand,
											efficacyD1Rand, efficacyD2Rand,
											D2LossRand);
						}
						else if (HPVDNA == 1 && CurrYear >= ImplementYR && (MnDCampaign == 0))
						{
							if ((HIVstage >= 5) || (HIVstage < 5 && AgeExact >= 30.0))
							{
								if (HPVDNAThermal == 0)
								{
									HPVScreenAlgorithm(
										ID, rea, ade, tts, res, ttC, CCd, SI,
										SII, SIII, SIV, SId, SIId, SIIId, SIVd);
								}
								else if (HPVDNAThermal == 1)
								{
									if (ThermalORPap < PropThermal)
									{
										HPV_ThermalScreenAlgorithm(
											ID, rea, ade, tts, res, ttC, CCd, SI,
											SII, SIII, SIV, SId, SIId, SIIId,
											SIVd);
									}
									else
									{
										ScreenAlgorithm(
											ID, rea, ade, tts, res, ttC, CCd, SI,
											SII, SIII, SIV, SId, SIId, SIIId,
											SIVd, acceptRand, efficacyD1Rand,
											efficacyD2Rand, D2LossRand);
									}
								}
							}
						}
						else if (HPVDNA == 1 && CurrYear >= ImplementYR && (MnDCampaign == 1))
						{
							if ((HIVstage >= 5) || (HIVstage < 5 && AgeExact >= 30.0))
							{
								MnDScreenAlgorithm(
									ID, rea, ade, tts, res, ttC, CCd, SI, SII,
									SIII, SIV, SId, SIId, SIIId, SIVd,
									acceptRand, efficacyD1Rand,
									efficacyD2Rand, D2LossRand);
							}
							else
							{
								ScreenAlgorithm(ID, rea, ade, tts, res, ttC,
												CCd, SI, SII, SIII, SIV, SId,
												SIId, SIIId, SIVd, acceptRand,
												efficacyD1Rand,
												efficacyD2Rand, D2LossRand);
							}
						}
						timePassed = 0;
					}
					else
					{
						timePassed += 1;
						ScreenCount += 1;
					}
				}
				if (InScreen == 0 &&
					(HIVstage > 0 ||
					 (HIVstage == 0 && AgeGroup > 3)))
				{
					if (scr < ScreenProb[zz * 4 + yy][CurrYear - StartYear] / 48.0)
					{
						if (rea < ScreenReason[zz * 4 + yy][CurrYear - StartYear])
						{
							reason = 0;
						}
						else
						{
							reason = 1;
						}
						// timetoScreen=0;
						if ((InvCaseAlgorithm == 1) && (acceptRand < HPVScaleUpProp[CurrYear - StartYear]))
						{
							HPVScreenAlgorithm_InvCase(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV, SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand);
						}
						else if (InvCaseAlgorithm == 1)
						{
							ScreenAlgorithm(ID, rea, ade, tts, res, ttC,
											CCd, SI, SII, SIII, SIV, SId,
											SIId, SIIId, SIVd, acceptRand,
											efficacyD1Rand, efficacyD2Rand,
											D2LossRand);
						}
						
						else if ((HPVDNA == 0 || CurrYear < ImplementYR))
						{
							ScreenAlgorithm(ID, rea, ade, tts, res, ttC,
											CCd, SI, SII, SIII, SIV, SId,
											SIId, SIIId, SIVd, acceptRand,
											efficacyD1Rand, efficacyD2Rand,
											D2LossRand);
						}
						
						else if ((HPVDNA == 1 && CurrYear >= ImplementYR) && (MnDCampaign == 0))
						{
							if ((HIVstage >= 5) || (HIVstage < 5 && AgeExact >= 30.0))
							{
								if (HPVDNAThermal == 0)
								{
									HPVScreenAlgorithm(
										ID, rea, ade, tts, res, ttC, CCd, SI,
										SII, SIII, SIV, SId, SIId, SIIId, SIVd);
								}
								else if (HPVDNAThermal == 1)
								{
									if (ThermalORPap < PropThermal)
									{
										HPV_ThermalScreenAlgorithm(
											ID, rea, ade, tts, res, ttC, CCd, SI,
											SII, SIII, SIV, SId, SIId, SIIId,
											SIVd);
									}
									else
									{
										ScreenAlgorithm(
											ID, rea, ade, tts, res, ttC, CCd, SI,
											SII, SIII, SIV, SId, SIId, SIIId,
											SIVd, acceptRand, efficacyD1Rand,
											efficacyD2Rand, D2LossRand);
									}
								}
							}
						}
						else if ((HPVDNA == 1 && CurrYear >= ImplementYR) && (MnDCampaign == 1))
						{
							if ((HIVstage >= 5) || (HIVstage < 5 && AgeExact >= 30.0))
							{
								MnDScreenAlgorithm(
									ID, rea, ade, tts, res, ttC, CCd, SI, SII,
									SIII, SIV, SId, SIId, SIIId, SIVd,
									acceptRand, efficacyD1Rand,
									efficacyD2Rand, D2LossRand);
							}
							else
							{
								ScreenAlgorithm(ID, rea, ade, tts, res, ttC,
												CCd, SI, SII, SIII, SIV, SId,
												SIId, SIIId, SIVd, acceptRand,
												efficacyD1Rand, efficacyD2Rand,
												D2LossRand);
							}
						}
						InScreen = 1;
						timePassed = 0;
					}
				}
			}

void Indiv::AdministerTherapeuticVaccine(int ID, double acceptRand, double efficacyD1Rand, double efficacyD2Rand, double D2LossRand)
{
	double adjustedTxVEfficacyD1[13]; 
	double adjustedTxVEfficacyD2[13];
	int zz;
	if (HIVstage == 0 ) {zz=0;}
	else if (HIVstage == 5){zz =1;}
	else {zz=2;}
	if (acceptRand < 0.9) // &&  ((CatchUpVacc == 1 && GotTxV == 0) || (CatchUpVacc == 0)))// - removed because can't track if they are a previous TxV recipient. (means increase wastage if TxV is prophylactic)
	{ 
		RSApop.NewTxV[zz*18 + AgeGroup][CurrYear - StartYear] += 1;
		GotTxV = 1;	
		if (CatchUpVacc == 1){
			RSApop.NewVACC[18 * SexInd + AgeGroup][CurrYear - StartYear] += 1;
			GotVacc = 1;
		}
		for (int xx = 0; xx < 13; xx++)
		{
			if (HPVstage[xx] > 2 && HPVstage[xx] <= 4) // CIN2 and 3
			{
				adjustedTxVEfficacyD1[xx] = TxVEfficacyD1CIN[xx];
				adjustedTxVEfficacyD2[xx] = TxVEfficacyD2CIN[xx];
			}
			else
			{
				adjustedTxVEfficacyD1[xx] = TxVEfficacyD1[xx];
				adjustedTxVEfficacyD2[xx] = TxVEfficacyD2[xx];
			}
			if (!(HIVstage == 0 || HIVstage == 5))
			{
				adjustedTxVEfficacyD1[xx] *= ReductionFactorHIV;
				adjustedTxVEfficacyD2[xx] *= ReductionFactorHIV;
			}
			if (HIVstage == 5)
			{
				adjustedTxVEfficacyD1[xx] *= ReductionFactorHIVART;
				adjustedTxVEfficacyD2[xx] *= ReductionFactorHIVART;
			}
			if ((HPVstage[xx] >= 1 && HPVstage[xx] <= 4) || (HPVstage[xx] == 6))
			{
				if (efficacyD1Rand < adjustedTxVEfficacyD1[xx])
				{
					HPVstageE[xx] = 0;
					if (CatchUpVacc == 1){
						VaccinationStatus[xx] = 1;
					}
				}
			}
		} 
		if (D2LossRand <= 0.7)
		{	RSApop.NewTxV[zz*18 + AgeGroup][CurrYear - StartYear] += 1;
			GotTxV = 2;		
			if (CatchUpVacc == 1 ){RSApop.NewVACC[18 * SexInd + AgeGroup][CurrYear - StartYear] += 1;
				GotVacc =1;
			}
			for (int xx = 0; xx < 13; xx++)
			{if (efficacyD2Rand < adjustedTxVEfficacyD2[xx])
				{HPVstageE[xx] = 0;
					if (CatchUpVacc == 1){
						VaccinationStatus[xx] = 1;
					}
				}
			}
		}
	}
}
	

void Indiv::ScreenAlgorithm(int ID, double rea,  double ade, double tts, double res, double ttC, double CCd, double SI, double SII, double SIII, double SIV, 
							double SId, double SIId, double SIIId, double SIVd, double acceptRand, double efficacyD1Rand, double efficacyD2Rand, double D2LossRand)
{
	int PrevResult, yy, zz;
	PrevResult = ScreenResult;
	//ofstream file("adequate.txt", std::ios::app);
	//file << CurrYear << " " << AgeExact << " " << HIVstage << " " << InScreen << " "<< timetoScreen  << endl;
	timetoScreen=0;
	//file.close();
	int SimCount2 = (CurrSim - 1)/IterationsPerPC;
	yy = 0;
	if (AgeGroup==6||AgeGroup==7) { yy = 1;}
	else if (AgeGroup==8||AgeGroup==9) { yy = 2;}
	else if (AgeGroup==10||AgeGroup==11) { yy = 3;}
	else if (AgeGroup>=12) { yy = 4;}
 
	if (HIVstage==0) {zz=0;}
	//else if(HIVstage==5||HIVstage==6) {zz=1;}
	else if(HIVstage==5) {zz=1;}
	else  {zz=2;}	

	RSApop.NewScreen[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
	/*if((PerfectSchedule==0 ) || CurrYear <ImplementYR){
		//if(HIVstage<5 && (InScreen==0 || ScreenCount>=10*48)){ RSApop.NewScreen[zz*18 + AgeGroup][CurrYear-StartYear] += 1;}
		//if(HIVstage>=5 &&  (InScreen==0 || ScreenCount>=3*48)){ RSApop.NewScreen[zz*18 + AgeGroup][CurrYear-StartYear] += 1;}
		if(HIVstage==0 && (InScreen==0 || ScreenCount>=10*48)){ RSApop.NewScreen[zz*18 + AgeGroup][CurrYear-StartYear] += 1;}
		if(HIVstage!=0 &&  (InScreen==0 || ScreenCount>=3*48)){ RSApop.NewScreen[zz*18 + AgeGroup][CurrYear-StartYear] += 1;}
		//if(HIVstage!=0 &&  (InScreen==0 || ScreenCount>=10*48)){ RSApop.NewScreen[zz*18 + AgeGroup][CurrYear-StartYear] += 1;}
	
	}
	if((PerfectSchedule==1 ) && CurrYear>=ImplementYR){
		if(repeat==0 ) {RSApop.NewScreen[zz*18 + AgeGroup][CurrYear-StartYear] += 1;}
	}*/

	//if(HIVstage==0 && ScreenCount>=10*48){ScreenCount=0;}
	//if(HIVstage!=0 &&  ScreenCount>=3*48){ScreenCount=0;}
	//if(HIVstage!=0 &&  ScreenCount>=10*48){ScreenCount=0;}

  	//ofstream file("GetScreened.txt", std::ios::app);
	//file <<CurrYear << " " << CurrSim << ID << " " <<  timePassed << " " << HIVstage << " " << AgeExact << " " << InScreen << " " << ScreenCount << " " << TrueStage<< endl;
	//file.close();

	//if(HIVstage<5 && ScreenCount>=10*48){ScreenCount=0;}
	//if(HIVstage>=5 &&  ScreenCount>=3*48){ScreenCount=0;}
	//when using the dumb coverage estimate:
	//if(repeat==0 && reason==0 ){ RSApop.NewScreen[zz*18 + AgeGroup][CurrYear-StartYear] += 1;}

	//is pap adequate?
	if(ade < PapAdequacy[CurrYear-StartYear]){
		//if yes, assign Pap result based on estimated sensitivities found in review of SA studies (summing studies together)
		//ofstream file1("adequate.txt", std::ios::app);
		//file1 <<CurrYear << " " << CurrSim << ID << " " <<  timePassed << " " << HIVstage << " " << AgeExact << " " << InScreen << " " << reason << " " << repeat << " " << TrueStage<< endl;
		//file1.close();
		//if(HIVstage>0){
			if(TrueStage==0){
				if(res < 0.954) { ScreenResult = 0;}
				else if(res<(0.954+0.034)){ ScreenResult = 1;}
				else { ScreenResult = 2;}
			}
			if(TrueStage==1){
				if(res < 0.487) { ScreenResult = 0;}
				else if(res < (0.487+0.35)) { ScreenResult = 1;}
				else { ScreenResult = 2;}
			}
			if(TrueStage==2||TrueStage==3){
				if(res < 0.278) { ScreenResult = 0;}
				else if (res < (0.278 + 0.18)) {ScreenResult = 1;}
				else { ScreenResult = 2;}
			}
			if(TrueStage==3 && ScreenResult == 2 && CCd<0.35) {  //CCd is sensitivity of Pap to pick up cancer
				DiagnosedCC = 1;
				//ofstream file("age.txt", std::ios::app);
				//file << CurrYear << " " << HIVstage << " " << AgeExact << endl;
				//file.close();	
				RSApop.NewDiagCancer[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
				if(AnyHPV(HPVstage, hpv1618, cc_un34)) {RSApop.NewDiagCancer1618[AgeGroup][CurrYear-StartYear] += 1;}
				if(HIVstage==5 ) {RSApop.NewDiagCancerART[CurrYear-StartYear] += 1;}
				for(int xx=0; xx<13; xx++) {
					if(HPVstage[xx]==5) {
						HPVstageE[xx]=11;
						RSApop.StageDiag[0][CurrYear-StartYear] += 1;
						RSApop.StageIdiag[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
						//DiagCCPost2000.out[CurrSim-1][0] += 1;
						//StageIdeath = static_cast<int> (48.0 * 126.5 * pow(-log(SI), 1.0/0.61));
						if(SId<0.192){StageIdeath = static_cast<int> (48.0 * 3.08 * pow(-log(SI), 1.0/1.23));}
				        else{StageIrecover = 8 ;}
					}
					if(HPVstage[xx]==8) {
						HPVstageE[xx]=12;
						RSApop.StageDiag[1][CurrYear-StartYear] += 1;
						RSApop.StageIIdiag[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
						//DiagCCPost2000.out[CurrSim-1][0] += 1;
						//StageIIdeath = static_cast<int> (48.0 * 16.28 * pow(-log(SII), 1.0/0.67));
						if(SIId<0.466){StageIIdeath = static_cast<int> (48.0 * 2.39 * pow(-log(SII), 1.0/1.17));}
				        else{StageIIrecover = 24 ;}
					}
					if(HPVstage[xx]==9) {
						HPVstageE[xx]=13;
						RSApop.StageDiag[2][CurrYear-StartYear] += 1;
						RSApop.StageIIIdiag[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
						//DiagCCPost2000.out[CurrSim-1][0] += 1;
						//StageIIIdeath = static_cast<int> (48.0 * 3.91 * pow(-log(SIII), 1.0/0.56));
						if(SIIId<0.715){StageIIIdeath = static_cast<int> (48.0 * 1.18 * pow(-log(SIII), 1.0/0.91));}
       					else{StageIIIrecover = 24 ;}
					}
					if(HPVstage[xx]==10) {
						HPVstageE[xx]=14;
						RSApop.StageDiag[3][CurrYear-StartYear] += 1;
						RSApop.StageIVdiag[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
						//DiagCCPost2000.out[CurrSim-1][0] += 1;
						//StageIVdeath = static_cast<int> (48.0 * 0.53 * pow(-log(SIV), 1.0/0.78));
						if(SIVd<0.931){StageIVdeath = static_cast<int> (48.0 * 0.46 * pow(-log(SIV), 1.0/0.9));}
				        else{StageIVrecover = 24 ;}
					}
				}
			}
			
		//ofstream file2("screenresult.txt", std::ios::app);
		//file2 <<CurrYear << " " << CurrSim << ID << " " <<  HIVstage << " " << AgeExact << " " << ScreenResult << " " <<  TrueStage<< endl;
		//file2.close();
		//Normal screen:
		if(ScreenResult == 0){
			if((WHOScreening==0 && PerfectSchedule==0)||CurrYear<ImplementYR){
				//if(HIVstage==5) {timetoScreen = 5.3 * pow(-log(tts),(1.0/0.78)) * 48;}
				//else if(AgeExact<50) { timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48; }
				if(HIVstage==5) {timetoScreen = 7.9 * pow(-log(tts),(1.0/1.0)) * 48;} //
				else if(AgeExact<50) { timetoScreen = 15.0 * pow(-log(tts),(1.0/1.0)) * 48; }
				else { timetoScreen = 200 * 48; }
				if(timetoScreen==0) {timetoScreen=1;}
			}
			if( PerfectSchedule==1 && CurrYear>=ImplementYR){
				if(HIVstage==5) {timetoScreen = 3 * 48;}
				else if(AgeExact<50) { timetoScreen = 10 * 48; }
				else { timetoScreen = 200 * 48; }
			}
			repeat=0;
		}
		//LSIL screen:
		if(ScreenResult == 1){
			if(PrevResult==1) { 
				//Refer to colposcopy 
				RSApop.GetReferred[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
				if(HIVstage==5) {
					if(ttC<AttendColposcopy[2][CurrYear-StartYear]){
						timetoCol = 24;
						timetoScreen=0;
					}
					else {
						timetoCol=0;
						if( PerfectSchedule==0||CurrYear<ImplementYR){
							//timetoScreen = 5.3 * pow(-log(tts),(1.0/0.78)) * 48;
							timetoScreen = 7.9 * pow(-log(tts),(1.0/1.0)) * 48;
							if(timetoScreen==0) {timetoScreen=1;}
						}	
					}
				}
				if(HIVstage>0 && HIVstage!=5) {
					if(ttC<AttendColposcopy[1][CurrYear-StartYear]){
						timetoCol = 24;
						timetoScreen=0;
					}
					else {
						timetoCol=0;
						if( PerfectSchedule==0||CurrYear<ImplementYR){
							//if(AgeExact<50) { timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48; }
							if(AgeExact<50) { timetoScreen = 15.0 * pow(-log(tts),(1.0/1.0)) * 48; }
							else { timetoScreen = 200 * 48; }
							if(timetoScreen==0) {timetoScreen=1;}
						}	
					}
				}
				if(HIVstage==0) {
					if(ttC<AttendColposcopy[0][CurrYear-StartYear]){
						timetoCol = 24;
						timetoScreen=0;
					}
					else {
						timetoCol=0;
						if( PerfectSchedule==0||CurrYear<ImplementYR){
							//if(AgeExact<50) { timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48; }
							if(AgeExact<50) { timetoScreen = 15.0 * pow(-log(tts),(1.0/1.0)) * 48; }
							else { timetoScreen = 200 * 48; }
							if(timetoScreen==0) {timetoScreen=1;}
						}
					}
				}
				repeat=0;
				//ofstream file("refer.txt", std::ios::app);
				//file << CurrYear << " " <<CurrSim<< ID << " " << HIVstage << " " << TrueStage << " " << AgeExact << " " << timetoCol << endl;
				//file.close();
			}
			else {
				//Supposed to repeat smear in year 
				//Derived from NHLS data 
				if( PerfectSchedule==0||CurrYear<ImplementYR){
					if(HIVstage==5) {timetoScreen = 4.5 * pow(-log(tts),(1.0/0.71)) * 48;}
					else  { timetoScreen = 9.0 * pow(-log(tts),(1.0/0.56)) * 48; }
					if(timetoScreen==0) {timetoScreen=1;}
				}
				if( PerfectSchedule==1 && CurrYear>=ImplementYR){timetoScreen = 48;}
				repeat=1;
			}
		}
		//HSIL/CC screen:
		if(ScreenResult == 2){
			//Refer to colposcopy 
			RSApop.GetReferred[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
				if(HIVstage==5) {
				if(ttC<AttendColposcopy[2][CurrYear-StartYear]){
					timetoCol = 24;
					timetoScreen=0;
				}
				else {
					timetoCol=0;
					if(PerfectSchedule==0||CurrYear<ImplementYR){
						//timetoScreen = 5.3 * pow(-log(tts),(1.0/0.78)) * 48;
						timetoScreen = 7.9 * pow(-log(tts),(1.0/1.0)) * 48;
						if(timetoScreen==0) {timetoScreen=1;}
					}
						
				}
			}
			if(HIVstage>0 && HIVstage!=5) {
				if(ttC<AttendColposcopy[1][CurrYear-StartYear]){
					timetoCol = 24;
					timetoScreen=0;
				}
				else {
					timetoCol=0;
					if( PerfectSchedule==0||CurrYear<ImplementYR){
						//if(AgeExact<50) { timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48; }
						if(AgeExact<50) { timetoScreen = 15.0 * pow(-log(tts),(1.0/1.0)) * 48; }
						else { timetoScreen = 200 * 48; }
						if(timetoScreen==0) {timetoScreen=1;}
					}
						
				}
			}
			if(HIVstage==0) {
				if(ttC<AttendColposcopy[0][CurrYear-StartYear]){
					timetoCol = 24;
					timetoScreen=0;
				}
				else {
					timetoCol=0;
					if( PerfectSchedule==0||CurrYear<ImplementYR){
						//if(AgeExact<50) { timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48; }
						if(AgeExact<50) { timetoScreen = 15.0 * pow(-log(tts),(1.0/1.0)) * 48; }
						else { timetoScreen = 200 * 48; }
						if(timetoScreen==0) {timetoScreen=1;}
					}
				}
			}
			//if(timetoCol>0) {ofstream file3("refer.txt", std::ios::app);
			//file3 << CurrYear << " " <<CurrSim<< ID << " " << HIVstage << " " << TrueStage << " " << AgeExact << " " << timetoCol << endl;
			//file3.close();}
			repeat=0;
		}
	}
	else { 
		//Supposed to repeat smear in 3 months 
		//derived from NHLS data - Weibull distr with scale 31.2  3-months and shape 0.57
		if( PerfectSchedule==0||CurrYear<ImplementYR){
			timetoScreen = 31.2 * pow(-log(tts),(1.0/0.57)) * 12;
			if(timetoScreen==0) {timetoScreen=1;}
		}
		if( PerfectSchedule==1 && CurrYear>=ImplementYR){timetoScreen = 12;}
		
		repeat=1;
	}	
	 if (TxVviaScreeningAlgorithm==1 && CurrYear >= 2030 && AliveInd ==1 ){
		AdministerTherapeuticVaccine( ID,  acceptRand, efficacyD1Rand, efficacyD2Rand, D2LossRand); 
	 }
}
void Indiv::WHOGetScreened(int ID, double rea, double scr, double ade, double tts, double res, double ttC, double CCd, 
									double SI, double SII, double SIII, double SIV, double clr, 
							double SId, double SIId, double SIIId, double SIVd, double acceptRand, double efficacyD1Rand, double efficacyD2Rand, double D2LossRand) 
{
	int iy;	
	iy = CurrYear - StartYear;

	int xx, yy, zz;
	//Get age that matches coverage age
	yy = 0;
	if (AgeGroup==6||AgeGroup==7) { yy = 1;}
	else if (AgeGroup==8||AgeGroup==9) { yy = 2;}
	else if (AgeGroup>=10) { yy = 3;}
		
	if(S7S11==0){
		if(HIVstage==5 || HIVstage==6) {zz=1;}
		else  {zz=0;}
	}
	else {
		if(HIVstage>0) {zz=1;}
		else  {zz=0;}
	}
	//ofstream file("GetScreened.txt", std::ios::app);
	//file <<CurrYear << " " << AgeExact << " " << HIVstage << " " << WHOcoverage[zz*4 + yy][CurrYear-StartYear] << " " << ScreenProb[zz*4 + yy][CurrYear-StartYear] << endl;
	//file.close();
	

	if((S7S11==0 || (S7S11==1 && HIVstage==0)) && scr<WHOcoverage[zz*4 + yy][CurrYear-StartYear]/(10.0*48.0)){
			if(AgeExact>=30.0 && AgeExact<40.0 && Scr35==0){
				Scr35=1;
				WHOScreenAlgorithm(ID, tts, ttC, clr, SI, SII, SIII, SIV, SId, SIId, SIIId, SIVd);
			}
			if(S5S6==1  && AgeExact>=40.0 && AgeExact<50.0 && Scr45==0){
				Scr45=1;
				WHOScreenAlgorithm(ID, tts, ttC, clr, SI, SII, SIII, SIV, SId, SIId, SIIId, SIVd);
			}
	}
	if(S7S11==1 && HIVstage>0 &&  scr<WHOcoverage[zz*4 + yy][CurrYear-StartYear]/(3.0*48.0)) {   
		if(AgeExact>=25.0 && AgeExact<28.0 && Scr25==0){
				Scr25=1;
				WHOScreenAlgorithm(ID, tts, ttC, clr, SI, SII, SIII, SIV, SId, SIId, SIIId, SIVd);
		}
		if(AgeExact>=28.0 && AgeExact<31.0 && Scr28==0){
				Scr28=1;
				WHOScreenAlgorithm(ID, tts, ttC, clr, SI, SII, SIII, SIV, SId, SIId, SIIId, SIVd);
		}
		if(AgeExact>=31.0 && AgeExact<34.0 && Scr31==0){
				Scr31=1;
				WHOScreenAlgorithm(ID, tts, ttC, clr, SI, SII, SIII, SIV, SId, SIId, SIIId, SIVd);
		}
		if(AgeExact>=34.0 && AgeExact<37.0 && Scr34==0){
				Scr34=1;
				WHOScreenAlgorithm(ID, tts, ttC, clr, SI, SII, SIII, SIV, SId, SIId, SIIId, SIVd);
		}
		if(AgeExact>=37.0 && AgeExact<40.0 && Scr37==0){
				Scr37=1;
				WHOScreenAlgorithm(ID, tts, ttC, clr, SI, SII, SIII, SIV, SId, SIId, SIIId, SIVd);
		}
		if(AgeExact>=40.0 && AgeExact<43.0 && Scr40==0){
				Scr40=1;
				WHOScreenAlgorithm(ID, tts, ttC, clr, SI, SII, SIII, SIV, SId, SIId, SIIId, SIVd);
		}
		if(AgeExact>=43.0 && AgeExact<46.0 && Scr43==0){
				Scr43=1;
				WHOScreenAlgorithm(ID, tts, ttC, clr, SI, SII, SIII, SIV, SId, SIId, SIIId, SIVd);
		}
		if(AgeExact>=46.0 && AgeExact<49.0 && Scr46==0){
				Scr46=1;
				WHOScreenAlgorithm(ID, tts, ttC, clr, SI, SII, SIII, SIV, SId, SIId, SIIId, SIVd);
		}
	}

	if(((S3S4==1 && S5S6==0 && S7S11==0 && (AgeExact<30.0||AgeExact>=40.0)) ||
		(S5S6==1 && S7S11==0 && (AgeExact<30.0||AgeExact>=50.0)) ||
		(S7S11==1 && ((HIVstage==0 && (AgeExact<30.0||AgeExact>=50.0))||(HIVstage>0 && (AgeExact<25.0||AgeExact>=49.0))))) &&
		Scr35==0 && Scr45==0 && Scr25==0 && Scr28==0 && Scr31==0 && Scr34==0 && Scr37==0 && Scr40==0 && Scr43==0 && Scr46==0 ){

		if(InScreen==1 && timetoCol==0  &&
			((CurrYear<2045&&AgeExact>=30)||(AgeExact<30)) && timetoScreen<5*48){
			if(timePassed > timetoScreen && timetoScreen>0 ) {
				ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV, SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand, efficacyD2Rand, D2LossRand);
				timePassed=0;
			}
			else {
				timePassed += 1;
				ScreenCount += 1; 
			} 
		}

		if(InScreen==0 && (HIVstage>0||(HIVstage==0 && AgeGroup>3))){
			if (scr < ScreenProb[zz*4 + yy][CurrYear-StartYear]/48.0){
				ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV, SId, SIId, SIIId, SIVd,acceptRand, efficacyD1Rand, efficacyD2Rand, D2LossRand);
				InScreen=1;			
				timePassed=0;			
			}
		}	
	}


}
void Indiv::WHOScreenAlgorithm(int ID,  double tts, double ttC, double clr, double SI, double SII, double SIII, double SIV, 
							double SId, double SIId, double SIIId, double SIVd)
{
	int  xx,  zz;
	
	if (HIVstage==0) {zz=0;}
	else if(HIVstage==5) {zz=1;} //||HIVstage==6
	else  {zz=2;}	

	RSApop.NewHPVScreen[zz*18 + AgeGroup][CurrYear-StartYear] += 1;

	//Assume sensitivity of HPV test as screen is 90% to pick up CIN2 and 94% to pick up CIN3+ (ttS)
	//Assume 10% loss-to-follow-up between screen and treatment. (ttC)
	//Assume 100% treatment efficacy.
	//Cari add: specificity 
	if(AnyHPV(HPVstage, allhpv, {4}) &&
		TrueStage < 3 && tts < 0.94 && ttC < 0.9){
			if(TrueStage==0){RSApop.NewUnnecessary[zz*18 + AgeGroup][CurrYear-StartYear] += 1; }
			RSApop.NewLLETZ[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
					if (clr < 0.15){
						for (xx = 0; xx < 13; xx++)	{
							if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4) {
								HPVstageE[xx] = 1;
								if(HPVstage[xx] == 4){
									WeibullCIN3[xx]=0;
									TimeinCIN3[xx]=0;
								}	
							}
						}
					}
					else{
						for (xx = 0; xx < 13; xx++){
							if (HPVstage[xx] == 1 ||HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4){
								HPVstageE[xx] = 0;
								if(HPVstage[xx] == 4){
									WeibullCIN3[xx]=0;
									TimeinCIN3[xx]=0;
								}	
							}
						}
					}	
	}
	else if(AnyHPV(HPVstage, allhpv, {1,2,3}) &&
		TrueStage < 3 && tts < 0.9 && ttC < 0.9){
			if(TrueStage==0){RSApop.NewUnnecessary[zz*18 + AgeGroup][CurrYear-StartYear] += 1; }
			RSApop.NewLLETZ[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
			if (clr < 0.15){
				for (xx = 0; xx < 13; xx++)	{
					if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4) {
						HPVstageE[xx] = 1;
						if(HPVstage[xx] == 4){
							WeibullCIN3[xx]=0;
							TimeinCIN3[xx]=0;
						}	
					}
				}
			}
			else{
				for (xx = 0; xx < 13; xx++){
					if (HPVstage[xx] == 1 ||HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4){
						HPVstageE[xx] = 0;
						if(HPVstage[xx] == 4){
							WeibullCIN3[xx]=0;
							TimeinCIN3[xx]=0;
						}	
					}
				}
			}
	}
	


	if(TrueStage==3 && tts < 0.94){
		DiagnosedCC=1;
		RSApop.NewDiagCancer[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
		if(AnyHPV(HPVstage, hpv1618, cc_un34)) {RSApop.NewDiagCancer1618[AgeGroup][CurrYear-StartYear] += 1;}
		if(HIVstage==5 ) {RSApop.NewDiagCancerART[CurrYear-StartYear] += 1;}
		for(int xx=0; xx<13; xx++) {
			if(HPVstage[xx]==5) {
				HPVstageE[xx]=11;
				RSApop.StageDiag[0][CurrYear-StartYear] += 1;
				RSApop.StageIdiag[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
				//StageIdeath = static_cast<int> (48.0 * 126.5 * pow(-log(SI), 1.0/0.61));
				if(SId<0.192){StageIdeath = static_cast<int> (48.0 * 3.08 * pow(-log(SI), 1.0/1.23));}
        		else{StageIrecover = 8 ;}

			}
			if(HPVstage[xx]==8) {
				HPVstageE[xx]=12;
				RSApop.StageDiag[1][CurrYear-StartYear] += 1;
				RSApop.StageIIdiag[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
				//StageIIdeath = static_cast<int> (48.0 * 16.28 * pow(-log(SII), 1.0/0.67));
				if(SIId<0.466){StageIIdeath = static_cast<int> (48.0 * 2.39 * pow(-log(SII), 1.0/1.17));}
        		else{StageIIrecover = 24 ;}

			}
			if(HPVstage[xx]==9) {
				HPVstageE[xx]=13;
				RSApop.StageDiag[2][CurrYear-StartYear] += 1;
				RSApop.StageIIIdiag[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
				//StageIIIdeath = static_cast<int> (48.0 * 3.91 * pow(-log(SIII), 1.0/0.56));
				if(SIIId<0.715){StageIIIdeath = static_cast<int> (48.0 * 1.18 * pow(-log(SIII), 1.0/0.91));}
				else{StageIIIrecover = 24 ;}

			}
			if(HPVstage[xx]==10) {
				HPVstageE[xx]=14;
				RSApop.StageDiag[3][CurrYear-StartYear] += 1;
				RSApop.StageIVdiag[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
				//StageIVdeath = static_cast<int> (48.0 * 0.53 * pow(-log(SIV), 1.0/0.78));
				if(SIVd<0.931){StageIVdeath = static_cast<int> (48.0 * 0.46 * pow(-log(SIV), 1.0/0.9));}
		        else{StageIVrecover = 24 ;}

			}
		}
	}
}

void Indiv::GetTreated(int ID, double res, double trt, double clr, double regr, double tts, double SI, double SII, double SIII, double SIV , 
							double SId, double SIId, double SIIId, double SIVd)
{
	int  xx,  zz;
	if (HIVstage==0) {zz=0;}
	else if(HIVstage==5) {zz=1;} //||HIVstage==6
	else  {zz=2;}	

	if (timePassed > timetoCol  ){
		RSApop.NewColposcopy[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
		
		if (TrueStage == 0){
			if(res<0.71){ ColResult=1; }
			else { ColResult = 0; }
		}
		if (TrueStage > 0 && TrueStage < 4){
			if (res < 0.91) { ColResult = 1;}
			else {  ColResult = 0; }
		}
		if(TrueStage==3){
			DiagnosedCC = 1;
			RSApop.NewDiagCancer[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
			if(AnyHPV(HPVstage, hpv1618, cc_un34)) {RSApop.NewDiagCancer1618[AgeGroup][CurrYear-StartYear] += 1;}
			if(HIVstage==5 ) {RSApop.NewDiagCancerART[CurrYear-StartYear] += 1;}
			for (xx = 0; xx < 13; xx++)	{
				if(HPVstage[xx]==5) {
					HPVstageE[xx]=11;
					RSApop.StageDiag[0][CurrYear-StartYear] += 1;
					RSApop.StageIdiag[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
					//DiagCCPost2000.out[CurrSim-1][0] += 1;
					//StageIdeath = static_cast<int> (48.0 * 126.5 * pow(-log(SI), 1.0/0.61));
					if(SId<0.192){StageIdeath = static_cast<int> (48.0 * 3.08 * pow(-log(SI), 1.0/1.23));}
        			else{StageIrecover = 8 ;}
				}
				if(HPVstage[xx]==8) {
					HPVstageE[xx]=12;
					RSApop.StageDiag[1][CurrYear-StartYear] += 1;
					RSApop.StageIIdiag[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
					//DiagCCPost2000.out[CurrSim-1][0] += 1;
					//StageIIdeath = static_cast<int> (48.0 * 16.28 * pow(-log(SII), 1.0/0.67));
					if(SIId<0.466){StageIIdeath = static_cast<int> (48.0 * 2.39 * pow(-log(SII), 1.0/1.17));}
        			else{StageIIrecover = 24 ;}
				}
				if(HPVstage[xx]==9) {
					HPVstageE[xx]=13;
					RSApop.StageDiag[2][CurrYear-StartYear] += 1;
					RSApop.StageIIIdiag[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
							//DiagCCPost2000.out[CurrSim-1][0] += 1;
					//StageIIIdeath = static_cast<int> (48.0 * 3.91 * pow(-log(SIII), 1.0/0.56));
					if(SIIId<0.715){StageIIIdeath = static_cast<int> (48.0 * 1.18 * pow(-log(SIII), 1.0/0.91));}
       				else{StageIIIrecover = 24 ;}
				}
				if(HPVstage[xx]==10) {
					HPVstageE[xx]=14;
					RSApop.StageDiag[3][CurrYear-StartYear] += 1;
					RSApop.StageIVdiag[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
							//DiagCCPost2000.out[CurrSim-1][0] += 1;
					//StageIVdeath = static_cast<int> (48.0 * 0.53 * pow(-log(SIV), 1.0/0.78));
					if(SIVd<0.931){StageIVdeath = static_cast<int> (48.0 * 0.46 * pow(-log(SIV), 1.0/0.9));}
        			else{StageIVrecover = 24 ;}
				}	
			}
		}
	
		if (ColResult == 1 && TrueStage<3){
			//Treat
			GetTreatment.out[AgeGroup][CurrYear-StartYear] += 1;
			if(TrueStage==0){RSApop.NewUnnecessary[zz*18 + AgeGroup][CurrYear-StartYear] += 1; }
			RSApop.NewLLETZ[zz*18 + AgeGroup][CurrYear-StartYear] += 1;
			
			if (HIVstage == 0){
				if (trt < 0.752){
					if (clr < 0.15){
						for (xx = 0; xx < 13; xx++)	{
							if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4) {
								HPVstageE[xx] = 1;
								if(HPVstage[xx] == 4){
									WeibullCIN3[xx]=0;
									TimeinCIN3[xx]=0;
								}	
							}
						}
					}
					else{
						for (xx = 0; xx < 13; xx++){
							if (HPVstage[xx] == 1 ||HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4){
								HPVstageE[xx] = 0;
								if(HPVstage[xx] == 4){
									WeibullCIN3[xx]=0;
									TimeinCIN3[xx]=0;
								}	
							}
						}
					}
				}
				else {
					if(regr<0.5){
						for (xx = 0; xx < 13; xx++)	{
							if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4){
								HPVstageE[xx] = 1;
								if(HPVstage[xx] == 4){
									WeibullCIN3[xx]=0;
									TimeinCIN3[xx]=0;
								}	
							}
						}	
					}
				}
			}
			else{
				if (trt < 0.4){
					if (clr < 0.15){
						for (xx = 0; xx < 13; xx++){
							if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4 ){
								HPVstageE[xx] = 1;
								if(HPVstage[xx] == 4){
									WeibullCIN3[xx]=0;
									TimeinCIN3[xx]=0;
								}
							}
						}
					}
					else{
						for (xx = 0; xx < 13; xx++){
							if (HPVstage[xx] == 1 ||HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4 ){
								HPVstageE[xx] = 0; 
								if(HPVstage[xx] == 4){
									WeibullCIN3[xx]=0;
									TimeinCIN3[xx]=0;
								}	
							}
						}
					}
				}
				else{
					if(regr<0.5){
						for (xx = 0; xx < 13; xx++)	{
							if (HPVstage[xx] == 1 || HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4) {
								HPVstageE[xx] = 1; 
								if(HPVstage[xx] == 4){
									WeibullCIN3[xx]=0;
									TimeinCIN3[xx]=0;
								}	
							}
						}
					}	
				}
			}	
		}
		
		timetoCol = 0;
		timePassed = 0;
		
		if (ColResult == 1 && TrueStage<3){
			if((WHOScreening==0 && PerfectSchedule==0)||CurrYear<ImplementYR){
					if(HIVstage==5) {timetoScreen = 4.5 * pow(-log(tts),(1.0/0.71)) * 48;}
					else  { timetoScreen = 9.0 * pow(-log(tts),(1.0/0.56)) * 48; }
					if(timetoScreen==0) {timetoScreen=1;}
			}
			if( PerfectSchedule==1 && CurrYear>=ImplementYR){timetoScreen = 48;}	
			repeat = 1;
			if(HPVDNAThermal==1){ HPVrepeat=1; }
		}
		if(ColResult==0){
			if((WHOScreening==0 &&  PerfectSchedule==0)||CurrYear<ImplementYR){
				//if (HIVstage == 5){timetoScreen = 5.3 * pow(-log(tts), (1.0 / 0.78)) * 48;}
				//else if(AgeExact<50) { timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48; }
				if (HIVstage == 5){timetoScreen = 7.9 * pow(-log(tts), (1.0 /1.0)) * 48;}
				else if(AgeExact<50) { timetoScreen = 15.0 * pow(-log(tts),(1.0/1.0)) * 48; }
				else { timetoScreen = 200 * 48; }
				
				if (timetoScreen == 0) {timetoScreen = 1;}
			}
			if(PerfectSchedule==1 && CurrYear>=ImplementYR){
				if (HIVstage == 5){timetoScreen = 3 * 48;}
				else if(AgeExact<50) { timetoScreen = 10 * 48; }
				else { timetoScreen = 200 * 48; }
			}
		}
		repeat = 0;
	}
	else { 	timePassed += 1;}
}
void ReadScreenData()
{
	int ia, iy;
	ifstream file;
	stringstream s;
	s << "ScreeningByYear.txt";
	string path = "./input/" + s.str();
	
	file.open(path.c_str());
	if (file.fail()) {
		cerr << "Could not open ScreeningByYear.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	for (ia = 0; ia<8; ia++){
		for (iy = 0; iy<136; iy++){
			file >> ScreenProb[ia][iy]; 
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<8; ia++){
		for (iy = 0; iy<136; iy++){
			file >> ScreenReason[ia][iy]; 
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<3; ia++){
		for (iy = 0; iy<136; iy++){
			file >> AttendColposcopy[ia][iy]; 
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<136; iy++){
			file >> PapAdequacy[iy]; 
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<136; iy++){
			file >> PropVaccinated[iy]; 
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<8; ia++){
		for (iy = 0; iy<136; iy++){
			file >> WHOcoverage[ia][iy]; 
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<136; iy++){
			file >> CampaignYear[iy]; 
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<136; iy++){
			file >> TAScaleUpRate[iy]; 
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<136; iy++){
			file >> HPVScaleUpProp[iy]; 
	}
}

void Pop::SaveNewScreen(const char* filout)
{
	int iy, is;
		ostringstream s;

		if (process_num >0){
			s << process_num << "_"  << filout;
		}
		else{
			s <<  filout;
		}
				
		string path = "./output/" + s.str();
		ofstream file(path.c_str()); // Converts s to a C string

	for (iy = 0; iy < 54; iy++){
		for (is = 0; is<136; is++){
			file << right << RSApop.ModelCoverage[iy][is] << "	"; //Pap smears 1, 2, 3
		}
		file << endl;
	}
	for (iy = 0; iy < 54; iy++){
		for (is = 0; is<136; is++){
			file << right << RSApop.ModelHPVCoverage[iy][is] << "	"; //HPV-DNA tests 4, 5, 6
		}
		file << endl;
	}
	for (iy = 0; iy < 54; iy++){
		for (is = 0; is<136; is++){
			file << right << RSApop.ModelColpCoverage[iy][is] << "	"; //Colposcopies 7, 8, 9
		}
		file << endl;
	}
	for (iy = 0; iy < 54; iy++){
		for (is = 0; is<136; is++){
			file << right << RSApop.ModelLLETZCoverage[iy][is] << "	"; //LLETZ performed 10, 11, 12
		}
		file << endl;
	}
	for (iy = 0; iy < 54; iy++){
		for (is = 0; is<136; is++){
			file << right << RSApop.ModelVATCoverage[iy][is] << "	"; //VAT 13, 14, 15 
		}
		file << endl;
	}
	for (iy = 0; iy < 54; iy++){
		for (is = 0; is<136; is++){
			file << right << RSApop.ModelThermalCoverage[iy][is] << "	"; //Thermal ablations 16, 17, 18
		}
		file << endl;
	}
	for (iy = 0; iy < 54; iy++){
		for (is = 0; is<136; is++){
			file << right << RSApop.ModelUnnecessaryCoverage[iy][is] << "	"; //Unnecessary treatments 19, 20, 21
		}
		file << endl;
	}
	for (iy = 0; iy < 36; iy++){
		for (is = 0; is<136; is++){
			file << right << RSApop.ModelVaccCoverage[iy][is] << "	"; //Catch-up vaccs for people on/initiating ART 22, 23
		}
		file << endl;
	}
	for (iy = 0; iy < 54; iy++){
		for (is = 0; is<136; is++){
			file << right << RSApop.ModelGetReferred[iy][is] << "	"; //were referred to treatment 24, 25, 26
		}
		file << endl;
	}
	for (iy = 0; iy < 54; iy++){
		for (is = 0; is<136; is++){
			file << right << RSApop.ModelTxVCoverage[iy][is] << "	"; // Therapeutic vaccines given 27, 28, 29
		}
		file << endl;
	}
	for (iy = 0; iy < 54; iy++){
		for (is = 0; is<136; is++){
			file << right << RSApop.ModelReflexCoverage[iy][is] << "	"; // Reflex cytology given 30, 31, 32
		}
		file << endl;
	}
	file.close();
}



void ReadLifeTimePartners()
{
	int ia, iy;
	ifstream file;
	stringstream s;
	s << "LifetimePartnersInit.txt";
	string path = "./input/" + s.str();
	
	file.open(path.c_str());
	for (ia = 0; ia<320; ia++){
		for (iy = 0; iy<2; iy++){
			file >> LTP[ia][iy];
		}
	}
	file.close();
}

void Pop::CalcModelCoverage()
{
	
	double TotPop[15], TotVacc[54];
	int ic, ig, xx, yy, zz;
	
	/*for(ig = 0; ig<15; ig++){
		TotPop[ig] = 0; 
	}

	int tpp = Register.size();

	for (ic = 0; ic<tpp; ic++){
		if(Register[ic].SexInd==1 && Register[ic].AliveInd == 1 && Register[ic].AgeGroup>=3){
			if (Register[ic].AgeGroup==3||Register[ic].AgeGroup==4||Register[ic].AgeGroup==5) { yy = 0;}
			//if (Register[ic].AgeGroup==4|| Register[ic].AgeGroup==5) { yy = 0;}
			else if (Register[ic].AgeGroup==6||Register[ic].AgeGroup==7) { yy = 1;}
			else if (Register[ic].AgeGroup==8||Register[ic].AgeGroup==9) { yy = 2;}
			else if (Register[ic].AgeGroup==10||Register[ic].AgeGroup==11) { yy = 3;}
			else if (Register[ic].AgeGroup>=12) { yy = 4;}
			
			if (Register[ic].HIVstage==0) {zz=0;}
			else if(Register[ic].HIVstage==5||Register[ic].HIVstage==6) {zz=1;}
			else  {zz=2;}

			TotPop[zz*5 + yy] += 1;
		}	
	}*/
	/*for(ig = 0; ig<54; ig++){
		TotVacc[ig] = 0; 
	}

	int tpp = Register.size();

	for (ic = 0; ic<tpp; ic++){
		if(Register[ic].SexInd==1 && Register[ic].AliveInd == 1  && CurrYear>=ImplementYR && 
		((WHOscenario==0 && Register[ic].AgeExact>=9)||(WHOscenario==1 && Register[ic].AgeExact>=10))){	
			if (Register[ic].HIVstage==0) {zz=0;}
			else if(Register[ic].HIVstage==5) {zz=1;}
			else  {zz=2;}

			if(Register[ic].GotVacc==1)	{TotVacc[zz*18 + Register[ic].AgeGroup] += 1;}
			//if(CurrYear==ImplementYR && Register[ic].GotVacc==1){
			//	cout <<ic << " " << Register[ic].AgeExact << " " << Register[ic].AgeGroup << " " << Register[ic].DOB << endl;}
		}	
	}*/
	for(ig = 0; ig<54; ig++){
		RSApop.ModelHPVCoverage[ig][CurrYear-StartYear] += RSApop.NewHPVScreen[ig][CurrYear-StartYear];
		RSApop.ModelCoverage[ig][CurrYear-StartYear] += RSApop.NewScreen[ig][CurrYear-StartYear];
		RSApop.ModelReflexCoverage[ig][CurrYear-StartYear] += RSApop.NewReflex[ig][CurrYear-StartYear];
		RSApop.ModelTxVCoverage[ig][CurrYear-StartYear] += RSApop.NewTxV[ig][CurrYear-StartYear];
		RSApop.ModelColpCoverage[ig][CurrYear-StartYear] += RSApop.NewColposcopy[ig][CurrYear-StartYear];
		RSApop.ModelLLETZCoverage[ig][CurrYear-StartYear] += RSApop.NewLLETZ[ig][CurrYear-StartYear];
		RSApop.ModelUnnecessaryCoverage[ig][CurrYear-StartYear] += RSApop.NewUnnecessary[ig][CurrYear-StartYear];
		RSApop.ModelVATCoverage[ig][CurrYear-StartYear] += RSApop.NewVAT[ig][CurrYear-StartYear];
		RSApop.ModelThermalCoverage[ig][CurrYear-StartYear] += RSApop.NewThermal[ig][CurrYear-StartYear];
		RSApop.ModelGetReferred[ig][CurrYear-StartYear] += RSApop.GetReferred[ig][CurrYear-StartYear];
	}
	for(ig = 0; ig<36; ig++){
		RSApop.ModelVaccCoverage[ig][CurrYear-StartYear] += RSApop.NewVACC[ig][CurrYear-StartYear];
	}
	
}	

void Pop::GetCurrHPVstage()
{
    int ic, xx;
	int tpp = Register.size();

    for (ic = 0; ic<tpp; ic++){
        if(Register[ic].AliveInd == 1 &&  Register[ic].VirginInd == 0){
			if(Register[ic].SexInd == 1){
				if(Register[ic].HIVstage==0){
					if(Register[ic].AgeGroup>=6 && Register[ic].AgeGroup<13){
						HPVtransitionCC.TotPop30to65neg[CurrYear-StartYear] += 1;
						if (Register[ic].TrueStage > 0 ){HPVtransitionCC.TotAB30to65neg[CurrYear-StartYear] += 1; }
						if (Register[ic].TrueStage > 1 ){HPVtransitionCC.TotHSIL30to65neg[CurrYear-StartYear] += 1; }   
					}
				}
				if(Register[ic].HIVstage==1||Register[ic].HIVstage==2||Register[ic].HIVstage==3||Register[ic].HIVstage==4||Register[ic].HIVstage==6){
					if (Register[ic].AgeExact>=18 && Register[ic].AgeExact<60 ){
						HPVtransitionCC.TotPop18to60noart[CurrYear-StartYear] += 1;
						if (Register[ic].TrueStage > 0 ){HPVtransitionCC.TotAB18to60noart[CurrYear-StartYear] += 1; }
						if (Register[ic].TrueStage > 1 ){HPVtransitionCC.TotHSIL18to60noart[CurrYear-StartYear] += 1; }   
					}
				}
				if(Register[ic].HIVstage==5){   
					if (Register[ic].AgeExact>=18 && Register[ic].AgeExact<66 ){
						HPVtransitionCC.TotPop18to65art[CurrYear-StartYear] += 1;
						if (Register[ic].TrueStage > 0 ){HPVtransitionCC.TotAB18to65art[CurrYear-StartYear] += 1; }
						if (Register[ic].TrueStage > 1 ){HPVtransitionCC.TotHSIL18to65art[CurrYear-StartYear] += 1; }   
					}
				}
				if(Register[ic].HIVstage>0 && Register[ic].HIVstage!=5){   
					if (Register[ic].AgeGroup>=6 && Register[ic].AgeGroup<13){
						HPVtransitionCC.TotPop30to65pos[CurrYear-StartYear] += 1;
						if (Register[ic].TrueStage > 0 ){HPVtransitionCC.TotAB30to65pos[CurrYear-StartYear] += 1; }
						if (Register[ic].TrueStage > 1 ){HPVtransitionCC.TotHSIL30to65pos[CurrYear-StartYear] += 1; }   
					}
				}
				if(Register[ic].HIVstage==5){   
					if (Register[ic].AgeGroup>=6 && Register[ic].AgeGroup<13){
						HPVtransitionCC.TotPop30to65art[CurrYear-StartYear] += 1;
						if (Register[ic].TrueStage > 0 ){HPVtransitionCC.TotAB30to65art[CurrYear-StartYear] += 1; }
						if (Register[ic].TrueStage > 1 ){HPVtransitionCC.TotHSIL30to65art[CurrYear-StartYear] += 1; }   
					}
				}
				if(Register[ic].HIVstage==0 && Register[ic].AgeExact>=15 && Register[ic].AgeExact<65){
					HPVtransitionCC.TotPop15to65negF[CurrYear-StartYear] += 1;
					if(Indiv::AnyHPV(Register[ic].HPVstage,  Register[ic].allhpv, {1})||Register[ic].TrueStage>0){
						HPVtransitionCC.TotHPV15to65negF[CurrYear-StartYear] += 1;
					}
					if (Register[ic].TrueStage > 1 ){HPVtransitionCC.TotHSILNEG[CurrYear-StartYear] += 1; }   
				}
				if(Register[ic].HIVstage>0 && Register[ic].AgeExact>=15 && Register[ic].AgeExact<65){
					HPVtransitionCC.TotPop15to65posF[CurrYear-StartYear] += 1;
					if(Indiv::AnyHPV(Register[ic].HPVstage,  Register[ic].allhpv, {1})||Register[ic].TrueStage>0){
						HPVtransitionCC.TotHPV15to65posF[CurrYear-StartYear] += 1;
					}
					if (Register[ic].TrueStage > 1 ){HPVtransitionCC.TotHSILPOS[CurrYear-StartYear] += 1; }   	
				}
				if(Register[ic].HIVstage==5 && Register[ic].AgeExact>=15 && Register[ic].AgeExact<65){
					HPVtransitionCC.TotPop15to65artF[CurrYear-StartYear] += 1;
					if(Indiv::AnyHPV(Register[ic].HPVstage,  Register[ic].allhpv, {1})||Register[ic].TrueStage>0){
						HPVtransitionCC.TotHPV15to65artF[CurrYear-StartYear] += 1;
					}
					if (Register[ic].TrueStage > 1 ){HPVtransitionCC.TotHSILART[CurrYear-StartYear] += 1; }   	
				}
				if(Register[ic].HIVstage==0 && Register[ic].AgeExact>=15){
					HPVtransitionCC.TotPop15upnegF[CurrYear-StartYear] += 1;
					if(Register[ic].TrueStage>2){
						HPVtransitionCC.TotCCNEG[CurrYear-StartYear] += 1;
					}
				}
				if(Register[ic].HIVstage>0 && Register[ic].AgeExact>=15){
					HPVtransitionCC.TotPop15upposF[CurrYear-StartYear] += 1;
					if(Register[ic].TrueStage>2){
						HPVtransitionCC.TotCCPOS[CurrYear-StartYear] += 1;
					} 	
				}
				if(Register[ic].HIVstage==5 && Register[ic].AgeExact>=15){
					HPVtransitionCC.TotPop15upartF[CurrYear-StartYear] += 1;
					if(Register[ic].TrueStage>2){
						HPVtransitionCC.TotCCART[CurrYear-StartYear] += 1;
					}
				}
			}
			if(Register[ic].SexInd == 0){
				if(Register[ic].HIVstage==0 && Register[ic].AgeExact>=15 && Register[ic].AgeExact<65){
					HPVtransitionCC.TotPop15to65negM[CurrYear-StartYear] += 1;
					if(Indiv::AnyHPV(Register[ic].HPVstage,  Register[ic].allhpv, {1})){
						HPVtransitionCC.TotHPV15to65negM[CurrYear-StartYear] += 1;
					}	
				}
				if(Register[ic].HIVstage>0 && Register[ic].AgeExact>=15 && Register[ic].AgeExact<65){
					HPVtransitionCC.TotPop15to65posM[CurrYear-StartYear] += 1;
					if(Indiv::AnyHPV(Register[ic].HPVstage,  Register[ic].allhpv, {1})){
						HPVtransitionCC.TotHPV15to65posM[CurrYear-StartYear] += 1;
					}					
				}
				if(Register[ic].HIVstage==5 && Register[ic].AgeExact>=15 && Register[ic].AgeExact<65){
					HPVtransitionCC.TotPop15to65artM[CurrYear-StartYear] += 1;
					if(Indiv::AnyHPV(Register[ic].HPVstage,  Register[ic].allhpv, {1})){
						HPVtransitionCC.TotHPV15to65artM[CurrYear-StartYear] += 1;
					}					
				}
			}
		}
    }
}
void Pop::SaveCancerCases(const char* filout)
{
	int ia, is;
	ostringstream s;

		if (process_num >0){
			s << process_num << "_"  << filout;
		}
		else{
			s <<  filout;
		}
				
		string path = "./output/" + s.str();
		ofstream file(path.c_str()); // Converts s to a C string

	for (ia = 0; ia < 54; ia++){
		for (is = 0; is<136; is++){
			file << right << NewDiagCancer[ia][is] << "	";
		}
		file << endl;
	}
	for (ia = 0; ia < 18; ia++){
		for (is = 0; is<136; is++){
			file << right << NewCancer[ia][is] << "	";
		}
		file << endl;
	}

	file.close();
	
}

void Indiv::MnDScreenAlgorithm(int ID, double rea, double ade, double tts,
                               double res, double ttC, double CCd, double SI,
                               double SII, double SIII, double SIV, double SId,
                               double SIId, double SIIId, double SIVd,
                               double acceptRand, double efficacyD1Rand,
                               double efficacyD2Rand, double D2LossRand) {
  int xx, yy, zz;
  int SimCount2 = (CurrSim - 1) / IterationsPerPC;
  if (HIVstage == 0) {
    zz = 0;
  } else if (HIVstage == 5) {
    zz = 1;
  } //||HIVstage==6
  else {
    zz = 2;
  }
  // before handling all cytology/treatment logic give TxV dose
	 if (TxVviaScreeningAlgorithm == 1 && CurrYear >= 2030) {
          AdministerTherapeuticVaccine(ID, acceptRand, efficacyD1Rand,
                                       efficacyD2Rand, D2LossRand);
      }
  RSApop.NewHPVScreen[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
  // assign time to screen based on HIV status and age
  if (ade < 0.95) {
    if (HPVstatus == 0) { // rescreen after 10 years or 3 years if HIV +
      if ((WHOScreening == 0 && PerfectSchedule == 0) ||
          CurrYear < ImplementYR) {
        // if(HIVstage==5) {timetoScreen = 5.3 * pow(-log(tts),(1.0/0.78)) *
        // 48;} else if(AgeExact<50) { timetoScreen = 15.0 *
        // pow(-log(tts),(1.0/0.83)) * 48; }
        if (HIVstage == 5) {
          timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
        } else if (AgeExact < 50) {
          timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
        } else {
          timetoScreen = 200 * 48;
        }
        if (timetoScreen == 0) {
          timetoScreen = 1;
        }
      }
      if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
        if (HIVstage == 5) {
          timetoScreen = 3 * 48;
        } else if (AgeExact < 50) {
          timetoScreen = 10 * 48;
        } else {
          timetoScreen = 200 * 48;
        }
      }
      repeat = 0;
      HPVrepeat = 0;
    } else {
		if (TxVtoHPVPos == 1 && HPVstatus == 1 && (CurrYear >= 2030)){
			AdministerTherapeuticVaccine(ID, acceptRand, efficacyD1Rand, efficacyD2Rand, D2LossRand);
		}
		 RSApop.NewReflex[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
		// HPVstatus == 1
        // reflex cytology means no need for acceptance rate drop off
      // assuming no genotyping, as per MnD,  screen the individual. 2 is
      // HSIL, 1 is LSIL, 0 is normal
      if (HIVstage == 0) {
        if (TrueStage == 0) {
          if (res < 0.732) {
            ScreenResult = 0;
          } else if (res < 0.732 + 0.146) {
            ScreenResult = 1;
          } else {
            ScreenResult = 2;
          }
        }
        if (TrueStage == 1) {
          if (res < 0.586) {
            ScreenResult = 0;
          } else if (res < 0.586 + 0.261) {
            ScreenResult = 1;
          } else {
            ScreenResult = 2;
          }
        }
        if (TrueStage >= 2 && TrueStage < 5) {
          if (res < 0.453) {
            ScreenResult = 0;
          } else if (res < 0.453 + 0.203) {
            ScreenResult = 1;
          } else {
            ScreenResult = 2;
          }
        }
      }
      if (HIVstage > 0) {
        if (TrueStage == 0) {
          if (res < 0.62) {
            ScreenResult = 0;
          } else if (res < 0.62 + 0.16) {
            ScreenResult = 1;
          } else {
            ScreenResult = 2;
          }
        }
        if (TrueStage == 1) {
          if (res < 0.545) {
            ScreenResult = 0;
          } else if (res < 0.545 + 0.15) {
            ScreenResult = 1;
          } else {
            ScreenResult = 2;
          }
        }
        if (TrueStage >= 2 && TrueStage < 5) {
          if (res < 0.209) {
            ScreenResult = 0;
          } else if (res < 0.209 + 0.157) {
            ScreenResult = 1;
          } else {
            ScreenResult = 2;
          }
        }
      }
      // track cancer diagnoses, stage at diagnosis for screen-detected cases,
      // and survival outcomes
      if (TrueStage == 3 && ScreenResult == 2 &&
          CCd < 0.35) { // CCd is sensitivity of Pap to pick up cancer
        DiagnosedCC = 1;
        RSApop.NewDiagCancer[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
        if (HIVstage == 5) {
          RSApop.NewDiagCancerART[CurrYear - StartYear] += 1;
        }
        for (int xx = 0; xx < 13; xx++) {
          if (HPVstage[xx] == 5) {
            HPVstageE[xx] = 11;
            RSApop.StageDiag[0][CurrYear - StartYear] += 1;
            RSApop.StageIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
            if (SId < 0.192) {
              StageIdeath =
                  static_cast<int>(48.0 * 3.08 * pow(-log(SI), 1.0 / 1.23));
            } else {
              StageIrecover = 8;
            }
          }
          if (HPVstage[xx] == 8) {
            HPVstageE[xx] = 12;
            RSApop.StageDiag[1][CurrYear - StartYear] += 1;
            RSApop.StageIIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
            if (SIId < 0.466) {
              StageIIdeath =
                  static_cast<int>(48.0 * 2.39 * pow(-log(SII), 1.0 / 1.17));
            } else {
              StageIIrecover = 24;
            }
          }
          if (HPVstage[xx] == 9) {
            HPVstageE[xx] = 13;
            RSApop.StageDiag[2][CurrYear - StartYear] += 1;
            RSApop.StageIIIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
            if (SIIId < 0.715) {
              StageIIIdeath =
                  static_cast<int>(48.0 * 1.18 * pow(-log(SIII), 1.0 / 0.91));
            } else {
              StageIIIrecover = 24;
            }
          }
          if (HPVstage[xx] == 10) {
            HPVstageE[xx] = 14;
            RSApop.StageDiag[3][CurrYear - StartYear] += 1;
            RSApop.StageIVdiag[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
            if (SIVd < 0.931) {
              StageIVdeath =
                  static_cast<int>(48.0 * 0.46 * pow(-log(SIV), 1.0 / 0.9));
            } else {
              StageIVrecover = 24;
            }
          }
        }
      }
      // Normal screen: rescreen in 1 year
      if (ScreenResult == 0) {
        // Supposed to repeat smear in year
        // Derived from NHLS data
        if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
          if (HIVstage == 5) {
            timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
          } else {
            timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
          }
          if (timetoScreen == 0) {
            timetoScreen = 1;
          }
        }
        if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
          timetoScreen = 48;
        }
        repeat = 1;
      }
      // LSIL screen: VIA/Assess suitability for treatment
      //  if suitable, TA, if not suitable or if TA not avaliable, refer to
      //  colposcopy
      else if (ScreenResult == 1) {
        if (CCd < 0.913) { // 91.3% suitable for ablation after VAT
          GetTreatment.out[AgeGroup][CurrYear - StartYear] += 1;
          if (TrueStage == 0) {
            RSApop.NewUnnecessary[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                1; 
          }
          RSApop.NewThermal[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
          if (HIVstage == 0) {
            if (res < 0.652) { // 0.652
              for (xx = 0; xx < 13; xx++) {
                if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                    HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                  HPVstageE[xx] = 0;
                  if (HPVstage[xx] == 4) {
                    WeibullCIN3[xx] = 0;
                    TimeinCIN3[xx] = 0;
                  }
                }
              }
            } else {
              for (xx = 0; xx < 13; xx++) {
                if (HPVstage[xx] == 2) {
                  HPVstageE[xx] = 1;
                }
                if (HPVstage[xx] == 3) {
                  HPVstageE[xx] = 2;
                }
                if (HPVstage[xx] == 4) {
                  HPVstageE[xx] = 3;
                  WeibullCIN3[xx] = 0;
                  TimeinCIN3[xx] = 0;
                }
              }
            }
          } else {
            if (res < 0.385) { //0.385
              for (xx = 0; xx < 13; xx++) {
                if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                    HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                  HPVstageE[xx] = 0;
                  if (HPVstage[xx] == 4) {
                    WeibullCIN3[xx] = 0;
                    TimeinCIN3[xx] = 0;
                  }
                }
              }
            } else {
              for (xx = 0; xx < 13; xx++) {
                if (HPVstage[xx] == 2) {
                  HPVstageE[xx] = 1;
                }
                if (HPVstage[xx] == 3) {
                  HPVstageE[xx] = 2;
                }
                if (HPVstage[xx] == 4) {
                  HPVstageE[xx] = 3;
                  WeibullCIN3[xx] = 0;
                  TimeinCIN3[xx] = 0;
                }
              }
            }
          }
          // Supposed to repeat smear in year
          if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
            if (HIVstage == 5) {
              timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
            } else {
              timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
            }
            if (timetoScreen == 0) {
              timetoScreen = 1;
            }
          }
          if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
            timetoScreen = 48;
          }
          repeat = 1;
        } else { // if not suitable for ablation, refer to colposcopy 
                 // Refer to colposcopy
          RSApop.GetReferred[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
          if (HIVstage == 5) {
            if (ttC < AttendColposcopy[2][CurrYear - StartYear]) {
              timetoCol = 24;
              timetoScreen = 0;
            } else {
              timetoCol = 0;
              if (PerfectSchedule == 0 ||
                  CurrYear < ImplementYR) { // timetoScreen = 5.3 *
                                            // pow(-log(tts),(1.0/0.78)) * 48;
                timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
                if (timetoScreen == 0) {
                  timetoScreen = 1;
                }
              }
            }
          }
          if (HIVstage > 0 && HIVstage != 5) {
            if (ttC < AttendColposcopy[1][CurrYear - StartYear]) {
              timetoCol = 24;
              timetoScreen = 0;
            } else {
              timetoCol = 0;
              if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                // if(AgeExact<50) { timetoScreen = 15.0 *
                // pow(-log(tts),(1.0/0.83)) * 48; }
                if (AgeExact < 50) {
                  timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                } else {
                  timetoScreen = 200 * 48;
                }
                if (timetoScreen == 0) {
                  timetoScreen = 1;
                }
              }
            }
          }
          if (HIVstage == 0) {
            if (ttC < AttendColposcopy[0][CurrYear - StartYear]) {
              timetoCol = 24;
              timetoScreen = 0;
            } else {
              timetoCol = 0;
              if (PerfectSchedule == 0 ||
                  CurrYear <
                      ImplementYR) { // if(AgeExact<50) {timetoScreen = 15.0 *
                                     // pow(-log(tts),(1.0/0.83)) * 48; }
                if (AgeExact < 50) {
                  timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                } else {
                  timetoScreen = 200 * 48;
                }
                if (timetoScreen == 0) {
                  timetoScreen = 1;
                }
              }
            }
          }
          repeat = 0;
          // ofstream file("refer.txt", std::ios::app);
          // file << CurrYear << " " <<CurrSim<< ID << " " << HIVstage << " "
          // << TrueStage << " " << AgeExact << " " << timetoCol << endl;
          // file.close();
          // Supposed to repeat smear in year
          // Derived from NHLS data
          if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
            if (HIVstage == 5) {
              timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
            } else {
              timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
            }
            if (timetoScreen == 0) {
              timetoScreen = 1;
            }
          }
          if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
            timetoScreen = 48;
          }
          repeat = 1;
        }
      }
      // HSIL/CC screen: Colposcopy and treatment (LLETZ/TA)
      else if (ScreenResult == 2 && TrueStage < 3) {
        RSApop.GetReferred[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
        if (HIVstage == 5) {
          if (ttC < AttendColposcopy[2][CurrYear - StartYear]) {
            timetoCol = 24;
            timetoScreen = 0;
          } else {
            timetoCol = 0;
            if (PerfectSchedule == 0 ||
                CurrYear < ImplementYR) { // timetoScreen = 5.3 *
                                          // pow(-log(tts),(1.0/0.78)) * 48;
              timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
              if (timetoScreen == 0) {
                timetoScreen = 1;
              }
            }
          }
        }
        if (HIVstage > 0 && HIVstage != 5) {
          if (ttC < AttendColposcopy[1][CurrYear - StartYear]) {
            timetoCol = 24;
            timetoScreen = 0;
          } else {
            timetoCol = 0;
            if (PerfectSchedule == 0 ||
                CurrYear <
                    ImplementYR) { // if(AgeExact<50) { timetoScreen = 15.0 *
                                   // pow(-log(tts),(1.0/0.83)) * 48; }
              if (AgeExact < 50) {
                timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
              } else {
                timetoScreen = 200 * 48;
              }
              if (timetoScreen == 0) {
                timetoScreen = 1;
              }
            }
          }
        }
        if (HIVstage == 0) {
          if (ttC < AttendColposcopy[0][CurrYear - StartYear]) {
            timetoCol = 24;
            timetoScreen = 0;
          } else {
            timetoCol = 0;
            if (PerfectSchedule == 0 ||
                CurrYear <
                    ImplementYR) { // if(AgeExact<50) { timetoScreen = 15.0 *
                                   // pow(-log(tts),(1.0/0.83)) * 48; }
              if (AgeExact < 50) {
                timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
              } else {
                timetoScreen = 200 * 48;
              }
              if (timetoScreen == 0) {
                timetoScreen = 1;
              }
            }
          }
        }
        // schedule post-treatment follow-up smear in 1 year
        // Supposed to repeat smear in year
        // Derived from NHLS data
        if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
          if (HIVstage == 5) {
            timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
          } else {
            timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
          }
          if (timetoScreen == 0) {
            timetoScreen = 1;
          }
        }
        if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
          timetoScreen = 48;
        }
        repeat = 1;
      }
     
    }
	 
  } else { // Supposed to repeat smear in 3 months
    // derived from NHLS data - Weibull distr with scale 31.2  3-months
    // and shape 0.57
    if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
      timetoScreen = 31.2 * pow(-log(tts), (1.0 / 0.57)) * 12;
      if (timetoScreen == 0) {
        timetoScreen = 1;
      }
    }
    if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
      timetoScreen = 12;
    }
    repeat = 1;
  }
}
void Indiv::HPVScreenAlgorithm_InvCase(
	int ID, double rea, double ade, double tts, double res, double ttC,
	double CCd, double SI, double SII, double SIII, double SIV, double SId,
	double SIId, double SIIId, double SIVd,
	double viaVis, double taRand)
{

	int xx, yy, zz;
	int SimCount2 = (CurrSim - 1) / IterationsPerPC;

	if (HIVstage == 0)
	{
		zz = 0;
	}
	else if (HIVstage == 5)
	{
		zz = 1;
	}
	else
	{
		zz = 2;
	}

	RSApop.NewHPVScreen[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;

	if (ade < 0.95)
	{
		if (HPVstatus == 0)
		{ // HPV negative
			if ((WHOScreening == 0 && PerfectSchedule == 0) || CurrYear < ImplementYR)
			{
				if (HIVstage == 5)
				{
					timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
				}
				else if (AgeExact < 50)
				{
					timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
				}
				else
				{
					timetoScreen = 200 * 48;
				}
				if (timetoScreen == 0)
				{
					timetoScreen = 1;
				}
			}
			if (PerfectSchedule == 1 && CurrYear >= ImplementYR)
			{
				if (HIVstage == 5)
				{
					timetoScreen = 5 * 48;
				}
				else if (AgeExact < 55)
				{
					timetoScreen = 10 * 48;
				}
				else
				{
					timetoScreen = 200 * 48;
				}
			}
			repeat = 0;
			HPVrepeat = 0;
		}

		else if (rea < 0.9)
		{ // 90% come back for results

			if (HPVGenotyping == 1)
			{

				if (HPVrepeat == 0)
				{

					// HPV 16/18/45 positive (not cancer stages 3/4)
					if (AnyHPV(HPVstage, hpv161845, {1, 2, 3, 4}) && TrueStage < 4)
					{

						// Thermal Ablation scale-up
						double TA_coverage = 0.0;
						int iy = CurrYear - StartYear;
						TA_coverage = TAScaleUpRate[iy];
						if (taRand < TA_coverage)
						{
							// Full pathway: VIA -- Thermal Ablation
							RSApop.NewVAT[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;

							if (viaVis < 0.95){ // 0.95 assummed based suitability for TA ** double check **
							if((TrueStage ==0 || TrueStage == 1)){ 
							 // no visible lesion -- Thermal Ablation

								GetTreatment.out[AgeGroup][CurrYear - StartYear] += 1;
								if (TrueStage == 0)
								{
									RSApop.NewUnnecessary[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
								}
								RSApop.NewThermal[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;

								if (HIVstage == 0)
								{
									if (res < 0.652)
									{
										for (xx = 0; xx < 13; xx++)
										{
											if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
												HPVstage[xx] == 3 || HPVstage[xx] == 4)
											{
												HPVstageE[xx] = 0;
												if (HPVstage[xx] == 4)
												{
													WeibullCIN3[xx] = 0;
													TimeinCIN3[xx] = 0;
												}
											}
										}
									}
									else
									{
										for (xx = 0; xx < 13; xx++)
										{
											if (HPVstage[xx] == 2)
											{
												HPVstageE[xx] = 1;
											}
											if (HPVstage[xx] == 3)
											{
												HPVstageE[xx] = 2;
											}
											if (HPVstage[xx] == 4)
											{
												HPVstageE[xx] = 3;
												WeibullCIN3[xx] = 0;
												TimeinCIN3[xx] = 0;
											}
										}
									}
								}
								else
								{
									if (res < 0.385)
									{
										for (xx = 0; xx < 13; xx++)
										{
											if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
												HPVstage[xx] == 3 || HPVstage[xx] == 4)
											{
												HPVstageE[xx] = 0;
												if (HPVstage[xx] == 4)
												{
													WeibullCIN3[xx] = 0;
													TimeinCIN3[xx] = 0;
												}
											}
										}
									}
									else
									{
										for (xx = 0; xx < 13; xx++)
										{
											if (HPVstage[xx] == 2)
											{
												HPVstageE[xx] = 1;
											}
											if (HPVstage[xx] == 3)
											{
												HPVstageE[xx] = 2;
											}
											if (HPVstage[xx] == 4)
											{
												HPVstageE[xx] = 3;
												WeibullCIN3[xx] = 0;
												TimeinCIN3[xx] = 0;
											}
										}
									}
								}

								if (PerfectSchedule == 0 || CurrYear < ImplementYR)
								{
									if (HIVstage == 5)
									{
										timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
									}
									else
									{
										timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
									}
									if (timetoScreen == 0)
									{
										timetoScreen = 1;
									}
								}
								if (PerfectSchedule == 1 && CurrYear >= ImplementYR)
								{
									timetoScreen = 48;
								}
								repeat = 1;
								HPVrepeat = 1;
							}

							else
							{ // lesion seen -- Colposcopy + LLETZ
								RSApop.GetReferred[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
								RSApop.NewLLETZ[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
								GetTreatment.out[AgeGroup][CurrYear - StartYear] += 1;
								if (TrueStage == 0)
								{
									RSApop.NewUnnecessary[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
								}

								if (HIVstage == 0)
								{
									if (res < 0.752)
									{
										if (ttC < 0.15)
										{
											for (xx = 0; xx < 13; xx++)
											{
												if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4)
												{
													HPVstageE[xx] = 1;
													if (HPVstage[xx] == 4)
													{
														WeibullCIN3[xx] = 0;
														TimeinCIN3[xx] = 0;
													}
												}
											}
										}
										else
										{
											for (xx = 0; xx < 13; xx++)
											{
												if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
													HPVstage[xx] == 3 || HPVstage[xx] == 4)
												{
													HPVstageE[xx] = 0;
													if (HPVstage[xx] == 4)
													{
														WeibullCIN3[xx] = 0;
														TimeinCIN3[xx] = 0;
													}
												}
											}
										}
									}
									else
									{
										if (CCd < 0.5)
										{
											for (xx = 0; xx < 13; xx++)
											{
												if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
													HPVstage[xx] == 3 || HPVstage[xx] == 4)
												{
													HPVstageE[xx] = 1;
													if (HPVstage[xx] == 4)
													{
														WeibullCIN3[xx] = 0;
														TimeinCIN3[xx] = 0;
													}
												}
											}
										}
									}
								}
								else
								{ // HIV positive
									if (res < 0.4)
									{
										if (ttC < 0.15)
										{
											for (xx = 0; xx < 13; xx++)
											{
												if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4)
												{
													HPVstageE[xx] = 1;
													if (HPVstage[xx] == 4)
													{
														WeibullCIN3[xx] = 0;
														TimeinCIN3[xx] = 0;
													}
												}
											}
										}
										else
										{
											for (xx = 0; xx < 13; xx++)
											{
												if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
													HPVstage[xx] == 3 || HPVstage[xx] == 4)
												{
													HPVstageE[xx] = 0;
													if (HPVstage[xx] == 4)
													{
														WeibullCIN3[xx] = 0;
														TimeinCIN3[xx] = 0;
													}
												}
											}
										}
									}
									else
									{
										if (CCd < 0.5)
										{
											for (xx = 0; xx < 13; xx++)
											{
												if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
													HPVstage[xx] == 3 || HPVstage[xx] == 4)
												{
													HPVstageE[xx] = 1;
													if (HPVstage[xx] == 4)
													{
														WeibullCIN3[xx] = 0;
														TimeinCIN3[xx] = 0;
													}
												}
											}
										}
									}
								}

								if (PerfectSchedule == 0 || CurrYear < ImplementYR)
								{
									if (HIVstage == 5)
									{
										timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
									}
									else
									{
										timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
									}
									if (timetoScreen == 0)
									{
										timetoScreen = 1;
									}
								}
								if (PerfectSchedule == 1 && CurrYear >= ImplementYR)
								{
									timetoScreen = 48;
								}
								repeat = 1;
								HPVrepeat = 1;
							}
							}
						}
						else
						{
							// Scale-up not active -- skip VIA, go straight to Colposcopy + LLETZ
							RSApop.GetReferred[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
							RSApop.NewLLETZ[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
							GetTreatment.out[AgeGroup][CurrYear - StartYear] += 1;
							if (TrueStage == 0)
							{
								RSApop.NewUnnecessary[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
							}

							if (HIVstage == 0)
							{
								if (res < 0.752)
								{
									if (ttC < 0.15)
									{
										for (xx = 0; xx < 13; xx++)
										{
											if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4)
											{
												HPVstageE[xx] = 1;
												if (HPVstage[xx] == 4)
												{
													WeibullCIN3[xx] = 0;
													TimeinCIN3[xx] = 0;
												}
											}
										}
									}
									else
									{
										for (xx = 0; xx < 13; xx++)
										{
											if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
												HPVstage[xx] == 3 || HPVstage[xx] == 4)
											{
												HPVstageE[xx] = 0;
												if (HPVstage[xx] == 4)
												{
													WeibullCIN3[xx] = 0;
													TimeinCIN3[xx] = 0;
												}
											}
										}
									}
								}
								else
								{
									if (CCd < 0.5)
									{
										for (xx = 0; xx < 13; xx++)
										{
											if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
												HPVstage[xx] == 3 || HPVstage[xx] == 4)
											{
												HPVstageE[xx] = 1;
												if (HPVstage[xx] == 4)
												{
													WeibullCIN3[xx] = 0;
													TimeinCIN3[xx] = 0;
												}
											}
										}
									}
								}
							}
							else
							{ // HIV positive
								if (res < 0.4)
								{
									if (ttC < 0.15)
									{
										for (xx = 0; xx < 13; xx++)
										{
											if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4)
											{
												HPVstageE[xx] = 1;
												if (HPVstage[xx] == 4)
												{
													WeibullCIN3[xx] = 0;
													TimeinCIN3[xx] = 0;
												}
											}
										}
									}
									else
									{
										for (xx = 0; xx < 13; xx++)
										{
											if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
												HPVstage[xx] == 3 || HPVstage[xx] == 4)
											{
												HPVstageE[xx] = 0;
												if (HPVstage[xx] == 4)
												{
													WeibullCIN3[xx] = 0;
													TimeinCIN3[xx] = 0;
												}
											}
										}
									}
								}
								else
								{
									if (CCd < 0.5)
									{
										for (xx = 0; xx < 13; xx++)
										{
											if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
												HPVstage[xx] == 3 || HPVstage[xx] == 4)
											{
												HPVstageE[xx] = 1;
												if (HPVstage[xx] == 4)
												{
													WeibullCIN3[xx] = 0;
													TimeinCIN3[xx] = 0;
												}
											}
										}
									}
								}
							}

							if (PerfectSchedule == 0 || CurrYear < ImplementYR)
							{
								if (HIVstage == 5)
								{
									timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
								}
								else
								{
									timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
								}
								if (timetoScreen == 0)
								{
									timetoScreen = 1;
								}
							}
							if (PerfectSchedule == 1 && CurrYear >= ImplementYR)
							{
								timetoScreen = 48;
							}
							repeat = 1;
							HPVrepeat = 1;
						}
					}
				
					// HPV 16/18/45 with cancer (stages 3/4)
					else if (AnyHPV(HPVstage, hpv161845, cc_un34))
					{
						DiagnosedCC = 1;
						RSApop.NewDiagCancer[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
						if (AnyHPV(HPVstage, hpv1618, cc_un34))
						{
							RSApop.NewDiagCancer1618[AgeGroup][CurrYear - StartYear] += 1;
						}
						if (HIVstage == 5)
						{
							RSApop.NewDiagCancerART[CurrYear - StartYear] += 1;
						}
						for (xx = 0; xx < 13; xx++)
						{
							if (HPVstage[xx] == 5)
							{
								HPVstageE[xx] = 11;
								RSApop.StageDiag[0][CurrYear - StartYear] += 1;
								RSApop.StageIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
								if (SId < 0.192)
								{
									StageIdeath = static_cast<int>(48.0 * 3.08 * pow(-log(SI), 1.0 / 1.23));
								}
								else
								{
									StageIrecover = 8;
								}
							}
							if (HPVstage[xx] == 8)
							{
								HPVstageE[xx] = 12;
								RSApop.StageDiag[1][CurrYear - StartYear] += 1;
								RSApop.StageIIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
								if (SIId < 0.466)
								{
									StageIIdeath = static_cast<int>(48.0 * 2.39 * pow(-log(SII), 1.0 / 1.17));
								}
								else
								{
									StageIIrecover = 24;
								}
							}
							if (HPVstage[xx] == 9)
							{
								HPVstageE[xx] = 13;
								RSApop.StageDiag[2][CurrYear - StartYear] += 1;
								RSApop.StageIIIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
								if (SIIId < 0.715)
								{
									StageIIIdeath = static_cast<int>(48.0 * 1.18 * pow(-log(SIII), 1.0 / 0.91));
								}
								else
								{
									StageIIIrecover = 24;
								}
							}
							if (HPVstage[xx] == 10)
							{
								HPVstageE[xx] = 14;
								RSApop.StageDiag[3][CurrYear - StartYear] += 1;
								RSApop.StageIVdiag[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
								if (SIVd < 0.931)
								{
									StageIVdeath = static_cast<int>(48.0 * 0.46 * pow(-log(SIV), 1.0 / 0.9));
								}
								else
								{
									StageIVrecover = 24;
								}
							}
						}
					}

					// Other high-risk types -- Reflex Cytology 
					else if (AnyHPV(HPVstage, {2, 3, 4, 5, 7, 8, 9, 10, 11, 12}, {1, 2, 3, 4, 5, 8, 9, 10}))
					{

						RSApop.NewReflex[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;

						// Simulate reflex cytology result
						// 0 = Normal, 1 = ASCUS/LSIL, 2 = HSIL/ASC-H
						int ScreenResult = 0;

						if (HIVstage == 0)
						{
							if (TrueStage == 0)
							{
								if (res < 0.732)
								{
									ScreenResult = 0;
								}
								else if (res < 0.732 + 0.146)
								{
									ScreenResult = 1;
								}
								else
								{
									ScreenResult = 2;
								}
							}
							else if (TrueStage == 1)
							{
								if (res < 0.586)
								{
									ScreenResult = 0;
								}
								else if (res < 0.586 + 0.261)
								{
									ScreenResult = 1;
								}
								else
								{
									ScreenResult = 2;
								}
							}
							else if (TrueStage >= 2 && TrueStage < 5)
							{
								if (res < 0.453)
								{
									ScreenResult = 0;
								}
								else if (res < 0.453 + 0.203)
								{
									ScreenResult = 1;
								}
								else
								{
									ScreenResult = 2;
								}
							}
						}
						else
						{ // HIV positive
							if (TrueStage == 0)
							{
								if (res < 0.62)
								{
									ScreenResult = 0;
								}
								else if (res < 0.62 + 0.16)
								{
									ScreenResult = 1;
								}
								else
								{
									ScreenResult = 2;
								}
							}
							else if (TrueStage == 1)
							{
								if (res < 0.545)
								{
									ScreenResult = 0;
								}
								else if (res < 0.545 + 0.15)
								{
									ScreenResult = 1;
								}
								else
								{
									ScreenResult = 2;
								}
							}
							else if (TrueStage >= 2 && TrueStage < 5)
							{
								if (res < 0.209)
								{
									ScreenResult = 0;
								}
								else if (res < 0.209 + 0.157)
								{
									ScreenResult = 1;
								}
								else
								{
									ScreenResult = 2;
								}
							}
						}

						// Cancer picked up by cytology
						if (TrueStage == 3 && ScreenResult == 2 && CCd < 0.35)
						{
							DiagnosedCC = 1;
							RSApop.NewDiagCancer[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
							if (HIVstage == 5)
							{
								RSApop.NewDiagCancerART[CurrYear - StartYear] += 1;
							}
							for (xx = 0; xx < 13; xx++)
							{
								if (HPVstage[xx] == 5)
								{
									HPVstageE[xx] = 11;
									RSApop.StageDiag[0][CurrYear - StartYear] += 1;
									RSApop.StageIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
									if (SId < 0.192)
									{
										StageIdeath = static_cast<int>(48.0 * 3.08 * pow(-log(SI), 1.0 / 1.23));
									}
									else
									{
										StageIrecover = 8;
									}
								}
								if (HPVstage[xx] == 8)
								{
									HPVstageE[xx] = 12;
									RSApop.StageDiag[1][CurrYear - StartYear] += 1;
									RSApop.StageIIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
									if (SIId < 0.466)
									{
										StageIIdeath = static_cast<int>(48.0 * 2.39 * pow(-log(SII), 1.0 / 1.17));
									}
									else
									{
										StageIIrecover = 24;
									}
								}
								if (HPVstage[xx] == 9)
								{
									HPVstageE[xx] = 13;
									RSApop.StageDiag[2][CurrYear - StartYear] += 1;
									RSApop.StageIIIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
									if (SIIId < 0.715)
									{
										StageIIIdeath = static_cast<int>(48.0 * 1.18 * pow(-log(SIII), 1.0 / 0.91));
									}
									else
									{
										StageIIIrecover = 24;
									}
								}
								if (HPVstage[xx] == 10)
								{
									HPVstageE[xx] = 14;
									RSApop.StageDiag[3][CurrYear - StartYear] += 1;
									RSApop.StageIVdiag[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
									if (SIVd < 0.931)
									{
										StageIVdeath = static_cast<int>(48.0 * 0.46 * pow(-log(SIV), 1.0 / 0.9));
									}
									else
									{
										StageIVrecover = 24;
									}
								}
							}
						}

						// ----- Cytology Normal -- HPV follow-up at 12 months -----
						else if (ScreenResult == 0)
						{
							if (PerfectSchedule == 0 || CurrYear < ImplementYR)
							{
								if (HIVstage == 5)
								{
									timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
								}
								else
								{
									timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
								}
								if (timetoScreen == 0)
								{
									timetoScreen = 1;
								}
							}
							if (PerfectSchedule == 1 && CurrYear >= ImplementYR)
							{
								timetoScreen = 48;
							}
							repeat = 1;
							HPVrepeat = 1;
						}

						//  ASCUS / LSIL -- VIA + Thermal Ablation (with scale-up) 
						else if (ScreenResult == 1)
						{

							double TA_coverage = 0.0;
							int iy = CurrYear - StartYear;
							if (iy >= 0 && iy < 136)
							{
								TA_coverage = TAScaleUpRate[iy];
							}

							if (taRand < TA_coverage)
							{
								// Full pathway available -- do VIA
								RSApop.NewVAT[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;

								if (viaVis < 0.65)
								{ // no visible lesion -- Thermal Ablation
									GetTreatment.out[AgeGroup][CurrYear - StartYear] += 1;
									if (TrueStage == 0)
									{
										RSApop.NewUnnecessary[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
									}
									RSApop.NewThermal[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;

									// Thermal efficacy (same numbers as 16/18/45 arm)
									if (HIVstage == 0)
									{
										if (res < 0.652)
										{
											for (xx = 0; xx < 13; xx++)
											{
												if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
													HPVstage[xx] == 3 || HPVstage[xx] == 4)
												{
													HPVstageE[xx] = 0;
													if (HPVstage[xx] == 4)
													{
														WeibullCIN3[xx] = 0;
														TimeinCIN3[xx] = 0;
													}
												}
											}
										}
										else
										{
											for (xx = 0; xx < 13; xx++)
											{
												if (HPVstage[xx] == 2)
													HPVstageE[xx] = 1;
												if (HPVstage[xx] == 3)
													HPVstageE[xx] = 2;
												if (HPVstage[xx] == 4)
												{
													HPVstageE[xx] = 3;
													WeibullCIN3[xx] = 0;
													TimeinCIN3[xx] = 0;
												}
											}
										}
									}
									else
									{
										if (res < 0.385)
										{
											for (xx = 0; xx < 13; xx++)
											{
												if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
													HPVstage[xx] == 3 || HPVstage[xx] == 4)
												{
													HPVstageE[xx] = 0;
													if (HPVstage[xx] == 4)
													{
														WeibullCIN3[xx] = 0;
														TimeinCIN3[xx] = 0;
													}
												}
											}
										}
										else
										{
											for (xx = 0; xx < 13; xx++)
											{
												if (HPVstage[xx] == 2)
													HPVstageE[xx] = 1;
												if (HPVstage[xx] == 3)
													HPVstageE[xx] = 2;
												if (HPVstage[xx] == 4)
												{
													HPVstageE[xx] = 3;
													WeibullCIN3[xx] = 0;
													TimeinCIN3[xx] = 0;
												}
											}
										}
									}

									// Follow-up HPV at 12 months
									if (PerfectSchedule == 0 || CurrYear < ImplementYR)
									{
										if (HIVstage == 5)
										{
											timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
										}
										else
										{
											timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
										}
										if (timetoScreen == 0)
											timetoScreen = 1;
									}
									if (PerfectSchedule == 1 && CurrYear >= ImplementYR)
									{
										timetoScreen = 48;
									}
									repeat = 1;
									HPVrepeat = 1;
								}
								else
								{
									// Visible lesion -- Colposcopy + LLETZ
									RSApop.GetReferred[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
									RSApop.NewLLETZ[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
									GetTreatment.out[AgeGroup][CurrYear - StartYear] += 1;
									if (TrueStage == 0)
									{
										RSApop.NewUnnecessary[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
									}

									// LLETZ efficacy block
									if (HIVstage == 0)
									{
										if (res < 0.752)
										{
											if (ttC < 0.15)
											{
												for (xx = 0; xx < 13; xx++)
												{
													if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4)
													{
														HPVstageE[xx] = 1;
														if (HPVstage[xx] == 4)
														{
															WeibullCIN3[xx] = 0;
															TimeinCIN3[xx] = 0;
														}
													}
												}
											}
											else
											{
												for (xx = 0; xx < 13; xx++)
												{
													if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
														HPVstage[xx] == 3 || HPVstage[xx] == 4)
													{
														HPVstageE[xx] = 0;
														if (HPVstage[xx] == 4)
														{
															WeibullCIN3[xx] = 0;
															TimeinCIN3[xx] = 0;
														}
													}
												}
											}
										}
										else
										{
											if (CCd < 0.5)
											{
												for (xx = 0; xx < 13; xx++)
												{
													if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
														HPVstage[xx] == 3 || HPVstage[xx] == 4)
													{
														HPVstageE[xx] = 1;
														if (HPVstage[xx] == 4)
														{
															WeibullCIN3[xx] = 0;
															TimeinCIN3[xx] = 0;
														}
													}
												}
											}
										}
									}
									else
									{
										if (res < 0.4)
										{
											if (ttC < 0.15)
											{
												for (xx = 0; xx < 13; xx++)
												{
													if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4)
													{
														HPVstageE[xx] = 1;
														if (HPVstage[xx] == 4)
														{
															WeibullCIN3[xx] = 0;
															TimeinCIN3[xx] = 0;
														}
													}
												}
											}
											else
											{
												for (xx = 0; xx < 13; xx++)
												{
													if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
														HPVstage[xx] == 3 || HPVstage[xx] == 4)
													{
														HPVstageE[xx] = 0;
														if (HPVstage[xx] == 4)
														{
															WeibullCIN3[xx] = 0;
															TimeinCIN3[xx] = 0;
														}
													}
												}
											}
										}
										else
										{
											if (CCd < 0.5)
											{
												for (xx = 0; xx < 13; xx++)
												{
													if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
														HPVstage[xx] == 3 || HPVstage[xx] == 4)
													{
														HPVstageE[xx] = 1;
														if (HPVstage[xx] == 4)
														{
															WeibullCIN3[xx] = 0;
															TimeinCIN3[xx] = 0;
														}
													}
												}
											}
										}
									}

									if (PerfectSchedule == 0 || CurrYear < ImplementYR)
									{
										if (HIVstage == 5)
										{
											timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
										}
										else
										{
											timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
										}
										if (timetoScreen == 0)
											timetoScreen = 1;
									}
									if (PerfectSchedule == 1 && CurrYear >= ImplementYR)
									{
										timetoScreen = 48;
									}
									repeat = 1;
									HPVrepeat = 1;
								}
							}
							else
							{
								// Scale-up not active -- skip VIA, go straight to Colposcopy + LLETZ
								RSApop.GetReferred[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
								RSApop.NewLLETZ[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
								GetTreatment.out[AgeGroup][CurrYear - StartYear] += 1;
								if (TrueStage == 0)
								{
									RSApop.NewUnnecessary[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
								}

								// LLETZ efficacy block (same as above)
								if (HIVstage == 0)
								{
									if (res < 0.752)
									{
										if (ttC < 0.15)
										{
											for (xx = 0; xx < 13; xx++)
											{
												if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4)
												{
													HPVstageE[xx] = 1;
													if (HPVstage[xx] == 4)
													{
														WeibullCIN3[xx] = 0;
														TimeinCIN3[xx] = 0;
													}
												}
											}
										}
										else
										{
											for (xx = 0; xx < 13; xx++)
											{
												if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
													HPVstage[xx] == 3 || HPVstage[xx] == 4)
												{
													HPVstageE[xx] = 0;
													if (HPVstage[xx] == 4)
													{
														WeibullCIN3[xx] = 0;
														TimeinCIN3[xx] = 0;
													}
												}
											}
										}
									}
									else
									{
										if (CCd < 0.5)
										{
											for (xx = 0; xx < 13; xx++)
											{
												if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
													HPVstage[xx] == 3 || HPVstage[xx] == 4)
												{
													HPVstageE[xx] = 1;
													if (HPVstage[xx] == 4)
													{
														WeibullCIN3[xx] = 0;
														TimeinCIN3[xx] = 0;
													}
												}
											}
										}
									}
								}
								else
								{
									if (res < 0.4)
									{
										if (ttC < 0.15)
										{
											for (xx = 0; xx < 13; xx++)
											{
												if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4)
												{
													HPVstageE[xx] = 1;
													if (HPVstage[xx] == 4)
													{
														WeibullCIN3[xx] = 0;
														TimeinCIN3[xx] = 0;
													}
												}
											}
										}
										else
										{
											for (xx = 0; xx < 13; xx++)
											{
												if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
													HPVstage[xx] == 3 || HPVstage[xx] == 4)
												{
													HPVstageE[xx] = 0;
													if (HPVstage[xx] == 4)
													{
														WeibullCIN3[xx] = 0;
														TimeinCIN3[xx] = 0;
													}
												}
											}
										}
									}
									else
									{
										if (CCd < 0.5)
										{
											for (xx = 0; xx < 13; xx++)
											{
												if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
													HPVstage[xx] == 3 || HPVstage[xx] == 4)
												{
													HPVstageE[xx] = 1;
													if (HPVstage[xx] == 4)
													{
														WeibullCIN3[xx] = 0;
														TimeinCIN3[xx] = 0;
													}
												}
											}
										}
									}
								}

								if (PerfectSchedule == 0 || CurrYear < ImplementYR)
								{
									if (HIVstage == 5)
									{
										timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
									}
									else
									{
										timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
									}
									if (timetoScreen == 0)
										timetoScreen = 1;
								}
								if (PerfectSchedule == 1 && CurrYear >= ImplementYR)
								{
									timetoScreen = 48;
								}
								repeat = 1;
								HPVrepeat = 1;
							}
						}

						//  HSIL / ASC-H --  Colposcopy + LLETZ 
						else if (ScreenResult == 2 && TrueStage < 3)
						{
							RSApop.GetReferred[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
							RSApop.NewLLETZ[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
							GetTreatment.out[AgeGroup][CurrYear - StartYear] += 1;
							if (TrueStage == 0)
							{
								RSApop.NewUnnecessary[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
							}

							// LLETZ efficacy block
							if (HIVstage == 0)
							{
								if (res < 0.752)
								{
									if (ttC < 0.15)
									{
										for (xx = 0; xx < 13; xx++)
										{
											if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4)
											{
												HPVstageE[xx] = 1;
												if (HPVstage[xx] == 4)
												{
													WeibullCIN3[xx] = 0;
													TimeinCIN3[xx] = 0;
												}
											}
										}
									}
									else
									{
										for (xx = 0; xx < 13; xx++)
										{
											if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
												HPVstage[xx] == 3 || HPVstage[xx] == 4)
											{
												HPVstageE[xx] = 0;
												if (HPVstage[xx] == 4)
												{
													WeibullCIN3[xx] = 0;
													TimeinCIN3[xx] = 0;
												}
											}
										}
									}
								}
								else
								{
									if (CCd < 0.5)
									{
										for (xx = 0; xx < 13; xx++)
										{
											if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
												HPVstage[xx] == 3 || HPVstage[xx] == 4)
											{
												HPVstageE[xx] = 1;
												if (HPVstage[xx] == 4)
												{
													WeibullCIN3[xx] = 0;
													TimeinCIN3[xx] = 0;
												}
											}
										}
									}
								}
							}
							else
							{
								if (res < 0.4)
								{
									if (ttC < 0.15)
									{
										for (xx = 0; xx < 13; xx++)
										{
											if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4)
											{
												HPVstageE[xx] = 1;
												if (HPVstage[xx] == 4)
												{
													WeibullCIN3[xx] = 0;
													TimeinCIN3[xx] = 0;
												}
											}
										}
									}
									else
									{
										for (xx = 0; xx < 13; xx++)
										{
											if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
												HPVstage[xx] == 3 || HPVstage[xx] == 4)
											{
												HPVstageE[xx] = 0;
												if (HPVstage[xx] == 4)
												{
													WeibullCIN3[xx] = 0;
													TimeinCIN3[xx] = 0;
												}
											}
										}
									}
								}
								else
								{
									if (CCd < 0.5)
									{
										for (xx = 0; xx < 13; xx++)
										{
											if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
												HPVstage[xx] == 3 || HPVstage[xx] == 4)
											{
												HPVstageE[xx] = 1;
												if (HPVstage[xx] == 4)
												{
													WeibullCIN3[xx] = 0;
													TimeinCIN3[xx] = 0;
												}
											}
										}
									}
								}
							}

							if (PerfectSchedule == 0 || CurrYear < ImplementYR)
							{
								if (HIVstage == 5)
								{
									timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
								}
								else
								{
									timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
								}
								if (timetoScreen == 0)
									timetoScreen = 1;
							}
							if (PerfectSchedule == 1 && CurrYear >= ImplementYR)
							{
								timetoScreen = 48;
							}
							repeat = 1;
							HPVrepeat = 1;
						}
					}
				}

				else { // HPVrepeat == 1  --  12-month post-follow up HPV test
					HPVrepeat = 0;
				
					if (HPVstatus == 0) {
						// HPV negative at 12-month follow-up -- return to normal long interval
						if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
							if (HIVstage == 5) {
								timetoScreen = 5 * 48;          // every 5 years
							} else {
								timetoScreen = 10 * 48;         // every 10 years
							}
						} else {
							if (HIVstage == 5) {
								timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
							} else if (AgeExact < 50) {
								timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
							} else {
								timetoScreen = 200 * 48;
							}
							if (timetoScreen == 0) timetoScreen = 1;
						}
						repeat = 0;
					}
					else {
						// Still HPV positive at 12-month follow-up -- treat immediately
						// (Colposcopy + LLETZ as per Figure 6)
						RSApop.GetReferred[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
						RSApop.NewLLETZ[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
						GetTreatment.out[AgeGroup][CurrYear - StartYear] += 1;
						if (TrueStage == 0) {
							RSApop.NewUnnecessary[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
						}
				
						// LLETZ efficacy block 
						if (HIVstage == 0) {
							if (res < 0.752) {
								if (ttC < 0.15) {
									for (xx = 0; xx < 13; xx++) {
										if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4) {
											HPVstageE[xx] = 1;
											if (HPVstage[xx] == 4) {
												WeibullCIN3[xx] = 0;
												TimeinCIN3[xx] = 0;
											}
										}
									}
								} else {
									for (xx = 0; xx < 13; xx++) {
										if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
											HPVstage[xx] == 3 || HPVstage[xx] == 4) {
											HPVstageE[xx] = 0;
											if (HPVstage[xx] == 4) {
												WeibullCIN3[xx] = 0;
												TimeinCIN3[xx] = 0;
											}
										}
									}
								}
							} else {
								if (CCd < 0.5) {
									for (xx = 0; xx < 13; xx++) {
										if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
											HPVstage[xx] == 3 || HPVstage[xx] == 4) {
											HPVstageE[xx] = 1;
											if (HPVstage[xx] == 4) {
												WeibullCIN3[xx] = 0;
												TimeinCIN3[xx] = 0;
											}
										}
									}
								}
							}
						} else { // HIV positive
							if (res < 0.4) {
								if (ttC < 0.15) {
									for (xx = 0; xx < 13; xx++) {
										if (HPVstage[xx] == 2 || HPVstage[xx] == 3 || HPVstage[xx] == 4) {
											HPVstageE[xx] = 1;
											if (HPVstage[xx] == 4) {
												WeibullCIN3[xx] = 0;
												TimeinCIN3[xx] = 0;
											}
										}
									}
								} else {
									for (xx = 0; xx < 13; xx++) {
										if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
											HPVstage[xx] == 3 || HPVstage[xx] == 4) {
											HPVstageE[xx] = 0;
											if (HPVstage[xx] == 4) {
												WeibullCIN3[xx] = 0;
												TimeinCIN3[xx] = 0;
											}
										}
									}
								}
							} else {
								if (CCd < 0.5) {
									for (xx = 0; xx < 13; xx++) {
										if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
											HPVstage[xx] == 3 || HPVstage[xx] == 4) {
											HPVstageE[xx] = 1;
											if (HPVstage[xx] == 4) {
												WeibullCIN3[xx] = 0;
												TimeinCIN3[xx] = 0;
											}
										}
									}
								}
							}
						}
				
						// After treatment, schedule another HPV in 12 months
						if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
							timetoScreen = 48;
						} else {
							if (HIVstage == 5) {
								timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
							} else {
								timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
							}
							if (timetoScreen == 0) timetoScreen = 1;
						}
						repeat = 1;
						HPVrepeat = 1;   // will come back again in 12 months
					}
				}
			}
		}
		else
		{ // did not come back for results
			if ((WHOScreening == 0 && PerfectSchedule == 0) || CurrYear < ImplementYR)
			{
				if (HIVstage == 5)
				{
					timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
				}
				else if (AgeExact < 50)
				{
					timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
				}
				else
				{
					timetoScreen = 200 * 48;
				}
				if (timetoScreen == 0)
				{
					timetoScreen = 1;
				}
			}
			if (PerfectSchedule == 1 && CurrYear >= ImplementYR)
			{
				if (HIVstage == 5)
				{
					timetoScreen = 5 * 48;
				}
				else if (AgeExact < 55)
				{
					timetoScreen = 10 * 48;
				}
				else
				{
					timetoScreen = 200 * 48;
				}
			}
			repeat = 0;
		}
	} else { // inadequate sample
		if (PerfectSchedule == 0 || CurrYear < ImplementYR)
		{
			timetoScreen = 31.2 * pow(-log(tts), (1.0 / 0.57)) * 12;
			if (timetoScreen == 0)
			{
				timetoScreen = 1;
			}
		}
		if (PerfectSchedule == 1 && CurrYear >= ImplementYR)
		{
			timetoScreen = 12;
		}
		repeat = 1;
	}
}

  void Indiv::HPVScreenAlgorithm(
      int ID, double rea, double ade, double tts, double res, double ttC,
      double CCd, double SI, double SII, double SIII, double SIV, double SId,
      double SIId, double SIIId, double SIVd) {
    int xx, yy, zz;
    int SimCount2 = (CurrSim - 1) / IterationsPerPC;
    // cout << "do this!" << endl;
    if (HIVstage == 0) {
      zz = 0;
    } else if (HIVstage == 5) {
      zz = 1;
    } //||HIVstage==6
    else {
      zz = 2;
    }

    RSApop.NewHPVScreen[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;

    if (ade < 0.95) {

      if (HPVstatus == 0) {
        if ((WHOScreening == 0 && PerfectSchedule == 0) ||
            CurrYear < ImplementYR) {
          // if(HIVstage==5) {timetoScreen = 5.3 * pow(-log(tts),(1.0/0.78)) *
          // 48;} else if(AgeExact<50) { timetoScreen = 15.0 *
          // pow(-log(tts),(1.0/0.83)) * 48; }
          if (HIVstage == 5) {
            timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
          } else if (AgeExact < 50) {
            timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
          } else {
            timetoScreen = 200 * 48;
          }
          if (timetoScreen == 0) {
            timetoScreen = 1;
          }
        }
        if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
          if (HIVstage == 5) {
            timetoScreen = 3 * 48;
          } else if (AgeExact < 50) {
            timetoScreen = 10 * 48;
          } else {
            timetoScreen = 200 * 48;
          }
        }
        repeat = 0;
        HPVrepeat = 0;
      }

      else if (rea < 0.9) { // 90% come back for results
        if (HPVGenotyping == 1) {
          if (HPVrepeat == 0) {
            // if test positive with HPV16/18/45 (but not cancer), treat
            if (AnyHPV(HPVstage, hpv161845, {1, 2, 3, 4}) && TrueStage < 3) {
              GetTreatment.out[AgeGroup][CurrYear - StartYear] += 1;
              if (TrueStage == 0) {
                RSApop.NewUnnecessary[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
              }
              RSApop.NewLLETZ[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
              if (HIVstage == 0) {
                if (res < 0.752) {
                  if (ttC < 0.15) {
                    for (xx = 0; xx < 13; xx++) {
                      if (HPVstage[xx] == 2 || HPVstage[xx] == 3 ||
                          HPVstage[xx] == 4) {
                        HPVstageE[xx] = 1;
                        if (HPVstage[xx] == 4) {
                          WeibullCIN3[xx] = 0;
                          TimeinCIN3[xx] = 0;
                        }
                      }
                    }
                  } else {
                    for (xx = 0; xx < 13; xx++) {
                      if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                          HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                        HPVstageE[xx] = 0;
                        if (HPVstage[xx] == 4) {
                          WeibullCIN3[xx] = 0;
                          TimeinCIN3[xx] = 0;
                        }
                      }
                    }
                  }
                } else {
                  if (CCd < 0.5) {
                    for (xx = 0; xx < 13; xx++) {
                      if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                          HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                        HPVstageE[xx] = 1;
                        if (HPVstage[xx] == 4) {
                          WeibullCIN3[xx] = 0;
                          TimeinCIN3[xx] = 0;
                        }
                      }
                    }
                  }
                }
              } else {
                if (res < 0.4) {
                  if (ttC < 0.15) {
                    for (xx = 0; xx < 13; xx++) {
                      if (HPVstage[xx] == 2 || HPVstage[xx] == 3 ||
                          HPVstage[xx] == 4) {
                        HPVstageE[xx] = 1;
                        if (HPVstage[xx] == 4) {
                          WeibullCIN3[xx] = 0;
                          TimeinCIN3[xx] = 0;
                        }
                      }
                    }
                  } else {
                    for (xx = 0; xx < 13; xx++) {
                      if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                          HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                        HPVstageE[xx] = 0;
                        if (HPVstage[xx] == 4) {
                          WeibullCIN3[xx] = 0;
                          TimeinCIN3[xx] = 0;
                        }
                      }
                    }
                  }
                } else {
                  if (CCd < 0.5) {
                    for (xx = 0; xx < 13; xx++) {
                      if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                          HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                        HPVstageE[xx] = 1;
                        if (HPVstage[xx] == 4) {
                          WeibullCIN3[xx] = 0;
                          TimeinCIN3[xx] = 0;
                        }
                      }
                    }
                  }
                }
              }

              if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                if (HIVstage == 5) {
                  timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
                } else {
                  timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
                }
                if (timetoScreen == 0) {
                  timetoScreen = 1;
                }
              }
              if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
                timetoScreen = 48;
              }
              repeat = 1;
            }
            // if test positive with HPV16/18/45 (cancer), 'diagnose'
            else if (AnyHPV(HPVstage, hpv161845, cc_un34)) {
              DiagnosedCC = 1;
              RSApop.NewDiagCancer[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                  1;
              if (AnyHPV(HPVstage, hpv1618, cc_un34)) {
                RSApop.NewDiagCancer1618[AgeGroup][CurrYear - StartYear] += 1;
              }
              if (HIVstage == 5) {
                RSApop.NewDiagCancerART[CurrYear - StartYear] += 1;
              }
              for (xx = 0; xx < 13; xx++) {
                if (HPVstage[xx] == 5) {
                  HPVstageE[xx] = 11;
                  RSApop.StageDiag[0][CurrYear - StartYear] += 1;
                  RSApop.StageIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                      1;
                  if (SId < 0.192) {
                    StageIdeath = static_cast<int>(48.0 * 3.08 *
                                                   pow(-log(SI), 1.0 / 1.23));
                  } else {
                    StageIrecover = 8;
                  }
                }
                if (HPVstage[xx] == 8) {
                  HPVstageE[xx] = 12;
                  RSApop.StageDiag[1][CurrYear - StartYear] += 1;
                  RSApop
                      .StageIIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                      1;
                  if (SIId < 0.466) {
                    StageIIdeath = static_cast<int>(48.0 * 2.39 *
                                                    pow(-log(SII), 1.0 / 1.17));
                  } else {
                    StageIIrecover = 24;
                  }
                }
                if (HPVstage[xx] == 9) {
                  HPVstageE[xx] = 13;
                  RSApop.StageDiag[2][CurrYear - StartYear] += 1;
                  RSApop
                      .StageIIIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                      1;
                  if (SIIId < 0.715) {
                    StageIIIdeath = static_cast<int>(
                        48.0 * 1.18 * pow(-log(SIII), 1.0 / 0.91));
                  } else {
                    StageIIIrecover = 24;
                  }
                }
                if (HPVstage[xx] == 10) {
                  HPVstageE[xx] = 14;
                  RSApop.StageDiag[3][CurrYear - StartYear] += 1;
                  RSApop
                      .StageIVdiag[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                      1;
                  if (SIVd < 0.931) {
                    StageIVdeath = static_cast<int>(48.0 * 0.46 *
                                                    pow(-log(SIV), 1.0 / 0.9));
                  } else {
                    StageIVrecover = 24;
                  }
                }
              }
            }

            // if test positive for 31/33/35/52/58, triage, treat
            else if (AnyHPV(HPVstage, {2, 3, 4, 8, 10},
                            {1, 2, 3, 4, 5, 8, 9, 10})) {
              if (PapTRIAGE == 0) {
                // Refer to colposcopy
                RSApop.GetReferred[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                    1;
                if (HIVstage == 5) {
                  if (ttC < AttendColposcopy[2][CurrYear - StartYear]) {
                    timetoCol = 24;
                    timetoScreen = 0;
                  } else {
                    timetoCol = 0;
                    if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                      // timetoScreen = 5.3 * pow(-log(tts),(1.0/0.78)) * 48;
                      timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
                      if (timetoScreen == 0) {
                        timetoScreen = 1;
                      }
                      repeat = 1;
                    }
                  }
                }
                if (HIVstage > 0 && HIVstage != 5) {
                  if (ttC < AttendColposcopy[1][CurrYear - StartYear]) {
                    timetoCol = 24;
                    timetoScreen = 0;
                  } else {
                    timetoCol = 0;
                    if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                      // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48;
                      timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                      if (timetoScreen == 0) {
                        timetoScreen = 1;
                      }
                      repeat = 1;
                    }
                  }
                }
                if (HIVstage == 0) {
                  if (ttC < AttendColposcopy[0][CurrYear - StartYear]) {
                    timetoCol = 24;
                    timetoScreen = 0;
                  } else {
                    timetoCol = 0;
                    if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                      // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48;
                      timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                      if (timetoScreen == 0) {
                        timetoScreen = 1;
                      }
                      repeat = 1;
                    }
                  }
                }
              } else {
                RSApop.NewScreen[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
                if (ade < PapAdequacy[CurrYear - StartYear]) {
                  if (HIVstage == 0) {
                    if (TrueStage < 2) {
                      if (res < 0.75) {
                        ScreenResult = 0;
                      } else {
                        ScreenResult = 2;
                      }
                    }
                    if (TrueStage >= 2 && TrueStage < 5) {
                      if (res < 0.72) {
                        ScreenResult = 2;
                      } else {
                        ScreenResult = 0;
                      }
                    }
                  }
                  if (HIVstage > 0) {
                    if (TrueStage < 2) {
                      if (res < 0.44) {
                        ScreenResult = 0;
                      } else {
                        ScreenResult = 2;
                      }
                    }
                    if (TrueStage >= 2 && TrueStage < 5) {
                      if (res < 0.92) {
                        ScreenResult = 2;
                      } else {
                        ScreenResult = 0;
                      }
                    }
                  }

                  if (TrueStage == 3 && ScreenResult == 2 &&
                      CCd <
                          0.35) { // CCd is sensitivity of Pap to pick up cancer
                    DiagnosedCC = 1;
                    RSApop.NewDiagCancer[zz * 18 + AgeGroup]
                                        [CurrYear - StartYear] += 1;
                    if (HIVstage == 5) {
                      RSApop.NewDiagCancerART[CurrYear - StartYear] += 1;
                    }
                    for (int xx = 0; xx < 13; xx++) {
                      if (HPVstage[xx] == 5) {
                        HPVstageE[xx] = 11;
                        RSApop.StageDiag[0][CurrYear - StartYear] += 1;
                        RSApop.StageIdiag[zz * 18 + AgeGroup]
                                         [CurrYear - StartYear] += 1;
                        if (SId < 0.192) {
                          StageIdeath = static_cast<int>(
                              48.0 * 3.08 * pow(-log(SI), 1.0 / 1.23));
                        } else {
                          StageIrecover = 8;
                        }
                      }
                      if (HPVstage[xx] == 8) {
                        HPVstageE[xx] = 12;
                        RSApop.StageDiag[1][CurrYear - StartYear] += 1;
                        RSApop.StageIIdiag[zz * 18 + AgeGroup]
                                          [CurrYear - StartYear] += 1;
                        if (SIId < 0.466) {
                          StageIIdeath = static_cast<int>(
                              48.0 * 2.39 * pow(-log(SII), 1.0 / 1.17));
                        } else {
                          StageIIrecover = 24;
                        }
                      }
                      if (HPVstage[xx] == 9) {
                        HPVstageE[xx] = 13;
                        RSApop.StageDiag[2][CurrYear - StartYear] += 1;
                        RSApop.StageIIIdiag[zz * 18 + AgeGroup]
                                           [CurrYear - StartYear] += 1;
                        if (SIIId < 0.715) {
                          StageIIIdeath = static_cast<int>(
                              48.0 * 1.18 * pow(-log(SIII), 1.0 / 0.91));
                        } else {
                          StageIIIrecover = 24;
                        }
                      }
                      if (HPVstage[xx] == 10) {
                        HPVstageE[xx] = 14;
                        RSApop.StageDiag[3][CurrYear - StartYear] += 1;
                        RSApop.StageIVdiag[zz * 18 + AgeGroup]
                                          [CurrYear - StartYear] += 1;
                        if (SIVd < 0.931) {
                          StageIVdeath = static_cast<int>(
                              48.0 * 0.46 * pow(-log(SIV), 1.0 / 0.9));
                        } else {
                          StageIVrecover = 24;
                        }
                      }
                    }
                  }

                  // Normal screen:
                  if (ScreenResult == 0) {
                    if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                      // if(HIVstage==5) {timetoScreen = 5.3 *
                      // pow(-log(tts),(1.0/0.78)) * 48;} else if(AgeExact<50) {
                      // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48; }
                      if (HIVstage == 5) {
                        timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
                      } else if (AgeExact < 50) {
                        timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                      } else {
                        timetoScreen = 200 * 48;
                      }
                      if (timetoScreen == 0) {
                        timetoScreen = 1;
                      }
                    }
                    if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
                      if (HIVstage == 5) {
                        timetoScreen = 3 * 48;
                      } else if (AgeExact < 50) {
                        timetoScreen = 10 * 48;
                      } else {
                        timetoScreen = 200 * 48;
                      }
                    }
                    repeat = 0;
                  }
                  // HSIL/CC screen:
                  if (ScreenResult == 2 && TrueStage < 3 && SI < 0.9) {
                    GetTreatment.out[AgeGroup][CurrYear - StartYear] += 1;
                    if (TrueStage == 0) {
                      RSApop.NewUnnecessary[zz * 18 + AgeGroup]
                                           [CurrYear - StartYear] += 1;
                    }
                    RSApop.NewLLETZ[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                        1;
                    if (HIVstage == 0) {
                      if (SII < 0.752) {
                        if (SIII < 0.15) {
                          for (xx = 0; xx < 13; xx++) {
                            if (HPVstage[xx] == 2 || HPVstage[xx] == 3 ||
                                HPVstage[xx] == 4) {
                              HPVstageE[xx] = 1;
                              if (HPVstage[xx] == 4) {
                                WeibullCIN3[xx] = 0;
                                TimeinCIN3[xx] = 0;
                              }
                            }
                          }
                        } else {
                          for (xx = 0; xx < 13; xx++) {
                            if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                                HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                              HPVstageE[xx] = 0;
                              if (HPVstage[xx] == 4) {
                                WeibullCIN3[xx] = 0;
                                TimeinCIN3[xx] = 0;
                              }
                            }
                          }
                        }
                      } else {
                        if (SIV < 0.5) {
                          for (xx = 0; xx < 13; xx++) {
                            if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                                HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                              HPVstageE[xx] = 1;
                              if (HPVstage[xx] == 4) {
                                WeibullCIN3[xx] = 0;
                                TimeinCIN3[xx] = 0;
                              }
                            }
                          }
                        }
                      }
                    } else {
                      if (SII < 0.4) {
                        if (SIII < 0.15) {
                          for (xx = 0; xx < 13; xx++) {
                            if (HPVstage[xx] == 2 || HPVstage[xx] == 3 ||
                                HPVstage[xx] == 4) {
                              HPVstageE[xx] = 1;
                              if (HPVstage[xx] == 4) {
                                WeibullCIN3[xx] = 0;
                                TimeinCIN3[xx] = 0;
                              }
                            }
                          }
                        } else {
                          for (xx = 0; xx < 13; xx++) {
                            if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                                HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                              HPVstageE[xx] = 0;
                              if (HPVstage[xx] == 4) {
                                WeibullCIN3[xx] = 0;
                                TimeinCIN3[xx] = 0;
                              }
                            }
                          }
                        }
                      } else {
                        if (SIV < 0.5) {
                          for (xx = 0; xx < 13; xx++) {
                            if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                                HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                              HPVstageE[xx] = 1;
                              if (HPVstage[xx] == 4) {
                                WeibullCIN3[xx] = 0;
                                TimeinCIN3[xx] = 0;
                              }
                            }
                          }
                        }
                      }
                    }

                    if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                      if (HIVstage == 5) {
                        timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
                      } else {
                        timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
                      }
                      if (timetoScreen == 0) {
                        timetoScreen = 1;
                      }
                    }
                    if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
                      timetoScreen = 48;
                    }
                    repeat = 1;
                  }
                } else {
                  // Supposed to repeat smear in 3 months
                  // derived from NHLS data - Weibull distr with scale 31.2
                  // 3-months and shape 0.57
                  if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                    timetoScreen = 31.2 * pow(-log(tts), (1.0 / 0.57)) * 12;
                    if (timetoScreen == 0) {
                      timetoScreen = 1;
                    }
                  }
                  if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
                    timetoScreen = 12;
                  }

                  repeat = 1;
                }
              }
            }

            // if other type, ask to come back in year
            else if (AnyHPV(HPVstage, {5, 7, 9, 11, 12},
                            {1, 2, 3, 4, 5, 8, 9, 10})) {
              if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                if (HIVstage == 5) {
                  timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
                } else {
                  timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
                }
                if (timetoScreen == 0) {
                  timetoScreen = 1;
                }
              }
              if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
                timetoScreen = 48;
              }
              repeat = 1;
              HPVrepeat = 1;
            }
          } else {
            HPVrepeat = 0;
            RSApop.GetReferred[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
            if (HIVstage == 5) {
              if (ttC < AttendColposcopy[2][CurrYear - StartYear]) {
                timetoCol = 24;
                timetoScreen = 0;
              } else {
                timetoCol = 0;
                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  // timetoScreen = 5.3 * pow(-log(tts),(1.0/0.78)) * 48;
                  timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                  repeat = 1;
                }
              }
            }
            if (HIVstage > 0 && HIVstage != 5) {
              if (ttC < AttendColposcopy[1][CurrYear - StartYear]) {
                timetoCol = 24;
                timetoScreen = 0;
              } else {
                timetoCol = 0;
                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48;
                  timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                  repeat = 1;
                }
              }
            }
            if (HIVstage == 0) {
              if (ttC < AttendColposcopy[0][CurrYear - StartYear]) {
                timetoCol = 24;
                timetoScreen = 0;
              } else {
                timetoCol = 0;
                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48;
                  timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                  repeat = 1;
                }
              }
            }
          }
        } else {
          if (PapTRIAGE == 0) {
            // Refer to colposcopy
            RSApop.GetReferred[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
            if (HIVstage == 5) {
              if (ttC < AttendColposcopy[2][CurrYear - StartYear]) {
                timetoCol = 24;
                timetoScreen = 0;
              } else {
                timetoCol = 0;
                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  // timetoScreen = 5.3 * pow(-log(tts),(1.0/0.78)) * 48;
                  timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                  repeat = 1;
                }
              }
            }
            if (HIVstage > 0 && HIVstage != 5) {
              if (ttC < AttendColposcopy[1][CurrYear - StartYear]) {
                timetoCol = 24;
                timetoScreen = 0;
              } else {
                timetoCol = 0;
                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48;
                  timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                  repeat = 1;
                }
              }
            }
            if (HIVstage == 0) {
              if (ttC < AttendColposcopy[0][CurrYear - StartYear]) {
                timetoCol = 24;
                timetoScreen = 0;
              } else {
                timetoCol = 0;
                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48;
                  timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                  repeat = 1;
                }
              }
            }
          } else {
            RSApop.NewScreen[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
            if (ade < PapAdequacy[CurrYear - StartYear]) {
              if (HIVstage == 0) {
                if (TrueStage < 2) {
                  if (res < 0.75) {
                    ScreenResult = 0;
                  } else {
                    ScreenResult = 2;
                  }
                }
                if (TrueStage >= 2 && TrueStage < 5) {
                  if (res < 0.72) {
                    ScreenResult = 2;
                  } else {
                    ScreenResult = 0;
                  }
                }
              }
              if (HIVstage > 0) {
                if (TrueStage < 2) {
                  if (res < 0.44) {
                    ScreenResult = 0;
                  } else {
                    ScreenResult = 2;
                  }
                }
                if (TrueStage >= 2 && TrueStage < 5) {
                  if (res < 0.92) {
                    ScreenResult = 2;
                  } else {
                    ScreenResult = 0;
                  }
                }
              }

			  if (TrueStage == 3 && ScreenResult == 2 &&
				  CCd < 0.35)
			  { // CCd is sensitivity of Pap to pick up cancer
				  DiagnosedCC = 1;
				  RSApop.NewDiagCancer[zz * 18 + AgeGroup][CurrYear - StartYear] +=
					  1;
				  if (HIVstage == 5)
				  {
					  RSApop.NewDiagCancerART[CurrYear - StartYear] += 1;
				  }
				  for (int xx = 0; xx < 13; xx++)
				  {
					  if (HPVstage[xx] == 5)
					  {
						  HPVstageE[xx] = 11;
						  RSApop.StageDiag[0][CurrYear - StartYear] += 1;
						  RSApop.StageIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
						  if (SId < 0.192)
						  {
							  StageIdeath = static_cast<int>(48.0 * 3.08 *
															 pow(-log(SI), 1.0 / 1.23));
						  }
						  else
						  {
							  StageIrecover = 8;
						  }
					  }
					  if (HPVstage[xx] == 8)
					  {
						  HPVstageE[xx] = 12;
						  RSApop.StageDiag[1][CurrYear - StartYear] += 1;
						  RSApop.StageIIdiag[zz * 18 + AgeGroup]
											[CurrYear - StartYear] += 1;
						  if (SIId < 0.466)
						  {
							  StageIIdeath = static_cast<int>(
								  48.0 * 2.39 * pow(-log(SII), 1.0 / 1.17));
						  }
						  else
						  {
							  StageIIrecover = 24;
						  }
					  }
					  if (HPVstage[xx] == 9)
					  {
						  HPVstageE[xx] = 13;
						  RSApop.StageDiag[2][CurrYear - StartYear] += 1;
						  RSApop.StageIIIdiag[zz * 18 + AgeGroup]
											 [CurrYear - StartYear] += 1;
						  if (SIIId < 0.715)
						  {
							  StageIIIdeath = static_cast<int>(
								  48.0 * 1.18 * pow(-log(SIII), 1.0 / 0.91));
						  }
						  else
						  {
							  StageIIIrecover = 24;
						  }
					  }
					  if (HPVstage[xx] == 10)
					  {
						  HPVstageE[xx] = 14;
						  RSApop.StageDiag[3][CurrYear - StartYear] += 1;
						  RSApop.StageIVdiag[zz * 18 + AgeGroup]
											[CurrYear - StartYear] += 1;
						  if (SIVd < 0.931)
						  {
							  StageIVdeath = static_cast<int>(
								  48.0 * 0.46 * pow(-log(SIV), 1.0 / 0.9));
						  }
						  else
						  {
							  StageIVrecover = 24;
						  }
					  }
				  }
			  }

			  // Normal screen:
              if (ScreenResult == 0) {
                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  // if(HIVstage==5) {timetoScreen = 5.3 *
                  // pow(-log(tts),(1.0/0.78)) * 48;} else if(AgeExact<50) {
                  // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48; }
                  if (HIVstage == 5) {
                    timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
                  } else if (AgeExact < 50) {
                    timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                  } else {
                    timetoScreen = 200 * 48;
                  }
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                }
                if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
                  if (HIVstage == 5) {
                    timetoScreen = 3 * 48;
                  } else if (AgeExact < 50) {
                    timetoScreen = 10 * 48;
                  } else {
                    timetoScreen = 200 * 48;
                  }
                }
                repeat = 0;
              }
              // HSIL/CC screen:
              if (ScreenResult == 2 && TrueStage < 3 && SI < 0.9) {
                GetTreatment.out[AgeGroup][CurrYear - StartYear] += 1;
                if (TrueStage == 0) {
                  RSApop.NewUnnecessary[zz * 18 + AgeGroup]
                                       [CurrYear - StartYear] += 1;
                }
                RSApop.NewLLETZ[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
                if (HIVstage == 0) {
                  if (SII < 0.752) {
                    if (SIII < 0.15) {
                      for (xx = 0; xx < 13; xx++) {
                        if (HPVstage[xx] == 2 || HPVstage[xx] == 3 ||
                            HPVstage[xx] == 4) {
                          HPVstageE[xx] = 1;
                          if (HPVstage[xx] == 4) {
                            WeibullCIN3[xx] = 0;
                            TimeinCIN3[xx] = 0;
                          }
                        }
                      }
                    } else {
                      for (xx = 0; xx < 13; xx++) {
                        if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                            HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                          HPVstageE[xx] = 0;
                          if (HPVstage[xx] == 4) {
                            WeibullCIN3[xx] = 0;
                            TimeinCIN3[xx] = 0;
                          }
                        }
                      }
                    }
                  } else {
                    if (SIV < 0.5) {
                      for (xx = 0; xx < 13; xx++) {
                        if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                            HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                          HPVstageE[xx] = 1;
                          if (HPVstage[xx] == 4) {
                            WeibullCIN3[xx] = 0;
                            TimeinCIN3[xx] = 0;
                          }
                        }
                      }
                    }
                  }
                } else {
                  if (SII < 0.4) {
                    if (SIII < 0.15) {
                      for (xx = 0; xx < 13; xx++) {
                        if (HPVstage[xx] == 2 || HPVstage[xx] == 3 ||
                            HPVstage[xx] == 4) {
                          HPVstageE[xx] = 1;
                          if (HPVstage[xx] == 4) {
                            WeibullCIN3[xx] = 0;
                            TimeinCIN3[xx] = 0;
                          }
                        }
                      }
                    } else {
                      for (xx = 0; xx < 13; xx++) {
                        if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                            HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                          HPVstageE[xx] = 0;
                          if (HPVstage[xx] == 4) {
                            WeibullCIN3[xx] = 0;
                            TimeinCIN3[xx] = 0;
                          }
                        }
                      }
                    }
                  } else {
                    if (SIV < 0.5) {
                      for (xx = 0; xx < 13; xx++) {
                        if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                            HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                          HPVstageE[xx] = 1;
                          if (HPVstage[xx] == 4) {
                            WeibullCIN3[xx] = 0;
                            TimeinCIN3[xx] = 0;
                          }
                        }
                      }
                    }
                  }
                }

                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  if (HIVstage == 5) {
                    timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
                  } else {
                    timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
                  }
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                }
                if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
                  timetoScreen = 48;
                }
                repeat = 1;
              }
            } else {
              // Supposed to repeat smear in 3 months
              // derived from NHLS data - Weibull distr with scale 31.2 3-months
              // and shape 0.57
              if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                timetoScreen = 31.2 * pow(-log(tts), (1.0 / 0.57)) * 12;
                if (timetoScreen == 0) {
                  timetoScreen = 1;
                }
              }
              if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
                timetoScreen = 12;
              }

              repeat = 1;
            }
          }
        }
      } else { // if don't come for results, assume normal screening resumes
        if ((WHOScreening == 0 && PerfectSchedule == 0) ||
            CurrYear < ImplementYR) {
          // if(HIVstage==5) {timetoScreen = 5.3 * pow(-log(tts),(1.0/0.78)) *
          // 48;} else if(AgeExact<50) { timetoScreen = 15.0 *
          // pow(-log(tts),(1.0/0.83)) * 48; }
          if (HIVstage == 5) {
            timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
          } else if (AgeExact < 50) {
            timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
          } else {
            timetoScreen = 200 * 48;
          }
          if (timetoScreen == 0) {
            timetoScreen = 1;
          }
        }
        if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
          if (HIVstage == 5) {
            timetoScreen = 3 * 48;
          } else if (AgeExact < 50) {
            timetoScreen = 10 * 48;
          } else {
            timetoScreen = 200 * 48;
          }
        }
        repeat = 0;
      }
    } else {
      // Supposed to repeat smear in 3 months
      // derived from NHLS data - Weibull distr with scale 31.2  3-months and
      // shape 0.57
      if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
        timetoScreen = 31.2 * pow(-log(tts), (1.0 / 0.57)) * 12;
        if (timetoScreen == 0) {
          timetoScreen = 1;
        }
      }
      if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
        timetoScreen = 12;
      }
      repeat = 1;
    }
  }
  void Indiv::HPV_ThermalScreenAlgorithm(
      int ID, double rea, double ade, double tts, double res, double ttC,
      double CCd, double SI, double SII, double SIII, double SIV, double SId,
      double SIId, double SIIId, double SIVd) {
    int xx, yy, zz;
    int SimCount2 = (CurrSim - 1) / IterationsPerPC;
    if (HIVstage == 0) {
      zz = 0;
    } else if (HIVstage == 5) {
      zz = 1;
    } //||HIVstage==6
    else {
      zz = 2;
    }

    RSApop.NewHPVScreen[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;

    if (ade < 0.95) {
      if (HPVstatus == 0) { // if HPV negative, resume normal screening
                            // intervals
        if ((WHOScreening == 0 && PerfectSchedule == 0) ||
            CurrYear < ImplementYR) {
          if (HIVstage == 5) {
            timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
          } else if (AgeExact < 50) {
            timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
          } else {
            timetoScreen = 200 * 48;
          }
          if (timetoScreen == 0) {
            timetoScreen = 1;
          }
        }
        if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
          if (HIVstage == 5) {
            timetoScreen = 3 * 48;
          } else if (AgeExact < 50) {
            timetoScreen = 10 * 48;
          } else {
            timetoScreen = 200 * 48;
          }
        }
        repeat = 0;
        HPVrepeat = 0;
      } else if (rea < 0.9) { // 90% come back for results
        if (HPVGenotyping == 1) {
          if (HPVrepeat == 0) {
            // if test positive with P1&2: HPV16/18/45 (but not cancer)
            if (AnyHPV(HPVstage, hpv161845, {1, 2, 3, 4}) && TrueStage < 3) {
              RSApop.NewVAT[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
              if (CCd < 0.913) { // 91.3% suitable for ablation after VAT
                GetTreatment.out[AgeGroup][CurrYear - StartYear] += 1;
                if (TrueStage == 0) {
                  RSApop.NewUnnecessary[zz * 18 + AgeGroup]
                                       [CurrYear - StartYear] += 1;
                }
                RSApop.NewThermal[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                    1;

                if (HIVstage == 0) {
                  if (res < 0.652) { 
                    for (xx = 0; xx < 13; xx++) {
                      if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                          HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                        HPVstageE[xx] = 0;
                        if (HPVstage[xx] == 4) {
                          WeibullCIN3[xx] = 0;
                          TimeinCIN3[xx] = 0;
                        }
                      }
                    }
                  } else {
                    for (xx = 0; xx < 13; xx++) {
                      if (HPVstage[xx] == 2) {
                        HPVstageE[xx] = 1;
                      }
                      if (HPVstage[xx] == 3) {
                        HPVstageE[xx] = 2;
                      }
                      if (HPVstage[xx] == 4) {
                        HPVstageE[xx] = 3;
                        WeibullCIN3[xx] = 0;
                        TimeinCIN3[xx] = 0;
                      }
                    }
                  }
                } else {
                  if (res < 0.385) {// 0.385
                    for (xx = 0; xx < 13; xx++) {
                      if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                          HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                        HPVstageE[xx] = 0;
                        if (HPVstage[xx] == 4) {
                          WeibullCIN3[xx] = 0;
                          TimeinCIN3[xx] = 0;
                        }
                      }
                    }
                  } else {
                    for (xx = 0; xx < 13; xx++) {
                      if (HPVstage[xx] == 2) {
                        HPVstageE[xx] = 1;
                      }
                      if (HPVstage[xx] == 3) {
                        HPVstageE[xx] = 2;
                      }
                      if (HPVstage[xx] == 4) {
                        HPVstageE[xx] = 3;
                        WeibullCIN3[xx] = 0;
                        TimeinCIN3[xx] = 0;
                      }
                    }
                  }
                }
                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  if (HIVstage == 5) {
                    timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
                  } else {
                    timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
                  }
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                }
                if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
                  timetoScreen = 48;
                }
                repeat = 1;
                HPVrepeat = 1;
              } else { // rest get referred to colposcopy
                RSApop.GetReferred[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                    1;
                if (HIVstage == 5) {
                  if (ttC < AttendColposcopy[2][CurrYear - StartYear]) {
                    timetoCol = 24;
                    timetoScreen = 0;
                  } else {
                    timetoCol = 0;
                    if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                      // timetoScreen = 5.3 * pow(-log(tts),(1.0/0.78)) * 48;
                      timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
                      if (timetoScreen == 0) {
                        timetoScreen = 1;
                      }
                      repeat = 1;
                    }
                  }
                }
                if (HIVstage > 0 && HIVstage != 5) {
                  if (ttC < AttendColposcopy[1][CurrYear - StartYear]) {
                    timetoCol = 24;
                    timetoScreen = 0;
                  } else {
                    timetoCol = 0;
                    if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                      // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48;
                      timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                      if (timetoScreen == 0) {
                        timetoScreen = 1;
                      }
                      repeat = 1;
                    }
                  }
                }
                if (HIVstage == 0) {
                  if (ttC < AttendColposcopy[0][CurrYear - StartYear]) {
                    timetoCol = 24;
                    timetoScreen = 0;
                  } else {
                    timetoCol = 0;
                    if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                      // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48;
                      timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                      if (timetoScreen == 0) {
                        timetoScreen = 1;
                      }
                      repeat = 1;
                    }
                  }
                }
              }
            } else if (AnyHPV(HPVstage, hpv161845, cc_un34)) {
              DiagnosedCC = 1;
              RSApop.NewDiagCancer[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                  1;
              if (AnyHPV(HPVstage, hpv1618, cc_un34)) {
                RSApop.NewDiagCancer1618[AgeGroup][CurrYear - StartYear] += 1;
              }
              if (HIVstage == 5) {
                RSApop.NewDiagCancerART[CurrYear - StartYear] += 1;
              }
              for (xx = 0; xx < 13; xx++) {
                if (HPVstage[xx] == 5) {
                  HPVstageE[xx] = 11;
                  RSApop.StageDiag[0][CurrYear - StartYear] += 1;
                  RSApop.StageIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                      1;
                  if (SId < 0.192) {
                    StageIdeath = static_cast<int>(48.0 * 3.08 *
                                                   pow(-log(SI), 1.0 / 1.23));
                  } else {
                    StageIrecover = 8;
                  }
                }
                if (HPVstage[xx] == 8) {
                  HPVstageE[xx] = 12;
                  RSApop.StageDiag[1][CurrYear - StartYear] += 1;
                  RSApop
                      .StageIIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                      1;
                  if (SIId < 0.466) {
                    StageIIdeath = static_cast<int>(48.0 * 2.39 *
                                                    pow(-log(SII), 1.0 / 1.17));
                  } else {
                    StageIIrecover = 24;
                  }
                }
                if (HPVstage[xx] == 9) {
                  HPVstageE[xx] = 13;
                  RSApop.StageDiag[2][CurrYear - StartYear] += 1;
                  RSApop
                      .StageIIIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                      1;
                  if (SIIId < 0.715) {
                    StageIIIdeath = static_cast<int>(
                        48.0 * 1.18 * pow(-log(SIII), 1.0 / 0.91));
                  } else {
                    StageIIIrecover = 24;
                  }
                }
                if (HPVstage[xx] == 10) {
                  HPVstageE[xx] = 14;
                  RSApop.StageDiag[3][CurrYear - StartYear] += 1;
                  RSApop
                      .StageIVdiag[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                      1;
                  if (SIVd < 0.931) {
                    StageIVdeath = static_cast<int>(48.0 * 0.46 *
                                                    pow(-log(SIV), 1.0 / 0.9));
                  } else {
                    StageIVrecover = 24;
                  }
                }
              }
            }

            // if test positive with P3: HPV31,33, 35, 52, 58
            else if (AnyHPV(HPVstage, {2, 3, 4, 8, 10},
                            {1, 2, 3, 4, 5, 8, 9, 10})) {
              if (Portal3 == 0) { // Follow-up after 1 year
                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  if (HIVstage == 5) {
                    timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
                  } else {
                    timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
                  }
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                }
                if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
                  timetoScreen = 48;
                }
                repeat = 1;
                HPVrepeat = 1;
              }
              if (Portal3 == 1) { // screen-and-treat
                RSApop.NewVAT[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
                if (TrueStage < 3) {
                  if (CCd < 0.913) { // 91.3% suitable for ablation after VAT
                    GetTreatment.out[AgeGroup][CurrYear - StartYear] += 1;
                    if (TrueStage == 0) {
                      RSApop.NewUnnecessary[zz * 18 + AgeGroup]
                                           [CurrYear - StartYear] += 1;
                    }
                    RSApop.NewThermal[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
                    if (HIVstage == 0) {
                      if (res < 0.689) {
                        for (xx = 0; xx < 13; xx++) {
                          if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                              HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                            HPVstageE[xx] = 0;
                            if (HPVstage[xx] == 4) {
                              WeibullCIN3[xx] = 0;
                              TimeinCIN3[xx] = 0;
                            }
                          }
                        }
                      } else {
                        for (xx = 0; xx < 13; xx++) {
                          if (HPVstage[xx] == 2) {
                            HPVstageE[xx] = 1;
                          }
                          if (HPVstage[xx] == 3) {
                            HPVstageE[xx] = 2;
                          }
                          if (HPVstage[xx] == 4) {
                            HPVstageE[xx] = 3;
                            WeibullCIN3[xx] = 0;
                            TimeinCIN3[xx] = 0;
                          }
                        }
                      }
                    } else {
                      if (res < 0.585) {
                        for (xx = 0; xx < 13; xx++) {
                          if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                              HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                            HPVstageE[xx] = 0;
                            if (HPVstage[xx] == 4) {
                              WeibullCIN3[xx] = 0;
                              TimeinCIN3[xx] = 0;
                            }
                          }
                        }
                      } else {
                        for (xx = 0; xx < 13; xx++) {
                          if (HPVstage[xx] == 2) {
                            HPVstageE[xx] = 1;
                          }
                          if (HPVstage[xx] == 3) {
                            HPVstageE[xx] = 2;
                          }
                          if (HPVstage[xx] == 4) {
                            HPVstageE[xx] = 3;
                            WeibullCIN3[xx] = 0;
                            TimeinCIN3[xx] = 0;
                          }
                        }
                      }
                    }
                    if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                      if (HIVstage == 5) {
                        timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
                      } else {
                        timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
                      }
                      if (timetoScreen == 0) {
                        timetoScreen = 1;
                      }
                    }
                    if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
                      timetoScreen = 48;
                    }
                    repeat = 1;
                    HPVrepeat = 1;
                  } else { // rest get referred to colposcopy
                    RSApop.GetReferred[zz * 18 + AgeGroup]
                                      [CurrYear - StartYear] += 1;
                    if (HIVstage == 5) {
                      if (ttC < AttendColposcopy[2][CurrYear - StartYear]) {
                        timetoCol = 24;
                        timetoScreen = 0;
                      } else {
                        timetoCol = 0;
                        if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                          // timetoScreen = 5.3 * pow(-log(tts),(1.0/0.78)) *
                          // 48;
                          timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
                          if (timetoScreen == 0) {
                            timetoScreen = 1;
                          }
                          repeat = 1;
                        }
                      }
                    }
                    if (HIVstage > 0 && HIVstage != 5) {
                      if (ttC < AttendColposcopy[1][CurrYear - StartYear]) {
                        timetoCol = 24;
                        timetoScreen = 0;
                      } else {
                        timetoCol = 0;
                        if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                          // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) *
                          // 48;
                          timetoScreen =
                              15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                          if (timetoScreen == 0) {
                            timetoScreen = 1;
                          }
                          repeat = 1;
                        }
                      }
                    }
                    if (HIVstage == 0) {
                      if (ttC < AttendColposcopy[0][CurrYear - StartYear]) {
                        timetoCol = 24;
                        timetoScreen = 0;
                      } else {
                        timetoCol = 0;
                        if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                          // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) *
                          // 48;
                          timetoScreen =
                              15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                          if (timetoScreen == 0) {
                            timetoScreen = 1;
                          }
                          repeat = 1;
                        }
                      }
                    }
                  }
                } else if (TrueStage == 3) {
                  DiagnosedCC = 1;
                  RSApop.NewDiagCancer[zz * 18 + AgeGroup]
                                      [CurrYear - StartYear] += 1;
                  if (HIVstage == 5) {
                    RSApop.NewDiagCancerART[CurrYear - StartYear] += 1;
                  }
                  for (xx = 0; xx < 13; xx++) {
                    if (HPVstage[xx] == 5) {
                      HPVstageE[xx] = 11;
                      RSApop.StageDiag[0][CurrYear - StartYear] += 1;
                      RSApop.StageIdiag[zz * 18 + AgeGroup]
                                       [CurrYear - StartYear] += 1;
                      if (SId < 0.192) {
                        StageIdeath = static_cast<int>(
                            48.0 * 3.08 * pow(-log(SI), 1.0 / 1.23));
                      } else {
                        StageIrecover = 8;
                      }
                    }
                    if (HPVstage[xx] == 8) {
                      HPVstageE[xx] = 12;
                      RSApop.StageDiag[1][CurrYear - StartYear] += 1;
                      RSApop.StageIIdiag[zz * 18 + AgeGroup]
                                        [CurrYear - StartYear] += 1;
                      if (SIId < 0.466) {
                        StageIIdeath = static_cast<int>(
                            48.0 * 2.39 * pow(-log(SII), 1.0 / 1.17));
                      } else {
                        StageIIrecover = 24;
                      }
                    }
                    if (HPVstage[xx] == 9) {
                      HPVstageE[xx] = 13;
                      RSApop.StageDiag[2][CurrYear - StartYear] += 1;
                      RSApop.StageIIIdiag[zz * 18 + AgeGroup]
                                         [CurrYear - StartYear] += 1;
                      if (SIIId < 0.715) {
                        StageIIIdeath = static_cast<int>(
                            48.0 * 1.18 * pow(-log(SIII), 1.0 / 0.91));
                      } else {
                        StageIIIrecover = 24;
                      }
                    }
                    if (HPVstage[xx] == 10) {
                      HPVstageE[xx] = 14;
                      RSApop.StageDiag[3][CurrYear - StartYear] += 1;
                      RSApop.StageIVdiag[zz * 18 + AgeGroup]
                                        [CurrYear - StartYear] += 1;
                      if (SIVd < 0.931) {
                        StageIVdeath = static_cast<int>(
                            48.0 * 0.46 * pow(-log(SIV), 1.0 / 0.9));
                      } else {
                        StageIVrecover = 24;
                      }
                    }
                  }
                }
              }
            }
            // if test positive with P4&5: HPV39, 51, 56, 59, 68
            else if (AnyHPV(HPVstage, {5, 7, 9, 11, 12},
                            {1, 2, 3, 4, 5, 8, 9, 10})) {
              if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                if (HIVstage == 5) {
                  timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
                } else {
                  timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
                }
                if (timetoScreen == 0) {
                  timetoScreen = 1;
                }
              }
              if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
                timetoScreen = 48;
              }
              repeat = 1;
              HPVrepeat = 1;
            }
          } else {
            HPVrepeat = 0;
            RSApop.GetReferred[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
            if (HIVstage == 5) {
              if (ttC < AttendColposcopy[2][CurrYear - StartYear]) {
                timetoCol = 24;
                timetoScreen = 0;
              } else {
                timetoCol = 0;
                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  // timetoScreen = 5.3 * pow(-log(tts),(1.0/0.78)) * 48;
                  timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                  repeat = 1;
                }
              }
            }
            if (HIVstage > 0 && HIVstage != 5) {
              if (ttC < AttendColposcopy[1][CurrYear - StartYear]) {
                timetoCol = 24;
                timetoScreen = 0;
              } else {
                timetoCol = 0;
                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48;
                  timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                  repeat = 1;
                }
              }
            }
            if (HIVstage == 0) {
              if (ttC < AttendColposcopy[0][CurrYear - StartYear]) {
                timetoCol = 24;
                timetoScreen = 0;
              } else {
                timetoCol = 0;
                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48;
                  timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                  repeat = 1;
                }
              }
            }
          }
        } else {
          if (HPVrepeat == 0) {
            if (TrueStage < 3) {
              RSApop.NewVAT[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
              if (CCd < 0.913) { // 91.3% suitable for ablation after VAT
                GetTreatment.out[AgeGroup][CurrYear - StartYear] += 1;
                if (TrueStage == 0) {
                  RSApop.NewUnnecessary[zz * 18 + AgeGroup]
                                       [CurrYear - StartYear] += 1;
                }
                RSApop.NewThermal[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                    1;
                if (HIVstage == 0) {
                  if (res < 0.689) {
                    for (xx = 0; xx < 13; xx++) {
                      if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                          HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                        HPVstageE[xx] = 0;
                        if (HPVstage[xx] == 4) {
                          WeibullCIN3[xx] = 0;
                          TimeinCIN3[xx] = 0;
                        }
                      }
                    }
                  } else {
                    for (xx = 0; xx < 13; xx++) {
                      if (HPVstage[xx] == 2) {
                        HPVstageE[xx] = 1;
                      }
                      if (HPVstage[xx] == 3) {
                        HPVstageE[xx] = 2;
                      }
                      if (HPVstage[xx] == 4) {
                        HPVstageE[xx] = 3;
                        WeibullCIN3[xx] = 0;
                        TimeinCIN3[xx] = 0;
                      }
                    }
                  }
                } else {
                  if (res < 0.585) {
                    for (xx = 0; xx < 13; xx++) {
                      if (HPVstage[xx] == 1 || HPVstage[xx] == 2 ||
                          HPVstage[xx] == 3 || HPVstage[xx] == 4) {
                        HPVstageE[xx] = 0;
                        if (HPVstage[xx] == 4) {
                          WeibullCIN3[xx] = 0;
                          TimeinCIN3[xx] = 0;
                        }
                      }
                    }
                  } else {
                    for (xx = 0; xx < 13; xx++) {
                      if (HPVstage[xx] == 2) {
                        HPVstageE[xx] = 1;
                      }
                      if (HPVstage[xx] == 3) {
                        HPVstageE[xx] = 2;
                      }
                      if (HPVstage[xx] == 4) {
                        HPVstageE[xx] = 3;
                        WeibullCIN3[xx] = 0;
                        TimeinCIN3[xx] = 0;
                      }
                    }
                  }
                }
                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  if (HIVstage == 5) {
                    timetoScreen = 4.5 * pow(-log(tts), (1.0 / 0.71)) * 48;
                  } else {
                    timetoScreen = 9.0 * pow(-log(tts), (1.0 / 0.56)) * 48;
                  }
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                }
                if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
                  timetoScreen = 48;
                }
                repeat = 1;
                HPVrepeat = 1;
              } else { // rest get referred to colposcopy
                RSApop.GetReferred[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                    1;
                if (HIVstage == 5) {
                  if (ttC < AttendColposcopy[2][CurrYear - StartYear]) {
                    timetoCol = 24;
                    timetoScreen = 0;
                  } else {
                    timetoCol = 0;
                    if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                      // timetoScreen = 5.3 * pow(-log(tts),(1.0/0.78)) * 48;
                      timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
                      if (timetoScreen == 0) {
                        timetoScreen = 1;
                      }
                      repeat = 1;
                    }
                  }
                }
                if (HIVstage > 0 && HIVstage != 5) {
                  if (ttC < AttendColposcopy[1][CurrYear - StartYear]) {
                    timetoCol = 24;
                    timetoScreen = 0;
                  } else {
                    timetoCol = 0;
                    if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                      // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48;
                      timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                      if (timetoScreen == 0) {
                        timetoScreen = 1;
                      }
                      repeat = 1;
                    }
                  }
                }
                if (HIVstage == 0) {
                  if (ttC < AttendColposcopy[0][CurrYear - StartYear]) {
                    timetoCol = 24;
                    timetoScreen = 0;
                  } else {
                    timetoCol = 0;
                    if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                      // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48;
                      timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                      if (timetoScreen == 0) {
                        timetoScreen = 1;
                      }
                      repeat = 1;
                    }
                  }
                }
              }
            } else if (TrueStage == 3) {
              DiagnosedCC = 1;
              RSApop.NewDiagCancer[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                  1;
              if (AnyHPV(HPVstage, hpv1618, cc_un34)) {
                RSApop.NewDiagCancer1618[AgeGroup][CurrYear - StartYear] += 1;
              }
              if (HIVstage == 5) {
                RSApop.NewDiagCancerART[CurrYear - StartYear] += 1;
              }
              for (xx = 0; xx < 13; xx++) {
                if (HPVstage[xx] == 5) {
                  HPVstageE[xx] = 11;
                  RSApop.StageDiag[0][CurrYear - StartYear] += 1;
                  RSApop.StageIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                      1;
                  if (SId < 0.192) {
                    StageIdeath = static_cast<int>(48.0 * 3.08 *
                                                   pow(-log(SI), 1.0 / 1.23));
                  } else {
                    StageIrecover = 8;
                  }
                }
                if (HPVstage[xx] == 8) {
                  HPVstageE[xx] = 12;
                  RSApop.StageDiag[1][CurrYear - StartYear] += 1;
                  RSApop
                      .StageIIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                      1;
                  if (SIId < 0.466) {
                    StageIIdeath = static_cast<int>(48.0 * 2.39 *
                                                    pow(-log(SII), 1.0 / 1.17));
                  } else {
                    StageIIrecover = 24;
                  }
                }
                if (HPVstage[xx] == 9) {
                  HPVstageE[xx] = 13;
                  RSApop.StageDiag[2][CurrYear - StartYear] += 1;
                  RSApop
                      .StageIIIdiag[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                      1;
                  if (SIIId < 0.715) {
                    StageIIIdeath = static_cast<int>(
                        48.0 * 1.18 * pow(-log(SIII), 1.0 / 0.91));
                  } else {
                    StageIIIrecover = 24;
                  }
                }
                if (HPVstage[xx] == 10) {
                  HPVstageE[xx] = 14;
                  RSApop.StageDiag[3][CurrYear - StartYear] += 1;
                  RSApop
                      .StageIVdiag[zz * 18 + AgeGroup][CurrYear - StartYear] +=
                      1;
                  if (SIVd < 0.931) {
                    StageIVdeath = static_cast<int>(48.0 * 0.46 *
                                                    pow(-log(SIV), 1.0 / 0.9));
                  } else {
                    StageIVrecover = 24;
                  }
                }
              }
            }
          } else {
            HPVrepeat = 0;
            RSApop.GetReferred[zz * 18 + AgeGroup][CurrYear - StartYear] += 1;
            if (HIVstage == 5) {
              if (ttC < AttendColposcopy[2][CurrYear - StartYear]) {
                timetoCol = 24;
                timetoScreen = 0;
              } else {
                timetoCol = 0;
                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  // timetoScreen = 5.3 * pow(-log(tts),(1.0/0.78)) * 48;
                  timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                  repeat = 1;
                }
              }
            }
            if (HIVstage > 0 && HIVstage != 5) {
              if (ttC < AttendColposcopy[1][CurrYear - StartYear]) {
                timetoCol = 24;
                timetoScreen = 0;
              } else {
                timetoCol = 0;
                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48;
                  timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                  repeat = 1;
                }
              }
            }
            if (HIVstage == 0) {
              if (ttC < AttendColposcopy[0][CurrYear - StartYear]) {
                timetoCol = 24;
                timetoScreen = 0;
              } else {
                timetoCol = 0;
                if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
                  // timetoScreen = 15.0 * pow(-log(tts),(1.0/0.83)) * 48;
                  timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
                  if (timetoScreen == 0) {
                    timetoScreen = 1;
                  }
                  repeat = 1;
                }
              }
            }
          }
        }
      } else { // if didn't come for results, resume normal screening schedule
        if ((WHOScreening == 0 && PerfectSchedule == 0) ||
            CurrYear < ImplementYR) {
          // if(HIVstage==5) {timetoScreen = 5.3 * pow(-log(tts),(1.0/0.78)) *
          // 48;} else if(AgeExact<50) { timetoScreen = 15.0 *
          // pow(-log(tts),(1.0/0.83)) * 48; }
          if (HIVstage == 5) {
            timetoScreen = 7.9 * pow(-log(tts), (1.0 / 1.0)) * 48;
          } else if (AgeExact < 50) {
            timetoScreen = 15.0 * pow(-log(tts), (1.0 / 1.0)) * 48;
          } else {
            timetoScreen = 200 * 48;
          }
          if (timetoScreen == 0) {
            timetoScreen = 1;
          }
        }
        if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
          if (HIVstage == 5) {
            timetoScreen = 3 * 48;
          } else if (AgeExact < 50) {
            timetoScreen = 10 * 48;
          } else {
            timetoScreen = 200 * 48;
          }
        }
        repeat = 0;
      }
    } else {
      // If inadequate sample, repeat HPV in 3 months
      // derived from NHLS data - Weibull distr with scale 31.2  3-months and
      // shape 0.57
      if (PerfectSchedule == 0 || CurrYear < ImplementYR) {
        timetoScreen = 31.2 * pow(-log(tts), (1.0 / 0.57)) * 12;
        if (timetoScreen == 0) {
          timetoScreen = 1;
        }
      }
      if (PerfectSchedule == 1 && CurrYear >= ImplementYR) {
        timetoScreen = 12;
      }
      repeat = 1;
    }
  }

  void Indiv::PerfectGetScreened(
      int ID, double rea, double scr, double ade, double tts, double res,
      double ttC, double CCd, double SI, double SII, double SIII, double SIV,
      double SId, double SIId, double SIIId, double SIVd, double acceptRand,
      double efficacyD1Rand, double efficacyD2Rand, double D2LossRand) {
    int iy;
    iy = CurrYear - StartYear;

    int xx, yy, zz;
    // Get age that matches coverage age
    yy = 0;
    if (AgeGroup == 6 || AgeGroup == 7) {
      yy = 1;
    } else if (AgeGroup == 8 || AgeGroup == 9) {
      yy = 2;
    } else if (AgeGroup >= 10) {
      yy = 3;
    }

    if (HIVstage == 5 || HIVstage == 6) {
      zz = 1;
    } else {
      zz = 0;
    }

    if ((timetoScreen > 0 && timetoScreen < 96) &&
        timetoCol == 0) { // give it two years for those who were already
                          // scheduled for a screen to be screened
      if (timePassed > timetoScreen) {
        timetoScreen = 0;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        } else if (HPVDNA == 1 && CurrYear >= ImplementYR &&
                   ((HIVstage >= 5 && AgeExact >= 25.0 && AgeExact < 60.0) ||
                    (HIVstage < 5 && AgeExact >= 30.0 && AgeExact < 60.0))) {
          if (HPVDNAThermal == 0) {
            HPVScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                               SIV, SId, SIId, SIIId, SIVd);
          } else {
            if (ThermalORPap < PropThermal) {
              HPV_ThermalScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI,
                                         SII, SIII, SIV, SId, SIId, SIIId,
                                         SIVd);
            } else {
              ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                              SIV, SId, SIId, SIIId, SIVd, acceptRand,
                              efficacyD1Rand, efficacyD2Rand, D2LossRand);
            }
          }
        }
        timePassed = 0;
      } else {
        timePassed += 1;
        ScreenCount += 1;
      }
    }

    else if (HIVstage < 5 &&
             scr < WHOcoverage[zz * 4 + yy][CurrYear - StartYear] /
                       (48.0 * 10.0)) { // otherwise, schedule screens

      if (AgeExact >= 30.0 && AgeExact < 40.0 && Scr30 == 0) {
        Scr30 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        } else if (HPVDNA == 1 && CurrYear >= ImplementYR) {
          if (HPVDNAThermal == 0) {
            HPVScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                               SIV, SId, SIId, SIIId, SIVd);
          } else {
            if (ThermalORPap < PropThermal) {
              HPV_ThermalScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI,
                                         SII, SIII, SIV, SId, SIId, SIIId,
                                         SIVd);
            } else {
              ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                              SIV, SId, SIId, SIIId, SIVd, acceptRand,
                              efficacyD1Rand, efficacyD2Rand, D2LossRand);
            }
          }
        }
      }
      if (AgeExact >= 40.0 && AgeExact < 50.0 && Scr40 == 0) {
        Scr40 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        } else if (HPVDNA == 1 && CurrYear >= ImplementYR) {
          if (HPVDNAThermal == 0) {
            HPVScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                               SIV, SId, SIId, SIIId, SIVd);
          } else {
            if (ThermalORPap < PropThermal) {
              HPV_ThermalScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI,
                                         SII, SIII, SIV, SId, SIId, SIIId,
                                         SIVd);
            } else {
              ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                              SIV, SId, SIId, SIIId, SIVd, acceptRand,
                              efficacyD1Rand, efficacyD2Rand, D2LossRand);
            }
          }
        }
      }
      if (AgeExact >= 50.0 && AgeExact < 60.0 && Scr50 == 0) {
        Scr50 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        } else if (HPVDNA == 1 && CurrYear >= ImplementYR) {
          if (HPVDNAThermal == 0) {
            HPVScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                               SIV, SId, SIId, SIIId, SIVd);
          } else {
            if (ThermalORPap < PropThermal) {
              HPV_ThermalScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI,
                                         SII, SIII, SIV, SId, SIId, SIIId,
                                         SIVd);
            } else {
              ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                              SIV, SId, SIId, SIIId, SIVd, acceptRand,
                              efficacyD1Rand, efficacyD2Rand, D2LossRand);
            }
          }
        }
      }
    } else if ((HIVstage == 5 || HIVstage == 6) &&
               scr < WHOcoverage[zz * 4 + yy][CurrYear - StartYear] /
                         (3.0 * 48.0)) {
      if (AgeExact >= 15.0 && AgeExact < 18.0 && Scr16 == 0) {
        Scr16 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        }
      }
      if (AgeExact >= 18.0 && AgeExact < 21.0 && Scr19 == 0) {
        Scr19 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        }
      }
      if (AgeExact >= 21.0 && AgeExact < 24.0 && Scr22 == 0) {
        Scr22 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        }
      }
      if (AgeExact >= 24.0 && AgeExact < 27.0 && Scr25 == 0) {
        Scr25 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR || AgeExact == 24.0) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        } else if (AgeExact >= 25.0 && HPVDNA == 1 && CurrYear >= ImplementYR) {
          if (HPVDNAThermal == 0) {
            HPVScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                               SIV, SId, SIId, SIIId, SIVd);
          } else {
            if (ThermalORPap < PropThermal) {
              HPV_ThermalScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI,
                                         SII, SIII, SIV, SId, SIId, SIIId,
                                         SIVd);
            } else {
              ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                              SIV, SId, SIId, SIIId, SIVd, acceptRand,
                              efficacyD1Rand, efficacyD2Rand, D2LossRand);
            }
          }
        }
      }
      if (AgeExact >= 27.0 && AgeExact < 30.0 && Scr28 == 0) {
        Scr28 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        } else if (HPVDNA == 1 && CurrYear >= ImplementYR) {
          if (HPVDNAThermal == 0) {
            HPVScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                               SIV, SId, SIId, SIIId, SIVd);
          } else {
            if (ThermalORPap < PropThermal) {
              HPV_ThermalScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI,
                                         SII, SIII, SIV, SId, SIId, SIIId,
                                         SIVd);
            } else {
              ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                              SIV, SId, SIId, SIIId, SIVd, acceptRand,
                              efficacyD1Rand, efficacyD2Rand, D2LossRand);
            }
          }
        }
      }
      if (AgeExact >= 30.0 && AgeExact < 33.0 && Scr31 == 0) {
        Scr31 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        } else if (HPVDNA == 1 && CurrYear >= ImplementYR) {
          if (HPVDNAThermal == 0) {
            HPVScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                               SIV, SId, SIId, SIIId, SIVd);
          } else {
            if (ThermalORPap < PropThermal) {
              HPV_ThermalScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI,
                                         SII, SIII, SIV, SId, SIId, SIIId,
                                         SIVd);
            } else {
              ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                              SIV, SId, SIId, SIIId, SIVd, acceptRand,
                              efficacyD1Rand, efficacyD2Rand, D2LossRand);
            }
          }
        }
      }
      if (AgeExact >= 33.0 && AgeExact < 36.0 && Scr34 == 0) {
        Scr34 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        } else if (HPVDNA == 1 && CurrYear >= ImplementYR) {
          if (HPVDNAThermal == 0) {
            HPVScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                               SIV, SId, SIId, SIIId, SIVd);
          } else {
            if (ThermalORPap < PropThermal) {
              HPV_ThermalScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI,
                                         SII, SIII, SIV, SId, SIId, SIIId,
                                         SIVd);
            } else {
              ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                              SIV, SId, SIId, SIIId, SIVd, acceptRand,
                              efficacyD1Rand, efficacyD2Rand, D2LossRand);
            }
          }
        }
      }
      if (AgeExact >= 36.0 && AgeExact < 39.0 && Scr37 == 0) {
        Scr37 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        } else if (HPVDNA == 1 && CurrYear >= ImplementYR) {
          if (HPVDNAThermal == 0) {
            HPVScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                               SIV, SId, SIId, SIIId, SIVd);
          } else {
            if (ThermalORPap < PropThermal) {
              HPV_ThermalScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI,
                                         SII, SIII, SIV, SId, SIId, SIIId,
                                         SIVd);
            } else {
              ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                              SIV, SId, SIId, SIIId, SIVd, acceptRand,
                              efficacyD1Rand, efficacyD2Rand, D2LossRand);
            }
          }
        }
      }
      if (AgeExact >= 39.0 && AgeExact < 42.0 && Scr40 == 0) {
        Scr40 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        } else if (HPVDNA == 1 && CurrYear >= ImplementYR) {
          if (HPVDNAThermal == 0) {
            HPVScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                               SIV, SId, SIId, SIIId, SIVd);
          } else {
            if (ThermalORPap < PropThermal) {
              HPV_ThermalScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI,
                                         SII, SIII, SIV, SId, SIId, SIIId,
                                         SIVd);
            } else {
              ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                              SIV, SId, SIId, SIIId, SIVd, acceptRand,
                              efficacyD1Rand, efficacyD2Rand, D2LossRand);
            }
          }
        }
      }
      if (AgeExact >= 42.0 && AgeExact < 45.0 && Scr43 == 0) {
        Scr43 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        } else if (HPVDNA == 1 && CurrYear >= ImplementYR) {
          if (HPVDNAThermal == 0) {
            HPVScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                               SIV, SId, SIId, SIIId, SIVd);
          } else {
            if (ThermalORPap < PropThermal) {
              HPV_ThermalScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI,
                                         SII, SIII, SIV, SId, SIId, SIIId,
                                         SIVd);
            } else {
              ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                              SIV, SId, SIId, SIIId, SIVd, acceptRand,
                              efficacyD1Rand, efficacyD2Rand, D2LossRand);
            }
          }
        }
      }
      if (AgeExact >= 45.0 && AgeExact < 48.0 && Scr46 == 0) {
        Scr46 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        } else if (HPVDNA == 1 && CurrYear >= ImplementYR) {
          if (HPVDNAThermal == 0) {
            HPVScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                               SIV, SId, SIId, SIIId, SIVd);
          } else {
            if (ThermalORPap < PropThermal) {
              HPV_ThermalScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI,
                                         SII, SIII, SIV, SId, SIId, SIIId,
                                         SIVd);
            } else {
              ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                              SIV, SId, SIId, SIIId, SIVd, acceptRand,
                              efficacyD1Rand, efficacyD2Rand, D2LossRand);
            }
          }
        }
      }
      if (AgeExact >= 48.0 && AgeExact < 51.0 && Scr49 == 0) {
        Scr49 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        } else if (HPVDNA == 1 && CurrYear >= ImplementYR) {
          if (HPVDNAThermal == 0) {
            HPVScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                               SIV, SId, SIId, SIIId, SIVd);
          } else {
            if (ThermalORPap < PropThermal) {
              HPV_ThermalScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI,
                                         SII, SIII, SIV, SId, SIId, SIIId,
                                         SIVd);
            } else {
              ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                              SIV, SId, SIId, SIIId, SIVd, acceptRand,
                              efficacyD1Rand, efficacyD2Rand, D2LossRand);
            }
          }
        }
      }
      if (AgeExact >= 51.0 && AgeExact < 54.0 && Scr52 == 0) {
        Scr52 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        } else if (HPVDNA == 1 && CurrYear >= ImplementYR) {
          if (HPVDNAThermal == 0) {
            HPVScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                               SIV, SId, SIId, SIIId, SIVd);
          } else {
            if (ThermalORPap < PropThermal) {
              HPV_ThermalScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI,
                                         SII, SIII, SIV, SId, SIId, SIIId,
                                         SIVd);
            } else {
              ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                              SIV, SId, SIId, SIIId, SIVd, acceptRand,
                              efficacyD1Rand, efficacyD2Rand, D2LossRand);
            }
          }
        }
      }
      if (AgeExact >= 54.0 && AgeExact < 57.0 && Scr55 == 0) {
        Scr55 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        } else if (HPVDNA == 1 && CurrYear >= ImplementYR) {
          if (HPVDNAThermal == 0) {
            HPVScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                               SIV, SId, SIId, SIIId, SIVd);
          } else {
            if (ThermalORPap < PropThermal) {
              HPV_ThermalScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI,
                                         SII, SIII, SIV, SId, SIId, SIIId,
                                         SIVd);
            } else {
              ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                              SIV, SId, SIId, SIIId, SIVd, acceptRand,
                              efficacyD1Rand, efficacyD2Rand, D2LossRand);
            }
          }
        }
      }
      if (AgeExact >= 57.0 && AgeExact < 60.0 && Scr58 == 0) {
        Scr58 = 1;
        if (HPVDNA == 0 || CurrYear < ImplementYR) {
          ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII, SIV,
                          SId, SIId, SIIId, SIVd, acceptRand, efficacyD1Rand,
                          efficacyD2Rand, D2LossRand);
        } else if (HPVDNA == 1 && CurrYear >= ImplementYR) {
          if (HPVDNAThermal == 0) {
            HPVScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                               SIV, SId, SIId, SIIId, SIVd);
          } else {
            if (ThermalORPap < PropThermal) {
              HPV_ThermalScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI,
                                         SII, SIII, SIV, SId, SIId, SIIId,
                                         SIVd);
            } else {
              ScreenAlgorithm(ID, rea, ade, tts, res, ttC, CCd, SI, SII, SIII,
                              SIV, SId, SIId, SIIId, SIVd, acceptRand,
                              efficacyD1Rand, efficacyD2Rand, D2LossRand);
            }
          }
        }
      }
    }
  }
  void Pop::SaveCancerDeaths(const char *filout) {
    int ia, is;
    ostringstream s;

    if (process_num > 0) {
      s << process_num << "_" << filout;
    } else {
      s << filout;
    }

    string path = "./output/" + s.str();
    ofstream file(path.c_str()); // Converts s to a C string

    for (ia = 0; ia < 54; ia++) {
      for (is = 0; is < 136; is++) {
        file << right << NewCancerDeathHIV[ia][is] << "	";
      }
      file << endl;
    }
    for (ia = 0; ia < 18; ia++) {
      for (is = 0; is < 136; is++) {
        file << right << NewDiagCancerDeath[ia][is] << "	";
      }
      file << endl;
    }
    for (ia = 0; ia < 18; ia++) {
      for (is = 0; is < 136; is++) {
        file << right << NewCancerDeath[ia][is] << "	";
      }
      file << endl;
    }
    for (ia = 0; ia < 18; ia++) {
      for (is = 0; is < 136; is++) {
        file << right << HIVDeath[ia][is] << "	";
      }
      file << endl;
    }
    for (ia = 0; ia < 18; ia++) {
      for (is = 0; is < 136; is++) {
        file << right << HIVDeathM[ia][is] << "	";
      }
      file << endl;
    }

    file.close();
  }

  void Pop::SaveWeeksInStage(const char *filout) {
    int iy, is;

    ostringstream s;

    if (process_num > 0) {
      s << process_num << "_" << filout;
    } else {
      s << filout;
    }

    string path = "./output/" + s.str();
    ofstream file(path.c_str()); // Converts s to a C string

    for (iy = 0; iy < 54; iy++) {
      for (is = 0; is < 101; is++) {
        file << right << RSApop.WeeksInStageI[iy][is] << "	";
      }
      file << endl;
    }
    for (iy = 0; iy < 54; iy++) {
      for (is = 0; is < 101; is++) {
        file << right << RSApop.WeeksInStageII[iy][is] << "	";
      }
      file << endl;
    }
    for (iy = 0; iy < 54; iy++) {
      for (is = 0; is < 101; is++) {
        file << right << RSApop.WeeksInStageIII[iy][is] << "	";
      }
      file << endl;
    }
    for (iy = 0; iy < 54; iy++) {
      for (is = 0; is < 101; is++) {
        file << right << RSApop.WeeksInStageIV[iy][is] << "	";
      }
      file << endl;
    }
    file.close();
  }

  void Pop::SavePopPyramid9(const char *filout) {
    int ia, is;
    ostringstream s;

    if (process_num > 0) {
      s << process_num << "_" << filout;
    } else {
      s << filout;
    }

    string path = "./output/" + s.str();
    ofstream file(path.c_str()); // Converts s to a C string

    for (ia = 0; ia < 6; ia++) {
      for (is = 0; is < 136; is++) {
        file << right << PopPyramid9[ia][is] << "	";
      }
      file << endl;
    }

    file.close();
  }
  void Pop::SaveStagediag(const char *filout) {
    int iy, is;

    ostringstream s;

    if (process_num > 0) {
      s << process_num << "_" << filout;
    } else {
      s << filout;
    }

    string path = "./output/" + s.str();
    ofstream file(path.c_str()); // Converts s to a C string

    for (iy = 0; iy < 54; iy++) {
      for (is = 0; is < 136; is++) {
        file << right << RSApop.StageIdiag[iy][is] << "	";
      }
      file << endl;
    }
    for (iy = 0; iy < 54; iy++) {
      for (is = 0; is < 136; is++) {
        file << right << RSApop.StageIIdiag[iy][is] << "	";
      }
      file << endl;
    }
    for (iy = 0; iy < 54; iy++) {
      for (is = 0; is < 136; is++) {
        file << right << RSApop.StageIIIdiag[iy][is] << "	";
      }
      file << endl;
    }
    for (iy = 0; iy < 54; iy++) {
      for (is = 0; is < 136; is++) {
        file << right << RSApop.StageIVdiag[iy][is] << "	";
      }
      file << endl;
    }
    file.close();
  }
  void ReadCCStrategies() {
    std::ifstream f("CCPreventionStrategiesJ.json");
    json data = json::parse(f);

    data.at("ImplementYR").get_to(ImplementYR);
    data.at("HPVvacc").get_to(HPVvacc);
    data.at("BOYSvacc").get_to(BOYSvacc);
    data.at("RoutineScreening").get_to(RoutineScreening);
    data.at("CCcalib").get_to(CCcalib);
    data.at("PerfectSchedule").get_to(PerfectSchedule);
    data.at("HPVDNA").get_to(HPVDNA);
    data.at("MnDCampaign").get_to(MnDCampaign);
    data.at("InvCaseAlgorithm").get_to(InvCaseAlgorithm);
    data.at("HPVDNAThermal").get_to(HPVDNAThermal);
    data.at("PropThermal").get_to(PropThermal);
    data.at("HPVGenotyping").get_to(HPVGenotyping);
    data.at("PapTRIAGE").get_to(PapTRIAGE);
    data.at("Portal3").get_to(Portal3);
    data.at("WHOscenario").get_to(WHOscenario);
    data.at("WHOvacc").get_to(WHOvacc);
    data.at("WHOScreening").get_to(WHOScreening);
    data.at("S3S4").get_to(S3S4);
    data.at("S5S6").get_to(S5S6);
    data.at("S7S11").get_to(S7S11);
    data.at("CatchUpVaccHIV").get_to(CatchUpVaccHIV);
    data.at("CatchUpVacc").get_to(CatchUpVacc);
    data.at("CatchUpAgeMIN").get_to(CatchUpAgeMIN);
    data.at("CatchUpAgeMAX").get_to(CatchUpAgeMAX);
    data.at("CatchUpCoverage").get_to(CatchUpCoverage);
    data.at("VaccineWane").get_to(VaccineWane);
    data.at("VaccDur").get_to(VaccDur);
    data.at("AdministerMassTxV").get_to(AdministerMassTxV);
    data.at("MassTxVAgeMIN").get_to(MassTxVAgeMIN);
    data.at("MassTxVAgeMAX").get_to(MassTxVAgeMAX);
    data.at("AdministerMassTxVtoART").get_to(AdministerMassTxVtoART);
    data.at("MassTxVtoARTAgeMIN").get_to(MassTxVtoARTAgeMIN);
    data.at("MassTxVtoARTAgeMAX").get_to(MassTxVtoARTAgeMAX);
    data.at("ReductionFactorHIV").get_to(ReductionFactorHIV);
    data.at("ReductionFactorHIVART").get_to(ReductionFactorHIVART);
    data.at("TxVviaScreeningAlgorithm").get_to(TxVviaScreeningAlgorithm);
	data.at("TxVtoHPVPos").get_to(TxVtoHPVPos);
    data.at("UpdateStart").get_to(UpdateStart);
    data.at("CreateCohort").get_to(CreateCohort);
    data.at("UseMedians").get_to(UseMedians);
    data.at("OneType").get_to(OneType);
    data.at("WhichType").get_to(WhichType);
    data.at("targets").get_to(targets);
    data.at("GetMacD").get_to(GetMacD);
  }

  bool Indiv::AnyHPV(const int *XXX, const vector<int> &type_subset,
                     const vector<int> &stage_subset) {

    for (int type : type_subset) {
      for (int stage : stage_subset) {
        if (XXX[type] == stage) {
          return true;
        }
      }
    }
    return false;
  }