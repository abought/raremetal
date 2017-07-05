////////////////////////////////////////////////////////////////////// 
// TransformResiduals.h 
// (c) 2012-2013 Shuang Feng, Dajiang Liu, Goncalo Abecasis
// 
// This file is distributed as part of the RareMetalWorker source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile RareMetalWorker.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Wednesday November 28, 2012
// 

#ifndef __TRANSFORMRESIDULAS_H__
#define __TRANSFORMRESIDULAS_H__

#include "Kinship.h"
#include "KinshipX.h"
#include "KinshipEmp.h"
#include "MathNormal.h"
#include "MathVector.h"
#include "StringBasics.h"
#include "MathMatrix.h"
#include "MathCholesky.h"
#include "AutoFit.h"
#include <Eigen/Eigenvalues>
#include <savvy/savvy.hpp>

class FastTransform
{
	public:
		FastTransform();
		~FastTransform();

		static bool empKin;
		static bool pedKin;
		static bool mergedVCFID;
		static String readInEmp;
		static String readInEmpX;

		//# of useful families and individuals
		int families,persons,fixedEffects,analyzedFounders,numFounder,totalN,numFamily;
		double traitVar,traitMean;

		//Hash of sample ids from vcf file, if reading from vcf is needed
		StringIntHash sampleVCFIDHash,samplePEDIDHash,foundersHash;

		IntArray * pPheno;
		IntArray genotypedSamplePED, genotypedSampleVCF;
		StringArray founders;
		//Saved Matrices for likelihood calculation
		Eigen::MatrixXf transU; 
		Matrix transU_del;
		Matrix X,UX;
		Vector D,UY,Y; 
		//Matrix UDX,UD,inv; //here D is D^1/2
		//Vector UDY; //D is D^1/2

		//Reformat LL* into UDU*U
		void CleanUp()
		{
		   //X.Dimension(0,0);
		   UX.Dimension(0,0);
		};

		int GetCovarianceMatrix(AutoFit & engine, Pedigree &ped, int f, int count,Matrix & omega);
		void ReformCholesky(Matrix & L, Vector & D);
		int GetTraitID(Pedigree & ped, const char * name);
		//These functions are transformation work horses
		void GetFixedEffectsCount(Pedigree & ped, bool useCovariates);
		void Prepare(Pedigree & ped, int traitNum,bool useCovariates,FILE * log,bool shortVersion);
		void SubSet(Matrix & kin,Pedigree & ped,Matrix & tmp);
		void SelectSamplesPED(Pedigree & ped,bool useCovariates);
	  template<typename VecType>
		void SelectSamplesVCF(const String& path, Pedigree & ped, bool useCovariates);
		void ScreenSampleID(Pedigree & ped, bool useCovariates);

		void TransformPedkin(Pedigree & ped, int traitNum,bool useCovaraites,FILE * log);
		void TransformPedkinSepX(Pedigree & ped, int traitNum,bool useCovaraites,FILE * log);
		void TransformEmpkin(Pedigree & ped, int traitNum,bool useCovariates,KinshipEmp & kin_emp,Matrix & input,FILE * log);
		void TransformPedkinX(AutoFit & engine, Pedigree &ped);
		void TransformEmpkinX(Matrix & covMatrix);
		//void FinalizeTransform(double sigma_g2,double sigma_e2);
		int EigenDecompose(Matrix & matrix_in, Matrix & U, Vector & D);
		bool FinalizeProducts();
};



// 07/08/16 update: read vcf once instead of multiple times
template<typename VecType>
void FastTransform::SelectSamplesVCF(const String& path, Pedigree & ped,bool useCovariates)
{
	// printf("Selecting useful samples ...\n");
	savvy::reader reader(path.c_str());
	VecType record;


	int numSamples = reader.samples_end() - reader.samples_begin();
	totalN=numSamples;

	StringIntHash VCFID;
	for(auto s = reader.samples_begin(); s != reader.samples_end(); ++s) {
		VCFID.SetInteger(s->c_str(), s - reader.samples_begin());
	}

	// now see if each sample is genotyped in at least one site
	std::map<int, bool> g1; // store the samples that haven't found a genotyped site
	for(int s=0; s<numSamples; s++)
		g1[s] = 0;

	while(reader >> record) {
		if (g1.empty())
			break;

		std::map<int, bool> to_remove;
		for(std::map<int, bool>::iterator t=g1.begin(); t!=g1.end(); t++) {
			int s = t->first;
			if(std::is_same<savvy::dense_dosage_vector<float>, VecType>::value) {
				float geno = record[s];
				if(std::isnan(geno))
					to_remove[s] = 1;
			}
			else { // gt
				int numGTs = savvy::get_ploidy(reader, record);
				bool bad = 0;
				for(int j = 0; j < numGTs; j++) {
					if(std::isnan(record[s * numGTs + j])) {
						bad = 1;
						break;
					}
				}
				if (!bad)
					to_remove[s] = 1;
			}
		}
		for(std::map<int, bool>::iterator pm=to_remove.begin(); pm!=to_remove.end(); pm++)
			g1.erase(pm->first);
	}

	// now g1 has index of samples that are not genotyped
	for(int i=0; i<ped.count; i++) {
		bool phenotyped =false;
		for(int tr=0;tr<ped.traitNames.Length();tr++) {
			if(ped[i].isPhenotyped(tr)) {
				phenotyped=true;
				break;
			}
		}
		if (!phenotyped)
			continue;
		int s;
		String sample;
		if(mergedVCFID) {
			sample = ped[i].famid+"_"+ped[i].pid;
			s = VCFID.Integer(sample);
		}
		else
			s = VCFID.Integer(ped[i].pid);
		if (s==-1)
			continue;
		if (g1.find(s) != g1.end())
			continue;
		if(mergedVCFID) {
			sampleVCFIDHash.SetInteger(sample,s);
			samplePEDIDHash.SetInteger(sample,i);
		}
		else {
			sampleVCFIDHash.SetInteger(ped[i].pid,s);
			samplePEDIDHash.SetInteger(ped[i].pid,i);
		}
		genotypedSampleVCF.Push(s);
	}
}

#endif
