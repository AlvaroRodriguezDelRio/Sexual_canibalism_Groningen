// Alvaro Rodríguez del Río
// C++ script for simulating the evolution of allele frequency in a haploid spider population with two alleles (mono for monogyny and bi fot bigyny) for the "mating strategy" gene
// In this simulations I keep the male and female population size constant by replacing each female and male by a new one everytime there is a death event


#include <cassert>
#include <string>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>

using namespace std;

const double Tmax = 100000000;
const double sexRatio = 5;  //males/females
const int numberFemales = 1000;
const int numberMales = (int) numberFemales * sexRatio;
const double intrinsicMortality = 0.1; // intrinsic mortality of males
const double spermTransferByMonogynists = 3; //Monogynist males can transfer more sperm than biginist males, which transfer 1 sperm load
const double travelMortality = 0.1; //When looking for a female, males may die because of extrinsic mortality
const int NumberOfPedipalps = 2; // Spiders have two pedipalps. Monogynist males use both with the same female, whereas byginist males use one per female
const double decreaseProbMatingByMonogyny = 0; //Monogynists may decrease the probability of the female to remate
const int numberOfEvents = 2;// this is an event based model. This 2 events are 1. male dies and 2. copulation
double z = 0.5; //proportion of mono/bi sperm
const int windowSize = 20; // the z value is calculated by meassuring the sperm of the 20 last females which laid eggs
const double probabilityFemaleDiesAfterCopulation = 0.20;
const double mutationProbability = 0.001;

//functions which may be used in the class definition, I did not use all of them in this script, but I am using them in others

double GetValueNormal (double mean, double sd){//gets a random value from a normal distribution with mean and sd provided
  mt19937_64 rng;
  normal_distribution <double> Normal (mean,sd);
  std::chrono::high_resolution_clock::time_point tp = std::chrono::high_resolution_clock::now();
  unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
  rng.seed(seed);
  return Normal(rng);
}

int getValueFromUniformDistrinution(int begining, int ending){
  mt19937_64 rng;
  uniform_int_distribution<int> Uniform (begining,ending);
  std::chrono::high_resolution_clock::time_point tp = std::chrono::high_resolution_clock::now();
  unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
  rng.seed(seed);
  int chosen = Uniform(rng);
  return chosen;
}

int getValueBernuilliDistribution(double factor){
  mt19937_64 rng;
  bernoulli_distribution Berm (factor);
  std::chrono::high_resolution_clock::time_point tp = std::chrono::high_resolution_clock::now();
  unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
  rng.seed(seed);
  int chosen = Berm(rng);
  return chosen;
}

double getValueExponentialDistribution(double factor){
  mt19937_64 rng;
  exponential_distribution<double> exponential(factor);
  std::chrono::high_resolution_clock::time_point tp = std::chrono::high_resolution_clock::now();
  unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
  rng.seed(seed);
  double chosen = exponential(rng);
  return chosen;
}

//creation of the variable type "Allele"

enum Allele {mono,bi};

//Class declaration

class male{
  public:
    male(){allele = mono;matingsLeft = NumberOfPedipalps;life = getValueExponentialDistribution(intrinsicMortality);};
    Allele getAllele (){return allele;};
    void setAllele (Allele GetAllele){allele = GetAllele;};
    int getMatingsLeft(){return matingsLeft;};
    void decreaseMatingsLeft(int matings){matingsLeft-=matings;};
    double getLife(){return life;};//males have a life taken from an exp distribution
    int getWhenDidIAppear(){return whenDidIappear;};//they are given a "birthday when appearing in the population"
    void setwhenDidIappear(int when){whenDidIappear = when;};
  private:
    Allele allele;
    int matingsLeft;
    double life;
    int whenDidIappear;
};

class female{
  public:
    female(){probabilityMate = 1;amountMonoSperm = 0;amountBiSperm = 0;matedByMono=false;};
    double getProbToMate(){return probabilityMate;};
    void decreaseProbToMate(double decreaseMatingPro){if (probabilityMate>0){probabilityMate -= decreaseMatingPro;}};
    int getAmountMonoSperm(){return amountMonoSperm;};
    void AddMonoSperm(int amount){amountMonoSperm+=amount;};// at the beginnign I was using a vector for string all the sperm, but I thought this was a simpler implementation
    int getAmountBiSperm(){return amountBiSperm;};
    void AddBiSperm(int amount){amountBiSperm+=amount;};
    bool getMatedByMono (){return matedByMono;};
    void changeMatedByMono (bool change){matedByMono = change;};
  private:
    double probabilityMate;
    int amountMonoSperm;
    int amountBiSperm;
    double life;
    bool matedByMono;
};

//creates a population, in this case, whith 50% mono males and 50% bi males
void setFirstPopulation(vector<male> &malePopulation,vector<female> &femalePopulation){
  for (int i = 0; i < malePopulation.size();++i){
    int whichAllele = getValueBernuilliDistribution(0.5);
    if (whichAllele == 0){
      malePopulation[i].setAllele(mono);
    }
    else{
      malePopulation[i].setAllele(bi);
    }
    malePopulation[i].setwhenDidIappear(0);
  }
}

//incorporates a new male into the population
void NewMale(vector <male> &malePopulation, int whichMale, vector <female> femalesWhoLaidEggs, int timeMeasurement){
  male in;
  if (femalesWhoLaidEggs.size() > 0){
    int whichFemale = getValueFromUniformDistrinution(0,femalesWhoLaidEggs.size()-1);
    int monoAlleles = femalesWhoLaidEggs[whichFemale].getAmountMonoSperm();
    int biAlleles = femalesWhoLaidEggs[whichFemale].getAmountBiSperm();
    double zFemale = (double)monoAlleles/(monoAlleles+biAlleles);
    int whichAllele = getValueBernuilliDistribution(zFemale);
    if (whichAllele == 1){
      in.setAllele(mono);
      int IsThereMutation = getValueBernuilliDistribution(mutationProbability);
      if (IsThereMutation == 1){
        in.setAllele(bi);
      }
    }
    else{
      in.setAllele(bi);
      int IsThereMutation = getValueBernuilliDistribution(mutationProbability);
      if (IsThereMutation == 1){
        in.setAllele(mono);
      }
    }
  }
  else{// if there are not available females for laying eggs I randomly set the allele for the new male
    int whichAllele = getValueBernuilliDistribution(0.5);
    if (whichAllele == 1){
      in.setAllele(mono);
    }
    else{
      in.setAllele(bi);
    }
  }
  in.setwhenDidIappear(timeMeasurement);
  malePopulation[whichMale] = in;//substitutes the old male position
}

int main(){

  vector<female> femalePopulation(numberFemales);
  vector<male> malePopulation(numberMales);

  setFirstPopulation(malePopulation,femalePopulation);

  vector<female> femalesWhoLaidEggs(0);

  //I usually prefer to print results to the screen and redirect the output to a text file instead of creating files within the script

  cout << "#sex ratio: " << sexRatio << " increased Fecundity By Cannibalism: " << spermTransferByMonogynists << ", Travel Mortality = "<< travelMortality << "intrinsic mortality = "<< intrinsicMortality <<" Decreased prob to mate by monogyny = " << decreaseProbMatingByMonogyny << endl;
  cout << "proportionMono,intrinsicMortality,travelMortality,sexRatio,spermTransferByMonogynists,decreaseProbMatingByMonogyny" << endl;

  int print = 0;//variable which controls when to print results

  int timeMeasurement = 0;//the time scale was set so that evertime that, n average, each male mated once, the time increased by 1
  int numberMaleTrials = 0;//controls the number of times males tried to copulate in order to control the timeMeasurement

  for (int i = 1; i < Tmax; ++i){

    //time control
    bool resetnumberMaleTrials = false;
    if (numberMaleTrials==numberFemales){
      timeMeasurement++;
      resetnumberMaleTrials = true;
    }
    if(resetnumberMaleTrials){
      numberMaleTrials=0;
    }

    //events
    int event = getValueFromUniformDistrinution(1,numberOfEvents);
    if (event == 1){//male dies
      bool dies = false;
      int whichMale = getValueFromUniformDistrinution(0,numberMales-1);
      if (timeMeasurement - malePopulation[whichMale].getWhenDidIAppear() > malePopulation[whichMale].getLife()){//dies because he is old
        NewMale(malePopulation,whichMale,femalesWhoLaidEggs,timeMeasurement);
      }
    }
    else if (event ==2){//male tries to copulates a female
      numberMaleTrials++;
      int whichMale = getValueFromUniformDistrinution(0,numberMales-1);
      double doesHeDieWhileSearching = getValueBernuilliDistribution(travelMortality);
      if (doesHeDieWhileSearching == 1){//dies because of travel mortality
        NewMale(malePopulation, whichMale, femalesWhoLaidEggs,timeMeasurement);
        continue;
      }
      int whichFemale = getValueFromUniformDistrinution(0,numberFemales-1);
      int DoesTheFemaleWantToMate = getValueBernuilliDistribution(femalePopulation[whichFemale].getProbToMate());
      if (DoesTheFemaleWantToMate == 1){//copulation happens
        Allele maleAllele = malePopulation[whichMale].getAllele();
        if (maleAllele == mono){//if the male is mono
          femalePopulation[whichFemale].AddMonoSperm(spermTransferByMonogynists);
          malePopulation[whichMale].decreaseMatingsLeft(NumberOfPedipalps);
          if (!femalePopulation[whichFemale].getMatedByMono()){//only the first mono male decreases the probability of the female to remate
            femalePopulation[whichFemale].decreaseProbToMate(decreaseProbMatingByMonogyny);
          }
          femalePopulation[whichFemale].changeMatedByMono(true);
        }
        else{//if the male is bi
          femalePopulation[whichFemale].AddBiSperm(1);
          malePopulation[whichMale].decreaseMatingsLeft(1);
        }
        if (malePopulation[whichMale].getMatingsLeft()==0){//dies because he does not have any mating left
          NewMale(malePopulation,whichMale,femalesWhoLaidEggs,timeMeasurement);
        }
        int doesSheDie = getValueBernuilliDistribution(probabilityFemaleDiesAfterCopulation);//Evaluates if the female dies after copulation
        if (doesSheDie==1){//if she dies s
          femalesWhoLaidEggs.push_back(femalePopulation[whichFemale]);
          if (femalesWhoLaidEggs.size() > windowSize){//change female laying eggs to fit the size
              femalesWhoLaidEggs.erase(femalesWhoLaidEggs.begin());
              int ammountMono = 0;
              int ammountBi = 0;
              for (int j = 0; j < femalesWhoLaidEggs.size(); j++){
                ammountMono = ammountMono + femalesWhoLaidEggs[j].getAmountMonoSperm();
                ammountBi = ammountBi + femalesWhoLaidEggs[j].getAmountBiSperm();
              }
              double proportion = (double)ammountMono/(ammountMono+ammountBi);
              z = proportion;
              female in;
              femalePopulation[whichFemale] = in;
          }
        }
      }
    }

    //prints results to screen every 5000 steps
    print++;
    if (print >= 5000){
      cout  << z << "," << intrinsicMortality << "," << travelMortality << "," << sexRatio << "," << spermTransferByMonogynists << "," << decreaseProbMatingByMonogyny << endl;
      print=0;
    }
  }
}
