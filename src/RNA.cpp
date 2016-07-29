
//Class RNA -- Wrapper for RNAstructure for use in object-oriented applications
#include "RNA.h"
#include "structure.h"
#include "arrayclass.h"
#include "forceclass.h"
#include "rna_library.h"
#include "stackclass.h"
#include "stackstruct.h"
#include "algorithm.h"
#include "outputconstraints.h"
#include "boltzmann.h"
#include "alltrace.h"
#include <iostream>

const float epsilon = 1e-6; // a small number for a tolerance in comparing floats
//constructor where user provides a string with the sequence
RNA::RNA(const char sequence[], const bool IsRNA) : Thermodynamics(IsRNA) {
  int i;
  //allocate ct
  ct = new structure();
  //Specify the sequence length based on the string length
  //ct->numofbases = (short) strlen(sequence);
  ct->allocate((int) strlen(sequence));//allocate the space required in the arrays for the sequence
  //Now store the sequence information
  for (i = 1; i <= ct->GetSequenceLength(); i++) {
    if (sequence[i - 1] == 'A' || sequence[i - 1] == 'a') ct->numseq[i] = 1;
    else if (sequence[i - 1] == 'C' || sequence[i - 1] == 'c') ct->numseq[i] = 2;
    else if (sequence[i - 1] == 'G' || sequence[i - 1] == 'g') ct->numseq[i] = 3;
    else if (sequence[i - 1] == 'U' || sequence[i - 1] == 'u' || sequence[i - 1] == 'T' || sequence[i - 1] == 't')
      ct->numseq[i] = 4;
    else ct->numseq[i] = 0;
    ct->nucs[i] = sequence[i - 1];
    ct->hnumber[i] = i;
  }
  //These should not be located here: (DHM 5/3/2011)
  // Experimental bonuses array
  //ct->EX = new double *[ct->GetSequenceLength()+1];
  //for(i=0;i<ct->GetSequenceLength()+1;i++) {
  //  ct->EX[i] = new double[ct->GetSequenceLength()+1];
  //}
  //for (i=0;i<ct->GetSequenceLength()+1;i++){
  //  for (j=0;j<ct->GetSequenceLength()+1;j++){
  //   ct->EX[i][j] = 0.0;
  // }
  //   }
  //set error status to zero
  ErrorCode = 0;
  //Do not report progress by default:
  progress = NULL;
}

//constructor where user provides a string with a filename
//	The flag type indicates the file type: type = 1 => ct file, type = 2 => .seq file, type = 3 => .pfs file.
//	The fconstructor saves an error code in ErrorCode:
//	0 = no error, 1 = file not found
//	2 = error opening file.
//  GetErrorCode provides public access to this errorcode.
RNA::RNA(const char filename[], const int type, const bool IsRNA) : Thermodynamics(IsRNA) {
  //allocate ct
  ct = new structure();
  //Do not report progress by default:
  progress = NULL;
  ErrorCode = FileReader(filename, type);
  return;
}

//Default constructor.
RNA::RNA(const bool IsRNA) : Thermodynamics(IsRNA) {
  //allocate the underlying structure class and nothing more.
  //User is then required to propogate the sequence and structural information manually.
  ct = new structure();
  //set the number of nucleotides to zero because no sequence has been read
  //ct->numofbases = 0;
  //Drawing coordinates have not been determined.
  ErrorCode = 0;
  //Do not report progress by default:
  progress = NULL;
}

//Return the value of ErrorCode
int RNA::GetErrorCode() {
  return ErrorCode;
}

//Return a c string that describes errors from GetErrorCode and other errors.
char* RNA::GetErrorMessage(const int error) {
  if (error == 0) return "No Error.\n";
  else if (error == 1) return "Input file not found.\n";
  else if (error == 2) return "Error opening file.\n";
  else if (error == 3) return "Structure number out of range.\n";
  else if (error == 4) return "Nucleotide number out of range.\n";
  else if (error == 5)
    return "Error reading thermodynamic parameters.\nPlease set environment variable DATAPATH to the location of the thermodynamic parameters.\n";
  else if (error == 6) return "This would form a pseudoknot and is not allowed.\n";
  else if (error == 7) return "This pair is non-canonical and is therefore not allowed.\n";
  else if (error == 8) return "Too many restraints specified.\n";
  else if (error == 9) return "This nucleotide already under a conflicting constraint.\n";
  else if (error == 10) return "There are no structures to write to file.\n";
  else if (error == 11) return "Nucleotide is not a U.\n";
  else if (error == 12) return "Maximum pairing distance is too short.\n";
  else if (error == 13) return "Error reading constraint file.\n";
  else if (error == 14) return "A traceback error occurred.\n";
  else if (error == 15) return "No partition function data is available.\n";
  else if (error == 16) return "Wrong save file version used or file format not recognized.\n";
  else if (error == 17)
    return "This function cannot be performed unless a save file (.sav) was correctly loaded by the RNA constructor.\n";
  else if (error == 18) return "This threshold is too low to generate valide secondary structures.\n";
  else if (error == 19)
    return "The structure coordinates have not been determined, use DetermineDrawingCoordinates() to calculate the coordinates.\n";
  else if (error == 20) return "No sequence has been read.\n";
  else if (error == 21) return "Probabilities summed to greater than 1 in stochastic traceback.\n";
  else if (error == 22) return "Programming error.  Incorrect file type passed to constructor.\n";
  else if (error == 23) return "There are no structures present.\n";
  else if (error == 24) return "Too few iterations.  There must be at least one iteration.\n";
  else if (error == 25) return "Index is not a multiple of 10.\n";
  else if (error == 26) return "k, the equilibrium constant, needs to be greater than or equal to 0.\n";
  else if (error == 27) return "Lyngso O(N^3) internal loop search is not compatible with a parallel calculation\n";
  else return "Unknown Error\n";
}

//Return a string that describes errors from GetErrorCode and other errors.
//This uses GetErrorMessage to acrually get the errors.
std::string RNA::GetErrorMessageString(const int error) {
  std::string temp;
  temp = GetErrorMessage(error);
  return temp;
}

void RNA::ResetError() {
  this->ErrorCode = 0;
}

//User specifies a base pair between i and j in structure # structurenumber.
int RNA::SpecifyPair(const int i, const int j, const int structurenumber) {
  //start with error checking:
  if (i < 0 || i > ct->GetSequenceLength() || j < 0 || j > ct->GetSequenceLength()) return 4;
  else if (structurenumber < 1) return 3;
  //also keep track of the maximum number of structures for which pairs have been specified
  if (structurenumber > ct->GetNumberofStructures()) {
    //Add one structure for each position between those available and the one being specified:
    for (int index = ct->GetNumberofStructures() + 1; index <= structurenumber; ++index) ct->AddStructure();
  }
  //now register the pair:
  ct->SetPair(i, j, structurenumber);
  return 0;
}

//Break a pair that i is involved in
int RNA::RemoveBasePair(const int i, const int structurenumber) {
  //start with error checking:
  if (i < 0 || i > ct->GetSequenceLength()) return 4;
  else if (structurenumber < 1 || structurenumber > ct->GetNumberofStructures()) return 3;
  //Call the function for this in the underlying structure class.
  ct->RemovePair(i, structurenumber);
  //return that there was no error
  return 0;
}

//remove all pairs in structure # structurenumber.
//Also, roll back the number of specified structures if this is the last specified structure.
int RNA::RemovePairs(const int structurenumber) {
  //do some error checking
  if (structurenumber > ct->GetNumberofStructures() || structurenumber < 1) return 5;
  //decrement the number of structures, if appropriate, i.e. this is the last structure
  if (structurenumber == ct->GetNumberofStructures()) {
    ct->RemoveLastStructure();
    return 0;
  }
  //otherwise, clean the selected structure of pairs:
  ct->CleanStructure(structurenumber);
  return 0;
}

//Calculate and return the folding free energy change for structure number structurenumber.
double RNA::CalculateFreeEnergy(const int structurenumber, const bool UseSimpleMBLoopRules) {
  //Do some simple error checking
  if (structurenumber < 1 || structurenumber > ct->GetNumberofStructures()) return 0.0;
  if (!energyread) {
    //The thermodynamic data tables have not yet been read
    if (ReadThermodynamic() != 0) {
      ErrorCode = 5;//record an error
      return 0.0;//return 0.0 if a problem occurs
    }
    else ErrorCode = 0;//record that there was no error.
  }
  else ErrorCode = 0;//Set the error code to zero because no errors were encountered.
  efn2(data, ct, structurenumber, UseSimpleMBLoopRules);
  //conversion factor is set in defines.h.  Free energies are multiplied by this factor internally so that integer math can be used.
  return (((double) ct->GetEnergy(structurenumber) / conversionfactor));
}

//Write the details on the energy caclulation for all structures.
int RNA::WriteThermodynamicDetails(const char filename[], const bool UseSimpleMBLoopRules) {
  if (!energyread) {
    //The thermodynamic data tables have not yet been read
    if (ReadThermodynamic() != 0) return 5;//return non-zero if a problem occurs
  }
  efn2(data, ct, 0, UseSimpleMBLoopRules, filename);
  return 0;
}

#ifndef DYNALIGN_II

//Predict the secondary structure by free energy minimization.
//Also generate subooptimal solutions using a heuristic.
int RNA::FoldSingleStrand(const float percent, const int maximumstructures, const int window, const char savefile[],
                          const int maxinternalloopsize, bool mfeonly, bool simple_iloops) {
  char* savefilename;
  int percenti;
  int tracebackstatus;
  //check to make sure that a sequence has been read
  if (ct->GetSequenceLength() == 0) return 20;
  if (!energyread) {
    //The thermodynamic data tables have not been read and need to be read now.
    if (ReadThermodynamic() != 0) return 5;//return non-zero if a problem occurs
  }
  //savefile will be passed to the function dynamic for structure prediction.
  //Dynamic expects a null pointer if no file is to be created, so set savefilename to null if savefile is an empty string.
  if (!strcmp(savefile, "")) savefilename = NULL;
  else {
    savefilename = new char[((int) strlen(savefile)) + 1];
    strcpy(savefilename, savefile);
  }
  //right now, dynamic requires an integer specification of percent energy change.
  //FoldSingleStrand takes this is a float to provide the opportunity to reform this in the future.
  //For now, cast percent as an integer to pass to dynamic.
  percenti = (int) percent;
  //Predict the secondary structures.
  tracebackstatus = dynamic(ct, data, maximumstructures, percenti, window, progress, false, savefilename,
                            maxinternalloopsize, mfeonly, simple_iloops);
  //Clean up the memory use.
  delete[] savefilename;
  if (tracebackstatus != 0) return 14;//This indicates a traceback error.
  else return 0;
}

#else
#endif

// Predict the lowest free energy secondary structure and generate all suboptimal structures.
int RNA::GenerateAllSuboptimalStructures(const float percent, const double deltaG) {
  //check to make sure that a sequence has been read
  if (ct->GetSequenceLength() == 0) return 20;
  if (!energyread) {
    //The thermodynamic data tables have not been read and need to be read now.
    if (ReadThermodynamic() != 0) return 5;//return non-zero if a problem occurs
  }
  //Call the alltrace function to do the work:
  alltrace(ct, data, ((short) percent), ((short) (deltaG * conversionfactor)), progress, NULL);
  return 0;
}

//Force a nucleotide to be double stranded (base paired).
//Return an integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified, 9 = same nucleotide in conflicting restraint).
int RNA::ForceDoubleStranded(const int i) {
  int index;
  //check to make sure that a sequence has been read
  if (ct->GetSequenceLength() == 0) return 20;
  //Check that nucleotide is valid
  if (i < 1 || i > ct->GetSequenceLength()) return 4;//i is out of range
  //Check for conflicting restraints (anything specifying i to be unpaired):
  for (index = 0; index < ct->GetNumberofSingles(); index++) {
    if (ct->GetSingle(index) == i) return 9;
  }
  ct->AddDouble(i);
  return 0;
}

//Function to specify a nucleotide, i, that is a U in a GU pair
//Returns an integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified, 9 = same nucleotide in conflicting restraint, 11 = nucleotide not U).
int RNA::ForceFMNCleavage(const int i) {
  int index;
  //check to make sure that a sequence has been read
  if (ct->GetSequenceLength() == 0) return 20;
  //Check that nucleotide is valid
  if (i < 1 || i > ct->GetSequenceLength()) return 4;//i is out of range
  //Check to make sure the nucleotide is U.
  if (ct->numseq[i] != 4) return 11;
  //Check for conflicting restraints (anything specifying i to be unpaired):
  for (index = 0; index < ct->GetNumberofSingles(); index++) {
    if (ct->GetSingle(index) == i) return 9;
  }
  //Check for a nucleotide already forced to be in a pair that is not a GU pair.
  for (index = 0; index < ct->GetNumberofPairs(); index++) {
    if (i == ct->GetPair5(index) && ct->numseq[ct->GetPair3(index)] != 3) return 9;
    else if (i == ct->GetPair3(index) && ct->numseq[ct->GetPair5(index)] != 3) return 9;
  }
  ct->AddGUPair(i);
  return 0;
}

//Specify the maximum distance allowed between paired nucleotides in subsequent structure prediction.
//return An integer that indicates an error code (0 = no error, 12 = too long or too short).
int RNA::ForceMaximumPairingDistance(const int distance) {
  //check to make sure that a sequence has been read
  if (ct->GetSequenceLength() == 0) return 20;
  if (distance < minloop + 1) return 12;
  else {
    ct->SetPairingDistance(distance);
    return 0;
  }
}

//Indicate a nucleotide that is accessible to chemical modification.
//Returns an integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified).
int RNA::ForceModification(const int i) {
  //check to make sure that a sequence has been read
  if (ct->GetSequenceLength() == 0) return 20;
  //Check that nucleotide is valid
  if (i < 1 || i > ct->GetSequenceLength()) return 4;//i is out of range
  //Go ahead and record the constraint.
  ct->AddModified(i);
  return 0;
}

//Force a base pair between nucleotides i and j.
//Returns an error code: (0 = no error, 4 = nucleotide out of range, 6 = pseudoknot formation, 7 = non-canonical pair, 8 = too many restraints specified, 9 = same nucleotide in conflicting restraint).
int RNA::ForcePair(const int i, const int j) {
  bool allowedpairs[6][6] = {{false, false, false, false, false, false},
                             {false, false, false, false, true,  false},
                             {false, false, false, true,  false, false},
                             {false, false, true,  false, true,  false},
                             {false, true,  false, true,  false, false},
                             {false, false, false, false, false, false}};
  int index;
  int locali, localj;
  //First perform the error checking:
  //check to make sure that a sequence has been read
  if (ct->GetSequenceLength() == 0) return 20;
  //Note: In structure, forced pairs run between index of 1 and a maximum of maxforce-1.
  //if (ct->npair==(maxforce-1)) return 8;//This means there are too many pair constraints.
  if (i < 1 || i > ct->GetSequenceLength()) return 4;//i is out of range
  if (j < 1 || j > ct->GetSequenceLength()) return 4;//j is out of range
  if (!allowedpairs[ct->numseq[i]][ct->numseq[j]]) return 7;//non-canonical pair
  //sort indexes from 5' to 3':
  if (i > j) {
    locali = j;
    localj = i;
  }
  else {
    locali = i;
    localj = j;
  }
  //check for pseudoknots with any other forced pair or the same nucleotide forced into two pairs:
  for (index = 0; index < ct->GetNumberofPairs(); index++) {
    if (locali < ct->GetPair5(index) && ct->GetPair5(index) < localj && localj < ct->GetPair3(index))
      return 6;//a pseudoknot
    if (locali == ct->GetPair5(index) || locali == ct->GetPair3(index) || localj == ct->GetPair5(index) ||
        localj == ct->GetPair3(index))
      return 9;//i or j is in another forced pair
  }
  //now check for other conflicting restraints:
  for (index = 0; index < ct->GetNumberofForbiddenPairs(); index++) {
    if (ct->GetForbiddenPair5(index) == locali && ct->GetForbiddenPair3(index) == localj)
      return 9;//The pair was forbidden.
  }
  for (index = 0; index < ct->GetNumberofSingles(); index++) {
    if (ct->GetSingle(index) == locali || ct->GetSingle(index) == localj)
      return 9;//i or j was previously forced single-stranded.
  }
  //Now register the restraint because the error checking was clear or errors.
  ct->AddPair(locali, localj);
  return 0;
}

//Prohibit a pair between two nucleotides in subsequent structure prediction.
//Returns an integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified, 9 = nucleotide in conflicting restraint).
int RNA::ForceProhibitPair(const int i, const int j) {
  int index, locali, localj;
  //First perform the error checking:
  //Note: In structure, forced pairs run between index of 0 and a maximum of maxforce-1.
  //check to make sure that a sequence has been read
  if (ct->GetSequenceLength() == 0) return 20;
  if (i < 1 || i > ct->GetSequenceLength()) return 4;//i is out of range
  if (j < 1 || j > ct->GetSequenceLength()) return 4;//j is out of range
  //sort indexes from 5' to 3':
  if (i > j) {
    locali = j;
    localj = i;
  }
  else {
    locali = i;
    localj = j;
  }
  //check to make sure this pair hasn't been forced:
  for (index = 0; index < ct->GetNumberofPairs(); index++) {
    if (locali == ct->GetPair5(index) && localj == ct->GetPair3(index)) return 9;//i or j is in a forced pair
  }
  //Now register the restraint because the error checking was clear or errors.
  ct->AddForbiddenPair(locali, localj);
  return 0;
}

//Force a nucleotide to be single stranded in subsequent structure prediction.
//An integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified, 9 = same nucleotide in conflicting restraint).
int RNA::ForceSingleStranded(const int i) {
  int index;
  //check to make sure that a sequence has been read
  if (ct->GetSequenceLength() == 0) return 20;
  //Check that nucleotide is valid
  if (i < 1 || i > ct->GetSequenceLength()) return 4;//i is out of range
  //Check for conflicting constraints; anything forcing a nucleotide to be paired.
  for (index = 0; index < ct->GetNumberofPairs(); index++) {//check all the forced pairs
    if (i == ct->GetPair5(index) || i == ct->GetPair3(index)) return 9;//i is in a forced pair
  }
  for (index = 0; index < ct->GetNumberofDoubles(); index++) {//check all the force doubles
    if (ct->GetDouble(index) == i) return 9;
  }
  for (index = 0; index < ct->GetNumberofGU(); index++) {//check all the force FMN
    if (ct->GetGUpair(index) == i) return 9;
  }
  //Register the constraint:
  ct->AddSingle(i);
  return 0;
}

//Return a nucleotide that is forced double stranded.
int RNA::GetForcedDoubleStranded(const int constraintnumber) {
  //First make sure the constraintnumber is valid.
  if (constraintnumber < 0 || constraintnumber >= ct->GetNumberofDoubles()) return 0;
  //Now return the constraint.
  return ct->GetDouble(constraintnumber);
}

//Return a nucleotide that is accessible to FMN cleavage.
int RNA::GetForcedFMNCleavage(const int constraintnumber) {
  //First make sure the constraintnumber is valid.
  if (constraintnumber < 0 || constraintnumber >= ct->GetNumberofGU()) return 0;
  //Now return the constraint.
  return ct->GetGUpair(constraintnumber);
}

//Return a nucleotide that is accessible to modification.
int RNA::GetForcedModification(const int constraintnumber) {
  //First make sure the constraintnumber is valid.
  if (constraintnumber < 0 || constraintnumber >= ct->GetNumberofModified()) return 0;
  //Now return the constraint.
  return ct->GetModified(constraintnumber);//note that the underlying ct indexes from 1 to ndbl.
}

//Return a nucleotide in a forced pair.
//fiveprime determines if the nucleotide is the five prime or the three prime nucleotide in the constraint.  true = five prime nucleotide.
int RNA::GetForcedPair(const int constraintnumber, const bool fiveprime) {
  //First make sure the constraintnumber is valid.
  if (constraintnumber < 0 || constraintnumber >= ct->GetNumberofPairs()) return 0;
  //Now return the constraint.
  if (fiveprime) return ct->GetPair5(constraintnumber);
  else return ct->GetPair3(constraintnumber);
}

//Return a nucleotide in a prohibited pair.
//fiveprime determines if the nucleotide is the five prime or the three prime nucleotide in the constraint.  true = five prime nucleotide.
int RNA::GetForcedProhibitedPair(const int constraintnumber, const bool fiveprime) {
  //First make sure the constraintnumber is valid.
  if (constraintnumber < 0 || constraintnumber >= ct->GetNumberofForbiddenPairs()) return 0;
  //Now return the constraint.
  if (fiveprime) return ct->GetForbiddenPair5(constraintnumber);
  else return ct->GetForbiddenPair3(constraintnumber);
}

//Return a nucleotide that is forced single stranded.
int RNA::GetForcedSingleStranded(const int constraintnumber) {
  //First make sure the constraintnumber is valid.
  if (constraintnumber < 0 || constraintnumber >= ct->GetNumberofSingles()) return 0;
  //Now return the constraint.
  return ct->GetSingle(constraintnumber);//note that the underlying ct indexes from 1 to ndbl.
}

//Return the maximum pairing distance.
//Return an integer that indicates the maximum distance allowed between paired nucleotides, where -1 indicates that the maximum distance is not set.
int RNA::GetMaximumPairingDistance() {
  if (ct->DistanceLimited()) return ct->GetPairingDistanceLimit();
  else return -1;
}

//Return the number of nucletides forced to be paired.
int RNA::GetNumberOfForcedDoubleStranded() {
  return ct->GetNumberofDoubles();
}

// Add an experimental bonus to a pair of nucleotides
//void RNA::SetExperimentalBonus(const int i, const int j, const double bonus){
//	ct->EX[i][j] = bonus;
//}
//!Return the number of nucleotides accessible to FMN cleavage.
int RNA::GetNumberOfForcedFMNCleavages() {
  return ct->GetNumberofGU();
}

//!Return the number of nucleotides accessible to chemical modification.
int RNA::GetNumberOfForcedModifications() {
  return ct->GetNumberofModified();
}

//!Return the number of forced base pairs.
int RNA::GetNumberOfForcedPairs() {
  return ct->GetNumberofPairs();
}

//!Return the number of prohibited base pairs.
int RNA::GetNumberOfForcedProhibitedPairs() {
  return ct->GetNumberofForbiddenPairs();
}

//!Return the number of nucleotides that are not allowed to pair.
int RNA::GetNumberOfForcedSingleStranded() {
  return ct->GetNumberofSingles();
}

//Read a set of folding constraints to disk in a plain text file.
//filename is a c string that is the file name to be read.
//Returns an integer that indicates an error code (0 = no error, 1 = file not found, 13 = error reading constraint file).
int RNA::ReadConstraints(const char filename[]) {
  FILE* check;
  //check that the file exists.
  if ((check = fopen(filename, "r")) == NULL) {
    //the file is not found
    fclose(check);
    return 1;
  }
  fclose(check);
  //Now read the constraints
  if (readconstraints(filename, ct)) return 0;
  else return 13;
}

//Read SHAPE data to constrain structure prediction on subsequent structure predictions.
//filename is a c string that indicates a file that contains SHAPE data.
//IsPseudoEnergy indicates whether this is the pseudo folding free energy constraint (the preferred method).  This defaults to true.
//parameter1 is the slope when IsPseudoEnergy=true and is a threshold above which nucleotides are forced single stranded otherwise.
//parameter2 is the intercept when IsPseudoEnergy=true and is a threshold above which a nucleotide is considered chemically modified otherwise.
//modifier is the type of chemical modification probe that was used (currently accepted values are SHAPE, diffSHAPE, DMS, and CMCT). Defaults to SHAPE.
//Returns an integer that indicates an error code (0 = no error, 1 = input file not found).
int RNA::ReadSHAPE(const char filename[], const double parameter1, const double parameter2, std::string modifier,
                   const bool IsPseudoEnergy) {
  FILE* check;
  //check that the SHAPE input file exists
  if ((check = fopen(filename, "r")) == NULL) {
    //the file is not found
    fclose(check);
    return 1;
  }
  fclose(check);
  if (IsPseudoEnergy) {
    //This is the pseudo energy version
    ct->SHAPEslope = parameter1 * conversionfactor;//register the slope in tenths of kcal/mol
    ct->SHAPEintercept = parameter2 * conversionfactor;//register the intercept in tenths of a kcal/mol
    ct->ReadSHAPE(filename, modifier);//call ReadSHAPE() to read the file and determine pseudo energies
  }
  else {
    ct->ReadSHAPE(filename, (float) parameter1,
                  (float) parameter2);//call ReadSHAPE() with parameters to parse thresholds
  }
  return 0;
}

//Read SHAPE data to constrain structure prediction on subsequent structure predictions.
//filename is a c string that indicates a file that contains SHAPE data.
//parameter1 is the slope.
//parameter2 is the intercept.
//Returns an integer that indicates an error code (0 = no error, 1 = input file not found).
int RNA::ReadExperimentalPairBonus(const char filename[], double const experimentalOffset,
                                   double const experimentalScaling) {
  FILE* check;
  //check that the SHAPE input file exists
  if (strlen(filename) > 0) {
    if ((check = fopen(filename, "r")) == NULL) {
      //the file is not found
      fclose(check);
      return 1;
    }
    fclose(check);
  }
  ct->ReadExperimentalPairBonus(filename, experimentalOffset, experimentalScaling);
  return 0;
}

//Read SHAPE data to constrain structure prediction on subsequent structure predictions - overloaded version for including single-stranded SHAPE.
//filename is a c string that indicates a file that contains SHAPE data.
//parameter1 is the double-stranded slope.
//parameter2 is the double-stranded intercept.
//modifier is the type of chemical modification probe that was used (currently accepted values are SHAPE, DMS, and CMCT). Defaults to SHAPE.
//ssm is the single-stranded slope.
//ssb in the single-stranded intercept.
//Returns an integer that indicates an error code (0 = no error, 1 = input file not found).
int RNA::ReadSHAPE(const char filename[], const double parameter1, const double parameter2, const double ssm,
                   const double ssb, std::string modifier) {
  FILE* check;
  //check that the SHAPE input file exists
  if ((check = fopen(filename, "r")) == NULL) {
    //the file is not found
    fclose(check);
    return 1;
  }
  fclose(check);
  ct->SHAPEslope = parameter1 * conversionfactor;//register the slope in tenths of kcal/mol
  ct->SHAPEintercept = parameter2 * conversionfactor;//register the intercept in tenths of a kcal/mol
  ct->SHAPEslope_ss = ssm * conversionfactor;//register the slope in tenths of kcal/mol
  ct->SHAPEintercept_ss = ssb * conversionfactor;//register the intercept in tenths of a kcal/mol
  ct->ReadSHAPE(filename, modifier);//call ReadSHAPE() to read the file and determine pseudo energies
  return 0;
}

//Read Double Strand Offset
int RNA::ReadDSO(const char filename[]) {
  FILE* check;
  //check that the SHAPE input file exists
  if ((check = fopen(filename, "r")) == NULL) {
    //the file is not found
    fclose(check);
    return 1;
  }
  ct->ReadOffset(NULL, filename);
  return 0;
}

//Read Single Strand Offset
int RNA::ReadSSO(const char filename[]) {
  FILE* check;
  //check that the SHAPE input file exists
  if ((check = fopen(filename, "r")) == NULL) {
    //the file is not found
    fclose(check);
    return 1;
  }
  ct->ReadOffset(filename, NULL);
  return 0;
}

//Remove all previously defined constraints.
void RNA::RemoveConstraints() {
  ct->RemoveConstraints();
  ct->min_gu = 0;
  ct->min_g_or_u = 0;
  ct->nneighbors = 0;
  ct->nregion = 0;
  ct->nmicroarray = 0;
}

//Add extrinsic restraints for partition function calculations.
int RNA::SetExtrinsic(int i, int j, double k) {
  int locali, localj;
  //First do the error checking:
  //check the indices
  if (i < 1 || i > ct->GetSequenceLength() || j < 1 || j > ct->GetSequenceLength()) return 4;
  //make sure the equilibrium constant is not less than zero.
  if (k < 0) return 26;
  //Now past error checking:
  //sort indexes from 5' to 3':
  if (i > j) {
    locali = j;
    localj = i;
  }
  else {
    locali = i;
    localj = j;
  }
  if (ct->constant == NULL) {
    //allocate the space needed in structure to store the constants
    ct->allocateconstant();
  }
  ct->constant[localj][locali] = k;
  return 0;
}

//Write the set of folding constraints to disk.
int RNA::WriteConstraints(const char filename[]) {
  outputconstraints(filename, ct);
  return 0;
}

//Specify a comment for inclusion in subsequently written .ct files.
int RNA::AddComment(const char comment[], const int structurenumber) {
  string label;
  //start with error checking:
  if (structurenumber < 1 || structurenumber > ct->GetNumberofStructures()) return 3;
  //now register the comment (at the end of existing string, but before the newline there by default):
  label = ct->GetCtLabel(structurenumber);
  //Remove the existing newline character, if it exists:
  if (label.length() > 0) if (label[label.length() - 1] == '\n') label.erase(label.length() - 1);
  //Add the comment
  label += comment;
  //Add back the newline
  label += "\n";
  //add a newline at the end of the comment -- required when ct is written
  ct->SetCtLabel(label, structurenumber);
  return 0;
}

//Write a ct file of the structures
int RNA::WriteCt(const char filename[], bool append) {
  if (ct->GetNumberofStructures() > 0) {
    ct->ctout(filename, append);
    return 0;
  }
  else return 10; //an error code
}

//Write a ct file of the structures
int RNA::WriteDotBracket(const char filename[]) {
  if (ct->GetNumberofStructures() > 0) {
    ct->writedotbracket(filename);
    return 0;
  }
  else return 10; //an error code
}

// Report if there are any pseudoknots in a structure.
bool RNA::ContainsPseudoknot(const int structurenumber) {
  int i, j;
  //make sure structurenumber is a valid structure
  if (structurenumber < 1 || structurenumber > ct->GetNumberofStructures()) {
    ErrorCode = 3;
    return false;
  }
  else {
    //passed error trapping:
    //check all nucs for pairs
    for (i = 1; i < ct->GetSequenceLength(); i++) {
      if (ct->GetPair(i, structurenumber) > i) {
        //found pair, check for crossing pairs
        for (j = i + 1; j < ct->GetPair(i, structurenumber); j++) {
          if (ct->GetPair(j, structurenumber) > j) {
            if (ct->GetPair(j, structurenumber) > ct->GetPair(i, structurenumber)) {
              return true;
            }
          }
        }
      }
    }
    return false;//no pseudoknot was found
  }
}

//Get the folding free energy change for a predicted structure.
double RNA::GetFreeEnergy(const int structurenumber) {
  //make sure structurenumber is a valid structure
  if (structurenumber < 1 || structurenumber > ct->GetNumberofStructures()) {
    ErrorCode = 3;
    return 0.0;
  }
  else {
    //error trapping complete
    return (((double) ct->GetEnergy(structurenumber)) / ((double) conversionfactor));
  }
}

// Get the nucleotide to which the specified nucleotide is paired.
int RNA::GetPair(const int i, const int structurenumber) {
  //check to make sure i is a valid nucleotide
  if (i < 1 || i > ct->GetSequenceLength()) {
    ErrorCode = 4;
    return 0;
  }
    //make sure there is structure information for this RNA
  else if (ct->GetNumberofStructures() == 0) {
    ErrorCode = 23;
    return 0;
  }
    //now make sure structurenumber is a valid structure
  else if (structurenumber < 1 || structurenumber > ct->GetNumberofStructures()) {
    ErrorCode = 3;
    return 0;
  }
  else {
    //error trapping complete
    return ct->GetPair(i, structurenumber);
  }
}

//Get the total number of specified structures
int RNA::GetStructureNumber() {
  return ct->GetNumberofStructures();
}

//Provide the comment from the ct file as a string.
std::string RNA::GetCommentString(const int structurenumber) {
  std::string temp;
  //Add some code for backwards compatibility:
  //In the past, all labels were associated with structures.  Now, labels can come with the sequence.
  //If there are no structures, return the sequence label:
  if (ct->GetNumberofStructures() == 0) {
    temp = ct->GetSequenceLabel();
    return temp;
  }
  //start with error checking:
  if (structurenumber < 1 || structurenumber > ct->GetNumberofStructures()) {
    //The request is for a structure that is out of range
    ErrorCode = 3;
    temp = "";
    return temp;
  }
  temp = ct->GetCtLabel(structurenumber);
  return temp;
}

//Get the identity of nucleotide i.
char RNA::GetNucleotide(const int i) {
  //check to make sure that a sequence has been read
  if (ct->GetSequenceLength() == 0) {
    ErrorCode = 20;
    return '-';
  }
  //Check that nucleotide is valid
  if (i < 1 || i > ct->GetSequenceLength()) {
    ErrorCode = 4;//i is out of range
    return '-';
  }
  return ct->nucs[i];
}

//Get the total length of the sequence
int RNA::GetSequenceLength() {
  return ct->GetSequenceLength();
}

//Return the type of backbone (true = RNA, false = DNA).
bool RNA::GetBackboneType() {
  return isrna;
}

//Access the underlying structure class.
structure* RNA::GetStructure() {
  return ct;
}

//Provide a TProgressDialog for following calculation progress.
//A TProgressDialog class has a public function void update(int percent) that indicates the progress of a long calculation.
void RNA::SetProgress(TProgressDialog& Progress) {
  progress = &Progress;
  return;
}

//Provide a means to stop using a TProgressDialog.
//StopProgress tells the RNA class to no longer follow progress.  This should be called if the TProgressDialog is deleted, so that this class does not make reference to it.
void RNA::StopProgress() {
  progress = NULL;
  return;
}

TProgressDialog* RNA::GetProgress() {
  return progress;
}

RNA::~RNA() {
  delete ct;//delete the structure
}

//This is a protected function for handling file input.
int RNA::FileReader(const char filename[], const int type) {
  FILE* check;
  short vers;
  //double scaling;
  int i;
  const char* errMsg; //used in catch blocks
  // RMW 2015-03-12: try/catch to prevent any unexpected errors from crashing the program. At least show an error message.
  // (For example, I found that passing the wrong type of file to openct caused an unhandled memory allocation exception.)
  try {
    if ((check = fopen(filename, "r")) == NULL) {
      //the file is not found
      return 1;
    }
    else {
      fclose(check);
      //open the file based on type:
      if (type == 1) {
        //type indicates a ct file
        long linenumber = ct->openct(filename);
        //if linenumber==0, there was no error
        if (linenumber == 0) return 0;
          //Otherwise, there was an error opening the file
        else return 2;
      }
      else if (type == 2) {
        //type indicates a .seq file
        //ct->numofstructures=0;
        if (ct->openseq(filename) == 1) return 0;//no error
        else return 2;//File open error
      }
      return 22;
    }
  } catch (std::exception* ex) {
    errMsg = ex->what();
  } catch (int err) {
    char result[20];
    sprintf(result, "Error number %d", err);
    errMsg = result;
  } catch (...) {
    errMsg = "(unknown)";
  }
  std::cerr << "Unexpected file read exception: " << errMsg << endl;
  return 2; //we only reach this line if there was an exception.
}
//The following should not be included for compilations for Windows:
#ifndef _WINDOWS_GUI

//A global function for error reporting
void errmsg(int err, int erri) {
  if (err == 30) {
    std::cout << "End Reached at traceback #" << erri << "\n";
    return;
  }
  if (err == 100) {
    std::cout << "error # " << erri;
    return;
  }
  switch (err) {
    case 1:
      std::cout << "Could not allocate enough memory";
      break;
    case 2:
      std::cout << "Too many possible base pairs";
      break;
    case 3:
      std::cout << "Too many helixes in multibranch loop";
    case 4:
      std::cout << "Too many structures in CT file";
    default:
      std::cout << "Unknown error";
  }
//std::cin >> err;
  return;
}

#endif
