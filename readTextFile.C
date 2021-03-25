#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include "TLorentzVector.h"
#include <fstream>
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"
#include <fstream>
#include <typeinfo>
#include <sys/stat.h> 
#include <sys/types.h>
#include <sstream>

using namespace std;

void readTextFile(){
    std::string inputTextFile = "run_01_tag_1_banner_0.151356GeV_beam2.9GeV_v8.txt";
    ifstream myfile;
    myfile.open(inputTextFile, ios::in);
    std::string str; 
    while (std::getline(myfile, str))
    {   
        if (str.find("Integrated weight") != std::string::npos){
            std::cout << str << std::endl;
            istringstream ss(str);
 
            string hash, word1, word2, pb, colon, number; // for storing each word
            ss >> hash >> word1 >> word2 >> pb >> colon >> number;
            std::cout << number << std::endl;
        }
    }
}
