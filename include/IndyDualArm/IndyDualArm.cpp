#include "IndyDualArm.h"

/* STEP PC
#include "../../include/ModernRobotics/ModernRobotics.h"
#include "../../include/json/json/json.h"
#pragma comment(lib, "jsoncpp.lib")
#include "../eigen/Eigen/Dense"
*/


/* X86 */
#include <json/json.h>
#pragma comment(lib, "jsoncpp.lib")
#include <ModernRobotics.h>

#include <RelativeModernRobotics.h>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
using namespace std;
using namespace Eigen;
IndyDualArm::IndyDualArm(const std::string mr_json_path){
	cout<<"INDYDUALARM LOAD"<<endl;
	this->actuated_joint_num = 12;
	this->MRSetup(mr_json_path);
	HOMEPOS_r << -0.3129 ,-1.1004, -1.2474, -0.7510, -0.7714 , 0.0695;
	HOMEPOS_l << 0.3129 , 1.1004,  1.2474 , 0.7510,  0.7714 ,-0.0695;
	HOMEPOS_rel= get_q_rel(HOMEPOS_r, HOMEPOS_l );

}

bool ReadFromFile_(const std::string filename, char* buffer, int len){
  FILE* r = fopen(filename.c_str(),"rb");
  if (NULL == r)
       return false;
  size_t fileSize = fread(buffer, 1, len, r);
  fclose(r);
  return true;

}
bool ReadMRData_(const std::string filename,Json::Value &rootr){
	//cout<<"START ReadMRData_"<<endl;
	const int BufferLength = 102400;
	char readBuffer[BufferLength] = {0,};
	if (false == ReadFromFile_(filename, readBuffer, BufferLength)) {
		cout<<"Failed"<<endl;
		return -1;
	}
	std::string config_doc = readBuffer;

	Json::Reader reader;
	bool parsingSuccessful = reader.parse(config_doc,rootr);
	if ( !parsingSuccessful ) { 
		std::cout << "Failed to parse configuration\n" << reader.getFormatedErrorMessages(); 
		return -1;
		
	}
}


void IndyDualArm::MRSetup(const std::string mr_json_path){
	//cout<<"START MRSetup"<<endl;
	Json::Value rootr;
	bool ret = ReadMRData_(mr_json_path,rootr);

//	cout<<"MR Setup 1"<<endl;
	for(int i =0;i<6 ; i++){
		for(int j =0;j<this->actuated_joint_num;j++){
			this->relSlist(i,j) = rootr["relSlist"][i][j].asDouble();
			this->relBlist(i,j) = rootr["relBlist"][i][j].asDouble();
		}
	}
	//cout<<"MR Setup 2"<<endl;
	for(int i =0;i<6 ; i++){
		for(int j =0;j<this->actuated_joint_num/2;j++){
			this->Slist(i,j) = rootr["Slist"][i][j].asDouble();
			this->Blist(i,j) = rootr["Blist"][i][j].asDouble();
		}
	}	
	//cout<<"MR Setup 3"<<endl;	
//	printMatrix(this->relSlist,"relSlist");
	//printMatrix(this->relBlist,"relBlist");
	
	for(int i = 0;i< rootr["Mlist"].size(); i++){
		MatrixXd M = MatrixXd::Identity(4,4);
		for(int j = 0;j< rootr["Mlist"][0].size(); j++){
			for(int k = 0;k< rootr["Mlist"][0][0].size(); k++){
				M(j,k) = rootr["Mlist"][i][j][k].asDouble();
			}
		}
		char str[50];
		//sprintf(str,"M%d%d",i,i+1);
		
		this->Mlist.push_back(M);
		//printMatrix(M,str);
	}
	for(int i = 0;i< rootr["Glist"].size(); i++){
		MatrixXd G = MatrixXd::Identity(6,6);
		for(int j = 0;j< rootr["Glist"][0].size(); j++){
			for(int k = 0;k< rootr["Glist"][0][0].size(); k++){
				G(j,k) = rootr["Glist"][i][j][k].asDouble();
			}
		}
		char str[50];
		//sprintf(str,"G%d",i);
		
		this->Glist.push_back(G);
		//printMatrix(G,str);
	}	

	for (int i = 0;i<4;i++){
		for (int j = 0;j<4;j++){
			this->relM(i,j) = rootr["relM"][i][j].asDouble();
		}
	}
	for (int i = 0;i<4;i++){
		for (int j = 0;j<4;j++){
			this->M(i,j) = rootr["M"][i][j].asDouble();
		}
	}	
	for (int i = 0;i<4;i++){
		for (int j = 0;j<4;j++){
			this->Tbr(i,j) = rootr["Tbr"][i][j].asDouble();
		}
	}	
	for (int i = 0;i<4;i++){
		for (int j = 0;j<4;j++){
			this->Tbl(i,j) = rootr["Tbl"][i][j].asDouble();
		}
	}			
	//printMatrix(this->relM,"relM");



	
}
IndyDualArm::~IndyDualArm(){
	
}
