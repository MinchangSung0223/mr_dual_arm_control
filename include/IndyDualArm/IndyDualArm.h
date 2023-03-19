#ifndef INDYDUALARM_SETUP_H
#define INDYDUALARM_SETUP_H

/* STEP PC
#include "../../include/ModernRobotics/ModernRobotics.h"
#include "../eigen/Eigen/Dense"
*/

#include <RelativeModernRobotics.h>
#include <ModernRobotics.h>

#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>      // std::setprecision

using namespace Eigen;
using namespace std;
using namespace relmr;
class IndyDualArm
{
	int actuated_joint_num;
public:

	mr::ScrewList Slist;
	mr::ScrewList Blist;
	mr::SE3 M;
	
	relmr::ScrewList relSlist;
	relmr::ScrewList relBlist;
	mr::JVec HOMEPOS_r=mr::JVec::Zero();
	mr::JVec HOMEPOS_l=mr::JVec::Zero();
	relmr::JVec HOMEPOS_rel=relmr::JVec::Zero();
	SE3 relM;				
	SE3 Tbr;				
	SE3 Tbl;		

	vector<mr::Matrix6d> Glist;	
	vector<mr::SE3> Mlist;			


	IndyDualArm(const std::string mr_json_path);
	relmr::JVec get_q_rel(const mr::JVec& q_r,const mr::JVec& q_l ){
		relmr::JVec q_rel =  relmr::JVec::Zero();
		for(int i = 0;i<JOINTNUM;i++){
			q_rel[i] = -q_r[JOINTNUM-i-1];
			q_rel[i+JOINTNUM] = q_l[i];
		}
		return q_rel;
	};
	void get_q_r_q_l(const relmr::JVec& q_rel, mr::JVec& q_r,mr::JVec& q_l ){
		for(int i = 0;i<JOINTNUM;i++){
			q_r[JOINTNUM-i-1] = -q_rel[i];
			q_l[i] = q_rel[i+JOINTNUM];
		}
	};	
	void MRSetup(const std::string mr_json_path);
	void print_mr_setting(){
		cout << u8"┌─────────────────────Slist────────────────────┐\n";
		cout<< std::fixed << std::setprecision(4) <<this->Slist<<endl;
		cout << u8"└──────────────────────────────────────────────┘"<<endl;	
		cout << u8"┌─────────────────────Blist────────────────────┐\n";
		cout<< std::fixed << std::setprecision(4) <<this->Blist<<endl;
		cout << u8"└──────────────────────────────────────────────┘"<<endl;	
		cout << u8"┌───────────────M──────────────┐\n";
		cout<< std::fixed << std::setprecision(4) <<this->M<<endl;
		cout << u8"└──────────────────────────────┘"<<endl;		
		cout << u8"┌─────────────────────────────────────────RelativeSlist────────────────────────────────────────┐\n";
		cout<< std::fixed << std::setprecision(4) <<this->relSlist<<endl;
		cout << u8"└──────────────────────────────────────────────────────────────────────────────────────────────┘"<<endl;		
		cout << u8"┌─────────────────────────────────────────RelativeBlist────────────────────────────────────────┐\n";
		cout<< std::fixed << std::setprecision(4) <<this->relBlist<<endl;
		cout << u8"└──────────────────────────────────────────────────────────────────────────────────────────────┘"<<endl;			
		cout << u8"┌──────────────relM────────────┐\n";
		cout<< std::fixed << std::setprecision(4) <<this->relM<<endl;
		cout << u8"└──────────────────────────────┘"<<endl;			
		cout << u8"┌───────────────Tbr────────────┐\n";
		cout<< std::fixed << std::setprecision(4) <<this->Tbr<<endl;
		cout << u8"└──────────────────────────────┘"<<endl;			
		cout << u8"┌───────────────Tbl────────────┐\n";
		cout<< std::fixed << std::setprecision(4) <<this->Tbl<<endl;
		cout << u8"└──────────────────────────────┘"<<endl;			
	};
	virtual ~IndyDualArm();

};
#endif  //INDYDUALARM_SETUP_H
