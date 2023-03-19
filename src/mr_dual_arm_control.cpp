#include "rclcpp/rclcpp.hpp"
#include "sensor_msgs/msg/joint_state.hpp"
#include "std_msgs/msg/bool.hpp"
#include "geometry_msgs/msg/pose_stamped.hpp"
#include "std_msgs/msg/header.hpp"
#include "nav_msgs/msg/path.hpp" 
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>      // std::setprecision
#include <ModernRobotics.h>
#include <RelativeModernRobotics.h>
#include <random>

#include <modern_robotics.h>
#include <IndyDualArm.h>
#include <chrono>
#include <iostream>

using namespace Eigen;
using namespace std;


class JointStatePublisherNode : public rclcpp::Node
{
public:
  JointStatePublisherNode()
  : Node("joint_state_publisher")
  {
    Tbr= dualarm->Tbr;
	  q_r  = dualarm->HOMEPOS_r;
	  q_l  = dualarm->HOMEPOS_l;
	  q_rel = dualarm->HOMEPOS_rel;    
    joint_state_publisher_ = create_publisher<sensor_msgs::msg::JointState>("joint_states", 10);
    joint_state_timer_ = create_wall_timer(std::chrono::microseconds(1000), std::bind(&JointStatePublisherNode::control, this));
  }

private:
  rclcpp::TimerBase::SharedPtr joint_state_timer_;
  rclcpp::Publisher<sensor_msgs::msg::JointState>::SharedPtr joint_state_publisher_;

  //DualArm Setup
  IndyDualArm *dualarm = new IndyDualArm("/home/sung/ros2_humble/src/mr_dual_arm_control/MR_dualarm_info.json");

	mr::JVec q_r=mr::JVec::Zero();
	mr::JVec dq_r=mr::JVec::Zero();
	mr::JVec ddq_r=mr::JVec::Zero();
	mr::JVec q_l=mr::JVec::Zero();
	mr::JVec dq_l=mr::JVec::Zero();
	mr::JVec ddq_l=mr::JVec::Zero();
	relmr::JVec q_rel = relmr::JVec::Zero();
	relmr::JVec dq_rel = relmr::JVec::Zero();
	relmr::JVec ddq_rel = relmr::JVec::Zero();

  static unsigned int print_count;
  int CONTROL_RATE = 1000;
  int DISPLAY_HZ = 30;
  double Tf = 5;
  double dt = 1.0/CONTROL_RATE;
  double t = 0;
	mr::SE3 Xd = mr::SE3::Identity();
  mr::SE3 Tbr =mr::SE3::Identity();
  void control()
  {
    Xd<<-1,0,0,0.1*sin(2*M_PI*t/2.0),
        0,1,0,0.1*cos(2*M_PI*t/2.0),
        0,0,-1,0.1,
        0,0,0,1;
    mr::SE3 Treef = mr::FKinSpace(dualarm->M,dualarm->Slist,q_r);
    Treef = Tbr*Treef;
    mr::SE3 Trel =  relmr::FKinSpace(dualarm->relM,dualarm->relSlist,q_rel);
    mr::SE3 Tleef = Treef*Trel;
    relmr::Jacobian Jb_rel=relmr::JacobianBody(dualarm->relBlist,q_rel);
    relmr::pinvJacobian pinvJb_rel = Jb_rel.transpose()*(Jb_rel*Jb_rel.transpose()).inverse();
    mr::Vector6d Xe  =mr::se3ToVec( mr::MatrixLog6(mr::TransInv(Trel)*Xd));
    dq_rel =1000.0*pinvJb_rel*Xe;

    relmr::EulerStep(q_rel,dq_rel,ddq_rel,dt);
    dualarm->get_q_r_q_l(q_rel,q_r,q_l );

    t=  t+dt;



    //ROS TOPIC
    auto joint_state_msg = std::make_shared<sensor_msgs::msg::JointState>();
    joint_state_msg->header.stamp = this->now();
    joint_state_msg->name = {"l_joint_0", "l_joint_1","l_joint_2","l_joint_3","l_joint_4","l_joint_5","r_joint_0", "r_joint_1","r_joint_2","r_joint_3","r_joint_4","r_joint_5"};
    joint_state_msg->position = {q_l[0],q_l[1],q_l[2],q_l[3],q_l[4],q_l[5],q_r[0],q_r[1],q_r[2],q_r[3],q_r[4],q_r[5]};

    if(++print_count>=10){
      joint_state_publisher_->publish(*joint_state_msg);
      print_count = 0;
    }

  }
};
unsigned int JointStatePublisherNode::print_count = 0;




int main(int argc, char *argv[])
{
  rclcpp::init(argc, argv);
  auto node = std::make_shared<JointStatePublisherNode>();
  rclcpp::WallRate loop_rate(1000);
  while (rclcpp::ok()) {
    rclcpp::spin_some(node);
    loop_rate.sleep();
  }
  rclcpp::shutdown();
  return 0;
}
