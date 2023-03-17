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
#include <random>

#include <modern_robotics.h>
#include <Indy7.h>
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
    // Initial Setup
    q<<0,0,-M_PI/2.0,0,-M_PI/2.0,0;
    q_start = q;
    q_end<<1,1,1,1,1,1;
    g<< 0, 0, -9.8;
    double scale = 10.0;
    Xstart = FKinSpace(indy7->M, indy7->Slist,q );
    SO3 Rend = TransToR(Xstart);
    Vector3d p_temp;
    p_temp << 0.1, 0.1, 0.1;
		Vector3d pend = TransToP(Xstart)+p_temp;
    Xend = RpToTrans(Rend, pend);
    for(int i=0; i<n; ++i)
    {
      switch(i)
      {
        case 0:
        case 1:
          Kp(i,i) = 70*scale;
          Kv(i,i) = 55*scale;
          break;
        case 2:
          Kp(i,i) = 40*scale;
          Kv(i,i) = 30*scale;
          break;
        case 3:
        case 4:
          Kp(i,i) = 25*scale;
          Kv(i,i) = 15*scale;
          break;
        case 5:
          Kp(i,i) = 18*scale;
          Kv(i,i) = 3*scale;
          break;
      }
    }    
    //

    joint_state_publisher_ = create_publisher<sensor_msgs::msg::JointState>("joint_states", 10);
    task_error_publisher_ = create_publisher<sensor_msgs::msg::JointState>("task_err", 10);
    joint_state_timer_ = create_wall_timer(std::chrono::microseconds(1000), std::bind(&JointStatePublisherNode::control, this));
    goal_pose_subscriber_ = create_subscription<geometry_msgs::msg::PoseStamped>(
        "goal_pose", 10, std::bind(&JointStatePublisherNode::goal_pose_callback, this, std::placeholders::_1));    
    pose_publisher_ = create_publisher<geometry_msgs::msg::PoseStamped>("new_pose", 10);
    path_publisher_ = create_publisher<nav_msgs::msg::Path>("path", 10); // Path 메시지 Publisher 생성

  }

private:
  rclcpp::TimerBase::SharedPtr joint_state_timer_;
  
  rclcpp::Publisher<sensor_msgs::msg::JointState>::SharedPtr joint_state_publisher_;
  rclcpp::Publisher<sensor_msgs::msg::JointState>::SharedPtr task_error_publisher_;
  rclcpp::Subscription<std_msgs::msg::Bool>::SharedPtr bool_subscriber_;
  rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr goal_pose_subscriber_;
  rclcpp::Publisher<geometry_msgs::msg::PoseStamped>::SharedPtr pose_publisher_;
  rclcpp::Publisher<nav_msgs::msg::Path>::SharedPtr path_publisher_; 

  //Indy Setup
  Indy7 *indy7 = new Indy7("/home/sung/ros2_humble/src/mr_dual_arm_control/MR_info.json");
  static unsigned int print_count;
  double t = 0;
  int CONTROL_RATE = 1000;
  int DISPLAY_PATH_NUM = 50;
  int n = 6;
  double Tf = 5;
  double dt = 1.0/CONTROL_RATE;
  JVec q=JVec::Zero();
  JVec q_start =JVec::Zero();
  JVec q_end =JVec::Zero();
	JVec dq=JVec::Zero();
	JVec ddq=JVec::Zero();
	Vector6d Ftip = Vector6d::Zero();
	Vector3d g= Vector3d::Zero();
	Eigen::Matrix<double, 6, 6> Kp =Matrix6d::Identity();
  Eigen::Matrix<double, 6, 6> Kv =Matrix6d::Identity();
  JVec q_desired=JVec::Zero();
  JVec dq_desired=JVec::Zero();
  JVec ddq_desired=JVec::Zero();
  JVec taulist = JVec::Zero();  
  SE3 Xstart = SE3::Identity();
  SE3 Xend = SE3::Identity();
  Vector6d Vstart = Vector6d::Zero();
  Vector6d Vend = Vector6d::Zero();
  Vector6d Vdotstart = Vector6d::Zero();
  Vector6d Vdotend = Vector6d::Zero();
  bool flag = 0;
  void goal_pose_callback(const geometry_msgs::msg::PoseStamped::SharedPtr msg)
  {
    cout<<"NEW POSE"<<endl;
    //Random  Params Setting
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_z(0.2, 0.7);
    std::uniform_real_distribution<> dis_r(-M_PI/8.0, M_PI/8.0);
    std::uniform_real_distribution<> dis_p(-M_PI/8.0, M_PI/8.0);
    std::uniform_real_distribution<> dis_y(-M_PI/8.0, M_PI/8.0);

    bool has_nan1 = q.array().isNaN().any();
    bool has_nan2 = dq.array().isNaN().any();
    bool has_nan3 = ddq.array().isNaN().any();
    if(has_nan1 ||  has_nan2 || has_nan3 ){
       q<<0,0,-M_PI/2.0,0,-M_PI/2.0,0;
       dq<<0,0,0,0,0,0;
       ddq<<0,0,0,0,0,0;
       Vstart<<0,0,0,0, 0,0;
       Vdotstart<<0,0,0,0,0,0;
    }    
    Jacobian Jb = JacobianBody(indy7->Blist, q);
    Jacobian Ja = AnalyticJacobianBody(indy7->M,indy7->Blist, q);
    Jacobian dJa= dAnalyticJacobianBody(indy7->M,indy7->Blist, q ,dq);

    Vstart = Ja*dq;
    Vdotstart = Ja*ddq + dJa*dq;
    cout<<"q : "<<q.transpose()<<endl;
    cout<<"dq : "<<dq.transpose()<<endl;
    cout<<"ddq : "<<ddq.transpose()<<endl;
    cout<<"Vstart : "<<Vstart.transpose()<<endl;
    cout<<"Vend : "<<Vend.transpose()<<endl;
    cout<<"Vdotstart : "<<Vdotstart.transpose()<<endl;
    cout<<"Vdotend : "<<Vdotend.transpose()<<endl;
    //GET goal_pose -> position and Rotation matrix
    Eigen::Vector3d position(msg->pose.position.x, msg->pose.position.y, msg->pose.position.z);
    Eigen::Quaterniond orientation(msg->pose.orientation.w, msg->pose.orientation.x, msg->pose.orientation.y, msg->pose.orientation.z);
    
    SO3 Rend = orientation.toRotationMatrix();      
    // Roffset = FKinspace(M,Slist,q0)
    SO3 Roffset;
    Roffset <<-1,0, 0,
               0,1,0,
              0,0,-1;
    //RESETTING TRAJECTORY 
    flag= 1;
    t= 0;
    Xstart = FKinSpace(indy7->M, indy7->Slist,q );
    Vector3d pend;
    pend << msg->pose.position.x, msg->pose.position.y, dis_z(gen);
    SO3 Rx;
    double r = dis_r(gen);
    double p = dis_p(gen);
    double y = dis_y(gen);
    Rx << 1,0,0,
          0,cos(r),-sin(r),
          0,sin(r),cos(r);
    SO3 Ry;
    Ry << cos(p),0,sin(p),
          0,1,0,
          -sin(p),0,cos(p);
    SO3 Rz;
    Rz << cos(y),-sin(y),0,
          sin(y),cos(y),0,
          0,0,1;
    Xend = RpToTrans(Rend*Roffset*Rz*Ry*Rx, pend);
    
    // new_pose publish
    SO3 new_R = Rend*Roffset*Rz*Ry*Rx;
    Quaterniond q(new_R);
    auto new_pose_msg = std::make_shared<geometry_msgs::msg::PoseStamped>();
    new_pose_msg->header.stamp = this->now();
    new_pose_msg->header.frame_id = "world";
    new_pose_msg->pose.position.x = pend[0];
    new_pose_msg->pose.position.y = pend[1];
    new_pose_msg->pose.position.z = pend[2];
    new_pose_msg->pose.orientation.w = q.w();
    new_pose_msg->pose.orientation.x = q.x();
    new_pose_msg->pose.orientation.y = q.y();
    new_pose_msg->pose.orientation.z = q.z();
    pose_publisher_->publish(*new_pose_msg);
    
    //Path Publish
    auto path_msg  = std::make_shared<nav_msgs::msg::Path>();
    path_msg->header.frame_id = "world";
    double timegap = Tf/DISPLAY_PATH_NUM/1.0;
    for(int i =0;i<DISPLAY_PATH_NUM;i++){
          double t_ = timegap*i;
          SE3 Xd= EulerT(Xstart, Xend, Vstart,Vend,Vdotstart,Vdotend, t_, Tf) ;
          SO3 Rd = TransToR(Xd);
          Vector3d pd = TransToP(Xd);   
          Quaterniond qd(Rd);
          auto temp_pose_msg = std::make_shared<geometry_msgs::msg::PoseStamped>();
          temp_pose_msg->pose.position.x = pd[0];
          temp_pose_msg->pose.position.y = pd[1];
          temp_pose_msg->pose.position.z = pd[2];
          temp_pose_msg->pose.orientation.w = qd.w();
          temp_pose_msg->pose.orientation.x = qd.x();
          temp_pose_msg->pose.orientation.y = qd.y();
          temp_pose_msg->pose.orientation.z = qd.z();
          path_msg->poses.push_back(*temp_pose_msg);
    }
    path_publisher_->publish(*path_msg);
    cout<<"PATH PUBLISH"<<endl;

  }

  void control()
  {
    if(flag==0)return;
    //SE3 Xd= CartesianT(Xstart, Xend, t, Tf, 5) ;
    SE3 Xd= EulerT(Xstart, Xend, Vstart,Vend,Vdotstart,Vdotend, t, Tf) ;
     SO3 Rd = TransToR(Xd);
     Vector3d pd = TransToP(Xd);
     Vector6d Vd= CartesianVel(Xstart, Xend, t, Tf, 5) ;
     Vector6d dVd= CartesianAcc(Xstart, Xend, t, Tf, 5) ;
    double lambda = 0.01;
     Jacobian Jb = JacobianBody(indy7->Blist, q);
     Jacobian Ja = AnalyticJacobianBody(indy7->M,indy7->Blist, q);
     Matrix6d Identity6d = Matrix6d::Identity();
     pinvJacobian Jb_pinv = Jb.transpose() * (Jb * Jb.transpose() + lambda * lambda *Identity6d ).inverse();
     pinvJacobian Ja_pinv = Ja.transpose() * (Ja* Ja.transpose() + lambda * lambda *Identity6d ).inverse();
    Jacobian dJb= dJacobianBody(indy7->Blist, q ,dq);
		Jacobian dJa= dAnalyticJacobianBody(indy7->M,indy7->Blist, q ,dq);

     SE3 Tsb = FKinBody(indy7->M,indy7->Blist,q);
     SO3 Rsb = TransToR(Tsb);
     Vector3d psb = TransToP(Tsb);
     Vector3d rot_err = so3ToVec(MatrixLog3(Rsb.transpose()*Rd));
     Vector3d pos_err = pd-psb;
     Vector6d Xe;
     Xe<< rot_err[0],rot_err[1],rot_err[2],pos_err[0],pos_err[1],pos_err[2];

     
     dq= Ja_pinv*Vd+1000.0*Ja_pinv*Xe;
	   ddq = Ja_pinv*dVd-Ja_pinv*dJa*dq;


  /*
    // INVERSE DYNAMICS 
      JVec taugrav = InverseDynamics(q,dq,ddq,g,Ftip,indy7->Mlist,indy7->Glist,indy7->Slist);
      MassMat Mmat = MassMatrix(q, indy7->Mlist,indy7->Glist,indy7->Slist);
      JVec C =  VelQuadraticForces(q,dq,indy7->Mlist,indy7->Glist,indy7->Slist);
      JVec qddot_ref = ddq_desired+Kv*(dq_desired-dq) + Kp*(q_desired- q);
      JVec tauCTM  =Mmat*qddot_ref + taugrav;
      taulist = tauCTM;
    // FOWRAD DYNAMIC SIMULATION
      JVec ddq = ForwardDynamics(q,dq, taulist,
                    g, Ftip, indy7->Mlist,
                    indy7->Glist, indy7->Slist);
  */                  
		EulerStep(q, dq, ddq, dt);
    bool has_nan = q.array().isNaN().any();
    if(has_nan){
       q<<0,0,-M_PI/2.0,0,-M_PI/2.0,0;
       Vstart<<0,0,0,0, 0,0;
       Vdotstart<<0,0,0,0,0,0;
       path_publisher_ = create_publisher<nav_msgs::msg::Path>("path", 10); 
    }
    JVec wrap_q =wrapToPI(q);
    // PUBLISH
    auto joint_state_msg = std::make_shared<sensor_msgs::msg::JointState>();
    joint_state_msg->header.stamp = this->now();
    joint_state_msg->name = {"l_joint_0", "l_joint_1","l_joint_2","l_joint_3","l_joint_4","l_joint_5","r_joint_0", "r_joint_1","r_joint_2","r_joint_3","r_joint_4","r_joint_5"};
    joint_state_msg->position = {sin(2*3.141592*t/10.0),sin(2*3.141592*t/10.0),sin(2*3.141592*t/10.0),sin(2*3.141592*t/10.0),sin(2*3.141592*t/10.0),sin(2*3.141592*t/10.0)
                                ,sin(2*3.141592*t/10.0),sin(2*3.141592*t/10.0),sin(2*3.141592*t/10.0),sin(2*3.141592*t/10.0),sin(2*3.141592*t/10.0),sin(2*3.141592*t/10.0)};

    if(++print_count>=15){
        auto task_err_msg = std::make_shared<sensor_msgs::msg::JointState>();
        task_err_msg->header.stamp = this->now();
        task_err_msg->name = {"wx", "wy","wz","vx","vy","vz"};
        task_err_msg->position = {Xe[0],Xe[1],Xe[2],Xe[3],Xe[4],Xe[5]};
        task_error_publisher_->publish(*task_err_msg);
        joint_state_publisher_->publish(*joint_state_msg);        
      // cout<<"t : "<< t << " - dq : "<<(dq).transpose()<<endl;
      print_count = 0;
    }
    t=  t+dt;
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
