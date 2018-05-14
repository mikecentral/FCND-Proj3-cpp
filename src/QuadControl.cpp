#include "Common.h"
#include "QuadControl.h"

#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

void QuadControl::Init()
{
  BaseController::Init();

  // variables needed for integral control
  integratedAltitudeError = 0;
    
#ifndef __PX4_NUTTX
  // Load params from simulator parameter system
  ParamsHandle config = SimpleConfig::GetInstance();
   
  // Load parameters (default to 0)
  kpPosXY = config->Get(_config+".kpPosXY", 0);
  kpPosZ = config->Get(_config + ".kpPosZ", 0);
  KiPosZ = config->Get(_config + ".KiPosZ", 0);
     
  kpVelXY = config->Get(_config + ".kpVelXY", 0);
  kpVelZ = config->Get(_config + ".kpVelZ", 0);

  kpBank = config->Get(_config + ".kpBank", 0);
  kpYaw = config->Get(_config + ".kpYaw", 0);

  kpPQR = config->Get(_config + ".kpPQR", V3F());

  maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
  maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
  maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
  maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

  maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

  minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
  maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);
#else
  // load params from PX4 parameter system
  //TODO
  param_get(param_find("MC_PITCH_P"), &Kp_bank);
  param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{
  // Convert a desired 3-axis moment and collective thrust command to 
  //   individual motor thrust commands
  // INPUTS: 
  //   desCollectiveThrust: desired collective thrust [N]
  //   desMoment: desired rotation moment about each axis [N m]
  // OUTPUT:
  //   set class member variable cmd (class variable for graphing) where
  //   cmd.desiredThrustsN[0..3]: motor commands, in [N]

  // HINTS: 
  // - you can access parts of desMoment via e.g. desMoment.x
  // You'll need the arm length parameter L, and the drag/thrust ratio kappa

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  float length = L / sqrtf(2.f);

  // F_tot = F0 + F1 + F2 + F3
  // tau_x = (F0 - F1 + F2 - F3) * length     // This is Roll
  // tau_y = (F0 + F1 - F2 - F3) * length     // This is Pitch
  // tau_z = (-F0 + F1 + F2 - F3) * kappa     // This is Yaw

  float A = collThrustCmd;
  float B = momentCmd.x / length;
  float C = momentCmd.y / length;
  float D = momentCmd.z / kappa;

  // I used WolframAlpha to solve the system of equations to get:

  float N0 = (A + B + C - D) / 4.f; // front left
  float N1 = (A - B + C + D) / 4.f; // front right
  float N2 = (A + B - C + D) / 4.f; // rear left
  float N3 = (A - B - C - D) / 4.f; // rear right

  // constrain each motor's thrust bu min and max motor thrust values 

  cmd.desiredThrustsN[0] = CONSTRAIN(N0, minMotorThrust, maxMotorThrust);
  cmd.desiredThrustsN[1] = CONSTRAIN(N1, minMotorThrust, maxMotorThrust);
  cmd.desiredThrustsN[2] = CONSTRAIN(N2, minMotorThrust, maxMotorThrust);
  cmd.desiredThrustsN[3] = CONSTRAIN(N3, minMotorThrust, maxMotorThrust);

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return cmd;
}

V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
  // Calculate a desired 3-axis moment given a desired and current body rate
  // INPUTS: 
  //   pqrCmd: desired body rates [rad/s]
  //   pqr: current or estimated body rates [rad/s]
  // OUTPUT:
  //   return a V3F containing the desired moments for each of the 3 axes

  // HINTS: 
  //  - you can use V3Fs just like scalars: V3F a(1,1,1), b(2,3,4), c; c=a-b;
  //  - you'll need parameters for moments of inertia Ixx, Iyy, Izz
  //  - you'll also need the gain parameter kpPQR (it's a V3F)

  V3F momentCmd;


  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  V3F MOI;
  
  MOI.x = Ixx;
  MOI.y = Iyy;
  MOI.z = Izz;

  momentCmd = MOI * kpPQR * (pqrCmd - pqr);

  // Limit the moments to the Max Torque value
  // if np.linalg.norm(momentCmd) > MAX_TORQUE:
  //	moments_cmd = momentCmd * MAX_TORQUE / np.linalg.norm(moments_cmd)

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return momentCmd;
}

// returns a desired roll and pitch rate 
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
  // Calculate a desired pitch and roll angle rates based on a desired global
  //   lateral acceleration, the current attitude of the quad, and desired
  //   collective thrust command
  // INPUTS: 
  //   accelCmd: desired acceleration in global XY coordinates [m/s2]
  //   attitude: current or estimated attitude of the vehicle
  //   collThrustCmd: desired collective thrust of the quad [N]
  // OUTPUT:
  //   return a V3F containing the desired pitch and roll rates. The Z
  //     element of the V3F should be left at its default value (0)

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the roll/pitch gain kpBank
  //  - collThrustCmd is a force in Newtons! You'll likely want to convert it to acceleration first

  V3F pqrCmd;
  Mat3x3F R = attitude.RotationMatrix_IwrtB();

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  float b_x_dot;
  float b_y_dot;
  float R13_cmd;
  float R23_cmd;
  float R11 = R(0, 0);
  float R12 = R(0, 1);
  float R13 = R(0, 2);
  float R21 = R(1, 0);
  float R22 = R(1, 1);
  float R23 = R(1, 2);
  float R33 = R(2, 2);

  // From lesson 14.16 we know that x_dot_dot = c * R13 and y_dot_dot = c * R23 where c is thrust_cmd / mass
  // R13 is - sin(pitch) and R23 is sin(roll)*cos(pitch)

  float c = -collThrustCmd / mass;

  if (collThrustCmd > 0.0) {
	  R13_cmd = accelCmd.x / c;
	  R23_cmd = accelCmd.y / c;
	  
	  // limit the tilt angles using the max_tilt value
	  R13_cmd = CONSTRAIN(R13_cmd, -maxTiltAngle, maxTiltAngle);
	  R23_cmd = CONSTRAIN(R23_cmd, -maxTiltAngle, maxTiltAngle);

	  b_x_dot = kpBank * (R13_cmd - R13);
	  b_y_dot = kpBank * (R23_cmd - R23);

	  pqrCmd.x = (1 / R33) * (R21 * b_x_dot - R11 * b_y_dot);
	  pqrCmd.y = (1 / R33) * (R22 * b_x_dot - R12 * b_y_dot);
  }
  else {
	  // If thrust is negative or = 0 then set pitch and roll rates to zero
	  pqrCmd.x = 0.0;
	  pqrCmd.y = 0.0;
  }

  pqrCmd.z = 0;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return pqrCmd;
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude, float accelZCmd, float dt)
{
  // Calculate desired quad thrust based on altitude setpoint, actual altitude,
  //   vertical velocity setpoint, actual vertical velocity, and a vertical 
  //   acceleration feed-forward command
  // INPUTS: 
  //   posZCmd, velZCmd: desired vertical position and velocity in NED [m]
  //   posZ, velZ: current vertical position and velocity in NED [m]
  //   accelZCmd: feed-forward vertical acceleration in NED [m/s2]
  //   dt: the time step of the measurements [seconds]
  // OUTPUT:
  //   return a collective thrust command in [N]

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the gain parameters kpPosZ and kpVelZ
  //  - maxAscentRate and maxDescentRate are maximum vertical speeds. Note they're both >=0!
  //  - make sure to return a force, not an acceleration
  //  - remember that for an upright quad in NED, thrust should be HIGHER if the desired Z acceleration is LOWER

  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  float thrust = 0;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  
  float z_err = posZCmd - posZ;

  integratedAltitudeError = integratedAltitudeError + z_err * dt;

  float b_z = R(2,2);   // This is matrix element R33

  //float vel_z_des = kpPosZ * z_err + velZCmd;  // this is the desired velocity

 // vel_z_des = CONSTRAIN(vel_z_des, -maxAscentRate, maxDescentRate);  // need to limit the desired vertical velocity by ascent / decent rates

  float vel_z_des = kpPosZ * z_err + velZCmd;
  vel_z_des = CONSTRAIN(vel_z_des, -maxAscentRate, maxDescentRate);

  float i_term = KiPosZ * integratedAltitudeError;

  float u_1_bar = kpVelZ*(vel_z_des-velZ) + i_term + accelZCmd;   // this is the desired vertical acceleration
  //float u_1_bar = kpVelZ * (kpPosZ * z_err - velZ) + i_term + accelZCmd;   // this is the desired vertical acceleration

  float c = (u_1_bar - CONST_GRAVITY) / b_z;   // this is net acceleration in the drone's z-direction, also accounting for gravity

  thrust = -c * mass;  // need positive thrust to move in negative z direction

  /////////////////////////////// END STUDENT CODE ////////////////////////////
  
  return thrust;
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelFF)
{
  // Calculate a desired horizontal acceleration based on 
  //  desired lateral position/velocity/acceleration and current pose
  // INPUTS: 
  //   posCmd: desired position, in NED [m]
  //   velCmd: desired velocity, in NED [m/s]
  //   pos: current position, NED [m]
  //   vel: current velocity, NED [m/s]
  //   accelCmd: desired acceleration, NED [m/s2]
  // OUTPUT:
  //   return a V3F with desired horizontal accelerations. 
  //     the Z component should be 0
  // HINTS: 
  //  - use CONSTRAIN to cap a value
  //  - use the gain parameters kpPosXY and kpVelXY
  //  - make sure you cap the horizontal velocity and acceleration
  //    to maxSpeedXY and maxAccelXY

  // make sure we don't have any incoming z-component
  accelFF.z = 0;
  velCmd.z = 0;
  posCmd.z = pos.z;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  float speed_cmd = sqrtf(velCmd.x*velCmd.x + velCmd.y*velCmd.y);  // calculate the speed being commanded

  if (speed_cmd > maxSpeedXY) // if the commanded speed is too high then reduce
  {
	  velCmd = velCmd * maxSpeedXY / speed_cmd;
  }

  V3F pos_err = posCmd - pos;
  V3F vel_err = velCmd - vel;

  V3F p_term_xy = kpPosXY * pos_err;
  V3F d_term_xy = kpVelXY * vel_err;

  V3F accelCmd = p_term_xy + d_term_xy + accelFF;
  
  accelCmd.z = 0;

  float ac_cmd = sqrtf(accelCmd.x*accelCmd.x + accelCmd.y*accelCmd.y);  // calculate the magnitude of the acceleration being commanded

  if (ac_cmd > maxAccelXY) // if the commanded acceleration is too high then reduce
  {
	  accelCmd = accelCmd * maxAccelXY / ac_cmd;
  }

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return accelCmd;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{
  // Calculate a desired yaw rate to control yaw to yawCmd
  // INPUTS: 
  //   yawCmd: commanded yaw [rad]
  //   yaw: current yaw [rad]
  // OUTPUT:
  //   return a desired yaw rate [rad/s]
  // HINTS: 
  //  - use fmodf(foo,b) to constrain float foo to range [0,b]
  //  - use the yaw control gain parameter kpYaw

  float yawRateCmd=0;
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  // Since yaw is decoupled from the other directions, we only need a P controller
  
  float Pi = 3.14159f;
  float twoPi = 2.f * Pi;

  yawCmd = fmodf(yawCmd, twoPi);  // constrain yaw to the range(0, 2pi)
  yaw = fmodf(yaw, twoPi);

  float yaw_err = yawCmd - yaw;

  // We have a choice on which way to rotate the drone to get to the desired yaw angle
  // And should pick the direction(CW or CCW) that requires the smaller rotation

  if (yaw_err > Pi)
  {
	  yaw_err = yaw_err - twoPi;
  }
  else if (yaw_err < -Pi)
  {
	  yaw_err = yaw_err + twoPi;
  }
  
  yawRateCmd = kpYaw * yaw_err;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return yawRateCmd;

}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
  curTrajPoint = GetNextTrajectoryPoint(simTime);

  float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);

  // reserve some thrust margin for angle control
  float thrustMargin = .1f*(maxMotorThrust - minMotorThrust);
  collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust+ thrustMargin)*4.f, (maxMotorThrust-thrustMargin)*4.f);
  
  V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);
  
  V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
  desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

  V3F desMoment = BodyRateControl(desOmega, estOmega);

  return GenerateMotorCommands(collThrustCmd, desMoment);
}
