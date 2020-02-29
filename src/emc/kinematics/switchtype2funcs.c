// this is a template file for switchkins_type==2 kinematics

// typical includes:
#include "kinematics.h"
#include "rtapi_math.h"

// Add for kins based on genserkins:
// #include "genserkins.h" //includes gomath,hal

//**********************************************************************
// static local variables and functions go here
//**********************************************************************

int switchtype2KinematicsForward(const double *joint,
                                 struct EmcPose * world,
                                 const KINEMATICS_FORWARD_FLAGS * fflags,
                                 KINEMATICS_INVERSE_FLAGS * iflags)
{
    // test linker inclusion of standard math functions:
    static volatile double tst=0;tst=sqrt(tst); // ensure -lm used

    return identityKinematicsForward(joint,world,fflags,iflags);
}

int switchtype2KinematicsInverse(const EmcPose * pos,
                                 double *joint,
                                 const KINEMATICS_INVERSE_FLAGS * iflags,
                                 KINEMATICS_FORWARD_FLAGS * fflags)
{
    return identityKinematicsForward(pos,joint,iflags,fflags);
}
