/********************************************************************
TEST: switchable kinematics identity or genserkins

1) genser kinematics provided by genserfuncs.c (shared with genserkins.c0
2) uses same pin names as genserkins
3) tesing mm configs: increase GO_REAL_EPSILON from 1e-7 to 1e-6

NOTE:
a) requires *exactly* 6 joints
b) identity assigments can use any of xyzabcuvw
   but should agree with [TRAJ]COORDINATES
   and may be confusing

www refs:

frame-larger-than:
https://www.mail-archive.com/emc-developers@lists.sourceforge.net/msg03790.html

angles:
https://www.mail-archive.com/emc-developers@lists.sourceforge.net/msg15285.html
*/

//----------------------------------------------------------------------
// genserKinematicsInverse() is 5104 with buster amd64 gcc 8.3.0-6
//#pragma GCC diagnostic error   "-Wframe-larger-than=6000"
  #pragma GCC diagnostic warning "-Wframe-larger-than=6000"

#include "rtapi.h"
#include "rtapi_app.h"
#include <rtapi_string.h>
#include "genserkins.h"
#include "motion.h"     // EMCMOT_MAX_JOINTS

//-7 is system defined -3 ok, -4 ok, -5 ok,-6 ok (mm system)
#undef  GO_REAL_EPSILON
#define GO_REAL_EPSILON (1e-6)

//*********************************************************************
struct switchdata {
    hal_bit_t *kinstype_is_0;
    hal_bit_t *kinstype_is_1;
    hal_bit_t *kinstype_is_2;
} *switchdata = 0;

static hal_u32_t  switchkins_type;

int kinematicsSwitchable() {return 1;}

int kinematicsSwitch(int new_kins_type)
{
    switchkins_type = new_kins_type;
    switch (switchkins_type) {
        case 0: rtapi_print_msg(RTAPI_MSG_INFO,
                "kinematicsSwitch:genserkins\n");
                *switchdata->kinstype_is_0 = 1;
                *switchdata->kinstype_is_1 = 0;
                *switchdata->kinstype_is_2 = 0;
                break;
        case 1: rtapi_print_msg(RTAPI_MSG_INFO,
                "kinematicsSwitch:Identity\n");
                *switchdata->kinstype_is_0 = 0;
                *switchdata->kinstype_is_1 = 1;
                *switchdata->kinstype_is_2 = 0;
                break;
        case 2: rtapi_print_msg(RTAPI_MSG_INFO,
                "kinematicsSwitch:switchtype2\n");
                *switchdata->kinstype_is_0 = 0;
                *switchdata->kinstype_is_1 = 0;
                *switchdata->kinstype_is_2 = 1;
                break;
       default: rtapi_print_msg(RTAPI_MSG_ERR,
                "kinematicsSwitch:BAD VALUE <%d>\n",
                switchkins_type);
                *switchdata->kinstype_is_1 = 0;
                *switchdata->kinstype_is_0 = 0;
                *switchdata->kinstype_is_2 = 0;
                return -1; // FAIL
    }

    return 0; // 0==> no error
} // kinematicsSwitch()

/* main function called by LinuxCNC for forward Kins */
int kinematicsForward(const double *joints,
                      EmcPose * pos,
                      const KINEMATICS_FORWARD_FLAGS * fflags,
                      KINEMATICS_INVERSE_FLAGS * iflags)
{
    switch (switchkins_type) {
       case 0: return      genserKinematicsForward(joints, pos, fflags, iflags);break;
       case 2: return switchtype2KinematicsForward(joints, pos, fflags, iflags);break;
      default: return    identityKinematicsForward(joints, pos, fflags, iflags);
    }

    return 0;
} // kinematicsForward()


int kinematicsInverse(const EmcPose * pos,
                      double *joints,
                      const KINEMATICS_INVERSE_FLAGS * iflags,
                      KINEMATICS_FORWARD_FLAGS * fflags)
{
    switch (switchkins_type) {
       case 0: return      genserKinematicsInverse(pos, joints, iflags, fflags);break;
       case 2: return switchtype2KinematicsInverse(pos, joints, iflags, fflags);break;
      default: return    identityKinematicsInverse(pos, joints, iflags, fflags);
    }

    return 0;
} // kinematicsInverse()

int kinematicsHome(EmcPose * world,
    double *joint,
    KINEMATICS_FORWARD_FLAGS * fflags, KINEMATICS_INVERSE_FLAGS * iflags)
{
    /* use joints, set world */
    return kinematicsForward(joint, world, fflags, iflags);
}

KINEMATICS_TYPE kinematicsType()
{
    return KINEMATICS_BOTH;
}

// support 6 joints (GENSER_MAX_JOINTS)
static char *coordinates = "XYZABC";
RTAPI_MP_STRING(coordinates, "Axes-to-joints-identity-ordering");

EXPORT_SYMBOL(kinematicsSwitchable);
EXPORT_SYMBOL(kinematicsSwitch);
EXPORT_SYMBOL(kinematicsType);
EXPORT_SYMBOL(kinematicsForward);
EXPORT_SYMBOL(kinematicsInverse);
MODULE_LICENSE("GPL");

static int comp_id;

int rtapi_app_main(void)
{
#define DISALLOW_DUPLICATES 0
    int res = 0;
    int njoints = (int)strlen(coordinates);

    if (njoints != GENSER_MAX_JOINTS) {
        rtapi_print_msg(RTAPI_MSG_ERR,
             "genser-switchkins requires exactly %d letters in coordinates=%s\n",
             GENSER_MAX_JOINTS,coordinates);
        return -1;
    }

    if (identityKinematicsSetup(coordinates,
                          GENSER_MAX_JOINTS,
                          DISALLOW_DUPLICATES)) goto error;

    comp_id = hal_init("genser-switchkins");
    if (comp_id < 0) return comp_id;

    if (genser_hal_setup(comp_id)) goto error;
    switchdata = hal_malloc(sizeof(struct switchdata));

    if (!switchdata) goto error;
    if((res=hal_pin_bit_new("kinstype.is-0", HAL_OUT,
           &(switchdata->kinstype_is_0), comp_id)) < 0)
        goto error;
    if((res=hal_pin_bit_new("kinstype.is-1", HAL_OUT,
            &(switchdata->kinstype_is_1), comp_id)) < 0)
        goto error;
    if((res=hal_pin_bit_new("kinstype.is-2", HAL_OUT,
            &(switchdata->kinstype_is_2), comp_id)) < 0)
        goto error;


    switchkins_type   = 0;            // startup with default (0) type
    kinematicsSwitch(switchkins_type);

    rtapi_print("genser-switchkins GO_REAL_EPSILON=%g\n",GO_REAL_EPSILON);
    hal_ready(comp_id);
    return 0;

  error:
    hal_exit(comp_id);
    return res;
} // rtapi_app_main()

void rtapi_app_exit(void)
{
    hal_exit(comp_id);
}
