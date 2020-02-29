/********************************************************************
* Description: trivkins.c
*   general trivkins for 3 axis Cartesian machine
*
*   Derived from a work by Fred Proctor & Will Shackleford
*
* License: GPL Version 2
*    
* Copyright (c) 2009 All rights reserved.
*
********************************************************************/

#include "motion.h"
#include "hal.h"
#include "rtapi.h"
#include "rtapi.h"      /* RTAPI realtime OS API */
#include "rtapi_app.h"  /* RTAPI realtime module decls */
#include "rtapi_math.h"
#include "rtapi_string.h"
#include "kinematics.h"


#define SET(f) pos->f = joints[i]

int kinematicsForward(const double *joints,
                      EmcPose * pos,
                      const KINEMATICS_FORWARD_FLAGS * fflags,
                      KINEMATICS_INVERSE_FLAGS * iflags)
{
    return identityKinematicsForward(joints, pos, fflags, iflags);
}

int kinematicsInverse(const EmcPose * pos,
                      double *joints,
                      const KINEMATICS_INVERSE_FLAGS * iflags,
                      KINEMATICS_FORWARD_FLAGS * fflags)
{
    return identityKinematicsInverse(pos, joints, iflags, fflags);
}

/* implemented for these kinematics as giving joints preference */
int kinematicsHome(EmcPose * world,
                   double *joint,
                   KINEMATICS_FORWARD_FLAGS * fflags,
                   KINEMATICS_INVERSE_FLAGS * iflags)
{
    *fflags = 0;
    *iflags = 0;

    return kinematicsForward(joint, world, fflags, iflags);
}

static KINEMATICS_TYPE ktype = -1;

KINEMATICS_TYPE kinematicsType()
{
    return ktype;
}

#define TRIVKINS_DEFAULT_COORDINATES "XYZABCUVW"
static char *coordinates = TRIVKINS_DEFAULT_COORDINATES;
RTAPI_MP_STRING(coordinates, "Existing Axes");

static char *kinstype = "1"; // use KINEMATICS_IDENTITY
RTAPI_MP_STRING(kinstype, "Kinematics Type (Identity,Both)");

KINS_NOT_SWITCHABLE
EXPORT_SYMBOL(kinematicsType);
EXPORT_SYMBOL(kinematicsForward);
EXPORT_SYMBOL(kinematicsInverse);
MODULE_LICENSE("GPL");

static int comp_id;

int rtapi_app_main(void) {
    #define ALLOW_DUPLICATES 1
#if 0
    static int axis_idx_for_jno[EMCMOT_MAX_JOINTS];
    if (map_coordinates_to_jnumbers(coordinates,
                                    EMCMOT_MAX_JOINTS,
                                    ALLOW_DUPLICATES,
                                    axis_idx_for_jno)) {
       return -1; //mapping failed
    }
    // identityKinematics uses same mapping:
#endif
    switch (*kinstype) {
      case 'b': case 'B': ktype = KINEMATICS_BOTH;         break;
      case 'f': case 'F': ktype = KINEMATICS_FORWARD_ONLY; break;
      case 'i': case 'I': ktype = KINEMATICS_INVERSE_ONLY; break;
      case '1': default:  ktype = KINEMATICS_IDENTITY;
    }

    if (identityKinematicsSetup(coordinates,
                                EMCMOT_MAX_JOINTS,
                                ALLOW_DUPLICATES)) {
       return -1; //setup failed
    }

    comp_id = hal_init("trivkins");
    if(comp_id < 0) return comp_id;

#if 0
    /* print message for unconventional ordering;
    **   a) duplicate coordinate letters
    **   b) letters not ordered by "XYZABCUVW" sequence
    **      (use kinstype=both works best for these)
    */
    {
        int jno,islathe,show=0;
        for (jno=0; jno<EMCMOT_MAX_JOINTS; jno++) {
            if (axis_idx_for_jno[jno] == -1) break; //fini
            if (axis_idx_for_jno[jno] != jno) { show++; } //not default order
        }
        islathe = !strcasecmp(coordinates,"xz"); // no show if simple lathe
        if (show && !islathe) {
            rtapi_print("\ntrivkins: coordinates:%s\n", coordinates);
            char *p="XYZABCUVW";
            for (jno=0; jno<EMCMOT_MAX_JOINTS; jno++) {
                if (axis_idx_for_jno[jno] == -1) break; //fini
                rtapi_print("   Joint %d ==> Axis %c\n",
                           jno,*(p+axis_idx_for_jno[jno]));
            }
            if (ktype != KINEMATICS_BOTH) {
                rtapi_print("trivkins: Recommend: kinstype=both\n");
            }
            rtapi_print("\n");
        }
    }
//duh
#endif

    hal_ready(comp_id);
    return 0;
}

void rtapi_app_exit(void) { hal_exit(comp_id); }
