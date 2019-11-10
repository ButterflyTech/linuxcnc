#include "rtapi.h"
#include "motion.h"
#include "rtapi_string.h"

/* Utility routines for kinematics modules
**
**
** map_coordinates_to_jnumbers()
**
** Map a string of coordinate letters to joint numbers sequentially.
** If allow_duplicates==1, a coordinate letter may be specified more
** than once to assign it to multiple joint numbers (the kinematics
** module must support such usage).
**
**   Default mapping if coordinates==NULL is:
**           X:0 Y:1 Z:2 A:3 B:4 C:5 U:6 V:7 W:8
**
**   Example coordinates-to-joints mappings:
**      coordinates=XYZ      X:0   Y:1   Z:2
**      coordinates=ZYX      Z:0   Y:1   X:2
**      coordinates=XYZZZZ   x:0   Y:1   Z:2,3,4,5
**      coordinates=XXYZ     X:0,1 Y:2   Z:3
*/
int map_coordinates_to_jnumbers(char *coordinates,
                                int  max_joints,
                                int  allow_duplicates,
                                int  axis_idx_for_jno[] ) //result
{
    char* errtag="map_coordinates_to_jnumbers: ERROR:\n  ";
    int  jno=0;
    int  found=0;
    int  dups[EMCMOT_MAX_AXIS] = {0};
    char *coords = coordinates;
    char coord_letter[] = {'X','Y','Z','A','B','C','U','V','W'};

    if ( (max_joints <= 0) || (max_joints > EMCMOT_MAX_JOINTS) ) {
        rtapi_print_msg(RTAPI_MSG_ERR,"%s bogus max_joints=%d\n",
          errtag,max_joints);
        return -1;
    }

    // init all axis_idx_for_jno[] to -1 ==> unspecified
    for(jno=0; jno<EMCMOT_MAX_JOINTS; jno++) { axis_idx_for_jno[jno] = -1; }

    if (coords == NULL) { coords = "XYZABCUVW"; }
    jno = 0; // begin: assign joint numbers at 0th coords position
    while (*coords) {
        found = 0;
        switch(*coords) {
          case 'x': case 'X': axis_idx_for_jno[jno]= 0;dups[0]++;found=1;break;
          case 'y': case 'Y': axis_idx_for_jno[jno]= 1;dups[1]++;found=1;break;
          case 'z': case 'Z': axis_idx_for_jno[jno]= 2;dups[2]++;found=1;break;
          case 'a': case 'A': axis_idx_for_jno[jno]= 3;dups[3]++;found=1;break;
          case 'b': case 'B': axis_idx_for_jno[jno]= 4;dups[4]++;found=1;break;
          case 'c': case 'C': axis_idx_for_jno[jno]= 5;dups[5]++;found=1;break;
          case 'u': case 'U': axis_idx_for_jno[jno]= 6;dups[6]++;found=1;break;
          case 'v': case 'V': axis_idx_for_jno[jno]= 7;dups[7]++;found=1;break;
          case 'w': case 'W': axis_idx_for_jno[jno]= 8;dups[8]++;found=1;break;
          case ' ': case '\t': coords++;continue; //whitespace
        }
        if (found) {
            coords++; // next coordinates letter
            jno++;    // next joint number
        } else {
            rtapi_print_msg(RTAPI_MSG_ERR,
              "%s Invalid character '%c' in coordinates '%s'\n",
              errtag,*coords,coordinates);
            return -1;
        }
        if (jno > max_joints) {
            rtapi_print_msg(RTAPI_MSG_ERR,
              "%s too many coordinates <%s> for max_joints=%d\n",
              errtag,coordinates,max_joints);
            return -1;
        }
    } // while

    if (!found) {
        rtapi_print_msg(RTAPI_MSG_ERR,"%s missing coordinates '%s'\n",
          errtag,coordinates);
        return -1;
    }
    if (!allow_duplicates) {
        int ano;
        for(ano=0; ano<EMCMOT_MAX_AXIS; ano++) {
            if (dups[ano] > 1) {
                rtapi_print_msg(RTAPI_MSG_ERR,
                "%s duplicates not allowed in coordinates=%s, letter=%c\n",
                errtag,coordinates,coord_letter[ano]);
                return -1;
            }
        }
    }
    return 0;
} //map_coordinates_to_jnumbers()

// IDENTITY kinematics implementation (local)
// joint number assignments when switched to identity kinematics
// are set by module coordinates= parameter,default ordering is:
#define DEFAULT_LETTER_TO_JOINT_MAP "XYZABC"
static int _JX = -1;
static int _JY = -1;
static int _JZ = -1;
static int _JA = -1;
static int _JB = -1;
static int _JC = -1;
static int _JU = -1;
static int _JV = -1;
static int _JW = -1;

static int identity_kinematics_inited = 0;

int identityKinematicsSetup(char *coordinates,
                      int  max_joints,
                      int  allow_duplicates
                     )
{
    static int axis_idx_for_jno[EMCMOT_MAX_JOINTS];
    int jno;
    int islathe;
    int show=0;

    if (map_coordinates_to_jnumbers(coordinates,
                                    max_joints,
                                    allow_duplicates,
                                    axis_idx_for_jno)) {
       return -1; //mapping failed
    }

    for (jno=0; jno<EMCMOT_MAX_JOINTS; jno++) {
      if (axis_idx_for_jno[jno] == 0) {_JX = jno;}
      if (axis_idx_for_jno[jno] == 1) {_JY = jno;}
      if (axis_idx_for_jno[jno] == 2) {_JZ = jno;}
      if (axis_idx_for_jno[jno] == 3) {_JA = jno;}
      if (axis_idx_for_jno[jno] == 4) {_JB = jno;}
      if (axis_idx_for_jno[jno] == 5) {_JC = jno;}
      if (axis_idx_for_jno[jno] == 6) {_JU = jno;}
      if (axis_idx_for_jno[jno] == 7) {_JV = jno;}
      if (axis_idx_for_jno[jno] == 8) {_JW = jno;}
    }

#if 0
    rtapi_print("\nidentity_kin_init coordinates=%s identity assignments:\n",
               coordinates);
    for (jno=0; jno<EMCMOT_MAX_JOINTS; jno++) {
        if (axis_idx_for_jno[jno] == -1) break; //fini
        rtapi_print("   Joint %d ==> Axis %c\n",
                   jno,*("XYZABCUVW"+axis_idx_for_jno[jno]));
    }
#endif

    /* print message for unconventional ordering;
    **   a) duplicate coordinate letters
    **   b) letters not ordered by "XYZABCUVW" sequence
    **      (use kinstype=both works best for these)
    */
    for (jno=0; jno<EMCMOT_MAX_JOINTS; jno++) {
        if (axis_idx_for_jno[jno] == -1) break; //fini
        if (axis_idx_for_jno[jno] != jno) { show++; } //not default order
    }
    islathe = !strcasecmp(coordinates,"xz"); // no show if simple lathe
    if (show && !islathe) {
        rtapi_print("\nidentityKinematicsSetup: coordinates:%s\n", coordinates);
        char *p="XYZABCUVW";
        for (jno=0; jno<EMCMOT_MAX_JOINTS; jno++) {
            if (axis_idx_for_jno[jno] == -1) break; //fini
            rtapi_print("   Joint %d ==> Axis %c\n",
                       jno,*(p+axis_idx_for_jno[jno]));
        }
        if (kinematicsType() != KINEMATICS_BOTH) {
            rtapi_print("identityKinematicsSetup: Recommend: kinstype=both\n");
        }
        rtapi_print("\n");
    }

    identity_kinematics_inited = 1;
    return 0;
} // identityKinematicsSetup()

int identityKinematicsForward(const double *joints,
                              EmcPose * pos,
                              const KINEMATICS_FORWARD_FLAGS * fflags,
                              KINEMATICS_INVERSE_FLAGS * iflags)
{
    if (!identity_kinematics_inited) {
        rtapi_print_msg(RTAPI_MSG_ERR,
            "identityKinematicsForward: not initialized\n");
        return -1;
    }

    if (_JX >= 0) pos->tran.x = joints[_JX];
    if (_JY >= 0) pos->tran.y = joints[_JY];
    if (_JZ >= 0) pos->tran.z = joints[_JZ];
    if (_JA >= 0) pos->a      = joints[_JA];
    if (_JB >= 0) pos->b      = joints[_JB];
    if (_JC >= 0) pos->c      = joints[_JC];
    if (_JU >= 0) pos->u      = joints[_JU];
    if (_JV >= 0) pos->v      = joints[_JV];
    if (_JW >= 0) pos->w      = joints[_JW];
    return 0;
} // identityKinematicsForward()

int identityKinematicsInverse(const EmcPose * pos,
                              double *joints,
                              const KINEMATICS_INVERSE_FLAGS * iflags,
                              KINEMATICS_FORWARD_FLAGS * fflags)
{
    if (!identity_kinematics_inited) {
        rtapi_print_msg(RTAPI_MSG_ERR,
            "identityKinematicsInverse: not initialized\n");
        return -1;
    }
    if (_JX >= 0) joints[_JX] = pos->tran.x;
    if (_JY >= 0) joints[_JY] = pos->tran.y;
    if (_JZ >= 0) joints[_JZ] = pos->tran.z;
    if (_JA >= 0) joints[_JA] = pos->a;
    if (_JB >= 0) joints[_JB] = pos->b;
    if (_JC >= 0) joints[_JC] = pos->c;
    if (_JU >= 0) joints[_JU] = pos->u;
    if (_JV >= 0) joints[_JV] = pos->v;
    if (_JW >= 0) joints[_JW] = pos->w;
    return 0;
} // identityKinematicsInverse()
