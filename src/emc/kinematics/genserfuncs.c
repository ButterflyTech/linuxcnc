/********************************************************************
* Description: genserfuncs.c
*   Kinematics for a generalised serial kinematics machine
*
*   Derived from a work by Fred Proctor,
*   changed to work with emc2 and HAL
*
* Adapting Author: Alex Joni
* License: GPL Version 2
* System: Linux
*    
*******************************************************************

  These are the forward and inverse kinematic functions for a general
  serial-link manipulator. Thanks to Herman Bruyninckx and John
  Hallam at http://www.roble.info/ for this.

  The functions are general enough to be configured for any serial
  configuration.  
  The kinematics use Denavit-Hartenberg definition for the joint and
  links. The DH definitions are the ones used by John J Craig in
  "Introduction to Robotics: Mechanics and Control"
  The parameters for the manipulator are defined by hal pins.
  Currently the type of the joints is hardcoded to ANGULAR, although 
  the kins support both ANGULAR and LINEAR axes.
  
  TODO:
    * make number of joints a loadtime parameter
    * add HAL pins for all settable parameters, including joint type: ANGULAR / LINEAR
*/

#include "rtapi_math.h"
#include "gotypes.h"    /* go_result, go_integer */
#include "gomath.h"     /* go_pose */
#include "genserkins.h" /* these decls */
#include "kinematics.h"
#include "hal.h"

#ifdef RTAPI
#include "rtapi.h"
//#include "rtapi_app.h"
#endif

struct haldata {
    hal_u32_t     *max_iterations;
    hal_u32_t     *last_iterations;
    hal_float_t   *a[GENSER_MAX_JOINTS];
    hal_float_t   *alpha[GENSER_MAX_JOINTS];
    hal_float_t   *d[GENSER_MAX_JOINTS];
    hal_s32_t     *unrotate[GENSER_MAX_JOINTS];
    genser_struct *kins;
    go_pose *pos; // used in various functions, we malloc it
                  // only once in rtapi_app_main
} *haldata = 0;

double j[GENSER_MAX_JOINTS];

#define A(i) (*(haldata->a[i]))
#define ALPHA(i) (*(haldata->alpha[i]))
#define D(i) (*(haldata->d[i]))

#define KINS_PTR (haldata->kins)

#if GENSER_MAX_JOINTS < 6
#error GENSER_MAX_JOINTS must be at least 6; fix genserkins.h
#endif
//---------------------------------------------------------------------------------------------------------------------------
double Ao = 0;
double Bo = 0;
double Co = 0;

double aw = 0;
double bw = 0;
double cw = 0;
double To[4][4];
double Tt[4][4];
double Tw[4][4];

int elbup = 1;
int wrflp = 0;
int wsing = 0;

double eps = 0.0001;

//=======================================================================
int fwABC2(double ux, double uy, double uz, double vx, double vy, 
           double vz, double wz)
{
// to Craig: p40 RPY(c,b,a) => (Rz,Ry,Rz)

    bw = atan2(-uz, sqrt(ux*ux + uy*uy) );
    if (fabs(fabs(bw) - M_PI/2.0) > eps) {
       aw = atan2(vz, wz);
       cw = atan2(uy, ux);
    } else if (bw > 0.0) {
       aw = atan2(vx, vy);
       bw = M_PI/2.0;
       cw = 0;
    } else {
       aw = -atan2(vx, vy);
       bw = -M_PI/2.0;
       cw = 0;
    }
    if (aw < 0) {aw = aw + 2*M_PI;}

    return 0;
}
//=======================================================================   
static void MatMult(double A[][4], double B[][4], double C[][4])
{
  int i, j, k;
  
  for (i=0; i<=3; ++i){
    for (j=0; j<=3; ++j){
      C[i][j] = 0;
      for (k=0; k<=3; ++k){
        C[i][j] = C[i][j] + A[i][k]*B[k][j];
      }
    }
  }
}
//=====================================================================
static void PoseMat(double X, double Y, double Z,
                    double A, double B, double C, double T[][4])
{
double cr, sr, cp, sp, cy, sy;
  
       sr = sin(A); cr = cos(A);
       sp = sin(B); cp = cos(B);
       sy = sin(C); cy = cos(C);
    
       T[0][0] = cp*cy;
       T[1][0] = cp*sy;
       T[2][0] = -sp;
       T[3][0] = 0;

       T[0][1] = sr*sp*cy - cr*sy;
       T[1][1] = sr*sp*sy + cr*cy;
       T[2][1] = sr*cp;
       T[3][1] = 0;

       T[0][2] = cr*sp*cy + sr*sy;
       T[1][2] = cr*sp*sy - sr*cy;
       T[2][2] = cr*cp;
       T[3][2] = 0;

       T[0][3] = X;
       T[1][3] = Y;
       T[2][3] = Z;
       T[3][3] = 1; 
}
//----------------------------------------------------------------------------------------------------------------------------------



int genser_kin_init(void) {
    genser_struct *genser = KINS_PTR;
    int t;

    /* init them all and make them revolute joints */
    /* FIXME: should allow LINEAR joints based on HAL param too */
    for (t = 0; t < GENSER_MAX_JOINTS; t++) {
        genser->links[t].u.dh.a = A(t);
        genser->links[t].u.dh.alpha = ALPHA(t);
        genser->links[t].u.dh.d = D(t);
        genser->links[t].u.dh.theta = 0;
        genser->links[t].type = GO_LINK_DH;
        genser->links[t].quantity = GO_QUANTITY_ANGLE;
    }

    /* set a select few to make it PUMA-like */
    // FIXME-AJ: make a hal pin, also set number of joints based on it
    genser->link_num = 6;

    return GO_RESULT_OK;
}

/* compute the forward jacobian function: 
   the jacobian is a linear aproximation of the kinematics function.
   It is calculated using derivation of the position transformation matrix, 
   and usually used for feeding velocities through it.
   It is analytically possible to calculate the inverse of the jacobian 
   (sometimes only the pseudoinverse) and to use that for the inverse kinematics.
*/
int compute_jfwd(go_link * link_params, 
                 int link_number, 
                 go_matrix * Jfwd, 
                 go_pose * T_L_0) 
{
    GO_MATRIX_DECLARE(Jv, Jvstg, 3, GENSER_MAX_JOINTS);
    GO_MATRIX_DECLARE(Jw, Jwstg, 3, GENSER_MAX_JOINTS);
    GO_MATRIX_DECLARE(R_i_ip1, R_i_ip1stg, 3, 3);
    GO_MATRIX_DECLARE(scratch, scratchstg, 3, GENSER_MAX_JOINTS);
    GO_MATRIX_DECLARE(R_inv, R_invstg, 3, 3);
    go_pose pose;
    go_quat quat;
    go_vector P_ip1_i[3];
    int row, col;

    /* init matrices to possibly smaller size */
    go_matrix_init(Jv, Jvstg, 3, link_number);
    go_matrix_init(Jw, Jwstg, 3, link_number);
    go_matrix_init(R_i_ip1, R_i_ip1stg, 3, 3);
    go_matrix_init(scratch, scratchstg, 3, link_number);
    go_matrix_init(R_inv, R_invstg, 3, 3);

    Jv.el[0][0] = 0, Jv.el[1][0] = 0, Jv.el[2][0] = (GO_QUANTITY_LENGTH == link_params[0].quantity ? 1 : 0);
    Jw.el[0][0] = 0, Jw.el[1][0] = 0, Jw.el[2][0] = (GO_QUANTITY_ANGLE == link_params[0].quantity ? 1 : 0);

    /* initialize inverse rotational transform */
    if (GO_LINK_DH == link_params[0].type) {
        go_dh_pose_convert(&link_params[0].u.dh, &pose);
    } else if (GO_LINK_PP == link_params[0].type) {
        pose = link_params[0].u.pp.pose;
    } else {
        return GO_RESULT_IMPL_ERROR;
    }

    *T_L_0 = pose;

    for (col = 1; col < link_number; col++) {
        /* T_ip1_i */
        if (GO_LINK_DH == link_params[col].type) {
            go_dh_pose_convert(&link_params[col].u.dh, &pose);
        } else if (GO_LINK_PP == link_params[col].type) {
            pose = link_params[col].u.pp.pose;
        } else {
            return GO_RESULT_IMPL_ERROR;
    }

        go_cart_vector_convert(&pose.tran, P_ip1_i);
        go_quat_inv(&pose.rot, &quat);
        go_quat_matrix_convert(&quat, &R_i_ip1);

        /* Jv */
        go_matrix_vector_cross(&Jw, P_ip1_i, &scratch);
        go_matrix_matrix_add(&Jv, &scratch, &scratch);
        go_matrix_matrix_mult(&R_i_ip1, &scratch, &Jv);
        Jv.el[0][col] = 0, Jv.el[1][col] = 0, Jv.el[2][col] = (GO_QUANTITY_LENGTH == link_params[col].quantity ? 1 : 0);
        /* Jw */
        go_matrix_matrix_mult(&R_i_ip1, &Jw, &Jw);
        Jw.el[0][col] = 0, Jw.el[1][col] = 0, Jw.el[2][col] = (GO_QUANTITY_ANGLE == link_params[col].quantity ? 1 : 0);
        if (GO_LINK_DH == link_params[col].type) {
            go_dh_pose_convert(&link_params[col].u.dh, &pose);
        } else if (GO_LINK_PP == link_params[col].type) {
            pose = link_params[col].u.pp.pose;
        } else {
            return GO_RESULT_IMPL_ERROR;
        }
        go_pose_pose_mult(T_L_0, &pose, T_L_0);
    }

    /* rotate back into {0} frame */
    go_quat_matrix_convert(&T_L_0->rot, &R_inv);
    go_matrix_matrix_mult(&R_inv, &Jv, &Jv);
    go_matrix_matrix_mult(&R_inv, &Jw, &Jw);

    /* put Jv atop Jw in J */
    for (row = 0; row < 6; row++) {
        for (col = 0; col < link_number; col++) {
            if (row < 3) {
                Jfwd->el[row][col] = Jv.el[row][col];
            } else {
                Jfwd->el[row][col] = Jw.el[row - 3][col];
            }
        }
    }

    return GO_RESULT_OK;
}

/* compute the inverse of the jacobian matrix */
int compute_jinv(go_matrix * Jfwd, go_matrix * Jinv)
{
    int retval;
    GO_MATRIX_DECLARE(JT, JTstg, GENSER_MAX_JOINTS, 6);

    /* compute inverse, or pseudo-inverse */
    if (Jfwd->rows == Jfwd->cols) {
        retval = go_matrix_inv(Jfwd, Jinv);
        if (GO_RESULT_OK != retval)
            return retval;
    } else if (Jfwd->rows < Jfwd->cols) {
        /* underdetermined, optimize on smallest sum of square of speeds */
        /* JT(JJT)inv */
        GO_MATRIX_DECLARE(JJT, JJTstg, 6, 6);

        go_matrix_init(JT, JTstg, Jfwd->cols, Jfwd->rows);
        go_matrix_init(JJT, JJTstg, Jfwd->rows, Jfwd->rows);
        go_matrix_transpose(Jfwd, &JT);
        go_matrix_matrix_mult(Jfwd, &JT, &JJT);
        retval = go_matrix_inv(&JJT, &JJT);
        if (GO_RESULT_OK != retval)
            return retval;
        go_matrix_matrix_mult(&JT, &JJT, Jinv);
    } else {
        /* overdetermined, do least-squares best fit */
        /* (JTJ)invJT */
        GO_MATRIX_DECLARE(JTJ, JTJstg, GENSER_MAX_JOINTS, GENSER_MAX_JOINTS);

        go_matrix_init(JT, JTstg, Jfwd->cols, Jfwd->rows);
        go_matrix_init(JTJ, JTJstg, Jfwd->cols, Jfwd->cols);
        go_matrix_transpose(Jfwd, &JT);
        go_matrix_matrix_mult(&JT, Jfwd, &JTJ);
        retval = go_matrix_inv(&JTJ, &JTJ);
        if (GO_RESULT_OK != retval)
            return retval;
        go_matrix_matrix_mult(&JTJ, &JT, Jinv);
    }

    return GO_RESULT_OK;
}

int genser_kin_jac_inv(void *kins,
    const go_pose * pos,
    const go_screw * vel, const go_real * joints, go_real * jointvels)
{
    genser_struct *genser = (genser_struct *) kins;
    GO_MATRIX_DECLARE(Jfwd, Jfwd_stg, 6, GENSER_MAX_JOINTS);
    GO_MATRIX_DECLARE(Jinv, Jinv_stg, GENSER_MAX_JOINTS, 6);
    go_pose T_L_0;
    go_link linkout[GENSER_MAX_JOINTS];
    go_real vw[6];
    int link;
    int retval;

    go_matrix_init(Jfwd, Jfwd_stg, 6, genser->link_num);
    go_matrix_init(Jinv, Jinv_stg, GENSER_MAX_JOINTS, 6);

    for (link = 0; link < genser->link_num; link++) {
        retval =
            go_link_joint_set(&genser->links[link], joints[link],
            &linkout[link]);
        if (GO_RESULT_OK != retval)
            return retval;
    }
    retval = compute_jfwd(linkout, genser->link_num, &Jfwd, &T_L_0);
    if (GO_RESULT_OK != retval)
        return retval;
    retval = compute_jinv(&Jfwd, &Jinv);
    if (GO_RESULT_OK != retval)
        return retval;

    vw[0] = vel->v.x;
    vw[1] = vel->v.y;
    vw[2] = vel->v.z;
    vw[3] = vel->w.x;
    vw[4] = vel->w.y;
    vw[5] = vel->w.z;

    return go_matrix_vector_mult(&Jinv, vw, jointvels);
}

int genser_kin_jac_fwd(void *kins,
    const go_real * joints,
    const go_real * jointvels, const go_pose * pos, go_screw * vel)
{
    genser_struct *genser = (genser_struct *) kins;
    GO_MATRIX_DECLARE(Jfwd, Jfwd_stg, 6, GENSER_MAX_JOINTS);
    go_pose T_L_0;
    go_link linkout[GENSER_MAX_JOINTS];
    go_real vw[6];
    int link;
    int retval;

    go_matrix_init(Jfwd, Jfwd_stg, 6, genser->link_num);

    for (link = 0; link < genser->link_num; link++) {
        retval =
            go_link_joint_set(&genser->links[link], joints[link],
            &linkout[link]);
        if (GO_RESULT_OK != retval)
            return retval;
    }

    retval = compute_jfwd(linkout, genser->link_num, &Jfwd, &T_L_0);
    if (GO_RESULT_OK != retval)
        return retval;

    go_matrix_vector_mult(&Jfwd, jointvels, vw);
    vel->v.x = vw[0];
    vel->v.y = vw[1];
    vel->v.z = vw[2];
    vel->w.x = vw[3];
    vel->w.y = vw[4];
    vel->w.z = vw[5];

    return GO_RESULT_OK;
}

/* main function called by emc2 for forward Kins */
int genserKinematicsForward(const double *joint, 
                            EmcPose * world, 
                            const KINEMATICS_FORWARD_FLAGS * fflags, 
                            KINEMATICS_INVERSE_FLAGS * iflags) {

    go_pose *pos;
    go_rpy rpy;
    go_real jcopy[GENSER_MAX_JOINTS]; // will hold the radian conversion of joints
    int ret = 0;
    int i, changed=0;
    
    for (i=0; i< 6; i++)  {
        // FIXME - debug hack
        if (!GO_ROT_CLOSE(j[i],joint[i])) changed = 1;
        // convert to radians to pass to genser_kin_fwd
        jcopy[i] = joint[i] * PM_PI / 180;
        if ((i) && *(haldata->unrotate[i]))
            jcopy[i] -= *(haldata->unrotate[i])*jcopy[i-1];
    }
    
    if (changed) {
        for (i=0; i< 6; i++)
            j[i] = joint[i];
            // rtapi_print("genserKinematicsForward(joints: %f %f %f %f %f %f)\n",
            //joint[0],joint[1],joint[2],joint[3],joint[4],joint[5]);
    }
    // AJ: convert from emc2 coords (XYZABC - which are actually rpy euler
    // angles)
    // to go angles (quaternions)
    pos = haldata->pos;
    rpy.y = world->c * PM_PI / 180;
    rpy.p = world->b * PM_PI / 180;
    rpy.r = world->a * PM_PI / 180;

    go_rpy_quat_convert(&rpy, &pos->rot);
    pos->tran.x = world->tran.x;
    pos->tran.y = world->tran.y;
    pos->tran.z = world->tran.z;

    // pos will be the world location
    // jcopy: joitn position in radians
    ret = genser_kin_fwd(KINS_PTR, jcopy, pos);
    if (ret < 0)
        return ret;

    // AJ: convert back to emc2 coords
    ret = go_quat_rpy_convert(&pos->rot, &rpy);
    if (ret < 0)
        return ret;
    world->tran.x = pos->tran.x;
    world->tran.y = pos->tran.y;
    world->tran.z = pos->tran.z;
    world->a = rpy.r * 180 / PM_PI;
    world->b = rpy.p * 180 / PM_PI;
    world->c = rpy.y * 180 / PM_PI;

    if (changed) {
// rtapi_print("genserKinematicsForward(world: %f %f %f %f %f %f)\n", world->tran.x, world->tran.y, world->tran.z, world->a, world->b, world->c);
    }
    return 0;
}

int genser_kin_fwd(void *kins, const go_real * joints, go_pose * pos)
{
    genser_struct *genser = kins;
    go_link linkout[GENSER_MAX_JOINTS];

    int link;
    int retval;

    genser_kin_init();

    for (link = 0; link < genser->link_num; link++) {
        retval = go_link_joint_set(&genser->links[link], joints[link], &linkout[link]);
        if (GO_RESULT_OK != retval)
            return retval;
    }

    retval = go_link_pose_build(linkout, genser->link_num, pos);
    if (GO_RESULT_OK != retval)
        return retval;

    return GO_RESULT_OK;
}

int genserKinematicsInverse(const EmcPose * world,
                            double *joints,
                            const KINEMATICS_INVERSE_FLAGS * iflags,
                            KINEMATICS_FORWARD_FLAGS * fflags)
{

    genser_struct *genser = KINS_PTR;
    GO_MATRIX_DECLARE(Jfwd, Jfwd_stg, 6, GENSER_MAX_JOINTS);
    GO_MATRIX_DECLARE(Jinv, Jinv_stg, GENSER_MAX_JOINTS, 6);
    go_pose T_L_0;
    go_real dvw[6];
    go_real jest[GENSER_MAX_JOINTS];
    go_real dj[GENSER_MAX_JOINTS];
    go_pose pest, pestinv, Tdelta; // pos = converted pose from EmcPose
    go_rpy rpy;
    go_rvec rvec;
    go_cart cart;
    go_link linkout[GENSER_MAX_JOINTS];
    int link;
    int smalls;
    int retval;

    // rtapi_print("kineInverse(joints: %f %f %f %f %f %f)\n",
    //      joints[0],joints[1],joints[2],joints[3],joints[4],joints[5]);
    // rtapi_print("kineInverse(world: %f %f %f %f %f %f)\n",
    //      world->tran.x, world->tran.y, world->tran.z, world->a, world->b, world->c);

//    genser_kin_init();
    
    // FIXME-AJ: rpy or zyx ?
    rpy.y = world->c * PM_PI / 180;
    rpy.p = world->b * PM_PI / 180;
    rpy.r = world->a * PM_PI / 180;

    go_rpy_quat_convert(&rpy, &haldata->pos->rot);
    haldata->pos->tran.x = world->tran.x;
    haldata->pos->tran.y = world->tran.y;
    haldata->pos->tran.z = world->tran.z;

    go_matrix_init(Jfwd, Jfwd_stg, 6, genser->link_num);
    go_matrix_init(Jinv, Jinv_stg, genser->link_num, 6);

    /* jest[] is a copy of joints[], which is the joint estimate */
    for (link = 0; link < genser->link_num; link++) {
        // jest, and the rest of joint related calcs are in radians
        jest[link] = joints[link] * (PM_PI / 180);
    }

    for (genser->iterations = 0;
         genser->iterations < *haldata->max_iterations;
         genser->iterations++) {
         *(haldata->last_iterations) = genser->iterations;
        /* update the Jacobians */
        for (link = 0; link < genser->link_num; link++) {
            go_link_joint_set(&genser->links[link], jest[link], &linkout[link]);
        }
        retval = compute_jfwd(linkout, genser->link_num, &Jfwd, &T_L_0);
        if (GO_RESULT_OK != retval) {
            rtapi_print("ERR kI - compute_jfwd (joints: %f %f %f %f %f %f), (iterations=%d)\n",
                 joints[0],joints[1],joints[2],joints[3],joints[4],joints[5], genser->iterations);
            return retval;
        }
        retval = compute_jinv(&Jfwd, &Jinv);
        if (GO_RESULT_OK != retval) {
            rtapi_print("ERR kI - compute_jinv (joints: %f %f %f %f %f %f), (iterations=%d)\n",
                 joints[0],joints[1],joints[2],joints[3],joints[4],joints[5], genser->iterations);
            return retval;
        }

        /* pest is the resulting pose estimate given joint estimate */
        genser_kin_fwd(KINS_PTR, jest, &pest);
        //printf("jest: %f %f %f %f %f %f\n",jest[0],jest[1],jest[2],jest[3],jest[4],jest[5]);
        /* pestinv is its inverse */
        go_pose_inv(&pest, &pestinv);
        /*
            Tdelta is the incremental pose from pest to pos, such that

            0        L         0
            . pest *  Tdelta =  pos, or
            L        L         L

            L         L          0
            .Tdelta =  pestinv *  pos
            L         0          L
        */
        go_pose_pose_mult(&pestinv, haldata->pos, &Tdelta);

        /*
            We need Tdelta in 0 frame, not pest frame, so rotate it
            back. Since it's effectively a velocity, we just rotate it, and
            don't translate it.
        */

        /* first rotate the translation differential */
        go_quat_cart_mult(&pest.rot, &Tdelta.tran, &cart);
        dvw[0] = cart.x;
        dvw[1] = cart.y;
        dvw[2] = cart.z;

        /* to rotate the rotation differential, convert it to a
           velocity screw and rotate that */
        go_quat_rvec_convert(&Tdelta.rot, &rvec);
        cart.x = rvec.x;
        cart.y = rvec.y;
        cart.z = rvec.z;
        go_quat_cart_mult(&pest.rot, &cart, &cart);
        dvw[3] = cart.x;
        dvw[4] = cart.y;
        dvw[5] = cart.z;

        /* push the Cartesian velocity vector through the inverse Jacobian */
        go_matrix_vector_mult(&Jinv, dvw, dj);

        /* check for small joint increments, if so we're done */
        for (link = 0, smalls = 0; link < genser->link_num; link++) {
            if (GO_QUANTITY_LENGTH == linkout[link].quantity) {
            if (GO_TRAN_SMALL(dj[link]))
                smalls++;
            } else {
                if (GO_ROT_SMALL(dj[link]))
                    smalls++;
            }
        }
        if (smalls == genser->link_num) {
            /* converged, copy jest[] out */
            for (link = 0; link < genser->link_num; link++) {
                // convert from radians back to angles
                joints[link] = jest[link] * 180 / PM_PI;
                if ((link) && *(haldata->unrotate[link]))
                    joints[link] += *(haldata->unrotate[link]) * joints[link-1];
            }
            //rtapi_print("DONEkineInverse(joints: %f %f %f %f %f %f), (iterations=%d)\n",
            //     joints[0],joints[1],joints[2],joints[3],joints[4],joints[5], genser->iterations);
            return GO_RESULT_OK;
        }
        /* else keep iterating */
        for (link = 0; link < genser->link_num; link++) {
            jest[link] += dj[link]; //still in radians
        }
    } /* for (iterations) */

    rtapi_print("ERRkineInverse(joints: %f %f %f %f %f %f), (iterations=%d)\n",
         joints[0],joints[1],joints[2],joints[3],joints[4],joints[5], genser->iterations);
    return GO_RESULT_ERROR;
}

// ---------------------------------------------------------------------------------------------------------------

int gensertoolKinematicsInverse(const EmcPose * pos,
		               double *joint,
		               const KINEMATICS_INVERSE_FLAGS * iflags,
		               KINEMATICS_FORWARD_FLAGS * fflags)
{      
    double a1 = 85;
    double a2 = 380;
    double a3 = 100;
    double d1 = 350;
    double d2 = 0;
    double d4 = 425;
    double d6 = 85;
    int tmode = 1;

    double Xt = pos->tran.x;
    double Yt = pos->tran.y;
    double Zt = pos->tran.z;
    double At = pos->a/180*M_PI;
    double Bt = pos->b/180*M_PI;
    double Ct = pos->c/180*M_PI;

    double th1, th2, th3, th4, th5, th6;
    double ux, uy, uz, vx, vy, vz, wx, wy, wz, px, py, pz;
    double r, k1, k2, qx, qy, qz;
    double c1, s1, c2, s2, c3, s3, c4, s4, c5, s5, c23, s23;
    double v114, v124, v214, v224, v113, v313, v323;
    double v111, v131, v311, v331, v411, v431;
    static double th4old = 0.0;
// 
    int n1 = 1; // shoulder on the right
    int n2 = 1; // elbow up
    int n4 = 1; // wrist not flipped

    joint[6] = To[0][3];  joint[7] = To[1][3]; joint[8] = To[2][3]; 

    if (tmode) {
//               tool coordinates
       PoseMat(Xt-To[0][3], Yt-To[1][3], Zt-To[2][3],
               At-Ao, Bt-Bo, Ct-Co, Tt);
       MatMult(To, Tt, Tw);
    } else {
//               world coordinates
       PoseMat(Xt, Yt, Zt, At, Bt, Ct, Tw);
    }

      ux = Tw[0][0]; vx = Tw[0][1]; wx = Tw[0][2]; qx = Tw[0][3];
      uy = Tw[1][0]; vy = Tw[1][1]; wy = Tw[1][2]; qy = Tw[1][3];
      uz = Tw[2][0]; vz = Tw[2][1]; wz = Tw[2][2]; qz = Tw[2][3];

//  wrist position --------------------------------------------------
    px = qx - d6*wx;
    py = qy - d6*wy;
    pz = qz - d6*wz;

//  solve for th1 ---------------------------------------------------
    r = sqrt(px*px + py*py);
    if (r < d2) {
//      'ERROR:--------- point not reachable' 
        return 1;
    }
    k1 = atan2(py, px);
    k2 = asin(d2/r);
    if (n1 == 1) { th1 = k1 + k2;}
    else { th1 = k1 - k2 + M_PI;}
    c1 = cos(th1);  s1 = sin(th1);

//  solve for th2 ---------------------------------------------------
    v114 = px*c1 + py*s1 - a1;
    v124 = pz - d1;
    r = sqrt(v114*v114 + v124*v124);
    k1 = (a2*a2 - d4*d4 - a3*a3 + v114*v114 + v124*v124)/(2*a2*r);
    if (fabs(k1) > 1) {
//       'ERROR:--------- point not reachable'; 
         return 2;
    }
    k2 = acos(k1);
    if (elbup == 1) {n2 = 1;}
    else {n2 = -1;}
    th2 = atan2(v124, v114) + n2*k2;
    c2 = cos(th2); s2 = sin(th2);

//  solve for  th3  -----------------------------------------------
    v214 =  c2*v114 + s2*v124 - a2;
    v224 = -s2*v114 + c2*v124;
    th3 = -atan2(a3, d4) + atan2(v214, -v224);
    c3 = cos(th3); s3 = sin(th3);

//  solve for  th4  -----------------------------------------------
    c23 = cos(th2+th3); s23 = sin(th2+th3);
    v113 = c1*wx + s1*wy;
    v313 = c23*v113 + s23*wz;
    v323 = s1*wx - c1*wy;

    if ((fabs(v323) < eps) && (fabs(v313) < eps)){ th4 = 0;}
    else {th4 = atan2(n4*v323, n4*v313);}
//         take care of singularities and map for continuity   
    if ((fabs(v323) < eps) && (v313 < eps)) {th4 = th4old;}
    if ((v323 > eps) && (v313 < eps)) {th4 = th4 - 2*M_PI;}
    if ((fabs(v113) < eps) && (fabs(v313) < eps) &&
        (fabs(v323) < eps) ) {th4 = th4old;}
    th4old = th4;
//  
    c4 = cos(th4); s4 = sin(th4);

//  solve for  th5  -------------------------------------------------
    k1 = c4*v313 + s4*v323;
    k2 = s23*v113 - c23*wz;
    th5 = atan2(k1, k2);
//  
    c5 = cos(th5); s5 = sin(th5);

//  solve for th6  --------------------------------------------------
    v111 = c1*ux + s1*uy;
    v131 = s1*ux - c1*uy;
    v311 = c23*v111 + s23*uz;
    v331 = s23*v111 - c23*uz;
    v411 =  c4*v311 + s4*v131;
    v431 = -s4*v311 + c4*v131;
// 
    k1 = v431;  
    k2 = c5*v411 - s5*v331;
    th6 = atan2(k1, k2);
//  
// convert to degrees ------------------------------------------------
   joint[0] = th1*180/PM_PI;
   joint[1] = th2*180/PM_PI;
   joint[2] = th3*180/PM_PI;
   joint[3] = th4*180/PM_PI;
   joint[4] = th5*180/PM_PI;
   joint[5] = th6*180/PM_PI;

   return 0;

}

//-------------------------------------------------------------------------------------------------------------------------------
/*
  Extras, not callable using go_kin_ wrapper but if you know you have
  linked in these kinematics, go ahead and call these for your ad hoc
  purposes.
*/

int genser_kin_inv_iterations(genser_struct * genser)
{
    return genser->iterations;
}

int genser_kin_inv_set_max_iterations(int i)
{
    if (i <= 0) return GO_RESULT_ERROR;
    *haldata->max_iterations = i;
    return GO_RESULT_OK;
}

int genser_kin_inv_get_max_iterations()
{
    return *haldata->max_iterations;
}

int genser_hal_setup(int comp_id) {
    int i,res;
    haldata = hal_malloc(sizeof(struct haldata));
    if (!haldata) goto error;

    for (i = 0; i < GENSER_MAX_JOINTS; i++) {
        if ((res =
            hal_pin_float_newf(HAL_IN, &(haldata->a[i]), comp_id,
               "genserkins.A-%d", i)) < 0) goto error;
        *(haldata->a[i])=0;
        if ((res =
            hal_pin_float_newf(HAL_IN, &(haldata->alpha[i]), comp_id,
               "genserkins.ALPHA-%d", i)) < 0) goto error;
        *(haldata->alpha[i])=0;
        if ((res =
            hal_pin_float_newf(HAL_IN, &(haldata->d[i]), comp_id,
               "genserkins.D-%d", i)) < 0) goto error;
        *(haldata->d[i])=0;
        if ((res =
           hal_pin_s32_newf(HAL_IN, &(haldata->unrotate[i]), comp_id,
              "genserkins.unrotate-%d", i)) < 0) goto error;
        *haldata->unrotate[i]=0;
    }
    if ((res=
       hal_pin_u32_newf(HAL_OUT, &(haldata->last_iterations), comp_id,
          "genserkins.last-iterations")) < 0) goto error;


    KINS_PTR = hal_malloc(sizeof(genser_struct));
    haldata->pos = (go_pose *) hal_malloc(sizeof(go_pose));
    if (KINS_PTR == NULL)
        goto error;
    if (haldata->pos == NULL)
        goto error;
    if ((res=
        hal_pin_u32_newf(HAL_IN, &haldata->max_iterations, comp_id,
          "genserkins.max-iterations")) < 0) goto error;

    *haldata->max_iterations = GENSER_DEFAULT_MAX_ITERATIONS;


    A(0) = DEFAULT_A1;
    A(1) = DEFAULT_A2;
    A(2) = DEFAULT_A3;
    A(3) = DEFAULT_A4;
    A(4) = DEFAULT_A5;
    A(5) = DEFAULT_A6;
    ALPHA(0) = DEFAULT_ALPHA1;
    ALPHA(1) = DEFAULT_ALPHA2;
    ALPHA(2) = DEFAULT_ALPHA3;
    ALPHA(3) = DEFAULT_ALPHA4;
    ALPHA(4) = DEFAULT_ALPHA5;
    ALPHA(5) = DEFAULT_ALPHA6;
    D(0) = DEFAULT_D1;
    D(1) = DEFAULT_D2;
    D(2) = DEFAULT_D3;
    D(3) = DEFAULT_D4;
    D(4) = DEFAULT_D5;
    D(5) = DEFAULT_D6;

    return 0;

error:
    return -1;
} // genser_hal_setup()
