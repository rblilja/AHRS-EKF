/**
 *  Copyright 2020 Robin Lilja <robin.lilja@gmail.com>
 *
 *  Created on: Jun 26, 2015
 *      Author: Robin Lilja
 * 		License: MIT License (see LICENSE file)
 */

#ifndef _AHRS_EKF_H_
#define _AHRS_EKF_H_

#include <cmath>        // std::sqrt()

#define AHRS_UNROLLED   ///< Makes the filter to unroll the matrix calculations

namespace AHRS {

static constexpr double AHRS_SAFE_EPSILON = 1e-9;   ///< Make sure we avoid division by zero

/** 
 * Attitude Heading Reference System implemented by an Extended Kalman Filter. 
 */
template <typename T>
class AHRS {

public:

	/**
	 * Constructor.
	 * 
	 * \param Q_attitude attitude quaternion process covariance
	 * \param Q_bias gyro bias process covariance
	 * \param R_acc accelerometer measurement covariance
	 * \param R_mag magnetometer measurement covariance
	 */
	AHRS<T>(const T &Q_attitude, const T &Q_bias, const T &R_acc, const T &R_mag) {

		// Make sure initial quaternion is a valid
		m_q0 = 1.0;
		m_q1 = m_q2 = m_q3 = 0.0;

		// Assume no bias at start
		m_bwx = m_bwy = m_bwz = 0.0;

		m_Qa = Q_attitude;
		m_Qb = Q_bias;

		m_Ra = R_acc;
		m_Rm = R_mag;
	};

	/**
	 * Time update (propagation) of the state estimate.
	 * \param wx x-axis angular rate [rad/s]
	 * \param wy y-axis angular rate [rad/s] 
	 * \param wz z-axis angular rate [rad/s] 
	 * \param dt time delta [s]
	 */
	void propagate(const T &wx, const T &wy, const T &wz, const T &dt);

	/**
	 * Accelerometer measurement update (provides roll and pitch correction).
	 * \param ax x-axis accelerometer measurement [m.s-2]
	 * \param ay y-axis accelerometer measurement [m.s-2]
	 * \param az z-axis accelerometer measurement [m.s-2]
	 */
	void accelerometer_update(T ax, T ay, T az);

	/**
	 * Magnetometer measurement update (provides yaw correction).
	 * \param Bx x-axis magnetometer measurement [uT]
	 * \param By y-axis magnetometer measurement [uT]
	 * \param Bz z-axis magnetometer measurement [uT]
	 */
	void magnetometer_update(T Bx, T By, T Bz);

	/**
	 * Get angular rates (sensor frame).
	 * \param rolldot roll rate [rad/s].
	 * \param pitchdot pitch rate [rad/s].
	 * \param yawdot yaw rate [rad/s].
	 */
	void getRate(T* rolldot, T* pitchdot, T* yawdot) { *rolldot = m_wx; *pitchdot = m_wy; *yawdot = m_wz; };

	void getBias(T* bwx, T* bwy, T* bwz) { *bwx = m_bwx; *bwy = m_bwy; *bwz = m_bwz; };

	void getAttitude(T* roll, T* pitch, T* yaw) { T q[4] = {m_q0, m_q1, m_q2, m_q3}; q2euler(q, roll, pitch, yaw); };

	void getAttitude(T q[4]) { q[0] = m_q0; q[1] = m_q1; q[2] = m_q2; q[3] = m_q3; };

	//void getAttitude(float dcm[3][3]) { float q[4] = {this->q0, this->q1, this->q2, this->q3}; q2dcm(q, dcm);  }

	/**
	 * Attitude state quaternion (q = w + ix + jy + kz).
	 */
	T m_q0, m_q1, m_q2, m_q3;

	/**
	 * True (bias compensated) angular rate (sensor frame) [rad/s].
	 */
	T m_wx, m_wy, m_wz;

	/**
	 * Gyro bias (sensor frame) [rad/s].
	 */
	T m_bwx, m_bwy, m_bwz;

	/**
	 * Process noise covariance of attitude quaternion.
	 * This one more or less sets how much our model deviates from reality.
	 */
	T m_Qa;

	/**
	 * Process noise covariance of gyro bias.
	 * Gyro bias process noise is kept near zero since we don't expect it to vary much.
	 */
	T m_Qb;

	/**
	 * Accelerometer measurement noise covariance.
	 */
	T m_Ra;

	/**
	 * Magnetometer measurement noise covariance.
	 */
	T m_Rm;

private:

#ifndef PI
	static constexpr T PI = T(3.141592653589793);
#endif

	static constexpr T PI_HALF = T(PI / 2.0);

	/**
	 * State estimate covariance matrix 'P'.
	 */
	T m_P[7][7] =
	{
		{ 100,     0.0,      0.0,     0.0,	0.0,	0.0,	0.0 },
		{ 0.0,     100,      0.0,     0.0,	0.0,	0.0,	0.0 },
		{ 0.0,     0.0,      100,     0.0,	0.0,	0.0,	0.0 },
		{ 0.0,     0.0,      0.0,     100,	0.0,	0.0,	0.0 },
		{ 0.0,     0.0,      0.0,     0.0,	1.0,	0.0,	0.0 },
		{ 0.0,     0.0,      0.0,     0.0,	0.0,	1.0,	0.0 },
		{ 0.0,     0.0,      0.0,     0.0,	0.0,	0.0,	1.0 }
	};

	/**
	 * Convert quaternion to Euler angles.
	 * @param q quaternion (q = w + ix + jy + kz).
	 * @param roll [rad].
	 * @param pitch [rad].
	 * @param yaw [rad].
	 */
	void q2euler(const T q[4], T* roll, T* pitch, T* yaw);

	/**
	 * Convert quaternion to DCM.
	 * @param q quternion (q = w + ix + jy + kz).
	 * @param dcm DCM (dcm[row][column]).
	 */
	void q2dcm(const T q[4], T dcm[3][3]);

	/**
	 * Convert DCM to Euler angles.
	 * @param dcm DCM (dcm[row][column]).
	 * @param roll [rad].
	 * @param pitch [rad].
	 * @param yaw [rad].
	 */
	void dcm2euler(const T dcm[3][3], T* roll, T* pitch, T* yaw);
};

template <typename T>
inline void AHRS<T>::q2euler(const T q[4], T* roll, T* pitch, T* yaw) {

    T dcm[3][3];

    q2dcm(q, dcm);
    dcm2euler(dcm, roll, pitch, yaw);
}

template <typename T>
inline void AHRS<T>::q2dcm(const T q[4], T dcm[3][3]) {

	// ZYX rotational order DCM

    // http://www.mathworks.se/help/aeroblks/quaternionstorotationangles.html

    const T q0q0 = q[0] * q[0];
    const T q1q1 = q[1] * q[1];
    const T q2q2 = q[2] * q[2];
    const T q3q3 = q[3] * q[3];

    const T q0q1 = q[0] * q[1];
    const T q0q2 = q[0] * q[2];
    const T q0q3 = q[0] * q[3];
    const T q1q2 = q[1] * q[2];
    const T q1q3 = q[1] * q[3];
    const T q2q3 = q[2] * q[3];

    dcm[0][0] = q0q0 + q1q1 - q2q2 - q3q3;
    dcm[0][1] = 2.0  * (q1q2 + q0q3);
    dcm[0][2] = 2.0  * (q1q3 - q0q2);

    dcm[1][0] = 2.0  * (q1q2 - q0q3);
    dcm[1][1] = q0q0 - q1q1 + q2q2 - q3q3;
    dcm[1][2] = 2.0  * (q2q3 + q0q1);

    dcm[2][0] = 2.0  * (q1q3 + q0q2);
    dcm[2][1] = 2.0  * (q2q3 - q0q1);
    dcm[2][2] = q0q0 - q1q1 - q2q2 + q3q3;
}

template <typename T>
inline void AHRS<T>::dcm2euler(const T dcm[3][3], T* roll, T* pitch, T* yaw) {

    // ZYX rotation order

    // http://www.mathworks.se/help/aeroblks/quaternionstorotationangles.html

    // Roll = phi
    // Pitch = theta
    // Yaw = psi

    static constexpr T EPSILON = 1.0e-6;

    *pitch = std::asin(-dcm[0][2]);

    // If pitch equals +/- PI/2 we have a gimbal lock situation

	// No gimbal lock
	if ( std::abs(std::abs(*pitch) - PI_HALF) > EPSILON ) {

        *roll = std::atan2(dcm[1][2], dcm[2][2]);
        *yaw = std::atan2(dcm[0][1], dcm[0][0]);
    }
    // Gimbal lock (here roll and yaw give the same rotation)
    else {

        // Assume roll to be zero
        *roll = 0.0;

        // Pitch near PI/2
        if (*pitch > 0) { *yaw = -1.0 * std::atan2(-dcm[1][2], dcm[1][1]); }

        // Pitch near -PI/2
        if (*pitch < 0) { *yaw = std::atan2(-dcm[1][2], dcm[1][1]); }
    }
}

} /* namespace */

#endif /* _AHRS_EKF_H_ */
