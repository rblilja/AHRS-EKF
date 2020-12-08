/**
 *  Copyright 2020 Robin Lilja <robin.lilja@gmail.com>
 *
 *  Created on: Jun 26, 2015
 *      Author: Robin Lilja
 * 		License: MIT License (see LICENSE file)
 */

#include "ahrs.h"

namespace AHRS {

template <typename T>
void AHRS<T>::propagate(const T &wx, const T &wy, const T &wz, const T &dt) {

	// State vector is defined as x = [q0 q1 q2 q3 bwx bwy bwz]'

	// Body rates are represented as quaternion time derivative 'qdot' by the following
	// d(quaternion)/dt = 1/2 * [  0 -wx -wy -wz ][q0]
	// 							[ wx   0  wz -wy ][q1]
  	// 							[ wy -wz   0  wx ][q2]
  	// 							[ wz  wy -wx   0 ][q3]
	//
	
	// If the sampling rate is sufficiently high, the quaternion dynamics behave in a
	// quasi-linear fashion since, with small timesteps, the integration steps propagate the quaternions
	// only small deviations away from the unit sphere, making the error in linearization minimal.

	// Remove bias from rate measurements
	m_wx = wx - m_bwx;
	m_wy = wy - m_bwy;
	m_wz = wz - m_bwz;

	// Calculate the linearized sate model 'F_k-1 = d(xdot)/dx' around the current estimated state 'xhat'
	//
	// Repeated arithmetics for the linearized state model
	const T _wx2 = 0.5 * m_wx * dt;
	const T _wy2 = 0.5 * m_wy * dt;
	const T _wz2 = 0.5 * m_wz * dt;
	//
	const T _q02 = 0.5 * m_q0 * dt;
	const T _q12 = 0.5 * m_q1 * dt;
	const T _q22 = 0.5 * m_q2 * dt;
	const T _q32 = 0.5 * m_q3 * dt;
	//
	// Discretization by (I + F*dt) is included in the 'F_k-1' matrix for convenience. Note that angular rates has been compensated for bias already.
	const T F[7][7] =
	{
		{ 1.0, -_wx2, -_wy2, -_wz2, _q12, _q22, _q32 },
		{ _wx2,	1.0, _wz2, -_wy2, -_q02, _q32, -_q22 },
		{ _wy2,	-_wz2, 1.0, _wx2, -_q32, -_q02, _q12 },
		{ _wz2,	_wy2, -_wx2, 1.0, _q22, -_q12, -_q02 },
		{ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 }
	};

	T FP[7];		// Intermediate array holding the i:th row of the 'F_k-1 * P_k-1|k-1' multiplication product
	T P_k[7][7];	// The new state estimate covariance matrix 'P_k|k-1'

	// Calculate the state estimate covariance 'P_k|k-1 = F_k-1 * P_k-1|k-1 * F'_k-1 + Q_k-1'
	//
	// Begin with the multiplications
	for (uint8_t i=0; i<7; ++i) {

		// Calculate i:th row of the 'F_k-1 * P_k-1|k-1' multiplication product
		for (uint8_t j=0; j<7; ++j) {

#ifdef AHRS_UNROLLED

			// Exploit the fact that 'F_k-1' contains quite a few zeros, and some ones on the diagonal
			switch (i) {
				case 0: FP[j] =         m_P[0][j] + F[i][1]*m_P[1][j] + F[i][2]*m_P[2][j] + F[i][3]*m_P[3][j] + F[i][4]*m_P[4][j] + F[i][5]*m_P[5][j] + F[i][6]*m_P[6][j]; break;
				case 1: FP[j] = F[i][0]*m_P[0][j] +         m_P[1][j] + F[i][2]*m_P[2][j] + F[i][3]*m_P[3][j] + F[i][4]*m_P[4][j] + F[i][5]*m_P[5][j] + F[i][6]*m_P[6][j]; break;
				case 2: FP[j] = F[i][0]*m_P[0][j] + F[i][1]*m_P[1][j] +         m_P[2][j] + F[i][3]*m_P[3][j] + F[i][4]*m_P[4][j] + F[i][5]*m_P[5][j] + F[i][6]*m_P[6][j]; break;
				case 3: FP[j] = F[i][0]*m_P[0][j] + F[i][1]*m_P[1][j] + F[i][2]*m_P[2][j] +         m_P[3][j] + F[i][4]*m_P[4][j] + F[i][5]*m_P[5][j] + F[i][6]*m_P[6][j]; break;
				case 4: FP[j] =                                                                         		m_P[4][j];                                     			   break;
				case 5: FP[j] =                                                                                          		  	m_P[5][j];                   		   break;
				case 6: FP[j] =                                                                                                             			m_P[6][j];
			}
#elif
			// Without exploitation of zeros the following unrolled calculation shall be used
			FP[j] = F[i][0]*m_P[0][j] + F[i][1]*m_P[1][j] + F[i][2]*m_P[2][j] + F[i][3]*m_P[3][j] + F[i][4]*m_P[4][j] + F[i][5]*m_P[5][j] + F[i][6]*m_P[6][j];
#endif
		}

		// Calculate the i:th row of 'P_k|k-1' by multiplying the intermediate i:th row with the transpose of 'F_k-1'
		for (uint8_t j=0; j<7; ++j) {

#ifdef AHRS_UNROLLED

			// As for 'F_k-1' we exploit the fact that its transpose contains zeros and ones as well
			// Note that the transpose is not necessary to calculate explicitly, instead of multiplying the i:th row with each column of 'F'_k-1' (the transpose) we
			// simply multiply with the rows of 'F_k-1' instead
			switch (j) {
				case 0: P_k[i][j] = FP[0]         + FP[1]*F[j][1] + FP[2]*F[j][2] + FP[3]*F[j][3] + FP[4]*F[j][4] + FP[5]*F[j][5] + FP[6]*F[j][6]; break;
				case 1: P_k[i][j] = FP[0]*F[j][0] + FP[1]         + FP[2]*F[j][2] + FP[3]*F[j][3] + FP[4]*F[j][4] + FP[5]*F[j][5] + FP[6]*F[j][6]; break;
				case 2: P_k[i][j] = FP[0]*F[j][0] + FP[1]*F[j][1] + FP[2]         + FP[3]*F[j][3] + FP[4]*F[j][4] + FP[5]*F[j][5] + FP[6]*F[j][6]; break;
				case 3: P_k[i][j] = FP[0]*F[j][0] + FP[1]*F[j][1] + FP[2]*F[j][2] + FP[3]         + FP[4]*F[j][4] + FP[5]*F[j][5] + FP[6]*F[j][6]; break;
				case 4: P_k[i][j] =                                                                 FP[4];		                                   break;
				case 5: P_k[i][j] =                                                                                 FP[5];		                   break;
				case 6: P_k[i][j] =                                                                                                 FP[6];
			}
#elif
			// Without exploitation of zeros the following unrolled calculation shall be used
			P_k[i][j] = FP[0]*F[j][0] + FP[1]*F[j][1] + FP[2]*F[j][2] + FP[3]*F[j][3] + FP[4]*F[j][4] + FP[5]*F[j][5] + FP[6]*F[j][6];
#endif
		}
	}
	//
	// Finalize 'P_k|k-1' by adding process noise, for the sake of computational simplicity only a diagonal noise matrix 'Q_k-1' is considered
	//
	// Some repeated arithmetics
	T _Qa = 0.5 * dt * m_Qa;
	T _Qb = 0.5 * dt * m_Qb;
	//
	P_k[0][0] += _Qa;
	P_k[1][1] += _Qa;
	P_k[2][2] += _Qa;
	P_k[3][3] += _Qa;
	P_k[4][4] += _Qb;
	P_k[5][5] += _Qb;
	P_k[6][6] += _Qb;

	// Take a copy of the state estimate covariance matrix
	memcpy(m_P, P_k, sizeof(T)*49);

	// Calculate the state derivative 'xdot', effectively being the rate quaternion 'qdot'
	//
	T qdot0, qdot1, qdot2, qdot3;
	//
	qdot0 = 0.5 * (-m_q1 * m_wx - m_q2 * m_wy - m_q3 * m_wz );
	qdot1 = 0.5 * ( m_q0 * m_wx + m_q2 * m_wz - m_q3 * m_wy );
	qdot2 = 0.5 * ( m_q0 * m_wy - m_q1 * m_wz + m_q3 * m_wx );
	qdot3 = 0.5 * ( m_q0 * m_wz + m_q1 * m_wy - m_q2 * m_wx );

	// Propagate the 'xhat' state estimate 'dt' seconds from time instance 'k-1' to 'k'
	//
	// Euler forward integration
	m_q0 += qdot0 * dt;
	m_q1 += qdot1 * dt;
	m_q2 += qdot2 * dt;
	m_q3 += qdot3 * dt;
	//
	// Normalize quaternion
	T reciprocal = 1.0 / std::sqrt(m_q0 * m_q0 + m_q1 * m_q1 + m_q2 * m_q2 + m_q3 * m_q3);
	m_q0 *= reciprocal;
	m_q1 *= reciprocal;
	m_q2 *= reciprocal;
	m_q3 *= reciprocal;
}

template <typename T>
void AHRS<T>::accelerometer_update(T ax, T ay, T az) {

	// Normalize accelerometer measurement
	const T reciprocal = 1.0 / std::sqrt(ax * ax + ay * ay + az * az);
	ax *= reciprocal;
	ay *= reciprocal;
	az *= reciprocal;

	// Repeated arithmetics used the filter update
	const T _2q0 = 2.0 * m_q0;
	const T _2q1 = 2.0 * m_q1;
	const T _2q2 = 2.0 * m_q2;
	const T _2q3 = 2.0 * m_q3;
	//
	const T _q0q0 = m_q0 * m_q0;
	const T _q1q1 = m_q1 * m_q1;
	const T _q2q2 = m_q2 * m_q2;
	const T _q3q3 = m_q3 * m_q3;
	//
	const T _2q1q0 = 2.0 * (m_q1 * m_q0);
	//const T _2q1q2 = 2.0 * (m_q1 * m_q2);
	const T _2q1q3 = 2.0 * (m_q1 * m_q3);
	//
	const T _2q2q0 = 2.0 * (m_q2 * m_q0);
	//const T _2q3q0 = 2.0 * (m_q3 * m_q0);
	//
	const T _2q2q3 = 2.0 * (m_q2 * m_q3);

	// 'H = d(xhat)/dx'
	const T H[6][7] =
	{
		{ -_2q2,	_2q3,		-_2q0,		_2q1,		0.0,	0.0,	0.0 },
		{ _2q1,		_2q0,		_2q3,		_2q2,		0.0,	0.0,	0.0 },
		{ _2q0,		-_2q1,		-_2q2,		_2q3,		0.0,	0.0,	0.0 },
		{ _2q0,		_2q1,		-_2q2,		-_2q3,		0.0,	0.0, 	0.0 },
		{ -_2q3,	_2q2,		_2q1,		-_2q0,		0.0,	0.0, 	0.0 },
		{ _2q2,		_2q3,		_2q0,		_2q1,		0.0,	0.0, 	0.0 }
	};

	// // Magnetometer measurement is present
	// if ((sensors & MAGNETOMETER)  == MAGNETOMETER) {

	// 	// Transform measured body frame magnetic vector into Earth (NED) frame, the z-axis is neglected
	// 	// since the magnetic "heading" vector in the Earth (NED) frame is always in the xy-plane i.e. north and east directions
	// 	//
	// 	// Be = DCM' * Bb
	// 	T Bex = (_q0q0 + _q1q1 - _q2q2 - _q3q3)*Bx + (_2q1q2 - _2q3q0)*By + (_2q1q3 + _2q2q0)*Bz;
	// 	T Bey = (_2q1q2 + _2q3q0)*Bx + (_q0q0 - _q1q1 + _q2q2 - _q3q3)*By + (_2q2q3 - _2q1q0)*Bz;

	// 	// Normalize the magnetic vector, it will now be a unit vector in the horizontal plane of the NED frame
	// 	reciprocal = 1.0 / sqrt(Bex * Bex + Bey * Bey);
	// 	Bex *= reciprocal;
	// 	Bey *= reciprocal;

	// 	// Transform the vector back to the body frame
	// 	//
	// 	// Bb = DCM * Be
	// 	Bx = (_q0q0 + _q1q1 - _q2q2 - _q3q3)*Bex + (_2q1q2 + _2q3q0)*Bey;
	// 	By = (_2q1q2 - _2q3q0)*Bex + (_q0q0 - _q1q1 + _q2q2 - _q3q3)*Bey;
	// 	Bz = (_2q1q3 + _2q2q0)*Bex + (_2q2q3 - _2q1q0)*Bey;
	// }

	// Innovation vector 'y' being the difference between actual observation vector 'z' and
	// recreation of observation vector 'zhat' from the current state estimate 'xhat'
	//
	// The following rotation matrix (DCM) transform a vector in the Earth (NED) coordinate frame into the body frame of reference
	//
	// [ (q0*q0 + q1*q1 - q2*q2 - q3*q3)	2*(q1*q2 + q3*q0)					2*(q1*q3 - q2*q0)		    	]
	// [ 2*(q1*q2 - q3*q0)    				(q0*q0 - q1*q1 + q2*q2 - q3*q3) 	2*(q2*q3 + q1*q0)				]
	// [ 2*(q1*q3 + q2*q0)    				2*(q2*q3 - q1*q0)					(q0*q0 - q1*q1 - q2*q2 + q3*q3)	]
	//
	// In the Earth (NED) frame two reference vectors exist, one being the acceleration vector [0 0 1]' pointing at Earth's center
	// and the second one being the magnetic field vector [1 0 0]' pointing at Earth's magnetic north.
	//
	// Our task here is to recreate those vectors from the current estimated state and express them in the body frame of reference.
	//
	// When the rotation matrix is used for transforming the NED gravitation and magnetic vectors into the body frame, the former only receives
	// contribution from the third column, and the latter only from the first column.
	//
	const T y[3] =
	{
		ax - (_2q1q3 - _2q2q0),
		ay - (_2q2q3 + _2q1q0),
		az - (_q0q0 - _q1q1 - _q2q2 + _q3q3)
		// Bx - (_q0q0 + _q1q1 - _q2q2 - _q3q3),
		// By - (_2q1q2 - _2q3q0),
		// Bz - (_2q1q3 + _2q2q0)
	};

	T S;		// The innovation covariance 'S_k'
	T Sinv;		// The inverse of the innovation covariance 'S_k'
	T K[7];		// The Kalman gain 'K_k'
	T HP[7];	// Used for calculating both 'S_k' and 'K_k'

	// Sequential update, see references: 
    // 1) KALMAN FILTERING Theory and Practice Using MATLAB, 4.2.2, p. 141, 3:rd
    // 2) http://www.anuncommonlab.com/articles/how-kalman-filters-work/part2.html#speed
	for (uint8_t m=0; m<3; ++m) {

		// Calculate the m:th row of the 'H_k * P_k|k-1' multiplication product
		//
		// Iterate the columns of 'P_k|k-1'
		for (uint8_t i=0; i<7; ++i) {

			HP[i] = H[m][0]*m_P[0][i] + H[m][1]*m_P[1][i] + H[m][2]*m_P[2][i] + H[m][3]*m_P[3][i];
		}

		// Calculate the innovation covariance 'S_k = H_k * P_k|k-1 * H'_k + R_k'
		S = HP[0]*H[m][0] + HP[1]*H[m][1] + HP[2]*H[m][2] + HP[3]*H[m][3] + m_Ra;
		//
		// Multiplying is faster than dividing (so divide once)
		Sinv = 1.0 / (S + AHRS_SAFE_EPSILON);

		// Calculate the Kalman gain 'K_k = (H_k * P_k|k-1) / S
		for (uint8_t i=0; i<7; ++i) {

			K[i] = HP[i]*Sinv;
		}

		// Update the state estimate covariance 'P_k|k = P_k|k-1 - K_k * H_k * P_k|k-1'
		for (uint8_t i=0; i<7; ++i) {

			// Exploit mirrored values over the diagonal
			for (uint8_t j=i; j<7; ++j) {

				m_P[i][j] = m_P[j][i] = m_P[i][j] - K[i]*HP[j];
			}
		}

		// Update the state estimate 'xhat_k|k = xhat_k|k-1 + K_k * y_k'
		m_q0 += K[0] * y[m];
		m_q1 += K[1] * y[m];
		m_q2 += K[2] * y[m];
		m_q3 += K[3] * y[m];
		m_bwx += K[4] * y[m];
		m_bwy += K[5] * y[m];
		m_bwz += K[6] * y[m];
	}
}

template <typename T>
void AHRS<T>::magnetometer_update(T Bx, T By, T Bz) {

	// Repeated arithmetics used the filter update
	const T _2q0 = 2.0 * m_q0;
	const T _2q1 = 2.0 * m_q1;
	const T _2q2 = 2.0 * m_q2;
	const T _2q3 = 2.0 * m_q3;
	//
	const T _q0q0 = m_q0 * m_q0;
	const T _q1q1 = m_q1 * m_q1;
	const T _q2q2 = m_q2 * m_q2;
	const T _q3q3 = m_q3 * m_q3;
	//
	const T _2q1q0 = 2.0 * (m_q1 * m_q0);
	const T _2q1q2 = 2.0 * (m_q1 * m_q2);
	const T _2q1q3 = 2.0 * (m_q1 * m_q3);
	//
	const T _2q2q0 = 2.0 * (m_q2 * m_q0);
	const T _2q3q0 = 2.0 * (m_q3 * m_q0);
	//
	const T _2q2q3 = 2.0 * (m_q2 * m_q3);

	// 'H = d(xhat)/dx'
	const T H[6][7] =
	{
		{ -_2q2,	_2q3,		-_2q0,		_2q1,		0.0,	0.0,	0.0 },
		{ _2q1,		_2q0,		_2q3,		_2q2,		0.0,	0.0,	0.0 },
		{ _2q0,		-_2q1,		-_2q2,		_2q3,		0.0,	0.0,	0.0 },
		{ _2q0,		_2q1,		-_2q2,		-_2q3,		0.0,	0.0, 	0.0 },
		{ -_2q3,	_2q2,		_2q1,		-_2q0,		0.0,	0.0, 	0.0 },
		{ _2q2,		_2q3,		_2q0,		_2q1,		0.0,	0.0, 	0.0 }
	};

	// Transform measured body frame magnetic vector into Earth (NED) frame, the z-axis is neglected
	// since the magnetic heading vector in the Earth (NED) frame is always in the xy-plane i.e. north and east directions
	//
	// Be = DCM' * Bb
	T Bex = (_q0q0 + _q1q1 - _q2q2 - _q3q3) * Bx + (_2q1q2 - _2q3q0) 			   * By + (_2q1q3 + _2q2q0) * Bz;
	T Bey = (_2q1q2 + _2q3q0)				* Bx + (_q0q0 - _q1q1 + _q2q2 - _q3q3) * By + (_2q2q3 - _2q1q0) * Bz;

	// Normalize the magnetic vector, it will now be a unit vector in the horizontal plane of the NED frame
	T reciprocal = 1.0 / sqrt(Bex * Bex + Bey * Bey);
	Bex *= reciprocal;
	Bey *= reciprocal;

	// Transform the vector back to the body frame
	//
	// Bb = DCM * Be
	Bx = (_q0q0 + _q1q1 - _q2q2 - _q3q3) * Bex + (_2q1q2 + _2q3q0) 				 * Bey;
	By = (_2q1q2 - _2q3q0)			     * Bex + (_q0q0 - _q1q1 + _q2q2 - _q3q3) * Bey;
	Bz = (_2q1q3 + _2q2q0)				 * Bex + (_2q2q3 - _2q1q0)				 * Bey;

	// Normalize accelerometer measurement
	reciprocal = 1.0 / std::sqrt(Bx * Bx + By * By + Bz * Bz);
	Bx *= reciprocal;
	By *= reciprocal;
	Bz *= reciprocal;
	
	// Innovation vector 'y' being the difference between actual observation vector 'z' and
	// recreation of observation vector 'zhat' from the current state estimate 'xhat'
	//
	// The following rotation matrix (DCM) transform a vector in the Earth (NED) coordinate frame into the body frame of reference
	//
	// [ (q0*q0 + q1*q1 - q2*q2 - q3*q3)	2*(q1*q2 + q3*q0)					2*(q1*q3 - q2*q0)		    	]
	// [ 2*(q1*q2 - q3*q0)    				(q0*q0 - q1*q1 + q2*q2 - q3*q3) 	2*(q2*q3 + q1*q0)				]
	// [ 2*(q1*q3 + q2*q0)    				2*(q2*q3 - q1*q0)					(q0*q0 - q1*q1 - q2*q2 + q3*q3)	]
	//
	// In the Earth (NED) frame two reference vectors exist, one being the acceleration vector [0 0 1]' pointing at Earth's center
	// and the second one being the magnetic field vector [1 0 0]' pointing at Earth's magnetic north.
	//
	// Our task here is to recreate those vectors from the current estimated state and express them in the body frame of reference.
	//
	// When the rotation matrix is used for transforming the NED gravitation and magnetic vectors into the body frame, the former only receives
	// contribution from the third column, and the latter only from the first column.
	//
	const T y[6] =
	{
		0,
		0,
		0,
		Bx - (_q0q0 + _q1q1 - _q2q2 - _q3q3),
		By - (_2q1q2 - _2q3q0),
		Bz - (_2q1q3 + _2q2q0)
	};

	T S;		// The innovation covariance 'S_k'
	T Sinv;		// The inverse of the innovation covariance 'S_k'
	T K[7];		// The Kalman gain 'K_k'
	T HP[7];	// Used for calculating both 'S_k' and 'K_k'

	// Sequential update, see references: 
    // 1) KALMAN FILTERING Theory and Practice Using MATLAB, 4.2.2, p. 141, 3:rd
    // 2) http://www.anuncommonlab.com/articles/how-kalman-filters-work/part2.html#speed
	for (uint8_t m=3; m<6; ++m) {

		// Calculate the m:th row of the 'H_k * P_k|k-1' multiplication product
		//
		// Iterate the columns of 'P_k|k-1'
		for (uint8_t i=0; i<7; ++i) {

			HP[i] = H[m][0]*m_P[0][i] + H[m][1]*m_P[1][i] + H[m][2]*m_P[2][i] + H[m][3]*m_P[3][i];
		}

		// Calculate the innovation covariance 'S_k = H_k * P_k|k-1 * H'_k + R_k'
		S = HP[0]*H[m][0] + HP[1]*H[m][1] + HP[2]*H[m][2] + HP[3]*H[m][3] + m_Ra;
		//
		// Multiplying is faster than dividing (so divide once)
		Sinv = 1.0 / (S + AHRS_SAFE_EPSILON);

		// Calculate the Kalman gain 'K_k = (H_k * P_k|k-1) / S
		for (uint8_t i=0; i<7; ++i) {

			K[i] = HP[i]*Sinv;
		}

		// Update the state estimate covariance 'P_k|k = P_k|k-1 - K_k * H_k * P_k|k-1'
		for (uint8_t i=0; i<7; ++i) {

			// Exploit mirrored values over the diagonal
			for (uint8_t j=i; j<7; ++j) {

				m_P[i][j] = m_P[j][i] = m_P[i][j] - K[i]*HP[j];
			}
		}

		// Update the state estimate 'xhat_k|k = xhat_k|k-1 + K_k * y_k'
		m_q0 += K[0] * y[m];
		m_q1 += K[1] * y[m];
		m_q2 += K[2] * y[m];
		m_q3 += K[3] * y[m];
		m_bwx += K[4] * y[m];
		m_bwy += K[5] * y[m];
		m_bwz += K[6] * y[m];
	}
}

// Explicit declarations
template class AHRS<float>;
template class AHRS<double>;

} /* namespace */

//	const float DECLINATION = 5.0f*M_PI/180.0f;
//
//	const float cos = cosf(DECLINATION);
//	const float sin = sinf(DECLINATION);
//
//	const float H[6][7] =
//	{
//		{ -_2q2,				_2q3,						-_2q0,						_2q1,						0.0,	0.0,	0.0 },
//		{ _2q1,					_2q0,						_2q3,						_2q2,						0.0,	0.0,	0.0 },
//		{ 0.0,					-2.0*_2q1,					-2.0*_2q2,					0.0,						0.0,	0.0,	0.0 },
//		{ -_2q3*sin,			-_2q2*sin,					-2.0*_2q2*cos - _2q1*sin,	-2.0*_2q3*cos - _2q0*sin,	0.0,	0.0, 	0.0 },
//		{ -_2q3*cos,			_2q2*cos + 2.0*_2q1*sin,	_2q1*cos,					-_2q0*cos + 2.0*_2q3*sin,	0.0,	0.0, 	0.0 },
//		{ _2q2*cos + _2q1*sin,	_2q3*cos + _2q0*sin,		_2q0*cos - _2q3*sin,		_2q1*cos - _2q2*sin,		0.0,	0.0, 	0.0 }
//	};
// Magnetic declination included
//	const float y[6] =
//	{
//		ax - (2.0f * (q1*q3 - q2*q0)),
//		ay - (2.0f * (q2*q3 + q1*q0)),
//		az - (1.0f - 2.0f * (q1*q1 + q2*q2)),
//
//		Bx - (1.0f - 2.0f * (q2*q2 + q3*q3))*cos - (2.0f * (q1*q2 - q3*q0))*sin,
//		By - (2.0f * (q1*q2 - q3*q0))*cos - (1.0f - 2.0f * (q1*q1 + q2*q2))*sin,
//		Bz - (2.0f * (q1*q3 + q2*q0))*cos - (2.0f * (q2*q3 + q1*q0))*sin
//	};
