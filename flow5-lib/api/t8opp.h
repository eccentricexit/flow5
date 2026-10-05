/****************************************************************************

    flow5 application
    Copyright (C) 2025 André Deperrois 
    
    This file is part of flow5.

    flow5 is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    flow5 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with flow5.
    If not, see <https://www.gnu.org/licenses/>.


*****************************************************************************/


#pragma once



/**
 * An operating point of a type-8 polar: the angle of attack and the sideslip (degrees), the speed (m/s) and, optionally, the steady body rates
 * p, q and r (rad/s) of the plane about its CoG, in the stability axes of the derivatives (x forward along the projection of the flight path on the
 * plane of symmetry, y to the right, z down): p positive right wing down, q positive nose up, r positive nose to the right. The rates default to 0,
 * a straight flight point. With rates the point is the steady rotating state of the plane (a turn, a pull-up): see PlaneTask::T123458Loop().
 */
struct T8Opp
{
        T8Opp() : m_bActive(true), m_Alpha(0), m_Beta(0), m_Vinf(1.0), m_P(0), m_Q(0), m_R(0)
        {
        }

        T8Opp(bool bActive, double alpha, double beta, double vinf, double p=0.0, double q=0.0, double r=0.0) :
            m_bActive(bActive), m_Alpha(alpha), m_Beta(beta), m_Vinf(vinf), m_P(p), m_Q(q), m_R(r)
        {
        }

        bool isActive() const {return m_bActive;}
        double alpha() const {return m_Alpha;}
        double beta() const {return m_Beta;}
        double Vinf() const {return m_Vinf;}
        double p() const {return m_P;}   /**< roll rate, rad/s, stability axes, positive right wing down */
        double q() const {return m_Q;}   /**< pitch rate, rad/s, stability axes, positive nose up */
        double r() const {return m_R;}   /**< yaw rate, rad/s, stability axes, positive nose to the right */
        bool hasRates() const {return m_P!=0.0 || m_Q!=0.0 || m_R!=0.0;}

        bool m_bActive;
        double m_Alpha;
        double m_Beta;
        double m_Vinf;
        double m_P;
        double m_Q;
        double m_R;
};
