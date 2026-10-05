/****************************************************************************
 *
 * 	flow5 application
 *
 * 	Copyright (C) 2025 André Deperrois 
 *
 * 	
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

 *
 *****************************************************************************/


#pragma once

#include <vector>

#include <vector3d.h>

class WingXfl;

struct FL5LIB_EXPORT SpanDistribs
{
    public:
        SpanDistribs();
        void resizeResults(int NStation);
        void resizeGeometry(int NStation);
        void setGeometry(WingXfl const*pWing);


        void clearGeometry();

        double stripLift(int m, double qDyn) const;
        double stripArea(int m) const {return m_StripArea.at(m);}

        void initializeToZero();

        int nStations() const {return int(m_Cl.size());}

        /** |V|/QInf at the strip: the local speed over the free-stream speed, 1 unless the point rotates (a type-8 point with body rates) */
        double speedRatio(int m) const {return m<int(m_VOnset.size()) ? m_VOnset.at(m).norm() : 1.0;}
        /** The unit vector along the local velocity at the strip, or `freestream` unless the point rotates */
        Vector3d flowDirection(int m, Vector3d const &freestream) const {return m<int(m_VOnset.size()) ? m_VOnset.at(m)*(1.0/m_VOnset.at(m).norm()) : freestream;}

    public:
        std::vector<double> m_Ai;            /**< the induced angles, in degrees */
        std::vector<double> m_Alpha_0;       /**< the zero-lift angle at the span station */
        std::vector<double> m_Cl;            /**< the lift coefficient on the strips */
        std::vector<double> m_ICd;           /**< the induced drag coefficient on the strips */
        std::vector<double> m_PCd;           /**< the viscous drag coefficient on the strips */
        std::vector<double> m_Re;            /**< the Reynolds number on the strips */
        std::vector<double> m_XTrTop;        /**< the upper transition location on the strips */
        std::vector<double> m_XTrBot;        /**< the lower transition location on the strips */
        std::vector<double> m_CmViscous;     /**< the pitching moment of viscous drag on the strips, w.r.t CoG, normalized by the strip's chord and area  */
        std::vector<double> m_CmPressure;    /**< the pitching moment of the pressure forces on the strips, w.r.t CoG, normalized by the strip's chord and area  */
        std::vector<double> m_CmC4;          /**< the pitching moment coefficient on the strips w.r.t. the chord's quarter point, normalized by the strip's chord and area */
        std::vector<double> m_XCPSpanRel;    /**< the relative position of the strip's center of pressure on the strips as a % of the local chord length*/
        std::vector<double> m_XCPSpanAbs;    /**< the absolute position of the strip's center of pressure pos on the strips */
        std::vector<double> m_BendingMoment; /**< the bending moment on the strips */
        std::vector<double> m_VTwist;        /**< the virtual twist in viscous loops */
        std::vector<double> m_Gamma;         /**< the circulation on the strip */
        std::vector<bool> m_bConverged;      /**< true if the local viscous interpolation or OTF calculation has converged */
        std::vector<Vector3d> m_Vd;          /**< the downwash vector at span stations in m/s. The downwash is calculated at the mid wake point, i.e. where the induced drag is evaluated. */
        std::vector<Vector3d> m_F;           /**< the force vector at span stations, in N and in body axes */
        std::vector<Vector3d> m_FInduced;    /**< the Trefftz-plane force on each strip, in N and in body axes, which m_F leaves out: m_F + m_FInduced is the strip's far-field force. Zero with a vorton wake, where m_F includes it. Set by the triangle methods only. */
        std::vector<Vector3d> m_FPressure;   /**< the sum of the panel forces on each strip, in N and in body axes: the near-field counterpart of m_F (panel pressures; for VLM, the panels' vortex forces). The tip patches' forces go to the nearest strip, so that the strips sum to the wing's Fsum. */


        std::vector<double> m_Chord;         /**< the chord on the strips */
        std::vector<double> m_Offset;        /**< the offset at the span stations */
        std::vector<double> m_Twist;         /**< the twist at the span stations */
        std::vector<double> m_StripArea;     /**< the area of each chordwise strip */
        std::vector<double> m_StripPos;       /**< the span positions of the stations */
        std::vector<Vector3d> m_PtC4;        /**< the quarter chord points */
        std::vector<Vector3d> m_VOnset;      /**< the velocity of the air relative to the strip's trailing panel (free stream plus the plane's rotation), over the free-stream speed, with body rates only: then m_Re, m_Cl and the viscous drag are those of the local speed. Empty otherwise. */
        std::vector<Vector3d> m_PtLE;        /**< the midpoints of the strips' leading edges, on the panel mesh. Set by the triangle methods for thin surfaces only. */
};

