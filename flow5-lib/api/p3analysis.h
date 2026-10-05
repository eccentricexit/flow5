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



#include <panel3.h>
#include <panelanalysis.h>
#include <trimesh.h>
#include <vector3d.h>


class StabDerivatives;


class FL5LIB_EXPORT P3Analysis : public PanelAnalysis
{
    friend class  Task3d;
    friend class  PlaneTask;
    friend class  BoatTask;

    public:
        P3Analysis();
        virtual ~P3Analysis() = default;
        virtual bool isTriUniMethod() const =0;
        virtual bool isTriLinMethod() const =0;

        void makeInfluenceMatrix() override;
        virtual void makeMatrixBlock(int iBlock) = 0;

        bool initializeAnalysis(const Polar3d *pPolar3d, int nRHS) override;

        void setTriMesh(TriMesh const &trimesh);

        bool isQuadMethod()     const {return false;}
        bool isTriangleMethod() const {return true;}

//        virtual void makeUnitRHSBlock(int iBlock, double *rhs, Vector3d const &Vel, bool bRot) override = 0;

        void releasePanelArrays();

        bool hasWarning()  const {return m_bWarning;}

        int nPanels()      const override {return int(m_Panel3.size());}
        int nWakePanels()  const {return int(m_WakePanel3.size());}
        int nNodes() const {return m_pRefTriMesh->nodeCount();}
        int matSize() const override;

        Panel *panel(int p) override {if(p>=0 && p<nPanels()) return m_Panel3.data()+p; else return nullptr;}
        Panel const *panelAt(int p) const override {if(p>=0 && p<nPanels()) return m_Panel3.data()+p; else return nullptr;}
        Panel3 const*panel3(int p) const {if(p>=0 && p<nPanels()) return m_Panel3.data()+p; else return nullptr;}
        Panel3 const& panel3At(int i) {return m_Panel3.at(i);}
        Panel3 const& wakePanelAt(int iw) {return m_WakePanel3.at(iw);}
        std::vector<Panel3> const &panels() {return m_Panel3;}
        std::vector<Panel3> const &wakePanels() {return m_WakePanel3;}
        void makeNegatingVortices(std::vector<Vortex> &negvortices) override;

        void savePanels() override;
        void restorePanels() override;
        bool computeTrimmedConditions(double mass, Vector3d const &CoG, double &alphaeq, double &u0, bool bFuseMi) override;
        int makeWakePanels(Vector3d const &WindDirection, bool bVortonWake) override;
        bool hasThickSurfaces() const override;
        bool wakeRowExcluded(int i3, int k3) const;
        bool yawedWakeProduct(double theta, std::vector<double> const &mu, std::vector<double> &product) override;
        void scaleResultsToSpeed(double ratio) override;
        void makeMu(int qrhs) override;
        virtual void getDoubletInfluence(Vector3d const &C, Panel3 const &p3, Vector3d *V, double *phi, double coreradius, bool bUseRFF) const;

        void combineLocalVelocities(double alpha, double beta, std::vector<Vector3d> &VLocal) const override;


        void forces(const double *Mu3, const double *Sigma3, double alpha, double beta, const Vector3d &CoG, bool bFuseMi, std::vector<Vector3d> const &VInf, Vector3d &Force, Vector3d &Moment) override;
        void sheetVorticityLoad(double const *Mu3, std::vector<Vector3d> const &VInf, Vector3d const &CoG, double G4pi, Vector3d &Force, Vector3d &Moment) const;
        void inducedForce(int nPanel3, double QInf, double alpha, double beta, int pos3, Vector3d &ForceBodyAxes, SpanDistribs &distribFF) const override;
        void makeTrefftzStrips(double const *mu3, Vector3d const &wakedir, std::vector<TrefftzStrip> &strips) const;
        void trefftzDrag(int nPanel3, double QInf, double alpha, double beta, int pos3, Vector3d &Drag, SpanDistribs &distribFF, bool bTrefftzPlane) const override;

        double getPotential(Vector3d const &C, const double *mu, const double *sigma) const;

        void getVelocityVector(Vector3d const &C, double const *Mu, double const *Sigma, Vector3d &VT, double coreradius, bool bWakeOnly, bool bMultiThread=true) const override;

        void velocityVectorBlock(int iBlock, Vector3d const &C, Vector3d *VT, const VelBlockData &data) const;
        void getFarFieldVelocity(const Vector3d &C, const std::vector<Panel3> &panel3, const double *Mu, Vector3d &VT, double coreradius) const;
        void getDebugPotential(Vector3d const &C, bool bSelf, const double *Mu, const double *Sigma, double &phi, bool bSource=true, bool bDoublet=true, bool bWake=true) const;

        void makePanelDoubletSurfaceVelocity(int p, const double *Mu, Vector3d &VLocal);

        void getSourceInfluence(Polar3d const *pPolar3d, Vector3d const &C, const Panel3 &p3, bool bSelf, Vector3d *V, double *phi) const;
        void makeNodeAverage(int iNode, const std::vector<double> &muPanel, std::vector<double> &muLin) const;

        void makeStripAreas();

        int allocateRHS3(int nRHS);

        void moments(double *Mu3, double alpha, double beta, Vector3d CoG, double *VInf, Vector3d &Force, Vector3d &Moment);

        void sumPanelForces(double *Cp, Vector3d &F);

        bool getZeroMomentAngle(Vector3d const &CoG, double &alphaeq, bool bFuseMi);

        double computeCm(const Vector3d &CoG, double Alpha, bool bFuseMi);
        double stripArea(const Panel3 &p3, bool bThinSurfaces) const;
        int nextTopTrailingPanelIndex(const Panel3 &p3) const;
        Vector3d stripQuarterPoint(Panel3 const &p3) const;

        // Steady rotating states, for the type-8 points with body rates (thin triangle methods; see PlaneTask::T123458Loop()). A state is the onset flow V + (r-CoG) x Omega
        // (V the velocity of the air relative to the plane's CoG, Omega its rotation vector, both in geometry axes), solved with the unit solutions of u, v, w, p, q and r.
        void makeStateSolution(Vector3d const &V, Vector3d const &Omega, Vector3d const &CoG, std::vector<double> &sol) const;
        void makeOnsetField(Vector3d const &V, Vector3d const &Omega, Vector3d const &CoG, std::vector<Vector3d> &field) const;
        void setStateDoublets(std::vector<double> const &sol);
        void makeStateCp(std::vector<double> const &sol, std::vector<Vector3d> const &field, double QInf, double alpha);
        bool makeYawedWakeCorrection(std::vector<double> const &sol0, std::vector<double> &dsol);
        void stateFarField(std::vector<double> const &sol, std::vector<Vector3d> const &field, double alpha, Vector3d const &CoG, bool bFuseMi, bool bSideslipModel,
                           Vector3d &Force, Vector3d &Moment);
        /** The local onset velocities (m/s) of a steady rotating point, per panel, that inducedForce() applies at the strips; empty for a point without rates. */
        void setOnset(std::vector<Vector3d> const &field) {m_Onset=field;}
        void clearOnset() {m_Onset.clear();}

        void midWakePoint(const Panel3 *pWakePanel, Vector3d &midleft, Vector3d &midright) const;
        void trailingWakePoint(const Panel3 *pWakePanel, Vector3d &left, Vector3d &right) const;
        void trailingWakePanels(const Panel3 *pWakePanel, Panel3 &p3WU, Panel3 &p3WD) const;



    protected:
        std::vector<Vector3d> m_Onset;              /**< the velocity of the air relative to each panel's centroid, m/s, of the steady rotating point being solved (setOnset()); empty for a point without rates */
        std::vector<char> m_RowSkip;                /**< the panels left out of the yawed wake's product, set by yawedWakeProduct() */
        std::vector<int> m_Component;               /**< the connected piece (wing) of each panel, set by yawedWakeProduct() */
        double m_WakeClearance = 0.0;               /**< the distance from a wake column within which another wing's panels are left out of the yawed wake's product, set by yawedWakeProduct() */
        TriMesh const *m_pRefTriMesh;
        std::vector<Panel3> m_Panel3;               /**< the panel array for the currently loaded plane */
        std::vector<Panel3> m_refPanel3;
        std::vector<Panel3> m_WakePanel3;           /**< the wake panel array for the currently loaded plane */
        std::vector<Panel3> m_refWakePanel3;

        std::vector<double> m_uRHSVertex, m_vRHSVertex, m_wRHSVertex; /** The unit doublet densities at the triangle's nodes. */
        std::vector<double> m_pRHSVertex, m_qRHSVertex, m_rRHSVertex; /** The unit doublet densities at the triangle's nodes. */

};


