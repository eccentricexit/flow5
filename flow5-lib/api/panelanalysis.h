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


#include <vorton.h>
#include <vortex.h>
#include <aeroforces.h>
#include <spandistribs.h>
#include <utils.h>
#include <enums_objects.h>

class Polar3d;
class Panel;
class StabDerivatives;
class Vortex;


struct VelBlockData
{
    double coreradius = 0.0;
    bool bWakeOnly = false;
    double const *Mu = nullptr;
    double const *Sigma = nullptr;
    xfl::enumAnalysisMethod m_AnalysisMethod = xfl::NOMETHOD;
};


/** A trailing strip as seen in the Trefftz plane: the bound vortex runs from A to B with circulation Gamma (N/rho/V), the trailing legs
 *  at A and B run along the wake direction (the mesh x axis), so that in the plane they are two point vortices of opposite sign. */
struct TrefftzStrip
{
    int iPanel = -1;        /**< the index of the strip's trailing panel */
    Vector3d A, B;          /**< the ends of the bound vortex; projected in the Trefftz plane by projectTrefftzStrips() */
    double Gamma = 0.0;
};


class FL5LIB_EXPORT PanelAnalysis
{
    friend class  Task3d;
    friend class  PlaneTask;
    friend class  BoatTask;

    public:
        PanelAnalysis();
        virtual ~PanelAnalysis() = default;

        void setAnalysisStatus(xfl::enumAnalysisStatus status) {m_AnalysisStatus=status;}
        void cancelAnalysis() {m_AnalysisStatus=xfl::CANCELLED;}
        bool isCancelled() const {return m_AnalysisStatus==xfl::CANCELLED;}
        bool isRunning()   const {return m_AnalysisStatus==xfl::RUNNING;}
        bool isPending()   const {return m_AnalysisStatus==xfl::PENDING;}
        bool isFinished()  const {return m_AnalysisStatus==xfl::FINISHED || m_AnalysisStatus==xfl::CANCELLED;}

        bool allocateMatrix(int N);
        void releaseMatrix();

        virtual Panel *panel(int p) = 0;
        virtual Panel const *panelAt(int p) const = 0;
        virtual bool initializeAnalysis(Polar3d const *pPolar3d, int nRHS);
        virtual int  matSize() const = 0;
        virtual int  makeWakePanels(Vector3d const &WindDirection, bool bVortonWake) = 0;
        virtual void scaleResultsToSpeed(double ratio) = 0;
        virtual void savePanels() = 0;
        virtual void restorePanels() = 0;
        virtual void makeInfluenceMatrix() = 0;
        virtual void makeUnitRHSBlock(int iBlock) = 0;
        virtual void makeRHSBlock(int iBlock, double *RHS, std::vector<Vector3d> const &VField, Vector3d const*normals) const = 0;
        virtual void makeUnitDoubletStrengths(double alpha, double beta) = 0;
        virtual void makeWakeMatrixBlock(int iBlock) = 0;
        /** Sets product = (dA_wake/dtheta).mu: the change of the wake's rows of the influence matrix when the wake, which leaves the trailing edges along the
         *  body x axis, is yawed by theta (rad, positive towards +y), times the doublet strengths mu. Triangle methods only; returns false if the method has no such matrix. */
        /** true if some surface is a closed thick body (triangle methods): its sideslip derivatives keep the older model, see computeTranslationDerivatives() */
        virtual bool hasThickSurfaces() const {return false;}
        virtual bool yawedWakeProduct(double theta, std::vector<double> const &mu, std::vector<double> &product) {(void)theta; (void)mu; (void)product; return false;}

        virtual void makeMu(int qrhs) = 0;
        virtual void makeLocalVelocities(std::vector<double> const &uRHS, std::vector<double> const &vRHS, std::vector<double> const &wRHS,
                                         std::vector<Vector3d> &uVLocal, std::vector<Vector3d> &vVLocal, std::vector<Vector3d> &wVLocal,
                                         Vector3d const &WindDirection) const = 0;
        virtual void combineLocalVelocities(double alpha, double beta, std::vector<Vector3d> &VLocal) const = 0;
        virtual void computeOnBodyCp(std::vector<Vector3d> const &VInf, std::vector<Vector3d> const &VLocal, std::vector<double>&Cp) const = 0;
        virtual bool computeTrimmedConditions(double mass, Vector3d const &CoG, double &alphaeq, double &u0, bool bFuseMi) = 0;
        virtual void forces(double const *Mu, double const *Sigma, double alpha, double beta, Vector3d const &CoG, bool bFuseMi, std::vector<Vector3d> const &VInf, Vector3d &Force, Vector3d &Moment) = 0;
        virtual void inducedForce(int nPanel3, double QInf, double alpha, double beta, int pos, Vector3d &ForceBodyAxes, SpanDistribs &SpanResFF) const = 0;
        /** If bTrefftzPlane, the velocity induced at each strip is evaluated analytically in the Trefftz plane, from the trailing legs of all the surfaces,
 *  else it is sampled halfway down the finite wake. */
        virtual void trefftzDrag(int nPanel3, double QInf, double alpha, double beta, int pos, Vector3d &Drag, SpanDistribs &SpanResFF, bool bTrefftzPlane) const = 0;
        virtual int  nPanels() const = 0;
        virtual void makeVertexDoubletDensities(std::vector<double> const &muPanel, std::vector<double> &muNode) const {(void)muPanel; (void)muNode;} //dummy virtual method to enable a call to P3analysis subclass...

        virtual void makeVortons(double dl, double const *mu3Vertex, int pos3, int nPanel3, int nStations, int nVtn0,
                                 std::vector<Vorton> &vortons, std::vector<Vortex> &vortexneg) const = 0;
        virtual void makeRHSVWVelocities(std::vector<Vector3d>& VPW, bool bVLM=false);
        virtual void makeRHSVWVelocitiesBlock(int iFirst, int iLast, bool bVLM, Vector3d *Vpanel);
        virtual void makeNegatingVortices(std::vector<Vortex> &negvortices) = 0;

        virtual void getVelocityVector(Vector3d const &C, double const *Mu, double const *Sigma, Vector3d &VT, double coreradius, bool bWakeOnly, bool bMultiThread=true) const = 0;

        void makeUnitRHSVectors();
        void makeWakeContribution();

        void makeSourceStrengths(Vector3d const &WindDirection);
        void makeSourceStrengths(std::vector<Vector3d> const &WindDirection);

        static void projectTrefftzStrips(Vector3d const &wakedir, std::vector<TrefftzStrip> &strips);
        static Vector3d trefftzVelocity(Vector3d const &pt, Vector3d const &wakedir, std::vector<TrefftzStrip> const &strips);

        void vortonDrag(double alpha, double beta, double QInf, int n0, int nStations, Vector3d &Drag, SpanDistribs &SpanResFF) const;

        void setVortons(const std::vector<std::vector<Vorton> > &vortons);

        void combineUnitRHS(std::vector<double> &RHS, const Vector3d &VInf, const Vector3d &Omega);
        void makeRHS(const std::vector<Vector3d> &VField, std::vector<double> &RHS, const Vector3d *normals);

        void addWakeContribution();
        virtual void computeAddedMass(bool bConsumeMatrix=false);
        /** The added-mass tensor over the density of the closed thick surfaces (computeAddedMass()): 36 values, row-major, the generalised coordinates the
         *  translations along x, y, z and the rotations about x, y, z through addedMassReference(), geometry axes; m^3 between translations, m^4 between a
         *  translation and a rotation, m^5 between rotations. Empty where it is not derived: addedMassReason() says why. */
        std::vector<double> const &addedMass() const {return m_AddedMass;}
        Vector3d const &addedMassReference() const {return m_AddedMassRef;}
        std::string const &addedMassReason() const {return m_AddedMassReason;}
        double addedMassAsymmetry() const {return m_AddedMassAsymmetry;}
        /** The thickness' part of the impulse (P3Analysis::makeThicknessImpulse()), 36 values row-major as addedMass(), not symmetric: row i the impulse's
         *  component, column j the unit motion. Empty where it is not derived: thicknessReason() says why. It is what addThicknessImpulse() puts in the derivatives. */
        std::vector<double> const &thicknessMass() const {return m_ThicknessMass;}
        std::string const &thicknessReason() const {return m_AddedMass.empty() ? m_AddedMassReason : m_ThicknessReason;}
        void addThicknessImpulse(double alphaeq, double u0, Vector3d const &CoG, StabDerivatives &SD) const;
        void thicknessImpulseLoad(Vector3d const &V, Vector3d const &Omega, Vector3d const &CoG, Vector3d &F, Vector3d &G) const;
        bool thicknessMassAbout(Vector3d const &CoG, double Mc[6][6]) const;
        bool pressuresCarryThickness() const;
        void takeAddedMass(PanelAnalysis const *pFrom, std::string const &reason);
        void computeStabilityDerivatives(  double alphaeq, double u0, Vector3d const &CoG, bool bFuseMi, StabDerivatives &SD, Vector3d &Force0, Vector3d &Moment0);
        void computeTranslationDerivatives(double alphaeq, double u0, Vector3d const &CoG, bool bFuseMi, StabDerivatives &SD, Vector3d &Force0, Vector3d &Moment0);
        void computeAngularDerivatives(    double alphaeq, double u0, Vector3d const &CoG, bool bFuseMi, StabDerivatives &SD);

        void traceStdLog(const std::string &str) const;
        inline double sourceStrength(const Vector3d &normal, Vector3d const &Velocity) const {return -1.0/4.0/PI * Velocity.dot(normal);}
        void releasePanelArrays();

        Polar3d const *polar3d() const {return m_pPolar3d;}

        int nVortonRows() const {return int(m_Vorton.size());}
        void clearVortons() {m_Vorton.clear();}
        void getVortonVelocity(Vector3d const &C, double vtncorelength, Vector3d &VelVtn, bool bMultiThread=false) const;
        void getVortonRowVelocity(int iRow, Vector3d const &C, double vtncorelength, Vector3d *VelVtn) const;
        void getVortonVelocityGradient(Vector3d const &C, double *G) const;


        virtual void testResults(double alpha, double beta, double QInf) const = 0;

        void setFarField(bool b) {m_bFarField=b;}
        Vector3d const &farFieldForce() const {return m_ForceFF;}
        Vector3d const &referenceFarFieldForce() const {return m_ForceFF0;}
        Vector3d const &farFieldMoment() const {return m_MomentFF;}
        Vector3d const &referenceFarFieldMoment() const {return m_MomentFF0;}

        static void setMaxThreadCount(int maxthreads) {s_MaxThreads=maxthreads;}
        static void setDoublePrecision(bool bDouble) {s_bDoublePrecision=bDouble;}
        static bool bDoublePrecision() {return s_bDoublePrecision;}

        static void clearDebugPts() {s_DebugPts.clear(); s_DebugVecs.clear();}


    protected:
        bool LUfactorize();
        virtual void backSubUnitRHS(double *uRHS, double *vRHS, double*wRHS, double *pRHS, double *qRHS, double*rRHS);
        bool backSubRHS(std::vector<double> &RHS);
        static bool solveDense(std::vector<double> &A, int n, std::vector<double> &B, int nrhs);

        /** Adds v to the entry (row, col) of the influence matrix, or, in the matrix-vector mode of yawedWakeProduct(), v times mu[col] to entry row of the product. Rows are never shared between threads. */
        void addWakeEntry(int row, int col, int N, double v)
        {
            if(m_pWakeProduct)         (*m_pWakeProduct)[unsigned(row)] += v * m_pWakeMu[col];
            else if(s_bDoublePrecision) m_aijd[unsigned(row*N+col)] += v;
            else                        m_aijf[unsigned(row*N+col)] += float(v);
        }

    protected:

        mutable std::string m_ErrorLog;

        Polar3d const *m_pPolar3d;

        std::vector<double> m_aijd;  /**< the matrix of panel influences - double precision; std::vector is limited to 2 GB and is unusable*/
        std::vector<float>  m_aijf;  /**< the matrix of panel influences - single precision; std::vector is limited to 2 GB and is unusable*/
        std::vector<int>    m_ipiv;  /** the array of pivot indices for the LAPACK LU solver */

        std::vector<double> m_AddedMass;      /**< the added-mass tensor over the density, 6 x 6 row-major, set by computeAddedMass(); empty where it is not derived */
        Vector3d m_AddedMassRef;              /**< the point the tensor's rotations are about: the polar's CoG, as in the unit right-hand sides */
        std::string m_AddedMassReason;        /**< why the tensor is not derived, when it is empty */
        double m_AddedMassAsymmetry = 0.0;    /**< the largest |m_ij - m_ji| over the largest |m_ij| of the collocated tensor before it is symmetrised */
        std::vector<double> m_ThicknessMass;  /**< the thickness' part of the impulse, 6 x 6 row-major, over the density, about m_AddedMassRef; empty where it is not derived */
        std::string m_ThicknessReason;        /**< why m_ThicknessMass is not derived, when the added masses are */


        // unit RHS for the 6 motion d.o.f
        std::vector<double> m_uRHS, m_vRHS, m_wRHS;
        std::vector<double> m_pRHS, m_qRHS, m_rRHS;

        // for control derivatives

        std::vector<double> m_cRHS;

        std::vector<Vector3d> m_uVLocal;              /**< the array of unit velocity vectors for aoa=0,  in local coordinates */
        std::vector<Vector3d> m_vVLocal;              /**< the array of unit velocity vectors for beta=PI/2, in local coordinates */
        std::vector<Vector3d> m_wVLocal;              /**< the array of unit velocity vectors for aoa=PI/2, in local coordinates */

        Vector3d m_ForceFF;                  /**< The sum of the far-field strip forces, at the mid-wake point as the force of the polar, set by forces() with m_MomentFF; forces() returns the force at the trailing wake point */
        Vector3d m_ForceFF0;                 /**< The far-field force of the reference state of computeTranslationDerivatives(), for the control derivatives */
        Vector3d m_MomentFF;                 /**< The moment about the CoG of the far-field strip forces, set by forces(); for VLM2 the moment of the near-field load of the vortex lattice (P4Analysis::boundVortexLoad()) */
        Vector3d m_MomentFF0;                /**< m_MomentFF of the reference state of computeTranslationDerivatives(), for the control derivatives of VLM2 (PlaneTask::computeControlDerivatives()) */
        Vector3d m_Omega;                    /**< The body rotation of the state forces() evaluates, set by computeAngularDerivatives() and zero otherwise: the onset at a point r is V0 + (r - CoG) x m_Omega. VLM2's near-field load takes it at each vortex segment's own midpoint (P4Analysis::boundVortexLoad()), where the onset of VInf is the one at the panel's centre */
        Vector3d m_MomentLocal;              /**< The moment of each strip's own chordwise load (panel pressures and, on thin surfaces, the leading-edge force) about the strip's quarter-chord point, summed over the strips, set by forces() with m_MomentFF: the load that m_MomentFF, which applies each strip's force at its quarter-chord point, does not have */
        Vector3d m_MomentLift;               /**< The part of m_MomentLocal that the rate derivatives p and r take (computeAngularDerivatives(), triangle methods): the moment about the strips' quarter-chord points of their load normal to the mean surface, the panel pressures and, on thin strips, the normal part of the leading-edge force at the leading edge, thick strips cell by cell (P3Analysis::thickStripLiftPressureForce()); the in-plane leading-edge suction stays at the quarter-chord point (P3Analysis::forces(), rule 2). Set by P3Analysis::forces() with m_MomentFF */

        std::vector<double> m_Cp;                /**< The array of pressure coefficients on the panels. 1 value/panel in the case of the quad methods, 3 values/panel in the case of the triangular methods */
        std::vector<double> m_Mu;                /**< The array of doublet strengths, or vortex circulations, associated to the panels. 1 value/panel in the case of the quad methods, 3 values/panel in the case of the triangular methods */
        std::vector<double> m_Sigma;             /**< The array of resulting source strengths of the analysis. 1 value/panel for all methods. */

        bool m_bCancel; /** to interrupt the matrix solver only; */

        bool m_bSequence;           /**< true if the calculation is should be performed for a range of aoa */
        bool m_bWarning;     /**< true if one the OpPoints could not be properly interpolated */
        bool m_bMatrixError;
        bool m_bSideslipFarField = false;   /**< true while computeTranslationDerivatives() evaluates the sideslip (v) derivatives: the far-field force has its Trefftz plane perpendicular to the free stream, along which the wake leaves */
        bool m_bSideslipVorticity;   /**< true while computeTranslationDerivatives() evaluates the sideslip (v) derivatives: forces() then puts the free-stream Kutta-Joukowski force of the doublet sheet's vorticity, chordwise as well as spanwise, at the vortex lines' own positions in m_MomentFF, instead of lumping it at the strips' quarter-chord points (thick strips: their lifting vorticity on their mean surface); VLM2 adds the free-stream force of its rings' chordwise legs to the near-field load of its bound vortices */
        std::vector<double> *m_pWakeProduct;   /**< if not null, makeWakeMatrixBlock() does not touch the matrix: it adds each wake entry times m_pWakeMu to this vector (see yawedWakeProduct()) */
        double const *m_pWakeMu;
        bool m_bRollRate = false;    /**< true while forces() evaluates the roll-rate states of computeAngularDerivatives(): P4Analysis::boundVortexLoad() then loads the chordwise legs with the onset velocity (the side-edge suction of the rolling wing); not in the pitch and yaw rates, where the same load reads against the curved-flow measurements (see there) */
        bool m_bRotationSuction;     /**< true while forces() evaluates a rotating state (computeAngularDerivatives(), P3Analysis::stateSuctionChange()): the strips' leading-edge suction is then directed perpendicular to the leading edge (P3Analysis::stripSuctionCorrection(), flow5#28), the change kept where the strip's force sits (the quarter-chord point in m_MomentFF) until the near-field moment (flow5#29); in a state without rotation the strips' far-field force is exact in total (momentum) and is left as it is. VLM2 returns the near-field force of its bound vortices in a rotating state (P4Analysis::forces()) */
        bool m_bFarField;            /**< true if forces() should also evaluate the far-field force m_ForceFF and moments m_MomentFF, m_MomentLocal and m_MomentLift; it is costly, and only the stability and control derivatives (computeTranslationDerivatives(), computeAngularDerivatives(), PlaneTask::computeControlDerivatives()) need it */

        xfl::enumAnalysisStatus m_AnalysisStatus;
        int m_nBlocks;              /** the number of row blocks for multithreading */

        int m_nStations;          /**< the number of chordwise strips,
                                   which is also the number of panels or stations in the spanwise direction,
                                   which is also the number of wake columns*/


        std::vector<std::vector<Vorton>> m_Vorton; /** The array of vorton rows. Vortons are organized in rows. Each row is located in a crossflow plane. The number of vortons is variable for each row, due to vortex stretching and vorton redistribution. */
        std::vector<Vortex> m_VortexNeg;    /** The array of negating vortices at the trailing edge of the trailing wake panel of each wake column. cf. Willis 2005 fig. 3*/


        static bool s_bDoublePrecision;
        static int s_MaxThreads;


    public:
        static std::vector<Vector3d> s_DebugPts;
        static std::vector<Vector3d> s_DebugVecs;
};





