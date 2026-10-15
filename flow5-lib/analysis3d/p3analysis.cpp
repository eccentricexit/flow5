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



#include <iostream>
#include <format>
#include <thread>
#include <array>


#include <p3analysis.h>

#include <mathelem.h>
#include <matrix.h>
#include <gaussquadrature.h>
#include <polar3d.h>
#include <objects2d.h>
#include <stabderivatives.h>
#include <vortex.h>
#include <vorton.h>



P3Analysis::P3Analysis() : PanelAnalysis()
{
    m_nStations = 0;
    m_pRefTriMesh = nullptr;
}


/** The Analysis has finished, the results have been stored, and the arrays are not needed anymore. */
void P3Analysis::releasePanelArrays()
{
    PanelAnalysis::releasePanelArrays();
    m_Mu.clear();
    m_Cp.clear();
    m_Sigma.clear();
    m_uRHSVertex.clear();
    m_vRHSVertex.clear();
    m_wRHSVertex.clear();
    m_pRHSVertex.clear();
    m_qRHSVertex.clear();
    m_rRHSVertex.clear();
    m_Panel3.clear();
    m_WakePanel3.clear();
}


void P3Analysis::setTriMesh(TriMesh const &trimesh)
{
    m_pRefTriMesh = &trimesh;
    m_refPanel3 = trimesh.panels();

    m_Panel3 = trimesh.panels();
    m_WakePanel3 = trimesh.wakePanels();
}


bool P3Analysis::initializeAnalysis(Polar3d const *pPolar3d, int nRHS)
{
    PanelAnalysis::initializeAnalysis(pPolar3d, nRHS);

    int N = nPanels();

    if(m_pPolar3d->isTriLinearMethod()) N*=3;

    if(nRHS>0)
    {
        if(!allocateMatrix(N)) return false;

        int memsize =     allocateRHS3(nRHS);
        std::string strange = std::format("Allocating {:2f} Mb for {:d} RHS vectors", double(memsize)/1024.0/1024.0, nRHS);
        traceStdLog(strange+"\n\n");
    }

    m_bWarning  = false;

    return true;
}


int P3Analysis::allocateRHS3(int nRHS)
{
    if(!m_pPolar3d) return 0;

    unsigned int memsize = 0;

    int N = m_pPolar3d->isTriUniformMethod() ? nPanels() : nPanels()*3;

    // the new elements are initialized with a default-constructed value,
    // but for primitive types like int and double, as well as for pointer types,
    // the C++ language doesn't specify any initialization;
    // in those cases, Qt's containers automatically initialize the value to 0.

    // nRHS independent
    m_uRHSVertex.resize(nPanels()*3);
    m_vRHSVertex.resize(nPanels()*3);
    m_wRHSVertex.resize(nPanels()*3);

    memsize += 7 * nPanels()*3 * sizeof(double);

    m_uRHS.resize(N);
    m_vRHS.resize(N);
    m_wRHS.resize(N);
    m_pRHS.resize(N);
    m_qRHS.resize(N);
    m_rRHS.resize(N);
    m_cRHS.resize(N);

    memsize += 14 * nPanels() * sizeof(double);

    m_uVLocal.resize(nPanels()*3); //the new elements are initialized with a default-constructed value
    m_vVLocal.resize(nPanels()*3); //the new elements are initialized with a default-constructed value
    m_wVLocal.resize(nPanels()*3);
    memsize += 2 * nPanels()*3 * sizeof(Vector3d);

    // nRHS dependent

    m_Cp.resize(nPanels()*3);
    m_Mu.resize(nPanels()*3);     // also independent of polarType
    memsize +=  2*nPanels()*3* int(sizeof(double));

    int sizeCst = nPanels() * nRHS;
    m_Sigma.resize(sizeCst);
    memsize +=  3*sizeCst* int(sizeof(double));

    return memsize;
}


/**
 * @brief Debug only, get the potential at any given field point
 */
void P3Analysis::getDebugPotential(Vector3d const &C, bool bSelf, double const *Mu, double const *Sigma,
                                   double &phi, bool bSource, bool bDoublet, bool bWake) const
{
    Vector3d *Vs = nullptr;
    Vector3d *Vd = nullptr;

    double phiD[]={0,0,0};
    double phiS=0, sign=1.0;
    phi = 0.0;
    for (int i3=0; i3<nPanels(); i3++)
    {
        Panel3 const & p3 = m_Panel3[i3];
        if(bSource)
        {
            if(!p3.isMidPanel() && !p3.isWakePanel())
            {
                getSourceInfluence(m_pPolar3d, C, p3, bSelf, Vs, &phiS);
                phi += phiS* Sigma[i3];
            }
        }

        if(bDoublet)
        {
            getDoubletInfluence(C, p3, Vd, phiD, 0.0, true);

            if(m_pPolar3d->isTriUniformMethod()) phi += (phiD[0]+phiD[1]+phiD[2]) * Mu[i3]; // weight 1.0 for constant density
            else                                 phi += (phiD[0]*Mu[3*i3] + phiD[1]*Mu[3*i3+1] + phiD[2]*Mu[3*i3+2]);
        }

        if(bWake)
        {
            // Is the panel shedding a wake?
            if(p3.isTrailing())
            {
                //If so, we need to add the contribution of the wake column shedded by this panel
                if(p3.isBotPanel()) sign=-1.0; else sign=1.0;

                assert(m_Panel3.at(i3).iWake()>=0);

                Panel3 const * p3w = &m_WakePanel3[p3.iWake()];
                do
                {
                    // do not use RFF approximation for wake panels
                    getDoubletInfluence(C, *p3w, Vd, phiD, 0.0, false);

                    if(m_pPolar3d->isTriUniformMethod())
                    {
                        phi += (phiD[0]+phiD[1]+phiD[2]) * Mu[i3]*sign; // weight 1.0 for constant density
                    }
                    else if(m_pPolar3d->isTriLinearMethod())
                    {
                        // whether p3 is on the left or right wing, node 1 is its left trailing node and node 2 is its right trailing node
                        // if p3w is a left  wake panel, node 0 and 2 are left  side and node 1 is right side
                        if(p3w->isLeftPanel())  phi += ((phiD[0]+phiD[2])*Mu[3*i3+1] +  phiD[1]         *Mu[3*i3+2]) *sign;
                        // if p3w is a right wake panel, node 0 and 1 are right side and node 2 is left side
                        else                        phi += ( phiD[2]         *Mu[3*i3+1] + (phiD[0]+phiD[1])*Mu[3*i3+2]) *sign;
                    }
                    //is there another wake panel downstream?
                    if(p3w->m_iPD>=0) p3w = m_WakePanel3.data() + p3w->m_iPD;
                    else              p3w = nullptr;
                }
                while(p3w);
            }
        }
    }
}


void P3Analysis::getVelocityVector(Vector3d const &C,
                                   double const *Mu, double const *Sigma,
                                   Vector3d &VT, double coreradius, bool bWakeOnly, bool bMultiThread) const
{
    VelBlockData data;
    data.coreradius = coreradius;
    data.bWakeOnly = bWakeOnly;
    data.Mu = Mu;
    data.Sigma = Sigma;


    std::vector<Vector3d> VBlock(m_nBlocks);

    if(bMultiThread)
    {
        std::vector<std::thread> threads;

        for(int iBlock=0; iBlock<m_nBlocks; iBlock++)
        {
            threads.push_back(std::thread(&P3Analysis::velocityVectorBlock, this, iBlock, C, &VBlock[iBlock], data));
        }

        for(int iBlock=0; iBlock<m_nBlocks; iBlock++)
        {
            threads[iBlock].join();
        }
//        std::cout << "P3Analysis::getVelocityVector joined all " << m_nBlocks << " threads" <<std::endl;
    }
    else
    {
        for(int iBlock=0; iBlock<m_nBlocks; iBlock++)
        {
            velocityVectorBlock(iBlock, C, &VBlock[iBlock], data);
        }
    }

    VT.reset(); // total induced velocity
    for(int ib=0; ib<m_nBlocks; ib++) VT += VBlock[ib];

    if(m_pPolar3d->bVortonWake())
    {
        double vtncorelength = m_pPolar3d->vortonCoreSize()*m_pPolar3d->referenceChordLength();

        Vector3d VVtn;
        getVortonVelocity(C, vtncorelength, VVtn, bMultiThread);
        VT += VVtn;
    }
}


void P3Analysis::velocityVectorBlock(int iBlock, Vector3d const &C, Vector3d *VT, VelBlockData const &data) const
{
    int blockSize = int(nPanels()/m_nBlocks) +1;
    int iStart = iBlock*blockSize;
    int maxRows = nPanels();
    int iMax = std::min(iStart+blockSize, maxRows);

    Vector3d Vd[3], Vs;
    double sign=0;

    for (int i3=iStart; i3<iMax; i3++)
    {
        Panel3 const &p3 = m_Panel3.at(i3);

        if(!data.bWakeOnly)
        {
            if(data.Sigma && fabs(data.Sigma[i3])>0.0)
            {
                getSourceInfluence(m_pPolar3d, C, p3, C.isSame(p3.CoG()), &Vs, nullptr);
                *VT += Vs * data.Sigma[i3];
            }

            getDoubletInfluence(C, p3, Vd, nullptr, data.coreradius, true);
            VT->x += Vd[0].x*data.Mu[3*i3+0] + Vd[1].x*data.Mu[3*i3+1] + Vd[2].x*data.Mu[3*i3+2];
            VT->y += Vd[0].y*data.Mu[3*i3+0] + Vd[1].y*data.Mu[3*i3+1] + Vd[2].y*data.Mu[3*i3+2];
            VT->z += Vd[0].z*data.Mu[3*i3+0] + Vd[1].z*data.Mu[3*i3+1] + Vd[2].z*data.Mu[3*i3+2];

        }

        // Is the panel shedding a wake?
        if(p3.isTrailing())
        {
            //If so, we need to add the contribution of the wake column shedded by this panel
            if(p3.isBotPanel()) sign=-1.0; else sign=1.0;
            if(p3.iWake()<0) continue; // requesting the velocity before the wake has been set - usually the consequence of a repaint signal for streamlines

            //in the case of wake roll-up, only the panels upstream of row nWakeRow
//            int iRow=0;
            Panel3 const *p3w = &m_WakePanel3[p3.iWake()];

            // whether p3 is on the left or right wing, node 1 is its left trailing node and node 2 is its right trailing node
            //Mu3[3*i3+1] is the doublet density at the wing's panel left trailing node, and Mu3[3*i3+3] at the right trailing node
            double mu3left  = data.Mu[3*i3+1];
            double mu3right = data.Mu[3*i3+2];

            while(p3w)
            {
                // do not use RFF approximation for wake panels?
                getDoubletInfluence(C, *p3w, Vd, nullptr, data.coreradius, false);

                if(p3w->isLeftPanel())
                {
                    VT->x += ((Vd[0].x+Vd[2].x)*mu3left +  Vd[1].x       *mu3right) *sign;
                    VT->y += ((Vd[0].y+Vd[2].y)*mu3left +  Vd[1].y       *mu3right) *sign;
                    VT->z += ((Vd[0].z+Vd[2].z)*mu3left +  Vd[1].z       *mu3right) *sign;
                    // if p3w is a left wake panel, node 0 and 2 are left side and node 1 is right side - cf. TriMesh::makeWakePanels()
                }
                else
                {
                    // if p3w is a right wake panel, node 2 is left side and node 0 and 1 are right side
                    VT->x += ( Vd[2].x       *mu3left + (Vd[0].x+Vd[1].x)*mu3right) *sign;
                    VT->y += ( Vd[2].y       *mu3left + (Vd[0].y+Vd[1].y)*mu3right) *sign;
                    VT->z += ( Vd[2].z       *mu3left + (Vd[0].z+Vd[1].z)*mu3right) *sign;
                }
                // is there another wake panel downstream?
                if(p3w->m_iPD>=0) p3w = m_WakePanel3.data() + p3w->m_iPD;
                else              p3w = nullptr;

//                iRow++;
            }
        }
        if(isCancelled()) return;
    }
}


/**
 * Returns the perturbation velocity vector far downstream using a line vortex model for the wake
 * irrespective of the analysis method.
 */
void P3Analysis::getFarFieldVelocity(const Vector3d &C, const std::vector<Panel3> &panel3, const double *Mu,
                                     Vector3d &VT, double coreradius) const
{
    if(isCancelled()) return;

    VT.set(0.0,0.0,0.0);

    Vector3d VL, VR, A, B;
    double sign = 1.0;
    double fardist = m_pPolar3d->TrefftzDistance();
    for (int i3=0; i3<int(panel3.size()); i3++)
    {
        if(isCancelled()) return;
        Panel3 const &p3 = panel3.at(i3);
        if(p3.isTrailing())
        {
            sign = p3.isBotPanel() ? -1.0 : 1.0;
            A.set(p3.leftTrailingNode());
            A.x += fardist;
            VL = vortexInducedVelocity(p3.leftTrailingNode(), A, C, coreradius);
            VT +=  VL*(+Mu[3*i3+1] * sign);

            B.set(p3.rightTrailingNode());
            B.x += fardist;
            VR = vortexInducedVelocity(p3.rightTrailingNode(), B, C, coreradius);
            // circulations of each seg are oppposite
            VT +=  VR*(-Mu[3*i3+2] * sign);
        }
    }
}


/**
* Returns the perturbation potential at a given point, induced by the distribution of source and doublet/circulation strengths.
* @param C the point where the influence is to be evaluated
* @param Mu a pointer to the array of doublet strength or vortex circulations
* @param sigma a pointer to the array of source strengths
* @param VT the resulting perturbation velocity
* @param nWakeRow in the case of wake roll-up, only the panels upstream of row nWakeRow will be taken into account
*/
double P3Analysis::getPotential(Vector3d const &C, double const *mu, double const*sigma) const
{
    double phiT=0, phiSource=0;
    double phiBasis[]{0,0,0};
    double sign=0;
    for (int i3=0; i3<nPanels(); i3++)
    {
        Panel3 const &p3 = m_Panel3.at(i3);
        getSourceInfluence(m_pPolar3d, C, p3, false, nullptr, &phiSource);
        phiT += phiSource * sigma[i3];

        getDoubletInfluence(C, p3, nullptr, phiBasis, 0.0, true);
//        strange += std::format("   {:9g}", phiBasis[0]+phiBasis[1]+phiBasis[2]);
        phiT += phiBasis[0]*mu[3*i3] + phiBasis[1]*mu[3*i3+1] + phiBasis[2]*mu[3*i3+2];

        // Is the panel shedding a wake?
        if(p3.isTrailing())
        {
            //If so, we need to add the contribution of the wake column shedded by this panel
            if(p3.isBotPanel()) sign=-1.0; else sign=1.0;
//            assert(m_Panel3.at(i3).m_iWake>=0);

            //in the case of wake roll-up, only the panels upstream of row nWakeRow
//            int iRow=0;
            Panel3 const *p3w = &m_WakePanel3[p3.iWake()];
            while(p3w )
            {
                // do not use RFF approximation for wake panels
                getDoubletInfluence(C, *p3w, nullptr, phiBasis, 0.0, false);

                // whether p3 is on the left or right wing, node 1 is its left trailing node and node 2 is its right trailing node
                // if p3w is a left  wake panel, node 0 and 2 are left  side and node 1 is right side
                if(p3w->isLeftPanel())
                    phiT += ((phiBasis[0]+phiBasis[2])*mu[3*i3+1] + phiBasis[1] *mu[3*i3+2]) *sign;
                // if p3w is a right wake panel, node 0 and 1 are right side and node 2 is right side
                else
                    phiT += ( phiBasis[2]*mu[3*i3+1] + (phiBasis[0]+phiBasis[1])*mu[3*i3+2]) *sign;

                //is there another wake panel downstream?
                if(p3w->m_iPD>=0) p3w = m_WakePanel3.data() + p3w->m_iPD;
                else              p3w = nullptr;

 //               iRow++;
            }
        }
    }
//qDebug("{:s}", strange.toStdString().c_str());
    return phiT;
}


/**
* Returns the influence at point C of a uniform source distribution on the panel pPanel
* The panel is necessarily located on a thick surface, else the source strength is zero
* @param C the point where the influence is to be evaluated
* @param pPanel a pointer to the Panel with the doublet strength
* @param V the perturbation velocity at point C
* @param phi the potential at point C
*/
void P3Analysis::getSourceInfluence(Polar3d const *pPolar3d, Vector3d const &C, Panel3 const &p3, bool bSelf, Vector3d *V, double *phi) const
{
    if(phi) p3.sourcePotential(C, bSelf, *phi);
    if(V)
        p3.sourceN4023Velocity(C, bSelf, *V, Vortex::coreRadius());

    if(pPolar3d->bGroundEffect() || pPolar3d->bFreeSurfaceEffect())
    {
        double coef = pPolar3d->bGroundEffect() ? 1.0 : -1.0;

        // To enforce the condition Vz=0 on the ground, add the contribution
        // of a symmetric panel oppposiwte the ground surface, so that the z
        // contributions cancel each other
        // To avoid building a symmetric panel, calculate the influence of this panel
        // at a symmetric point below the surface, and add its z-opposite contribution

        Vector3d VG, CG;
        double phiG;
        CG.set(C.x, C.y, -C.z-2.0*m_pPolar3d->groundHeight());
        if(phi)
        {
            p3.sourcePotential(CG, false, phiG);
            *phi += phiG * coef;
        }
        if(V)
        {
            p3.sourceVelocity(CG, false, VG);
            V->x += VG.x * coef;
            V->y += VG.y * coef;
            V->z -= VG.z * coef;
        }
    }
}


/**
 * @brief For a uniform strength analysis, makes the doublet density at the nodes from the values at the panel's center of any two arrays
 */
void P3Analysis::makeNodeAverage(int iNode, std::vector<double> const &muPanel, std::vector<double> &muLin) const
{
    Node const &nd = m_pRefTriMesh->nodeAt(iNode);
    if(!nd.triangleCount()) return;

    //average the density over the panels that have been listed
    double muNode = 0.0;
    for(int i3=0; i3<nd.triangleCount(); i3++)
    {
        muNode += muPanel.at(nd.triangleIndex(i3));
    }
    assert(nd.triangleCount()>0);
    muNode /= nd.triangleCount();

    //store the result
    for(int i3=0; i3<nd.triangleCount(); i3++)
    {
        int idx3 = nd.triangleIndex(i3);
        Panel3 const &p3 = m_Panel3.at(idx3);
        for(int k=0; k<3; k++)
        {
            if(p3.nodeIndex(k)==iNode)
            {
                muLin[3*idx3+k]    = muNode;
            }
        }
    }
}


/**
 * This method calculates the Cp coefficient on a panel based on the distribution of doublet strengths.
 * Finds the plane best fitting the doublet values at the three vertices
 * and calculates it slope.
 * The local velocity is the slope of the plane.
 * @param p the index of the panel for which the calculation is performed
 * @param Mu the array of doublet strength values
 * @param VLocal a reference to the calculated local velocity on the panel
 */
void P3Analysis::makePanelDoubletSurfaceVelocity(int p, double const *Mu, Vector3d &VLocal)
{
    Vector3d const* sl = m_Panel3.at(p).m_Sl; // use local reference frame

    assert(p*3+2<3*nPanels());

    //find the gradient of steepest descent from the three values at the vertices
    double mu0 = Mu[p*3];
    double mu1 = Mu[p*3+1];
    double mu2 = Mu[p*3+2];

    // the three points which define the plane are A(s0.x, s0.y, µ0), B(s1.x, s1.y, µ1), C(s2.x, s2.y, µ2)
    // the plane equation is ax+by+cµ+d=0, or as a function: µ(x,y) = -(ax+by+d)/c
    double a = (sl[1].y-sl[0].y) * (mu2-mu0)         - (sl[2].y-sl[0].y)*(mu1-mu0);
    double b = (sl[2].x-sl[0].x) * (mu1-mu0)         - (sl[1].x-sl[0].x)*(mu2-mu0);
    double c = (sl[1].x-sl[0].x) * (sl[2].y-sl[0].y) - (sl[2].x-sl[0].x)*(sl[1].y-sl[0].y);
//  double d = -(a*sl[0].x + b*sl[0].y + c*mu0);  // Unused

    // the gradient is (dµ/dx, dµ/dy) =  (-a/c, -b/c)
    VLocal.x = -4.0*PI * (-a/c);
    VLocal.y = -4.0*PI * (-b/c);
    VLocal.z = 0.0;
//    VLocal.z = VinfL.z +4.0*PI * m_Sigma[p]; // is zero
}


double P3Analysis::stripArea(Panel3 const &p3, bool bThinSurfaces) const
{
    double striparea = 0.0;
    //sum panel areas of bottom strip
    int index = p3.index();
    do
    {
        Panel3 const &p3k = m_Panel3.at(index);
        striparea += p3k.area();
        if(p3k.m_iPU>=0) index = p3k.m_iPU;
        else             break;
    }
    while (index>=0); //  = while(true)
    if(!bThinSurfaces) striparea /= 2;

    return striparea;
}


/**
 * Returns an estimate of the strip's quarter-chord point, used as the application point
 * of the strip's far-field force when summing moments.
 * The leading edge point is the midpoint of the strip's leading edge, found by walking the strip
 * upstream from the trailing panel via the iPU panel links;
 * for thick surfaces the walk stops at the leading edge i.e. where the chain reaches the top surface.
 * The last triangle of the walk is the upper triangle of the leading cell, see Panel3::leadingEdgeMidpoint().
 */
Vector3d P3Analysis::stripQuarterPoint(Panel3 const &p3) const
{
    // the trailing nodes have indexes 1 & 2
    Vector3d ptTE = (p3.vertexAt(1) + p3.vertexAt(2)) * 0.5;

    Panel3 const *pLast = &p3;
    while(pLast->m_iPU>=0)
    {
        Panel3 const &p3k = m_Panel3.at(pLast->m_iPU);
        if(p3k.isTopPanel()) break;
        pLast = &p3k;
    }

    Vector3d ptLE = pLast->leadingEdgeMidpoint();

    return ptLE + (ptTE-ptLE)*0.25;
}


/**
 * The frame of a strip's leading-edge suction (flow5#28), from its trailing panel (the mid panel of a thin strip, the bottom one of a thick strip):
 * n the normal of the mean surface at the leading edge (thin: the leading panel's; thick: the difference of the normals of the top and bottom panels
 * that meet at the leading edge), e = n x x the span direction in that plane, t = e x n the chordwise direction, d = n x l_LE the direction in that
 * plane perpendicular to the strip's leading-edge segment l_LE, and ptLE the midpoint of that segment. False for a degenerate strip.
 */
bool P3Analysis::stripSuctionFrame(Panel3 const &p3, Vector3d &n, Vector3d &e, Vector3d &t, Vector3d &d, Vector3d &ptLE) const
{
    Panel3 const *pLast = &p3;
    while(pLast->m_iPU>=0)
    {
        Panel3 const &p3k = m_Panel3.at(pLast->m_iPU);
        if(p3k.isTopPanel()) break;
        pLast = &p3k;
    }

    // the leading nodes of the last panel, see Panel3::leadingEdgeMidpoint()
    int iA=0, iB=1;
    if(pLast->isMidPanel())
    {
        if(pLast->isLeftWingPanel()) {iA=0; iB=2;}
        else                         {iA=1; iB=2;}
        n = pLast->normal();
    }
    else
    {
        if(pLast->m_iPU<0) return false;
        n = m_Panel3.at(pLast->m_iPU).normal() - pLast->normal();
    }
    Vector3d lLE = pLast->vertexAt(iB) - pLast->vertexAt(iA);
    ptLE = (pLast->vertexAt(iA) + pLast->vertexAt(iB)) * 0.5;

    if(n.norm()<1.0e-9) return false;
    n.normalize();
    e = n * Vector3d(1.0, 0.0, 0.0);
    if(e.norm()<1.0e-6) return false;     // a surface facing the x axis has no leading edge to speak of
    e.normalize();
    t = e * n;
    d = n * lLE;
    if(d.norm()<1.0e-9) return false;
    d.normalize();
    return fabs(d.dot(t))>1.0e-6;
}


/**
 * The change of a strip's force F (N) that directs its leading-edge suction perpendicular to the leading edge (flow5#28).
 *
 * In thin-wing theory the force on a lifting surface is the pressure, normal to the surface, plus the leading-edge suction, which lies in the plane of
 * the surface at the leading edge and is perpendicular to the leading edge (Polhamus, NASA TN D-3767, 1966; Lan, J. Aircraft 11(9), 1974; Katz and
 * Plotkin, Low-Speed Aerodynamics, 2nd ed., 2001, sec. 12.4). The strip's Kutta-Joukowski force rho Gamma V x l_TE has the right lift and chordwise
 * force but is perpendicular to the trailing edge: where the trailing edge is swept differently from the leading edge (a tapered or elliptic
 * planform, a twisted wing) its component along the span is wrong, and so are the side force of a rolling or yawing wing and its yawing moment.
 *
 * The strip's force is kept along n (the normal force) and t (the chordwise force); along e it becomes that of the panel pressures FP (N), which
 * follow the real surface, its twist and dihedral included, plus that of the leading-edge force F - FP, whose component in the plane (the suction)
 * is made perpendicular to the leading edge with its chordwise component unchanged:
 *     F.e  ->  FP.e + ((F-FP).t) (d.e)/(d.t)
 * Nothing is chosen along the chord. Thick strips take the same with FP the pressure force of their lift (thickStripLiftPressureForce(): the panel forces
 * projected on the mean surface's normal, cell by cell), so that their suction, which their pressures carry in the plane of the mean surface, is directed
 * as a thin strip's. With the whole pressure force instead (flow5#30) a thick strip's span component took the thickness' own pressures as well: the
 * fluid's impulse turned by a yaw rate, -Omega x (M_a V) (Kirchhoff; Lamb, Hydrodynamics, sec. 124), which the thin sheet and AVL have not, and only the
 * part of it on the strips' spanwise slopes, not on the tips (flow5#31: CYr of an elliptic NACA 0012 wing -0.0057 against -0.0009 thin, of the
 * rectangular one the thin's). Returns the change, along e.
 * forces() applies it in rotating states only (m_bRotationSuction): without rotation the strips' far-field force is exact in total.
 */
Vector3d P3Analysis::stripSuctionCorrection(Panel3 const &p3, Vector3d const &F, Vector3d const &FP) const
{
    Vector3d n, e, t, d, ptLE;
    if(!stripSuctionFrame(p3, n, e, t, d, ptLE)) return Vector3d();
    double r = d.dot(e)/d.dot(t);
    double Fe = FP.dot(e) + (F-FP).dot(t)*r;
    return e * (Fe - F.dot(e));
}


/**
 * The change that the strips' leading-edge suction (stripSuctionCorrection()) makes to the far-field force (N) and to the moment (N.m, about the CoG)
 * of a steady rotating state of a polar point (PlaneTask::solveRatePoint()): the change in the state with the body rates, solution sol and onset
 * field, less the change in the state of the same free stream V without rates. A state without rotation keeps its far-field force (see forces(),
 * rule 1), so only the part due to the rotation is added. The moment's change is placed as the rate derivatives place it (forces(), rule 2): its
 * rolling and yawing parts are those of m_MomentFF, the change at the strip's quarter-chord point, its pitching part that of m_MomentFF +
 * m_MomentLocal, as Mq; so a point with small rates differs from the point without by the rate derivatives times the rates. sol and sol0 are the
 * plain solutions, without the yawed wake's change that PlaneTask::applySideslipModel() adds: the two changes are superposed (forces(), rule 3),
 * exact to first order in beta and the rates.
 */
void P3Analysis::stateSuctionChange(std::vector<double> const &sol, std::vector<Vector3d> const &field, Vector3d const &V, double alpha, Vector3d const &CoG,
                                    bool bFuseMi, Vector3d &dForce, Vector3d &dMoment)
{
    dForce.set(0.0, 0.0, 0.0);
    dMoment.set(0.0, 0.0, 0.0);
    std::vector<double> sol0;
    makeStateSolution(V, Vector3d(), CoG, sol0);
    std::vector<Vector3d> field0(nPanels(), V);

    bool bFarField = m_bFarField, bRotation = m_bRotationSuction;
    m_bFarField = true;
    for(int k=0; k<2; k++)
    {
        std::vector<double> const &x = k==0 ? sol : sol0;
        std::vector<Vector3d> const &f = k==0 ? field : field0;
        std::vector<double> none(x.size(), 0.0);
        std::vector<double> mu3(3*nPanels());
        if(m_pPolar3d->isTriUniformMethod()) makeVertexDoubletDensities(x, mu3);
        else                                 std::copy(x.begin(), x.end(), mu3.begin());
        std::vector<double> sigma(nPanels(), 0.0);   // thin surfaces: no sources (rotating points are refused on thick surfaces)
        makeLocalVelocities(x, none, none, m_uVLocal, m_vVLocal, m_wVLocal, objects::windDirection(alpha, 0.0));
        computeOnBodyCp(f, m_uVLocal, m_Cp);

        Vector3d F, M;
        m_bRotationSuction = false;
        forces(mu3.data(), sigma.data(), alpha, 0.0, CoG, bFuseMi, f, F, M);
        Vector3d FF0 = m_ForceFF, MFF0 = m_MomentFF, ML0 = m_MomentLocal;
        m_bRotationSuction = true;
        forces(mu3.data(), sigma.data(), alpha, 0.0, CoG, bFuseMi, f, F, M);
        double sign = k==0 ? 1.0 : -1.0;
        Vector3d dMFF = m_MomentFF - MFF0, dML = m_MomentLocal - ML0;
        dForce  += (m_ForceFF - FF0) * sign;
        dMoment += Vector3d(dMFF.x, dMFF.y + dML.y, dMFF.z) * sign;
    }
    m_bFarField = bFarField;
    m_bRotationSuction = bRotation;
}


int P3Analysis::nextTopTrailingPanelIndex(Panel3 const &p3) const
{
    if(!p3.isBotPanel()) return -1;

/*    int index = p3.index();
    do
    {
        Panel3 const &p3k = m_Panel3.at(index);
        if(p3k.m_iPU>=0) index = p3k.m_iPU;
        else             return index;
    }
    while (index>=0); //  = while(true)

    return -1;*/

    return p3.oppositeIndex();
}


void P3Analysis::trailingWakePanels(const Panel3 *pWakePanel, Panel3 &p3WU, Panel3 &p3WD) const
{
    do
    {
        if(pWakePanel->m_iPD<0)
        {
            break;
        }
        else pWakePanel = m_WakePanel3.data() + pWakePanel->m_iPD;
    }
    while (pWakePanel);

    if(pWakePanel)
    {
        p3WD = m_WakePanel3[pWakePanel->index()];
        p3WU = m_WakePanel3[pWakePanel->m_iPU];
    }
}

/**
 * Returns the center point of the of the wake column to which the panel belongs
 * Assumes that pWakePanel points to the leading panel of the wake column
 */
void P3Analysis::midWakePoint(Panel3 const*pWakePanel, Vector3d &midleft, Vector3d &midright) const
{
    Vector3d LA, LB, TA, TB;
    if(!pWakePanel) return;
    LA = pWakePanel->vertexAt(2);
    LB = pWakePanel->vertexAt(1);

    do
    {
        if(pWakePanel->m_iPD<0)
        {
            TA.set(pWakePanel->leftTrailingNode());
            TB.set(pWakePanel->rightTrailingNode());
            break;
        }
        else pWakePanel = m_WakePanel3.data() + pWakePanel->m_iPD;
    }
    while (pWakePanel);

    midleft  = (LA+TA)*0.5;
    midright = (LB+TB)*0.5;
}


void P3Analysis::trailingWakePoint(Panel3 const*pWakePanel, Vector3d &left, Vector3d &right) const
{
    do
    {
        if(pWakePanel->m_iPD<0)
        {
            left = pWakePanel->leftTrailingNode();
            right = pWakePanel->rightTrailingNode();
            break;
        }
        else pWakePanel = m_WakePanel3.data() + pWakePanel->m_iPD;
    }
    while (pWakePanel);
}


int P3Analysis::makeWakePanels(const Vector3d &WindDirection, bool bVortonWake)
{
    if(!bVortonWake)
    {
        Vector3d pt;
        m_pRefTriMesh->getLastTrailingPoint(pt);
        TriMesh::makeWakePanels(m_Panel3,
                                m_pPolar3d->NXWakePanel4(), m_pPolar3d->wakePanelFactor(), pt.x + m_pPolar3d->wakeLength(),
                                WindDirection, m_WakePanel3, m_nStations, true);
    }
    else
    {
        TriMesh::makeWakePanels(m_Panel3,
                                m_pPolar3d->NXBufferWakePanels(), 1.0, m_pPolar3d->bufferWakeLength(),
                                WindDirection, m_WakePanel3, m_nStations, false);
    }

    return nWakePanels();
}


/**
 * Whether the panel i3 is left out of the product of yawedWakeProduct() for the wake column shed by the trailing panel k3: a panel of another wing within
 * half the reference chord of the sheet of that column and downstream of it. The sheet's lateral velocity jumps across it, so that such a point
 * (a fin crossing the wake of the wing) has no derivative with respect to the yaw of the wake that a mesh can resolve: it keeps the unyawed wake.
 */
bool P3Analysis::wakeRowExcluded(int i3, int k3) const
{
    if(m_Component.empty() || m_Component[i3]==m_Component[k3]) return false;
    Panel3 const &p3k = m_Panel3.at(k3);
    Vector3d Tl = p3k.leftTrailingNode(), Tr = p3k.rightTrailingNode();
    Vector3d e = Tr - Tl;
    double L = e.norm();
    if(L<=0.0) return false;
    e *= 1.0/L;
    Vector3d nw = e * Vector3d(1.0, 0.0, 0.0);
    double nn = nw.norm();
    if(nn<1e-9) return false;
    nw *= 1.0/nn;
    Vector3d d = m_Panel3.at(i3).CoG() - Tl;
    double c = m_WakeClearance;
    return d.x>0.0 && fabs(d.dot(nw))<c && d.dot(e)>-c && d.dot(e)<L+c;
}


bool P3Analysis::hasThickSurfaces() const
{
    for(Panel3 const &p3 : m_Panel3)
        if(p3.isTrailing() && p3.isBotPanel()) return true;
    return false;
}


/**
 * product = (dA_wake/dtheta).mu, by a central difference of the wake's rows of the influence matrix between the wake yawed by +theta and by -theta
 * (it leaves each trailing edge along (cos theta, sin theta, 0)), applied to the doublet strengths mu without building the matrix.
 * The wake panels of the analysis are restored on return.
 */
bool P3Analysis::yawedWakeProduct(double theta, std::vector<double> const &mu, std::vector<double> &product)
{
    int N = matSize();
    if(int(mu.size())!=N || theta<=0.0) return false;

    // the wings of the mesh: panels sharing a node are one piece
    {
        std::vector<int> parent(m_pRefTriMesh->nodeCount()+1);
        for(unsigned i=0; i<parent.size(); i++) parent[i] = int(i);
        auto find = [&](int a){while(parent[a]!=a){parent[a]=parent[parent[a]]; a=parent[a];} return a;};
        for(Panel3 const &p3 : m_Panel3)
            for(int k=1; k<3; k++)
            {
                int a = find(p3.nodeIndex(0)), b = find(p3.nodeIndex(k));
                if(a!=b) parent[a] = b;
            }
        m_Component.resize(nPanels());
        for(int i=0; i<nPanels(); i++) m_Component[i] = find(m_Panel3.at(i).nodeIndex(0));
    }
    m_WakeClearance = 0.5 * m_pPolar3d->referenceChordLength();

    m_RowSkip.assign(nPanels(), 0);
    {
        std::vector<int> te;
        for(int k=0; k<nPanels(); k++) if(m_Panel3.at(k).isTrailing() && (m_Panel3.at(k).isBotPanel() || m_Panel3.at(k).isMidPanel())) te.push_back(k);
        for(int i=0; i<nPanels(); i++)
            for(int k : te)
                if(wakeRowExcluded(i, k)) {m_RowSkip[i] = 1; break;}
        int nSkip = 0;
        for(char c : m_RowSkip) nSkip += c;
        if(nSkip>0)
            traceStdLog(std::format("      sideslip derivatives: {:d} panels lie within half a reference chord of another surface's wake sheet and keep the unyawed wake\n", nSkip));
    }

    std::vector<Panel3> wake0 = m_WakePanel3;
    std::vector<double> pp(N, 0.0), pm(N, 0.0);
    bool bOk = true;
    for(int s=0; s<2 && bOk; s++)
    {
        double th = s==0 ? theta : -theta;
        std::vector<double> &out = s==0 ? pp : pm;
        makeWakePanels(Vector3d(cos(th), sin(th), 0.0), false);
        m_pWakeProduct = &out;
        m_pWakeMu = mu.data();
        makeWakeContribution();
        m_pWakeProduct = nullptr;
        m_pWakeMu = nullptr;
        bOk = !m_bMatrixError;
    }
    m_WakePanel3 = wake0;
    if(!bOk) return false;

    product.resize(N);
    for(int i=0; i<N; i++) product[i] = (pp[i]-pm[i])/(2.0*theta);
    return true;
}


int P3Analysis::matSize() const
{
    if(!m_pPolar3d) return nPanels();
    if(m_pPolar3d->isTriUniformMethod()) return nPanels();
    if(m_pPolar3d->isTriLinearMethod())  return nPanels()*3;
    return nPanels();
}


void P3Analysis::savePanels()
{
    m_refPanel3     = m_Panel3;
    m_refWakePanel3 = m_WakePanel3;
}


void P3Analysis::restorePanels()
{
    m_Panel3     = m_refPanel3;
    m_WakePanel3 = m_refWakePanel3;
}


void P3Analysis::scaleResultsToSpeed(double ratio)
{
    //______________________________________________________________________________________
    // Scale RHS and Sigma i.a.w. speeds (so far we have unit doublet and source strengths)
    for(int pp=0; pp<3*nPanels(); pp++)
    {
        m_Mu[pp] *= ratio;
    }
    for(int pp=0; pp<nPanels(); pp++)
    {
        m_Sigma[pp] *= ratio;
    }
}


/**
 * Evaluates the cross-flow forces by application of the Kutta-Joukowski theorem
 * Drela § 5.7.
 * The induced drag is evaluated separately in the Trefftz plane or in the vorton wake
 * */
void P3Analysis::inducedForce(int nPanel3, double QInf, double alpha, double beta,
                              int pos3, Vector3d &ForceBodyAxes, SpanDistribs &SpanResFF) const
{
    Vector3d VInf;
    double gLeft=0, gRight=0, gMid=0;
    //dynamic pressure, kg/m³

    //   Define wind axes
    Vector3d winddir = objects::windDirection(alpha, beta);
    VInf = winddir * QInf;

    double q = 0.5 * m_pPolar3d->density() * QInf * QInf;
    Vector3d StripForce; // body axes

    double const *mu3 = m_Mu.data();
    int m=0;
    for(int i3=0; i3<nPanel3; i3++)
    {
        int index = i3+pos3;
        Panel3 const &p3 = m_Panel3.at(index);
        Vector3d const &surfacenormal = p3.surfaceNormal();

        if(p3.isTrailing() && (p3.isBotPanel() || p3.isMidPanel()))
        {
            if(p3.isMidPanel())
            {
                int idxM = p3.index();
                gLeft  = -mu3[3*idxM+1] *4.0*PI;
                gRight = -mu3[3*idxM+2] *4.0*PI;
            }
            else if(p3.isBotPanel())
            {
                int idxB = p3.index();
                int idxU = nextTopTrailingPanelIndex(p3);
                assert(idxU>=0);

                // the trailing nodes have indexes 1 & 2
                gLeft  = (-mu3[3*idxU+1] + mu3[3*idxB+1])*4.0*PI;
                gRight = (-mu3[3*idxU+2] + mu3[3*idxB+2])*4.0*PI;
            }

            gMid = (gLeft+gRight)/2.0;

            SpanResFF.m_Gamma[m] = gMid;

            Vector3d vtx = p3.trailingVortex();
            // a steady rotating point: the strip's force is the one of the local onset flow at its trailing panel, as in forces() (m_Onset, m/s)
            Vector3d const &VStrip = m_Onset.empty() ? VInf : m_Onset.at(index);
            StripForce = VStrip * vtx;
            StripForce *= gMid * m_pPolar3d->density();     // N
            StripForce *=  1.0/q;      // N/q

            SpanResFF.m_Cl[m]  = StripForce.dot(surfacenormal) /SpanResFF.stripArea(m);
            SpanResFF.m_F[m]   = StripForce * q;                        // N, body axes
            if(!m_Onset.empty())
            {
                // the strip's Cl is the one of its own dynamic pressure, which is what its section sees (the station Cl of the stall check, XFoil's Cl and Re)
                if(int(SpanResFF.m_VOnset.size())!=int(SpanResFF.m_F.size())) SpanResFF.m_VOnset.assign(SpanResFF.m_F.size(), Vector3d());
                SpanResFF.m_VOnset[m] = VStrip * (1.0/QInf);
                SpanResFF.m_Cl[m] /= SpanResFF.m_VOnset[m].dot(SpanResFF.m_VOnset[m]);
            }

            m++;
        }
    }

    for(unsigned int m=0; m<SpanResFF.m_F.size(); m++)
        ForceBodyAxes += SpanResFF.m_F.at(m);
    ForceBodyAxes *= 1.0/q;    // N/q, body axes
}


/**
 * The trailing strips of every surface of the plane as the point vortices of the Trefftz plane perpendicular to wakedir:
 * the strip's circulation is the one of the polar's induced force (mid panel, or bottom and top trailing panels on thick surfaces).
 * Used by trefftzDrag() and by forces(), so that the polar and the derivatives share one far-field velocity.
 */
void P3Analysis::makeTrefftzStrips(double const *mu3, Vector3d const &wakedir, std::vector<TrefftzStrip> &strips) const
{
    strips.clear();
    for(int index=0; index<nPanels(); index++)
    {
        Panel3 const &p3 = m_Panel3.at(index);
        if(!p3.isTrailing() || !(p3.isBotPanel() || p3.isMidPanel())) continue;

        double gamma=0;
        if(p3.isMidPanel())
        {
            int idxM = p3.index();
            gamma = -(mu3[3*idxM+1]+mu3[3*idxM+2])/2.0 *4.0*PI;
        }
        else
        {
            int idxB = p3.index();
            int idxU = nextTopTrailingPanelIndex(p3);
            assert(idxU>=0);
            gamma = (-(mu3[3*idxU+1]+mu3[3*idxU+2]) + (mu3[3*idxB+1]+mu3[3*idxB+2]))/2.0 *4.0*PI;
        }
        TrefftzStrip s;
        s.iPanel = index;
        s.A = p3.leftTrailingNode();
        s.B = p3.rightTrailingNode();
        s.Gamma = gamma;
        strips.push_back(s);
    }
    projectTrefftzStrips(wakedir, strips);
}


void P3Analysis::trefftzDrag(int nPanel3, double QInf, double alpha, double beta, int pos3,
                             Vector3d &Drag, SpanDistribs &SpanResFF, bool bTrefftzPlane) const
{
    Vector3d left, right;
    Vector3d Wg_m;
    Vector3d midwakept;
    Vector3d ForceBodyAxes; // Strip and global forces, body axes
    double gLeft(0), gRight(0), gMid(0);

    //   Define wind axis
    Vector3d winddir = objects::windDirection(alpha, beta);
    // The wake sheet leaves the trailing edges along the x axis of the mesh (the free stream in a control polar, which turns the mesh):
    // the Trefftz plane is perpendicular to it, so that the induced velocity is normal to the sheet and the induced force along x
    Vector3d wakedir(1.0, 0.0, 0.0);

    //dynamic pressure, kg/m³
    double qDyn = 0.5 * m_pPolar3d->density() * QInf * QInf;

    Vector3d StripForce; // body axes
    Vector3d theforce;

    //dynamic pressure, kg/m³

    double const *mu3    = m_Mu.data();
    double const *sigma3 = m_Sigma.data();

//    clearDebugPts();

    if(SpanResFF.m_FInduced.size()!=SpanResFF.m_F.size()) SpanResFF.m_FInduced.assign(SpanResFF.m_F.size(), Vector3d());

    // The Trefftz plane is perpendicular to the wake sheet (the mesh x axis). The trailing legs of every strip of every surface, projected in it, are point vortices:
    // the finite length of the wake and the direction of the wake sheet do not enter.
    std::vector<TrefftzStrip> strips;
    if(bTrefftzPlane) makeTrefftzStrips(mu3, wakedir, strips);

    // Note: parallelization fails, incompatibility with std::vectors of SpanDistribs
    int m=0;
    for(int i3=0; i3<nPanel3; i3++)
    {
        int index = i3+pos3;
        Panel3 const &p3 = m_Panel3.at(index);
        Vector3d const &surfacenormal = p3.surfaceNormal();
        StripForce.set(0,0,0);
        if(p3.isTrailing() && (p3.isBotPanel() || p3.isMidPanel()))
        {
            if(bTrefftzPlane)
            {
                // the velocity at the strip's mid-point, in the plane; half of it is the velocity at the bound vortex
                Vector3d pt = (p3.leftTrailingNode() + p3.rightTrailingNode())/2.0;
                pt -= wakedir * pt.dot(wakedir);
                Wg_m = trefftzVelocity(pt, wakedir, strips) * 0.5;
            }
            else
            {
            // get the last triangle of the wake column
            assert(p3.iWake()>=0 && p3.iWake()<nWakePanels());
            Panel3 const *p3W = m_WakePanel3.data() + p3.iWake();

            // evaluate the induced drag at the half of the wake's length to avoid end effects
            // exact position is not significant, result is not affected by panel side singularities

            midWakePoint(p3W, left, right);
            midwakept.set((left + right)/2.0);

            getVelocityVector(midwakept, mu3, sigma3, Wg_m, 0.0001, true);

            Wg_m *= 0.5;
            }

//s_DebugPts.push_back(mid);
//s_DebugVecs.push_back(Wg_m);

            if(p3.isMidPanel())
            {
                int idxM = p3.index();
                gLeft  = -mu3[3*idxM+1] *4.0*PI;
                gRight = -mu3[3*idxM+2] *4.0*PI;
                gMid   = (gLeft+gRight)/2.0;
            }
            else if(p3.isBotPanel())
            {
                int idxB = p3.index();
                int idxU = nextTopTrailingPanelIndex(p3);
                assert(idxU>=0);

                // the trailing nodes have indexes 1 & 2
                gLeft  = (-mu3[3*idxU+1] + mu3[3*idxB+1])*4.0*PI;
                gRight = (-mu3[3*idxU+2] + mu3[3*idxB+2])*4.0*PI;
                gMid = (gLeft+gRight)/2.0;
            }

            theforce = Wg_m * p3.trailingVortex() * gMid;

            StripForce += theforce;

            StripForce *= m_pPolar3d->density() / qDyn;      // N/q
            ForceBodyAxes += StripForce;                     // N/q
            SpanResFF.m_FInduced[m] = StripForce * qDyn;     // N, body axes

            SpanResFF.m_Vd[m]  = Wg_m;
            SpanResFF.m_Ai[m]  = atan2(Wg_m.dot(surfacenormal), QInf)*180.0/PI;
            SpanResFF.m_ICd[m] = StripForce.dot(winddir) /SpanResFF.stripArea(m);

            m++;
        }
    }

    Drag.set(ForceBodyAxes);    // N/q, body axes
}


/**
 * +1 or -1: the sense in which a ring of positive doublet density, counter-clockwise about the panel's normal, circulates, in the sense of the strip circulation
 * Gamma of forces(), whose force is V x Gamma T: the trailing edge of the panel, traversed counter-clockwise, runs along (+1) or against (-1) T.
 */
static double ringOrientation(Panel3 const &p3, Vector3d const &T)
{
    Vector3d a = p3.vertexAt(1) - p3.vertexAt(0);
    Vector3d b = p3.vertexAt(2) - p3.vertexAt(0);
    bool bCCW = (a*b).dot(p3.normal()) >= 0.0;
    Vector3d te = bCCW ? p3.vertexAt(2) - p3.vertexAt(1) : p3.vertexAt(1) - p3.vertexAt(2);
    return te.dot(T) >= 0.0 ? 1.0 : -1.0;
}


/**
 * The free-stream Kutta-Joukowski force of the thin sheets' vorticity, at the positions of the vorticity (uav#79).
 *
 * A doublet density mu on a sheet is the vorticity n x grad(mu) over the triangles plus, along every edge, the vortex line of circulation (jump of mu), and
 * along the outer edges the line of circulation mu itself. For the nodal densities Mu3 (linear on each triangle, 3 values per panel) each triangle is then a
 * ring of circulation G.mu(s) along its edges, G = 4 pi per unit doublet density with the sign of its strip, closed by the area vorticity -G n x grad(mu):
 * its total force is zero and its moment is not. Summed over the triangles the jumps along shared edges are left, and the edges between the strips
 * are the chordwise vorticity d(mu)/dy that the strip force V x Gamma l of forces() cannot see, since it puts the whole force at the trailing edge, at
 * the quarter point. The trailing edges are left out: the wake continues the doublet density there, so that no vortex line closes it. The force of the
 * sheets is then the free-stream part of the strips' far-field force (Stokes), to rounding; only its distribution differs, hence the moment.
 *
 * Adds the moment about the CoG (N.m) to Moment and sets Force (N), a check. Thin strips only; the Trefftz-plane part stays on the strips.
 */
void P3Analysis::sheetVorticityLoad(double const *Mu3, std::vector<Vector3d> const &VInf, Vector3d const &CoG, double G4pi, Vector3d &Force, Vector3d &Moment) const
{
    Force.set(0.0, 0.0, 0.0);
    double G = G4pi * m_pPolar3d->density();
    for(int i3=0; i3<nPanels(); i3++)
    {
        Panel3 const &p3 = m_Panel3.at(i3);
        if(!p3.isWingPanel() || p3.isFusePanel()) continue;
        Vector3d const &V = VInf.at(i3);
        Vector3d n = p3.normal();

        // the vertices counter-clockwise about the normal
        int vi[3] = {0, 1, 2};
        {
            Vector3d a = p3.vertexAt(1) - p3.vertexAt(0);
            Vector3d b = p3.vertexAt(2) - p3.vertexAt(0);
            if((a*b).dot(n) < 0.0) {vi[1]=2; vi[2]=1;}
        }
        Vector3d x[3];
        double mu[3];
        for(int k=0; k<3; k++)
        {
            x[k]  = p3.vertexAt(vi[k]);
            mu[k] = Mu3[3*i3+vi[k]];
        }
        double area = ((x[1]-x[0]) * (x[2]-x[0])).dot(n) / 2.0;
        if(area<=0.0) continue;

        // the area vorticity, at the centroid
        Vector3d grad;
        for(int k=0; k<3; k++) grad += (n * (x[(k+2)%3] - x[(k+1)%3])) * mu[k];
        grad *= 1.0/(2.0*area);
        Vector3d centroid = (x[0]+x[1]+x[2])/3.0;
        Vector3d Fa = (V * (n * grad)) * (-G * area);
        Moment += (centroid - CoG) * Fa;
        Force  += Fa;

        // the lines along the edges, but the trailing edge
        for(int k=0; k<3; k++)
        {
            int k1 = (k+1)%3;
            if(p3.isTrailing() && vi[k]+vi[k1]==3) continue;     // vertices 1 and 2 are the trailing nodes
            Vector3d l = x[k1] - x[k];
            double m0 = (mu[k]+mu[k1])/2.0;
            Vector3d lever = (x[k]-CoG)*m0 + l*((mu[k]+2.0*mu[k1])/6.0);   // integral of mu (x-CoG) along the edge
            Vector3d dF = V * l;
            Moment += (lever * dF) * G;
            Force  += dF * (m0 * G);
        }
    }
}


/**
 * The panels of the thick strip whose bottom trailing panel is i3, in the order of the iPU links (bottom surface from the trailing edge, leading edge, top surface
 * to the trailing edge), and its number of chordwise cells a side NX: cell c (panels 2c and 2c+1) of the bottom faces cell 2 NX - 1 - c of the top. False if the
 * chain is not 2 NX cells of two triangles, the first NX on the bottom.
 */
bool P3Analysis::thickStripChain(int i3, std::vector<int> &chain, int &NX) const
{
    chain.clear();
    NX = 0;
    if(i3<0 || i3>=nPanels() || !m_Panel3.at(i3).isBotPanel()) return false;
    for(int k=i3; k>=0; k=m_Panel3.at(k).m_iPU)
    {
        chain.push_back(k);
        if(int(chain.size())>nPanels()) return false;
    }
    if(chain.size()%4!=0 || chain.size()<4) return false;
    NX = int(chain.size())/4;
    for(int j=0; j<int(chain.size()); j++)
        if(m_Panel3.at(chain[j]).isBotPanel() != (j<2*NX)) return false;
    return true;
}


/**
 * The pressure force of the lift of a thick strip (N/rho): per chordwise cell, the panel forces of the bottom cell and of the top cell facing it, panelForce (N/rho
 * per panel), projected on the cell's mean-surface normal, the difference of the top's and the bottom's area-weighted normals. On a thin strip the panel
 * force is that of the pressure jump across the sheet, normal to it; on a thick strip the projection keeps that part, the jump between facing points, and leaves
 * out the in-plane part of the pressures: the leading-edge suction, which the rotating states direct as on a thin strip (stripSuctionCorrection()), and the
 * thickness' own pressures (uav#80). False if the strip's chain is not regular (the caller then takes the whole strip's pressure force).
 */
bool P3Analysis::thickStripLiftPressureForce(int i3, std::vector<Vector3d> const &panelForce, Vector3d &FP) const
{
    FP.set(0.0, 0.0, 0.0);
    std::vector<int> chain;
    int NX = 0;
    if(!thickStripChain(i3, chain, NX)) return false;
    for(int c=0; c<NX; c++)
    {
        int cb[2] = {chain[2*c], chain[2*c+1]};
        int ct[2] = {chain[2*(2*NX-1-c)], chain[2*(2*NX-1-c)+1]};
        Vector3d n, f;
        for(int k=0; k<2; k++)
        {
            Panel3 const &pb = m_Panel3.at(cb[k]);
            Panel3 const &pt = m_Panel3.at(ct[k]);
            n += pt.normal()*pt.area() - pb.normal()*pb.area();
            f += panelForce.at(cb[k]) + panelForce.at(ct[k]);
        }
        if(n.norm()<1.0e-12) return false;
        n.normalize();
        FP += n * f.dot(n);
    }
    return true;
}


/**
 * The moment about QP (N.m/rho) of a thick strip's load normal to its chord plane, cell by cell (m_MomentLift, uav#198): per chordwise cell, the panel forces
 * of the bottom cell and of the top cell facing it, panelForce (N/rho per panel), projected on the normal of the strip's chord plane and placed at the
 * area-weighted centre of the four panels. The chord plane holds the chord, from the trailing edge's midpoint to the quarter-chord point QP, and the span along the
 * trailing edge; it is the plane of the leading-edge suction and of the thickness' in-plane pressures, which stay at QP (P3Analysis::forces(), rule 2).
 * Not the cells' own mean-surface normal: at the nose both faces look forward and the mesh's leading-edge node need not lie on the mean surface, so that
 * normal tilts towards the suction (by 0.05 on the leading-edge cells of a 45 deg NACA 0012 wing, whose leading-edge node is 0.4 mm above the chord) and
 * takes a part of it in. False if the strip's chain is not regular.
 */
bool P3Analysis::thickStripLiftMoment(int i3, std::vector<Vector3d> const &panelForce, Vector3d const &QP, Vector3d &M) const
{
    M.set(0.0, 0.0, 0.0);
    std::vector<int> chain;
    int NX = 0;
    if(!thickStripChain(i3, chain, NX)) return false;
    Panel3 const &p3 = m_Panel3.at(i3);
    Vector3d span = p3.rightTrailingNode() - p3.leftTrailingNode();
    Vector3d chord = QP - (p3.leftTrailingNode() + p3.rightTrailingNode())*0.5;   // its sign does not matter: n enters twice
    Vector3d n = chord * span;
    if(n.norm()<1.0e-12) return false;
    n.normalize();
    for(int c=0; c<NX; c++)
    {
        int cell[4] = {chain[2*c], chain[2*c+1], chain[2*(2*NX-1-c)], chain[2*(2*NX-1-c)+1]};
        Vector3d f, cc;
        double area = 0.0;
        for(int j : cell)
        {
            Panel3 const &pj = m_Panel3.at(j);
            f += panelForce.at(j);
            cc += pj.CoG()*pj.area();
            area += pj.area();
        }
        if(area<=0.0) return false;
        M += (cc/area - QP) * (n * f.dot(n));
    }
    return true;
}


/**
 * The free-stream Kutta-Joukowski force of the lifting vorticity of the thick strips, on their mean surface (uav#80).
 *
 * A thick strip carries the doublet density mu on its bottom and top surfaces. Its vorticity is the lift's and the thickness': the part of mu common to both
 * surfaces, with the sources, makes the thickness, which the free stream alone does not load (on a closed surface the force on that vorticity is that of the
 * mean velocity across the sheet, not of the free stream). The jump mu_top - mu_bottom between facing points is the lift's: it is the doublet density of the
 * mean surface that, with the thickness' sources, has the same far field (thin-wing theory: the lifting problem is the antisymmetric part, Katz and Plotkin,
 * Low-Speed Aerodynamics, 2nd ed., 2001, sec. 4.4 and 12.1). At the trailing edge it is the strip's circulation, at the leading edge 0.
 *
 * Each strip is walked from its bottom trailing panel along the iPU links (bottom surface, leading edge, top surface); its chordwise stations are the pairs of nodes
 * that two consecutive cells share, the station s of the bottom facing the station 2 NX - s of the top. The mean surface is the bottom's triangles with each node
 * moved to the midpoint of itself and its facing node, its nodal density mu_top - mu_bottom (each node's value the mean of the strip's panels that hold it);
 * its load is that of sheetVorticityLoad() on a thin sheet: per triangle a ring of circulation 4 pi mu along its edges, closed by the area vorticity, the trailing
 * edge left out. The sense of each strip's rings is set by its circulation (forces(): the force rho Gamma V x l of the bottom trailing panel's vortex), so that the
 * total is the free-stream part of the strips' far-field force (Stokes): measured within 0.3 percent (the rectangular and the 45 deg swept wing), the nodal
 * densities being means over each strip's own panels.
 *
 * Adds the moment about the CoG (N.m) to Moment, sets Force (N), a check. irregular: the bottom trailing panels of the strips whose mesh is not the expected chain
 * of cells (a null triangle), which add nothing: forces() keeps their free-stream force at their quarter points.
 */
void P3Analysis::meanSheetVorticityLoad(double const *Mu3, std::vector<Vector3d> const &VInf, Vector3d const &CoG, Vector3d &Force, Vector3d &Moment,
                                        std::vector<int> &irregular) const
{
    struct MeanTri {Vector3d x[3]; double mu[3]; bool te[3]; Vector3d V; double G;};
    std::vector<MeanTri> tris;
    Force.set(0.0, 0.0, 0.0);
    irregular.clear();

    // the mean surface of one strip, from its bottom trailing panel i3; false if its mesh is not the expected chain of cells
    auto addStrip = [&](int i3) -> bool
    {
        Panel3 const &p3 = m_Panel3.at(i3);
        std::vector<int> chain;
        int NX = 0;
        if(!thickStripChain(i3, chain, NX)) return false;
        int nCells = 2*NX;

        // the stations: st[s] = {node on the side of the bottom trailing panel's left trailing node, the other}
        std::vector<std::array<int,2>> st(nCells+1);
        st[0] = {p3.leftTrailingNode().index(), p3.rightTrailingNode().index()};
        for(int c=0; c<nCells; c++)
        {
            Panel3 const &t1 = m_Panel3.at(chain[2*c]);
            Panel3 const &t2 = m_Panel3.at(chain[2*c+1]);
            auto nOld = [&](Panel3 const &t){int n=0; for(int k=0; k<3; k++) if(t.nodeIndex(k)==st[c][0] || t.nodeIndex(k)==st[c][1]) n++; return n;};
            Panel3 const *pOne=nullptr, *pTwo=nullptr;
            if     (nOld(t1)==1 && nOld(t2)==2) {pOne=&t1; pTwo=&t2;}
            else if(nOld(t1)==2 && nOld(t2)==1) {pOne=&t2; pTwo=&t1;}
            else return false;
            int diag=-1;
            for(int k=0; k<3; k++) if(pTwo->nodeIndex(k)!=st[c][0] && pTwo->nodeIndex(k)!=st[c][1]) diag = pTwo->nodeIndex(k);
            int xside=-1, other=-1;
            for(int k=0; k<3; k++)
            {
                int n = pOne->nodeIndex(k);
                if(n==st[c][0]) xside = 0;
                else if(n==st[c][1]) xside = 1;
                else if(n!=diag) other = n;
            }
            if(diag<0 || other<0 || xside<0 || !pOne->hasVertex(diag)) return false;
            st[c+1][xside]   = other;
            st[c+1][1-xside] = diag;
        }

        // the nodal densities, each the mean of the strip's panels that hold the node
        auto nodeMu = [&](int node)
        {
            double sum=0.0; int n=0;
            for(int j : chain)
                for(int k=0; k<3; k++)
                    if(m_Panel3.at(j).nodeIndex(k)==node) {sum += Mu3[3*j+k]; n++;}
            return n>0 ? sum/double(n) : 0.0;
        };
        auto nodePos = [&](int node) -> Vector3d
        {
            for(int j : chain)
                for(int k=0; k<3; k++)
                    if(m_Panel3.at(j).nodeIndex(k)==node) return m_Panel3.at(j).vertexAt(k);
            return Vector3d();
        };

        // the mean surface's nodes, station s = 0 (trailing edge) to NX (leading edge)
        std::vector<std::array<Vector3d,2>> mx(NX+1);
        std::vector<std::array<double,2>> mmu(NX+1);
        for(int s=0; s<=NX; s++)
            for(int side=0; side<2; side++)
            {
                int nb = st[s][side], nt = st[2*NX-s][side];
                mx[s][side]  = (nodePos(nb) + nodePos(nt)) * 0.5;
                mmu[s][side] = nb==nt ? 0.0 : nodeMu(nt) - nodeMu(nb);
            }

        // the bottom's triangles on the mean surface, oriented as the bottom's normal reversed (towards the top)
        double sense = 0.0;
        size_t first = tris.size();
        for(int j=0; j<2*NX; j++)
        {
            Panel3 const &pb = m_Panel3.at(chain[j]);
            int c = j/2;
            MeanTri t;
            for(int k=0; k<3; k++)
            {
                int n = pb.nodeIndex(k);
                int s=-1, side=-1;
                for(int ss=c; ss<=c+1; ss++)
                    for(int sd=0; sd<2; sd++)
                        if(st[ss][sd]==n) {s=ss; side=sd;}
                if(s<0) return false;
                t.x[k]  = mx[s][side];
                t.mu[k] = mmu[s][side];
                t.te[k] = s==0;
            }
            Vector3d nrm = (t.x[1]-t.x[0]) * (t.x[2]-t.x[0]);
            if(nrm.dot(pb.normal())>0.0)
            {
                std::swap(t.x[1], t.x[2]);
                std::swap(t.mu[1], t.mu[2]);
                std::swap(t.te[1], t.te[2]);
            }
            t.V = VInf.at(chain[j]);
            t.G = 0.0;
            for(int k=0; k<3; k++)
            {
                int k1 = (k+1)%3;
                if(t.te[k] && t.te[k1]) sense = (t.x[k1]-t.x[k]).dot(p3.trailingVortex())>=0.0 ? 1.0 : -1.0;
            }
            tris.push_back(t);
        }
        if(sense==0.0) return false;
        for(size_t k=first; k<tris.size(); k++) tris[k].G = sense * 4.0*PI * m_pPolar3d->density();
        return true;
    };

    for(int i3=0; i3<nPanels(); i3++)
    {
        Panel3 const &p3 = m_Panel3.at(i3);
        if(!p3.isTrailing() || !p3.isBotPanel()) continue;
        size_t n0 = tris.size();
        if(!addStrip(i3))
        {
            tris.resize(n0);
            irregular.push_back(i3);
        }
    }

    // the load of each triangle, as sheetVorticityLoad()
    for(MeanTri const &t : tris)
    {
        Vector3d n = (t.x[1]-t.x[0]) * (t.x[2]-t.x[0]);
        double area = n.norm()/2.0;
        if(area<=0.0) continue;
        n *= 1.0/(2.0*area);
        Vector3d grad;
        for(int k=0; k<3; k++) grad += (n * (t.x[(k+2)%3] - t.x[(k+1)%3])) * t.mu[k];
        grad *= 1.0/(2.0*area);
        Vector3d centroid = (t.x[0]+t.x[1]+t.x[2])/3.0;
        Vector3d Fa = (t.V * (n * grad)) * (-t.G * area);
        Moment += (centroid - CoG) * Fa;
        Force  += Fa;
        for(int k=0; k<3; k++)
        {
            int k1 = (k+1)%3;
            if(t.te[k] && t.te[k1]) continue;
            Vector3d l = t.x[k1] - t.x[k];
            double m0 = (t.mu[k]+t.mu[k1])/2.0;
            Vector3d lever = (t.x[k]-CoG)*m0 + l*((t.mu[k]+2.0*t.mu[k1])/6.0);
            Vector3d dF = t.V * l;
            Moment += (lever * dF) * t.G;
            Force  += dF * (m0 * t.G);
        }
    }
}


/**
* Calculates the forces using a far-field method.
* Calculates the moments by a near field method, i.e. direct summation on the panels.
* The downwash is evaluated far downstream (i.e. where the inluence of the bounded vortices is negligible)
*
* On thin surfaces the panel force is normal to the mean surface, so the panel pressures carry no leading-edge suction,
* and the load near the leading edge is under-integrated; the far-field force carries both. Strip by strip, the difference
* between the far-field force (the Kutta-Joukowski force with the local velocity plus the Trefftz-plane force) and the sum of the
* strip's panel forces is a force that belongs at the strip's leading edge: its moment is added to the on-body moment.
* Then the moment and the force agree in what they include, and the neutral point does not depend on the CG.
* Thick surfaces need no correction: their pressures carry the suction.
*
* The force returned is the sum of the strips' forces at their trailing wake points. If m_bFarField is set, forces() also sums the
* far-field strip forces, at the mid-wake point as the force of the polar, in m_ForceFF, and their moment about the CoG in m_MomentFF.
* A trailing wake point that lies in another surface's wake sheet (a stabiliser at the wing's height) does not see the far-field velocity:
* the derivatives taken from it are wrong, so the derivatives in the plane of symmetry (Xu, Zu, Xw, Zw, Xq, Zq) and the control
* derivatives take m_ForceFF. The sideslip (v) derivatives take it too, with the Trefftz plane perpendicular to the free stream, along which the wake of a
* sideslipped state leaves the trailing edges (m_bSideslipFarField); with the plane perpendicular to the body x axis the far-field velocity is not
* force-free under a lateral perturbation (issue #7).
*
* In the sideslip derivatives m_MomentFF puts the free-stream part of the strips' forces on the vorticity of the sheet instead of the quarter points
* (sheetVorticityLoad(), thin surfaces; uav#79): the strips' force V x Gamma l_TE is right as a total (Stokes), but lumps at the quarter point the chordwise
* vorticity d(mu)/dy between the strips, whose force in a sideslip carries the rolling and yawing moments. Thick strips put it on their lifting vorticity,
* the jump mu_top - mu_bottom on their mean surface (meanSheetVorticityLoad(), uav#80).
*
* In rotating states (m_bRotationSuction) each strip's leading-edge suction is directed perpendicular to the leading edge, in the forces and the moments:
* without rotation the strips' far-field force is exact in total, with it the near-field force is the reference (flow5#28).
*
* m_MomentFF applies each strip's far-field force at the strip's quarter-chord point, so it has no chordwise load. m_MomentLocal is the
* moment about these points of the strips' own chordwise loads: the panel pressures and, on thin surfaces, the leading-edge force. The
* pitch damping of a wing needs it, since its load moves along the chord in a pitch rate.
*/
void P3Analysis::forces(double const *Mu3, double const *Sigma3, double alpha, double beta, Vector3d const &CoG, bool bFuseMi,
                        std::vector<Vector3d> const &VInf, Vector3d &Force, Vector3d &Moment)
{
    if(!m_pPolar3d) return;

    double Cp=0.0, viscousDrag=0.0;
    double QInf=1.0;
    Vector3d left, mid, right;
    Vector3d WindDirection, WindNormal, PanelLeverArm, Wg;
    Vector3d Velocity, stripforce, viscousMoment, PanelForce;

    //   Define the wind axis
    WindNormal = objects::windNormal(alpha, beta);
    WindDirection = objects::windDirection(alpha, beta);

    Force.set( 0.0, 0.0, 0.0);
    Moment.set(0.0, 0.0, 0.0);
    m_ForceFF.set(0.0, 0.0, 0.0);
    m_MomentFF.set(0.0, 0.0, 0.0);
    m_MomentLocal.set(0.0, 0.0, 0.0);
    m_MomentLift.set(0.0, 0.0, 0.0);
    viscousDrag = 0.0;
    viscousMoment.set(0.0,0.0,0.0);

    // the strips' far-field forces (N), their quarter-chord and leading edge points, the sums of their panel forces (N/rho),
    // whether their leading-edge force is to be added (thin strips), and the strip that each panel belongs to (-1 for none: fuselage, tips)
    std::vector<Vector3d> stripFF, stripQP, stripLE, stripFP;
    std::vector<double> stripG;   // sign of the circulation of a counter-clockwise ring of a unit doublet density, in 4 pi units: the sideslip vorticity's orientation
    std::vector<bool> stripThin;
    std::vector<int> stripOfPanel(nPanels(), -1);
    std::vector<Vector3d> stripNF;   // the strips' near-field forces (N), which Force sums
    std::vector<int> stripTP;        // the strips' trailing panels (mid panel, or bottom panel of a thick strip)

    // The far-field velocity of the strips is the one of the polar's induced force (trefftzDrag): the analytic Trefftz plane perpendicular to the wake sheet,
    // from the trailing legs of every strip of every surface for this solution; half of it is the velocity at the bound vortex. It does not depend on the
    // length of the wake.
    bool bSheetLoad = m_bFarField && m_bSideslipVorticity;   // sideslip derivatives (uav#79; thick surfaces uav#80)
    std::vector<std::pair<int, Vector3d>> thickFreeStream;   // the thick strips' bottom trailing panels and the free-stream parts of their far-field forces (N), with bSheetLoad

    Vector3d wakedir(1.0, 0.0, 0.0);
    if(m_bSideslipFarField && !VInf.empty() && VInf.front().norm()>0.0) wakedir = VInf.front()/VInf.front().norm();   // the sideslip derivatives' wake leaves along the free stream
    std::vector<TrefftzStrip> trefftzStrips;
    makeTrefftzStrips(Mu3, wakedir, trefftzStrips);
    auto trefftzWash = [&](Panel3 const &p3)
    {
        Vector3d pt = (p3.leftTrailingNode() + p3.rightTrailingNode())/2.0;
        pt -= wakedir * pt.dot(wakedir);
        return trefftzVelocity(pt, wakedir, trefftzStrips) * 0.5;
    };

//    int m=0;

    for(int i3=0; i3<nPanels(); i3++)
    {
        Panel3 const &p3 = m_Panel3.at(i3);

        if(p3.isMidPanel())
        {
            if(p3.isTrailing())
            {
                // get the last triangle of the wake column
                assert(p3.iWake()>=0 && p3.iWake()<nWakePanels());
                Panel3 *p3W = m_WakePanel3.data() + p3.iWake();
                trailingWakePoint(p3W, left, right);
                mid = (left + right)/2.0;

                //    s_PointDbg.push_back(C);
                getVelocityVector(mid, Mu3, Sigma3, Wg, 0.00001, true, true);

                int idxM = p3.index();
                double GammaStrip = -(Mu3[3*idxM+1] + Mu3[3*idxM+2])/2.0 *4.0*PI;
                //qDebug("m={:3d}  %13.5f  %13.5f  %13.5f  %13.5f", m, GammaStrip, Wg.x, Wg.y, Wg.z);

                Wg += VInf.at(i3);

                Vector3d vortex = p3.trailingVortex();
                stripforce  =  Wg * vortex;
                stripforce *= GammaStrip * m_pPolar3d->density();     // N
                Force += stripforce;
                Vector3d nearforce = stripforce;

                // the strip's far-field force: the Kutta-Joukowski lift plus the
                // Trefftz-plane force, evaluated at the mid-wake point as in trefftzDrag;
                // the trailing-point velocity used for the force above is left out: the wake-induced
                // velocity is not force-free under a lateral perturbation (issue #7)
                // The thin strips need it too: their leading-edge force is the difference of this force and the panel forces
                {
                    Vector3d WgFF = trefftzWash(p3);
                    Vector3d wash = WgFF;
                    WgFF += VInf.at(i3);
                    stripforce  =  WgFF * vortex;
                    stripforce *= GammaStrip * m_pPolar3d->density();     // N
                    m_ForceFF += stripforce;
                    Vector3d quarterpoint = stripQuarterPoint(p3);
                    if(m_bFarField)
                    {
                        if(bSheetLoad)
                        {
                            // only the Trefftz-plane part is lumped at the quarter point: the free-stream part is put on the vortex lines of the sheet below
                            Vector3d inducedforce = wash * vortex;
                            inducedforce *= GammaStrip * m_pPolar3d->density();     // N
                            m_MomentFF += (quarterpoint - CoG) * inducedforce;   // N.m
                        }
                        else m_MomentFF += (quarterpoint - CoG) * stripforce;   // N.m
                    }

                    // walk the strip upstream, to its leading edge
                    int iStrip = int(stripFF.size());
                    Panel3 const *pLast = &p3;
                    stripOfPanel[i3] = iStrip;
                    while(pLast->m_iPU>=0)
                    {
                        stripOfPanel[pLast->m_iPU] = iStrip;
                        pLast = &m_Panel3.at(pLast->m_iPU);
                    }
                    stripFF.push_back(stripforce);
                    stripQP.push_back(quarterpoint);
                    stripLE.push_back(pLast->leadingEdgeMidpoint());
                    stripFP.push_back(Vector3d());
                    stripThin.push_back(true);
                    stripG.push_back(ringOrientation(p3, vortex) * 4.0*PI);
                    stripNF.push_back(nearforce);
                    stripTP.push_back(i3);
                }

 //               m++;
            }
        }
        else if(p3.isTrailing() && p3.isBotPanel())
        {
            int idxB = p3.index();
            int idxU = nextTopTrailingPanelIndex(p3);
            assert(idxU>=0);

            //Get the strip's lifting force

            // get the last triangle of the wake column
            Panel3 *p3W = m_WakePanel3.data() + p3.iWake();
            trailingWakePoint(p3W, left, right);
            mid = (left + right)/2.0;

            /** @todo used to be bWakeOnly=false up to beta07 */
            getVelocityVector(mid, Mu3, Sigma3, Wg, Vortex::coreRadius(), true, true);

            Wg += VInf.at(i3);

            // the trailing nodes have indexes 1 & 2
            double gU = (Mu3[3*idxU+1] + Mu3[3*idxU+2])/2.0;
            double gB = (Mu3[3*idxB+1] + Mu3[3*idxB+2])/2.0;
            double GammaStrip = (-gU + gB) *4.0*PI;

            Vector3d vortex = p3.trailingVortex()*(-1.0);
            stripforce  = vortex * Wg;
            stripforce *= GammaStrip * m_pPolar3d->density();     // N

            Force += stripforce;
            Vector3d nearforce = stripforce;

            // same far-field force as in the mid-panel branch above
            Vector3d quarterpoint = stripQuarterPoint(p3);
            stripforce.set(0.0, 0.0, 0.0);
            if(m_bFarField)
            {
                Vector3d WgFF = trefftzWash(p3);
                Vector3d wash = WgFF;
                WgFF += VInf.at(i3);
                stripforce  = vortex * WgFF;
                stripforce *= GammaStrip * m_pPolar3d->density();     // N
                m_ForceFF += stripforce;
                if(bSheetLoad)
                {
                    // as the thin strips: the Trefftz-plane part at the quarter point, the free-stream part on the mean surface's vorticity below
                    Vector3d inducedforce = vortex * wash;
                    inducedforce *= GammaStrip * m_pPolar3d->density();     // N
                    m_MomentFF += (quarterpoint - CoG) * inducedforce;   // N.m
                    thickFreeStream.push_back({i3, stripforce - inducedforce});
                }
                else m_MomentFF += (quarterpoint - CoG) * stripforce;   // N.m
            }

            // the panels of the strip: the chain of the iPU links from the bottom trailing panel runs along the bottom surface
            // and round the leading edge to the top trailing panel
            int iStrip = int(stripFF.size());
            Panel3 const *pLast = &p3;
            stripOfPanel[i3] = iStrip;
            while(pLast->m_iPU>=0)
            {
                stripOfPanel[pLast->m_iPU] = iStrip;
                pLast = &m_Panel3.at(pLast->m_iPU);
            }
            stripFF.push_back(stripforce);
            stripQP.push_back(quarterpoint);
            stripLE.push_back(quarterpoint);
            stripFP.push_back(Vector3d());
            stripThin.push_back(false);
            stripG.push_back(ringOrientation(m_Panel3.at(idxU), p3.trailingVortex()) * 4.0*PI);   // the ring of the top trailing panel
            stripNF.push_back(nearforce);
            stripTP.push_back(i3);
//            m++;
        }
    }

    // The free-stream Kutta-Joukowski force on the vorticity of the sheets, at the vortex lines' own positions (sideslip derivatives)
    if(bSheetLoad && !stripG.empty())
    {
        Vector3d sheetForce;
        if(!thickFreeStream.empty())
        {
            // thick strips: the lifting vorticity on the mean surface (uav#80); a strip whose mesh is not the regular chain of cells keeps its free-stream
            // force at its quarter point, as the older model
            std::vector<int> irregular;
            meanSheetVorticityLoad(Mu3, VInf, CoG, sheetForce, m_MomentFF, irregular);
            for(int i : irregular)
                for(auto const &fs : thickFreeStream)
                    if(fs.first==i) m_MomentFF += (stripQuarterPoint(m_Panel3.at(i)) - CoG) * fs.second;   // N.m
        }
        else sheetVorticityLoad(Mu3, VInf, CoG, stripG.front(), sheetForce, m_MomentFF);
    }

    //On-Body moment

    Moment.set(0.0,0.0,0.0);

    // the rotating states' suction rule and m_MomentLift need the panel forces of the thick strips (thickStripLiftPressureForce())
    bool bThickFF = false;
    if(m_bFarField)
        for(bool b : stripThin) if(!b) bThickFF = true;
    std::vector<Vector3d> panelForce(bThickFF ? nPanels() : 0);

    for(int i3=0; i3<nPanels(); i3++)
    {
        Panel3 const &p3 = m_Panel3.at(i3);
        if(!p3.isFusePanel() || bFuseMi)
        {
            Velocity = VInf.at(i3);
            QInf = Velocity.norm();
            Cp = (m_Cp[3*i3]+m_Cp[3*i3+1]+m_Cp[3*i3+2])/3.0;
            PanelForce = p3.normal() * (-Cp) * p3.area() *1/2.*QInf*QInf;      // Newtons/rho
            PanelLeverArm = p3.CoG() - CoG;
            Moment += PanelLeverArm * PanelForce;                     // N.m/rho
            if(bThickFF) panelForce[i3] = PanelForce;
            if(stripOfPanel[i3]>=0)
            {
                int iStrip = stripOfPanel[i3];
                stripFP[iStrip] += PanelForce;   // N/rho
                m_MomentLocal += (p3.CoG() - stripQP[iStrip]) * PanelForce * m_pPolar3d->density();   // N.m
                if(stripThin[iStrip]) m_MomentLift += (p3.CoG() - stripQP[iStrip]) * PanelForce * m_pPolar3d->density();   // N.m: a thin panel's force is normal to the sheet
            }
            // the far-field strip forces only cover the wings: the fuse keeps its on-body pressure moment
            if(m_bFarField && p3.isFusePanel()) m_MomentFF += PanelLeverArm * PanelForce * m_pPolar3d->density();  // N.m
        }
    }

    if(m_pPolar3d->isViscous())
    {
        Force += WindDirection * viscousDrag;
        Moment += viscousMoment;
    }

    // add fuse drag and extra drag
/*    double frictiondrag = m_pPolar3d->extraDragTotal();
    if(m_pPlane->fuse(0) && pWPolar->isFuseDragIncluded())
    {
        Re = m_pPlane->fuse(0)->length() * m_QInfMin / m_pPolar3d->viscosity();
        double fusedrag = pWPolar->fuseDragCoef(Re) * m_pPlane->fuse(0)->wettedArea() * m_pPlane->fuse(0)->formFactor();
        frictiondrag += fusedrag;
    }
    frictiondrag *= 1./2.*m_pPolar3d->density()*Vinc.norm()*Vinc.norm();   // N

    Force += WindDirection*frictiondrag;
*/
    Moment *= m_pPolar3d->density();                          // N.m

    // The leading-edge suction of every strip, thin or thick, perpendicular to the leading edge (stripSuctionCorrection(), flow5#28), in rotating
    // states only (m_bRotationSuction: the rate derivatives, the steady rotating points). Three rules, each with its basis and its limit:
    // 1. Rotation only. Basis: with a uniform onset (alpha, beta, controls) the sum of the strips' forces rho Gamma (V + w/2) x l_TE is the
    //    Trefftz-plane force, exact in total by momentum whatever the direction of each strip's segment; AVL's near-field and Trefftz-plane CY agree
    //    in a sideslip (-0.0024 and -0.0025 on a wing with 8 deg dihedral) and not in a roll. Applying the change there too moved CL and CD of a
    //    bell-twisted wing by -0.6 and -3.4 percent and CYb of the dihedral wing 12 percent away from AVL. Limit: momentum is exact only for a
    //    force-free wake; flow5's wake is the straight one of the polar. A rotating onset V + Omega x r is not uniform, the far field is no longer
    //    the total, and the near-field force of thin-wing theory is taken.
    // 2. The change sits where the strip's force sits in each moment: the quarter-chord point in m_MomentFF, the leading-edge force in the pitch
    //    damping's m_MomentLocal; the steady rotating points take their moment's change the same way (stateSuctionChange()), so that a point and
    //    the derivatives agree. Basis: m_MomentFF lumps every strip's far-field force there; moving the suction alone to the leading edge, its true
    //    point of action, turned Cnp/CL positive on swept wings of aspect ratio 1.34 (+0.012 and +0.032 against AVL's -0.018 and -0.058 at 100
    //    chordwise vortices, the wind tunnel's -0.38 and -0.56). Limit: Cnp and Cnr carry no leading-edge lever for the suction until the
    //    near-field moment (flow5#29) puts every load where it acts.
    // 3. Rotation with sideslip (a type-8 point with both): the sideslip's change (PlaneTask::applySideslipModel()) and the rotation's
    //    (stateSuctionChange()) are each made without the other and added. Limit: superposition, exact to first order; the cross term of beta and
    //    the rates is left out.
    if(m_bRotationSuction && m_bFarField)
    {
        for(unsigned int is=0; is<stripFF.size(); is++)
        {
            Panel3 const &p3 = m_Panel3.at(stripTP[is]);
            Vector3d FP = stripFP[is] * m_pPolar3d->density();   // N
            Vector3d FPlift;
            if(!stripThin[is] && thickStripLiftPressureForce(stripTP[is], panelForce, FPlift)) FP = FPlift * m_pPolar3d->density();   // uav#80
            Force += stripSuctionCorrection(p3, stripNF[is], FP);
            Vector3d dF = stripSuctionCorrection(p3, stripFF[is], FP);
            m_ForceFF   += dF;
            stripFF[is] += dF;
            m_MomentFF  += (stripQP[is] - CoG) * dF;   // N.m
        }
    }

    // the leading-edge force of the thin strips: what the far field carries and the panel forces do not
    for(unsigned int is=0; is<stripFF.size(); is++)
    {
        if(!stripThin[is]) continue;
        Vector3d leforce = stripFF[is] - stripFP[is] * m_pPolar3d->density();   // N
        Moment += (stripLE[is] - CoG) * leforce;   // N.m
        m_MomentLocal += (stripLE[is] - stripQP[is]) * leforce;   // N.m
    }

    // m_MomentLift (uav#198): the strips' load normal to the mean surface where it acts, the in-plane suction left at the quarter-chord point (rule 2 above).
    // Basis: the far-field strip forces at the quarter-chord points carry no chordwise lever, and a rate whose onset varies along the chord (a yaw rate on a
    // fin, as a pitch rate on a wing) loads the strip aft of it: a fin's Cnr read 0.84 of AVL's, its Cnp 1.08. The normal load is the pressure jump, at the
    // panels, and on thin strips the part the panels under-integrate near the leading edge, the normal part of the leading-edge force, at the leading edge,
    // as in Cm and Mq, the normal that of the suction's frame (stripSuctionFrame()), so the suction's change (rule 1), in that plane, does not move it.
    // Thick strips: the cells' load normal to the strip's chord plane, each at its cell (thickStripLiftMoment()), the in-plane pressures (suction and
    // thickness) at the quarter-chord point.
    if(m_bFarField)
    {
        for(unsigned int is=0; is<stripFF.size(); is++)
        {
            Panel3 const &p3 = m_Panel3.at(stripTP[is]);
            if(stripThin[is])
            {
                Vector3d n, e, t, d, ptLE;
                if(!stripSuctionFrame(p3, n, e, t, d, ptLE)) n = p3.normal();
                Vector3d leforce = stripFF[is] - stripFP[is] * m_pPolar3d->density();   // N
                m_MomentLift += (stripLE[is] - stripQP[is]) * (n * leforce.dot(n));   // N.m
            }
            else
            {
                Vector3d Mlift;
                if(thickStripLiftMoment(stripTP[is], panelForce, stripQP[is], Mlift))
                    m_MomentLift += Mlift * m_pPolar3d->density();   // N.m; an irregular chain adds nothing: its load stays at the quarter-chord point
            }
        }
    }
}


/**
 * Calculates the forces using a far-field method.
 * Calculates the moments by a near field method, i.e. direct summation on the panels.
 * Downwash is evaluated at a distance 100 times the span downstream (i.e. infinite)
 * @param Mu a pointer to the array of doublet strengths or vortex circulations
 * @param Sigma a pointer to the array of source strengths
 * @param *VInf a pointer to the array of the velocity vectors on the panels
 * @param Force the resulting force vector
 * @param Moment the resulting moment vector
 */
void P3Analysis::moments(double *Mu3, double alpha, double beta, Vector3d CoG,
                         double *VInf, Vector3d &Force, Vector3d &Moment)
{
    if(!m_pPolar3d) return;

    double Cp=0.0, viscousDrag=0.0;
    double QInf=0.0;
    Vector3d WindDirection, WindNormal, PanelLeverArm, Wg;
    Vector3d Velocity, viscousMoment, PanelForce;

    //   Define the wind axis
    WindNormal = objects::windNormal(alpha, beta);
    WindDirection = objects::windDirection(alpha, beta);

    Force.set( 0.0, 0.0, 0.0);
    Moment.set(0.0, 0.0, 0.0);
    viscousDrag = 0.0;
    viscousMoment.set(0.0,0.0,0.0);


    //On-Body moment
    Vector3d Vp3;
    Moment.set(0.0,0.0,0.0);

    for(int i3=0; i3<nPanels(); i3++)
    {
        Panel3 *p3 = m_Panel3.data()+i3;
        Velocity.x = *(VInf               +i3);
        Velocity.y = *(VInf +   nPanels() +i3);
        Velocity.z = *(VInf + 2*nPanels() +i3);
        QInf = Velocity.norm();
        makePanelDoubletSurfaceVelocity(i3, Mu3, Vp3);

        Vp3 += p3->globalToLocal(Velocity);
        Cp = 1.0 - (Vp3.x*Vp3.x + Vp3.y*Vp3.y + Vp3.z*Vp3.z)/QInf/QInf;

        PanelForce = p3->normal() * (-Cp) * p3->area() *1/2.*QInf*QInf;      // Newtons/rho
        PanelLeverArm = p3->CoG() - CoG;

        Force += PanelForce;
        Moment += PanelLeverArm * PanelForce;                     // N.m/rho
    }

    if(m_pPolar3d->isViscous())
    {
        Force += WindDirection * viscousDrag;
        Moment += viscousMoment;
    }

    Force  *= m_pPolar3d->density();                          // N
    Moment *= m_pPolar3d->density();                          // N.m
}


/**
 * Returns the geometric pitching moment coefficient for the specified angle of attack
 * The effect of the viscous drag is not included.
 * This is the moment of forces(), i.e. with the leading-edge force of the thin surfaces, so that the equilibrium angle
 * is the one where the operating point's Cm is zero.
 * @param Alpha the aoa for which Cm is calculated
 */
double P3Analysis::computeCm(Vector3d const &CoG, double Alpha, bool bFuseMi)
{
    // Define the wind axis
    double beta=0;
    double cosa = cos(Alpha*PI/180.0);
    double sina = sin(Alpha*PI/180.0);
    Vector3d WindDirection(cosa, 0.0, sina);

    std::vector<Vector3d> VLocal(nPanels()*3);
    std::vector<Vector3d> VInf(nPanels(), WindDirection);
    combineLocalVelocities(Alpha, beta, VLocal);
    computeOnBodyCp(VInf, VLocal, m_Cp);

    // forces() needs the doublet and source strengths for this angle
    makeSourceStrengths(WindDirection);
    makeUnitDoubletStrengths(Alpha, beta);

    Vector3d Force, Moment;
    forces(m_Mu.data(), m_Sigma.data(), Alpha, beta, CoG, bFuseMi, VInf, Force, Moment);

    return 2.0 * Moment.y;      // forces() has the dynamic pressure of a unit speed, N.m
}


/**
* Performs the summation of on-body forces to calculate the total lift and drag forces
* @param Cp a pointer to the array of previously calculated Cp coefficients
* @param Alpha the aoa for which this calculation is performed
* @param Lift the resulting lift force
* @param Drag the resulting drag force
*/
void P3Analysis::sumPanelForces(double *Cp, Vector3d &F)
{
    for(int i3=0; i3<nPanels(); i3++)
    {
        F += m_Panel3.at(i3).normal() * (-Cp[3*i3]) * m_Panel3.at(i3).area();
    }
}


void P3Analysis::makeMu(int )
{
/*    for(int i3=0; i3<nPanels(); i3++)
    {
        m_Mu[qrhs*nPanels() +i3] = m_uRHS[i3];
    }*/
}

#define CM_ITER_MAX 50
/**
* Finds the zero-pitching-moment aoa such that Cm=0.
* Proceeds by iteration between -PI/4 and PI/4
* @return true if an equlibrium angle was found false otherwise.
*/
bool P3Analysis::getZeroMomentAngle(Vector3d const &CoG, double &alphaeq, bool bFuseMi)
{
    double tmp=0;
    double eps = 1.e-7;

    int iter = 0;
    double a0 = -PI/4.0;
    double a1 =  PI/4.0;

    double a = 0.0;
    double Cm0 = computeCm(CoG, a0*180.0/PI, bFuseMi);
    double Cm1 = computeCm(CoG, a1*180.0/PI, bFuseMi);
    double Cm = 1.0;

    //are there two initial values of opposite signs?
    int itermax = 50;
    while(Cm0*Cm1>0.0 && iter <=itermax)
    {
        a0 *=0.9;
        a1 *=0.9;
        Cm0 = computeCm(CoG, a0*180.0/PI, bFuseMi);
        Cm1 = computeCm(CoG, a1*180.0/PI, bFuseMi);
//        qDebug(" iter={:3d},  {:11g}   {:11g},  {:11g}   {:11g}", iter, a0*180.0/PI, a1*180.0/PI, Cm0, Cm1);
        iter++;
        if(isCancelled()) break;
    }
    if(iter>=itermax || isCancelled()) return false;

    iter = 0;

    //Cm0 and Cm1 are of opposite sign
    if(Cm0>Cm1)
    {
        tmp = Cm1;
        Cm1 = Cm0;
        Cm0 = tmp;
        tmp = a0;
        a0  = a1;
        a1  = tmp;
    }

    while (fabs(Cm)>eps && iter<=CM_ITER_MAX)
    {
        a = a0 - (a1-a0) * Cm0/(Cm1-Cm0);
        Cm = computeCm(CoG, a*180.0/PI, bFuseMi);
        if(Cm>0.0)
        {
            a1  = a;
            Cm1 = Cm;
        }
        else
        {
            a0  = a;
            Cm0 = Cm;
        }
        iter++;
        if(isCancelled()) break;
    }

    if(iter>=CM_ITER_MAX || isCancelled()) return false;

    computeCm(CoG, a*180.0/PI, bFuseMi);
    alphaeq = a*180.0/PI;

    return true;
}


bool P3Analysis::computeTrimmedConditions(double mass, Vector3d const &CoG, double &alphaeq, double &u0, bool bFuseMi)
{
    Vector3d VInf = objects::windDirection(alphaeq, 0.0);
    Vector3d WindNormal = objects::windNormal(alphaeq, 0.0);

    if(isCancelled()) return false;

    if(!getZeroMomentAngle(CoG, alphaeq, bFuseMi))
    {
        traceStdLog("      no zero-moment angle found...     \n");
        return false;
    }


    makeSourceStrengths(objects::windDirection(alphaeq, 0.0));
    if (isCancelled()) return true;

    //reconstruct doublet strengths from unit cosine and sine vectors
    makeUnitDoubletStrengths(alphaeq, 0.0);
    if(isCancelled()) return false;

    //______________________________________________________________________________________
    // Calculate the trimmed conditions for this control setting and calculated Alpha_eq
/*        phi = bank angle
        V   = sqrt(2 m g / rho S CL cos(phi))    (airspeed)
        R   = V^2 / g tan(phi)        (turn radius, positive for right turn)
        W   = V / R                   (turn rate, positive for right turn)
        p   = 0                       (roll rate, zero for steady turn)
        q   = W sin(phi)              (pitch rate, positive nose upward)
        r   = W cos(phi)              (yaw rate, positive for right turn) */

    // so far we have a unit Vortex Strength
    // find the speeds which will create a lift equal to the weight

    u0 = 1.0;

    //extra drag is unused when calculating lifting velocity
//displayArray(m_Mu);

    Vector3d Force, Moment;
    std::vector<Vector3d> VField(nPanels(), VInf);
    forces(m_Mu.data(), m_Sigma.data(), alphaeq, 0.0, CoG, bFuseMi, VField, Force, Moment);

    double Lift   = Force.dot(WindNormal);        //N/rho ; bank effect not included
//    VerticalCl = Lift*2.0/m_pWPolar->referenceArea() * cos(phi)/m_pWPolar->density();
    if(Lift<=PRECISION)
    {
        u0 = -100.0;
        std::string strong;
        strong = "        Found a negative lift for " + ALPHAstr;
        strong = std::format("=%.5f", alphaeq) + DEGstr;
        strong = ".... skipping the angle...\n";
        m_bWarning = true;
        traceStdLog("\n"+strong);
        return false;
    }
    else
    {
        u0 = sqrt(9.81 * mass / Force.z);

        if(fabs(m_pPolar3d->phi())>PRECISION)
        {
        }
    }
    return true;
}


void P3Analysis::makeInfluenceMatrix()
{
    m_bMatrixError = false;

    bool bMultiThread = true;

    if(bMultiThread)
    {
        std::vector<std::thread> threads;

        for(int iBlock=0; iBlock<m_nBlocks; iBlock++)
        {
            threads.push_back(std::thread(&P3Analysis::makeMatrixBlock, this, iBlock));
        }
        for(int iBlock=0; iBlock<m_nBlocks; iBlock++)
        {
            threads[iBlock].join();
        }
//        std::cout << "P3Analysis::makeInfluenceMatrix joined all " << m_nBlocks << " threads" <<std::endl;

    }
    else
    {
        for(int iBlock=0; iBlock<m_nBlocks; iBlock++)
        {
            if(m_bMatrixError) break; // no point to continue
            makeMatrixBlock(iBlock);
        }
    }


    if(m_bMatrixError)
    {
        traceStdLog("      Error building the influence matrix\n");
    }
}


/**
 * UNUSED
 * Makes the array of negating vortices at the downstream end of the wake panels.
 * Used to cancel the effect of the transverse vortex between the last wake panel and the array of vortons.
 * The circulation of the vortices is stored in the SpanDistrib results
 * Cf. Willis 2005 fig.3
 */
void P3Analysis::makeNegatingVortices(std::vector<Vortex> &negvortices)
{
/*    s_DebugPts.clear();
    s_DebugVecs.clear();*/
    Vector3d left, right;
    negvortices.clear();
    for(int iw=0; iw<nWakePanels(); iw++)
    {
        Panel3 const &p3w = m_WakePanel3.at(iw);
        if(p3w.m_iPU<0)
        {
            // start of the wake column
            trailingWakePoint(&p3w, left, right);
            negvortices.push_back({left, right});
/*            s_DebugPts.push_back(left);
            s_DebugVecs.push_back(negvortices.back().segment());*/
        }
    }
}


void P3Analysis::getDoubletInfluence(Vector3d const &C, Panel3 const &p3, Vector3d *V, double *phi,
                                     double , bool bUseRFF) const
{
    if(!phi && !V) return;

    if(phi) p3.doubletBasisPotential(C, p3.isPointInTriangle(C), phi, bUseRFF);
    if(V)
    {
        if(m_pPolar3d->isTriLinearMethod())
            p3.doubletBasisVelocity(C, V, bUseRFF);
        else
        {
            // faster
            p3.doubletN4023Velocity(C, false, *V, Vortex::coreRadius(), bUseRFF);
            V[0] = V[1] = V[2] = V[0]/3.0;
        }
    }

    if(m_pPolar3d->bHPlane())
    {
        double coef = m_pPolar3d->bGroundEffect() ? 1.0 : -1.0;

        double phiG[]{0,0,0};
        Vector3d VG[3];
        Vector3d CG(C.x, C.y, -C.z-2.0*m_pPolar3d->groundHeight());

        if(phi) p3.doubletBasisPotential(CG, false, phiG, bUseRFF);
        if(V)
        {
            if(m_pPolar3d->isTriLinearMethod())
                p3.doubletBasisVelocity(CG, VG, bUseRFF);
            else
            {
                // faster
                p3.doubletN4023Velocity(CG, false, *VG, Vortex::coreRadius(), bUseRFF);
                VG[0] = VG[1] = VG[2] = VG[0]/3.0;
            }
        }


        for(int ig=0; ig<3; ig++)
        {
            if(phi) phi[ig] += phiG[ig] * coef;
            if(V)
            {
                V[ig].x += VG[ig].x * coef;
                V[ig].y += VG[ig].y * coef;
                V[ig].z -= VG[ig].z * coef;
            }
        }
    }
}


void P3Analysis::combineLocalVelocities(double alpha, double beta, std::vector<Vector3d> &VLocal) const
{
    VLocal.resize(nPanels()*3);
    double cosa = cos(alpha*PI/180.0);
    double sina = sin(alpha*PI/180.0);
    double cosb = cos(-beta*PI/180.0); //change of beta sign introduced in v7.24 to be consistent with AVL
    double sinb = sin(-beta*PI/180.0); //change of beta sign introduced in v7.24 to be consistent with AVL
    for(int i=0; i<int(VLocal.size()); i++)
    {
        VLocal[i].x = cosa*cosb* m_uVLocal.at(i).x + sinb*m_vVLocal.at(i).x + sina*cosb* m_wVLocal.at(i).x;
        VLocal[i].y = cosa*cosb* m_uVLocal.at(i).y + sinb*m_vVLocal.at(i).y + sina*cosb* m_wVLocal.at(i).y;
        VLocal[i].z = cosa*cosb* m_uVLocal.at(i).z + sinb*m_vVLocal.at(i).z + sina*cosb* m_wVLocal.at(i).z;
    }

//    for(int i=0; i<VLocal.size(); i++)        qDebug("  {:11g}  {:11g}  {:11g}", VLocal.at(i).x, VLocal.at(i).y, VLocal.at(i).z);
}


/**
 * The doublet solution (in the space of the right-hand sides: one value per panel for TRIUNIFORM, per vertex for TRILINEAR) of the steady state in which the air
 * flows past the plane's CoG at V (m/s) while the plane turns at Omega (rad/s), both in geometry axes: the onset at a point r is V + (r-CoG) x Omega,
 * the apparent air velocity of a plane that rotates (computeAngularDerivatives() builds its p, q and r perturbations from the same unit solutions).
 * The unit rotation solutions take their lever arm about the polar's CoG: a CoG elsewhere is a uniform velocity, which goes with V.
 */
void P3Analysis::makeStateSolution(Vector3d const &V, Vector3d const &Omega, Vector3d const &CoG, std::vector<double> &sol) const
{
    Vector3d Vc = V + (m_pPolar3d->CoG() - CoG) * Omega;
    sol.resize(m_uRHS.size());
    for(unsigned int i=0; i<sol.size(); i++)
    {
        sol[i]  = Vc.x    * m_uRHS[i] + Vc.y    * m_vRHS[i] + Vc.z    * m_wRHS[i];
        sol[i] += Omega.x * m_pRHS[i] + Omega.y * m_qRHS[i] + Omega.z * m_rRHS[i];
    }
}


/** The onset velocity at each panel's centroid, V + (r-CoG) x Omega, in m/s, geometry axes */
void P3Analysis::makeOnsetField(Vector3d const &V, Vector3d const &Omega, Vector3d const &CoG, std::vector<Vector3d> &field) const
{
    field.resize(nPanels());
    for(int i3=0; i3<nPanels(); i3++)
        field[i3] = V + (m_Panel3.at(i3).CoG() - CoG) * Omega;
}


/** Sets the doublet strengths of the analysis (3 per panel) from a solution in the space of the right-hand sides */
void P3Analysis::setStateDoublets(std::vector<double> const &sol)
{
    if(m_pPolar3d->isTriUniformMethod()) makeVertexDoubletDensities(sol, m_Mu);
    else                                 std::copy(sol.begin(), sol.end(), m_Mu.begin());
}


/**
 * The on-body Cp of a steady state (m_Cp), per panel the local onset flow's: Cp = 1-V^2/Vlocal^2. Rescaled to the free-stream dynamic pressure of QInf, so that the
 * panel forces come out in N/q like those of a point without rates (a thin panel's force, the jump in velocity times the onset, does not depend on the
 * reference of the pressure).
 */
void P3Analysis::makeStateCp(std::vector<double> const &sol, std::vector<Vector3d> const &field, double QInf, double alpha)
{
    std::vector<double> none(sol.size(), 0.0);
    makeLocalVelocities(sol, none, none, m_uVLocal, m_vVLocal, m_wVLocal, objects::windDirection(alpha, 0.0));
    computeOnBodyCp(field, m_uVLocal, m_Cp);
    for(int i3=0; i3<nPanels(); i3++)
    {
        double ratio2 = field.at(i3).dot(field.at(i3))/(QInf*QInf);
        for(int k=0; k<3; k++) m_Cp[3*i3+k] *= ratio2;
    }
}


/**
 * dsol = -A^-1 (dA_wake/dtheta) sol0: the first-order change of the doublet solution of the state sol0 when the wake, which leaves the trailing edges along the body x axis,
 * is yawed by theta (rad, positive towards +y) as it is along a sideslipped free stream: mu(theta) = mu_linear + theta dsol. This is the wake part of the
 * sideslip model of the derivatives (uav#79, computeTranslationDerivatives()). Returns false if the wake's matrix product could not be made.
 */
bool P3Analysis::makeYawedWakeCorrection(std::vector<double> const &sol0, std::vector<double> &dsol)
{
    if(!yawedWakeProduct(1.0e-2, sol0, dsol)) return false;
    for(double &d : dsol) d = -d;
    return backSubRHS(dsol);
}


/**
 * The far-field force (N) and the moment about the CoG (N.m) of a state, as the derivatives evaluate them: the moment of the far-field strip forces lumped at the strips'
 * quarter-chord points plus the strips' own chordwise loads, which is the pressure moment with the leading-edge force that the panels miss (forces(): Moment).
 * With bSideslipModel the state is evaluated by the sideslip model of the v derivatives (uav#79): the Trefftz plane of the far field is perpendicular to the onset
 * flow, whose direction is that of the first panel's, and the free-stream force of the doublet sheet is put on the sheet's vorticity instead of at the quarter points
 * (the chordwise loads stay). pOnBodyMoment, if given, receives the on-body moment of forces() (N.m: the pressures, with the thin strips' leading-edge force), which
 * is the moment of a point of a thick plane (PlaneTask::applySideslipModel()). It sets m_Cp.
 */
void P3Analysis::stateFarField(std::vector<double> const &sol, std::vector<Vector3d> const &field, double alpha, Vector3d const &CoG, bool bFuseMi, bool bSideslipModel,
                               Vector3d &Force, Vector3d &Moment, Vector3d *pOnBodyMoment)
{
    std::vector<double> none(sol.size(), 0.0);
    std::vector<double> mu3(3*nPanels());
    if(m_pPolar3d->isTriUniformMethod()) makeVertexDoubletDensities(sol, mu3);
    else                                 std::copy(sol.begin(), sol.end(), mu3.begin());
    std::vector<double> sigma(nPanels(), 0.0);   // thin surfaces: no sources

    makeLocalVelocities(sol, none, none, m_uVLocal, m_vVLocal, m_wVLocal, objects::windDirection(alpha, 0.0));
    computeOnBodyCp(field, m_uVLocal, m_Cp);

    bool bFarField = m_bFarField;
    m_bFarField = true;
    m_bSideslipFarField = bSideslipModel;
    m_bSideslipVorticity = bSideslipModel;
    Vector3d dummyForce, dummyMoment;
    forces(mu3.data(), sigma.data(), alpha, 0.0, CoG, bFuseMi, field, dummyForce, dummyMoment);
    m_bFarField = bFarField;
    m_bSideslipFarField = false;
    m_bSideslipVorticity = false;

    Force  = m_ForceFF;
    Moment = m_MomentFF + m_MomentLocal;
    if(pOnBodyMoment) *pOnBodyMoment = dummyMoment;
}
