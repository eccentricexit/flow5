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


#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#include <objects_global.h>

struct EigenValues
{
    public:
        EigenValues()
        {
            reset();
        }

        void reset()
        {
            m_EV[0] = std::complex<double>(0.0,0.0);
            m_EV[1] = std::complex<double>(0.0,0.0);
            m_EV[2] = std::complex<double>(0.0,0.0);
            m_EV[3] = std::complex<double>(0.0,0.0);
            m_EV[4] = std::complex<double>(0.0,0.0);
            m_EV[5] = std::complex<double>(0.0,0.0);
            m_EV[6] = std::complex<double>(0.0,0.0);
            m_EV[7] = std::complex<double>(0.0,0.0);

            m_PhugoidDamping = 0.0;
            m_PhugoidFrequency = 0.0;
            m_RollDampingT2 = 0.0;
            m_ShortPeriodDamping = 0.0;
            m_ShortPeriodFrequency = 0.0;
            m_DutchRollDamping = 0.0;
            m_DutchRollFrequency = 0.0;
            m_SpiralDampingT2 = 0.0;

            for(int i=0; i<8; i++)
                for(int j=0; j<4; j++) m_Vec[i][j] = std::complex<double>(0.0,0.0);
            m_bVectors = false;
        }

        /**
         * Sets the figures of the modes from the eigenvalues. The modes are identified by their character rather than by their
         * index, since the eigenvalues are sorted on their real parts only, and not by the number of complex pairs: a small model's
         * short period is often overdamped (two real roots) while its phugoid stays an oscillation.
         * Longitudinal: the phugoid is the candidate (a complex pair as one unit, or two real roots) whose eigenvectors carry the
         * speed, the larger mean share |u|/(|u|+|w|) (state order u, w, q, theta), and the short period is the other two roots. Without
         * eigenvectors (setVectors), the phugoid is the candidate whose fastest root is the slower.
         * Lateral: the Dutch roll is the complex pair of the lateral roots; with two pairs, the one with the smaller roll share
         * |p|/(|p|+|r|) (state order v, p, r, phi), or the faster pair without eigenvectors. The roll mode is the most negative real
         * root and the spiral the real root closest to zero.
         * A mode that is not oscillatory (a short period or a phugoid split in two real roots, for instance) has zero frequency and
         * damping, and is not replaced by another root.
         */
        void computeModes()
        {
            double omegaN(0), omega1(0), zeta(0);
            double pi = 3.141592654;

            std::vector<std::complex<double>> pairs;
            std::vector<double> reals;
            std::vector<int> pairAt, realAt;    // where the roots are in m_EV

            m_ShortPeriodFrequency = m_ShortPeriodDamping = 0.0;
            m_PhugoidFrequency = m_PhugoidDamping = 0.0;
            splitRoots(m_EV, 4, pairs, reals, &pairAt, &realAt);
            if(2*pairs.size()+reals.size()==4)
            {
                // the units: a complex pair counts once, then the real roots
                int nPairs = int(pairs.size());
                std::vector<std::complex<double>> units(pairs);
                std::vector<int> at(pairAt);
                for(unsigned int i=0; i<reals.size(); i++)
                {
                    units.push_back(reals.at(i));
                    at.push_back(realAt.at(i));
                }
                // the candidates for the phugoid: a pair, or two real roots
                std::vector<std::vector<int>> cands;
                for(int i=0; i<nPairs; i++) cands.push_back({i});
                for(int i=nPairs; i<int(units.size()); i++)
                    for(int j=i+1; j<int(units.size()); j++) cands.push_back({i, j});

                std::vector<double> share;
                bool bShare = m_bVectors;
                for(unsigned int i=0; i<units.size() && bShare; i++)
                {
                    share.push_back(vectorShare(at.at(i), 0, 1));
                    if(share.back()<0.0) bShare = false;
                }

                int best = 0;
                double bestScore(0);
                for(unsigned int ic=0; ic<cands.size(); ic++)
                {
                    double score(0), fastest(0);
                    for(int i : cands.at(ic))
                    {
                        if(bShare) score += share.at(i)/double(cands.at(ic).size());
                        else       fastest = std::max(fastest, std::abs(units.at(i)));
                    }
                    if(!bShare) score = -fastest;
                    if(ic==0 || score>bestScore)
                    {
                        best = ic;
                        bestScore = score;
                    }
                }

                for(int i=0; i<nPairs; i++) // a real root has no frequency
                {
                    bool bPhugoid = std::find(cands.at(best).begin(), cands.at(best).end(), i) != cands.at(best).end();
                    objects::modeProperties(units.at(i), omegaN, omega1, zeta);
                    if(bPhugoid)
                    {
                        m_PhugoidFrequency = omegaN/2.0/pi;
                        m_PhugoidDamping   = zeta;
                    }
                    else
                    {
                        m_ShortPeriodFrequency = omegaN/2.0/pi;
                        m_ShortPeriodDamping   = zeta;
                    }
                }
            }

            splitRoots(m_EV+4, 4, pairs, reals, &pairAt);
            m_DutchRollFrequency = m_DutchRollDamping = 0.0;
            m_RollDampingT2 = m_SpiralDampingT2 = 0.0;
            if(pairs.size()>0)
            {
                int iDutch = 0;   // the faster pair
                if(pairs.size()>1 && m_bVectors)
                {
                    double smallest(2.0);
                    for(unsigned int i=0; i<pairs.size(); i++)
                    {
                        double s = vectorShare(4+pairAt.at(i), 1, 2);
                        if(s<0.0) {iDutch=0; break;}
                        if(s<smallest) {smallest=s; iDutch=i;}
                    }
                }
                objects::modeProperties(pairs.at(iDutch), omegaN, omega1, zeta);
                m_DutchRollFrequency = omegaN/2.0/pi;
                m_DutchRollDamping   = zeta;
            }
            std::sort(reals.begin(), reals.end());
            if(reals.size()>0) m_RollDampingT2 = log(2.0)/fabs(reals.front());
            if(reals.size()>1)
            {
                double spiral = reals.at(1);
                for(unsigned int i=2; i<reals.size(); i++) if(fabs(reals.at(i))<fabs(spiral)) spiral = reals.at(i);
                m_SpiralDampingT2 = log(2.0)/fabs(spiral);
            }
        }


        /** Gives the eigenvectors of the eight eigenvalues (the longitudinal then the lateral, 4 components each) for computeModes to name the modes by */
        void setVectors(std::complex<double> const v[8][4])
        {
            for(int i=0; i<8; i++)
                for(int j=0; j<4; j++) m_Vec[i][j] = v[i][j];
            m_bVectors = true;
        }


        /** Splits n eigenvalues into the complex conjugate pairs, one representative with a positive imaginary part for each and sorted by decreasing frequency, and the real roots.
         *  The optional pairAt and realAt receive the index in ev of a root of each pair and of each real root, in the same order. */
        static void splitRoots(std::complex<double> const *ev, int n, std::vector<std::complex<double>> &pairs, std::vector<double> &reals,
                               std::vector<int> *pairAt=nullptr, std::vector<int> *realAt=nullptr)
        {
            pairs.clear();
            reals.clear();
            std::vector<int> pairIdx, realIdx;
            std::vector<bool> used(n, false);
            for(int i=0; i<n; i++)
            {
                if(used.at(i)) continue;
                used[i] = true;
                if(fabs(ev[i].imag())>1.e-15)
                {
                    for(int j=i+1; j<n; j++)
                    {
                        if(!used.at(j) && std::abs(std::conj(ev[j])-ev[i])<1.e-6*std::max(1.0, std::abs(ev[i])))
                        {
                            used[j] = true;
                            break;
                        }
                    }
                    pairs.push_back({ev[i].real(), fabs(ev[i].imag())});
                    pairIdx.push_back(i);
                }
                else
                {
                    reals.push_back(ev[i].real());
                    realIdx.push_back(i);
                }
            }
            std::vector<int> order(pairs.size());
            for(unsigned int i=0; i<order.size(); i++) order[i] = i;
            std::stable_sort(order.begin(), order.end(), [&pairs](int a, int b){return pairs.at(a).imag()>pairs.at(b).imag();});
            std::vector<std::complex<double>> sorted;
            if(pairAt) pairAt->clear();
            for(int i : order)
            {
                sorted.push_back(pairs.at(i));
                if(pairAt) pairAt->push_back(pairIdx.at(i));
            }
            pairs = sorted;
            if(realAt) *realAt = realIdx;
        }


    private:
        /** The share |v[num]| / (|v[num]| + |v[den]|) of the eigenvector of root i, or -1 if both components are null (no eigenvector) */
        double vectorShare(int i, int num, int den) const
        {
            double a = std::abs(m_Vec[i][num]);
            double d = std::abs(m_Vec[i][den]);
            if(!(a+d>0.0)) return -1.0;
            return a/(a+d);
        }


    public:
        double m_PhugoidFrequency;        /**< the phugoid's frequency, as a result of stability analysis only */
        double m_PhugoidDamping;          /**< the phugoid's damping factor, as a result of stability analysis only */
        double m_RollDampingT2;           /**< the time to double or half for the damping of the roll-damping mode, as a result of stability analysis only */
        double m_ShortPeriodDamping;      /**< the damping of the short period mode, as a result of stability analysis only */
        double m_ShortPeriodFrequency;    /**< the frequency of the short period mode, as a result of stability analysis only */
        double m_DutchRollDamping;        /**< the damping of the Dutch roll mode, as a result of stability analysis only */
        double m_DutchRollFrequency;      /**< the frequency of the Dutch roll mode, as a result of stability analysis only */
        double m_SpiralDampingT2;         /**< the time to double or half for the damping of the spiral mode, as a result of stability analysis only >*/

        std::complex<double> m_EV[8];

        std::complex<double> m_Vec[8][4];  /**< the eigenvectors of the eight eigenvalues, if known; not saved with the polar, so modes read from a file are named by time scale */
        bool m_bVectors;                   /**< true if m_Vec holds the eigenvectors */

};
