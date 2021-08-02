/******************************************************************************
 * Copyright (C) 2015-2017 Darren Strash <strash@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

// local includes
#include "modified.h"
#include "branch_and_reduce_algorithm.h"

// system includes
#include <vector>

modified::modified(int const _add, std::vector<int> &_removed, std::vector<int> &_vs, std::vector<std::vector<int>> &newAdj, branch_and_reduce_algorithm *_pAlg)
: add(_add)
, pAlg(_pAlg) 
{
    removed.swap(_removed);
    vs.swap(_vs);
    oldAdj.resize(vs.size());
    pAlg->crt += add;
    for (int i = 0; i < static_cast<int>(removed.size()); i++) pAlg->vRestore[--(pAlg->rn)] = -1;
    for (int v : removed) {
        assert(pAlg->x[v] < 0);
        pAlg->x[v] = 2;
    }
    for (int i = 0; i < static_cast<int>(vs.size()); i++) {
        oldAdj[i].swap(pAlg->adj[vs[i]]);
        pAlg->adj[vs[i]].swap(newAdj[i]);
    }    
}

modified::modified(std::vector<int> &_removed, std::vector<int> &_vs, branch_and_reduce_algorithm *_pAlg)
: add(0)
, pAlg(_pAlg)
{
    removed.swap(_removed);
    vs.swap(_vs);
}

void modified::restore() {
    pAlg->crt -= add;
    pAlg->rn += removed.size();
    for (int v : removed) pAlg->x[v] = -1;
    for (int i = 0; i < static_cast<int>(vs.size()); i++) {
        pAlg->adj[vs[i]] = oldAdj[i];
        int inV = pAlg->in[vs[i]], outV = pAlg->out[vs[i]];
        for (int u : pAlg->adj[vs[i]]) {
            if (u == inV) inV = -1;
            if (u == outV) outV = -1;
        }
        if (inV >= 0) {
            pAlg->out[pAlg->in[vs[i]]] = -1;
            pAlg->in[vs[i]] = -1;
        }
        if (outV >= 0) {
            pAlg->in[pAlg->out[vs[i]]] = -1;
            pAlg->out[vs[i]] = -1;
        }
    }
}

fold::fold(int const add, std::vector<int> &_removed, std::vector<int> &_vs, std::vector<std::vector<int>> &newAdj, branch_and_reduce_algorithm *_pAlg)
    : modified(add, _removed, _vs, newAdj, _pAlg)
    { 
        if (pAlg->USE_DEPENDENCY_CHECKING) {
            foldingMapAdd = 0;

            for (int v : removed) {
                if (pAlg->packingMap[v].size() != 0) {
                    pAlg->foldingMap[vs[0]].push_back(v);
                    foldingMapAdd += 1;
                }
            }        
        }
    }

void fold::restore() {
    modified::restore();

    if (pAlg->USE_DEPENDENCY_CHECKING) {
        for (int i = 0; i < foldingMapAdd; i++){
            pAlg->foldingMap[vs[0]].pop_back();
        }        
    }
}

alternative::alternative(int const add, std::vector<int> &_removed, std::vector<int> &_vs, std::vector<std::vector<int>> &newAdj, branch_and_reduce_algorithm *_pAlg, int k)
    : modified(add, _removed, _vs, newAdj, _pAlg)
    {
        this->k = k;

        if (pAlg->USE_DEPENDENCY_CHECKING) {
            foldingMapAdd = 0;

            for (int v : removed) {
                if (pAlg->packingMap[v].size() != 0) {
                    foldingMapAdd += 1;
                    for (int vv : vs) {
                        pAlg->foldingMap[vv].emplace_back(v);
                    }
                }
            }        
        }
    }

void alternative::restore() {
    modified::restore();

    if (pAlg->USE_DEPENDENCY_CHECKING) {
        for (int i = 0; i < foldingMapAdd; i++) {
            for (int vv : vs) {
                pAlg->foldingMap[vv].pop_back();
            }
        }  
    }
}