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
    : add(_add), pAlg(_pAlg)
{
    removed.swap(_removed);
    vs.swap(_vs);
    oldAdj.resize(vs.size());
    pAlg->crt += add;
    for (int i = 0; i < static_cast<int>(removed.size()); i++)
        pAlg->vRestore[--(pAlg->rn)] = -1;
    for (int v : removed)
    {
        assert(pAlg->x[v] < 0);
        pAlg->x[v] = 2;
    }

    if (pAlg->bc_index_built)
    {
        std::vector<int>& mapping = pAlg->nodeMapping;

        // disable removed vertices in bc_index
        for (int v : removed)
        {
            for (int u : pAlg->adj[v])
            {
                if (pAlg->x[u] < 0)
                {
                    pAlg->bc_index->DeleteEdge(mapping[v], mapping[u]);
                    pAlg->bc_index->DeleteEdge(mapping[u], mapping[v]);

                    modStack.emplace_back(std::make_pair(mapping[v], mapping[u]), 1);
                    modStack.emplace_back(std::make_pair(mapping[u], mapping[v]), 1);
                }
            }
        }

        // remove deleted edges
        for (int v : vs)
        {
            for (int u : pAlg->adj[v])
            {
                if (pAlg->x[u] < 0)
                {
                    pAlg->bc_index->DeleteEdge(mapping[v], mapping[u]);
                    pAlg->bc_index->DeleteEdge(mapping[u], mapping[v]);

                    modStack.emplace_back(std::make_pair(mapping[v], mapping[u]), 1);
                    modStack.emplace_back(std::make_pair(mapping[u], mapping[v]), 1);
                }
            }
        }

        // add new edges
        for (int i = 0; i < static_cast<int>(vs.size()); i++)
        {
            for (int u : newAdj[i])
            {
                pAlg->bc_index->InsertEdge(mapping[vs[i]], mapping[u]);
                pAlg->bc_index->InsertEdge(mapping[u], mapping[vs[i]]);

                modStack.emplace_back(std::make_pair(mapping[vs[i]], mapping[u]), 2);
                modStack.emplace_back(std::make_pair(mapping[u], mapping[vs[i]]), 2);
            }
        }
    }

    for (int i = 0; i < static_cast<int>(vs.size()); i++)
    {
        oldAdj[i].swap(pAlg->adj[vs[i]]);
        pAlg->adj[vs[i]].swap(newAdj[i]);
    }
}

modified::modified(std::vector<int> &_removed, std::vector<int> &_vs, branch_and_reduce_algorithm *_pAlg)
    : add(0), pAlg(_pAlg)
{
    removed.swap(_removed);
    vs.swap(_vs);
}

void modified::restore()
{
    pAlg->crt -= add;
    pAlg->rn += removed.size();
    for (int v : removed)
        pAlg->x[v] = -1;

    if (pAlg->bc_index_built)
    {
        for (int i = 0; i < modStack.size(); i++)
        {
            std::pair<int, int> edge = modStack.back().first;
            int mode = modStack.back().second;
            modStack.pop_back();
             
            
            if (mode == 1)
                pAlg->bc_index->InsertEdge(edge.first, edge.second);
            else if (mode == 2)
                pAlg->bc_index->DeleteEdge(edge.first, edge.second);
        }
    }

    for (int i = 0; i < static_cast<int>(vs.size()); i++)
    {
        pAlg->adj[vs[i]] = oldAdj[i];
        int inV = pAlg->in[vs[i]], outV = pAlg->out[vs[i]];
        for (int u : pAlg->adj[vs[i]])
        {
            if (u == inV)
                inV = -1;
            if (u == outV)
                outV = -1;
        }
        if (inV >= 0)
        {
            pAlg->out[pAlg->in[vs[i]]] = -1;
            pAlg->in[vs[i]] = -1;
        }
        if (outV >= 0)
        {
            pAlg->in[pAlg->out[vs[i]]] = -1;
            pAlg->out[vs[i]] = -1;
        }
    }
}
