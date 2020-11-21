/******************************************************************************
 * Copyright (C) 2015-2017 Darren Strash <strash@kit.edu>
 * Copyright (C) 2019 Demian Hespe <hespe@kit.edu>
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

#include "branch_and_reduce_algorithm.h"
#include "fast_set.h"
#include "modified.h"

#include <stack>
#include <vector>
#include <set>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <algorithm> // sort()
#include <deque>
#include <chrono>
#include <list>
#include <set>

#include <cstring>
#include "../Metis/include/metis.h"

using namespace std;

int branch_and_reduce_algorithm::REDUCTION = 3;
int branch_and_reduce_algorithm::LOWER_BOUND = 4;
int branch_and_reduce_algorithm::BRANCHING = 2;
bool branch_and_reduce_algorithm::outputLP = false;
long branch_and_reduce_algorithm::nBranchings = 0;
int branch_and_reduce_algorithm::debug = 0;

int branch_and_reduce_algorithm::EXTRA_DECOMP = 0;
int branch_and_reduce_algorithm::ND_LEVEL = 32;

long branch_and_reduce_algorithm::TUNING_PARAM1 = 7;
double branch_and_reduce_algorithm::TUNING_PARAM2 = 0.1;
long branch_and_reduce_algorithm::TUNING_PARAM3 = 10;

long branch_and_reduce_algorithm::defaultBranchings = 0;
bool branch_and_reduce_algorithm::defaultBranch = false;
long branch_and_reduce_algorithm::defaultPicks = 0;
long branch_and_reduce_algorithm::stratPicks = 0;
long branch_and_reduce_algorithm::nDecomps = 0;
long branch_and_reduce_algorithm::prunes = 0;

long branch_and_reduce_algorithm::ro = 1;

inline int getReverseIndex(std::vector<std::vector<int>> &adj, int source, int target)
{
    for (int i = 0; i < adj[source].size(); i++)
        if (adj[source][i] == target)
            return i;
}

branch_and_reduce_algorithm::branch_and_reduce_algorithm(std::vector<std::vector<int>> &_adj, int const _N)
    : adj(), n(_adj.size()), used(n * 2), domination_graph(n), chainLength(n), chains(n), chain_vec(n)
{
    SHRINK = 0.5;
    depth = 0;
    maxDepth = 10;
    rootDepth = -1; // invalid value

    n = _adj.size();
    adj.swap(_adj);

    N = _N;
    opt = n;
    y.resize(N, 0);
    for (int i = 0; i < n; i++)
        y[i] = 1;
    for (int i = n; i < N; i++)
        y[i] = 2;
    crt = 0;
    x.resize(N, 0);
    for (int i = 0; i < n; i++)
        x[i] = -1;
    for (int i = n; i < N; i++)
        x[i] = 2;
    rn = n;
    in.resize(n, -1);
    out.resize(n, -1);
    lb = -1; // invalid value

    vRestore.resize(n, 0);

    que.resize(n * 2, 0);
    level.resize(n * 2, 0);
    iter.resize(n * 2, 0);

    modTmp.resize(n, 0);

    modifiedN = 0;
    modifieds.resize(N, shared_ptr<modified>());

    // MODIFICATIONS
    s = t = -1;
}

int branch_and_reduce_algorithm::deg(int v)
{
    assert(x[v] < 0);
    int deg = 0;
    for (int u : adj[v])
        if (x[u] < 0)
            deg++;

    return deg;
}

void branch_and_reduce_algorithm::set(int v, int a)
{
    assert(x[v] < 0);
    crt += a;
    x[v] = a;
    vRestore[--rn] = v;

    if (a == 0)
    {
        for (int u : adj[v])
            if (x[u] < 0)
            {
                x[u] = 1;
                crt++;
                vRestore[--rn] = u;
            }
    }
}

// methods that modify the graph

void branch_and_reduce_algorithm::compute_fold(std::vector<int> const &S, std::vector<int> const &NS)
{
    assert(NS.size() == S.size() + 1);

    // remove all vertices but the 1. neighbour
    std::vector<int> removed(S.size() * 2);
    for (unsigned int i = 0; i < S.size(); i++)
        removed[i] = S[i];
    for (unsigned int i = 0; i < S.size(); i++)
        removed[S.size() + i] = NS[1 + i];

    // s = first neighbour
    int s = NS[0];
    used.clear();
    for (int v : S)
        used.add(v);
    std::vector<int> &tmp = modTmp;
    int p = 0;
    for (int v : NS)
    {
        assert(!used.get(v));
        for (int u : adj[v])
            if (x[u] < 0 && used.add(u))
            {
                tmp[p++] = u;
            }
    }
    // tmp contains N(N(S))\S

    // set edges of s
    std::vector<std::vector<int>> newAdj(p + 1);
    {
        std::vector<int> copyOfTmp(tmp.begin(), tmp.begin() + p);
        newAdj[0].swap(copyOfTmp);
    }
    std::sort(newAdj[0].begin(), newAdj[0].end());
    std::vector<int> vs(p + 1);
    vs[0] = s;
    used.clear();
    for (int v : S)
        used.add(v);
    for (int v : NS)
        used.add(v);

    // set edges of vertices in tmp
    for (unsigned int i = 0; i < newAdj[0].size(); i++)
    {
        int v = newAdj[0][i];
        p = 0;
        bool add = false;
        for (int u : adj[v])
            if (x[u] < 0 && !used.get(u))
            {
                if (!add && s < u)
                {
                    tmp[p++] = s;
                    add = true;
                }
                tmp[p++] = u;
            }
        if (!add)
            tmp[p++] = s;
        vs[1 + i] = v;

        {
            std::vector<int> copyOfTmp(tmp.begin(), tmp.begin() + p);
            newAdj[i + 1].swap(copyOfTmp);
        }
    }

    // test
    if (ro == 1)
    {
        for (int i = 0; i < nd_order.size(); i++)
        {
            if (std::find(removed.begin(), removed.end(), nd_order[i]) != removed.end())
                nd_order[i] = s;
        }
    }

    modifieds[modifiedN++] = make_shared<fold>(fold(S.size(), removed, vs, newAdj, this));
}

void branch_and_reduce_algorithm::compute_alternative(std::vector<int> const &A, std::vector<int> const &B)
{
    assert(A.size() == B.size());
    used.clear();
    for (int b : B)
        for (int u : adj[b])
            if (x[u] < 0)
                used.add(u);
    for (int a : A)
        for (int u : adj[a])
            if (x[u] < 0 && used.get(u))
                set(u, 1);
    NodeID p = 0, q = 0;
    std::vector<int> &tmp = modTmp;
    used.clear();
    for (int b : B)
        used.add(b);
    for (int a : A)
        for (int u : adj[a])
            if (x[u] < 0 && used.add(u))
                tmp[p++] = u;
    std::vector<int> A2(tmp.begin(), tmp.begin() + p);
    std::sort(A2.begin(), A2.end());
    p = 0;
    used.clear();
    for (int a : A)
        used.add(a);
    for (int b : B)
        for (int u : adj[b])
            if (x[u] < 0 && used.add(u))
                tmp[p++] = u;
    std::vector<int> B2(tmp.begin(), tmp.begin() + p);
    std::sort(B2.begin(), B2.end());
    std::vector<int> removed(A.size() + B.size());
    for (unsigned int i = 0; i < A.size(); i++)
        removed[i] = A[i];
    for (unsigned int i = 0; i < B.size(); i++)
        removed[A.size() + i] = B[i];
    std::vector<int> vs(A2.size() + B2.size());
    for (unsigned int i = 0; i < A2.size(); i++)
        vs[i] = A2[i];
    for (unsigned int i = 0; i < B2.size(); i++)
        vs[A2.size() + i] = B2[i];
    std::vector<std::vector<int>> newAdj(vs.size());
    used.clear();
    for (int a : A)
        used.add(a);
    for (int b : B)
        used.add(b);
    for (unsigned int i = 0; i < vs.size(); i++)
    {
        unsigned int v = (i < A2.size()) ? A2[i] : B2[i - A2.size()];
        std::vector<int> const &C = (i < A2.size()) ? B2 : A2;
        p = q = 0;
        for (int u : adj[v])
            if (x[u] < 0 && !used.get(u))
            {
                while (q < C.size() && C[q] <= u)
                {
                    if (used.get(C[q]))
                        q++;
                    else
                        tmp[p++] = C[q++];
                }
                if (p == 0 || tmp[p - 1] != u)
                    tmp[p++] = u;
            }
        while (q < C.size())
        {
            if (used.get(C[q]))
                q++;
            else
                tmp[p++] = C[q++];
        }
        {
            std::vector<int> copyOfTmp(tmp.begin(), tmp.begin() + p);
            newAdj[i].swap(copyOfTmp);
        }
    }

    modifieds[modifiedN++] = make_shared<alternative>(alternative(removed.size() / 2, removed, vs, newAdj, this, A2.size()));
}

void branch_and_reduce_algorithm::restore(int n)
{
    while (rn < n)
    {
        int v = vRestore[rn];
        if (v >= 0)
        {
            crt -= x[v];
            x[v] = -1;
            rn++;
        }
        else
        {
            modifieds[--modifiedN]->restore();
            modifieds[modifiedN] = shared_ptr<modified>();
        }
    }
}

void branch_and_reduce_algorithm::reverse()
{
    for (int i = modifiedN - 1; i >= 0; i--)
    {
        modifieds[i]->reverse(y);
    }
}

// lower bounds

int branch_and_reduce_algorithm::lpLowerBound()
{
    return crt + (rn + 1) / 2;
}

int branch_and_reduce_algorithm::cycleLowerBound()
{
    int lb = crt;
    std::vector<int> &id = iter;
    for (int i = 0; i < n; i++)
        id[i] = -1;
    std::vector<int> &pos = que;
    std::vector<int> &S = level;
    std::vector<int> &S2 = modTmp;
    for (int i = 0; i < n; i++)
        if (x[i] < 0 && id[i] < 0)
        {
            int v = i;
            int size = 0;
            do
            {
                assert(id[v] < 0);
                id[v] = i;
                v = out[v];
                pos[v] = size;
                S[size++] = v;
            } while (v != i);
            bool clique = true;
            for (int j = 0; j < size; j++)
            {
                v = S[j];
                int num = 0;
                for (int u : adj[v])
                    if (x[u] < 0 && id[u] == id[v])
                        num++;
                if (num != size - 1)
                {
                    clique = false;
                    break;
                }
            }
            if (clique)
            {
                lb += size - 1;
            }
            else
            {
                while (size >= 6)
                {
                    int minSize = size, s = 0, t = size;
                    for (int j = 0; j < size; j++)
                    {
                        used.clear();
                        v = S[j];
                        for (int u : adj[v])
                            if (x[u] < 0 && id[u] == id[v])
                            {
                                used.add(u);
                            }
                        v = S[(j + 1) % size];
                        for (int u : adj[v])
                            if (x[u] < 0 && id[u] == id[v])
                            {
                                if (used.get(S[(pos[u] + 1) % size]))
                                {
                                    int size2 = (pos[u] - j + size) % size;
                                    if (minSize > size2 && size2 % 2 != 0)
                                    {
                                        minSize = size2;
                                        s = (j + 1) % size;
                                        t = (pos[u] + 1) % size;
                                    }
                                }
                            }
                    }
                    if (minSize == size)
                        break;
                    int p = 0;
                    for (int j = t; j != s; j = (j + 1) % size)
                    {
                        S2[p++] = S[j];
                    }
                    for (int j = s; j != t; j = (j + 1) % size)
                    {
                        id[S[j]] = n;
                    }

                    S.swap(S2);

                    size -= minSize;
                    assert(size == p);
                    assert(minSize > 1);
                    lb += (minSize + 1) / 2;
                    for (int j = 0; j < size; j++)
                        pos[S[j]] = j;
                }
                assert(size > 1);
                lb += (size + 1) / 2;
            }
        }

    if (static_cast<int>(level.size()) != n * 2)
    {
        level.swap(modTmp);
    }
    return lb;
}

int branch_and_reduce_algorithm::cliqueLowerBound()
{
    int need = crt;
    std::vector<long long> ls(rn, 0);
    int k = 0;
    for (int i = 0; i < n; i++)
        if (x[i] < 0)
            ls[k++] = ((long long)deg(i)) << 32 | i;
    std::sort(ls.begin(), ls.end());
    std::vector<int> &clique = que;
    std::vector<int> &size = level;
    std::vector<int> &tmp = iter;
    used.clear();
    for (int i = 0; i < rn; i++)
    {
        int v = (int)ls[i];
        int to = v, max = 0;
        for (int u : adj[v])
            if (x[u] < 0 && used.get(u))
                tmp[clique[u]] = 0;
        for (int u : adj[v])
            if (x[u] < 0 && used.get(u))
            {
                int c = clique[u];
                tmp[c]++;
                if (tmp[c] == size[c] && max < size[c])
                {
                    to = c;
                    max = size[c];
                }
            }
        clique[v] = to;
        if (to != v)
        {
            size[to]++;
            need++;
        }
        else
        {
            size[v] = 1;
        }
        used.add(v);
    }
    return need;
}

int branch_and_reduce_algorithm::lowerBound()
{
    int type = 0, tmp;
    if (lb < crt)
    {
        lb = crt;
        type = 1;
    }
    if (LOWER_BOUND == 1 || LOWER_BOUND == 4)
    {
        tmp = cliqueLowerBound();
        if (lb < tmp)
        {
            lb = tmp;
            type = 4;
        }
    }
    if (LOWER_BOUND == 2 || LOWER_BOUND == 4)
    {
        tmp = lpLowerBound();
        if (lb < tmp)
        {
            lb = tmp;
            type = 2;
        }
    }
    if (LOWER_BOUND == 3 || LOWER_BOUND == 4)
    {
        tmp = cycleLowerBound();
        if (lb < tmp)
        {
            lb = tmp;
            type = 3;
        }
    }
    if (debug >= 2 && depth <= maxDepth)
        fprintf(stderr, "%slb: %d (%d), %d\n", debugString().c_str(), lb, type, opt);
    return lb;
}

// helper for lpReduction
bool branch_and_reduce_algorithm::dinicDFS(int v)
{
    while (iter[v] >= 0)
    {
        int u = adj[v][iter[v]--], w = in[u];
        if (x[u] >= 0)
            continue;
        if (w < 0 || (level[v] < level[w] && iter[w] >= 0 && dinicDFS(w)))
        {
            in[u] = v;
            out[v] = u;
            return true;
        }
    }
    return false;
}

// helper for lpReduction
void branch_and_reduce_algorithm::updateLP()
{
#if 1
    for (int v = 0; v < n; v++)
        if (out[v] >= 0 && ((x[v] < 0) ^ (x[out[v]] < 0)))
        {
            in[out[v]] = -1;
            out[v] = -1;
        }
    for (;;)
    {
        used.clear();
        int qs = 0, qt = 0;
        for (int v = 0; v < n; v++)
            if (x[v] < 0 && out[v] < 0)
            {
                level[v] = 0;
                used.add(v);
                que[qt++] = v;
            }
        bool ok = false;
        while (qs < qt)
        {
            int v = que[qs++];
            iter[v] = adj[v].size() - 1;
            for (int u : adj[v])
                if (x[u] < 0 && used.add(n + u))
                {
                    int w = in[u];
                    if (w < 0)
                        ok = true;
                    else
                    {
                        level[w] = level[v] + 1;
                        used.add(w);
                        que[qt++] = w;
                    }
                }
        }
        if (!ok)
            break;
        for (int v = n - 1; v >= 0; v--)
            if (x[v] < 0 && out[v] < 0)
            {
                dinicDFS(v);
            }
    }

#endif // 0
}

// reductions

bool branch_and_reduce_algorithm::lpReduction()
{
    int oldn = rn;
    updateLP();

#if 1
    for (int v = 0; v < n; v++)
    {
        if (x[v] < 0 && used.get(v) && !used.get(n + v))
            set(v, 0);
    }
    used.clear();
    int p = 0;
    iter.assign(iter.size(), 0);
    for (int s = 0; s < n; s++)
        if (x[s] < 0 && used.add(s))
        {
            int qt = 0;
            que[qt] = s;
            while (qt >= 0)
            {
                int v = que[qt], u = -1;
                if (v < n)
                {
                    while (iter[v] < static_cast<int>(adj[v].size()))
                    {
                        u = n + adj[v][iter[v]++];
                        if (x[u - n] < 0 && used.add(u))
                        {
                            break;
                        }
                        u = -1;
                    }
                }
                else if (used.add(in[v - n]))
                {
                    u = in[v - n];
                }
                if (u >= 0)
                {
                    que[++qt] = u;
                }
                else
                {
                    level[p++] = v;
                    qt--;
                }
            }
        }
    used.clear();
    for (int i = p - 1; i >= 0; i--)
        if (used.add(level[i]))
        {
            int v = level[i];
            int qs = 0, qt = 0;
            que[qt++] = v;
            bool ok = true;
            while (qs < qt)
            {
                v = que[qs++];
                if (used.get(v >= n ? (v - n) : (v + n)))
                    ok = false;
                if (v >= n)
                {
                    for (int u : adj[v - n])
                        if (x[u] < 0 && used.add(u))
                        {
                            que[qt++] = u;
                        }
                }
                else if (used.add(n + out[v]))
                {
                    que[qt++] = n + out[v];
                }
            }
            ok = false;
            if (ok)
            {
                for (int j = 0; j < qt; j++)
                {
                    v = que[j];
                    if (v >= n)
                        set(v - n, 0);
                }
            }
        }

#endif // 0
    if (debug >= 3 && depth <= maxDepth && oldn != rn)
        fprintf(stderr, "%sLP: %d -> %d\n", debugString().c_str(), oldn, rn);
    return oldn != rn;
}

bool branch_and_reduce_algorithm::deg1Reduction()
{
    int oldn = rn;
    std::vector<int> &deg = iter;
    int qt = 0;
    used.clear();
    for (int v = 0; v < n; v++)
        if (x[v] < 0)
        {
            deg[v] = n == rn ? adj[v].size() : this->deg(v);
            if (deg[v] <= 1)
            {
                que[qt++] = v;
                used.add(v);
            }
        }
    while (qt > 0)
    {
        int v = que[--qt];
        if (x[v] >= 0)
            continue;
        assert(deg[v] <= 1);
        for (int u : adj[v])
            if (x[u] < 0)
            {
                for (int w : adj[u])
                    if (x[w] < 0)
                    {
                        deg[w]--;
                        if (deg[w] <= 1 && used.add(w))
                            que[qt++] = w;
                    }
            }
        set(v, 0);
    }
    if (debug >= 3 && depth <= maxDepth && oldn != rn)
        fprintf(stderr, "%sdeg1: %d -> %d\n", debugString().c_str(), oldn, rn);
    return oldn != rn;
}

bool branch_and_reduce_algorithm::dominateReduction()
{
    int oldn = rn;
#if 1
    for (int v = 0; v < n; v++)
        if (x[v] < 0)
        {
            used.clear();
            used.add(v);
            for (int u : adj[v])
                if (x[u] < 0)
                    used.add(u);
            for (int u : adj[v])
                if (x[u] < 0)
                {
                    int cnt = 0;
                    int vtx = -1;
                    for (int w : adj[u])
                    {
                        if (x[w] < 0 && !used.get(w))
                        {
                            cnt++;
                            vtx = w;
                        }
                        if (cnt >= 2)
                            goto loop;
                    }
                    if (cnt == 1)
                    {
                        domin_vtcs.push_back(vtx);
                    }
                    else
                    {
                        set(v, 1);
                        break;
                    }
                loop:;
                }
        }
    if (debug >= 3 && depth <= maxDepth && oldn != rn)
        fprintf(stderr, "%sdominate: %d -> %d\n", debugString().c_str(), oldn, rn);
#endif // 0
    return oldn != rn;
}

bool branch_and_reduce_algorithm::almost_dominated()
{
    bool found = false;
    for (int v = 0; v < n; v++)
        if (x[v] < 0)
        {
            used.clear();
            used.add(v);
            for (int u : adj[v])
                if (x[u] < 0)
                    used.add(u);
            for (int u : adj[v])
                if (x[u] < 0)
                {
                    int cnt = 0;
                    int vtx = -1;
                    for (int w : adj[u])
                    {
                        if (x[w] < 0 && !used.get(w))
                        {
                            cnt++;
                            vtx = w;
                        }
                        if (cnt >= 2)
                            goto loop;
                    }
                    found = true;
                    domin_vtcs.push_back(vtx);
                loop:;
                }
        }
    return found;
}

bool branch_and_reduce_algorithm::fold2Reduction()
{
    int oldn = rn;
    std::vector<int> &tmp = level;
    for (int v = 0; v < n; v++)
        if (x[v] < 0)
        {
            int p = 0;
            for (int u : adj[v])
                if (x[u] < 0)
                {
                    tmp[p++] = u;
                    if (p > 2)
                        goto loop;
                }
            if (p < 2)
                continue;
            for (int u : adj[tmp[0]])
                if (u == tmp[1]) // triangle
                {
                    set(v, 0);
                    goto loop;
                }
            {
                std::vector<int> copyOfTmp(tmp.begin(), tmp.begin() + 2);
                compute_fold(std::vector<int>{v}, copyOfTmp);
            }
        loop:;
        }
    if (debug >= 3 && depth <= maxDepth && oldn != rn)
        fprintf(stderr, "%sfold2: %d -> %d\n", debugString().c_str(), oldn, rn);
    return oldn != rn;
}

bool branch_and_reduce_algorithm::twinReduction()
{
    int oldn = rn;
    std::vector<int> &vUsed = iter;
    int uid = 0;
    std::vector<int> NS(3, 0);
    for (int i = 0; i < n; i++)
        vUsed[i] = 0;
    for (int v = 0; v < n; v++)
        if (x[v] < 0 && deg(v) == 3)
        {
            int p = 0;
            for (int u : adj[v])
                if (x[u] < 0)
                {
                    NS[p++] = u;
                    uid++;
                    for (int w : adj[u])
                        if (x[w] < 0 && w != v)
                        {
                            if (p == 1)
                                vUsed[w] = uid;
                            else if (vUsed[w] == uid - 1)
                            {
                                vUsed[w]++;
                                if (p == 3 && deg(w) == 3)
                                {
                                    uid++;
                                    for (int z : NS)
                                        vUsed[z] = uid;
                                    bool ind = true;
                                    for (int z : NS)
                                        for (int a : adj[z])
                                            if (x[a] < 0 && vUsed[a] == uid)
                                                ind = false;
                                    if (ind)
                                    {
                                        compute_fold(std::vector<int>{v, w}, NS);
                                    }
                                    else
                                    {
                                        set(v, 0);
                                        set(w, 0);
                                    }
                                    //cout << "twin: " << nBranchings << endl;
                                    goto loop;
                                }
                                else if (BRANCHING == 15 && p == 3 && deg(w) == 4)
                                {
                                    for (int z : adj[w])
                                    {
                                        if (x[z] < 0 && vUsed[z] == 0)
                                            twin_vtcs.push_back(z);
                                    }
                                }
                            }
                        }
                }
        loop:;
        }
    if (debug >= 3 && depth <= maxDepth && oldn != rn)
        fprintf(stderr, "%stwin: %d -> %d\n", debugString().c_str(), oldn, rn);
    return oldn != rn;
}

bool branch_and_reduce_algorithm::funnelReduction()
{
    if (BRANCHING == 16)
        return funnelReduction_a();

    int oldn = rn;
    for (int v = 0; v < n; v++)
        if (x[v] < 0)
        {
            used.clear();
            std::vector<int> &tmp = level;
            int p = 0;
            for (int u : adj[v])
                if (x[u] < 0 && used.add(u))
                {
                    tmp[p++] = u;
                }
            if (p <= 1)
            {
                set(v, 0);
                continue;
            }
            int u1 = -1;
            for (int i = 0; i < p; i++)
            {
                int d = 0;
                for (int u : adj[tmp[i]])
                    if (x[u] < 0 && used.get(u))
                        d++;
                if (d + 1 < p)
                {
                    u1 = tmp[i];
                    break;
                }
            }
            if (u1 < 0)
            {
                set(v, 0);
                continue;
            }
            else
            {
                std::vector<int> &id = iter;
                for (int i = 0; i < p; i++)
                    id[tmp[i]] = -1;
                for (int u : adj[u1])
                    if (x[u] < 0)
                        id[u] = 0;
                int u2 = -1;
                for (int i = 0; i < p; i++)
                    if (tmp[i] != u1 && id[tmp[i]] < 0)
                    {
                        u2 = tmp[i];
                        break;
                    }
                assert(u2 >= 0);
                used.remove(u1);
                used.remove(u2);
                int d1 = 0, d2 = 0;
                for (int w : adj[u1])
                    if (x[w] < 0 && used.get(w))
                        d1++;
                for (int w : adj[u2])
                    if (x[w] < 0 && used.get(w))
                        d2++;
                if (d1 < p - 2 && d2 < p - 2)
                    continue;
                for (int i = 0; i < p; i++)
                {
                    int u = tmp[i];
                    if (u == u1 || u == u2)
                        continue;
                    int d = 0;
                    for (int w : adj[u])
                        if (x[w] < 0 && used.get(w))
                            d++;
                    if (d < p - 3)
                    {
                        goto loop;
                    }
                }
                int u = (d1 == p - 2) ? u2 : u1;
                std::vector<int> const v1{v};
                std::vector<int> const v2{u};
                compute_alternative(v1, v2);
            }
        loop:;
        }
    if (debug >= 3 && depth <= maxDepth && oldn != rn)
        fprintf(stderr, "%sfunnel: %d -> %d\n", debugString().c_str(), oldn, rn);
    return oldn != rn;
}

bool branch_and_reduce_algorithm::checkFunnel(int v)
{
    used.clear();
    std::vector<int> &tmp = level;
    int p = 0;
    for (int u : adj[v])
        if (x[u] < 0 && used.add(u))
        {
            tmp[p++] = u;
        }
    if (p <= 1)
        return false;
    int u1 = -1;
    for (int i = 0; i < p; i++)
    {
        int d = 0;
        for (int u : adj[tmp[i]])
            if (x[u] < 0 && used.get(u))
                d++;
        if (d + 1 < p)
        {
            u1 = tmp[i];
            break;
        }
    }
    if (u1 < 0)
    {
        return false;
    }
    else
    {
        std::vector<int> &id = iter;
        for (int i = 0; i < p; i++)
            id[tmp[i]] = -1;
        for (int u : adj[u1])
            if (x[u] < 0)
                id[u] = 0;
        int u2 = -1;
        for (int i = 0; i < p; i++)
            if (tmp[i] != u1 && id[tmp[i]] < 0)
            {
                u2 = tmp[i];
                break;
            }
        assert(u2 >= 0);
        used.remove(u1);
        used.remove(u2);
        int d1 = 0, d2 = 0;
        for (int w : adj[u1])
            if (x[w] < 0 && used.get(w))
                d1++;
        for (int w : adj[u2])
            if (x[w] < 0 && used.get(w))
                d2++;
        if (d1 < p - 2 && d2 < p - 2)
            return false;
        for (int i = 0; i < p; i++)
        {
            int u = tmp[i];
            if (u == u1 || u == u2)
                continue;
            int d = 0;
            for (int w : adj[u])
                if (x[w] < 0 && used.get(w))
                    d++;
            if (d < p - 3)
                return false;
        }
        return true;
    }
}

bool branch_and_reduce_algorithm::funnelReduction_a()
{
    int oldn = rn;
    for (int v = 0; v < n; v++)
        if (x[v] < 0)
        {
            used.clear();
            std::vector<int> &tmp = level;
            int p = 0;
            for (int u : adj[v])
                if (x[u] < 0 && used.add(u))
                {
                    tmp[p++] = u;
                }
            if (p <= 1)
            {
                set(v, 0);
                continue;
            }
            int u1 = -1;
            for (int i = 0; i < p; i++)
            {
                int d = 0;
                for (int u : adj[tmp[i]])
                    if (x[u] < 0 && used.get(u))
                        d++;
                if (d + 1 < p)
                {
                    u1 = tmp[i];
                    break;
                }
            }
            if (u1 < 0)
            {
                set(v, 0);
                continue;
            }
            else
            {
                std::vector<int> &id = iter;
                for (int i = 0; i < p; i++)
                    id[tmp[i]] = -1;
                for (int u : adj[u1])
                    if (x[u] < 0)
                        id[u] = 0;
                int u2 = -1;
                for (int i = 0; i < p; i++)
                    if (tmp[i] != u1 && id[tmp[i]] < 0)
                    {
                        u2 = tmp[i];
                        break;
                    }
                assert(u2 >= 0);
                used.remove(u1);
                used.remove(u2);
                int d1 = 0, d2 = 0;
                for (int w : adj[u1])
                    if (x[w] < 0 && used.get(w))
                        d1++;
                for (int w : adj[u2])
                    if (x[w] < 0 && used.get(w))
                        d2++;
                if (!(d1 < p - 2 && d2 < p - 2))
                {
                    for (int i = 0; i < p; i++)
                    {
                        int u = tmp[i];
                        if (u == u1 || u == u2)
                            continue;
                        int d = 0;
                        for (int w : adj[u])
                            if (x[w] < 0 && used.get(w))
                                d++;
                        if (d < p - 3)
                        {
                            goto loop;
                        }
                    }
                    int u = (d1 == p - 2) ? u2 : u1;
                    std::vector<int> const v1{v};
                    std::vector<int> const v2{u};
                    compute_alternative(v1, v2);
                }
                else
                {
                    int cnt = 0;
                    for (int i = 0; i < p; i++)
                    {
                        int d = 0;
                        for (int w : adj[tmp[i]])
                        {
                            if (x[w] < 0 && used.get(w))
                                d++;
                        }
                        if (d < p - 3)
                            cnt++;

                        if (cnt >= 3)
                            goto loop;
                    }

                    for (int i : adj[v])
                    {
                        if (x[i] < 0)
                        {
                            x[i] = 3;
                            if (checkFunnel(v))
                            {
                                funnel_vtcs.push_back(i);
                            }

                            x[i] = -1;
                        }
                    }
                }
            }
        loop:;
        }
    if (debug >= 3 && depth <= maxDepth && oldn != rn)
        fprintf(stderr, "%sfunnel: %d -> %d\n", debugString().c_str(), oldn, rn);
    return oldn != rn;
}

bool branch_and_reduce_algorithm::funnelReduction_b()
{
    int oldn = rn;
    for (int v = 0; v < n; v++)
        if (x[v] < 0)
        {
            used.clear();
            std::vector<int> &tmp = level;
            int p = 0;
            for (int u : adj[v])
                if (x[u] < 0 && used.add(u))
                {
                    tmp[p++] = u;
                }
            if (p <= 1) // deg 1 vtx.
            {
                set(v, 0);
                continue;
            }
            int u1 = -1;
            for (int i = 0; i < p; i++)
            {
                int d = 0;
                for (int u : adj[tmp[i]])
                    if (x[u] < 0 && used.get(u))
                        d++;
                if (d + 1 < p)
                {
                    u1 = tmp[i];
                    break;
                }
            }
            if (u1 < 0) // d >= p-1 for all neighbors => v and neighbours induce a clique
            {
                set(v, 0);
                continue;
            }
            else
            {
                std::vector<int> &id = iter;
                for (int i = 0; i < p; i++)
                    id[tmp[i]] = -1;
                for (int u : adj[u1])
                    if (x[u] < 0)
                        id[u] = 0;
                int u2 = -1;
                for (int i = 0; i < p; i++)
                    if (tmp[i] != u1 && id[tmp[i]] < 0)
                    {
                        u2 = tmp[i];
                        break;
                    }
                assert(u2 >= 0);
                used.remove(u1);
                used.remove(u2);
                int d1 = 0, d2 = 0;
                for (int w : adj[u1])
                    if (x[w] < 0 && used.get(w))
                        d1++;
                for (int w : adj[u2])
                    if (x[w] < 0 && used.get(w))
                        d2++;
                if (d1 < p - 2 && d2 < p - 2)
                    continue;

                int cnt = 0;
                int kv = -1;
                for (int i = 0; i < p; i++)
                {
                    int u = tmp[i];
                    if (u == u1 || u == u2)
                        continue;
                    int d = 0;
                    for (int w : adj[u])
                        if (x[w] < 0 && used.get(w))
                            d++;
                    if (d < p - 3)
                    {
                        if (cnt >= 1)
                            goto loop;
                        kv = u;
                        cnt++;
                    }
                }

                if (cnt == 1)
                {
                    funnel_vtcs.push_back(kv);
                    goto loop;
                }

                int u = (d1 == p - 2) ? u2 : u1;
                std::vector<int> const v1{v};
                std::vector<int> const v2{u};
                compute_alternative(v1, v2);
                //std::cout << "funnel: " << nBranchings << endl;
            }
        loop:;
        }
    if (debug >= 3 && depth <= maxDepth && oldn != rn)
        fprintf(stderr, "%sfunnel: %d -> %d\n", debugString().c_str(), oldn, rn);
    return oldn != rn;
}

bool branch_and_reduce_algorithm::deskReduction()
{
    int oldn = rn;
#if 1
    std::vector<int> &tmp = level;
    std::vector<int> &nv = iter;
    for (int i = 0; i < n; i++)
        nv[i] = -1;
    for (int v = 0; v < n; v++)
        if (x[v] < 0)
        {
            int d = 0;
            for (int u : adj[v])
                if (x[u] < 0)
                {
                    tmp[d++] = u;
                    nv[u] = v;
                    if (d > 4)
                        break;
                }
            if (d == 3 || d == 4)
            {
                int d2 = 0;
                for (int i = 0; i < d; i++)
                {
                    int a = deg(tmp[i]);
                    if (a == 3 || a == 4)
                        tmp[d2++] = tmp[i];
                }
                for (int i = 0; i < d2; i++)
                {
                    int u1 = tmp[i];
                    int sB1 = 0;
                    used.clear();
                    for (int w : adj[u1])
                        if (x[w] < 0 && w != v)
                        {
                            used.add(w);
                            sB1++;
                        }
                    for (int j = i + 1; j < d2; j++)
                    {
                        int u2 = tmp[j];
                        if (used.get(u2))
                            continue;
                        int sB2 = 0;
                        for (int w : adj[u2])
                            if (x[w] < 0 && w != v && !used.get(w))
                                sB2++;
                        if (sB1 + sB2 <= 3)
                        {
                            for (int w : adj[u2])
                                if (x[w] < 0 && used.get(w) && nv[w] != v)
                                {
                                    int d3 = deg(w);
                                    if (d3 == 3 || d3 == 4)
                                    {
                                        int sA = d - 2;
                                        for (int z : adj[w])
                                            if (x[z] < 0 && z != u1 && z != u2 && nv[z] != v)
                                            {
                                                sA++;
                                            }
                                        if (sA <= 2)
                                        {
                                            compute_alternative(std::vector<int>{v, w}, std::vector<int>{u1, u2});
                                            //std::cout << "desk: " << nBranchings << endl;
                                            goto loop;
                                        }
                                    }
                                }
                        } // if
                    }
                }
            }
        loop:;
        }
    if (debug >= 3 && depth <= maxDepth && oldn != rn)
        fprintf(stderr, "%sdesk: %d -> %d\n", debugString().c_str(), oldn, rn);
#endif // 0
    return oldn != rn;
}

bool branch_and_reduce_algorithm::unconfinedReduction()
{
    if (BRANCHING == 18)
        return unconfinedReduction_a();

    int oldn = rn;
#if 1
    std::vector<int> &NS = level;
    std::vector<int> &deg = iter;
    for (int v = 0; v < n; v++)
        if (x[v] < 0)
        {
            used.clear();
            used.add(v);
            int p = 1, size = 0;
            for (int u : adj[v])
                if (x[u] < 0)
                {
                    used.add(u);
                    NS[size++] = u;
                    deg[u] = 1;
                }
            bool ok = false;

            while (!ok)
            {
                ok = true;
                for (int i = 0; i < size; i++)
                {
                    int const u = NS[i];
                    if (deg[u] != 1)
                    {
                        continue;
                    }
                    int z = -1;
                    for (int const w : adj[u])
                        if (x[w] < 0 && !used.get(w))
                        {
                            if (z >= 0)
                            {
                                z = -2;
                                break;
                            }
                            z = w;
                        }
                    if (z == -1)
                    {
                        if (REDUCTION >= 3)
                        {
                            std::vector<int> &qs = que;
                            int q = 0;
                            qs[q++] = 1;
                            for (int w : adj[v])
                                if (x[w] < 0)
                                    qs[q++] = w;
                            std::vector<int> copyOfqs(qs.begin(), qs.begin() + q);
                            packing.emplace_back(std::move(copyOfqs));
                        }
                        set(v, 1);
                        goto whileloopend;
                    }
                    else if (z >= 0)
                    {
                        ok = false;
                        used.add(z);
                        p++;
                        for (int w : adj[z])
                            if (x[w] < 0)
                            {
                                if (used.add(w))
                                {
                                    NS[size++] = w;
                                    deg[w] = 1;
                                }
                                else
                                {
                                    deg[w]++;
                                }
                            }
                    }
                }
            }
        whileloopend:
            if (x[v] < 0 && p >= 2)
            {
                used.clear();
                for (int i = 0; i < size; i++)
                    used.add(NS[i]);
                std::vector<int> &vs = que;
                for (int i = 0; i < size; i++)
                {
                    vs[i] = vs[n + i] = -1;
                    int u = NS[i];
                    if (deg[u] != 2)
                        continue;
                    int v1 = -1, v2 = -1;
                    for (int w : adj[u])
                        if (x[w] < 0 && !used.get(w))
                        {
                            if (v1 < 0)
                                v1 = w;
                            else if (v2 < 0)
                                v2 = w;
                            else
                            {
                                v1 = v2 = -1;
                                break;
                            }
                        }
                    if (v1 > v2)
                    {
                        int t = v1;
                        v1 = v2;
                        v2 = t;
                    }
                    vs[i] = v1;
                    vs[n + i] = v2;
                }
                for (int i = 0; i < size; i++)
                    if (vs[i] >= 0 && vs[n + i] >= 0)
                    {
                        int u = NS[i];
                        used.clear();
                        for (int w : adj[u])
                            if (x[w] < 0)
                                used.add(w);
                        for (int j = i + 1; j < size; j++)
                            if (vs[i] == vs[j] && vs[n + i] == vs[n + j] && !used.get(NS[j]))
                            {
                                if (REDUCTION >= 3)
                                {
                                    std::vector<int> &qs = que;
                                    int q = 0;
                                    qs[q++] = 1;
                                    for (int w : adj[v])
                                        if (x[w] < 0)
                                            qs[q++] = w;
                                    std::vector<int> copyOfqs(qs.begin(), qs.begin() + q);
                                    packing.emplace_back(std::move(copyOfqs));
                                }
                                set(v, 1);
                                goto forloopend;
                            }
                    }
            forloopend:;
            }
        }
    if (debug >= 3 && depth <= maxDepth && oldn != rn)
        fprintf(stderr, "%sunconfined: %d -> %d\n", debugString().c_str(), oldn, rn);
#endif //0
    return oldn != rn;
}

bool branch_and_reduce_algorithm::unconfinedReduction_a()
{
    int oldn = rn;
#if 1
    std::vector<int> &NS = level;
    std::vector<int> &deg = iter;
    for (int v = 0; v < n; v++)
        if (x[v] < 0)
        {
            used.clear();
            used.add(v); // add v to set S
            int p = 1, size = 0;
            for (int u : adj[v]) // add N(v) to set S
                if (x[u] < 0)
                {
                    used.add(u);
                    NS[size++] = u;
                    deg[u] = 1;
                }
            bool ok = false;

            while (!ok) //
            {
                ok = true;
                std::vector<int> extends;

                for (int i = 0; i < size; i++)
                {
                    int const u = NS[i];
                    if (deg[u] != 1)
                    {
                        continue;
                    }
                    int z = -1;
                    for (int const w : adj[u])
                        if (x[w] < 0 && !used.get(w)) // w is in N(u) but not N(S)
                        {
                            if (z >= 0)
                            {
                                z = -2; // |N(u)\N(S)| >= 2
                                break;
                            }
                            z = w;
                        }
                    if (z == -1) // u is child with N(u)\N(S) = 0
                    {
                        if (REDUCTION >= 3)
                        {
                            std::vector<int> &qs = que;
                            int q = 0;
                            qs[q++] = 1;
                            for (int w : adj[v])
                                if (x[w] < 0)
                                    qs[q++] = w;
                            std::vector<int> copyOfqs(qs.begin(), qs.begin() + q);
                            packing.emplace_back(std::move(copyOfqs));
                        }
                        set(v, 1);
                        goto whileloopend;
                    }
                    else if (z >= 0) // u is extending child
                    {
                        extends.push_back(z);
                    }
                }

                if (extends.size() > 0)
                {
                    if (extends.size() == 1)
                    {
                        unconf_vtcs.push_back(extends[0]);
                    }
                    int z = extends[0];
                    ok = false;
                    used.add(z);
                    p++;
                    for (int w : adj[z])
                        if (x[w] < 0)
                        {
                            if (used.add(w))
                            {
                                NS[size++] = w;
                                deg[w] = 1;
                            }
                            else
                            {
                                deg[w]++;
                            }
                        }
                }
            }
        whileloopend:
            if (x[v] < 0 && p >= 2)
            {
                used.clear();
                for (int i = 0; i < size; i++)
                    used.add(NS[i]);
                std::vector<int> &vs = que;
                for (int i = 0; i < size; i++)
                {
                    vs[i] = vs[n + i] = -1;
                    int u = NS[i];
                    if (deg[u] != 2)
                        continue;
                    int v1 = -1, v2 = -1;
                    for (int w : adj[u])
                        if (x[w] < 0 && !used.get(w))
                        {
                            if (v1 < 0)
                                v1 = w;
                            else if (v2 < 0)
                                v2 = w;
                            else
                            {
                                v1 = v2 = -1;
                                break;
                            }
                        }
                    if (v1 > v2)
                    {
                        int t = v1;
                        v1 = v2;
                        v2 = t;
                    }
                    vs[i] = v1;
                    vs[n + i] = v2;
                }
                for (int i = 0; i < size; i++)
                    if (vs[i] >= 0 && vs[n + i] >= 0)
                    {
                        int u = NS[i];
                        used.clear();
                        for (int w : adj[u])
                            if (x[w] < 0)
                                used.add(w);
                        for (int j = i + 1; j < size; j++)
                            if (vs[i] == vs[j] && vs[n + i] == vs[n + j] && !used.get(NS[j]))
                            {
                                if (REDUCTION >= 3)
                                {
                                    std::vector<int> &qs = que;
                                    int q = 0;
                                    qs[q++] = 1;
                                    for (int w : adj[v])
                                        if (x[w] < 0)
                                            qs[q++] = w;
                                    std::vector<int> copyOfqs(qs.begin(), qs.begin() + q);
                                    packing.emplace_back(std::move(copyOfqs));
                                }
                                set(v, 1);
                                goto forloopend;
                            }
                    }
            forloopend:;
            }
        }
    if (debug >= 3 && depth <= maxDepth && oldn != rn)
        fprintf(stderr, "%sunconfined: %d -> %d\n", debugString().c_str(), oldn, rn);
#endif //0
    return oldn != rn;
}

int branch_and_reduce_algorithm::packingReduction()
{
    int oldn = rn;
#if 1
    std::vector<int> x2(x);
    int a = -1;
    for (unsigned int pi = 0; pi < packing.size(); ++pi)
    {
        std::vector<int> &ps = packing[pi];
        if (a != rn)
        {
            for (int j = 0; j < N; j++)
                x2[j] = x[j];
            for (int j = modifiedN - 1; j >= 0; j--)
            {
                modifieds[j]->reverse(x2);
            }
            a = rn;
        }
        int max = ps.size() - 1 - ps[0], sum = 0, size = 0;
        std::vector<int> &S = level;
        for (unsigned int j = 1; j < ps.size(); j++)
        {
            int v = ps[j];
            if (x2[v] < 0)
                S[size++] = v;
            if (x2[v] == 1)
                sum++;
        }
        if (sum > max)
        {
            return -1;
        }
        else if (sum == max && size > 0)
        {
            std::vector<int> &count = iter;
            used.clear();
            for (int j = 0; j < size; j++)
            {
                used.add(S[j]);
                count[S[j]] = -1;
            }
            for (int j = 0; j < size; j++)
            {
                for (int u : adj[S[j]])
                    if (x[u] < 0)
                    {
                        if (used.add(u))
                        {
                            count[u] = 1;
                        }
                        else if (count[u] < 0)
                        {
                            return -1;
                        }
                        else
                        {
                            count[u]++;
                        }
                    }
            }
            for (int j = 0; j < size; j++)
            {
                for (int u : adj[S[j]])
                    if (x[u] < 0 && count[u] == 1)
                    {
                        std::vector<int> &tmp = que;
                        int p = 0;
                        tmp[p++] = 1;
                        for (int w : adj[u])
                            if (x[w] < 0 && !used.get(w))
                            {
                                tmp[p++] = w;
                            }
                        std::vector<int> copyOfTmp(tmp.begin(), tmp.begin() + p);
                        packing.emplace_back(std::move(copyOfTmp));
                    }
            }
            for (int j = 0; j < size; j++)
            {
                if (S[j] == 1)
                    return -1;
                assert(x[S[j]] < 0);
                set(S[j], 0);
            }
        }
        else if (sum + size > max)
        {
            assert(size >= 2);
            used.clear();
            for (int j = 0; j < size; j++)
                used.add(S[j]);
            for (int v : adj[S[0]])
                if (x[v] < 0 && !used.get(v))
                {
                    int p = 0;
                    for (int u : adj[v])
                        if (used.get(u))
                            p++;
                    if (sum + p > max)
                    {
                        std::vector<int> &qs = que;
                        int q = 0;
                        qs[q++] = 2;
                        for (int u : adj[v])
                            if (x[u] < 0)
                                qs[q++] = u;
                        std::vector<int> copyOfqs(qs.begin(), qs.begin() + q);
                        packing.emplace_back(std::move(copyOfqs));
                        set(v, 1);
                        break;
                    }
                }
        }
    }
    if (debug >= 3 && depth <= maxDepth && oldn != rn)
        fprintf(stderr, "%spacking: %d -> %d\n", debugString().c_str(), oldn, rn);
#endif // 0
    return oldn != rn ? 1 : 0;
}

void branch_and_reduce_algorithm::branching(timer &t, double time_limit)
{
    int oldLB = lb;
    int v = -1, dv = 0;
    std::vector<int> &mirrors = que;
    int mirrorN = 0;

    if (BRANCHING == 0) // branch on random vertex
    {
        int p = rand() % rn;
        for (int i = 0; i < n; i++)
            if (x[i] < 0 && p-- == 0)
                v = i;
        dv = deg(v);
    }
    else if (BRANCHING == 1) // branch on vertex with min. deg.
    {
        dv = n + 1;
        for (int u = 0; u < n; u++)
            if (x[u] < 0)
            {
                int degree = deg(u);
                if (dv > degree)
                {
                    v = u;
                    dv = degree;
                }
            }
    }
    else if (BRANCHING == 2) // branch on vertex with max. deg. && min. edges among N(v)
    {
        v = get_max_deg_vtx();
        dv = deg(v);
    }
    else if (BRANCHING == 3) // articulation points
    {
        int cv = get_articulation_point();

        if (cv == -1)
        {
            v = get_max_deg_vtx();
            defaultPicks++;
            defaultBranch = true;
        }
        else
        {
            stratPicks++;
            v = cv;
        }

        dv = deg(v);
    }
    else if (BRANCHING == 4) // global mincut
    {
        v = get_mincut_vertex();

        if (v == -1)
            std::cout << "FAILED !!" << std::endl;

        dv = deg(v);
    }
    else if (BRANCHING == 5) // st cuts
    {
        while (!cut.empty() && x[cut.back()] != -1)
            cut.pop_back();

        if (cut.empty())
        {
            if (branch_t == 0)
            {
                get_stcut_vertices();
                branch_t = TUNING_PARAM3;
            }
            else
            {
                cut.push_back(get_max_deg_vtx());
                branch_t--;

                defaultPicks++;
                defaultBranch = true;
            }
        }

        v = cut.back();
        cut.pop_back();
        dv = deg(v);
    }
    else if (BRANCHING == 55) // st cuts
    {
        while (!cut.empty() && x[cut.back()] != -1)
            cut.pop_back();

        if (cut.empty())
        {
            if (branch_t == 0)
            {
                get_stcut_vertices();
                branch_t = TUNING_PARAM3;
            }
            else
            {
                cut.push_back(get_max_deg_vtx());
                branch_t--;

                defaultPicks++;
                defaultBranch = true;
            }
        }

        v = cut.back();
        cut.pop_back();
        dv = deg(v);
    }
    else if (BRANCHING == 6) // nested dissection
    {
        if (nd_computed == false)
        {
            nd_computed = true;
            compute_nd_order();
        }

        while (!nd_order.empty() && x[nd_order.back()] != -1)
            nd_order.pop_back();

        branch_t--;

        if (nd_order.empty())
        {
            //cout << "nd is empty!" << endl;
            nd_order.push_back(get_max_deg_vtx());
            defaultPicks++;
            defaultBranch = true;
        }
        else
            stratPicks++;

        v = nd_order.back();
        nd_order.pop_back();
        dv = deg(v);

        dv = deg(v);
    }
    else if (BRANCHING == 66) 
    {
        if (nd_computed == false && branch_t == 0)
        {
            nd_computed = true;
            compute_nd_order();
            branch_t = 20;
        }

        while (!nd_order.empty() && x[nd_order.back()] != -1)
            nd_order.pop_back();

        branch_t--;

        if (nd_order.empty())
        {
            //cout << "nd is empty!" << endl;
            nd_order.push_back(get_max_deg_vtx());
            defaultPicks++;
            defaultBranch = true;
        }
        else
            stratPicks++;

        v = nd_order.back();
        nd_order.pop_back();
        dv = deg(v);
    }

    else if (BRANCHING == 13) // max N2
    {
        v = max_nh_vtx();
        dv = deg(v);
    }
    else if (BRANCHING == 14) // almost dominated
    {
        almost_dominated();

        int vv = -1;
        int dvv = -1;
        for (int i = 0; i < domin_vtcs.size(); i++)
        {
            if (x[domin_vtcs[i]] < 0 && deg(domin_vtcs[i]) > dvv)
            {
                dvv = deg(domin_vtcs[i]);
                vv = domin_vtcs[i];
            }
        }

        if (vv == -1)
        {
            v = get_max_deg_vtx();
            defaultPicks++;
        }
        else
        {
            int vvv = get_max_deg_vtx();
            if (deg(vv) >= deg(vvv) - TUNING_PARAM1)
            {
                v = vv;
                stratPicks++;
            }
            else
                v = vvv;
        }

        domin_vtcs.clear();
        dv = deg(v);
    }
    else if (BRANCHING == 15) // almost twin
    {
        int vv = -1;
        int dvv = -1;
        for (int i = 0; i < twin_vtcs.size(); i++)
        {
            if (x[twin_vtcs[i]] < 0 && deg(twin_vtcs[i]) > dvv)
            {
                dvv = deg(twin_vtcs[i]);
                vv = twin_vtcs[i];
            }
        }

        if (vv == -1)
        {
            v = get_max_deg_vtx();
            defaultPicks++;
        }
        else
        {
            int vvv = get_max_deg_vtx();
            if (deg(vv) >= deg(vvv) - TUNING_PARAM1)
            {
                v = vv;
                stratPicks++;
            }
            else
                v = vvv;
        }

        twin_vtcs.clear();
        dv = deg(v);
    }
    else if (BRANCHING == 16) // funnel
    {
        int vv = -1;
        int dvv = -1;
        for (int i = 0; i < funnel_vtcs.size(); i++)
        {
            if (x[funnel_vtcs[i]] < 0 && deg(funnel_vtcs[i]) > dvv)
            {
                dvv = deg(funnel_vtcs[i]);
                vv = funnel_vtcs[i];
            }
        }

        if (vv == -1)
        {
            v = get_max_deg_vtx();
            defaultPicks++;
        }
        else
        {
            int vvv = get_max_deg_vtx();
            if (deg(vv) >= deg(vvv) - TUNING_PARAM1)
            {
                v = vv;
                stratPicks++;
            }
            else
                v = vvv;
        }

        funnel_vtcs.clear();
        dv = deg(v);
    }
    else if (BRANCHING == 17) // unconfined
    {
        int vv = -1;
        int dvv = -1;
        for (int i = 0; i < unconf_vtcs.size(); i++)
        {
            if (x[unconf_vtcs[i]] < 0 && deg(unconf_vtcs[i]) > dvv)
            {
                dvv = deg(unconf_vtcs[i]);
                vv = unconf_vtcs[i];
            }
        }

        if (vv == -1)
        {
            v = get_max_deg_vtx();
            defaultPicks++;
        }
        else
        {
            int vvv = get_max_deg_vtx();
            if (deg(vv) >= deg(vvv) - TUNING_PARAM1)
            {
                v = vv;
                stratPicks++;
            }
            else
                v = vvv;
        }

        unconf_vtcs.clear();
        dv = deg(v);
    }

    else if (BRANCHING == 18)
    {
        build_domination_graph();
        find_chains();
        calc_chain_vec();
        int maxdel = 0;
        pair<int, int> vec;
        v = -1;
        for (int i = 0; i < n; i++)
        {
            if (x[i] < 0 && chain_vec[i].first + chain_vec[i].second > maxdel)
            {
                vec = chain_vec[i];
                maxdel = chain_vec[i].first + chain_vec[i].second;
                v = i;
            }
        }

        dv = deg(v);
    }


    int crntBest = opt;

    std::vector<int> &ps = iter;
    for (int i = 0; i < n; i++)
        ps[i] = -2;
    used.clear();
    used.add(v);
    for (int u : adj[v])
        if (x[u] < 0)
        {
            used.add(u);
            ps[u] = -1;
        }
    for (int u : adj[v]) // find mirrors
        if (x[u] < 0)
        {
            for (int w : adj[u])
                if (x[w] < 0 && used.add(w))
                {
                    int c1 = dv;
                    for (int z : adj[w])
                        if (x[z] < 0 && ps[z] != -2)
                        {
                            ps[z] = w;
                            c1--;
                        }
                    bool ok = true;
                    for (int u2 : adj[v])
                        if (x[u2] < 0 && ps[u2] != w)
                        {
                            int c2 = 0;
                            for (int w2 : adj[u2])
                                if (x[w2] < 0 && ps[w2] == w)
                                    c2++;
                            if (c2 != c1 - 1)
                            {
                                ok = false;
                                break;
                            }
                        }
                    if (ok)
                        mirrors[mirrorN++] = w;
                }
        }
    int pn = rn;
    unsigned int oldP = packing.size(); // update packing constr.
    if (REDUCTION >= 3)
    {
        std::vector<int> &tmp = level;
        int p = 0;
        tmp[p++] = mirrorN > 0 ? 2 : 1;
        for (int u : adj[v])
            if (x[u] < 0)
                tmp[p++] = u;
        std::vector<int> copyOfTmp(tmp.begin(), tmp.begin() + p);
        packing.emplace_back(std::move(copyOfTmp));
    }
    set(v, 1);
    for (int i = 0; i < mirrorN; i++)
        set(mirrors[i], 1);
    if (debug >= 2 && depth <= maxDepth)
    {
        if (mirrorN > 0)
            fprintf(stderr, "%sbranchMirror (%d, %d): 1\n", debugString().c_str(), dv, mirrorN);
        else
            fprintf(stderr, "%sbranch (%d): 1\n", debugString().c_str(), dv);
    }
    depth++;
    rec(t, time_limit);
    while (packing.size() > oldP)
        packing.pop_back();
    lb = oldLB;
    depth--;
    restore(pn);

    // optimal branch order
    crntBest = opt;

    if (lb >= opt)
    {
        if (startingSolutionIsBest)
        {
            ++numBranchesPrunedByStartingSolution;
        }
        return;
    }
    nBranchings++;
    if (defaultBranch)
    {
        defaultBranchings++;
        defaultBranch = false;
    }

    if (mirrorN == 0)
    {
        used.clear();
        used.add(v);
        for (int u : adj[v])
            if (x[u] < 0)
                used.add(u);
        if (REDUCTION >= 3)
        {
            std::vector<int> ws(n, -1);
            for (int u : adj[v])
                if (x[u] < 0)
                {
                    std::vector<int> &tmp = level;
                    int p = 0;
                    tmp[p++] = 1;
                    for (int w : adj[u])
                        if (x[w] < 0 && !used.get(w))
                        {
                            tmp[p++] = w;
                            ws[w] = u;
                        }
                    assert(p >= 2);
                    for (int u2 : adj[tmp[1]])
                        if (x[u2] < 0 && used.get(u2) && u2 != u)
                        {
                            int c = 0;
                            for (int w : adj[u2])
                                if (x[w] < 0)
                                {
                                    if (ws[w] == u)
                                        c++;
                                    else if (w == u || !used.get(w))
                                    {
                                        c = -1;
                                        break;
                                    }
                                }
                            if (c == p - 1)
                            {
                                tmp[0] = 2;
                                break;
                            }
                        }
                    std::vector<int> copyOfTmp(tmp.begin(), tmp.begin() + p);
                    packing.emplace_back(std::move(copyOfTmp));
                }
        }
    }
    set(v, 0);
    if (debug >= 2 && depth <= maxDepth)
        fprintf(stderr, "%sbranch (%d): 0\n", debugString().c_str(), dv);
    depth++;
    rec(t, time_limit);
    while (packing.size() > oldP)
        packing.pop_back();
    lb = oldLB;
    depth--;
    restore(pn);
}

bool branch_and_reduce_algorithm::decompose(timer &t, double time_limit)
{
    std::vector<std::vector<int>> vss; // components
    {
        std::vector<int> &id = level;  // representative of the component
        std::vector<int> &size = iter; // size of the components
        int nC = 0;
        {
            for (int i = 0; i < n; i++)
                id[i] = -1;
            for (int s = 0; s < n; s++)
                if (x[s] < 0 && id[s] < 0) // new components
                {
                    nC++;
                    int qs = 0, qt = 0;
                    que[qt++] = s;
                    id[s] = s;
                    while (qs < qt) // DFS -> find component
                    {
                        int v = que[qs++];
                        for (int u : adj[v])
                            if (x[u] < 0 && id[u] < 0)
                            {
                                id[u] = s;
                                que[qt++] = u;
                            }
                    }
                    size[s] = qt;
                }
        }
        if (nC <= 1 && (n <= 100 || n * SHRINK < rn))
            return false;

        // count decomps
        nDecomps++;

        std::vector<long long> cs(nC, 0);
        {
            int p = 0;
            for (int i = 0; i < n; i++)
                if (x[i] < 0 && id[i] == i)
                {
                    cs[p++] = ((long long)(size[i])) << 32 | i;
                }
            std::sort(cs.begin(), cs.end());
        }
        vss.resize(nC);
        std::vector<int> qs(n, 0);
        {
            for (int i = 0; i < nC; i++)
            {
                vss[i].resize(size[(int)cs[i]]);
                qs[(int)cs[i]] = i;
            }
            std::vector<int> ps(nC);
            for (int i = 0; i < n; i++)
                if (x[i] < 0)
                {
                    int j = qs[id[i]];
                    vss[j][ps[j]++] = i;
                }
        }
        for (int i = 0; i < n; i++)
            id[i] = -1;
        for (unsigned int i = 0; i < vss.size(); i++)
        {
            std::vector<int> &vs = vss[i];
            std::vector<long long> ls(vs.size());
            for (unsigned int j = 0; j < vs.size(); j++)
                ls[j] = ((long long)(n - deg(vs[j]))) << 32 | vs[j];
            std::sort(ls.begin(), ls.end());
            for (unsigned int j = 0; j < vs.size(); j++)
                vs[j] = (int)ls[j];
        }
    }
    std::vector<int> x2(x);
    for (int i = modifiedN - 1; i >= 0; i--)
        modifieds[i]->reverse(x2);
    std::vector<int> size(vss.size());
    for (unsigned int i = 0; i < vss.size(); i++)
        size[i] = vss[i].size();
    std::vector<int> pos1(N, -1);
    std::vector<int> pos2(N, 0);
    std::vector<std::vector<int>> packingB;
    {
        for (unsigned int i = 0; i < vss.size(); i++)
        {
            for (unsigned int j = 0; j < vss[i].size(); j++)
            {
                pos1[vss[i][j]] = i;
                pos2[vss[i][j]] = j;
            }
        }
        std::vector<bool> need(N, false);
        for (std::vector<int> const &ps : packing)
        {

            int max = ps.size() - 1 - ps[0], sum = 0, count = 0;
            for (unsigned int j = 1; j < ps.size(); j++)
            {
                int v = ps[j];
                if (x2[v] < 0 || x2[v] == 2)
                {
                    count++;
                }
                if (x2[v] == 1)
                    sum++;
            }
            if (sum > max)
                return true;
            if (sum + count > max)
            {
                packingB.push_back(ps);
                for (unsigned int k = 1; k < ps.size(); k++)
                {
                    if (x2[ps[k]] == 2)
                        need[ps[k]] = true;
                }
            }
        }
        for (int i = 0; i < modifiedN; i++)
        {
            bool b = false;
            shared_ptr<modified> mod = modifieds[i];
            for (int v : mod->removed)
                if (need[v])
                    b = true;
            if (b)
            {
                if (dynamic_cast<fold *>(mod.get()) != nullptr)
                {
                    if (x2[mod->vs[0]] == 2)
                        need[mod->vs[0]] = true;
                }
                else
                {
                    for (int v : mod->vs)
                        if (x2[v] == 2)
                        {
                            need[v] = true;
                        }
                }
            }
        }
        for (int i = modifiedN - 1; i >= 0; i--)
        {
            shared_ptr<modified> mod = modifieds[i];
            bool b = false;
            for (int v : mod->removed)
                if (need[v])
                    b = true;
            if (b)
            {
                if (dynamic_cast<fold *>(mod.get()) != nullptr)
                {
                    for (int v : mod->removed)
                    {
                        assert(pos1[v] == -1);
                        pos1[v] = pos1[mod->vs[0]];
                        assert(pos1[v] >= 0);
                        pos2[v] = size[pos1[v]]++;
                    }
                }
                else
                {
                    int max = -1;
                    for (int v : mod->vs)
                        if (max < pos1[v])
                            max = pos1[v];
                    assert(max >= 0);
                    for (int v : mod->removed)
                    {
                        assert(pos1[v] == -1);
                        pos1[v] = max;
                        pos2[v] = size[pos1[v]]++;
                    }
                }
            }
        }
        for (int i = 0; i < n; i++)
        {
            if ((x2[i] == 0 || x2[i] == 1) && pos1[i] >= 0)
            {
                assert(false);
            }
        }
    }
    std::vector<branch_and_reduce_algorithm *> vcs(vss.size(), nullptr); // create new VCSolver instances
    {
        for (int i = 0; i < static_cast<int>(vss.size()); i++)
        {
            std::vector<int> &vs = vss[i];
            size[i] += 2;
            std::vector<std::vector<int>> adj2(vs.size()); // create new adj. list
            for (int j = 0; j < static_cast<int>(vs.size()); j++)
            {
                adj2[j].resize(deg(vs[j]), 0);
                int p = 0;
                for (int u : adj[vs[j]])
                    if (x[u] < 0)
                        adj2[j][p++] = pos2[u];
                assert(p == adj2[j].size());
                std::sort(adj2[j].begin(), adj2[j].end());
            }

            vcs[i] = new branch_and_reduce_algorithm(adj2, size[i]);

            // inherit nd branching order ypp
            std::vector<int> sub_nd_order(0);
            for (int i = 0; i < this->nd_order.size(); i++)
            {
                if (std::find(vs.begin(), vs.end(), this->nd_order[i]) != vs.end())
                {
                    sub_nd_order.push_back(pos2[this->nd_order[i]]);
                }
            }

            vcs[i]->nd_computed = this->nd_computed;
            vcs[i]->nd_order.swap(sub_nd_order);

            /*if (BRANCHING == 8 && bc_index_built)
            {
                std::vector<int> sub_nodeMapping(nodeMapping.size());
                for (int i = 0; i < vs.size(); i++)
                    sub_nodeMapping[i] = nodeMapping[vs[i]];

                vcs[i]->bc_index_built = this->bc_index_built;
                vcs[i]->bc_index = this->bc_index;
                vcs[i]->nodeMapping.swap(sub_nodeMapping);
            }*/

            for (unsigned int j = 0; j < vs.size(); j++)
            {
                if (in[vs[j]] >= 0 && pos1[in[vs[j]]] == i && pos2[in[vs[j]]] < static_cast<int>(vs.size()))
                {
                    vcs[i]->in[j] = pos2[in[vs[j]]];
                }
                if (out[vs[j]] >= 0 && pos1[out[vs[j]]] == i && pos2[out[vs[j]]] < static_cast<int>(vs.size()))
                {
                    vcs[i]->out[j] = pos2[out[vs[j]]];
                }
            }
            vcs[i]->x[vcs[i]->N - 2] = vcs[i]->y[vcs[i]->N - 2] = 0;
            vcs[i]->x[vcs[i]->N - 1] = vcs[i]->y[vcs[i]->N - 1] = 1;
        }
    }
    {
        for (unsigned int i = 0; i < packingB.size(); i++)
        {
            std::vector<int> &ps(packingB[i]);
            int maxID = -1;
            for (unsigned int j = 1; j < ps.size(); j++)
            {
                int v = ps[j];
                if (x2[v] < 0 || x2[v] == 2)
                {
                    maxID = std::max(maxID, pos1[v]);
                }
            }
            vcs[maxID]->packing.push_back(ps);
        }
    }
    {
        for (int i = 0; i < modifiedN; i++)
        {
            shared_ptr<modified> mod = modifieds[i];
            int p = pos1[mod->removed[0]];
            if (p >= 0)
            {
                vcs[p]->modifieds[vcs[p]->modifiedN++] = mod;
            }
        }
    }
    std::vector<std::vector<int>> vss2(vss.size());
    {
        for (unsigned int i = 0; i < vss.size(); i++)
            vss2[i].resize(vcs[i]->N - 2, 0);
        for (int i = 0; i < N; i++)
            if (pos1[i] >= 0)
                vss2[pos1[i]][pos2[i]] = i;
    }
    int sum = crt;
    // std::vector<bool> optChanged(vss.size(), false);

    for (int i = 0; i < static_cast<int>(vss.size()) && opt > sum; i++)
    {
        branch_and_reduce_algorithm *vc = vcs[i];
        {
            std::vector<std::vector<int>> packing2;
            for (std::vector<int> const &ps : vc->packing)
            {
                std::vector<int> &tmp = level;
                int p = 0;
                tmp[p++] = ps[0];
                for (unsigned int k = 1; k < ps.size(); k++)
                {
                    int v = ps[k];
                    if (pos1[v] == i)
                    {
                        tmp[p++] = pos2[v];
                    }
                    else
                    {
                        assert(x2[v] == 0 || x2[v] == 1);
                        if (x2[v] == 0)
                            tmp[0]--;
                    }
                }
                if (p - 1 < tmp[0])
                    return true;
                if (tmp[0] <= 0)
                    continue;
                std::vector<int> copyOfTmp(tmp.begin(), tmp.begin() + p);
                packing2.emplace_back(std::move(copyOfTmp));
            }
            (vc->packing).swap(packing2);
        }
        {
            for (int j = 0; j < vc->modifiedN; j++)
            {
                shared_ptr<modified> mod = vc->modifieds[j];
                std::vector<int> removed(mod->removed.size());
                for (unsigned int k = 0; k < removed.size(); k++)
                {
                    int v = mod->removed[k];
                    assert(pos1[v] == i);
                    removed[k] = pos2[v];
                }
                if (dynamic_cast<fold *>(mod.get()) != nullptr)
                {
                    std::vector<int> vs(1, 0);
                    int v = mod->vs[0];
                    if (pos1[v] == i)
                    {
                        vs[0] = pos2[v];
                    }
                    else
                    {
                        assert(x2[v] == 0 || x2[v] == 1);
                        vs[0] = vc->N - 2 + x2[v];
                    }
                    mod = make_shared<fold>(fold(removed, vs, this));
                }
                else
                {
                    std::vector<int> vs(mod->vs.size(), 0);
                    for (unsigned int k = 0; k < vs.size(); k++)
                    {
                        int v = mod->vs[k];
                        if (pos1[v] == i)
                        {
                            vs[k] = pos2[v];
                        }
                        else
                        {
                            assert(x2[v] == 0 || x2[v] == 1);
                            vs[k] = vc->N - 2 + x2[v];
                        }
                    }
                    mod = make_shared<alternative>(alternative(removed, vs, this, dynamic_cast<alternative *>(mod.get())->k));
                }
                vc->modifieds[j] = mod;
            }
        }
        vc->depth = depth + (vss.size() > 1 ? 1 : 0);
        if (debug >= 2 && depth <= maxDepth)
        {
            if (vss.size() == 1)
                fprintf(stderr, "%sshrink: %d -> %d (%d)\n", debugString().c_str(), n, vcs[i]->n, vcs[i]->N);
            else
                fprintf(stderr, "%sdecompose: %d (%d)\n", debugString().c_str(), vcs[i]->n, vcs[i]->N);
        }
        if (i + 1 == static_cast<int>(vss.size()))
        {
            vc->opt = std::min(vss[i].size(), static_cast<size_t>(opt - sum));
        }

        vc->reverse();
        for (int j = 0; j < vc->N; j++)
            assert(vc->y[j] == 0 || vc->y[j] == 1);

        // Map current optimal solution to CC
        // int optForMapping = 0;
        // for (unsigned int j = 0; j < vss[i].size(); j++) if(y[vss[i][j]] && vcs[i]->x[j] == -1) optForMapping++;
        // int currentOpt = 0;
        // for (unsigned int j = 0; j < vss[i].size(); j++) if(vcs[i]->y[j] && vcs[i]->x[j] == -1) currentOpt++;
        // if(currentOpt > optForMapping) {
        //     for (unsigned int j = 0; j < vss[i].size(); j++) if(vcs[i]->x[j] == -1) vcs[i]->y[j] = y[vss[i][j]];
        //     optForMapping = 0;
        //     // for (unsigned int j = 0; j < vss[i].size(); j++) if(vcs[i]->y[j]) optForMapping++;
        //     for (unsigned int j = 0; j < vcs[i]->N; j++) if(vcs[i]->y[j]) optForMapping++;
        //     vc->opt = optForMapping;
        //     vc->numBranchesPrunedByStartingSolution = 0;
        //     vc->startingSolutionIsBest = true;
        // }
        // // for(int v : vss[i]) {
        // //     assert(pos1[v] == i);
        // //     vc->y[pos2[i]] = y[i];
        // //     if(y[i] == 1) {
        // //         ++optForMapping;
        // //     }
        // // }
        // // // std::cout << "Setting opt for subproblem: " << optForMapping << std::endl;
        // vc->numBranchesPrunedByStartingSolution = 0;
        // vc->startingSolutionIsBest = startingSolutionIsBest;

        vc->solve(t, time_limit);
        sum += vc->opt;

        // if(vc->opt != optForMapping) {
        //     optChanged[i] = true;
        // }

        // numBranchesPrunedByStartingSolution += vc->numBranchesPrunedByStartingSolution;

        for (int j = 0; j < vc->N - 2; j++)
        {
            x2[vss2[i][j]] = vc->y[j];
            assert(vc->y[j] == 0 || vc->y[j] == 1);
        }
    }

    for (int i = 0; i < vcs.size(); i++)
    {
        branch_and_reduce_algorithm *pAlg = vcs[i];
        std::vector<int> &vs = vss[i];
    }

    if (opt > sum) // new best solution found -> set solution
    {
        if (debug >= 2 && rootDepth <= maxDepth)
            fprintf(stderr, "%sopt: %d -> %d\n", debugString().c_str(), opt, sum);
        opt = sum;
        y = x;
        startingSolutionIsBest = false;
        for (unsigned int i = 0; i < vss.size(); i++)
        {
            // if(!optChanged[i]) continue;
            for (unsigned int j = 0; j < vss[i].size(); j++)
                y[vss[i][j]] = vcs[i]->y[j];
        }

        reverse();
    }

    for (branch_and_reduce_algorithm *pAlg : vcs) // clean up solver
    {
        delete pAlg;
        pAlg = nullptr;
    }

    return true;
}

bool branch_and_reduce_algorithm::reduce()
{
    int oldn = rn;
    for (;;)
    {
        if (REDUCTION >= 0)
            deg1Reduction();
        // if (n > 100 && n * SHRINK >= rn && !outputLP && decompose()) return true;
        if (REDUCTION >= 0 && REDUCTION < 2 && dominateReduction())
            continue;

        if (REDUCTION >= 2 && unconfinedReduction())
            continue;
        if (REDUCTION >= 1 && lpReduction())
            continue;
        if (REDUCTION >= 3)
        {
            int r = packingReduction();
            if (r < 0)
                return true;
            if (r > 0)
                continue;
        }
        if (REDUCTION >= 1 && fold2Reduction())
            continue;
        if (REDUCTION >= 2 && twinReduction())
            continue;
        if (REDUCTION >= 2 && funnelReduction())
            continue;
        if (REDUCTION >= 2 && deskReduction())
            continue;
        break;
    }
    if (debug >= 2 && depth <= maxDepth && oldn != rn)
        fprintf(stderr, "%sreduce: %d -> %d\n", debugString().c_str(), oldn, rn);
    return false;
}

void branch_and_reduce_algorithm::rec(timer &t, double time_limit)
{
    if (t.elapsed() >= time_limit)
        return;
    if (REDUCTION < 3)
        assert(packing.size() == 0);

    if (EXTRA_DECOMP == 1)
    {
        if (decompose(t, time_limit)) // check for CC's
            return;
    }

    if (reduce()) // kernelization
        return;
    if (lowerBound() >= opt) // pruned by LowerBound
    {
        prunes++;
        if (startingSolutionIsBest && rn != 0)
        {
            ++numBranchesPrunedByStartingSolution;
        }
        return;
    }
    if (rn == 0) // G is empty
    {
        if (debug >= 2 && rootDepth <= maxDepth)
            fprintf(stderr, "%sopt: %d -> %d\n", debugString().c_str(), opt, crt);
        opt = crt;
        y = x;
        startingSolutionIsBest = false;
        reverse();
        return;
    }
    if (decompose(t, time_limit)) // check for CC's
        return;
    branching(t, time_limit); // branch
}

void branch_and_reduce_algorithm::addStartingSolution(std::vector<int> solution, int solutionSize)
{
    if (solution.size() != y.size())
    {
        cout << "ERROR: invalid solution std::vector!" << endl
             << flush;
    }
    for (int i = 0; i < solution.size(); ++i)
    {
        y[i] = solution[i];
    }
    opt = solutionSize;

    startingSolutionIsBest = true;
    numBranchesPrunedByStartingSolution = 0;
}

int branch_and_reduce_algorithm::solve(timer &t, double time_limit)
{
    if (t.elapsed() >= time_limit)
        return -1;    

    // PrintState();
    if (LOWER_BOUND >= 2 && REDUCTION <= 0 && !outputLP)
    {
        cerr << "LP/cycle lower bounds require LP reduction." << endl
             << flush;
        assert(0);
    }
    rootDepth = depth;
    if (outputLP)
    {
        if (REDUCTION < 0)
        {
            lpReduction();
        }
        else
        {
            reduce();
        }
        printf("%.1f\n", crt + rn / 2.0);
        return opt;
    }
    rec(t, time_limit);
    if (debug >= 2 && depth <= maxDepth)
        fprintf(stderr, "%sopt: %d\n", debugString().c_str(), opt);
    if (t.elapsed() >= time_limit)
        return -1;
    else
        return opt;
}

std::string branch_and_reduce_algorithm::debugString() const
{
    stringstream ins;
#ifdef PUT_TIME
    time_t rawtime;
    struct tm *timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    ins << std::put_time(timeinfo, "%T") << "  ";
#else
    std::locale::global(std::locale("ja_JP.utf8"));
    std::time_t t = std::time(NULL);
    char mbstr[100];
    if (std::strftime(mbstr, sizeof(mbstr), "%T", std::localtime(&t)))
    {
        std::cout << mbstr << '\n';
    }
#endif
    for (int i = 0; i < depth && i < maxDepth; ++i)
    {
        ins << " ";
    }
    return ins.str();
}

void branch_and_reduce_algorithm::PrintState() const
{
    cout << "State(" << this << "):" << endl
         << flush;
    cout << "adj=" << endl
         << flush;
    for (unsigned int j = 0; j < adj.size(); ++j)
    {
        cout << j << " : ";
        for (int const k : adj[j])
        {
            cout << k << " ";
        }
        cout << endl;
    }
    cout << "N  =" << N << endl
         << flush;
    cout << "in =";
    for (int const i : in)
    {
        cout << i << " ";
    }
    cout << endl
         << flush;
    cout << "out=";
    for (int const i : out)
    {
        cout << i << " ";
    }
    cout << endl
         << flush;
}

void branch_and_reduce_algorithm::reduce_graph()
{
    int oldn = rn;
    for (;;)
    {
        if (REDUCTION >= 0)
            deg1Reduction();
        if (REDUCTION >= 0 && REDUCTION < 2 && dominateReduction())
            continue;
        if (REDUCTION >= 2 && unconfinedReduction())
            continue;
        if (REDUCTION >= 1 && lpReduction())
            continue;
        if (REDUCTION >= 3)
        {
            int r = packingReduction();
            if (r < 0)
                return;
            if (r > 0)
                continue;
        }
        if (REDUCTION >= 1 && fold2Reduction())
            continue;
        if (REDUCTION >= 2 && twinReduction())
            continue;
        if (REDUCTION >= 2 && funnelReduction())
            continue;
        if (REDUCTION >= 2 && deskReduction())
            continue;
        break;
    }
    if (debug >= 2 && depth <= maxDepth && oldn != rn)
        fprintf(stderr, "%sreduce: %d -> %d\n", debugString().c_str(), oldn, rn);
    size_t low_degree_count(0);
    for (int v = 0; v < n; v++)
        if (x[v] < 0)
        {
            if (deg(v) <= 1)
            {
                low_degree_count++;
            }
        }
    cout << "There are " << low_degree_count << " degree 0 and 1 vertices left!" << endl
         << flush;
}

void branch_and_reduce_algorithm::initial_reduce_graph()
{
    reduce_graph();
    snapshotX = x;
    reductionSnapshotSize = modifieds.size();
}

size_t branch_and_reduce_algorithm::get_current_is_size() const
{

    std::vector<int> x2(x);

    for (int i = modifiedN - 1; i >= 0; i--)
    {
        modifieds[i]->reverse(x2);
    }

    size_t current_is_size(0);
    for (unsigned int i = 0; i < adj.size(); ++i)
    {
        if (x2[i] == 0)
        {
            current_is_size++;
        }
    }

    return current_is_size;
}

size_t branch_and_reduce_algorithm::get_current_is_size_with_folds() const
{

    size_t folded_vertex_count(0);
    size_t current_is_size(0);
    for (int const i : x)
    {
        if (i == 0)
            current_is_size++;
        if (i == 2)
            folded_vertex_count++;
    }

    return current_is_size + folded_vertex_count / 2;
}

bool branch_and_reduce_algorithm::folded_vertices_exist() const
{

    std::vector<int> x2(x);

    for (int i = modifiedN - 1; i >= 0; i--)
    {
        modifieds[i]->reverse(x2);
    }

    for (int const i : x2)
    {
        if (i == 2)
            return true;
    }

    return false;
}

std::vector<int> branch_and_reduce_algorithm::compute_maximal_is()
{

    int vertexToForceInIndependentSet(0);
    while (vertexToForceInIndependentSet != -1)
    {
        reduce_graph();

        vertexToForceInIndependentSet = -1;
        for (unsigned int i = 0; i < x.size(); ++i)
        {
            if (x[i] == -1)
            { // status not determined
                vertexToForceInIndependentSet = i;
                break;
            }
        }

        // add vertex to independent set
        if (vertexToForceInIndependentSet != -1)
        {
            set(vertexToForceInIndependentSet, 0);
        }
    }

    std::vector<int> x2(x);

    for (int i = modifiedN - 1; i >= 0; i--)
    {
        modifieds[i]->reverse(x2);
    }

    size_t current_is_size(0);
    for (unsigned int i = 0; i < adj.size(); ++i)
    {
        if (x2[i] == 0)
        {
            current_is_size++;
        }
    }

    return x2;
}

size_t branch_and_reduce_algorithm::compute_alternative_maximal_is_size()
{

    int vertexToForceInIndependentSet(0);
    while (vertexToForceInIndependentSet != -1)
    {
        reduce_graph();

        vertexToForceInIndependentSet = -1;
        for (int i = 0; i < static_cast<int>(x.size()); ++i)
        {
            if (x[i] == -1)
            { // status not determined
                vertexToForceInIndependentSet = i;
                break;
            }
        }

        // add vertex to independent set
        if (vertexToForceInIndependentSet != -1)
        {
            set(vertexToForceInIndependentSet, 0);
        }
    }

    size_t numberOfFoldedVertices(0);
    size_t sizeOfIS(0);
    for (int const i : x)
    {
        if (i == 0)
            sizeOfIS++;
        if (i == 2)
            numberOfFoldedVertices++;
    }

    return sizeOfIS + numberOfFoldedVertices / 2;
}

void branch_and_reduce_algorithm::restore_to_snapshot()
{

    while (modifiedN > reductionSnapshotSize)
    {
        modifieds[--modifiedN]->restore();
        modifieds[modifiedN] = shared_ptr<modified>();
    }

    x = snapshotX;
}

size_t branch_and_reduce_algorithm::number_of_nodes_remaining() const
{

    size_t node_count(0);
    for (int i : x)
        if (i == -1)
            node_count++;

    return node_count;
}

size_t branch_and_reduce_algorithm::number_of_edges_remaining() const
{

    size_t edge_count(0);
    for (int i = 0; i < adj.size(); i++)
    {
        if (x[i] == -1)
        {
            for (int j = 0; j < adj[i].size(); j++)
                if (x[adj[i][j]] == -1)
                    edge_count++;
        }
    }

    return edge_count / 2;
}

// void branch_and_reduce_algorithm::force_into_independent_set(std::vector<NodeID> const &nodes)
// {

//     for (NodeID const node : nodes)
//     {
//         assert(x[node] == -1); // should not have been assigned yet.
//         if (x[node] != -1)
//         {
//             cout << "ERROR: invalid vertex selected for independent set!" << endl
//                  << flush;
//         }
//         ////        cout << "Switching node " << node << " from value " << x[node] << " to 0" << endl;
//         set(node, 0);
//     }
// }

// input: std::vector mapping vertices in independent set to true

// WARNING: this is destructive, can only be applied once, and then the class is not
//          guaranteed to be a good state!
void branch_and_reduce_algorithm::extend_finer_is(std::vector<bool> &independent_set)
{

    assert(independent_set.size() == adj.size());
    assert(independent_set.size() == x.size());

    // first fill in temporary std::vector with independent set vertices from finer set
    size_t new_independent_set_vertices(0);
    for (size_t index = 0; index < independent_set.size(); ++index)
        if (independent_set[index])
        {
            assert(x[index] == -1); // it needs to be unassigned
            if (x[index] != -1)
            {
                cout << "ERROR: invalid vertex selected for independent set!" << endl
                     << flush;
            }
            set(index, 0); // add to independent set
            new_independent_set_vertices++;
        }

    std::vector<int> x2(x);

    // undo reductions
    for (int i = modifiedN - 1; i >= 0; i--)
    {
        modifieds[i]->reverse(x2);
    }

    // update full independent set
    for (unsigned int i = 0; i < adj.size(); ++i)
        if (x2[i] == 0)
        {
            independent_set[i] = true;
        }
}

void branch_and_reduce_algorithm::get_solved_is(std::vector<bool> &independent_set)
{
    for (unsigned int i = 0; i < y.size(); ++i)
    {
        if (y[i] == 0)
            independent_set[i] = true;
    }
}

// void branch_and_reduce_algorithm::convert_adj_lists(graph_access &G, std::vector<NodeID> &reverse_mapping) const
// {
//     // Number of nodes
//     unsigned int const node_count = number_of_nodes_remaining();
//     // Number of edges
//     int m = 0;

//     // Nodes -> Range
//     std::vector<NodeID> mapping(adj.size(), 10000000);
//     // Get number of edges and reorder nodes
//     unsigned int node_counter = 0;
//     for (NodeID node = 0; node < adj.size(); ++node)
//         if (x[node] < 0)
//         {
//             for (int const neighbor : adj[node])
//                 if (x[neighbor] < 0)
//                     m++;
//             mapping[node] = node_counter;
//             reverse_mapping[node_counter] = node;
//             node_counter++;
//         }

//     // Create the adjacency array
//     std::vector<int> xadj(node_count + 1);
//     std::vector<int> adjncy(m);
//     unsigned int adjncy_counter = 0;
//     for (unsigned int i = 0; i < node_count; ++i)
//     {
//         xadj[i] = adjncy_counter;
//         for (int const neighbor : adj[reverse_mapping[i]])
//         {
//             if (mapping[neighbor] == i)
//                 continue;
//             if (mapping[neighbor] == 10000000)
//                 continue;
//             adjncy[adjncy_counter++] = mapping[neighbor];
//         }
//         std::sort(std::begin(adjncy) + xadj[i], std::begin(adjncy) + adjncy_counter);
//     }
//     xadj[node_count] = adjncy_counter;

//     // Build the graph
//     G.build_from_metis(node_count, &xadj[0], &adjncy[0]);

// #if 0
//     std::vector<std::set<int>> neighbors;
//     neighbors.resize(G.number_of_nodes());
//     // verify the graph
//     forall_nodes(G, node) {
//         forall_out_edges(G, edge, node) {
//             NodeID neighbor = G.getEdgeTarget(edge);
//             neighbors[node].insert(neighbor);
//         } endfor
//     } endfor

//     for (int vertex = 0; vertex < node_count; ++vertex) {
//         size_t neighbors_in_subgraph(0);
//         for (int const neighbor : adj[reverse_mapping[vertex]]) {
//             bool const in_mapping(mapping[neighbor] != UINT_MAX);
//             if (in_mapping) {
//                 bool const in_graph(neighbors[vertex].find(mapping[neighbor]) != neighbors[vertex].end());
//                 if (in_graph) neighbors_in_subgraph++;
//             }
//         }

//         if (neighbors_in_subgraph != neighbors[vertex].size()) {
//             cout << "ERROR: subgraph verification failed" << endl << flush;
//         }
//     }
// #endif // DEBUG
// }

void branch_and_reduce_algorithm::convert_to_metis(int32_t *nNodes, std::vector<int32_t> &xadj, std::vector<int32_t> &adjncy, std::vector<int> &reverse_mapping)
{
    unsigned int const node_count = number_of_nodes_remaining();
    // Number of edges
    int m = 0;

    // Nodes -> Range
    std::vector<int> mapping(adj.size(), 10000000);
    reverse_mapping.resize(node_count, -1);
    // Get number of edges and reorder nodes
    unsigned int node_counter = 0;
    for (NodeID node = 0; node < adj.size(); ++node)
        if (x[node] < 0)
        {
            for (int const neighbor : adj[node])
                if (x[neighbor] < 0)
                    m++;
            mapping[node] = node_counter;
            reverse_mapping[node_counter] = node;
            node_counter++;
        }

    // Create the adjacency array
    xadj.resize(node_count + 1);
    adjncy.resize(m);
    unsigned int adjncy_counter = 0;
    for (unsigned int i = 0; i < node_count; ++i)
    {
        xadj[i] = adjncy_counter;
        for (int const neighbor : adj[reverse_mapping[i]])
        {
            if (mapping[neighbor] == i)
                continue;
            if (mapping[neighbor] == 10000000)
                continue;
            adjncy[adjncy_counter++] = mapping[neighbor];
        }
        std::sort(std::begin(adjncy) + xadj[i], std::begin(adjncy) + adjncy_counter);
    }
    xadj[node_count] = adjncy_counter;
    *nNodes = node_count;
}
/*
void branch_and_reduce_algorithm::convert_to_ga(std::shared_ptr<graph_access> G, std::vector<NodeID> &reverse_mapping, std::vector<NodeID> &mapping)
{
    // Number of nodes
    unsigned int const node_count = number_of_nodes_remaining();
    // Number of edges
    int m = 0;

    // Nodes -> Range
    mapping.resize(adj.size(), 10000000);

    // Get number of edges and reorder nodes
    unsigned int node_counter = 0;
    for (NodeID node = 0; node < adj.size(); ++node)
        if (x[node] < 0)
        {
            for (int const neighbor : adj[node])
                if (x[neighbor] < 0)
                    m++;
            mapping[node] = node_counter;
            reverse_mapping[node_counter] = node;
            node_counter++;
        }

    // Create the adjacency array
    std::vector<std::vector<int>> adja;
    adja.resize(node_counter);

    for (int i = 0; i < node_counter; i++)
        for (const int neighbor : adj[reverse_mapping[i]])
            if (x[neighbor] < 0)
                adja[i].emplace_back(mapping[neighbor]);

    // Build the graph
    G->build_from_adj(adja);
}
*/
void branch_and_reduce_algorithm::convert_to_adj(std::vector<std::vector<int>> &G, std::vector<int> &reverse_mapping, std::vector<int> &mapping)
{
    unsigned int const node_count = number_of_nodes_remaining();

    mapping.resize(adj.size(), -1);
    reverse_mapping.resize(node_count, -1);

    unsigned int node_counter = 0;
    for (int node = 0; node < adj.size(); ++node)
    {
        if (x[node] < 0)
        {
            mapping[node] = node_counter;
            reverse_mapping[node_counter] = node;
            node_counter++;
        }
    }
    G.resize(node_count);

    for (int node = 0; node < adj.size(); ++node)
    {
        if (x[node] < 0)
        {
            for (int neigh : adj[node])
            {
                if (x[neigh] < 0)
                    G[mapping[node]].push_back(mapping[neigh]);
            }
        }
    }
}

// NEW BRANCHING RULES

int branch_and_reduce_algorithm::get_articulation_point()
{
    get_articulation_points_iteratively();

    for (int i = 0; i < n; i++)
        if (articulation_points[i] == 1 && x[i] < 0)
            return i;

    return -1;
}

void branch_and_reduce_algorithm::get_articulation_points()
{
    current_dfs_num = 0;
    const int n = adj.size();

    visited.resize(0);
    minNr.resize(0);
    articulation_points.resize(0);

    visited.resize(n, -1);
    minNr.resize(n, -1);
    articulation_points.resize(n, 0);

    // Start DFS
    for (int i = 0; i < n; i++)
        if (x[i] < 0 && visited[i] < 0) // node active but not yet visited
            dfs_root(i);
}

void branch_and_reduce_algorithm::dfs_root(int v)
{
    minNr[v] = current_dfs_num;
    visited[v] = current_dfs_num++;

    int direct_children = 0;

    for (int u : adj[v])
    {
        if (x[u] < 0 && visited[u] < 0) // tree edge
        {
            direct_children++;

            dfs(u, v);
            minNr[v] = min(minNr[v], minNr[u]);

            if (minNr[u] >= visited[v]) // cut vertex found
                articulation_points[v] = 1;
        }
        // else if (x[u] < 0) // backwards edge
        //     minNr[v] = min(minNr[v], visited[u]);
    }

    if (direct_children < 2) // root is no cut vertex
        articulation_points[v] = 0;
}

void branch_and_reduce_algorithm::dfs(int v, int in)
{
    minNr[v] = visited[v] = current_dfs_num++;

    for (int u : adj[v])
    {
        if (x[u] < 0 && visited[u] < 0) // tree edge
        {
            dfs(u, v);
            minNr[v] = min(minNr[v], minNr[u]);

            if (minNr[u] >= visited[v]) // cut vertex found
                articulation_points[v] = 1;
        }
        else if (x[u] < 0 && u != in) // backwards edge
            minNr[v] = min(minNr[v], visited[u]);
    }
}

void branch_and_reduce_algorithm::get_articulation_points_iteratively() {
    current_dfs_num = 0;
    const int n = adj.size();

    visited = {};
    minNr = {};
    articulation_points = {};

    visited.resize(n, -1);
    minNr.resize(n, -1);
    articulation_points.resize(n, 0);

    for (int i = 0; i < n; i++)
        if (x[i] < 0 && visited[i] < 0) // node active but not yet visited
            dfs_iteratively(i);
}

void branch_and_reduce_algorithm::dfs_iteratively(int s) 
{
    dfs_stack.emplace(s, s);
    int child_cnt = -1;

    while (!dfs_stack.empty()) 
    {
        std::pair<int, int> e = dfs_stack.top();
        int v = e.first;
        int p = e.second;

        if (visited[v] < 0) // not yet visited
        {
            if (p == s) child_cnt++;

            visited[v] = current_dfs_num;
            minNr[v] = current_dfs_num++;

            for (int u : adj[v])
            {
                if (x[u] < 0 && visited[u] < 0) // tree edge
                    dfs_stack.emplace(u, v);
                else if (x[u] < 0 && u != p) // back edge
                    minNr[v] = min(minNr[v], visited[u]); 
            }
        }
        else
        {
            dfs_stack.pop();
            if (visited[v] <= n) // backtrack
            {
                visited[v] += (n + 1); // set v finished
                minNr[p] = min(minNr[p], minNr[v]);

                if (minNr[v] >= visited[p]) // articulation point found
                    articulation_points[p] = 1;
            }
        }     
    }

    if (child_cnt < 2) // root is no cut vertex
        articulation_points[s] = 0;  
    else articulation_points[s] = 1;
}

int branch_and_reduce_algorithm::get_mincut_vertex()
{
    // std::shared_ptr<graph_access> graph(new graph_access());
    // std::vector<NodeID> reverseMapping(number_of_nodes_remaining());
    // std::vector<NodeID> mapping(number_of_nodes_remaining());
    // convert_to_ga(graph, reverseMapping, mapping);

    // cut_algo.perform_minimum_cut(graph);

    // NodeID branchNode = -1;
    // int coutSize = 0;

    // for (int i = 0; i < graph->number_of_nodes(); i++)
    // {
    //     if (graph->getNodeInCut(i))
    //     {
    //         coutSize++;
    //         for (auto neig : graph->edges_of(i))
    //         {
    //             NodeID id = graph->getEdgeTarget(neig);
    //             if (!graph->getNodeInCut(id))
    //             {
    //                 branchNode = reverseMapping[id];
    //             }
    //         }
    //     }
    // }

    // return branchNode;
    return get_max_deg_vtx();
}

// void branch_and_reduce_algorithm::find_st_vtcs(std::shared_ptr<graph_access> graph, NodeID ss, NodeID tt)
// {
//     int v = -1;
//     deque<NodeID> queue;

//     if (ss == -1)
//     {
//         visited.resize(0);
//         visited.resize(graph->number_of_nodes(), 0);
//         queue.push_back(0);

//         while (!queue.empty())
//         {
//             NodeID node = queue[0];
//             queue.pop_front();
//             visited[node] = 1;
//             v = node;

//             for (EdgeID e : graph->edges_of(node))
//             {
//                 NodeID n_id = graph->getEdgeTarget(e);
//                 if (visited[n_id] == 0)
//                     queue.push_back(n_id);
//             }
//         }
//         s = v;
//     }
//     else
//         s = ss;
//     if (tt == -1)
//     {
//         visited.resize(0);
//         visited.resize(graph->number_of_nodes(), 0);
//         queue.clear();
//         queue.push_back(s);

//         while (!queue.empty())
//         {
//             NodeID node = queue[0];
//             queue.pop_front();
//             visited[node] = 1;
//             v = node;

//             for (EdgeID e : graph->edges_of(node))
//             {
//                 NodeID n_id = graph->getEdgeTarget(e);
//                 if (visited[n_id] == 0)
//                     queue.push_back(n_id);
//             }
//         }
//         t = v;
//     }
//     else
//         t = tt;
// }

void hc_karp_DFS(std::vector<std::vector<int>> & G, std::vector<int> &dist, std::vector<int> &matched, std::stack<int> & stack, int u, std::vector<int> &vis)
{
    vis[u] = 1;
    for (int t : G[u])
    {
        if (vis[t] == 0 && dist[t] == dist[u] + 1)
        {
            stack.push(t);
            vis[t] = 1;
            if (matched[t] == -1)
            {
                while (!stack.empty())
                {
                    int u = stack.top();
                    stack.pop();
                    int v = stack.top();
                    stack.pop();

                    matched[u] = v;
                    matched[v] = u;
                }
                return;
            }
            else
            {
                stack.push(matched[t]);
                hc_karp_DFS(G, dist, matched, stack, matched[t], vis);
            };
        }
    }
    if (!stack.empty())
    {
        stack.pop();
        stack.pop();
    }
}

void hc_karp(std::vector<std::vector<int>> & G, std::vector<int> & U, std::vector<int> & V, std::vector<int> &vc) 
{
    int n = U.size() + V.size();
    std::vector<int> matched(n, -1);

    bool aug = true;
    std::vector<int> dist(n, -1);

    // Repeat until no more augmenting paths
    while (aug)
    {
        aug = false;
        dist.clear();
        dist.resize(n, -1);

        std::queue<int> q;

        // push unmatched nodes of U into queue
        for (auto u : U)
        {
            if (matched[u] == -1)
            {
                q.push(u);
                dist[u] = 0;
            }
        }

        // BFS until augmenting path found
        while (!q.empty())
        {
            int v = q.front();
            q.pop();

            int dist_v = dist[v];

            for (int t : G[v])
            {
                if (dist[t] == -1)
                {
                    dist[t] = dist_v + 1;

                    if (matched[t] == -1)
                    {
                        aug = true;
                        break;
                    }
                    else
                    {
                        dist[matched[t]] = dist_v + 2;
                        q.push(matched[t]);
                    }
                }
            }
        }

        if (aug == false)
            break;

        // augment matching along path
        std::vector<int> vis(n, 0);
        for (auto u : U)
        {
            if (dist[u] == 0)
            {
                std::stack<int> stack;
                stack.push(u);
                hc_karp_DFS(G, dist, matched, stack, u, vis);
            }
        }
    }

    // get vc from matching
    for (int i = 0; i < matched.size(); i++)
    {
        if (std::find(U.begin(), U.end(), i) != U.end() && dist[i] == -1)
            vc.push_back(i);
        if (std::find(U.begin(), U.end(), i) == U.end() && dist[i] != -1)
            vc.push_back(i);
    }
}

void branch_and_reduce_algorithm::get_stcut_vertices()
{
    // choose s and t with max deg
    NodeID v1 = get_max_deg_vtx();
    x[v1] = 0;
    NodeID v2 = get_max_deg_vtx();
    x[v1] = -1;

    s = v1;
    t = v2;

    std::vector<int> res = {};
    std::vector<std::pair<int, int>> edges;

    push_relabel flow_algo(adj, x);
    int s2 = flow_algo.solve_max_flow_min_cut(rn, s, t, true, res, true, edges);

    std::vector<int> mapping(n, -1);
    std::vector<int> reverseMapping;
    std::vector<int> U, V;
    std::vector<std::vector<int>> biGraph;

    int id = 0;
    for (auto edge: edges)
    {
        int u = edge.first;
        int v = edge.second;        
        if (mapping[u] == -1)
        {
            mapping[u] = id;
            reverseMapping.push_back(u);
            U.push_back(id);
            biGraph.push_back(std::vector<int>());
            id++;
        }
        if (mapping[v] == -1)
        {
            mapping[v] = id;
            reverseMapping.push_back(v);
            V.push_back(id);
            biGraph.push_back(std::vector<int>());
            id++;
        }

        biGraph[mapping[u]].push_back(mapping[v]);
        biGraph[mapping[v]].push_back(mapping[u]);
    }

    std::vector<int> vc;
    hc_karp(biGraph, U, V, vc);
    for (int i = 0; i < vc.size(); i++)
    {
        vc[i] = reverseMapping[vc[i]];
    }

    cut.swap(vc);

    double perc = (double)res.size() / (double)number_of_nodes_remaining();
    if (cut.size() > TUNING_PARAM1 || perc < TUNING_PARAM2 || perc > (1.0 - TUNING_PARAM2))
    {
        // to big, use max. deg. vertex instead
        cut.clear();
        cut.emplace_back(get_max_deg_vtx());
        defaultBranch = true;
        defaultPicks++;
    }
    else
        stratPicks += cut.size();
}

int inline branch_and_reduce_algorithm::get_max_deg_vtx()
{
    int v, dv = -1;
    long long minE = 0;
    for (int u = 0; u < n; u++)
        if (x[u] < 0)
        {
            int degree = deg(u);
            if (dv > degree)
                continue;
            long long e = 0;
            used.clear();
            for (int w : adj[u])
                if (x[w] < 0)
                    used.add(w);
            for (int w : adj[u])
                if (x[w] < 0)
                {
                    for (int w2 : adj[w])
                        if (x[w2] < 0 && used.get(w2))
                            e++;
                }
            if (dv < degree || (dv == degree && minE > e))
            {
                dv = degree;
                minE = e;
                v = u;
            }
        }

    return v;
}

std::vector<std::vector<int>> branch_and_reduce_algorithm::get_nd_separators(int32_t *perm, int32_t *part_sizes, int32_t *sep_sizes, int n, int p, int32_t *weights)
{
    int level = log2(p);

    std::vector<std::vector<int>> seps(0);
    std::vector<int> t_sep(0);

    // extract top level separator
    int tsep_size = *(sep_sizes + p - 2);
    if (weights == NULL)
    {
        for (int i = n - tsep_size; i < n; i++)
        {
            t_sep.push_back(perm[i]);
        }
    }
    else
    {
        int pnt = n - 1;
        while (tsep_size > 0)
        {
            t_sep.push_back(perm[pnt]);
            tsep_size -= weights[perm[pnt]];
            pnt--;
        }
    }

    seps.push_back(t_sep);

    if (level == 0)
        return seps;

    // recursively extract separators of lower levels
    int lsize = 0;
    int rsize = 0;
    int32_t *lsize_arr = (int32_t *)malloc(sizeof(int32_t) * p / 2);
    int32_t *rsize_arr = (int32_t *)malloc(sizeof(int32_t) * p / 2);

    int32_t *l_ptr = lsize_arr;
    int32_t *r_ptr = rsize_arr;


    for (int i = level - 2; i >= 0; i--)
    {
        for (int j = 0; j < pow(2, i); j++)
        {
            lsize += *sep_sizes;
            *l_ptr = *sep_sizes;
            l_ptr++;
            sep_sizes++;
        }

        for (int j = 0; j < pow(2, i); j++)
        {
            rsize += *sep_sizes;
            *r_ptr = *sep_sizes;
            r_ptr++;
            sep_sizes++;
        }
    }

    for (int i = 0; i < p / 2; i++)
    {
        lsize += part_sizes[i];
        rsize += part_sizes[i + p / 2];
    }

    std::vector<std::vector<int>> seps_l = get_nd_separators(perm, part_sizes, lsize_arr, lsize, p / 2, weights);
    std::vector<std::vector<int>> seps_r = get_nd_separators(perm + lsize, part_sizes + p / 2, rsize_arr, rsize, p / 2, weights);

    for (int i = 0; i < seps_l.size(); i++)
    {
        seps.push_back(seps_l[i]);
        seps.push_back(seps_r[i]);
    }

    free(lsize_arr);
    free(rsize_arr);

    return seps;
}

void branch_and_reduce_algorithm::compute_nd_order()
{

    nd_threshold = std::max(30, (int)(0.1 * number_of_nodes_remaining()));

    int32_t n;
    std::vector<int32_t> xadj_v;
    std::vector<int32_t> adjncy_v;
    std::vector<int> rm;
    convert_to_metis(&n, xadj_v, adjncy_v, rm);

    int p = 4;

    if (TUNING_PARAM2 > 0)
        p = pow(2, TUNING_PARAM1);
    else
    {
        p = pow(2, ceil(log2(n)) - TUNING_PARAM1);
        if (p < 2)
            p = 2;
    }

    int32_t *xadj = (int32_t *)malloc(sizeof(int32_t) * xadj_v.size());
    int32_t *adjncy = (int32_t *)malloc(sizeof(int32_t) * adjncy_v.size());
    int32_t *perm = (int32_t *)malloc(sizeof(int32_t) * n);
    int32_t *iperm = (int32_t *)malloc(sizeof(int32_t) * n);
    int32_t *sizes = (int32_t *)malloc(sizeof(int32_t) * p * 2);

    //int32_t *degs = (int32_t *)malloc(sizeof(int32_t) * n);

    for (int i = 0; i < xadj_v.size(); i++)
        xadj[i] = xadj_v[i];

    for (int i = 0; i < adjncy_v.size(); i++)
        adjncy[i] = adjncy_v[i];

    int maxdeg = 0;
    // for (int i = 0; i < n; i++)
    // {
    //     if (deg(rm[i]) > maxdeg)
    //         maxdeg = deg(rm[i]);
    //     degs[i] = deg(rm[i]);
    // }

    // for (int i = 0; i < n; i++)
    // {
    //     degs[i] = maxdeg + 1 - degs[i];
    // }
    int r = METIS_NodeNDP(n, xadj, adjncy, NULL, p, NULL, perm, iperm, sizes);
    //int r = METIS_NodeNDP(n, xadj, adjncy, degs, p, NULL, perm, iperm, sizes);

    if (BRANCHING != 6666)
    {
        double thresh = 0;
        if (TUNING_PARAM3 > 0)
            thresh = TUNING_PARAM3;
        else if (TUNING_PARAM3 == 0)
            thresh = std::max(30, (int)ceil(log2(rn)));
        else
            thresh = std::max(30, (int)((-TUNING_PARAM3) * (rn / 100)));

        for (int i = p; i < 2 * p; i++)
        {
            if (sizes[i] > nd_threshold)
            {
                nd_computed = false;
                free(xadj);
                free(adjncy);
                free(perm);
                free(iperm);
                free(sizes);
                return;
            }
        }
    }
    std::vector<std::vector<int>> seps = get_nd_separators(perm, sizes, sizes + p, n, p, NULL);
    for (int i = seps.size() - 1; i >= 0; i--)
    {
        for (int j = 0; j < seps[i].size(); j++)
        {
            nd_order.push_back(rm[seps[i][j]]);
        }
    }
    //cout << "p: " << p << endl;
    //cout << *(sizes + 2 * p - 2) << endl;
    //cout << nd_order.size() << endl;

    free(xadj);
    free(adjncy);
    free(perm);
    free(iperm);
    free(sizes);

    //free(degs);
}

int branch_and_reduce_algorithm::max_nh_vtx()
{
    int v = -1;
    int maxN2 = 0;

    for (int i = 0; i < n; i++)
    {
        if (x[i] < 0)
        {
            used.clear();
            for (int u : adj[i])
            {
                if (x[u] < 0)
                    used.add(u);
            }

            int cnt = 0;
            for (int u : adj[i])
            {
                for (int v : adj[u])
                {
                    if (!used.get(v))
                    {
                        cnt++;
                        used.add(v);
                    }
                }
            }

            if (cnt > maxN2)
            {
                v = i;
                maxN2 = cnt;
            }
        }
    }

    return v;
}

void branch_and_reduce_algorithm::build_domination_graph()
{
    for (int i = 0; i < n; i++)
    {
        if (x[i] < 0)
            domination_graph[i].clear();
    }

    for (int v = 0; v < n; v++)
        if (x[v] < 0)
        {
            used.clear();
            used.add(v);
            for (int u : adj[v])
                if (x[u] < 0)
                    used.add(u);
            for (int u : adj[v])
                if (x[u] < 0)
                {
                    int cnt = 0;
                    int vtx = -1;
                    for (int w : adj[u])
                    {
                        if (x[w] < 0 && !used.get(w))
                        {
                            cnt++;
                            vtx = w;
                        }
                        if (cnt >= 2)
                            goto loop;
                    }
                    domination_graph[vtx].push_back(v);
                loop:;
                }
        }
}

void branch_and_reduce_algorithm::find_chains()
{
    for (int i = 0; i < n; i++)
    {
        if (x[i] >= 0)
            continue;

        chains[i].clear();
        int len = 1;

        std::vector<int> scanned(n, 0);
        std::queue<int> q;
        q.emplace(i);
        chains[i].push_back(i);
        scanned[i] = 1;

        while (!q.empty())
        {
            int v = q.front();
            q.pop();
            for (int u : domination_graph[v])
            {
                if (scanned[u] == 0)
                {
                    scanned[u] = 1;
                    q.push(u);
                    len++;

                    chains[i].push_back(u);
                }
            }
        }
        chainLength[i] = len;
    }
}

void branch_and_reduce_algorithm::calc_chain_vec()
{
    for (int i = 0; i < n; i++)
    {
        if (x[i] >= 0)
            chain_vec[i] = make_pair(0, 0);

        used.clear();
        int cnt = 0;
        for (int w : chains[i])
        {
            if (used.add(w))
                cnt++;
        }
        for (int u : adj[i])
        {
            if (x[u] < 0)
            {
                for (int w : chains[u])
                {
                    if (used.add(w))
                        cnt++;
                }
            }
        }
        chain_vec[i] = make_pair(chainLength[i], cnt);
    }
}
