/*
#*****************************************************************************
#       Copyright (C) 2009 Carlo Hamalainen <carlo.hamalainen@gmail.com>, 
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
*/


/* To compile this file:

Visit http://www.boost.org/ and download the Boost library. Unpack
it somewhere.

Unpack nauty24b7 in the current directory, then:

$ cd nauty24b7
$ ./configure
$ make
$ cd ..
$ g++ -I /path_to_boost_1_37_0/ nauty24b7/nauty.o nauty24b7/nautil.o nauty24b7/naututil.o nauty24b7/rng.o nauty24b7/naugraph.o nauty24b7/naugroup.o klatin.cpp

*/

#include "mpicxx.h"

#include <string>
#include <cassert>
#include <iostream>
#include <map>
#include <vector>
#include <stack>
#include <utility>
#include <ctime>
#include <boost/multi_array.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_set.hpp>
#include <boost/random.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include "nauty24b7/nauty.h"
#include "nauty24b7/naugroup.h"
#include "nauty24b7/naututil.h"



#if 0
// from   http://snipplr.com/view/5907/stdrandomshuffle-and-boostrandomhpp/
class Random
{
public:
boost::mt19937 gen;
boost::uniform_int<int> dst;
boost::variate_generator< boost::mt19937, boost::uniform_int<int> > rand;
Random( int N ):// call instance:
gen( static_cast<unsigned long>(std::time(0)) ), dst( 0, N ), rand( gen, dst ) {
}
std::ptrdiff_t operator()( std::ptrdiff_t arg ) {
return static_cast< std::ptrdiff_t >( rand() );
}
};

Random rnd(4);
#endif

#define ESTIMATE_SIZE 0
boost::mt19937 rng;
long _C;
float _D;


// A k-latin square L of side n is an n by n array of multisets.
// For each cell (r,c), L[r][c] is a map from a symbol s to the number
// of times that it occurs in the cell (r,c). If a symbol s does not
// appear in the cell (r,c) then we assume that 
// L[r][c].find(s) == L[r][c].end()
typedef boost::multi_array<std::map<int,int> , 2> array_type;



int nr_times_in_cell(array_type &K, int r, int c, int s)
{
    std::map<int,int>::iterator iter = K[r][c].find(s);

    if (iter == K[r][c].end())  return 0;
    return iter->second;
}

bool has_symbol(array_type &K, int r, int c, int s)
{
    std::map<int,int>::iterator iter = K[r][c].find(s);

    return iter != K[r][c].end();
}

int unique_symbol(array_type &K, int r, int c)
{
    // Only makes sense when K is a 1-latin square.
    std::map<int,int>::iterator iter = K[r][c].begin();

    assert(iter != K[r][c].end());
    assert(iter->second == 1);

    return iter->first;
}

void add_symbol(array_type &K, int r, int c, int s)
{
    std::map<int,int>::iterator iter = K[r][c].find(s);

    if (iter == K[r][c].end())
        K[r][c].insert(std::make_pair(s, 1));
    else
        K[r][c][s] += 1;
}

void remove_symbol(array_type &K, int r, int c, int s)
{
    std::map<int,int>::iterator iter = K[r][c].find(s);

    assert(iter != K[r][c].end());

    if (iter->second == 1)
        K[r][c].erase(iter);
    else
        iter->second -= 1;
}

bool operator==(const std::vector<int> &A, const std::vector<int>&B)
{
    assert(A.size() == B.size());

    for(unsigned int i = 0; i < A.size(); ++i) {
        if (A[i] != B[i]) return false;
    }

    return true;
}

bool operator>(const std::vector<int> &A, const std::vector<int>&B)
{
    assert(A.size() == B.size());

    for(unsigned int i = 0; i < A.size(); ++i) {
        if (A[i] == B[i]) continue;

        if (A[i] > B[i]) return true;

        if (A[i] < B[i]) return false;
    }

    return false; // A == B
}


void print_klatin(array_type &K, int n)
{
    for(int r = 0; r < n; ++r) {
        for(int c = 0; c < n; ++c) {
            std::cout << "(";

            for(std::map<int,int>::iterator i = K[r][c].begin(); i != K[r][c].end(); ++i) {
                for(int w = 0; w < i->second; ++w) {
                    std::cout << i->first;
                }
            }
            std::cout << ")";
        }

        std::cout << std::endl;
    }
}

bool is_klatin(array_type &K, int n, int k)
{
    int nr_times;

    for(int s = 0; s < n; ++s) {
        // check that s appears k times in each row
        for(int r = 0; r < n; ++r) {
            nr_times = 0;

            for(int c = 0; c < n; ++c) {
                std::map<int,int>::iterator iter = K[r][c].find(s);

                if (iter != K[r][c].end()) {
                    nr_times += iter->second;
                }
            }

            if (nr_times != k) return false;
        }
    }

    for(int s = 0; s < n; ++s) {
        // check that s appears k times in each column
        for(int c = 0; c < n; ++c) {
            nr_times = 0;

            for(int r = 0; r < n; ++r) {
                std::map<int,int>::iterator iter = K[r][c].find(s);

                if (iter != K[r][c].end()) {
                    nr_times += iter->second;
                }
            }

            if (nr_times != k) return false;
        }
    }

    return true;
}

bool is_partial_klatin(array_type &K, int n, int k)
{
    int nr_times;

    for(int s = 0; s < n; ++s) {
        // check that s appears at most k times in each row
        for(int r = 0; r < n; ++r) {
            nr_times = 0;

            for(int c = 0; c < n; ++c) {
                std::map<int,int>::iterator iter = K[r][c].find(s);

                if (iter != K[r][c].end()) {
                    nr_times += iter->second;
                }
            }

            if (nr_times > k) return false;
        }
    }

    for(int s = 0; s < n; ++s) {
        // check that s appears at most k times in each column
        for(int c = 0; c < n; ++c) {
            nr_times = 0;

            for(int r = 0; r < n; ++r) {
                std::map<int,int>::iterator iter = K[r][c].find(s);

                if (iter != K[r][c].end()) {
                    nr_times += iter->second;
                }
            }

            if (nr_times > k) return false;
        }
    }

    return true;
}

array_type conjugate(int n, int k, array_type &K, int *p)
{
    array_type C;
    C.resize(boost::extents[n][n]);

    int x[3], y[3];

    for(int r = 0; r < n; r++) {
        for(int c = 0; c < n; c++) {
            for(std::map<int,int>::iterator i = K[r][c].begin(); i != K[r][c].end(); ++i) {
                int s = i->first;

                for(int w = 0; w < i->second; ++w) {

                    x[0] = r;
                    x[1] = c;
                    x[2] = s;

                    y[0] = x[p[0]];
                    y[1] = x[p[1]];
                    y[2] = x[p[2]];

                    add_symbol(C, y[0], y[1], y[2]);
                }
            }
        }
    }

    return C;
}

// conjugates of a k-latin square.
std::vector<array_type> conjugates(int n, int k, array_type &K)
{
    int c1[] = {0, 1, 2};
    int c2[] = {0, 2, 1};
    int c3[] = {1, 0, 2};
    int c4[] = {1, 2, 0};
    int c5[] = {2, 0, 1};
    int c6[] = {2, 1, 0};

    std::vector<array_type> conjs;

    array_type conj1 = conjugate(n, k, K, c1);
    array_type conj2 = conjugate(n, k, K, c2);
    array_type conj3 = conjugate(n, k, K, c3);
    array_type conj4 = conjugate(n, k, K, c4);
    array_type conj5 = conjugate(n, k, K, c5);
    array_type conj6 = conjugate(n, k, K, c6);

    conjs.push_back(conj1);
    conjs.push_back(conj2);
    conjs.push_back(conj3);
    conjs.push_back(conj4);
    conjs.push_back(conj5);
    conjs.push_back(conj6);

    return conjs;
}

// This function only checks that a *change* to a partial k-latin square
// preserves the k-latin property.
bool is_partial_klatin(array_type &K, int n, int k, int changed_row, int changed_column)
{
    int nr_times;

    for(int s = 0; s < n; ++s) {
        // check that s appears <= k times in the row
        int r = changed_row;
        nr_times = 0;

        for(int c = 0; c < n; ++c) {
            std::map<int,int>::iterator iter = K[r][c].find(s);

            if (iter != K[r][c].end()) {
                nr_times += iter->second;
            }
        }

        if (nr_times > k) return false;
    }

    for(int s = 0; s < n; ++s) {
        // check that s appears <= k times in the column
        int c = changed_column;
        nr_times = 0;

        for(int r = 0; r < n; ++r) {
            std::map<int,int>::iterator iter = K[r][c].find(s);

            if (iter != K[r][c].end()) {
                nr_times += iter->second;
            }
        }

        if (nr_times > k) return false;
    }

    return true;
}

int sum_of_image(std::map<int, int> m)
{
    int sum = 0;

    for(std::map<int,int>::iterator i = m.begin(); i != m.end(); ++i)
        sum += i->second;

    return sum;
}

// We construct each cell in nondecreasing order. So we require
// t <= s for each t in K[r][c].
bool is_nondecreasing_change(array_type &K, int n, int r, int c, int s)
{
    for(int t = s + 1; t < n; ++t) {
        std::map<int,int>::iterator iter = K[r][c].find(t);

        if (iter != K[r][c].end()) return false;
    }

    return true;
}

#define MODE_EXTEND 0
#define MODE_SOLUTION 1
#define MODE_BACKTRACK 2
#define MODE_DONE 3
#define MODE_ADVANCE 4

class klatin_dfs {
    int mode;
    int n, k;
    std::stack<boost::tuple<int, int, std::stack<int> > > choice;
    array_type P;
    bool parent;

public:
    array_type K;
    bool simple;
    int force_row;
    bool randomise_search;

    void set_parent(array_type _P) {
        P.resize(boost::extents[n][n]);
        P = _P;
        parent = true;

        assert(!simple); // fixme Not implemented

        assert(K[0][0].size() == 0); // we always start from an empty solution

        std::vector<int> possible_symbols;

        std::stack<int> tmp_stack;
        for(int s = 0; s < n; ++s) {
            if (!has_symbol(P, 0, 0, s)) continue;

            possible_symbols.push_back(s);
        }

        assert(!possible_symbols.empty());

        if (randomise_search)
            std::random_shuffle(possible_symbols.begin(), possible_symbols.end()); // , rnd);
        for(std::vector<int>::iterator iter = possible_symbols.begin(); iter != possible_symbols.end(); iter++)
            tmp_stack.push(*iter);

        // wipe out the default first choice.
        while (!choice.empty()) choice.pop();

        choice.push(boost::make_tuple(0, 0, tmp_stack));
    }

    void set_partial_solution(array_type _P) {
        // wipe out the default first choice.
        while (!choice.empty()) choice.pop();

        K = _P;

        boost::tuple<int, int> unfilled = unfilled_cell(K, n, k);
        int r = unfilled.get<0>();
        int c = unfilled.get<1>();

        assert(r >= 0 && c >= 0);

        std::stack<int> tmp_stack;
        for(int s = 0; s < n; ++s) {
            if (simple) {
                if (nr_times_in_cell(K, r, c, s) == 0)
                    tmp_stack.push(s); 
            } else {
                if (nr_times_in_cell(K, r, c, s) < k)
                    tmp_stack.push(s); 
            }
        }

        assert(!tmp_stack.empty());

        choice.push(boost::make_tuple(r, c, tmp_stack));
    }

    klatin_dfs(int _n, int _k) {
        simple = false;
        parent = false;
        randomise_search = false;

        n = _n;
        k = _k;

        K.resize(boost::extents[n][n]);

        for(int r = 0; r < n; ++r) {
            for(int c = 0; c < n; ++c) {
                K[r][c].clear();
            }
        }

        mode = MODE_ADVANCE;

        std::stack<int> tmp_stack;

        tmp_stack.push(0);
        choice.push(boost::make_tuple(0, 0, tmp_stack));

        force_row = -1;
    }

    boost::tuple<int, int> unfilled_cell(array_type &K, int n, int k)
    {
        int r, c;
        int r_to_fill = -1, c_to_fill = -1;

        if (force_row < 0)
            r = 0;
        else
            r = force_row;

        while (r < n) {
            if (force_row >= 0 && r != force_row) break;

            c = 0;
            while (c < n) {
                if (sum_of_image(K[r][c]) < k) {
                    r_to_fill = r;
                    c_to_fill = c;
                    break;
                }

                ++c;
            }

            if (r_to_fill >= 0) break;

            ++r;
        }

        return boost::make_tuple(r_to_fill, c_to_fill);
    }

    bool dfs_stack()
    {
        int r, c, s;
        std::stack<int> tmp_stack;
        std::vector<int> possible_symbols;
        boost::tuple<int, int, std::stack<int> > some_choice;

        boost::tuple<int, int> unfilled;
        
        while(true) {
            switch(mode) {
            case MODE_EXTEND:
                // if we get here with nothing in choice then we must have
                // run out of alternatives (assumes initial seeding from
                // outside the loop).
                if (choice.empty()) {
                    mode = MODE_DONE;
                    continue;
                }

                unfilled = unfilled_cell(K, n, k);
                r = unfilled.get<0>();
                c = unfilled.get<1>();

                if (r < 0 && c < 0) {
                    mode = MODE_SOLUTION;
                    continue;
                }

                // Make a stack of alternatives for this cell.
                while (!tmp_stack.empty()) tmp_stack.pop();

                possible_symbols.clear();
                for(int _s = 0; _s < n; ++_s) {
                    //K[r][c] = _s;

                    if (!is_nondecreasing_change(K, n, r, c, _s))
                        continue;

                    if (simple && has_symbol(K, r, c, _s))
                        continue;

                    if (parent && !has_symbol(P, r, c, _s))
                        continue;

                    // We can't put too many _s symbols into this cell
                    // if we are finding a sublatin square.
                    if (parent && nr_times_in_cell(K, r, c, _s) == nr_times_in_cell(P, r, c, _s)) continue;

                    add_symbol(K, r, c, _s);

                    if (is_partial_klatin(K, n, k, r, c))
                        possible_symbols.push_back(_s);

                    remove_symbol(K, r, c, _s);
                }

                if (randomise_search)
                    std::random_shuffle(possible_symbols.begin(), possible_symbols.end()); // , rnd);
                for(std::vector<int>::iterator iter = possible_symbols.begin(); iter != possible_symbols.end(); iter++)
                    tmp_stack.push(*iter);

                // If we can't find something to fill this cell then we have
                // reached an impossible situation and have to backtrack.
                if (tmp_stack.empty()) {
                    mode = MODE_BACKTRACK;
                    continue;
                }

                // Otherwise we have some alternatives and we can advance
                // the solution.
                choice.push(boost::make_tuple(r, c, tmp_stack));
                mode = MODE_ADVANCE;
                continue;

                break;
            case MODE_ADVANCE:
                assert(!choice.empty());

                some_choice = choice.top();

                r = some_choice.get<0>(); 
                c = some_choice.get<1>(); 
                tmp_stack = some_choice.get<2>(); 
                s = tmp_stack.top();

                if (simple) assert(!has_symbol(K, r, c, s));

                if (parent) assert(has_symbol(P, r, c, s));

                // We can't put too many _s symbols into this cell
                // if we are finding a sublatin square.
                if (parent && nr_times_in_cell(K, r, c, s) == nr_times_in_cell(P, r, c, s)) continue;

                if (is_nondecreasing_change(K, n, r, c, s)) {
                    add_symbol(K, r, c, s);

                    if (is_partial_klatin(K, n, k, r, c)) {
                        mode = MODE_EXTEND;
                        continue;
                    } else {
                        // not a partial k-latin square so we have to
                        // undo this change.
                        mode = MODE_BACKTRACK;
                        continue;
                    }
                }

                // otherwise we have tried to extend the cell in a
                // noncanonical way. Put this noncanonical change in and
                // then let MODE_BACKTRACK undo it.
                add_symbol(K, r, c, s);
                mode = MODE_BACKTRACK;
                continue;

                break;
            case MODE_SOLUTION:
                //assert(is_klatin(K, n, k));
                //print_klatin(K, n);
                //std::cout << std::endl;

                mode = MODE_BACKTRACK;
                //continue;
                return true;

                break;
            case MODE_BACKTRACK:
                // undo the current solution and try to extend in another
                // way.

                if (choice.empty()) {
                    mode = MODE_DONE;
                    continue;
                }

                some_choice = choice.top();
            
                r = some_choice.get<0>(); 
                c = some_choice.get<1>(); 
                tmp_stack = some_choice.get<2>(); 

                remove_symbol(K, r, c, tmp_stack.top());
               
                tmp_stack.pop();

                if (tmp_stack.empty()) {
                    // we ran out of choices for this cell so we need to
                    // try something else
                    choice.pop();
                    mode = MODE_BACKTRACK;
                    continue;
                } else {
                    // there's still something to try here so adjust the top
                    // of the mains stack and continue.
                    choice.pop();
                    choice.push(boost::make_tuple(r, c, tmp_stack));

                    mode = MODE_ADVANCE;
                    continue;
                }
     
                break;
            case MODE_DONE:
                return false;
                break;
            default:
                assert(false);
            }
        }
    }
};

array_type klatin_subtract(array_type &K, array_type &J, int n)
{
    // return K - J
    array_type S;
    S.resize(boost::extents[n][n]);

    for(int r = 0; r < n; ++r) {
        for(int c = 0; c < n; ++c) {
            for(int s = 0; s < n; ++s) {
                std::map<int,int>::iterator iter_K = K[r][c].find(s);
                std::map<int,int>::iterator iter_J = J[r][c].find(s);

                int a, b;

                if (iter_K != K[r][c].end())    a = iter_K->second;
                else                            a = 0;

                if (iter_J != J[r][c].end())    b = iter_J->second;
                else                            b = 0;

                assert(a >= b);

                if (a == 0) continue;

                if (a > b)
                    S[r][c].insert(std::make_pair(s, a-b));
            }
        }
    }

    return S;
}

bool has_no_klatin_subsquare(int n, int k, array_type &K, int k_1)
{
    assert(k_1 < k);
    assert(k > 0);

    klatin_dfs d(n, k_1);
    d.randomise_search = true;
    d.set_parent(K);

    return d.dfs_stack() == false;
}

bool is_separable(array_type K, int n, int k, int k_1)
{
    // Is the k-latin square K of side n separable into a
    // k_1 and (k-k_1)-latin square?

    assert(k_1 < k);
    assert(k > 0);

    klatin_dfs d(n, k_1);
    d.randomise_search = true;
    d.set_parent(K);

    while (d.dfs_stack()) {
        array_type KJ = klatin_subtract(K, d.K, n);

        if (is_klatin(KJ, n, k - k_1)) {
            return true;
        }
    }

    return false;
}


bool is_reducible(int n, int k, array_type &K)
{
    // Is the k-latin square K of side n reducible? For this we have to
    // find a sub-latin square J in K such that K-J is a (k-1)-latin
    // square.

    return is_separable(K, n, k, 1);
}

#define t_cell_vertex   0
#define t_row_vertex    1
#define t_col_vertex    2
#define t_sym_vertex    3
#define t_R_vertex      4
#define t_C_vertex      5
#define t_S_vertex      6

std::string type_to_str(int t)
{
    switch(t) {
    case t_cell_vertex: return "cell";
    case t_row_vertex:  return "row";
    case t_col_vertex:  return "col";
    case t_sym_vertex:  return "sym";
    case t_R_vertex:    return "R";
    case t_C_vertex:    return "C";
    case t_S_vertex:    return "S";
    }
}

// global, argh
// DYNALLSTAT(permutation,perm,perm_sz);

static void userautom(int count, permutation *perm, int *orbits, int numorbits, int stabvertex, int n)
{
    return;

	// afound ^= permhash(perm,n,107651L,count);
    for(int i = 0; i < n; i++)
        std::cout << perm[i] << " ";
    std::cout << std::endl;
}


class nauty_wrap
{
    std::vector<int> vertex_types;
    
    void nauty_init(int nr_rows, int nr_columns, int k_latin, array_type &L, int skip_row, bool allow_conjugates)
    {
        // L is a k-latin rectangle of size nr_rows by nr_columns.

        assert(nr_rows > 0);
        assert(nr_columns > 0);
        assert(k_latin > 0);

        assert(skip_row < 0 || (skip_row >= 0 && skip_row < nr_rows));

        DYNALLSTAT(graph,g,g_sz);
        DYNALLSTAT(graph,g_canonical,g_canonical_sz);
        DYNALLSTAT(int,lab,lab_sz);
        DYNALLSTAT(int,ptn,ptn_sz);
        DYNALLSTAT(int,orbits,orbits_sz);
        static DEFAULTOPTIONS_GRAPH(options);
        statsblk(stats);

        int n,m,v;
        set *gv;

        /* The following cause nauty to call two procedures which
        accumulate the group generators found by nauty. */
        // options.userautomproc = groupautomproc;

        // Our own procedure for the automorphism elements.
        options.userautomproc = userautom;

        options.userlevelproc = grouplevelproc;

        // we have k_latin*nr_rows*nr_columns vertices v_l which represent the
        // entries of the k-latin square; nr_rows vertices r_i, nr_columns vertices
        // c_j, nr_columns vertices s_k, plus special R, C, S.
        //
        // If skip_row >= 0 then we don't include that row in the canonical
        // form.
        if (skip_row < 0)
            n = k_latin*nr_rows*nr_columns + nr_rows + nr_columns + nr_columns + 3;
        else
            n = k_latin*(nr_rows-1)*nr_columns + (nr_rows-1) + nr_columns + nr_columns + 3;

        m = (n + WORDSIZE - 1) / WORDSIZE;
        nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

        vertex_types.resize(n);

        int worksize = 50*m; // recommended in the nauty manual
        setword workspace[worksize];

        DYNALLOC2(graph,g,g_sz,m,n,"malloc");
        DYNALLOC2(graph,g_canonical,g_canonical_sz,m,n,"malloc");
        DYNALLOC1(int,lab,lab_sz,n,"malloc");
        DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
        DYNALLOC1(int,orbits,orbits_sz,n,"malloc");

        for(v = 0; v < n; v++) {
            gv = GRAPHROW(g,v,m);
            EMPTYSET(gv,m);
        }

        int row_OFFSET, col_OFFSET, sym_OFFSET;
        int R_OFFSET, C_OFFSET, S_OFFSET;

        if (skip_row < 0) {
            row_OFFSET = k_latin*nr_rows*nr_columns;
            col_OFFSET = k_latin*nr_rows*nr_columns + nr_rows;
            sym_OFFSET = k_latin*nr_rows*nr_columns + nr_rows + nr_columns;

            R_OFFSET = k_latin*nr_rows*nr_columns + nr_rows + nr_columns + nr_columns;
            C_OFFSET = R_OFFSET + 1;
            S_OFFSET = C_OFFSET + 1;
        } else {
            row_OFFSET = k_latin*(nr_rows-1)*nr_columns;
            col_OFFSET = k_latin*(nr_rows-1)*nr_columns + (nr_rows-1);
            sym_OFFSET = k_latin*(nr_rows-1)*nr_columns + (nr_rows-1) + nr_columns;

            R_OFFSET = k_latin*(nr_rows-1)*nr_columns + (nr_rows-1) + nr_columns + nr_columns;
            C_OFFSET = R_OFFSET + 1;
            S_OFFSET = C_OFFSET + 1;
        }

        int row, col, sym;

        v = 0;
        for(int _row = 0; _row < nr_rows; _row++) { // _row in L
            if (skip_row >= 0 && _row == skip_row) continue;
            if (skip_row >= 0 && _row < skip_row) row = _row;
            if (skip_row >= 0 && _row > skip_row) row = _row - 1;
            if (skip_row < 0) row = _row; // no skipping

            // row is the 'nauty' variable. ugh.

            for(col = 0; col < nr_columns; col++) {
                for(std::map<int,int>::iterator iter = L[_row][col].begin(); iter != L[_row][col].end(); ++iter) {
                    sym = iter->first;

                    for(int w = 0; w < iter->second; ++w) {
                        // say that v is attached to row, col, and sym.
                        gv = GRAPHROW(g,v,m);
                        ADDELEMENT(gv, row_OFFSET + row);
                        ADDELEMENT(gv, col_OFFSET + col);
                        ADDELEMENT(gv, sym_OFFSET + sym);

                        // say that row is attached to v
                        gv = GRAPHROW(g, row_OFFSET + row, m);
                        ADDELEMENT(gv, v);
                        //n std::cout << row_OFFSET + row << " <-> " << v << std::endl;

                        // say that col is attached to v
                        gv = GRAPHROW(g, col_OFFSET + col, m);
                        ADDELEMENT(gv, v);
                        //n std::cout << col_OFFSET + col << " <-> " << v << std::endl;

                        // say that sym is attached to v
                        gv = GRAPHROW(g, sym_OFFSET + sym, m);
                        ADDELEMENT(gv, v);
                        //n std::cout << sym_OFFSET + sym << " <-> " << v << std::endl;

                        v++;
                    }
                }
            }
        }

        for(int _row = 0; _row < nr_rows; _row++) { // _row in L
            if (skip_row >= 0 && _row == skip_row) continue;
            if (skip_row >= 0 && _row < skip_row) row = _row;
            if (skip_row >= 0 && _row > skip_row) row = _row - 1;
            if (skip_row < 0) row = _row; // no skipping

            // R -> row
            gv = GRAPHROW(g, R_OFFSET, m);
            ADDELEMENT(gv, row_OFFSET + row);
            //n std::cout << R_OFFSET << " <-> " << row_OFFSET + row << std::endl;

            // R <- row
            gv = GRAPHROW(g, row_OFFSET + row, m);
            ADDELEMENT(gv, R_OFFSET);
        }

        for(int col = 0; col < nr_columns; col++) {
            // C -> col
            gv = GRAPHROW(g, C_OFFSET, m);
            ADDELEMENT(gv, col_OFFSET + col);
            //n std::cout << C_OFFSET << " <-> " << col_OFFSET + col << std::endl;

            // C <- col
            gv = GRAPHROW(g, col_OFFSET + col, m);
            ADDELEMENT(gv, C_OFFSET);
        }

        for(int sym = 0; sym < nr_columns; sym++) { // same number of columns as symbols!
            // S -> sym
            gv = GRAPHROW(g, S_OFFSET, m);
            ADDELEMENT(gv, sym_OFFSET + sym);
            //n std::cout << S_OFFSET << " <-> " << sym_OFFSET + sym << std::endl;

            // S <- sym
            gv = GRAPHROW(g, sym_OFFSET + sym, m);
            ADDELEMENT(gv, S_OFFSET);
        }

        // Now colour each class. The arrays lab and ptn hold the 
        // vertex labels and partition, respectively. A sequence of
        // vertices lab[i], lab[i+1], ..., lab[i+l] describes a contiguous
        // sequence of vertices, while ptn[i] = 1 unless i is the last vertex
        // in a colour class.

        int k = 0;

        // cell vertices
        if (skip_row < 0) {
            for(v = 0; v < k_latin*nr_rows*nr_columns; v++) {
                lab[k] = v;
                ptn[k] = 1;
                vertex_types[k] = t_cell_vertex;
                k++;
            }
            ptn[k-1] = 0; // marks the last vertex of this colour class

            // r_i
            for(v = row_OFFSET; v < row_OFFSET + nr_rows; v++) {
                lab[k] = v;
                ptn[k] = 1;
                vertex_types[k] = t_row_vertex;
                k++;
            }
        } else {
            for(v = 0; v < k_latin*(nr_rows-1)*nr_columns; v++) {
                lab[k] = v;
                ptn[k] = 1;
                vertex_types[k] = t_cell_vertex;
                k++;
            }
            ptn[k-1] = 0; // marks the last vertex of this colour class

            // r_i
            for(v = row_OFFSET; v < row_OFFSET + (nr_rows-1); v++) {
                lab[k] = v;
                ptn[k] = 1;
                vertex_types[k] = t_row_vertex;
                k++;
            }
        }
        
        // c_i
        for(v = col_OFFSET; v < col_OFFSET + nr_columns; v++) {
            lab[k] = v;
            ptn[k] = 1;
            vertex_types[k] = t_col_vertex;
            k++;
        }

        // s_i
        for(v = sym_OFFSET; v < sym_OFFSET + nr_columns; v++) {
            lab[k] = v;
            ptn[k] = 1;
            vertex_types[k] = t_sym_vertex;
            k++;
        }
        ptn[k-1] = 0;

        // colour R, C, S depending on whether we allow conjugaes
        lab[k] = R_OFFSET;
        ptn[k] = 1;
        vertex_types[k] = t_R_vertex;
        k++;
        if (!allow_conjugates)
            ptn[k-1] = 0;

        lab[k] = C_OFFSET;
        ptn[k] = 1;
        vertex_types[k] = t_C_vertex;
        k++;
        if (!allow_conjugates)
            ptn[k-1] = 0;

        lab[k] = S_OFFSET;
        ptn[k] = 1;
        vertex_types[k] = t_S_vertex;
        k++;
        ptn[k-1] = 0;

        assert(k == n);

        // Compute the canonical labelling.
        options.getcanon = TRUE;

        // We supplied an initial partition and labelling
        // so tell nauty not to compute its own.
        options.defaultptn = FALSE;

        nauty(  g,      // the graph
                lab,    // vertex labels
                ptn,    // colour partition
                NULL,   // active == NULL, unused
                orbits, // orbits of the automorphism group (I don't use these)
                &options,   // my options to nauty
                &stats,     // statistics about nauty's computation (I don't use these)
                workspace,  // nauty's workspace
                worksize,   // size of workspace
                m,      // number of setwords 
                n,      // number of vertices of the graph
                g_canonical // the canonically labelled graph
        );

        // std::cout << "grpsize1: " << stats.grpsize1 << " grpsize2: " << stats.grpsize2 << std::endl;

        // std::cout << "g:" << std::endl;
        // putgraph(stdout, g, options.linelength, m, n);
        // putptn(stdout, lab, ptn, 0, options.linelength, n);
        // putcanon(stdout, canonlab, graph *canong, int linelength, int m, int n)
        // std::cout << "g_canonical:" << std::endl;
        // putgraph(stdout, g_canonical, options.linelength, m, n);

        // Let's have a look at the orbits output variable.
        // According to the nauty manual, orbits[i] is the least-numbered
        // vertex in the same orbit as the vertex i, under the action of the
        // automorphism group of the graph.
        // std::cout << "orbits:" << std::endl;
        // for(int x = 0; x < n; x++) {
            // std::cout << "vertex " << x << " is in the orbit of vertex " << orbits[x] << std::endl;
        // }
        // std::cout << "/orbits" << std::endl;

        // The lab variable now tells us how to relabel the graph to get
        // it into canonical form.
        
        // print_klatin(L, nr_rows);
#if 0
        std::cout << "lab:" << std::endl;
        for(int x = 0; x < n; x++) {
            std::cout << "lab[" << x << "] =  " << lab[x] << ". ";
            std::cout << "type " << type_to_str(vertex_types[x]) 
                      << " to type " << type_to_str(vertex_types[lab[x]]) << std::endl;
        }
        std::cout << "/lab" << std::endl;

        int first_row = row_OFFSET, last_row = row_OFFSET + nr_rows - 1;

        std::cout << "skip_row = " << skip_row << std::endl;
        if (skip_row < 0) {
        } else {
            last_row--;
        }
        std::cout << "first row = " << first_row << std::endl;
        std::cout << "last row = " << last_row << std::endl;

        // find the minimum-labeled row under the canonical labelling.
        // fixme deal with the skip_row issue here. ugh.
        int min_row;
#endif

        // The following loop finds extracts the canonical label and 
        // is equivalent to the output of putcanon from nauty.
        //
        //std::cout << "putcanon" << std::endl;
        //putcanon(stdout, lab, g_canonical, options.linelength, m, n);
        //std::cout << "/putcanon" << std::endl;

        clabel.clear();

        for(v = 0; v < n; v++) {
            gv = GRAPHROW(g_canonical,v,m);
            for(int w = 0; w < n; w++) {
                if (ISELEMENT(gv, w))
                    clabel.push_back(w);
            }	
            clabel.push_back(-1);
        }

        // return clabel;
    }

public:
    // The canonical label of the latin square
    std::vector<int> clabel;
   
    nauty_wrap(int nr_rows, int nr_columns, int k_latin, array_type &L, bool allow_conjugates)
    {
        nauty_init(nr_rows, nr_columns, k_latin, L, -1, allow_conjugates);
    }

    nauty_wrap(int nr_rows, int nr_columns, int k_latin, array_type &L, int skip_row, bool allow_conjugates)
    {
        nauty_init(nr_rows, nr_columns, k_latin, L, skip_row, allow_conjugates);
    }
};


bool is_simple(int n, array_type &K)
{
    for(int r = 0; r < n; r++) {
        for(int c = 0; c < n; c++) {
            for(std::map<int,int>::iterator iter = K[r][c].begin(); iter != K[r][c].end(); iter++) {
                if (iter->second > 1) return false;
            }
        }
    } 

    return true;
}

void clear_square(int n, array_type &K)
{
    for(int row = 0; row < n; row++) {
        for(int col = 0; col < n; col++) {
            K[row][col].clear();
        }
    } 
}

bool is_canonical_augmentation_last(int nr_rows, int nr_columns, int k_latin, array_type &L,
                               std::vector<int> &clabel_parent, bool allow_conjugates)
{
    assert(nr_rows == nr_columns);
    std::vector<array_type> L_conjugates = conjugates(nr_rows, k_latin, L);

    std::vector<int> max_clabel;

    for(std::vector<array_type>::iterator iter = L_conjugates.begin(); iter != L_conjugates.end(); iter++) {
        for(int skip_row = 0; skip_row < nr_rows; skip_row++) {
            nauty_wrap nw(nr_rows, nr_columns, k_latin, *iter, skip_row, allow_conjugates);

            if (max_clabel.size() == 0) {
                max_clabel.resize(nw.clabel.size());
                std::copy(nw.clabel.begin(), nw.clabel.end(), max_clabel.begin());
            } else {
                if (nw.clabel > max_clabel) {
                    max_clabel.resize(nw.clabel.size());
                    std::copy(nw.clabel.begin(), nw.clabel.end(), max_clabel.begin());
                }
            }

        }
    }

    return clabel_parent == max_clabel;
}

bool is_canonical_augmentation(int nr_rows, int nr_columns, int k_latin, array_type &L,
                               std::vector<int> &clabel_parent, bool allow_conjugates)
{
    //std::cout << "clabel parent:" << std::endl;
    //std::copy(clabel_parent.begin(), clabel_parent.end(), std::ostream_iterator<int>(std::cout, " " ));
    //std::cout << std::endl;

    // compute the canonical form of each possible parent and choose the
    // minimal one.

    std::vector<int> max_clabel;
    for(int skip_row = 0; skip_row < nr_rows; skip_row++) {
        // std::vector<int> clabel = canonical_form(nr_rows, nr_columns, k_latin, L, skip_row, allow_conjugates);
        nauty_wrap nw(nr_rows, nr_columns, k_latin, L, skip_row, allow_conjugates);

        if (max_clabel.size() == 0) {
            max_clabel.resize(nw.clabel.size());
            std::copy(nw.clabel.begin(), nw.clabel.end(), max_clabel.begin());
        } else {
            if (nw.clabel > max_clabel) {
                max_clabel.resize(nw.clabel.size());
                std::copy(nw.clabel.begin(), nw.clabel.end(), max_clabel.begin());
            }
        }

    }

    return clabel_parent == max_clabel;
}

int stat_count, stat_nr_reducible, stat_nr_separable, stat_nr_simple;
int leaf_count;
bool first_run;
int save_size;
int rank;
std::vector<array_type> dump;

// boost::unordered_set<std::vector<int> > dump;

// find all canonical augmentations of L, a k-latin square of side n
// with the first nr_rows filled in.
void canonaug(int n, int k, int nr_rows, array_type &L, bool _simple, bool allow_conjugates)
{
    std::vector<array_type> isomorphs;
    boost::unordered_set<std::vector<int> > iso_labels;

    if (first_run && nr_rows == save_size) {
        dump.push_back(L);
        return;
    }

    if (nr_rows == n) {
        stat_count++;

        for(int k_1 = 1; k_1 <= k/2; k_1++) {
            if (is_separable(L, n, k, k_1)) {
                stat_nr_separable++;
                if (k_1 == 1) stat_nr_reducible++;
                break;
            }
        }

        if (is_simple(n, L)) stat_nr_simple++;

        if (ESTIMATE_SIZE) {
            // print_klatin(L, n); std::cout << _D << std::endl;
        }

        leaf_count++;

        return; 
    }

    klatin_dfs d(n, k);
    // currently L has rows 0, 1, ..., nr_rows-1 filled in, so 
    // fill in the next row.
    d.force_row = nr_rows;
    d.simple = _simple;
    d.set_partial_solution(L);

    nauty_wrap nw_parent(nr_rows, n, k, L, allow_conjugates);

    while(d.dfs_stack()) {
        bool iscanon;

        if (nr_rows + 1 == n)
            iscanon = is_canonical_augmentation_last(nr_rows + 1, n, k, d.K, nw_parent.clabel, allow_conjugates);
        else 
            iscanon = is_canonical_augmentation(nr_rows + 1, n, k, d.K, nw_parent.clabel, allow_conjugates);

        if (iscanon) {
            nauty_wrap nw(nr_rows + 1, n, k, d.K, allow_conjugates);

            if (iso_labels.find(nw.clabel) == iso_labels.end()) {
                iso_labels.insert(nw.clabel);
                isomorphs.push_back(d.K);
            }
        }
    }

    if (ESTIMATE_SIZE && isomorphs.empty()) {
        return;
    }

    if (ESTIMATE_SIZE) {
        std::cout << _D << " " << isomorphs.size() << std::endl;

        boost::uniform_int<> six(0, isomorphs.size() - 1);
        boost::variate_generator<boost::mt19937&, boost::uniform_int<> > die(rng, six);

        _D *= isomorphs.size();

        // std::cout << die() << " :: " << isomorphs.size() << std::endl;

        canonaug(n, k, nr_rows + 1, isomorphs[die()], _simple, allow_conjugates); 
    } else {
        if (isomorphs.empty()) leaf_count++;

        for(std::vector<array_type>::iterator iter = isomorphs.begin(); iter != isomorphs.end(); iter++) {
            canonaug(n, k, nr_rows + 1, *iter, _simple, allow_conjugates); 
        }
    }
}

// std::map<std::vector<int>, array_type> single_rows(int n, int k, bool _simple, bool allow_conjugates)
std::vector<array_type> single_rows(int n, int k, bool _simple, bool allow_conjugates)
{
    boost::unordered_set<std::vector<int> > isomorph_labels;
    std::vector<array_type> isomorphs;

    klatin_dfs d(n, k);
    d.force_row = 0;
    d.simple = _simple;

    while(d.dfs_stack()) {
        // std::vector<int> clabel = canonical_form(1, n, k, d.K, allow_conjugates);
        nauty_wrap nw(1, n, k, d.K, allow_conjugates);

        if (isomorph_labels.find(nw.clabel) == isomorph_labels.end()) {
            // isomorphs.insert(std::make_pair(nw.clabel, d.K));
            isomorph_labels.insert(nw.clabel);
            isomorphs.push_back(d.K);
        }
    }

    return isomorphs;
}

void find_nonseparable(int n, int k)
{
    klatin_dfs d(n, k);
    d.randomise_search = true;
    while(d.dfs_stack()) {
        bool has_separation = false;

        for(int k_1 = 1; k_1 <= k/2; k_1++) {
            if (is_separable(d.K, n, k, k_1)) {
                has_separation = true;
                break;
            }
        }

        if (!has_separation) {
            print_klatin(d.K, n);
            for(int k_1 = 1; k_1 <= k/2; k_1++) {
                std::cout << "has_no_klatin_subsquare k_1 = " << k_1 << " result: " 
                << has_no_klatin_subsquare(n, k, d.K, k_1) << std::endl;;
            }
            exit(0);
        }
    }
}

// fixme: run a test on this for Nick

/*
I have conjectured that for order 3, a k-latin square is always
separable for k\geq 3, but it would be good to check this by computer
for some small values of k.
*/

int main(int argc, char **argv)
{
    MPI::Init(argc, argv);

    rank = MPI::COMM_WORLD.Get_rank();
    int size = MPI::COMM_WORLD.Get_size();

    if (argc != 5) {
        printf("Usage: ./klatin seed n save_size k\n");
        printf("\n");
        printf("seed:       random seed\n");
        printf("n:          order of latin squares to find\n");
        printf("save_size:  compute up to save_size rows before distributing to other nodes\n");
        printf("k:          find k-latin squares of order n\n");
        printf("\n");
        printf("example:    ./klatin 0 3 2 4\n");
        printf("            n k count nr_erodible nr_separable nr_simple \n");
        printf("            3 4 24 22 24 0\n");

        return 1;
    }

    rng.seed(atoi(argv[1]));

	int n = atoi(argv[2]);
	save_size = atoi(argv[3]);
	int k = atoi(argv[4]);

    assert(save_size >= 0);
    assert(ESTIMATE_SIZE || (!ESTIMATE_SIZE && save_size < n));

    bool find_only_simple = false;
    bool allow_conjugates = true;

    std::vector<array_type> singles = single_rows(n, k, find_only_simple, allow_conjugates);

    stat_count = 0;
    stat_nr_reducible = 0;
    stat_nr_separable = 0;
    stat_nr_simple = 0;

    int sum_D = 0;
    int nr_runs = 5;
    int nr_done = 0;
lazy:

    // Only use one processor to estimate the size of the search tree.
    if (ESTIMATE_SIZE && rank == 0) {
        save_size = n; // wipe over our command line parameter

        _D = singles.size();

        boost::uniform_int<> six(0, singles.size() - 1);
        boost::variate_generator<boost::mt19937&, boost::uniform_int<> > die(rng, six);

        canonaug(n, k, 1, singles[die()], find_only_simple, allow_conjugates); 

        // if (_D == -1) goto lazy;

        sum_D += _D;
        nr_done++;
        if (nr_done == nr_runs) {
            std::cout << "average: " << sum_D/(1.0*nr_runs) << std::endl;
        } else {
            goto lazy;
        }
    } else {
        leaf_count = 0;
        first_run = true;
        for(std::vector<array_type>::iterator iter = singles.begin(); iter != singles.end(); iter++) {
            canonaug(n, k, 1, *iter, find_only_simple, allow_conjugates); 
        }

        first_run = false;
        const int slice_size = dump.size() / size;

        assert(slice_size > 0);

        if (rank < size - 1) {
            for(int i = 0; i < slice_size; i++) {
                canonaug(n, k, save_size, dump[slice_size*rank + i], find_only_simple, allow_conjugates); 
            }
        } else {
            // last CPU has slightly more work if size does not divide
            // bitrade_output.size() evenly.
            for(int i = slice_size*(size-1); i < dump.size(); i++) {
                canonaug(n, k, save_size, dump[i], find_only_simple, allow_conjugates); 
            }
        }

        // stat_count, stat_nr_reducible, stat_nr_separable, stat_nr_simple;
        int sendcount = 4;
        int sendbuf[4];
        int recvbuf[size*4];
        int recvcounts[size];
        int displs[size];

        for (int i = 0; i < size; i++) {
            displs[i] = i*sendcount;
            recvcounts[i] = sendcount;
        }

        sendbuf[0] = stat_count;
        sendbuf[1] = stat_nr_reducible;
        sendbuf[2] = stat_nr_separable;
        sendbuf[3] = stat_nr_simple;

        MPI::COMM_WORLD.Gatherv(sendbuf, sendcount, MPI::INT, recvbuf, recvcounts, displs, MPI::INT, 0);

        if (rank == 0) {
            for(int u = 1; u < size; u++) {
                for(int i = 0; i < sendcount; i++) {
                    recvbuf[i] += recvbuf[u*sendcount + i];
                }
            }

            std::cout << "n k count nr_erodible nr_separable nr_simple " << std::endl;
            std::cout << n << " " << k << " " << recvbuf[0] << " " << recvbuf[1] << " "
                         << recvbuf[2] << " " << recvbuf[3] << std::endl;

            if (ESTIMATE_SIZE)
                std::cout << "leaf count: " << leaf_count << std::endl;
        }
    }

    MPI::Finalize();

    return 0;

#if 0
    int n = 4;
    int k = 4;

    klatin_dfs d(n, k);

    while(d.dfs_stack()) {
        if (is_separable(d.K, n, k, 2) && !is_reducible(d.K, n, k)) {
            std::cout << "klatin:" << std::endl;
            print_klatin(d.K, n);

            std::cout << "separable: " << is_separable(d.K, n, k, 2) << std::endl;
            std::cout << "reducible: " << is_reducible(d.K, n, k) << std::endl;
            std::cout << std::endl;
            break;
        }
    }

    //std::cout << std::endl;
    //std::cout << std::endl;
    //std::cout << std::endl;
    //std::cout << std::endl;

    //std::cout << "klatin:" << std::endl;
    //print_klatin(d.K, n);

    //is_reducible(d.K, n, k);
#endif


#if 0
    int n = 4;
    int k = 2;

    klatin_dfs d(n, k);

    d.force_row = 0;

    assert(d.dfs_stack());

    klatin_dfs e(n, k);
    e.force_row = 1;
    e.set_partial_solution(d.K);
    assert(e.dfs_stack());

    klatin_dfs f(n, k);
    f.force_row = 2;
    f.set_partial_solution(e.K);
    assert(f.dfs_stack());

    klatin_dfs g(n, k);
    g.force_row = 3;
    g.set_partial_solution(f.K);
    assert(g.dfs_stack());

    print_klatin(g.K, n);
    assert(is_klatin(g.K, n, k));


    return 0;

    std::vector<int> cform1 = canonical_form(n, n, k, d.K, n - 1);
    std::vector<int> cform2 = canonical_form(n - 1, n, k, d.K, -1);

    assert(vec_eq(cform1, cform2));


    //cform = canonical_form(n - 1, n, k, d.K, -1);
    //std::copy(cform.begin(), cform.end(), std::ostream_iterator<int>(std::cout, " " ));
    //std::cout << std::endl;
 
    return 0;
#endif
}





