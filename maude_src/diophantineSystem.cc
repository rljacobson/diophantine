/*

    This file is part of the Maude 3 interpreter.

    Copyright 1997-2003 SRI International, Menlo Park, CA 94025, USA.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.

*/

//
//      Implementation for class DiophantineSystem
//
#include "diophantineSystem.hh"

// #include "macros.hh"
// #include "vector.hh"

DiophantineSystem::DiophantineSystem(int estNrRows, int estNrColumns) {
    rows = std::vector<Row>();
    rows.reserve(estNrRows);
    columns = std::vector<int>();
    columns.reserve(estNrColumns);

    columnSum      = 0;
    maxColumnValue = 0;
    closed         = false;
    complex        = false;
    failed         = false;
}

void DiophantineSystem::insertRow(int coeff, int minSize, int maxSize) {
    Assert(!closed, "system closed");
    Assert(coeff > 0, "bad row coefficient " << coeff);
    Assert(minSize >= 0, "negative minSize " << minSize);
    Assert(minSize <= maxSize, "minSize > maxSize");
    int nrRows = rows.size();
    rows.resize(1 + nrRows);
    Row& r    = rows[nrRows];
    r.name    = nrRows;
    r.coeff   = coeff;
    r.minSize = minSize;
    r.maxSize = maxSize;
}

void DiophantineSystem::insertColumn(int value) {
    Assert(value > 0, "bad column value " << value);
    Assert(!closed, "system closed");
    columns.push_back(value);
    columnSum += value;
    if(value > maxColumnValue)
        maxColumnValue = value;
}

bool DiophantineSystem::solve() {
    bool findFirst = !closed;
    if(findFirst && !precompute())
        return false;
    Assert(!failed, "attempt to solve failed system");
#ifdef DIO_STATS
    bool r = complex ? solveComplex(findFirst) : solveSimple(findFirst);
    std::cout << (r ? "success\t" : "failure\n");
    return r;
#else
    return complex ? solveComplex(findFirst) : solveSimple(findFirst);
#endif
}

void DiophantineSystem::dumpInfo() const {
    std::cout << "row permute: [";
    for(auto n : rowPermute){
        std::cout << n << ", ";
    }
    std::cout << "]\nclosed: " << closed
        << "\ncolumn_sum: " << columnSum
        << "\ncomplex: " << complex
        << "\nfailed: " << failed
        << "\nmax_column_value: " << maxColumnValue
        << "\nRows:";

    for(auto row: rows){
        std::cout << "\n\tname: " << row.name
        << "\n\tcoeff: " << row.coeff
        << "\n\tminSize: " << row.minSize
        << "\n\tminProduct: " << row.minProduct
        << "\n\tminLeave: " << row.minLeave
        << "\n\tmaxSize: " << row.maxSize
        << "\n\tmaxProduct: " << row.maxProduct
        << "\n\tmaxLeave: " << row.maxLeave
        << "\n\tcurrentSize: " << row.currentSize
        << "\n\tcurrentMaxSize: " << row.currentMaxSize
        << "\n\tselection: [";
         for(auto sel: row.selection) {
             std::cout << "{" << sel.base << "} ";
             if (sel.extra != 0) {
                 std::cout << "+ " << sel.extra << " extra, " << sel.maxExtra << " maxExtra";
             }
             std::cout << ", ";
         }


        std::cout << "]\n\tsoluble: [";
        for(auto sol: row.soluble) {
             std::cout << "{" << sol.min << "," << sol.max << "} ";
        }
        std::cout << "]" << std::endl;
    }

    std::cout << "columns: [";
    for(auto n: columns) {
        std::cout << n << ", ";
    }
    std::cout << "]" << std::endl;

}

//
//	We sort rows in order of descending coefficients; splitting ties
//	in order of ascending maximum allowed sums
//
inline bool DiophantineSystem::rowLt(const Row& i, const Row& j) {
    int t = i.coeff - j.coeff;
    return (t != 0) ? t < 0 : (i.maxSize - j.maxSize) < 0;
}

//
//	Check for trivial failure, sort R, fill out rowPermute vector, compute
//	minLeave and maxLeave values and allocate and initialize selection vectors.
//	For complex system we also build solubility vectors and check each compontent
//	of C for trivial failure.
//
bool DiophantineSystem::precompute() {
    std::cout << "precompute" << std::endl;
    int nrRows = rows.size();
    Assert(nrRows > 0, "no rows");
    int nrColumns = columns.size();
    Assert(nrColumns > 0, "no columns");
    closed = true;

#ifdef DIO_STATS
    for(int i = 0; i < nrRows; i++) {
        Row& r = rows[i];
          std::cout << "row " << i << ' ' << r.coeff << ' ' << r.minSize << ' ' << r.maxSize << '\n';
    }
      std::cout << "columns: ";
    for(int i = 0; i < nrColumns; i++)
          std::cout << columns[i] << ' ';
      std::cout <<   std::endl;
#endif

    int sumOfMinProducts = 0;
    int sumOfMaxProducts = 0;
    for(int i = 0; i < nrRows; i++) {
        Row& r = rows[i];
        if(r.maxSize == UNBOUNDED)
            r.maxSize = columnSum;  // good substitute for infinity!
        r.minProduct = r.minSize * r.coeff;
        sumOfMinProducts += r.minProduct;
        r.maxProduct = r.maxSize * r.coeff;
        sumOfMaxProducts += r.maxProduct;
    }
    if(sumOfMinProducts > columnSum || sumOfMaxProducts < columnSum) {
        failed = true;
        std::cout << "Precompute failed #1\n";
        // dumpInfo();
        return false;
    }
    sort(rows.begin(), rows.end(), rowLt);
    rowPermute.resize(nrRows);
    int minTotal = 0;
    int maxTotal = 0;
    for(int i = nrRows - 1; i >= 0; i--) {
        Row& r             = rows[i];
        rowPermute[r.name] = i;
        r.minLeave         = minTotal;
        r.maxLeave         = maxTotal;
        Vector<Select>& s  = r.selection;
        s.resize(nrColumns);
        for(int j = 0; j < nrColumns; j++)
            s[j].base = 0;
        minTotal += r.minProduct;
        maxTotal += r.maxProduct;
    }
    if(rows[nrRows - 1].coeff > 1 || rows[nrRows - 1].maxSize < maxColumnValue) {
        //
        //	Complex case
        //
        buildSolubilityVectors();
        Vector<Soluble>& soluble = rows[0].soluble;
        for(int j = 0; j < nrColumns; j++) {
            if(soluble[columns[j]].min < 0) {
                failed = true;
                std::cout << "Precompute failed #2\n";
                // dumpInfo();
                return false;
            }
        }
        complex = true;
    }

    // dumpInfo();
    return true;
}

//
//	Function to build the solubility vectors discussed above using a dynamic
//	programming approach.
//
void DiophantineSystem::buildSolubilityVectors() {
    std::cout << "buildSolubilityVectors" << std::endl;

    int nrRows = rows.size();
    //
    //	Compute solubility vector for last row
    //
    {
        Row&             r = rows[nrRows - 1];
        Vector<Soluble>& s = r.soluble;
        s.resize(maxColumnValue + 1);
        int coeff = r.coeff;
        for(int j = 1; j <= maxColumnValue; j++)
            s[j].min = s[j].max = INSOLUBLE;

        int count = 0;
        for(int j = 0; j <= maxColumnValue && count <= r.maxSize; j += coeff)
            s[j].min = s[j].max = count++;
    }
    //
    //	Compute remaining vectors in descending order
    //
    for(int i = nrRows - 2; i >= 0; i--) {
        Vector<Soluble>& prev = rows[i + 1].soluble;
        Vector<Soluble>& next = rows[i].soluble;
        next.resize(maxColumnValue + 1);
        int coeff   = rows[i].coeff;
        int maxSize = rows[i].maxSize;
        for(int j = 0; j <= maxColumnValue; j++) {
            int t = j - coeff;
            if(t >= 0 && next[t].min != INSOLUBLE && next[t].min < maxSize) {
                next[j].min = (prev[j].min == INSOLUBLE) ? next[t].min + 1 : 0;
                if(next[t].max < maxSize)
                    next[j].max = next[t].max + 1;
                else {
                    int newMax = maxSize;
                    for(int k = j - maxSize * coeff; prev[k].min == INSOLUBLE; k += coeff)
                        --newMax;
                    Assert(newMax >= next[t].min + 1, "bad newMax");
                    next[j].max = newMax;
                }
            } else {
                next[j].min = next[j].max = ((prev[j].min == INSOLUBLE) ? INSOLUBLE : 0);
            }
        }
    }
}

//
//	For each initial segment of the unsolved portion of R we check that there
//	is a large enough sum of large enough elements in (what is left of) C to
//	rule out a certain kind of failure. Return false if the current partial
//	solution fails this test (and must therefore fail).
//
inline bool DiophantineSystem::viable(int row_idx) const {
    int nrRows                = rows.size();
    int nrColumns             = columns.size();
    int localSumOfMinProducts = 0;
    for(int i = row_idx; i <= nrRows - 2; i++)  // no need to consider last row
    {
        int t = rows[i].minProduct;
        if(t > 0) {
            localSumOfMinProducts += t;
            int lowerLimit     = rows[i].coeff;
            int localColumnSum = 0;
            for(int j = 0; j < nrColumns; j++) {
                int c = columns[j];
                if(c >= lowerLimit) {
                    localColumnSum += c;
                    if(localColumnSum >= localSumOfMinProducts)
                        goto okay;
                }
            }
            return false;
        }
    okay:;
    }
    return true;
}

/////////////////////////////////////////////////////////////////////////
//
//	Code for solving simple case
//
/////////////////////////////////////////////////////////////////////////

//
//	Find a selection from a multiset by undoing the previous selection until
//	the selected amount of some element can be increased by one (without
//	exceeding overall selection size). Then make up the size of the selection
//	by selecting the earliest elements available.
//
inline bool DiophantineSystem::Row::multisetSelect(Vector<int>& bag, bool findFirst) {
    std::cout << "Row::multisetSelect" << std::endl;
    int bagLength = bag.size();
    int undone;

    if(!findFirst) {
        if(currentSize > 0) {
            undone = 0;
            for(int j = 0; j < bagLength; j++) {
                Assert(selection[j].extra <= selection[j].maxExtra, "extra > maxExtra");
                int t = selection[j].extra;
                if(undone > 0 && t < selection[j].maxExtra) {
                    ++(selection[j].extra);
                    --undone;
                    bag[j] -= coeff;
                    goto forwards;
                }
                if(t > 0) {
                    selection[j].extra = 0;
                    undone += t;
                    bag[j] += t * coeff;
                }
            }
        }
        return false;
    } else {
        undone = currentSize;
    }

forwards:
    for(int j = 0; undone > 0; j++) {
        Assert(j < bagLength, "overran bag");
        int t = std::min(undone, selection[j].maxExtra);
        if(t > 0) {
            selection[j].extra = t;
            undone -= t;
            bag[j] -= t * coeff;
        }
    }
    return true;
}

//
//	Solve last row by allocating what is left.
//
inline void DiophantineSystem::solveLastRowSimple() {
    Vector<Select>& selection = rows[rows.size() - 1].selection;
    int             nrColumns = columns.size();
    for(int i = 0; i < nrColumns; i++)
        selection[i].extra = columns[i];
}

//
//	Solve non-last row by trying to find a next selection for it and increasing
//	the size of selection we are looking for if necessary. If we are looking for
//	a first solution we first have to generate the multiset and determine the
//	feasable range of selection sizes.
inline bool DiophantineSystem::solveRowSimple(int row_idx, bool findFirst) {
    std::cout << "solveRowSimple" << std::endl;
    Row& r = rows[row_idx];

    if(findFirst) {
        if(!viable(row_idx))
            return false;
        int columnTotal = 0;
        int maxSum      = 0;
        int coeff       = r.coeff;
        int nrColumns   = columns.size();
        for(int i = 0; i < nrColumns; i++) {
            r.selection[i].extra = 0;
            int t                = columns[i];
            columnTotal += t;
            if(t >= coeff) {
                t /= coeff;
                maxSum += t;
                r.selection[i].maxExtra = t;
            } else
                r.selection[i].maxExtra = 0;
        }
        int minSize = std::max(
            r.minSize,
            ceilingDivision(columnTotal - r.maxLeave, coeff)
            );
        int maxSize = std::min(
            std::min(maxSum, r.maxSize),
            floorDivision(columnTotal - r.minLeave, coeff)
            );
        if(minSize > maxSize)
            return false;
        r.currentSize    = minSize;
        r.currentMaxSize = maxSize;
    } else {
        if(r.multisetSelect(columns, false))
            return true;
        if(r.currentSize == r.currentMaxSize)
            return false;
        ++(r.currentSize);
    }
    return r.multisetSelect(columns, true);  // always succeeds
}

bool DiophantineSystem::solveSimple(bool findFirst) {
    std::cout << "solveSimple" << std::endl;
    if(rows.size() > 1) {
        int penultimate = rows.size() - 2;
        int i           = findFirst ? 0 : penultimate;
        for(;;) {
            findFirst = solveRowSimple(i, findFirst);
            if(findFirst) {
                if(i == penultimate)
                    break;
                ++i;
            } else {
                if(i == 0)
                    break;
                --i;
            }
        }
    }
    if(findFirst)
        solveLastRowSimple();
    else
        failed = true;
    return findFirst;
}

/////////////////////////////////////////////////////////////////////////
//
//	Code for solving complex case
//
/////////////////////////////////////////////////////////////////////////

//
//	Find a selection from a multiset by undoing the previous selection until
//	the selected amount of some element can be increased by one (without
//	exceeding overall selection size or violating solubility constraints).
//	Then make up the size of the selection by selecting the earliest elements
//	available (backtracking if this violates solubility constraints).
//
inline bool DiophantineSystem::Row::multisetComplex(Vector<int>& bag, Vector<Soluble>& soluble, bool findFirst) {
    std::cout << "Row::multisetComplex" << std::endl;
    int undone;

    int bagLength = bag.size();
    if(!findFirst) {
        if(currentSize > 0) {
            undone = 0;
        backtrack:
            for(int j = 0; j < bagLength; j++) {
                Assert(selection[j].extra <= selection[j].maxExtra, "extra > maxExtra");
                int t = selection[j].extra;
                if(undone > 0 && t < selection[j].maxExtra) {
                    int c = bag[j];
                    for(int e = 1; e <= undone; e++) {
                        Assert(t + e <= selection[j].maxExtra, "t + e > maxExtra");
                        c -= coeff;
                        if(soluble[c].min != INSOLUBLE) {
                            selection[j].extra = t + e;
                            bag[j]             = c;
                            undone -= e;
                            goto forwards;
                        }
                    }
                }
                if(t > 0) {
                    selection[j].extra = 0;
                    undone += t;
                    bag[j] += t * coeff;
                }
            }
        }
        return false;
    } else{
        undone = currentSize;
    }
forwards:
    for(int j = 0; undone > 0; j++) {
        Assert(j < bagLength, "overran bag");
        int t = selection[j].maxExtra;
        if(t <= undone) {
            if(t > 0) {
                selection[j].extra = t;
                undone -= t;
                bag[j] -= t * coeff;
            }
        } else {
            selection[j].extra = undone;
            bag[j] -= undone * coeff;
            undone = 0;
            if(soluble[bag[j]].min == INSOLUBLE)
                goto backtrack;
        }
    }
    return true;
}

//
//	Solve last row by allocating what is left divided by coefficient. We
//	actually avoid the explicit division by using the solubility vector.
//
inline void DiophantineSystem::solveLastRowComplex() {
    Vector<Select>&  selection = rows[rows.size() - 1].selection;
    Vector<Soluble>& soluble   = rows[rows.size() - 1].soluble;
    int              nrColumns = columns.size();
    for(int i = 0; i < nrColumns; i++) {
        int t = soluble[columns[i]].min;
        Assert(t != INSOLUBLE, "solubility bug");
        selection[i].extra = t;
    }
}

//
//	Solve non-last row by trying to find a next selection for it and increasing
//	the size of selection we are looking for if necessary. If we are looking for
//	a first solution we first have to generate the multiset and determine the
//	feasable range of selection sizes taking into account solubility constraints.
//	Because of solubility constaints we may end up with minimum values for
//	certain elements of M which are subtracted out of C just before searching
//	for the first solution and added back when the final failure occurs (just
//	before returning false).
//
inline bool DiophantineSystem::solveRowComplex(int rowNr, bool findFirst) {
    std::cout << "solveRowComplex" << std::endl;
    int              nrColumns = columns.size();
    Row&             r         = rows[rowNr];
    int              coeff     = r.coeff;
    Vector<Select>&  selection = r.selection;
    Vector<Soluble>& soluble2  = rows[rowNr + 1].soluble;

    if(findFirst) {
        if(!viable(rowNr))
            return false;
        Vector<Soluble>& soluble     = r.soluble;
        int              columnTotal = 0;
        int              maxSum      = 0;
        int              minSum      = 0;
        for(int i = 0; i < nrColumns; i++) {
            int t   = columns[i];
            int min = soluble[t].min;
            Assert(min != INSOLUBLE, "min insoluble");
            int max = soluble[t].max;
            Assert(max != INSOLUBLE, "max insoluble");
            Assert(min <= max, "min > max");

            selection[i].base     = min;
            selection[i].extra    = 0;
            selection[i].maxExtra = max - min;

            columnTotal += t;
            minSum += min;
            maxSum += max;
        }
        int minSize = std::max(std::max(minSum, r.minSize), ceilingDivision(columnTotal - r.maxLeave, coeff));
        int maxSize = std::min(std::min(maxSum, r.maxSize), floorDivision(columnTotal - r.minLeave, coeff));
        if(minSize > maxSize)
            return false;
        r.currentSize    = minSize - minSum;
        r.currentMaxSize = maxSize - minSum;
        for(int i = 0; i < nrColumns; i++) {
            if(selection[i].base > 0) {
                columns[i] -= selection[i].base * coeff;
                Assert(columns[i] >= 0, "value -ve");
            }
        }
    } else {
        if(r.multisetComplex(columns, soluble2, false))
            return true;
        ++(r.currentSize);
    }
    for(; r.currentSize <= r.currentMaxSize; r.currentSize++) {
        if(r.multisetComplex(columns, soluble2, true))
            return true;
    }
    for(int i = 0; i < nrColumns; i++) {
        if(selection[i].base > 0) {
            columns[i] += selection[i].base * coeff;
            Assert(columns[i] <= maxColumnValue, "value too big");
        }
    }
    return false;
}

bool DiophantineSystem::solveComplex(bool findFirst) {
    std::cout << "solveComplex" << std::endl;
    if(rows.size() > 1) {
        int penultimate = rows.size() - 2;
        int i           = findFirst ? 0 : penultimate;
        for(;;) {
            findFirst = solveRowComplex(i, findFirst);
            if(findFirst) {
                if(i == penultimate)
                    break;
                ++i;
            } else {
                if(i == 0)
                    break;
                --i;
            }
        }
    }
    if(findFirst)
        solveLastRowComplex();
    else
        failed = true;
    return findFirst;
}
