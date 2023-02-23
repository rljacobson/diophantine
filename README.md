# Linear Diophantine Equation Solver

This library implements the linear diophantine equation solver due to Steven Eker found in the Maude source code. It was
factored out of Mod, a RIIR of Maude. Note that this library currently only solves systems in which every number is
nonnegative. In fact, it assumes nonnegative values that fit in an `i32`. 

The next section is taken verbatim from comments in the Maude source.

# Algorithm Description

Based on:

> Steven Eker,
> "Single Elementary Associative-Commutative Matching",
> _Journal of Automated Reasoning_, pp35-51, 28(1), 2002.

Given an $n$-component vector of positive integers $R$ and an $m$-component vector
of positive integers $C$ a solution is an $n\times m$ matrix $M$ of natural numbers such
that
$R \cdot M = C$.
The intuition is that $M_{i,j}$ is the multiplicity of the $j$th constant assigned
to the $i$th variable in an single elementary AC or ACU matching problem.
For AC we are only interested in solutions which all rows (possibly all but
one in the case of extension variables) of $M$ have a non-zero sum.

We solve a slightly more general problem where minimum and maximum
values may be specified for the sum of each row of $M$.
The algorithm used here is optimized for moderately large problem instances.
We do not combine repeated entries in $R$ or make use of failure information
since these approaches to reducing the search space only pay for their overheads
for very large problem instances and would considerably complicate the code.

Rows are created by the member function
`insertRow(int coeff, int minSize, int maxSize)`
when `coeff` is the component of $R$ and `minSize` and `maxSize` are the
minimum and maximum values for the sum of the corresponding row of $M$.Columns are created by the member function
`insertColumn(int value)`
where `value` is the component of $C$.

To generate the first and successive solutions the member function `solve()` is
called; this returns `true` if a solution was generated and `false` if no more
solutions exist. If `true` was returned the components of $M$ can be extracted
using the member function `solution(int row, int column)`.

Note that it is an error to add more rows or columns after the first call to
`solve()`, or to try to examine a non-existent solution.

The basic approach is to sort $R$ into descending order and solve one row
at a time, backtracking whenever a dead end is detected. To solve a row
we consider $C$ as a multiset and compute a sub-multiset that we could
potentially use for the current row. We then systematically try selections
from this sub-multiset starting with the smallest ones.

The bottom level of the implementation uses two different algorithms for
solving a row depending on whether the system as a whole is "simple" or
"complex". A system is "*simple*" iff at least one element of $R$ is $1$ and
has a maximum allowable sum greater or equal that the largest column value.
Otherwise the system is "*complex*".
Since we sort $R$ into descending order, simple systems
have the property that any natural number can be expressed as a natural number
linear combination of any final segment of the sorted $R$. Thus we can rule
out one cause of failure for partial solutions. For complex systems we keep
a "solubility vector" which allows us to detect this kind of failure early
and prune the useless branches from the search.
