/*!
A solver for linear Diophantine systems associated with AC and ACU matching. The exposition below is from Maude.


# Example

The system of linear (algebraic) equations

```
2x + 3y +  z + 5w + 6u + 0v == 26,
 x + 2y + 3z + 4w + 5u + 2v == 28,
3x + 2y + 5z +  w + 7u + 3v == 32,
2x +  y + 4z + 3w + 5u +  v == 25,
5x + 3y + 2z + 4w + 8u + 5v == 41,
 x + 4y + 2z +  w + 3u + 4v == 26
```
has solution x = 1, y = 2, z = 2, w = 2, u = 1, v = 2. So we have

```
 ┌─ ─┐T ┌─           ─┐   ┌────┐T
 │ 1 │  │ 2 1 3 2 5 1 │   │ 26 │
 │ 2 │  │ 3 2 2 1 3 4 │   │ 28 │
 │ 2 │  │ 1 3 5 4 2 2 │   │ 32 │
 │ 2 │  │ 5 4 1 3 4 1 │ = │ 25 │ ,
 │ 1 │  │ 6 5 7 5 8 3 │   │ 41 │
 │ 2 │  │ 0 2 3 1 5 4 │   │ 26 │
 └─ ─┘  └─           ─┘   └────┘

   R   *       M        =   C
```

We solve an alternative problem in which R and C are given and M is solved for. We constrain the matrix M by giving values `min_j` and `max_j` such that the sum of values in row `j` has minimum value `min_j` and maximum value `max_j`.

In general, there may be multiple solutions. To generate solutions, call `System.solve()` until it returns false. When it returns true, the solution is extracted with `System.solution(row, column)`.


```rust

let mut system = DiophantineSystem::new(6, 6);
system.insert_row(1, 10, 20);
system.insert_row(2, 11, 19);
system.insert_row(2, 15, 20);
system.insert_row(2, 15, 20);
system.insert_row(1, 30, 38);
system.insert_row(2, 12, 16);
system.insert_column(26);
system.insert_column(28);
system.insert_column(32);
system.insert_column(25);
system.insert_column(41);
system.insert_column(26);

while system.solve() {
  println!("\nSolution:");
  for row in 0..6 {
    for col in 0..6 {
      print!("{}  ", system.solution(row, col));
    }
    println!();
  }
}

```

# Description from Maude

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
*/

mod system;
pub(crate) mod row;

pub use system::DiophantineSystem;

// TODO: Templatize integer types.


#[derive(Copy, Clone, Default, Debug)]
pub(crate) struct Select {
  pub(crate) base      : u32,	// base value for element of $M$ (0 for simple systems)
  pub(crate) extra     : u32,	// extra value representing current state of solution
  pub(crate) max_extra : u32,	// maximum for extra
}


///	In a complex system, for each row with coefficient $R_i$ and for each possible
///	column value $V$ we compute and store the minimum and maximum $K$ such that
/// $V - K*R_i$ can be expressed as a natural number linear combination over
///	$R_j$ for $j > i$, respecting the maximum allowable sums but not the minimum
///	allowable sums (since some other column may make up the minimum).
///	If no such (natural number) $K$ exists we store `min = max = INSOLUBLE`.
#[derive(Copy, Clone, Default, Debug)]
pub(crate) struct Soluble {
  pub(crate) min: i32,	// minimum assignment to row for given column value
  pub(crate) max: i32,	// maximum assignment to row for given column value
}

impl Soluble {
  /// A special value used as a marker in the `Soluble` struct.
  pub(crate) const INSOLUBLE: i32 = -1;
  /// A special instance
  pub(crate) const INSOLUBLE_STRUCT: Soluble = Soluble{
    min: Soluble::INSOLUBLE,
    max: Soluble::INSOLUBLE,
  };
}

// Miscellaneous utility functions

#[inline(always)]
pub(crate) fn ceiling_division(dividend: i32, divisor: i32) -> i32 {
    if divisor > 0 {
      if dividend >= 0 {
          (dividend + divisor - 1) / divisor
      } else {
          -((-dividend)  / divisor)
      }
    }
    else {
      assert!(divisor < 0);
      if dividend >= 0 {
        -(dividend / (-divisor))
      } else {
        ((-dividend) + (-divisor) - 1) / (-divisor)
      }
    }
}

#[inline(always)]
pub(crate) fn floor_division(dividend: i32, divisor: i32) -> i32 {
    if divisor > 0 {
      if dividend >= 0 {
          dividend / divisor
      } else {
          -((divisor - dividend - 1) / divisor)
      }
    }
    else {
      assert!(divisor < 0);
      if dividend >= 0 {
        -((dividend - divisor - 1) / (-divisor))
      } else {
        (-dividend) / (-divisor)
      }
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn system_solver_test() {

      let mut system = DiophantineSystem::new(6, 6);
      system.insert_row(1, 14, 14);
      system.insert_row(2, 15, 15);
      system.insert_row(2, 17, 17);
      system.insert_row(2, 18, 18);
      system.insert_row(1, 34, 34);
      system.insert_row(2, 15, 15);
      system.insert_column(26);
      system.insert_column(28);
      system.insert_column(32);
      system.insert_column(25);
      system.insert_column(41);
      system.insert_column(26);

      // println!("Solve: {}", system.solve());
      let solution_count = 8;
      for _ in 0..solution_count{
      // let mut counter = 0usize;
      // while system.solve() {
        system.solve();
        println!("\nSolution:");
        for row in 0..6 {
          for col in 0..6 {
            print!("{}  ", system.solution(row, col));
          }
          println!();
        }
        // system.dump_info();


        // This program can computer ~1 million solutions per second on an M1 Max Macbook Pro.
        // counter +=1;
        // if counter % 1_000_000 == 0 {
        //   println!("{}", counter);
        // }
      }

      println!("Done!")

    }
}
