/*!

Structure representing the system of Diophantine equations.

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

We solve an alternative problem in which R and C are given and M is solved for. We constrain the matrix M by giving values min_j and max_j such that the sum of values in row j has minimum value min_j and maximum value max_j.

In general, there may be multiple solutions. To generate solutions, call `System.solve()` until it returns false. When it returns true, the solution is extracted with `System.solution(row, column)`.


```rust

let mut system = DiophantineSystem::new(6, 6);
system.insert_row(1, 10, 20); // 14 = actual sum of row
system.insert_row(2, 11, 19); // 15 = actual sum of row
system.insert_row(2, 15, 20); // 17 = actual sum of row
system.insert_row(2, 15, 20); // 18 = actual sum of row
system.insert_row(1, 30, 38); // 34 = actual sum of row
system.insert_row(2, 12, 16); // 15 = actual sum of row
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

*/

use std::cmp::{max, min};

use crate::{row::Row, ceiling_division, floor_division, Soluble, Select};

const UNBOUNDED: u32 = u32::MAX;

pub struct DiophantineSystem {
  rows        : Vec<Row>,
  columns     : Vec<u32>,
  row_permute : Vec<u32>,

  column_sum        : u32,
  max_column_value  : u32,
  closed            : bool, // System is closed once we start solving
  complex           : bool,
  failed            : bool  // Set when failure detected
}


impl DiophantineSystem {

  pub fn new(row_count: usize, col_count: usize) -> Self {
    DiophantineSystem {
      rows              : Vec::with_capacity(row_count),
      columns           : Vec::with_capacity(col_count),
      row_permute       : Vec::new(),
      column_sum        : 0,
      max_column_value  : 0,
      closed            : false,
      complex           : false,
      failed            : false,
    }
  }

  #[inline(always)]
  pub fn solution(&self, r: usize, c: usize) -> u32 {
    assert!(self.closed, "solve() not called");
    assert!(!self.failed, "non-existent soluiton");

    let s = &self.rows[self.row_permute[r] as usize].selection[c as usize];
    return s.base + s.extra;
  }

  #[inline(always)]
  pub fn row_count(&self) -> usize {
    self.rows.len()
  }

  #[inline(always)]
  pub fn column_count(&self) -> usize {
    self.columns.len()
  }

  pub fn insert_row(&mut self, coeff: u32, min_size: u32, max_size: u32) {
    assert!(!self.closed);
    assert!(coeff > 0);
    // assert!(min_size >= 0);
    assert!(min_size <= max_size);

    let row_count = self.rows.len();
    let new_row = Row{
      name: row_count as u32,
      coeff,
      min_size,
      max_size,
      ..Default::default()
    };

    self.rows.push(new_row);
  }

  pub fn insert_column(&mut self, value: u32) {
    assert!(value > 0);
    assert!(!self.closed);

    self.columns.push(value);
    self.column_sum += value;
    if value > self.max_column_value {
        self.max_column_value = value;
    }
  }
  // Check for trivial failure, sort R, fill out row_permute vector, compute
  // min_leave and max_leave values and allocate and initialize selection vectors.
  // For complex system we also build solubility vectors and check each compontent
  // of C for trivial failure.
  fn precompute(&mut self) -> bool {
    assert!(self.rows.len() > 0);
    assert!(self.columns.len() > 0);

    self.closed = true;

    #[cfg(feature = "dio_stats")]
    {
      for i in 0..self.rows.len() {
        let r = &self.rows[i];
        println!("row {} {} {} {}", i, r.coeff, r.min_size, r.max_size);
      }
      print!("columns: ");
      for i in 0..self.columns.len() {
        print!("{} ", self.columns[i]);
      }
      println!();
    }

    let mut sum_of_min_products = 0;
    let mut sum_of_max_products = 0;

    for r in self.rows.iter_mut() {
      if r.max_size == UNBOUNDED {
        r.max_size = self.column_sum;
      }
      r.min_product = r.min_size * r.coeff;
      sum_of_min_products += r.min_product;
      r.max_product = r.max_size * r.coeff;
      sum_of_max_products += r.max_product;
    }

    if sum_of_min_products > self.column_sum
        || sum_of_max_products < self.column_sum
    {
      self.failed = true;
      return false;
    }

    self.rows.sort();
    self.row_permute.resize(self.rows.len(), 0);

    let mut min_total: u32 = 0;
    let mut max_total: u32 = 0;
    for (i, row) in self.rows.iter_mut().enumerate().rev() {
      self.row_permute[row.name as usize] = i as u32;
      row.min_leave = min_total as i32;
      row.max_leave = max_total as i32;
      row.selection.resize(self.columns.len(), Select::default());

      min_total += row.min_product;
      max_total += row.max_product;
    }

    if self.rows.last().unwrap().coeff > 1
        || self.rows.last().unwrap().max_size < self.max_column_value
    {
      // The complex case
      self.build_solubility_vectors();
      let soluble = &mut self.rows[0].soluble;

      for column in self.columns.iter() {
        if soluble[*column as usize].min < 0 {
          self.failed = true;
          return false;
        }
      }

      self.complex = true;
    }

    true
  }


  // Function to build the solubility vectors discussed in [README.md] using a dynamic
  // programming approach.
  fn build_solubility_vectors(&mut self) {
    // Compute solubility vector for last row
    {
      let r         : &mut Row          = self.rows.last_mut().unwrap();
      let s         : &mut Vec<Soluble> = &mut r.soluble;
      let coeff     : u32               = r.coeff;
      let mut count : u32               = 0;

      s.resize(self.max_column_value as usize + 1, Soluble::INSOLUBLE_STRUCT);

      for j in (0..=self.max_column_value).step_by(coeff as usize) {
        s[j as usize].min = count as i32;
        s[j as usize].max = count as i32;
        count += 1;
        if count > r.max_size {
          break;
        }
      }
    }

    // Compute remaining vectors in descending order
    for i in (0..=(self.rows.len() - 2)).rev() {
      let max_size  : u32 = self.rows[i].max_size;
      let coeff     : u32 = self.rows[i].coeff;

      // Get mutable access to two elements at once.
      let (lower, upper) = self.rows.split_at_mut(i + 1);
      let next: &mut Vec<Soluble> = &mut lower.last_mut().unwrap().soluble; // self.rows[row_idx].soluble;
      let prev: &mut Vec<Soluble> = &mut upper.first_mut().unwrap().soluble; // self.rows[row_idx + 1].soluble;

      next.resize(self.max_column_value as usize + 1, Soluble::INSOLUBLE_STRUCT);

      for j in 0..=self.max_column_value as usize {
        if let Some(t) = j.checked_sub(coeff as usize) {
          if next[t].min != Soluble::INSOLUBLE && (max_size == UNBOUNDED || next[t].min < max_size as i32) {
            next[j].min = match prev[j].min {
              Soluble::INSOLUBLE => next[t].min + 1,
              _ => 0,
            };

            if max_size == UNBOUNDED || next[t].max < max_size as i32 {
              next[j].max = next[t].max + 1;
            }
            else {
              let mut new_max: i32 = max_size as i32;

              for k in ((j - ((max_size * coeff) as usize))..j).step_by(coeff as usize) {
                if prev[k].min == Soluble::INSOLUBLE {
                  new_max -= 1;
                } else {
                  break;
                }
              }

              assert!(new_max >= next[t].min + 1);
              next[j].max = new_max;
            }

          } else {
            let v = match prev[j].min {
              Soluble::INSOLUBLE => Soluble::INSOLUBLE,
              _ => 0,
            };

            next[j].min = v;
            next[j].max = v;
          }
        } else {
          let v = match prev[j].min {
            Soluble::INSOLUBLE => Soluble::INSOLUBLE,
            _ => 0,
          };

          next[j].min = v;
          next[j].max = v;
        }

      }
    }
  }



  pub fn solve(&mut self) -> bool {
    let find_first = !self.closed;
    if find_first && !self.precompute() {
      return false;
    }

    assert!(!self.failed);

    #[cfg(feature = "dio_stats")]
    {
      let r = if self.complex {
        self.solve_complex(find_first)
      } else {
        self.solve_simple(find_first)
      };
      if r {
        print!("success\t");
      } else {
        print!("failure");
      }
      return r;
    }

    if self.complex {
      self.solve_complex(find_first)
    } else {
      self.solve_simple(find_first)
    }
  }

  /// For each initial segment of the unsolved portion of R we check that there
  /// is a large enough sum of large enough elements in (what is left of) C to
  /// rule out a certain kind of failure. Return false if the current partial
  /// solution fails this test (and must therefore fail).
  #[inline]
  fn viable(&self, row_idx: usize) -> bool {
    let mut local_sum_of_min_products = 0;


    'okay:
    for row in self.rows[row_idx .. (self.rows.len() - 1)].iter() {
      let t = row.min_product;

      if t > 0 {
        local_sum_of_min_products += t;
        let lower_limit = row.coeff;
        let mut local_column_sum = 0;
        for c in self.columns.iter() {
          if *c >= lower_limit {
            local_column_sum += *c;
            if local_column_sum >= local_sum_of_min_products {
              continue 'okay;
            }
          }
        } // end iterate over cols
        return false;
      }

    } // end iterate over rows

    true
  }


  // region  The Simple Case

  /// Solve last row by allocating what is left.
  #[inline]
  fn solve_last_row_simple(&mut self) {
    let selection = &mut self.rows.last_mut().unwrap().selection;

    for i in 0..self.columns.len() {
      selection[i].extra = self.columns[i];
    }
  }


  /// Solve non-last row by trying to find a next selection for it and increasing
  // the size of selection we are looking for if necessary. If we are looking for
  // a first solution we first have to generate the multiset and determine the
  // feasable range of selection sizes.
  #[inline]
  fn solve_row_simple(&mut self, row_idx: usize, find_first: bool) -> bool {

    if find_first {
      if ! self.viable(row_idx) {
        return false;
      }
      let     r             : &mut Row = &mut self.rows[row_idx];
      let mut column_total  : u32      = 0;
      let mut max_sum       : u32      = 0;
      let     coeff         : u32      = r.coeff;

      for i in 0..self.columns.len() {
        r.selection[i].extra = 0;
        let mut t: u32       = self.columns[i];

        column_total += t;

        if t > coeff {
          t /= coeff;
          max_sum += t;
          r.selection[i].max_extra = t;
        }
        else {
          r.selection[i].max_extra = 0;
        }
      }

      let min_size: u32 = max(
        r.min_size,
        ceiling_division(
          (column_total as i32 - r.max_leave) as i32,
          coeff as i32
        ) as u32
      );
      let max_size: u32 = min(
        min(
          max_sum,
          r.max_size
        ),
        floor_division(
          (column_total as i32 - r.min_leave) as i32,
          coeff as i32
        ) as u32
      );

      if min_size > max_size {
        return false;
      }

      r.current_size = min_size;
      r.current_max_size = max_size;
    }
    else {
      let r: &mut Row = &mut self.rows[row_idx];

      if r.multiset_select(&mut self.columns, false) {
        return true;
      }
      else if r.current_size == r.current_max_size {
        return false;
      }

      r.current_size += 1;
    }

    // Always succeeds
    return self.rows[row_idx].multiset_select(&mut self.columns, true);
  }


  /// Solves the simple case using the auxiliary functions `solve_row_simple(..)` and `solve_last_row_simple(..)`.
  fn solve_simple(&mut self, mut find_first: bool) -> bool {
    if self.rows.len() > 1 {
      let penultimate_idx = self.rows.len() - 1;
      let mut i = if find_first { 0 } else { penultimate_idx };

      loop {
        find_first = self.solve_row_simple(i, find_first);
        if find_first {
          if i == penultimate_idx {
            break;
          }
          i += 1;
        }
        else {
          if i == 0 {
            break;
          }
          i -= 1;
        }
      }
    }

    if find_first {
      self.solve_last_row_simple();
    }
    else {
      self.failed = true;
    }

    find_first
  }
  // endregion


  // region The Complex Case

  #[inline]
  fn solve_last_row_complex(&mut self) {
    let last_row_idx  : usize                   = self.rows.len() - 1;
    let r             : &mut Row                = &mut self.rows[last_row_idx];
    let selection     : &mut Vec<crate::Select> = &mut r.selection;
    let soluble       : &Vec<Soluble>           = &r.soluble;
    let nr_columns    : usize                   = self.columns.len();

    for i in 0..nr_columns {
        let t = soluble[self.columns[i] as usize].min;
        assert!(t != Soluble::INSOLUBLE, "solubility bug");
        selection[i].extra = t as u32;
    }
  }


  fn solve_row_complex(&mut self, row_idx: usize, find_first: bool) -> bool {
    if find_first {
      if !self.viable(row_idx) {
        return false;
      }

      let     row          : &mut Row = &mut self.rows[row_idx];
      let     coeff        : u32      = row.coeff;
      let mut column_total : i32    = 0;
      let mut max_sum      : i32      = 0;
      let mut min_sum      : i32      = 0;

      for i in 0..self.columns.len() {
        let t   : usize = self.columns[i] as usize;
        let min : i32   = row.soluble[t].min;
        let max : i32   = row.soluble[t].max;
        assert!(min != Soluble::INSOLUBLE, "min Soluble::INSOLUBLE");
        assert!(max != Soluble::INSOLUBLE, "max Soluble::INSOLUBLE");
        assert!(min <= max, "min > max");

        row.selection[i].base      = min as u32;
        row.selection[i].extra     = 0;
        row.selection[i].max_extra = (max - min) as u32;

        column_total += t as i32;
        min_sum      += min;
        max_sum      += max;
      }

      let min_size = max(
        max(min_sum, row.min_size as i32),
        ceiling_division((column_total - row.max_leave) as i32, coeff as i32),
      );
      let max_size = min(
        min(max_sum, row.max_size as i32),
        floor_division((column_total - row.min_leave) as i32, coeff as i32),
      );

      if min_size > max_size {
        return false;
      }

      row.current_size     = (min_size - min_sum) as u32; // The maxes above gaurantee this is positive.
      row.current_max_size = (max_size - min_sum) as u32; // The mins  above gaurantee this is positive.

      for i in 0..self.columns.len() {
        if row.selection[i].base > 0 {
          self.columns[i] -= row.selection[i].base * coeff;
          // assert!(self.columns[i] >= 0, "value -ve");
        }
      }
    } //else

    // Get mutable access to two elements at once.
    let (lower, upper) = self.rows.split_at_mut(row_idx + 1);
    let row          : &mut Row          = lower.last_mut().unwrap();              // self.rows[row_idx];
    let coeff        : u32               = row.coeff;
    let next_soluble : &mut Vec<Soluble> = &mut upper.first_mut().unwrap().soluble; // self.rows[row_idx + 1].soluble;

    // This is an else for the previous if, but we want the bindings r and next_soluble in the outer scope.
    if !find_first {
      if row.multiset_complex(&mut self.columns, next_soluble, false) {
        return true;
      }

      row.current_size += 1;
    }

    while row.current_size <= row.current_max_size {
      if row.multiset_complex(&mut self.columns, next_soluble, true) {
        return true;
      }

      row.current_size += 1;
    }

    for i in 0..self.columns.len() {
      if (&mut row.selection)[i].base > 0 {
        self.columns[i] += row.selection[i].base * coeff;
        assert!(
          self.columns[i] <= self.max_column_value,
          "value too big"
        );
      }
    }

    false
  }




  fn solve_complex(&mut self, mut find_first: bool) -> bool {
    if self.rows.len() > 1 {
      let penultimate = self.rows.len() - 2;
      let mut i = if find_first { 0 } else { penultimate };
      loop {
        find_first = self.solve_row_complex(i, find_first);
        if find_first {
          if i == penultimate {
            break;
          }
          i += 1;
        } else {
          if i == 0 {
            break;
          }
          i -= 1;
        }
      }
    }
    if find_first {
      self.solve_last_row_complex();
    } else {
      self.failed = true;
    }
    find_first
  }

  // endregion

}


