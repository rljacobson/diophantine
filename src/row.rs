/*!

A row of a system of linear Diophantine equations.

*/

use std::cmp::{min, Ordering};

use crate::{Select, Soluble};

///	Structure for each row. We have a pair of member functions to handle
///	making a selection from a multiset, both normally and in the presence
///	of solubility constraints on the non-selected part.
#[derive(Default)]
pub(crate) struct Row {
    pub(crate) name             : u32,          // original position of row
    pub(crate) coeff            : u32,          // coefficient
    pub(crate) min_size         : u32,          // minimum acceptable sum
    pub(crate) min_product      : u32,          // coeff * minSize
    pub(crate) min_leave        : i32,          // minimum sum that must be left for
                                                // remaining rows
    pub(crate) max_size         : u32,          // maximum acceptable sum
    pub(crate) max_product      : u32,          // coeff * maxSize
    pub(crate) max_leave        : i32,          // maximum sum that may be left for
                                                // remaining rows
    pub(crate) current_size     : u32,          // current size of selection from multiset
    pub(crate) current_max_size : u32,          // maximum size of selection from multiset
    pub(crate) selection        : Vec<Select>,  // vector of values selected for this row
    pub(crate) soluble          : Vec<Soluble>, // solubility vector (complex systems only)
}

impl Row {
    ///	Find a selection from a multiset by undoing the previous selection until
    ///	the selected amount of some element can be increased by one (without
    ///	exceeding overall selection size). Then make up the size of the selection
    ///	by selecting the earliest elements available.
    pub fn multiset_select(&mut self, bag: &mut Vec<u32>, find_first: bool) -> bool {
        let mut undone: u32 = 0;
        let mut forwards: bool = false; // A flag directing control flow.

        if !find_first {
            if self.current_size > 0 {
                undone = 0;

                for j in 0..bag.len() {
                    assert!(self.selection[j].extra <= self.selection[j].max_extra);
                    let t = self.selection[j].extra;

                    if undone > 0 && t < self.selection[j].max_extra {
                        self.selection[j].max_extra += 1;
                        undone -= 1;
                        bag[j] -= self.coeff;
                        // Go to forwards section.
                        forwards = true;
                        break;
                    }

                    if t > 0 {
                        self.selection[j].extra = 0;
                        undone += t;
                        bag[j] += t * self.coeff;
                    }
                }
            }
            // If we got here via the innermost break, we continue to the forwards section.
            if !forwards {
                return false;
            }
        } else {
            undone = self.current_size;
        }

        // Forwards //
        let mut j: usize = 0;
        while undone > 0 {
            assert!(j < bag.len());

            let t: u32 = min(undone, self.selection[j].max_extra);
            if t > 0 {
                self.selection[j].extra = t;
                undone -= t;
                bag[j] -= t * self.coeff;
            }

            j += 1;
        }

        return true;
    }


    /*
    fn multiset_complex_nonlocal(
        &mut self,
        bag: &mut Vec<u32>,
        soluble: &mut Vec<Soluble>,
        find_first: bool,
    ) -> bool {
        let mut undone: u32;
        let bag_length = bag.len();

        if !find_first {
            if self.current_size > 0 {
                undone = 0;
                'backtrack: for j in 0..bag_length {
                    assert!(self.selection[j].extra <= self.selection[j].max_extra);
                    let t = self.selection[j].extra;

                    if undone > 0 && t < self.selection[j].max_extra {
                        let c = bag[j];

                        for e in 1..=undone {
                            assert!((t + e) <= (self.selection[j].max_extra));
                            let idx = c - self.coeff;
                            if soluble[idx as usize].min != Soluble::INSOLUBLE {
                                self.selection[j].extra = t + e;
                                bag[j] = idx;
                                undone -= e;
                                continue 'forwards;
                            }
                        }
                    }
                    if t > 0 {
                        self.selection[j].extra = 0;
                        undone += t;
                        bag[j] += t * self.coeff;
                    }
                }
            }
            return false;
        } else {
            undone = self.current_size;
        }

        'forwards: for j in 0..bag_length {
            assert!(j < bag_length);
            let t = self.selection[j].max_extra;
            if t <= undone {
                if t > 0 {
                    self.selection[j].extra = t;
                    undone -= t;
                    bag[j] -= t * self.coeff;
                }
            } else {
                self.selection[j].extra = undone;
                bag[j] -= undone * self.coeff;
                undone = 0;
                if soluble[bag[j] as usize].min == Soluble::INSOLUBLE {
                    continue 'backtrack;
                }
            }
        }
        return true;
    }
 */

    /// Find a selection from a multiset by undoing the previous selection until
    /// the selected amount of some element can be increased by one (without
    /// exceeding overall selection size or violating solubility constraints).
    /// Then make up the size of the selection by selecting the earliest elements
    /// available (backtracking if this violates solubility constraints).
    pub(crate) fn multiset_complex(
        &mut self,
        bag: &mut Vec<u32>,
        soluble: &mut Vec<Soluble>,
        mut find_first: bool,
    ) -> bool {
      let mut undone: u32;
      let bag_length = bag.len();

      // The control flow here is bananas, because Maude uses `GOTO`, which is considered bad.

      if find_first {
        undone = self.current_size;
      } else {
        if self.current_size > 0 {
          undone = 0;
          // How to skip the forward block in this case? `!find_first` is true, so we use `find_first` as a flag.
        } else {
          // The case `!find_first && self.current_size == 0`:
          return false;
        }
      }

      // The labels on the loops refer to which code block will run next, although the labeling is confusing in the
      // sense that the control flow FWD->BCK uses `break 'backtrack`, while the flow BCK->FWD uses `continue 'forward`.
      'forward: loop{
        'backtrack: loop {
          // if !find_first, we want to skip the forward block, but we only do this ONCE, the first time through.
          if !find_first{
            // This is not semantically true, but it serves as a flag so that we don't skip the forward block the next
            // time through the loop.
            find_first = true;
            break 'backtrack;
          }

          // The FORWARED block //
          for j in 0..bag_length {
            assert!(j < bag_length);
            let t = self.selection[j].max_extra;
            if t <= undone {
              if t > 0 {
                self.selection[j].extra = t;
                undone -= t;
                bag[j] -= t * self.coeff;
              }
            } else {
              self.selection[j].extra = undone;
              bag[j] -= undone * self.coeff;
              undone = 0;
              if soluble[bag[j] as usize].min == Soluble::INSOLUBLE {
                // Jump to the second half of the outer loop, which contains the backtrack block.
                // Same as `goto BACKTRACK block`
                break 'backtrack;
              }
            }
          }
          // If we fall all the way through the forward block, we don't loop but rather return true.
          return true;
        }

        // The BACKTRACK block //
        for j in 0..bag_length {
          assert!(self.selection[j].extra <= self.selection[j].max_extra);
          let t = self.selection[j].extra;

          if undone > 0 && t < self.selection[j].max_extra {
            let mut c = bag[j];

            for e in 1..=undone {
              assert!((t + e) <= (self.selection[j].max_extra));
              c -= self.coeff;
              if soluble[c as usize].min != Soluble::INSOLUBLE {
                self.selection[j].extra = t + e;
                bag[j] = c;
                undone -= e;
                // Same as `goto FORWARD block`
                continue 'forward;
              }
            }
          }
          if t > 0 {
            self.selection[j].extra = 0;
            undone += t;
            bag[j] += t * self.coeff;
          }
        }
        // If we fall through backtrack we return false
        return false;
      }

      // Why do I get an unreachable warning for an unreachable expression?!
      // unreachable!("Impossible code path in multiset_complex.")

    }





}

// impl Default for Row {
//     fn default() -> Self {
//         Self {
//           name: Default::default(),
//           coeff: Default::default(),
//           min_size: Default::default(),
//           min_product: Default::default(),
//           min_leave: Default::default(),
//           max_size: Default::default(),
//           max_product: Default::default(),
//           max_leave: Default::default(),
//           current_size: Default::default(),
//           current_max_size: Default::default(),
//           selection: Default::default(),
//           soluble: Default::default()
//         }
//     }
// }


impl Ord for Row{
    fn cmp(&self, other: &Self) -> Ordering {
      let t: Ordering = self.coeff.cmp(&other.coeff);

      if t == Ordering::Equal {
        // Break ties with `max_size`.
        self.max_size.cmp(&other.max_size)
      } else {
        t
      }
    }
}

impl PartialOrd for Row {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for Row{}

impl PartialEq for Row {
    fn eq(&self, other: &Self) -> bool {
        self.coeff == other.coeff && self.max_size == other.max_size
    }
}
