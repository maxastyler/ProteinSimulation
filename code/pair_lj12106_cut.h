/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   # potential given by
   pair_style lj12106/cut cutoff
   pair_coeff type1 type2 c6 c10 c12 cutoff
   V(r) = c12/r^12 + c10/r^10 + c6/r^6
   F(r) = 12 c12/r^13 + 10 c10/r^11 + 6 c6/r^7
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lj12106/cut,PairLJ12106Cut)

#else

#ifndef LMP_PAIR_LJ12106_CUT_H
#define LMP_PAIR_LJ12106_CUT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJ12106Cut : public Pair {
 public:
  PairLJ12106Cut(class LAMMPS *);
  virtual ~PairLJ12106Cut();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

  void compute_inner();
  void compute_middle();
  void compute_outer(int, int);

 protected:
  double cut_global;
  double **cut;
  double **c6,**c10,**c12;
  double **lj1,**lj2,**lj3,**lj4,**lj5,**lj6,**offset;
  double *cut_respa;

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

*/
