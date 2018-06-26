/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _aspect_initial_temperature_ascii_data_vs_h
#define _aspect_initial_temperature_ascii_data_vs_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/initial_temperature/S40RTS_perturbation.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialTemperature
  {
    using namespace dealii;

    /**
     * A class that implements a prescribed temperature field determined from
     * a Vs model given by an AsciiDataVs input file.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class AsciiDataVs : public Utilities::AsciiDataInitial<dim>, public Interface<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        AsciiDataVs ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize ();

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataInitial<dim>::initialize;

        /**
         * Return vs  as a function of position. For the
         * current class, this function returns value from the text files.
         */
        double
        ascii_grid_vs (const Point<dim> &position) const;

        /**
         * Return the boundary temperature as a function of position. For the
         * current class, this function returns value from the text files.
         */
        double
        initial_temperature (const Point<dim> &position) const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
        * Allows us to call functions from S40RTS_perturbation.cc
        * Declares class
        */
        S40RTSPerturbation<dim> s40rts;

//        No longer need these as everything is called from S40RTS_perturbations
//        /**
//        * An enum to describe which method should be chosen to scale vs to density.
//        */
//        enum VsToDensityMethod
//		{
//        	file,
//			constant
//		};

//        /**
//         * Currently chosen source for vs to density scaling.
//         */
//        VsToDensityMethod vs_to_density_method;

//        /**
//         * The parameters below describe the perturbation of shear wave
//         * velocity into a temperatures perturbation. The first parameter is
//         * constant so far but could be made depth dependent as constraint by
//         * e.g. Forte, A.M. & Woodward, R.L., 1997. Seismic-geodynamic
//         * constraints on three- dimensional structure, vertical flow, and
//         * heat transfer in the mantle, J. Geophys. Res. 102 (B8),
//         * 17,981-17,994.
//         * The last parameter is a depth down to which heterogeneities are
//         * zeroed out.
//         */
//        double vs_to_density_constant;
//        double thermal_alpha;
//        double no_perturbation_depth;

//        /**
//         * This parameter gives the reference temperature, which will be
//         * perturbed. In the compressional case the background temperature
//         * will be the adiabat.
//         */
//        double reference_temperature;

//        /**
//         * Object containing the data profile.
//         */
//        aspect::Utilities::AsciiDataProfile<dim> profile;

//        /**
//         * The column index of the vs to density scaling in the data file
//         */
//        unsigned int vs_to_density_index;
    };
  }
}


#endif
