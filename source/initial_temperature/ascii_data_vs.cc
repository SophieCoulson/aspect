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


#include <aspect/global.h>
#include <aspect/initial_temperature/ascii_data_vs.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace InitialTemperature
  {
    template <int dim>
    AsciiDataVs<dim>::AsciiDataVs ()
    {}


    template <int dim>
    void
    AsciiDataVs<dim>::initialize ()
    {
      Utilities::AsciiDataInitial<dim>::initialize(1);
      if (s40rts.vs_to_density_method == s40rts.file)
              {
                s40rts.profile.initialize(this->get_mpi_communicator());
                s40rts.vs_to_density_index = s40rts.profile.get_column_index_from_name("vs_to_density");
              }
    }

    //Read in and return vs perturbation from ascii grid
    template <int dim>
    double
    AsciiDataVs<dim>::
    ascii_grid_vs (const Point<dim> &position) const
    {
      const double vs_perturbation = Utilities::AsciiDataInitial<dim>::get_data_component(position,0);
      return vs_perturbation;
    }


    //Read in ascii grid vs data above 700 km and S40RTS vs data below 700 km
    template <int dim>
    double
    AsciiDataVs<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      double vs_perturbation;
      if (this->get_geometry_model().depth(position) < 700000)
        {
          vs_perturbation = ascii_grid_vs(position);
        }
      else
        {
          vs_perturbation = s40rts.get_Vs(position);
        }

      // use either the user-input reference temperature as background temperature
      // (incompressible model) or the adiabatic temperature profile (compressible model)
      const double background_temperature = this->get_material_model().is_compressible() ?
    		  this->get_adiabatic_conditions().temperature(position) :
              s40rts.reference_temperature;

      // Get the vs to density conversion
      const double depth = this->get_geometry_model().depth(position);

            double vs_to_density = 0.0;
            if (s40rts.vs_to_density_method == s40rts.file)
              vs_to_density = s40rts.profile.get_data_component(Point<1>(depth), s40rts.vs_to_density_index);
            else if (s40rts.vs_to_density_method == s40rts.constant)
              vs_to_density = s40rts.vs_to_density_constant;
            else
              // we shouldn't get here but instead should already have been
              // kicked out by the assertion in the parse_parameters()
              // function
              Assert (false, ExcNotImplemented());

            // scale the perturbation in seismic velocity into a density perturbation
            // vs_to_density is an input parameter
            const double density_perturbation = vs_to_density * vs_perturbation;

            double temperature_perturbation;
            if (depth > s40rts.no_perturbation_depth)
              // scale the density perturbation into a temperature perturbation
              temperature_perturbation =  -1./s40rts.thermal_alpha * density_perturbation;
            else
              // set heterogeneity to zero down to a specified depth
              temperature_perturbation = 0.0;

            return background_temperature + temperature_perturbation;
    }




    template <int dim>
    void
    AsciiDataVs<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/test/",
                                                          "box_2d.txt");

        prm.declare_entry ("Reference temperature", "1600.0",
                           Patterns::Double (0),
                           "The reference temperature that is perturbed by the spherical "
                           "harmonic functions. Only used in incompressible models.");

// No longer need these as everything is called from S40RTS_perturbation
//        aspect::Utilities::AsciiDataProfile<dim>::declare_parameters(prm,
//                                                                     "$ASPECT_SOURCE_DIR/data/initial-temperature/S40RTS/",
//                                                                     "vs_to_density_Steinberger.txt",
//                                                                     "Ascii data vs to density model");
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    AsciiDataVs<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial temperature model");
      {
        Utilities::AsciiDataBase<dim>::parse_parameters(prm);
        s40rts.reference_temperature   = prm.get_double ("Reference temperature");

        //Parameters from AsciiDataProfile
        //profile.parse_parameters(prm,"Ascii data vs to density model");
      }
      prm.leave_subsection ();
      s40rts.initialize_simulator (this->get_simulator());

      // note: parse_parameters will call initialize for us
      s40rts.parse_parameters(prm);

    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(AsciiDataVs,
                                              "ascii data vs",
                                              "Implementation of a model in which the initial "
                                              "temperature is derived from files containing data "
                                              "in ascii format. Note the required format of the "
                                              "input data: The first lines may contain any number of comments "
                                              "if they begin with '#', but one of these lines needs to "
                                              "contain the number of grid points in each dimension as "
                                              "for example '# POINTS: 3 3'. "
                                              "The order of the data columns "
                                              "has to be `x', `y', 'Temperature [K]' in a 2d model and "
                                              " `x', `y', `z', 'Temperature [K]' in a 3d model, which means that "
                                              "there has to be a single column "
                                              "containing the temperature. "
                                              "Note that the data in the input "
                                              "files need to be sorted in a specific order: "
                                              "the first coordinate needs to ascend first, "
                                              "followed by the second and the third at last in order to "
                                              "assign the correct data to the prescribed coordinates. "
                                              "If you use a spherical model, "
                                              "then the data will still be handled as Cartesian, "
                                              "however the assumed grid changes. `x' will be replaced by "
                                              "the radial distance of the point to the bottom of the model, "
                                              "`y' by the azimuth angle and `z' by the polar angle measured "
                                              "positive from the north pole. The grid will be assumed to be "
                                              "a latitude-longitude grid. Note that the order "
                                              "of spherical coordinates is `r', `phi', `theta' "
                                              "and not `r', `theta', `phi', since this allows "
                                              "for dimension independent expressions.")
  }
}
