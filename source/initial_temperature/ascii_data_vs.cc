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


//Read in background temperature and calculate temperatures from Vs
    template <int dim>
    double
    AsciiDataVs<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      double vs;
      if (this->get_geometry_model().depth(position) < 700000)
        {
          vs = ascii_grid_vs(position);
        }
      else
        {
          vs = s40rts.get_Vs(position);
        }
//...convert vs into temperature, return that...

// use either the user-input reference temperature as background temperature
      // (incompressible model) or the adiabatic temperature profile (compressible model)
      // const double background_temperature = this->get_material_model().is_compressible() ?                                                        this->get_adiabatic_conditions().temperature(position) :
      // reference_temperature;
      const double background_temperature=1000;
      const double density_scaling=1;

      //const double ascii_grid_vs_perturbation=ascii_grid_vs (position);
      //const double S40RTS_vs_perturbation = s40rts.get_Vs (position);
      const double vs_perturbation = vs+1;
      //std::cout << vs_perturbation << std::endl;
      const double density_perturbation = vs_perturbation*density_scaling;
      const double temperature_scaling = 1;
      const double temperature_perturbation = density_perturbation*temperature_scaling*10;
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
        reference_temperature   = prm.get_double ("Reference temperature");
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
