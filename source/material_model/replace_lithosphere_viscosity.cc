/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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

#include <aspect/material_model/replace_lithosphere_viscosity.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <limits>

#include <array>
#include <utility>



namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    ReplaceLithosphereViscosity<dim>::ReplaceLithosphereViscosity ()
      :
      lab_depths(1, 1.0)
    {}

    template <int dim>
    void
    ReplaceLithosphereViscosity<dim>::initialize()
    {
      base_model->initialize();

      if (LAB_depth_source == File)
        {
          const std::string filename = data_directory+LAB_file_name;

          if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
          std::cout << "   Loading Ascii data lookup file " << filename << "." << std::endl;

          lab_depths.load_file(filename,this->get_mpi_communicator());
        }
    }

    template <int dim>
    double
    ReplaceLithosphereViscosity<dim>::ascii_lab (const Point<2> &position) const
    {
      const double lab = lab_depths.get_data(position,0);
      return lab;
    }

    template <int dim>
    void
    ReplaceLithosphereViscosity<dim>::evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                                               typename Interface<dim>::MaterialModelOutputs &out) const
    {
      base_model->evaluate(in,out);

      for (unsigned int i=0; i < in.position.size(); ++i)
        {
          const double depth = this->SimulatorAccess<dim>::get_geometry_model().depth(in.position[i]);

          if (LAB_depth_source == File)
            {
              //Get spherical coordinates for model
              std::array<double,dim> scoord      = Utilities::Coordinates::cartesian_to_spherical_coordinates(in.position[i]);
              const double phi = scoord[1];
              const double theta = scoord[2];
              const Point<2> phi_theta (phi, theta);

              //Get lab depth for specific phi and theta
              const double lab_depth = ascii_lab(phi_theta);

              if (depth <= lab_depth)
                out.viscosities[i] = lithosphere_viscosity;
              else
                out.viscosities[i] *= 1;
            }

          else if (LAB_depth_source == Value)
            {
              if (depth <= max_depth)
                out.viscosities[i] = lithosphere_viscosity;
              else
                out.viscosities[i] *= 1;
            }

          else
            {
              Assert( false, ExcMessage("Invalid method for depth specification method") );
            }
        }
    }

    template <int dim>
    void
    ReplaceLithosphereViscosity<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Replace lithosphere viscosity");
        {
          prm.declare_entry("Base model","simple",
                            Patterns::Selection(MaterialModel::get_valid_model_names_pattern<dim>()),
                            "The name of a material model that will be modified by a replacing"
                            "the viscosity in the lithosphere by a constant value. Valid values for this parameter "
                            "are the names of models that are also valid for the "
                            "``Material models/Model name'' parameter. See the documentation for "
                            "more information.");
          prm.declare_entry ("Lithosphere viscosity", "1600",
                             Patterns::Double (0),
                             "The viscosity within lithosphere, applied above"
                             "the maximum lithosphere depth.");
          prm.declare_entry ("Depth specification method", "Value",
                             Patterns::Selection("File|Value"),
                             "Method that is used to specify the depth of the lithosphere-asthenosphere boundary.");
          prm.declare_entry ("Maximum lithosphere depth", "200000.0",
                             Patterns::Double (0),"Units: m."
                             "The maximum depth of the lithosphere. The viscosity will be "
                             "taken from the base model below this depth.");
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/initial-temperature/lithosphere-mask/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text `$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the `data/' subdirectory of ASPECT. ");
          prm.declare_entry ("LAB depth filename",
                             "LAB_CAM2016.txt",
                             Patterns::FileName (),
                             "File from which the lithosphere-asthenosphere boundary depth data is read.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    ReplaceLithosphereViscosity<dim>::parse_parameters (ParameterHandler &prm)
    {
      AssertThrow (dim == 3,
                   ExcMessage ("The 'Replace lithosphere viscosity' material model "
                               "is only available for 3d computations."));

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Replace lithosphere viscosity");
        {
          AssertThrow( prm.get("Base model") != "replace lithosphere viscosity",
                       ExcMessage("You may not use ``replace lithosphere viscosity'' as the base model for "
                                  "a replace lithosphere viscosity model.") );

          // create the base model and initialize its SimulatorAccess base
          // class; it will get a chance to read its parameters below after we
          // leave the current section
          base_model.reset(create_material_model<dim>(prm.get("Base model")));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
            sim->initialize_simulator (this->get_simulator());

          lithosphere_viscosity   = prm.get_double ("Lithosphere viscosity");

          if ( prm.get("Depth specification method") == "File" )
            {
              LAB_depth_source = File;
              data_directory = Utilities::expand_ASPECT_SOURCE_DIR (prm.get("Data directory"));
              LAB_file_name = prm.get("LAB depth filename");
            }
          else if ( prm.get("Depth specification method") == "Value" )
            {
              LAB_depth_source = Value;
              max_depth = prm.get_double ("Maximum lithosphere depth");
            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      /* After parsing the parameters for replace lithosphere viscosity, it is essential to parse
      parameters related to the base model. */
      base_model->parse_parameters(prm);
      this->model_dependence = base_model->get_model_dependence();
    }

    template <int dim>
    bool
    ReplaceLithosphereViscosity<dim>::
    is_compressible () const
    {
      return base_model->is_compressible();
    }

    template <int dim>
    double
    ReplaceLithosphereViscosity<dim>::
    reference_viscosity() const
    {
      /* Return reference viscosity from base model*/
      return base_model->reference_viscosity();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ReplaceLithosphereViscosity,
                                   "replace lithosphere viscosity",
                                   "The ``replace lithosphere viscosity'' Material model sets viscosity to a "
                                   "prescribed constant above the lithosphere-asthenosphere boundary (specified by "
                                   "an ascii file or maximum lithosphere depth). Below the lithosphere-asthenosphere"
                                   "boundary the viscosity is taken from any of the other available material model. "
                                   "In other words, it is a ``compositing material model''."
                                   "\n"
                                   "Parameters related to the depth dependent model are read from a subsection "
                                   "``Material model/replace lithosphere vicsosity''. "
                                   "The user must specify a ``Base model'' from which other material properties are "
                                   "derived.  "
                                   "\n"
                                   "Note the required format of the input data file: The first lines may "
                                   "contain any number of comments if they begin with ‘#’, but one of these lines "
                                   "needs to contain the number of grid points in each dimension as for example"
                                   "‘# POINTS: 3 3’. For a spherical model, the order of the data columns has to be"
                                   "'phi', 'theta','depth (m)', where phi is the  azimuth angle and theta is the "
                                   "polar angle measured positive from the north pole.")
  }
}
