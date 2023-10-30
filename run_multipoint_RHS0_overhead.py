# avoid multithreading of numpy
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import time as time_package

from mpi4py import MPI

from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint
from openaerostruct.utils.constants import grav_constant

import openmdao.api as om


class Scalars2Vector(om.ExplicitComponent):
    """
    Concatenates scalars into a vector
    Inputs: scalar0, scalar1, scalar2, ..., scalar<nn>
    Output: vector = [scalar0, scalar1, scalar2, ..., scalar<nn>]
    """

    def initialize(self):
        self.options.declare('num_nodes', types=int, desc='length of the output vector')
        self.options.declare('units', types=str, default=None)

    def setup(self):
        nn = self.options['num_nodes']
        units = self.options['units']

        for i in range(nn):
            self.add_input('scalar' + str(i), shape=(1,), units=units)
        # END FOR

        self.add_output('vector', shape=(nn,), units=units)

        # partials
        for i in range(nn):
            self.declare_partials('vector', 'scalar' + str(i), rows=[i], cols=[0], val=1.)

    def compute(self, inputs, outputs):
        nn = self.options['num_nodes']
        vector = np.zeros(nn)
        for i in range(nn):
            vector[i] = inputs['scalar' + str(i)]
        # END FOR

        outputs['vector'] = vector


def main(num_points, CL_constraints='vector'):
    """
    Parameters
    ----------
    num_points : int
        number of analysis points
    CL_constraints : str
        'vector' or 'scalers'
        'vector' to constraint the CL vector, CL = [CL0, CL1, CL2, ...]
        'scalers' to constraint each CL, CL0, CL1, CL2, ... separately.
    """

    # angle of attack for each conditions
    alpha_all = np.linspace(1., 6., num_points)

    # Create a dictionary to store options about the surface
    mesh_dict = {"num_y": 5,
                "num_x": 2,
                "wing_type": "CRM",
                "symmetry": True,
                "num_twist_cp": 3,
                }

    mesh, twist_cp = generate_mesh(mesh_dict)

    surface = {
        # Wing definition
        "name": "wing",  # name of the surface
        "symmetry": True,  # if true, model one half of wing
        # reflected across the plane y = 0
        "S_ref_type": "wetted",  # how we compute the wing area,
        # can be 'wetted' or 'projected'
        "fem_model_type": "tube",
        "thickness_cp": 0.2 * np.ones(10),
        "twist_cp": twist_cp,
        "mesh": mesh,
        # Aerodynamic performance of the lifting surface at
        # an angle of attack of 0 (alpha=0).
        # These CL0 and CD0 values are added to the CL and CD
        # obtained from aerodynamic analysis of the surface to get
        # the total CL and CD.
        # These CL0 and CD0 values do not vary wrt alpha.
        "CL0": 0.0,  # CL of the surface at alpha=0
        "CD0": 0.015,  # CD of the surface at alpha=0
        # Airfoil properties for viscous drag calculation
        "k_lam": 0.05,  # percentage of chord with laminar
        # flow, used for viscous drag
        "t_over_c_cp": np.array([0.15]),  # thickness over chord ratio (NACA0015)
        "c_max_t": 0.303,  # chordwise location of maximum (NACA0015)
        # thickness
        "with_viscous": True,
        "with_wave": False,  # if true, compute wave drag
        # Structural values are based on aluminum 7075
        "E": 70.0e9,  # [Pa] Young's modulus of the spar
        "G": 30.0e9,  # [Pa] shear modulus of the spar
        "yield": 500.0e6 / 2.5,  # [Pa] yield stress divided by 2.5 for limiting case
        "mrho": 3.0e3,  # [kg/m^3] material density
        "fem_origin": 0.35,  # normalized chordwise location of the spar
        "wing_weight_ratio": 2.0,
        "struct_weight_relief": True,  # True to add the weight of the structure to the loads on the structure
        "distributed_fuel_weight": False,
        # Constraints
        "exact_failure_constraint": False,  # if false, use KS function
    }

    # Create the problem and assign the model group
    prob = om.Problem(reports=False)

    # Add problem information as an independent variables component
    indep_var_comp = om.IndepVarComp()
    indep_var_comp.add_output("v", val=248.136 * np.ones(num_points), units="m/s")  # velocity vector for each analysis point
    indep_var_comp.add_output("alpha", val=alpha_all, units="deg")   # AoA vector for each analysis point
    indep_var_comp.add_output("Mach_number", val=0.84)
    indep_var_comp.add_output("re", val=1.0e6, units="1/m")
    indep_var_comp.add_output("rho", val=0.38, units="kg/m**3")
    indep_var_comp.add_output("CT", val=grav_constant * 17.0e-6, units="1/s")
    indep_var_comp.add_output("R", val=11.165e6, units="m")
    indep_var_comp.add_output("W0", val=0.4 * 3e5, units="kg")
    indep_var_comp.add_output("speed_of_sound", val=295.4, units="m/s")
    indep_var_comp.add_output("load_factor", val=1.0)
    indep_var_comp.add_output("empty_cg", val=np.zeros((3)), units="m")

    prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

    # add wing geometry
    name = "wing"
    aerostruct_group = AerostructGeometry(surface=surface)
    prob.model.add_subsystem(name, aerostruct_group)

    # --- add analysis points ---

    # parallelize OAS analysis at each node
    parallel_group = prob.model.add_subsystem('parallel', om.ParallelGroup(), promotes=['*'])

    # loop over analysis points
    for i in range(num_points):
        point_name = "AS_point_{}".format(i)

        # monolithic problem
        # Create the aero point group and add it to the model
        AS_point = AerostructPoint(surfaces=[surface])
        parallel_group.add_subsystem(point_name, AS_point, promotes_inputs=["Mach_number", "re", "rho", "CT", "R", "W0", "speed_of_sound", "empty_cg", "load_factor"])

        # connect alpha and v
        prob.model.connect("alpha", point_name + ".alpha", src_indices=[i])
        prob.model.connect("v", point_name + ".v", src_indices=[i])

        # some intra-node OAS connections
        prob.model.connect(name + ".local_stiff_transformed", point_name + ".coupled." + name + ".local_stiff_transformed")
        prob.model.connect(name + ".element_mass", point_name + "." + "coupled." + name + ".element_mass")
        prob.model.connect(name + ".nodes", point_name + ".coupled." + name + ".nodes")
        prob.model.connect(name + ".mesh", point_name + ".coupled." + name + ".mesh")
        com_name = point_name + "." + name + "_perf"
        prob.model.connect(name + ".radius", com_name + ".radius")
        prob.model.connect(name + ".thickness", com_name + ".thickness")
        prob.model.connect(name + ".nodes", com_name + ".nodes")
        prob.model.connect(name + ".cg_location", point_name + "." + "total_perf." + name + "_cg_location")
        prob.model.connect(name + ".structural_mass", point_name + "." + "total_perf." + name + "_structural_mass")
        prob.model.connect(name + ".t_over_c", com_name + ".t_over_c")
        # END IF
    # END IF

    prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options["tol"] = 1e-6

    # Setup problem and add design variables, constraint, and objective
    prob.model.add_design_var("wing.twist_cp", lower=-10.0, upper=15.0)
    prob.model.add_design_var("wing.thickness_cp", lower=0.01, upper=0.5, scaler=1e2)
    prob.model.add_design_var("v", lower=200, upper=250, units="m/s")
    prob.model.add_constraint('alpha', units='deg', lower=-10, upper=10)

    # put all CD and CLs into a vector
    prob.model.add_subsystem('CD_vector', Scalars2Vector(num_nodes=num_points, units=None), promotes_outputs=[('vector', 'CD_vec')])
    prob.model.add_subsystem('CL_vector', Scalars2Vector(num_nodes=num_points, units=None), promotes_outputs=[('vector', 'CL_vec')])
    for i in range(num_points):
        prob.model.connect('AS_point_' + str(i) + '.CD', 'CD_vector.scalar' + str(i))
        prob.model.connect('AS_point_' + str(i) + '.CL', 'CL_vector.scalar' + str(i))

    # add objective (just random dummy objective)
    prob.model.add_objective('Mach_number', ref=0.1)

    # CL constraints
    if CL_constraints == 'vector':
        # 1) impose to the vector
        prob.model.add_constraint('CL_vec', equals=0.5, ref=0.1)
    elif CL_constraints == 'scalers':
        # 2) impose to each points
        for i in range(num_points):
            prob.model.add_constraint('AS_point_' + str(i) + '.CL', equals=0.5, ref=0.1)

    # setup problem
    prob.setup(mode='rev', check=False)

    # --- update solver settings for each OAS system ---
    for subsystem in prob.model.system_iter():
        if subsystem.name == 'coupled':
            # turn off nonlinear solver print
            subsystem.nonlinear_solver.options['iprint'] = 0

            # use Krylov solver for each OAS
            subsystem.linear_solver = om.PETScKrylov(assemble_jac=True, iprint=2, err_on_non_converge=True)
            subsystem.linear_solver.precon = om.LinearRunOnce(iprint=-1)
    # END FOR

    prob.run_model()
    prob.compute_totals()


if __name__ == '__main__':
    # compute the derivatives of the CL vector. This includes the RHS=0 overheads
    print('--- constraining CL vector ---')
    main(num_points=2, CL_constraints='vector')

    # now, compute the derivatives of each CL separately. This does not include the RHS=0 overheads
    print('--- constraining CL of each point ---')
    main(num_points=2, CL_constraints='scalers')
