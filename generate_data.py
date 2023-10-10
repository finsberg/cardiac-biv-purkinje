import numpy as np
import logging
import cardiac_geometries
from pathlib import Path
import dolfin
import ldrb
import meshio
from fractal_tree import generate_fractal_tree, FractalTreeParameters, Mesh

# Set the logging level to INFO

logging.basicConfig(level=logging.INFO)

# Create mesh
# -----------
# We will now create the mesh. If the mesh already exist we skip this step.
# Here we need to import `cardiac-geometries` and call the `biv_ellipsoid`
# function from the `gmsh` sub-package. We alse set the characteristic
# length to 0.1 which will generate a finer mesh that what is the default.

datadir = Path("data")
msh_file = datadir / "biv_ellipsoid.msh"
if not msh_file.is_file():
    # Generate mesh with cardiac geometries
    # Use `cardiac-geometries create-biv-ellipsoid --help to see all options`

    cardiac_geometries.create_biv_ellipsoid(
        datadir,
        char_length=0.2,  # Reduce this value to get a finer mesh
        center_lv_y=0.2,
        center_lv_z=0.0,
        a_endo_lv=5.0,
        b_endo_lv=2.2,
        c_endo_lv=2.2,
        a_epi_lv=6.0,
        b_epi_lv=3.0,
        c_epi_lv=3.0,
        center_rv_y=1.0,
        center_rv_z=0.0,
        a_endo_rv=6.0,
        b_endo_rv=2.5,
        c_endo_rv=2.7,
        a_epi_rv=8.0,
        b_epi_rv=5.5,
        c_epi_rv=4.0,
    )

# Create fibers
# -------------

# In this step we will create the microstructure, i,e the fiber
# sheet and sheet normal directions. To do this we use the
# Laplace-Dirichlet Rule Based (ldrb) algorithm from the
# ldrb python library (pip install ldrb)

fiber_file = datadir / "fiber.xdmf"
if not fiber_file.is_file():
    geo = cardiac_geometries.geometry.Geometry.from_folder(datadir)

    ldrb_markers = {
        "base": geo.markers["BASE"][0],
        "lv": geo.markers["ENDO_LV"][0],
        "rv": geo.markers["ENDO_RV"][0],
        "epi": geo.markers["EPI"][0],
    }

    system = ldrb.dolfin_ldrb(
        mesh=geo.mesh,
        fiber_space="P_1",
        ffun=geo.ffun,
        markers=ldrb_markers,
        alpha_endo_lv=60,
        alpha_epi_lv=-60,
        beta_endo_lv=0,
        beta_epi_lv=0,
    )

    for attr in ["fiber", "sheet", "sheet_normal"]:
        f = getattr(system, attr)
        with dolfin.XDMFFile(f"{datadir}/{attr}.xdmf") as xdmf:
            xdmf.write_checkpoint(f, attr, 0.0, dolfin.XDMFFile.Encoding.HDF5, False)


# Create purkinje system
# ----------------------

# Create purkinje system for LV using fractal tree algorithm

lv_file = datadir / "lv_tree.vtu"
if not lv_file.is_file():
    # Read gmsh file
    msh = meshio.read(msh_file)

    # Extract LV endocardium
    tag = msh.field_data["ENDO_LV"][0]
    inds = [i for i, x in enumerate(msh.cell_data["gmsh:physical"]) if x[0] == tag]
    connectivity = np.vstack([msh.cells[i].data for i in inds])

    # Define initial node for the purkinje system in the LV (use paraview to find this)
    init_node = [0, 2.34484, 0.19]

    # Now that we have the approximate coordinate, we need to find the closest node
    # in the mesh. We can do this by selecting the point with the shortest
    # distance to this node.
    verts = msh.points
    index = np.linalg.norm(np.subtract(verts, init_node), axis=1).argmin()

    # We create the mesh, and pass in the initial node
    mesh = Mesh(verts=verts, connectivity=connectivity, init_node=verts[index, :])

    # We also set the number of generations to 15 and we set the initial direction
    # to point in the positive $x$-direction which is the direction from the base
    # towards the apex (this can also be found by visualizing the mesh in Paraview).
    # We also specify the initial length to be 3. This is about the length of the
    # septal wall, so that the first branch will represent the His bundle.
    # We set the length to 0.25 which will be the length each successive branch.
    #
    # First we set a seed so that we know that the results are reproducible
    np.random.seed(1234)

    param = FractalTreeParameters(
        filename=datadir / lv_file.stem,
        init_length=7.0,
        N_it=12,
        length=0.5,
        initial_direction=np.array([1, 0, 0]),
    )

    # Next we create the Purkinje networks for the LV
    lv_tree = generate_fractal_tree(mesh, param)


# Create purkinje system for RV using fractal tree algorithm

rv_file = datadir / "rv_tree.vtu"
if not rv_file.is_file():
    # Read gmsh file
    msh = meshio.read(msh_file)

    tag = msh.field_data["ENDO_RV"][0]
    inds = [i for i, x in enumerate(msh.cell_data["gmsh:physical"]) if x[0] == tag]
    connectivity = np.vstack([msh.cells[i].data for i in inds])

    # and defining the initial node (again, we use Paraview to find this location).

    init_node = [0, 3.19, 0.19]
    verts = msh.points
    index = np.linalg.norm(np.subtract(verts, init_node), axis=1).argmin()
    mesh = Mesh(verts=verts, connectivity=connectivity, init_node=verts[index, :])

    # Set a new seed (note that if you don't like the generated system
    # you can try a different seed)
    np.random.seed(1234)

    param = FractalTreeParameters(
        filename=datadir / rv_file.stem,
        init_length=9.7,
        N_it=15,
        length=0.5,
        initial_direction=np.array([1, 0, 0]),
    )

    rv_tree = generate_fractal_tree(mesh, param)
