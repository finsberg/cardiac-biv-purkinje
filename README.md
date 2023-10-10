# Sandbox for creating idealized BiV geometry with fibers and purkinje system.
This repository contains code for generating data of an idealized bi-ventricular geometry with microstructure (fiber, sheets, and sheet normal directions) as well as a purkinje system.

It uses [cardiac-geometries](https://github.com/ComputationalPhysiology/cardiac_geometries) to create the idealized BiV geometry, [ldrb](https://github.com/finsberg/ldrb) to generate the fiber orientations, and [fractal-tree](https://github.com/finsberg/fractal-tree) to generate the Purkinje system.

## Install
In order to generate the mesh you need to install [gmsh](https://gmsh.info) and to run the the ldrb algorithm you would need [FEniCS](https://fenicsproject.org/download/archive/). These are both packages that can be quite complicated to install.
The easiest way to install all the necessary dependencies is therefore by using [docker](https://www.docker.com/get-started/)
To start a docker container that comes with gmsh and fenics pre-installed, run the following command

```
docker run --name cardiac -w /home/shared -v $PWD:/home/shared -it ghcr.io/scientificcomputing/fenics-gmsh:2023-08-16
```
This will start a new docker container and also share your current working directory (`$PWD`) with the container (at `/home/shared`). This means that any files you generate inside the container can be shared with your local machine an vica-versa.

Next, you need to install the rest of the dependencies. This can be done with `pip`
```
python3 -m pip install -r requirements.txt
```

## Generate data
To generate the mesh with correctly marked boundaries, the microstructure and the Purkinje system you can run `generate_data.py`
```
python3 generate_data.py
```

Once you have run this script a folder called `data` will be created that will contain mesh, boundaries, fibers and purkinje system.

The folder will contain the following files
```
data
├── biv_ellipsoid.msh       - Mesh generated from gmsh
├── ffun.h5                 - Data for facet markers
├── ffun.xdmf               - File for visualizing facet markers in Paraview
├── fiber.h5                - Data for fiber fields
├── fiber.xdmf              - File for visualizing fibers in Paraview
├── info.json               - Parameter used for generating mesh
├── lv_tree.vtu             - File for visualizing purkinje system in the LV in Paraview
├── lv_tree_endnodes.txt    - List of leaf nodes in the purkinje system in the LV
├── lv_tree_lines.txt       - Connectivity to create lines for purkinje system in the LV
├── lv_tree_xyz.txt         - Coordinates of all points in the purkinje system in the LV
├── markers.json            - List of markers for the facets
├── mesh.h5                 - Data of the mesh
├── mesh.xdmf               - File for visualizing mesh in Paraview
├── rv_tree.vtu             - File for visualizing purkinje system in the RV in Paraview
├── rv_tree_endnodes.txt    - List of leaf nodes in the purkinje system in the RV
├── rv_tree_lines.txt       - Connectivity to create lines for purkinje system in the RV
├── rv_tree_xyz.txt         - Coordinates of all points in the purkinje system in the RV
├── sheet.h5                - Data for sheet fields
├── sheet.xdmf              - File for visualizing sheets in Paraview
├── sheet_normal.h5         - Data for sheet normal fields
├── sheet_normal.xdmf       - File for visualizing sheet normals in Paraview
├── triangle_mesh.h5        - Surface mesh (essentially same as ffun)
└── triangle_mesh.xdmf      - Surface mesh (essentially same as ffun)
```