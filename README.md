# LungKinematics
A set of tools and analysis pipelines for lung and lobar kinematics.

CURRENT TASK:
Find the "optimal" 3D velocity fields for reducing the exact component along the fissure surface
that is a linear combination of a set of given input vector fields.

SUB TASKS:
- Optimization task [ ]
  - Generate a set of input vector fields for analysis [ ]
  - Assemble the covariance matrix for the inner product between the exact components of each vector field [ ]
  - Find the eigenvectors and values of this covariance matrix [ ]
  - Visualize vector fields associated with those eigenvectors [ ]
- Have the HHD algorithm output 

STRUCTURAL CHANGES:


What's required to make this more "independent"?
- Address External Dependencies
	- TecPlot Tools folder
	- Loading data

- Add missing analyses
	- FE simulation deformation analyses
	- Sliding trajectory analysis