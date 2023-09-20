# LungKinematics
A set of tools and analysis pipelines for lung and lobar kinematics.

CURRENT TASK:
Test different vector field reconstruction methods.

SUB TASKS:
- Implemented function for comparing 1-form and gradient based reconstruction. [X]
- Implemented least squares reconstruction as an "optimal" method. [X]
- Compared all three reconstructions for an exact planar vector field. [X]
- Compared all three reconstructions for a non-exact planar vector field. [X]
- Compared all three reconstructions for a lung vector field. [X]
- Saved previous version of HHD algorithm as HHD_GradientRecon.m [X]
- Least-squares reconstruction implemented into the main HHD algorithm. [X]

What's required to make this more "independent"?
- Address External Dependencies
	- TecPlot Tools folder
	- Loading data

- Add missing analyses
	- FE simulation deformation analyses
	- Sliding trajectory analysis