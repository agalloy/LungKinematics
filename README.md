# LungKinematics
A set of tools and analysis pipelines for lung and lobar kinematics.

Task In Progress:
Test different vector field reconstruction methods.
	- Completed:
		- Implemented function for comparing 1-form and gradient based reconstruction.
		- Implemented least squares reconstruction as an "optimal" method.
		- Compared all three reconstructions for an exact planar vector field.
		- Compared all three reconstructions for a non-exact planar vector field.
		- Compared all three reconstructions for a lung vector field.
		- Saved previous version of HHD algorithm as HHD_GradientRecon.m
		- Least-squares reconstruction implemented into the main HHD algorithm.

What's required to make this more "independent"?
- Address External Dependencies
	- TecPlot Tools folder
	- Loading data

- Add missing analyses
	- FE simulation deformation analyses
	- Sliding trajectory analysis