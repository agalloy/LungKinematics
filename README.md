# LungKinematics
A set of tools and analysis pipelines for lung and lobar kinematics.

CURRENT TASK:
1) Create a flat but reasonable lobar fissure for analyzing boundary effects on HHD results.

SUB TASKS:
- Flat fissure HHD [X]
  - Load a subjects mesh and identify lobar fissure [X]
  - "Flatten" the lobar fissure [X]
  - Generate an input vector field on the flattened lobar fissure [X]
  - Perform HHD [X]
  - Plot HHD Components[X]
  - Analyze Results [X]

STRUCTURAL CHANGES:
- Grid Sample now works in 2D
- Flat fissure was incorporated to Test_HHDSimpleLobeTransformations.m



What's required to make this more "independent"?
- Address External Dependencies
	- TecPlot Tools folder
	- Loading data

- Add missing analyses
	- FE simulation deformation analyses
	- Sliding trajectory analysis