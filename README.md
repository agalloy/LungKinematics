# LungKinematics
A set of tools and analysis pipelines for lung and lobar kinematics.

CURRENT TASK:
Examine the results of HHD on simple lung lobe transformations.

SUB TASKS:
- Implement a new analysis code that does the following [X]
	- Loads a lung mesh and lobar fissure [X]
	- Creates vector fields for various simple transformations on one of the lobes [X]
	- Restricts the vector field to the lobar fissure and applies HHD [X]
- Implement a function that identifies a contact interface between two surfaces [X]

STRUCTURAL CHANGES:


What's required to make this more "independent"?
- Address External Dependencies
	- TecPlot Tools folder
	- Loading data

- Add missing analyses
	- FE simulation deformation analyses
	- Sliding trajectory analysis