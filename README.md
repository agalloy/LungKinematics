# LungKinematics
A set of tools and analysis pipelines for lung and lobar kinematics.

CURRENT TASK:
1) Implement the Galerkin Hodge Star (which gives the inner product of Whitney basis 1-forms) from
   de Goes et al. "Discrete 2-Tensor Fields..."

SUB TASKS:
- Implement the Galerkin Hodge Star [X]
  - Implement an option to toggle between it and the diagonal Hodge star [ ] 
- Verify that it has the same Laplace operator as the diagonal Hodge star [X]
- Verify orthogonality between HHD components [ ]
  - In the non-enhanced HHD [X]
  - In the enhanced HHD [ ]

STRUCTURAL CHANGES:
- 

What's required to make this more "independent"?
- Address External Dependencies
	- TecPlot Tools folder
	- Loading data

- Add missing analyses
	- FE simulation deformation analyses
	- Sliding trajectory analysis