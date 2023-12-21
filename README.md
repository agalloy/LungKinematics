# LungKinematics
A set of tools and analysis pipelines for lung and lobar kinematics.

CURRENT TASK:
1) Develop a "gradient enhanced" HHD where the component of the vector field aligned with the (co)gradient of 
   square magnitude is reintroduced to the (co)exact component.

SUB TASKS:
- "Gradient enhanced" HHD [ ]
  - Perform normal HHD [ ]
  - Exact component [ ]
    - Compute the magnitude of the 1-form at each node [ ]
    - Take the gradient of this scalar field to get a new 1-form [ ]
    - Add back the components of the input aligned with the new 1-form to the exact component [ ]
  - Coexact component [ ]
    - Compute the magnitude of the 1-form at each face [ ]
    - Take the gradient of this scalar field to get a new 1-form [ ]
    - Add back the components of the input aligned with the new 1-form to the coexact component [ ]
  - Harmonic component is the remainder [ ]
  - Perform usual verification checks to ensure orthogonality and whatnot [ ]
  - Analyze results [ ]

STRUCTURAL CHANGES:
- HHD now takes in an optional options structure

What's required to make this more "independent"?
- Address External Dependencies
	- TecPlot Tools folder
	- Loading data

- Add missing analyses
	- FE simulation deformation analyses
	- Sliding trajectory analysis