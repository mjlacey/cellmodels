# Basic modelling of battery cells and cell properties

This repository contains a framework for constructing basic models of battery cells, and calculation/analysis of cell properties. That is, the models are constructed for particular materials, dimensions, component thicknesses, to estimate properties such as cell capacity, energy density, mass breakdown, and so on.

There are several types of cell models currently defined in the `cellmodel.jl` file:

* `CylindricalCell` - a model for a generalised cylindrical cell, with specified diameter and height

* `PrismaticCell_JellyRoll` - a model for a generalised prismatic cell, with one or more wound jelly roll with the axis of winding parallel to the 'width' axis, and considering the tabs at the side of the cell

* `PrismaticCell_CinnamonRoll` - a model for a generalised prismatic cell, with one or more wound jelly roll with the axis of winding perpendicular to the 'width' axis, and considering the tabs at the top of the cell

* `PrismaticCell_Stacked` - a model for a stacked prismatic cell, with one or more stacked electrode assemblies, considering the tabs at the side of the cell

* `PouchCell` - a model for a stacked pouch cell of a defined number of layers.

# Intentions and limitations

The code is intended for battery scientists, engineers and other interested persons as a useful guide for estimating hypothetical battery cell properties based on known materials. This can be useful for predicting the basic characteristics of as-yet-undeveloped battery chemistries, cell formats, or other design considerations.

This is not designed to be a perfectly accurate representation of a real battery cell. There are currently a number of approximations and other fudge factors (it does not consider small details), and it does not provide any estimations related to, e.g., power, lifetime.

It also has almost nothing in the way of error handling at present, so it will not provide any warnings against unrealistic estimations. Users of this code for 'real' work are urged to use careful judgement...

# Get started

You must have Julia. The current code is tested to work with v1.5.3. If you are a first time user of Julia, I would recommend starting with [JuliaPro](https://juliacomputing.com/products/juliapro/), which bundles Julia with the Juno IDE.

To run this code you will also need to install the following add-on packages: `Measurements`, `Parameters`, `QuadGK`, `DataFrames`.

Examples of the use of this code are given in the Jupyter notebooks included in this repository.

Contributions and suggestions for further development are welcome.
