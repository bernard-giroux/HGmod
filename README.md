# HGmod
## A hydrogeophysical synthetic model generator

HGmod is a computer program that builds on stochastic realizations of porosity fields to derive electrical conductivity, dielectric permittivity and hydraulic permeability models. The presence of clay, the influence of salinity as well as temperature of the fluid of imbibition are taken into account in the underlying formulations. The saturated and unsaturated zones are also considered through the application of a saturation profile on the porosity field. A micro-geometrical model is used to relate the porosity to the clay fraction. This model also used to derive an expression for the pore specific surface of the sand-clay mixture. The specific surface is subsequently used to compute the conduction at the pore surface when building electrical conductivity models. Dielectric permittivity fields are built by successive applications of either the Hanai-Bruggeman or Maxwell-Garnett mixing models, depending on the relative proportions of sand, clay, water or air. In addition, the dielectric permittivity of water and clay follow a Cole-Cole behavior. HGmod is therefore a versatile tool useful to generate synthetic datasets needed to anticipate the geophysical response under specific conditions and to study hydrogeophysical sensitivity or resolution analysis.

For more information, see the paper

- B. Giroux and M. Chouteau, 2008. A hydrogeophysical synthetic model generator, _Computers & Geosciences_, 34, 9, 1080-1092, 2008 [DOI: 10.1016/j.cageo.2007.11.006]
