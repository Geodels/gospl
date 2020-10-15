---
title: 'gospl: Global Scalable Paleo Landscape Evolution'
tags:
  - Python
  - global landscape evolution model
  - basin evolution and stratigraphy
authors:
 - name: Tristan Salles
   orcid: 0000-0001-6095-7689
   affiliation: "1"
 - name: Claire Mallard
   orcid: 0000-0003-2595-2414
   affiliation: "1"
 - name: Sabin Zahirovic
   orcid: 0000-0002-6751-4976
   affiliation: "1"
affiliations:
 - name: School of Geosciences, The University of Sydney, Australia
   index: 1
date: 20 October 2020
bibliography: paper.bib
---

**gospl** (short for *Global Scalable Paleo Landscape Evolution*) is an open source, GPL-licensed library providing a scalable parallelised Python-based numerical model to simulate landscapes and basins reconstruction at global scale.

# Introduction

The source-to-sink (S2S) concept quantifies the different components of sedimentary systems: from source areas, through dispersal systems, to deposition in a number of sedimentary sinks. When applied to ancient sedimentary systems, it can be used to make inferences or predictions about upstream control, or downstream evolution of paleo-landscapes and stratigraphic record [@Helland:16]. Such concept is of keen interest to Earth system scientists studying the role of atmospheric circulation on physical denudation, the influence of mantle convection on erosion and deposition patterns, the location and abundance of natural resources or the implications of morphological changes and catchments dynamics on the evolution of life.

# Statement of Need

Traditional techniques in characterising S2S have consisted in the development of physical experiments and the use of both modern and outcrop analogues to constrain the size, shape, and complexity of sedimentary bodies. Over the years, growing detailed global datasets have improved our general understanding of sedimentary systems and scaling relationships have been proposed to assess the morphology and connectivities between modern S2S segments (*e.g.* the catchment, continental shelf, continental slope and submarine fan) and infer their temporal variability [@Nyberg:18].

It is recognised that such approaches provide a useful, first-order assessment of S2S in ancient sedimentary regions where the complete sediment routing system is not preserved [@Hobley:11; @Bhattacharya:16]. However, they often oversimplify the temporal variability and complexity of the different forcing conditions influencing S2S systems.

The emergence of forward modelling of landscapes and sedimentary systems has proven to be a valuable  avenue to integrate S2S concepts into a process-based  framework that considers sedimentation, eustatic changes  and tectonic influences [@Granjeon:99; @Tucker:10; @Salles:11].  Since the '90s, many software have been designed to estimate long-term landscape evolution as well as sedimentary basins formation in response to various mechanisms such as tectonic or climatic forcing [@Braun:97; @Tucker:01; @Bianchi:15]. These models rely on a set of mathematical and physical expressions that simulate sediment erosion, transport and deposition and can reproduce the first order complexity of Earth' surface geomorphological and sedimentary evolution [@Braun:13; @Hobley:17, @Salles:18].

Yet, we are still missing a tool to evaluate global scale patterns of paleo-Earth surface and the interactions with both the atmosphere, the hydrosphere, the tectonics and mantle dynamics.  **gospl** is designed to address parts of this gap and to make the link between deep Earth, climate and Earth' surface. It pioneers innovative and efficient techniques to estimate and predict ancient landscapes and sedimentary successions at global scale and over deep time [@Salles:18b]. It can be used to test different hypothesis related to Earth past evolution and to characterise the role of several drivers such as precipitation, dynamic topography, sea-level on Earth landscape evolution and sedimentary basins formation.


# Package summary


**gospl** is able to simulate global-scale forward model of landscape evolution, dual-lithology (coarse and fine) sediment routing and stratigraphic history forced with deforming plate tectonics, paleotopographies and paleoclimate reconstructions. It relates the complexity of the triggers and responses of sedimentary processes from the complete sediment routing perspective accounting for different scenarii of plate motion, tectonic uplift/subsidence, climate, geodynamic and sedimentary conditions.

The following physical processes are considered:
- river incision using the stream power law based on a multiple flow direction approach and an implicit parallel implicit drainage area method [@Richardson:14],
- inland river deposition in depressions computed using priority-flood techniques [@Barnes:14],
- dual-lithology marine deposition at river mouth based on a diffusion algorithm,
- hillslope processes in both marine and inland areas [@Salles:19], and
- sediment compaction as stratigraphic layers geometry and properties changes [@Yuan:19].

As mentioned previously, **gospl** can be forced with spatially and temporally varying tectonics (horizontal and vertical displacements) and climatic forces (temporal and spatial precipitation changes and sea-level fluctuations).

Internally, **gospl** is mostly written in Python and takes advantage of PETSc solvers [@Balay:12] over parallel computing architectures using MPI [@Barney:12]. For a detailed description of the physical processes and implicit numerical scheme implementation, the user is invited to read the technical guide available within the documentation.

Additionally, the documentation contains some workflows and functions that helps, among others, to pre-process netCDF elevation and climate dataset, plate velocity files and to create unstructured spherical mesh used as input for **gospl**.  Outputs are generated as compressed hdf5 files [@hdf5] that can be simply visualise as a temporal series in Paraview [@Ahrens:05]. In addition, post-processing scripts are also provided to build the output as VTK [@Schroeder:06] or netCDF files [@Brown:93] and to extract some specific information from the results such as erosion deposition thicknesses, stratigraphic records or river sediment load over time.

Finally, we provide different approaches to get you started with **gospl** from Docker containers, to pip or conda. **gospl** depends on several libraries sometimes hard to install locally, so the recommended starting approach is with Docker to ensure compatibility of all dependencies. Eventually, pip and conda are possible considering all dependency requirements are met.



# References
