"""
Neutrino mass models for GLoBES.jl: particle models read from
Mathematica-style model files (fields, symmetries, Lagrangian), their
tree-level spectra and LesHouches files, and neutrino mass mechanisms (the
scotogenic loop and the type-I seesaw) with the Casas-Ibarra
parametrisation and lepton flavour violation.
"""
module NeutrinoMassModels

using LinearAlgebra, Printf
using GLoBES
using GLoBES:
    MExpr, MNumber, parse_mathematica, is_call, complex_if_negative, G_FERMI, flavour
using GLoBES: escape_latex, guess_symbol_latex, latex_number, symbol_latex, tex
import GLoBES: latex

include("particlemodel.jl")
include("spectrum.jl")
include("mechanisms.jl")
include("latex.jl")

export BSMModel,
    load_particle_model,
    check_model,
    add_gauge!,
    add_global!,
    add_fermion!,
    add_scalar!,
    add_lagrangian!,
    add_vevs!,
    add_parameter!,
    spectrum,
    read_leshouches,
    leshouches_values,
    write_leshouches,
    ScotogenicLoop,
    TypeISeesaw,
    neutrino_mass_matrix,
    neutrino_spectrum,
    casas_ibarra,
    orthogonal_matrix,
    br_lfv,
    scotogenic_masses

end
