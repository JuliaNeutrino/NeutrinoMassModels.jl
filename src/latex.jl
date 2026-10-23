# LaTeX output of particle models, spectra and the light neutrinos: methods
# of GLoBES.latex.

function latex(s::NeutrinoSpectrum)
    p = s.params
    rows = [
        ["\$m_1, m_2, m_3\$", "\$" * join(latex_number.(s.m), ", ") * "\$ eV"],
        ["ordering", string(s.ordering)],
        [
            "\$\\sin^2\\theta_{12}, \\sin^2\\theta_{13}, \\sin^2\\theta_{23}\$",
            "\$" *
            join(latex_number.((sin(p.θ12)^2, sin(p.θ13)^2, sin(p.θ23)^2)), ", ") *
            "\$",
        ],
        [
            "\$\\delta_{\\mathrm{CP}}, \\alpha_{21}, \\alpha_{31}\$",
            "\$" *
            join(latex_number.(rad2deg.((p.δ, s.majorana_phases...))), "^\\circ, ") *
            "^\\circ\$",
        ],
        [
            "\$\\Delta m^2_{21}, \\Delta m^2_{31}\$",
            "\$" * latex_number(p.Δm21) * ", " * latex_number(p.Δm31) * "\$ eV\$^2\$",
        ],
        [
            "\$\\sum m_i, m_\\beta, m_{\\beta\\beta}\$",
            "\$" * join(latex_number.((s.Σm, s.mβ, s.mββ)), ", ") * "\$ eV",
        ],
    ]
    return latex_table(["Light neutrinos", ""], rows)
end

function latex(c::ModelCheck)
    passed(c) &&
        isempty(c.warnings) &&
        return "All $(c.terms) terms of the Lagrangian are invariant under the symmetries of the model.\n"
    items = [
        ["\\item \\textbf{Problem:} " * escape_latex(p) for p in c.problems];
        ["\\item " * escape_latex(w) for w in c.warnings]
    ]
    return (passed(c) ? "All $(c.terms) terms of the Lagrangian are invariant.\n" : "") *
           "\\begin{itemize}\n" *
           join(items, "\n") *
           "\n\\end{itemize}\n"
end

function latex(s::Spectrum)
    rows = [
        ["\$$(particle_latex_guess(n))\$", latex_number(m) * " GeV"] for (n, m) in s.scalars
    ]
    for (n, ms) in s.fermions
        push!(
            rows,
            ["\$$(particle_latex_guess(n))\$", join(latex_number.(ms), ", ") * " GeV"],
        )
    end
    for (k, v) in s.tadpoles
        push!(rows, ["\$$(guess_symbol_latex(k))\$ (tadpole)", latex_number(v)])
    end
    return latex_table(
        ["State", "Tree-level mass"],
        rows;
        caption = "Spectrum of the $(escape_latex(s.model)) model.",
    )
end
particle_latex_guess(n::Symbol) = guess_symbol_latex(n)


# ---------------------------------------------------------------------------
#  Particle models
# ---------------------------------------------------------------------------

function bsm_latex_names(m::BSMModel)
    names = Dict{Symbol,String}()
    for (k, info) in m.parameter_info
        haskey(info, :LaTeX) && (names[k] = info[:LaTeX])
    end
    for (k, info) in m.particle_info
        l = get(info, :LaTeX, nothing)
        l isa String && (names[k] = l)
    end
    return names
end

"""
    latex(model::BSMModel; Fields=true, Lagrangian=true, Check=true, Parameters=true) -> String

The particle content (with charges), the Lagrangian with LaTeX names from
`particles.m` and `parameters.m`, and the result of [`check_model`](@ref).

```julia
latex(scoto; Check = false)   # model; parts switched off by name
```
"""
function latex(m::BSMModel; kwargs...)
    on(k) = get(kwargs, k, true)
    names = bsm_latex_names(m)
    io = IOBuffer()
    println(
        io,
        "\\section{",
        escape_latex(isempty(m.latex_name) ? m.name : m.latex_name),
        "}",
    )
    isempty(m.authors) || println(io, "Model by ", escape_latex(m.authors), ".\n")
    if on(:Fields)
        header = ["Field", "Generations", "Components"]
        for g in m.gauge
            push!(
                header,
                g.group === :U1 ? "\$U(1)_{\\mathrm{$(g.name)}}\$" : "\$SU($(g.N))\$",
            )
        end
        for s in m.global_symmetries
            push!(header, s.kind === :Z ? "\$\\mathbb{Z}_$(s.N)\$" : "\$U(1)\$")
        end
        rows = Vector{String}[]
        for f in m.fields
            comps = join((symbol_latex(c, names) for c in f.components), ", ")
            f.conjugated && (comps = "(" * comps * ")^\\dagger")
            push!(
                rows,
                [
                    "\$$(symbol_latex(f.name, names))\$ ($(f.kind))",
                    string(f.generations),
                    "\$$comps\$",
                    (
                        q isa Rational ? "\$" * tex(q, names, 0) * "\$" : "\$$(q)\$" for
                        q in f.charges
                    )...,
                ],
            )
        end
        print(
            io,
            latex_table(
                header,
                rows;
                caption = "Fields of the $(escape_latex(m.name)) model (left-handed Weyl fermions).",
            ),
        )
    end
    if on(:Lagrangian)
        println(io, "\\subsection*{Lagrangian}\n\\begin{align}")
        lines = String[]
        for t in m.lagrangian
            push!(
                lines,
                "  \\mathcal{L}_{\\mathrm{$(escape_latex(String(t.name)))}} &= " *
                latex_expression(t.expr; names) *
                (t.hc ? " + \\mathrm{h.c.}" : ""),
            )
        end
        println(io, join(lines, " \\\\\n"), "\n\\end{align}")
    end
    if on(:Check)
        println(io, "\\subsection*{Symmetries}")
        print(io, latex(check_model(m)))
    end
    return String(take!(io))
end
