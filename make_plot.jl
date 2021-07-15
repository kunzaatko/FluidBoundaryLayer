using CSV, DataFrames, Plots, LaTeXStrings

function make_plots()
    if "-h" in ARGS || "--help" in ARGS
        print("""
    usage:
        julia -i make_plot.jl {{csv_filename}} {{lambda_value}} {{file_for_saving (in svg format)}}

        The file for save is optional.
    """)
        return # this is for the possibility to quit out of the include before running other stuff when wanting only help
    end

    pyplot()
    c = DataFrame(CSV.File(ARGS[1]))
    lambda = ARGS[2]

    plot(
        Vector(c[:, 1]),
        Vector(c[:, 2]);
        label=L"y(x)",
        title=L"\lambda = %$lambda",
        leg=:topleft,
    )

    plot!(Vector(c[:, 1]), Vector(c[:, 3]); label=L"y'(x)")

    plot!(Vector(c[:, 1]), Vector(c[:, 4]); label=L"y''(x)", show=true)

    if length(ARGS) == 3
        tex_output_file = ARGS[3]
        pgfplotsx()

        plot(
            Vector(c[:, 1]),
            Vector(c[:, 2]);
            label=L"y(x)",
            title=L"\lambda = %$lambda",
            leg=:topleft,
        )

        plot!(Vector(c[:, 1]), Vector(c[:, 3]); label=L"y'(x)")

        plot!(Vector(c[:, 1]), Vector(c[:, 4]); label=L"y''(x)")
        savefig(tex_output_file * ".svg")
    end
    return
end
make_plots()
