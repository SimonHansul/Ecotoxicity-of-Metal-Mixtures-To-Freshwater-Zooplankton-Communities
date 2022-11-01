using Plots, StatsPlots, Plots.Measures, LaTeXStrings
using DocStringExtensions
import Plots:RecipesPipeline
import StatsPlots:_cycle, grouped_xy

@userplot RugPlot
@recipe function f(h::RugPlot)
    if length(h.args)>1
        error("Rugplot does not support multiple arguments.")
    end
    @series begin
        seriestype := :scatter
        markershape := :vline
        x, y = h.args[1], repeat([0], length(h.args[1]))
    end
end

@userplot LinePlot
@recipe function f(h::LinePlot; estimator=mean)
    let x, y, q, mask
        # if no quantile argument given at all => assume q=(0.05, 0.95)
        if length(h.args)>=3
            x, y,  q = h.args
        else
            x, y, = h.args
            q = (.05, .95)
        end

        mask = isfinite.(x) .& isfinite.(y) 
        x = x[mask]
        y = y[mask]
        
        x_ag = []
        y_mn = []
        y_lo = []
        y_hi = []

        for (j,xi) in enumerate(unique(x))
            mask = x.==xi
            push!(x_ag, xi)
            push!(y_mn, estimator(y[mask]))
            push!(y_lo, quantile(y[mask], q[1]))
            push!(y_hi, quantile(y[mask], q[2]))
        end

        @series begin
            seriestype --> :path
            ribbon := @. (y_mn - y_lo, y_hi - y_mn)
            x_ag, y_mn
        end
    end
end

@userplot GroupedLinePlot
@recipe function f(h::GroupedLinePlot; cmap=palette(:tab10).colors, estimator=mean)
    let x, y, group, q, mask
        # if no quantile argument given at all => assume q=(0.05, 0.95)
        if length(h.args)>=4
            x, y, group, q = h.args
        else
            x, y, group = h.args
            q = (.05, .95)
        end
        if typeof(group[1])<:Number
            mask = isfinite.(x) .& isfinite.(y) .& isfinite.(group) .& (ismissing.(x).==false) .& (ismissing.(y).==false) .& (ismissing.(group).==false)
        else
            mask = isfinite.(x) .& isfinite.(x) .& (ismissing.(x).==false) .& (ismissing.(y).==false)
        end
        x = x[mask]
        y = y[mask]
        group = group[mask]
        # aggregate values
        x_ag = []
        y_mn = []
        y_lo = []
        y_hi = []
        g_ag = []
        for (i,g) in enumerate(unique(group))
            for (j,xi) in enumerate(unique(x[group.==g]))
                mask = (group.==g).&(x.==xi)
                push!(x_ag, xi)
                push!(y_mn, estimator(y[mask]))
                push!(y_lo, quantile(y[mask], q[1]))
                push!(y_hi, quantile(y[mask], q[2]))
                push!(g_ag, g)
            end
        end
        for (i,g) in enumerate(unique(g_ag))
            @series begin
                mask = g_ag.==g
                seriestype --> :path
                ribbon := (y_mn[mask] .- y_lo[mask], y_hi[mask] .- y_mn[mask])
                x_ag[mask], y_mn[mask]
            end
        end
    end
end

notch_width(q2, q4, N) = 1.58 * (q4-q2)/sqrt(N)

"""
This series recipe is adopted from StatsPlots. 
It does the same as StatsPlots.boxplot, except that the extent of box and whiskers are passed on as quantiles. <br>
`Estimator` is a function to calculate the location of the middle line (mean by default).
"""
@recipe function f(
    ::Type{Val{:boksplot}},
    x,
    y,
    z;
    estimator=mean,
    q_box=(0.25, 0.75),
    q_whiskers=(0.05, 0.95),
    notch=false,
    outliers=true,
    whisker_width=:half
)
    # if only y is provided, then x will be UnitRange 1:size(y,2)
    if typeof(x) <: AbstractRange
        if step(x) == first(x) == 1
            x = plotattributes[:series_plotindex]
        else
            x = [getindex(x, plotattributes[:series_plotindex])]
        end
    end
    xsegs, ysegs = Segments(), Segments()
    texts = String[]
    glabels = sort(collect(unique(x)))
    warning = false
    outliers_x, outliers_y = zeros(0), zeros(0)
    bw = plotattributes[:bar_width]
    isnothing(bw) && (bw = 0.8)
    @assert whisker_width == :match || whisker_width == :half || whisker_width >= 0 "whisker_width must be :match, :half, or a positive number"
    ww = whisker_width == :match ? bw :
         whisker_width == :half ? bw / 2 :
         whisker_width
    for (i, glabel) in enumerate(glabels)
        # filter y
        values = y[filter(i -> StatsPlots._cycle(x, i) == glabel, 1:length(y))]

        # compute quantiles
        estimate = estimator(values) # formerly q3
        box_lo = quantile(values, q_box[1]) # formerly q2
        box_hi = quantile(values, q_box[2]) # formerly q4
        whisker_lo = quantile(values, q_whiskers[1]) # q1
        whisker_hi = quantile(values, q_whiskers[2]) # q5
        #q1, q2, q3, q4, q5 = quantile(values, range(0, stop = 1, length = 5))
        q1, q2, q3, q4, q5 = whisker_lo, box_lo, estimate, box_hi, whisker_hi
        
        # notch
        n = notch_width(q2, q4, length(values))

        # warn on inverted notches?
        if notch && !warning && ((q2 > (q3 - n)) || (q4 < (q3 + n)))
            @warn("Boxplot's notch went outside hinges. Set notch to false.")
            warning = true # Show the warning only one time
        end

        # make the shape
        center = Plots.discrete_value!(plotattributes[:subplot][:xaxis], glabel)[1]
        hw = 0.5_cycle(bw, i) # Box width
        HW = 0.5_cycle(ww, i) # Whisker width
        l, m, r = center - hw, center, center + hw
        lw, rw = center - HW, center + HW

        # internal nodes for notches
        L, R = center - 0.5 * hw, center + 0.5 * hw

        # outliers
        if q_whiskers != (0., 0.)  # if the range is 0.0, the whiskers will extend to the data
            #limit = whisker_range * (q4 - q2)
            limit = q5 - q1
            inside = Float64[]
            for value in values
                if (value < (q2 - limit)) || (value > (q4 + limit))
                    if outliers
                        push!(outliers_y, value)
                        push!(outliers_x, center)
                    end
                else
                    push!(inside, value)
                end
            end
            # change q1 and q5 to show outliers
            # using maximum and minimum values inside the limits
            q1, q5 = Plots.ignorenan_extrema(inside)
            q1, q5 = (min(q1, q2), max(q4, q5)) # whiskers cannot be inside the box
        end
        # Box
        push!(xsegs, m, lw, rw, m, m)       # lower T
        push!(ysegs, q1, q1, q1, q1, q2)    # lower T
        push!(
            texts,
            "Lower fence: $q1",
            "Lower fence: $q1",
            "Lower fence: $q1",
            "Lower fence: $q1",
            "Q1: $q2",
            "",
        )

        if notch
            push!(xsegs, r, r, R, L, l, l, r, r) # lower box
            push!(xsegs, r, r, l, l, L, R, r, r) # upper box

            push!(ysegs, q2, q3 - n, q3, q3, q3 - n, q2, q2, q3 - n) # lower box
            push!(
                texts,
                "Q1: $q2",
                "Median: $q3 ± $n",
                "Median: $q3 ± $n",
                "Median: $q3 ± $n",
                "Median: $q3 ± $n",
                "Q1: $q2",
                "Q1: $q2",
                "Median: $q3 ± $n",
                "",
            )

            push!(ysegs, q3 + n, q4, q4, q3 + n, q3, q3, q3 + n, q4) # upper box
            push!(
                texts,
                "Median: $q3 ± $n",
                "Q3: $q4",
                "Q3: $q4",
                "Median: $q3 ± $n",
                "Median: $q3 ± $n",
                "Median: $q3 ± $n",
                "Median: $q3 ± $n",
                "Q3: $q4",
                "",
            )
        else
            push!(xsegs, r, r, l, l, r, r)         # lower box
            push!(xsegs, r, l, l, r, r, m)         # upper box
            push!(ysegs, q2, q3, q3, q2, q2, q3)   # lower box
            push!(
                texts,
                "Q1: $q2",
                "Median: $q3",
                "Median: $q3",
                "Q1: $q2",
                "Q1: $q2",
                "Median: $q3",
                "",
            )
            push!(ysegs, q4, q4, q3, q3, q4, q4)   # upper box
            push!(texts, "Q3: $q4", "Q3: $q4", "Median: $q3", "Median: $q3", "Q3: $q4", "Q3: $q4", "")
        end

        push!(xsegs, m, lw, rw, m, m)             # upper T
        push!(ysegs, q5, q5, q5, q5, q4)          # upper T
        push!(
            texts,
            "Upper fence: $q5",
            "Upper fence: $q5",
            "Upper fence: $q5",
            "Upper fence: $q5",
            "Q3: $q4",
            "",
        )

    end

    if !Plots.isvertical(plotattributes)
        # We should draw the plot horizontally!
        xsegs, ysegs = ysegs, xsegs
        outliers_x, outliers_y = outliers_y, outliers_x

        # Now reset the orientation, so that the axes limits are set correctly.
        orientation := default(:orientation)
    end

    @series begin
        # To prevent linecolor equal to fillcolor (It makes the median visible)
        if plotattributes[:linecolor] == plotattributes[:fillcolor]
            plotattributes[:linecolor] = plotattributes[:markerstrokecolor]
        end
        primary := true
        seriestype := :shape
        x := xsegs.pts
        y := ysegs.pts
        ()
    end

    # Outliers
    if outliers && !isempty(outliers)
        @series begin
            primary := false
            seriestype := :scatter
            if get!(plotattributes, :markershape, :circle) == :none
                plotattributes[:markershape] = :circle
            end

            fillrange := nothing
            x := outliers_x
            y := outliers_y
            ()
        end
    end

    # Hover
    primary := false
    seriestype := :path
    marker := false
    hover := texts
    linewidth := 0
    x := xsegs.pts
    y := ysegs.pts
    ()
end

@userplot GroupedBoksplot

recipetype(::Val{:groupedboksplot}, args...) = GroupedBoksplot(args)

@recipe function f(g::GroupedBoksplot; spacing = 0.1)
    x, y = grouped_xy(g.args...)

    # extract xnums and set default bar width.
    # might need to set xticks as well
    ux = unique(x)
    x = if eltype(x) <: Number
        bar_width --> (0.8 * mean(diff(sort(ux))))
        float.(x)
    else
        bar_width --> 0.8
        xnums = [findfirst(isequal(xi), ux) for xi in x] .- 0.5
        xticks --> (eachindex(ux) .- 0.5, ux)
        xnums
    end

    # shift x values for each group
    group = get(plotattributes, :group, nothing)
    if group != nothing
        gb = RecipesPipeline._extract_group_attributes(group)
        labels, idxs = getfield(gb, 1), getfield(gb, 2)
        n = length(labels)
        bws = plotattributes[:bar_width] / n
        bar_width := bws * clamp(1 - spacing, 0, 1)
        for i in 1:n
            groupinds = idxs[i]
            Δx = _cycle(bws, i) * (i - (n + 1) / 2)
            x[groupinds] .+= Δx
        end
    end

    seriestype := :boksplot
    x, y
end


"""
Generate an empty plot object for padding.
$(TYPEDSIGNATURES)
"""
function pad(;kwargs...)
    return plot(;grid=false, axis=false, xticks=[], yticks=[], xaxis=false, kwargs...)
end

"""
Add left and bottom padding to existing plot. E.g., `plot(x,y) |> pad`.
$(TYPEDSIGNATURES)
"""
function pad(p::Plots.Plot; kwargs...)
    return plot(pad(), p, pad(;kwargs...), layout=@layout([a{0.005w} b{0.9w,0.995h}; c{0.005h}]))
end

"""
Plot parameter probabilities as histogram, using weights `w`.
$(TYPEDSIGNATURES)
"""
function probplot!(x, w; kwargs...)
    histogram!(x; lw=0, fillalpha=0.5, weights=Weights(w), normalize=:pdf, kwargs...)
end

"""
Generate legend for boksplot default settings.
$(TYPEDSIGNATURES)
"""
function mk_bksplt_leg(;annfontsize=22)

    x = rand(Normal(),1_000)
    global bksplt_leg = plot(
        repeat([1], length(x)),
        x, 
        seriestype=:boksplot,
        q_whiskers=(.05, .95),
        leg=false, 
        fillalpha=.25, 
        outliers=false,
        xaxis=false, yaxis=false, 
        xticks=[], yticks=[],
        size=(400,400), 
        rightmargin=45mm,
        bar_width=0.25, lw=2
        )
    plot!([(1.3,mean(x)), (1.2,mean(x))], arrow = arrow(:closed, 0.1), color=:black, lw=2)
    plot!([(1.3,quantile(x, 0.25)), (1.2,quantile(x, 0.25))], arrow = arrow(:closed, 0.001), lw=2, color=palette(:default)[1])
    plot!([(1.3,quantile(x, 0.75)), (1.2,quantile(x, 0.75))], arrow = arrow(:closed, 0.001), lw=2, color=palette(:default)[1])
    plot!([(1.3,quantile(x, 0.05)), (1.2,quantile(x, 0.05))], arrow = arrow(:closed, 0.001), lw=2, color=palette(:default)[2])
    plot!([(1.3,quantile(x, 0.95)), (1.2,quantile(x, 0.95))], arrow = arrow(:closed, 0.001), lw=2, color=palette(:default)[2])
    annotate!(1.525, mean(x), Plots.text("Mean", annfontsize))
    annotate!(1.525, quantile(x, .25), Plots.text("25th %ile", annfontsize), color=palette(:default)[1])
    annotate!(1.525, quantile(x, .75), Plots.text("75th %ile", annfontsize), color=palette(:default)[1])
    annotate!(1.525, quantile(x, .05), Plots.text("5th %ile", annfontsize), color=palette(:default)[2])
    annotate!(1.525, quantile(x, .95), Plots.text("95th %ile", annfontsize), color=palette(:default)[2])
    
end

"""
Legend for grouped violinplot. 
$(TYPEDSIGNATURES)
"""
function violinleg()
    n = 1_000_000
    p = violin(
        sample([1], n), rand(Normal(), n), side=:left, color=3, lw=0, fillalpha=.5, 
        xticks=[], yticks=[], grid=false, xaxis=false, yaxis=false, 
        title="Treamtent level \n\n  low    |    high"
        )
    violin!(
        sample([1], n), rand(Normal(), n), side=:right, color=4, lw=0, fillalpha=.5, 
    )
    return p
end