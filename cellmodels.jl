# Calculators for battery cell energy density.

# Using Measurements for the error/uncertainty handling and Parameters for the
# fancier structs. QuadGK for integration.
using Measurements, Parameters, QuadGK, DataFrames

### Material classes!

# Active material. Has a name, a specific capacity and an average potential.
# For now. Might add some more properties later.
# Values are unions of numbers and measurements (they can be numbers without
# uncertainties, or they can have uncertainties.)
@with_kw struct ActiveMaterial
    name::String
    spec_cap::Union{Number, Measurement} # mAh g-1
    avg_E::Union{Number, Measurement}
end

# Current collector. Has a name. For Al or Cu, can just give name as "Al" or "Cu"
# along with thickness, and density will be automatically included.
@with_kw struct CurrentCollector
    name::String
    thickness::Union{Number, Measurement} # cm
    density::Union{Nothing, Number, Measurement} = nothing # g/cm3

    function CurrentCollector(name, thickness)
        if name == "Al"
            density = 2.7
            new(name, thickness, density)
        elseif name == "Cu"
            density = 8.96
            new(name, thickness, density)
        end
    end

    # Mass is included apparently? I think I want this
    mass::Union{Number, Measurement} = thickness * density
end

# Electrode. This is kinda complicated.
# Start: electrode is defined as an active material composite, coated on both sides
# of a current collector.

# First need to define the composite, which has a thickness, an areal capacity (mAh/cm2), an active loading (mg/cm2), an active fraction (%), and a composite density.
@with_kw struct ElectrodeComposite
    active_material::ActiveMaterial
#    cc::CurrentCollector
    thickness::Union{Nothing, Number, Measurement} = nothing ## cm
    areal_cap::Union{Nothing, Number, Measurement} = nothing ## mAh cm-2
    active_load::Union{Nothing, Number, Measurement} = nothing ## mg cm-2
    active_frac::Union{Nothing, Number, Measurement} # dimensionless
    density::Union{Nothing, Number, Measurement} = nothing # g cm-3

    # There are two ways to define the electrode.
    function ElectrodeComposite(active_material, thickness, areal_cap,
        active_load, active_frac, density)
        # First, if you know areal capacity, active fraction and density, you can leave
        # thickness and active_load empty, and these will be calculated.
        if thickness == nothing && active_load == nothing
            thickness = areal_cap / (active_frac * active_material.spec_cap * density)
            active_load = areal_cap / active_material.spec_cap
            new(active_material, thickness, areal_cap, active_load, active_frac, density)
        # Otherwise, if you know thickness, loading and active fraction, you can calculate areal capacity and
        # composite density.
        elseif areal_cap == nothing && density == nothing
            areal_cap = active_material.spec_cap * active_load
            density = (active_load / active_frac) / thickness
            new(active_material, thickness, areal_cap, active_load, active_frac, density)
        else
            error("Information given does not seem to fit existing methods :(")
        end
    end
end

# An electrode is made up of a composite and a current collector.
@with_kw struct Electrode
    composite::ElectrodeComposite
    cc::CurrentCollector
end

# Separator: has a name, thickness, porosity, and density.
@with_kw struct Separator
    name::String
    thickness::Union{Number, Measurement}
    porosity::Union{Number, Measurement}
    density::Union{Number, Measurement}
end

# Electrolyte: has a salt, solvent, concentration, density.
@with_kw struct Electrolyte
    salt::String
    solvent::String
    concentration::Union{Nothing, Number, Measurement} = nothing

#   On second thoughts probably the vol_cap_ratio shouldn't be a property of the electrolyte.
#     It should probably be a property of the cell.
#    vol_cap_ratio::Union{Number, Measurement}
    density::Union{Nothing, Number, Measurement} = nothing
    saltmassfrac::Union{Nothing, Number, Measurement} = nothing

    function Electrolyte(salt, solvent, concentration, density, saltmassfrac)
        if salt == "LiPF6" && concentration != nothing
            saltmassfrac = 0.1222 * concentration
            if density == nothing && solvent == "EC:DEC"
                density = 0.7641 * saltmassfrac + 1.1299
            end
        end
        new(salt, solvent, concentration, density, saltmassfrac)
    end

end

## Cell models

# Cylindrical cell, which is just a spiral stack constrained in a cylinder.
@with_kw struct CylindricalCell

    # Has a name
    name::String

    # Has the following materials.
    positive::Electrode
    negative::Electrode
    separator::Separator
    electrolyte::Electrolyte

    # Has the additional following choice
    ecap_ratio::Union{Number, Measurement}

    # Has the following dimensions
    diameter::Union{Number, Measurement}
    height::Union{Number, Measurement}

    # And the can and cell geometry
    canthickness::Union{Number, Measurement}
    candensity::Union{Number, Measurement}
    # voiddiameter is the diameter of the centre of the spiral from where the spiral starts
    voiddiameter::Union{Number, Measurement}
    headspace::Union{Number, Measurement}
    extramass::Union{Number, Measurement}

    # These give the jellyroll the following properties
    jr_area::Union{Nothing, Number, Measurement} = nothing

    # And the cell some key properties
    capacity::Union{Nothing, Number, Measurement} = nothing
    energy::Union{Nothing, Number, Measurement} = nothing

    # Input is everything but jellyroll area, capacity, energy; which will be calculated.
    function CylindricalCell(name, positive, negative, separator, electrolyte,
        ecap_ratio, diameter, height, canthickness, candensity, voiddiameter,
        headspace, extramass, jr_area, capacity, energy)

        # First we calculate the stack thickness
        stackthickness = thickness(positive, negative, separator)

        # Then we calculate the number of turns in the jellyroll
        turns = ((diameter - (2 * canthickness) - stackthickness - voiddiameter) / 2) / stackthickness

        # Then we get our jellyroll length from the Archimedes spiral equation, and area from that times the height
        l, err = quadgk(θ -> sqrt(((voiddiameter / 2) + stackthickness * θ/2π)^2 + (stackthickness / 2π)^2), 0, turns * 2π)
        jr_area = l * (height - headspace - 2 * canthickness)

        # Now capacity of the cell. It's the minimum of the two areal capacities of the electrode, times the jellyroll area, * 2
        # Then / 1000 to give answer in Ah.
        capacity = min(positive.composite.areal_cap, negative.composite.areal_cap) * jr_area * 2 / 1000

        # Energy is capacity times the potential difference
        energy = (positive.composite.active_material.avg_E - negative.composite.active_material.avg_E) * capacity

        # Now build the cell.
        new(name, positive, negative, separator, electrolyte, ecap_ratio,
        diameter, height, canthickness, candensity, voiddiameter, headspace,
        extramass, jr_area, capacity, energy)
    end
end

# Prismatic cell. Jelly roll is wound vertically (i.e. spiral axis is parallel to top plate and bus bars are at the sides.)
@with_kw struct PrismaticCell_Jellyroll

    # Has a name
    name::String

    # Has the following materials.
    positive::Electrode
    negative::Electrode
    separator::Separator
    electrolyte::Electrolyte

    # Has the additional following choice
    ecap_ratio::Union{Number, Measurement}

    # Has the following dimensions
    height::Union{Number, Measurement}
    width::Union{Number, Measurement}
    depth::Union{Number, Measurement}

    # And the can and cell geometry
    canthickness::Union{Number, Measurement}
    candensity::Union{Number, Measurement}
    headspace::Union{Number, Measurement}
    termclearance::Union{Number, Measurement}
    nrolls::Number
    topplatemass::Union{Number, Measurement}
    extramass::Union{Number, Measurement}

    # These give the jellyroll the following properties
    jr_area::Union{Nothing, Number, Measurement} = nothing

    # And the cell some key properties
    capacity::Union{Nothing, Number, Measurement} = nothing
    energy::Union{Nothing, Number, Measurement} = nothing

    function PrismaticCell_Jellyroll(name, positive, negative, separator,
        electrolyte, ecap_ratio, height, width, depth, canthickness,
        candensity, headspace, termclearance, nrolls, topplatemass,
        extramass, jr_area, capacity, energy)

        # First we calculate the stack thickness
        stackthickness = thickness(positive, negative, separator)

        # The jellyroll radius is calculated as follows. There are nrolls in the cell.
        # The 'box' that nrolls rolls have to fit in is depth minus twice the can thickness, divided by the number of rolls.
        rollthickness = (depth - (2 * canthickness)) / nrolls

        # Then the turns is the roll thickness minus one stack thickness, divided by two, then divided by the stack thickness.
        turns = ((rollthickness - stackthickness) / 2) / stackthickness

        # Then we get our "turn length" from the Archimedes spiral equation, starting from zero.
        turnl, err = quadgk(θ -> sqrt((stackthickness * θ/2π)^2 + (stackthickness / 2π)^2), 0, turns * 2π)

        # Jelly roll area is the length of the spiral, plus the flat parts, times the width and the number of rolls.
        jr_area = (turnl + ((height - headspace - rollthickness) * 2 * floor(turns))) * (width - (termclearance * 2)) * nrolls

        # Now capacity of the cell. It's the minimum of the two areal capacities of the electrode, times the jellyroll area, * 2
        # Then / 1000 to give answer in Ah.
        capacity = min(positive.composite.areal_cap, negative.composite.areal_cap) * jr_area * 2 / 1000

        # Energy is capacity times the potential difference
        energy = (positive.composite.active_material.avg_E - negative.composite.active_material.avg_E) * capacity

        # Now build the cell.
        new(name, positive, negative, separator, electrolyte, ecap_ratio,
        height, width, depth, canthickness, candensity, headspace,
        termclearance, nrolls, topplatemass, extramass, jr_area, capacity,
        energy)
    end
end

# PrismaticCell - "Cinnamon roll" is wound "horizontally", i.e., spiral axis is perpendicular to top plate and bus bars sit on the top.
@with_kw struct PrismaticCell_Cinnamonroll

    # Has a name
    name::String

    # Has the following materials.
    positive::Electrode
    negative::Electrode
    separator::Separator
    electrolyte::Electrolyte

    # Has the additional following choice
    ecap_ratio::Union{Number, Measurement}

    # Has the following dimensions
    height::Union{Number, Measurement}
    width::Union{Number, Measurement}
    depth::Union{Number, Measurement}

    # And the can and cell geometry
    canthickness::Union{Number, Measurement}
    candensity::Union{Number, Measurement}
#    headspace::Union{Number, Measurement} # No headspace, just term clearance.
    termclearance::Union{Number, Measurement}
    nrolls::Number
    topplatemass::Union{Number, Measurement}
    extramass::Union{Number, Measurement}

    # These give the jellyroll the following properties
    jr_area::Union{Nothing, Number, Measurement} = nothing

    # And the cell some key properties
    capacity::Union{Nothing, Number, Measurement} = nothing
    energy::Union{Nothing, Number, Measurement} = nothing

    function PrismaticCell_Cinnamonroll(name, positive, negative, separator,
        electrolyte, ecap_ratio, height, width, depth, canthickness,
        candensity, termclearance, nrolls, topplatemass, extramass,
        jr_area, capacity, energy)

        # First we calculate the stack thickness
        stackthickness = thickness(positive, negative, separator)

        # The jellyroll radius is calculated as follows. There are nrolls in the cell.
        # The 'box' that nrolls rolls have to fit in is depth minus twice the can thickness, divided by the number of rolls.
        rollthickness = (depth - (2 * canthickness)) / nrolls

        # Then the turns is the roll thickness minus one stack thickness, divided by two, then divided by the stack thickness.
        turns = ((rollthickness - stackthickness) / 2) / stackthickness

        # Then we get our "turn length" from the Archimedes spiral equation, starting from zero.
        turnl, err = quadgk(θ -> sqrt((stackthickness * θ/2π)^2 + (stackthickness / 2π)^2), 0, turns * 2π)

        # Jelly roll area is the length of the spiral, plus the flat parts, times the width and the number of rolls.
        jr_area = (turnl + ((width - rollthickness) * 2 * floor(turns))) * (height - termclearance) * nrolls

        # Now capacity of the cell. It's the minimum of the two areal capacities of the electrode, times the jellyroll area, * 2
        # Then / 1000 to give answer in Ah.
        capacity = min(positive.composite.areal_cap, negative.composite.areal_cap) * jr_area * 2 / 1000

        # Energy is capacity times the potential difference
        energy = (positive.composite.active_material.avg_E - negative.composite.active_material.avg_E) * capacity

        # Now build the cell.
        new(name, positive, negative, separator, electrolyte, ecap_ratio,
        height, width, depth, canthickness, candensity, termclearance,
        nrolls, topplatemass, extramass, jr_area, capacity, energy)
    end
end

# Pouch cell. Rectangular stacks of electrodes. Not constrained by can, so number of layers is defined. Tabs are added.
@with_kw struct PouchCell

    # Has a name
    name::String

    # Has the following materials.
    positive::Electrode
    negative::Electrode
    separator::Separator
    electrolyte::Electrolyte

    # Has the additional following choice
    ecap_ratio::Union{Number, Measurement}

    # Has the following dimensions
    height::Union{Number, Measurement}
    width::Union{Number, Measurement}
#    depth::Union{Number, Measurement}
    nlayers::Union{Number, Measurement}

    # And the pouch and cell geometry
    pouchthickness::Union{Number, Measurement}
    pouchdensity::Union{Number, Measurement}
#    headspace::Union{Number, Measurement} # No headspace, just term clearance.
    pouchclearance::Union{Number, Measurement}

    termh::Union{Number, Measurement}
    termw::Union{Number, Measurement}
    termt::Union{Number, Measurement} = 0.05

    termdenspos::Union{Number, Measurement} = 2.7
    termdensneg::Union{Number, Measurement} = 8.9

    extramass::Union{Number, Measurement}

    # These give the jellyroll the following properties
    jr_area::Union{Nothing, Number, Measurement} = nothing

    # And the cell some key properties
    capacity::Union{Nothing, Number, Measurement} = nothing
    energy::Union{Nothing, Number, Measurement} = nothing

    function PouchCell(name, positive, negative, separator, electrolyte, ecap_ratio, height, width, nlayers, pouchthickness, pouchdensity, pouchclearance, termh, termw, termt, termdenspos, termdensneg, extramass, jr_area, capacity, energy)

        # Total stack thickness is the stack thickness * nlayers, plus extra negative.
        #tstackthickness = thickness(positive, negative, separator) + thickness(negative)

        # Total active area is nlayers * w* h.
        jr_area = nlayers * width * height

        # Capacity
        capacity = min(positive.composite.areal_cap, negative.composite.areal_cap) * jr_area * 2 / 1000

        # Energy is capacity times the potential difference
        energy = (positive.composite.active_material.avg_E - negative.composite.active_material.avg_E) * capacity

        new(name, positive, negative, separator, electrolyte, ecap_ratio, height, width, nlayers, pouchthickness, pouchdensity, pouchclearance, termh, termw, termt, termdenspos, termdensneg, extramass, jr_area, capacity, energy)

    end
end

## Overwriting show() for these cell types. Made one for each type because I don't know how to make it accept all four in the same function...

function Base.show(io::IO, cell::CylindricalCell)
    println(io, cell.name, "\n", cell.positive.composite.active_material.name, " cathode @ ", cell.positive.composite.areal_cap, " mAh/cm2\n",
    cell.capacity, " Ah, ", round(cell.energy, digits = 2), " Wh\n", mass(cell), " g\n", gravimetric_energy(cell), " Wh/kg\n" ,volumetric_energy(cell), " Wh/L")
end

function Base.show(io::IO, cell::PouchCell)
    println(io, cell.name, "\n", cell.positive.composite.active_material.name, " cathode @ ", cell.positive.composite.areal_cap, " mAh/cm2\n",
    cell.capacity, " Ah, ", round(cell.energy, digits = 2), " Wh\n", mass(cell), " g\n", gravimetric_energy(cell), " Wh/kg\n" ,volumetric_energy(cell), " Wh/L")
end

function Base.show(io::IO, cell::PrismaticCell_Jellyroll)
    println(io, cell.name, "\n", cell.positive.composite.active_material.name, " cathode @ ", cell.positive.composite.areal_cap, " mAh/cm2\n",
    cell.capacity, " Ah, ", round(cell.energy, digits = 2), " Wh\n", mass(cell), " g\n", gravimetric_energy(cell), " Wh/kg\n" ,volumetric_energy(cell), " Wh/L")
end

function Base.show(io::IO, cell::PrismaticCell_Cinnamonroll)
    println(io, cell.name, "\n", cell.positive.composite.active_material.name, " cathode @ ", cell.positive.composite.areal_cap, " mAh/cm2\n",
    cell.capacity, " Ah, ", round(cell.energy, digits = 2), " Wh\n", mass(cell), " g\n", gravimetric_energy(cell), " Wh/kg\n" ,volumetric_energy(cell), " Wh/L")
end

## Thickness calculation methods

# Thickness of a composite. Simple because it is already specified in the class
function thickness(composite::ElectrodeComposite)
    composite.thickness
end

# Thickness of a separator. Same as for composite.
function thickness(separator::Separator)
    separator.thickness
end

# Current collector, same again.
function thickness(cc::CurrentCollector)
    cc.thickness
end

# Thickness of an electrode. Twice the composite thickness + thickness of the current collector.
function thickness(electrode::Electrode)
    (2 * electrode.composite.thickness) + electrode.cc.thickness
end

# Stack thickness. Two current collectors, twice each of positive and negative composite, twice separator thickness.
function thickness(positive::Electrode, negative::Electrode, separator::Separator)
    positive.cc.thickness + negative.cc.thickness + (2 * positive.composite.thickness) + (2 * negative.composite.thickness) + (2 * separator.thickness)
end

# Thickness of a pouch cell. Pouch cell is not constrained by can, so thickness is given by the number of layers + the pouch itself.
function thickness(cell::PouchCell)

    # Thickness is nlayers * stackthickness, plus 1 extra negative electrode, plus twice the pouch thickness itself.
    return (cell.nlayers * thickness(cell.positive, cell.negative, cell.separator)) + thickness(cell.negative) + (2 * cell.pouchthickness)

end



## Mass calculation methods

# Composite. Just active loading divided by active fraction.
function mass(composite::ElectrodeComposite)
    composite.active_load / composite.active_frac
end

# Mass of a separator
function mass(separator::Separator)
    separator.thickness * separator.density
end

# Mass of a current collector
function mass(cc::CurrentCollector)
    cc.thickness * cc.density
end

# Mass of an electrode. Twice composite plus current collector
function mass(electrode::Electrode)
    (2 * mass(electrode.composite)) + mass(electrode.cc)
end

# Mass of the stack. Two electrodes + two separators.
function mass(positive::Electrode, negative::Electrode, separator::Separator)
    mass(positive) + mass(negative) + (2 * mass(separator))
end

# Mass of cylindrical cell.
function mass(cell::CylindricalCell)
    # Mass of a cell is the mass of the can plus the mass of the jellyroll.
    # volume of the can itself is:
    canvol = (((π * (cell.diameter/2)^2) - (π * (cell.diameter/2 - cell.canthickness)^2)) * cell.height) +
        (2π * (cell.diameter/2)^2 * cell.canthickness)

    # Can weight is given by:
    canwt = (canvol * cell.candensity) + cell.extramass

    #Jellyroll mass is the stack mass times the area.
    jr_mass = mass(cell.positive, cell.negative, cell.separator) * cell.jr_area

    # Now add the electrolyte mass.
    e_mass = cell.ecap_ratio * cell.capacity  * cell.electrolyte.density

    # Return everything
    return canwt + jr_mass + e_mass
end

function massbreakdown(cell::CylindricalCell)

    # Initialise the data frame
    res = DataFrame(component = String[], mass = Number[], percentage = Measurement[])

    # Positive current collector
    poscc = mass(cell.positive.cc) * cell.jr_area
    push!(res, ["+ve cc", poscc, 100 * poscc / mass(cell)])

    # Positive composite
    poscomp = mass(cell.positive.composite) * cell.jr_area * 2
    push!(res, ["+ve composite", poscomp, 100 * poscomp / mass(cell)])

    # Negative current collector
    negcc = mass(cell.negative.cc) * cell.jr_area
    push!(res, ["-ve cc", negcc, 100 * negcc / mass(cell)])

    # Negative composite
    negcomp = mass(cell.negative.composite) * cell.jr_area * 2
    push!(res, ["-ve composite", negcomp, 100 * negcomp / mass(cell)])

    # Separator
    sep = mass(cell.separator) * cell.jr_area * 2
    push!(res, ["separator", sep, 100 * sep / mass(cell)])

    # Electrolyte
    elyte = cell.ecap_ratio * cell.capacity * cell.electrolyte.density
    push!(res, ["electrolyte", elyte, 100 * elyte / mass(cell)])

    # Packaging; from the mass() function
    canvol = (((π * (cell.diameter/2)^2) - (π * (cell.diameter/2 - cell.canthickness)^2)) * cell.height) +
        (2π * (cell.diameter/2)^2 * cell.canthickness)

    # Can weight is given by:
    canwt = (canvol * cell.candensity) + cell.extramass
    push!(res, ["packaging", canwt, 100 * canwt / mass(cell)])

    # and extramass
#    push!(res, ["est. extra mass", cell.extramass, 100 * cell.extramass / mass(cell)])

    return res

end

# Mass of prismatic cells - method is the same for each type.
function mass(cell::Union{PrismaticCell_Jellyroll, PrismaticCell_Cinnamonroll})
    # Mass of the cell is the mass of the can plus the mass of the jellyroll.
    # Volume of the can is:
    canvol = ((cell.width * cell.depth) + ((cell.height - 2 * cell.canthickness) * (cell.width - 2 * cell.canthickness) * 2) + ((cell.height - 2 * cell.canthickness) * cell.depth * 2)) * cell.canthickness

    # Can plus top plate is
    canwt = (canvol * cell.candensity) + cell.topplatemass + cell.extramass

    #Jellyroll mass is the stack mass times the area.
    jr_mass = mass(cell.positive, cell.negative, cell.separator) * cell.jr_area

    # Now add the electrolyte mass.
    e_mass = cell.ecap_ratio * cell.capacity  * cell.electrolyte.density

    # Return everything
    return canwt + jr_mass + e_mass
end

# Mass of a pouch cell.
function mass(cell::PouchCell)

    # First we need the mass of the "jellyroll"
    jr_mass = (mass(cell.positive, cell.negative, cell.separator) * cell.jr_area) +
    (mass(cell.negative) * cell.width * cell.height)

    # Now add the electrolyte mass.
    e_mass = cell.ecap_ratio * cell.capacity  * cell.electrolyte.density

    # Then we need the mass of the pouch

    pouchmass = cell.pouchthickness * (cell.width + cell.pouchclearance) * (cell.height + cell.pouchclearance) * cell.pouchdensity * 2

    termmass = ((cell.termh + cell.pouchclearance) * cell.termw * cell.termt * cell.termdenspos) + ((cell.termh + cell.pouchclearance) * cell.termw * cell.termt * cell.termdensneg)

    return jr_mass + e_mass + pouchmass + termmass + cell.extramass

end

function massbreakdown(cell::PouchCell)

    # Initialise the data frame
    res = DataFrame(component = String[], mass = Number[], percentage = Measurement[])

    # Positive current collector
    poscc = mass(cell.positive.cc) * cell.jr_area
    push!(res, ["+ve cc", poscc, 100 * poscc / mass(cell)])

    # Positive composite
    poscomp = mass(cell.positive.composite) * cell.jr_area * 2
    push!(res, ["+ve composite", poscomp, 100 * poscomp / mass(cell)])

    # Negative current collector
    negcc = mass(cell.negative.cc) * cell.jr_area + (mass(cell.negative.cc) * cell.width * cell.height)
    push!(res, ["-ve cc", negcc, 100 * negcc / mass(cell)])

    # Negative composite
    negcomp = mass(cell.negative.composite) * cell.jr_area * 2 + (mass(cell.negative.composite) * cell.width * cell.height * 2)
    push!(res, ["-ve composite", negcomp, 100 * negcomp / mass(cell)])

    # Separator
    sep = mass(cell.separator) * cell.jr_area * 2
    push!(res, ["separator", sep, 100 * sep / mass(cell)])

    # Electrolyte
    elyte = cell.ecap_ratio * cell.capacity * cell.electrolyte.density
    push!(res, ["electrolyte", elyte, 100 * elyte / mass(cell)])

    # Packaging
    pouchmass = cell.pouchthickness * (cell.width + cell.pouchclearance) * (cell.height + cell.pouchclearance) * cell.pouchdensity * 2
    termmass = ((cell.termh + cell.pouchclearance) * cell.termw * cell.termt * cell.termdenspos) + ((cell.termh + cell.pouchclearance) * cell.termw * cell.termt * cell.termdensneg)
    packaging = pouchmass + termmass + cell.extramass
    push!(res, ["packaging", packaging, 100 * packaging / mass(cell)])

    # and extramass
    #push!(res, ["est. extra mass", cell.extramass, 100 * cell.extramass / mass(cell)])

    return res

end
## Energy density calculation methods

function gravimetric_energy(cell::Union{CylindricalCell, PrismaticCell_Jellyroll, PrismaticCell_Cinnamonroll, PouchCell})
    # It's just the energy divided by the mass, converted to Wh/kg.
    return 1000 * cell.energy / mass(cell)
end

function volumetric_energy(cell::CylindricalCell)
    # Need the volume of the cell
    cellvol = π * (cell.diameter / 2)^2 * cell.height # cm^3
    return 1000 * cell.energy / cellvol # Wh/L
end

function volumetric_energy(cell::Union{PrismaticCell_Jellyroll, PrismaticCell_Cinnamonroll})
    cellvol = cell.width * cell.height * cell.depth
    return 1000 * cell.energy / cellvol
end

function volumetric_energy(cell::PouchCell)
    cellvol = (cell.width + cell.pouchclearance) * (cell.height + cell.pouchclearance + cell.termh) * thickness(cell)

    return 1000 * cell.energy / cellvol
end
