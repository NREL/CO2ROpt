

#OVERVIEW
#----------------
#Any JuMP model that describes an optimization problem must have four parts:

#Model Object,
#Variables,
#Objective,
#Constraints.

#PACKAGE IMPORT
#----------------
println("Initializing...")

using JuMP  # Need to say it whenever we use JuMP

using GLPK

using CSV # Import data

using DataFrames # For results export to CSV

using Plots

#=
electrolyzer_tea.jl is a function used to evaluate the TEA of the electrochemical
system, including separations, stack cost, and other operational parameters.
This function takes a series of inputs specific to an electrochemical pathway,
as shown below:

function electrolyzer_TEA(
    flow, product flow in [kg/year]
    product_phase, phase of product ["Gas" or "Liquid"]
    product_mw, molecular weight of product [g/mol]
    product_density, product density [kg/m^3]
    electrons_product, [# e-]
    electrons_co2, [#e-/co2]
    system_potential, whole cell potential [V]
    selectivity, selectivity, or faradaic efficiency of CO2R system to product [%]
    hours_year, hours in a given year [hours/year]
    current_density, [A/cm^2]
    single_pass_conversion, single pass conversion of CO2 to product [%]
    electrolyzer_cost, electrolyzer capital cost [$/m^2]
    electricity_price, [$/kWh]
    water_price, process water cost [$/kg]
    water_rate, rate of water consusumption [kg water/kg product]
    crf capital recovery factor [1/years]
    )
=#
include("./electrolyzer_tea.jl")

#=
The electrolyzer tea function operates identically to the electrolyzer tea function
discussed above, however it returns an array of results for costs rather
than a single value. This allows for a more detailed analysis of costs
included in the TEA. Input values are identical, output values are as follows:

    return:
    electrolyzer_capex + electrolyzer_BOP, electrolyzer costs [$/year]
    PSA_capex + PSA_opex + distil_capex + distil_opex, separations costs [$/year]
    battery_capex, battery energy storage costs [$/year]
    maintenance_cost, system maintenance costs [$/year]
    water_opex annual process water costs [$/year]
=#
include("./electrolyzer_tea_array.jl")


__precompile__()

# Define function to perform optimization routine - decreases runtime
function co2_optimizer()

    #DATA IMPORT
    #--------------------
    input_filename = "input_data_full.csv"
    input_filepath = joinpath(@__DIR__, input_filename)

    # Import point source and point sink data
    data = DataFrame!(CSV.File(input_filepath))



    #=
    Option to exclude specific point types from either sink or
    source data in optimizer. To use, fill the array with the proper name
    of the point source to be excluded, e.g. ["COAL", "GAS"]
    =#
    #--------------------
    source_exclusions =  ["COAL","GAS"]#["None"]#["None"]#["None"]
    sink_exclusions = ["None"]#["Import","Export","Acetic Acid Plant", "Petroleum Terminal"]#["Ethane Crackers"]#["Import","Export"]#["None"]

    for i = 1:length(source_exclusions)

        data = data[data[!,8].!= source_exclusions[i] ,:]

    end

    for i = 1:length(sink_exclusions)

        data = data[data[!,8].!= sink_exclusions[i] ,:]


    end


    println("Source and sink points to be considered:")
    println(size(data,1))





    #VARIABLES and generic functions
    #---------

    # A variable is modeled using @defVar(name of the model object, variable name and bound, variable type)
    # Bound can be lower bound, upper bound or both. If no variable type is defined, then it is treated as
    #real. For binary variable write Bin and for integer use Int.


    # Sink/source point datasets
    #_____________________________________________________________________________

    #= Isolate consumption point data by including only points where
    demand for the molecule is not equal to zero =#
    sinks_C2H4 = data[data[:,4].!= 0  ,:]
    sinks_HCOO = data[data[:,5].!= 0  ,:]
    sinks_CO = data[data[:,6].!= 0  ,:]

    #= Isolate CO2 point source data by including only points where
    feedstock supply is not equal to zero =#
    source_pts = data[data[:,7].!= 0  ,:]
    sink_pts = data[data[:,7].== 0  ,:]

    # Number of chemical products considered
    num_species = 3

    # Find number of CO2 point sources
    num_source_pts = size(data[data[:,7].!= 0  ,:],1)

    # Find number of sink points
    num_sink_pts = size(data[data[:,7].== 0  ,:],1)

    # Find total number of points
    num_pts = size(data,1)

    # Faradays constant
    F = 96485 # [C/mol]

    # Sink/source distance evaluation datasets
    # Goal: build an array showing distance between sink and source points
    #_____________________________________________________________________________

    # Array of lat/long for sources
    lat_source = source_pts[:, 2]
    long_source = source_pts[:, 3]

    # Array of lat/long for sinks
    lat_sink = sink_pts[:, 2]
    long_sink = sink_pts[:, 3]


    # Molecular weights
    #_____________________________________________________________________________
    C_mw = 12.011 # [g/mol]
    CO2_mw = 44 # [g/mol]
    C2H4_mw = 28.05 # [g/mol]
    HCOO_mw = 45.017 # [g/mol]
    CO_mw = 28.01 # [g/mol]
    H2O_mw = 18 # [g/mol]

    product_mw = [C2H4_mw , HCOO_mw, CO_mw]


    # Molecule carbons
    #_____________________________________________________________________________
    C2H4_carbons = 2 # [carbons/molecule]
    HCOO_carbons = 1 # [carbons/molecule]
    CO_carbons = 1 # [carbons/molecule]

    carbons = [C2H4_carbons , HCOO_carbons, CO_carbons]


    # Electrons per mol product
    #_____________________________________________________________________________
    C2H4_electrons_product = 12  # [mol e- / mol product]
    HCOO_electrons_product = 2 # [mol e- / mol product]
    CO_electrons_product = 2 # [mol e- / mol product]
    H2O_electrons = 4 # [mol e- / mol product]

    electrons_product = [C2H4_electrons_product , HCOO_electrons_product, CO_electrons_product]

    # Product state
    #_____________________________________________________________________________
    C2H4_phase = "Gas"
    HCOO_phase = "Liquid"
    CO_phase = "Gas"

    product_phase = [C2H4_phase, HCOO_phase, CO_phase]

    # Product densities
    #_____________________________________________________________________________
    C2H4_density = 1.18 # [kg/m^3]
    HCOO_density = 1221 # [kg/m^3]
    CO_density = 1.14 # [kg/m^3]

    product_density = [C2H4_density , HCOO_density, CO_density]

    # Electrons per mol CO2
    #_____________________________________________________________________________
    C2H4_electrons_co2 = 6  # [mol e- / mol co2]
    HCOO_electrons_co2 = 2 # [mol e- / mol co2]
    CO_electrons_co2 = 2 # [mol e- / mol co2]

    electrons_co2 = [C2H4_electrons_co2 , HCOO_electrons_co2, CO_electrons_co2]

    # Current densities
    #_____________________________________________________________________________
    C2H4_cd = 1.55 # [A/cm^2] 10.1126/science.aay4217
    HCOO_cd = 0.5 # [A/cm^2] 10.1021/acsenergylett.0c00860
    CO_cd = 0.35 # [A/cm^2] http://xlink.rsc.org/?DOI=C9EE02410G

    current_density = [C2H4_cd , HCOO_cd, CO_cd]

    # single pass conversion
    #_____________________________________________________________________________
    C2H4_spc = 0.5 # [%] assumed
    HCOO_spc = 0.5 # [%] assumed
    CO_spc = 0.5 # [%] assumed

    single_pass_conversion = [C2H4_spc , HCOO_spc, CO_spc]


    # Selectivities
    #___________________________________________________________________________
    C2H4_selectivity = 0.6 # 10.1126/science.aay4217
    HCOO_selectivity = 0.93 #  10.1021/acsenergylett.0c00860
    CO_selectivity = 0.95  #  http://xlink.rsc.org/?DOI=C9EE02410G

    selectivity = [C2H4_selectivity , HCOO_selectivity, CO_selectivity]

    # Potentials
    #_____________________________________________________________________________
    C2H4_system_potential = 3.9 #[V] 10.1126/science.aay4217 taken from fig s33, durability study
    HCOO_system_potential = 5.8 #[V] 10.1021/acsenergylett.0c00860 table 1
    CO_system_potential = 3.0 #[V] 10.1021/acs.iecr.7b03514

    system_potential = [C2H4_system_potential , HCOO_system_potential, CO_system_potential]

    # Market Values
    #_____________________________________________________________________________
    C2H4_value = 0.555 # [$/kg] C2H4 value from PEP yearbook and: https://www.spglobal.com/platts/en/market-insights/podcasts/focus/040720-north-america-lng-exports-price-plunge
    HCOO_value = 0.735 # [$/kg] CO value https://pubs.acs.org/doi/10.1021/acs.iecr.7b03514
    CO_value = 0.6 # [$/kg] CO value https://pubs.acs.org/doi/10.1021/acs.iecr.7b03514

    value = [C2H4_value , HCOO_value, CO_value]


    # Market Demand
    #_____________________________________________________________________________
    C2H4_demand = sink_pts[:,4]  # consumption at each point in [kg/year]
    HCOO_demand = sink_pts[: ,5] # consumption at each point in [kg/year]
    CO_demand = sink_pts[: ,6] # consumption at each point in [kg/year]

    demand = hcat(C2H4_demand, HCOO_demand, CO_demand)

    # Water rate of consumption
    #_____________________________________________________________________________
    C2H4_water = 1.28 #[kg H2O/kg product]
    HCOO_water = 0.391 #[kg H2O/kg product]
    CO_water = 0.6426 #[kg H2O/kg product]

    water_rate = [C2H4_water, HCOO_water, CO_water]


    #_____________________________________________________________________________
    capacity = data[data[:,7].!= 0  ,7] * 1e3 * (C_mw / CO2_mw) # Feedstock carbon at each point in [kg C/year]

    transport_A = 1e-5 # Best fit line for CO2 transportation costs from IPCC (see Excel sheet)
    transport_B = -0.0002 # Best fit line for CO2 transportation costs from IPCC (see Excel sheet)

    #CO2 flue gas capture cost for each location [$/kg C]
    CO2_cap_cost = data[data[:,7].!= 0  ,12] / 1e3 * (C_mw / CO2_mw) #[$/kg C]

    # Spatially distinct electricity prices in [$/kWh] (sink pts)
    electricity_price = data[data[:,7].!= 0  ,15]

    # Economic assumptions
    #_____________________________________________________________________________

    hours_year = 8000 # [hours of operation annually]
    electrolyzer_cost = 10000 # [$/m^2]
    interest_rate = 0.08 # [%]
    facility_lifetime = 15 #[years]
    water_price = 8.929e-4 # [$/kg] https://www.energy.gov/sites/prod/files/2017/10/f38/water_wastewater_escalation_rate_study.pdf
    crf = (interest_rate * (1 + interest_rate)^facility_lifetime) / ((1 + interest_rate)^facility_lifetime - 1)

    # Define an array with rows = number of sources, columns = number of sinks, and sheets = number of products
    D = zeros(Float16, num_source_pts, num_sink_pts, num_species)

    # Distance constraint evaluation.
    #_____________________________________________________________________________
    #=
    Define an array with rows = source points, columns = sink points,
    sheets = number of species. Default values are zeros. The loop below changes
    the value to 1 if distance between sink/source is <1,000 km.
    This is then incorporated in the flow variable definition as an upper
    bound on the variable.
    =#
    up_bound = zeros(num_source_pts,num_sink_pts,num_species)
    for i in 1:num_source_pts
        for j in 1:num_sink_pts
            for k in 1:num_species
                pt1 = [long_source[i]; lat_source[i]]
                pt2 = [long_sink[j]; lat_sink[j]]

                D[i, j, k] = (abs(pt1[1] - pt2[1]) + abs(pt1[2] - pt2[2])) / 1e3

                if (abs(source_pts[i,3] - sink_pts[j,3]) + abs(source_pts[i,2] - sink_pts[j,2])) / 1e3 < 1000
                    up_bound[i, j, k] = 1.0

                end
            end
        end
    end

    #MODEL CONSTRUCTION
    #--------------------

    # Preparing an optimization model
    eco2RR = Model(GLPK.Optimizer)


    #OPTIMIZATION VARIABLES
    #--------------------
    #Initialize flow array [kg/year]
    @variable(eco2RR, flow[i = 1:num_source_pts , j = 1:num_sink_pts , k = 1:num_species] >= 0, upper_bound = up_bound[i, j, k].*1e10) # distance constraint
    #@variable(eco2RR, flow[i = 1:num_source_pts , j = 1:num_sink_pts , k = 1:num_species] >= 0)

    #OPTIMIZATION CONSTRAINTS
    #--------------------

    # Ship no more than plant capacity, in units of [kg carbon]
    @constraint(eco2RR, capacity_con[i in 1:num_source_pts],
    + sum(flow[i, : , 1]) * carbons[1] * C_mw / product_mw[1] / selectivity[1]
    + sum(flow[i, : , 2]) * carbons[2] * C_mw / product_mw[2] / selectivity[2]
    + sum(flow[i, : , 3]) * carbons[3] * C_mw / product_mw[3] / selectivity[3]
    <= capacity[i])

    #=
    Constrain total electrochemical system power consumption at a source point
    to less than 300MW. Provides a realistic constraint on system size

    Modifying this equality only applies when the system is feedstock constrained
    i.e. total demand is greater than carbon feedstock available

    Two key options are available:

    == equality - model required to export 300MW worth at each CO2 point.

    <= equality - model can provide anything up to 300MW system at a CO2 point source.

    =#

    @constraint(eco2RR, elec_con[i in 1:num_source_pts],
    sum(
    + sum(flow[i, :, 1]) / hours_year / 3600 * 1e3 / product_mw[1] * electrons_product[1] * F * system_potential[1] / selectivity[1] / 1e6
    + sum(flow[i, :, 2]) / hours_year / 3600 * 1e3 / product_mw[2] * electrons_product[2] * F * system_potential[2] / selectivity[2] / 1e6
    + sum(flow[i, :, 3]) / hours_year / 3600 * 1e3 / product_mw[3] * electrons_product[3] * F * system_potential[3] / selectivity[3] / 1e6
    )
    <= 500)

    #=
    Choose how the model attempts to meet demand for different molecules.
    This constraint can cause the model to be unable to find a SOLUTION
    if feedstock inputs are not sufficient to meet demand.

    Modifying this equality only applies when the system is demand constrained
    i.e. total demand is less than carbon feedstock available

    Two key options are available:

    == equality - model required to meet demand at each point. In this scenario,
        if feedstock inputs are not sufficient then optimizer will not find a
        solution

    <= equality - model can provide anything up to market demand at a point.
        Using this equality allows the optimizer to only make flows it deems
        economically sustainable (positive NPV).

    =#
    @constraint(eco2RR, demand_con[j in 1:num_sink_pts, k in 1:num_species],
    sum(flow[:, j , k]) == demand[j,k])

    println("Constraints Defined")

    #OPTIMIZATION OBJECTIVE FUNCTION
    #--------------------
    println("Beginning Optimization...")

    @objective(
        eco2RR, # Designate model name
        Max, # Tell JuMP to minimize/optimize
        sum(
        (

        - (D[i, j, k] * transport_A + transport_B) * flow[i, j, k] # Transportation cost
        - electrolyzer_TEA(
            flow[i, j, k],
            product_phase[k],
            product_mw[k],
            product_density[k],
            electrons_product[k],
            electrons_co2[k],
            system_potential[k],
            selectivity[k],
            hours_year,
            current_density[k],
            single_pass_conversion[k],
            electrolyzer_cost,
            electricity_price[i],
            water_price,
            water_rate[k],
            crf) # System TEA
        - flow[i, j, k] * carbons[k] * C_mw / product_mw[k] * CO2_cap_cost[i] # Carbon capture cost
        - flow[i, j, k] / hours_year / 3600 * 1e3 / product_mw[k] * electrons_product[k] * F * system_potential[k] / selectivity[k] / 1e3 * hours_year * electricity_price[i] # Electrochemical system electricity cost
        + flow[i, j, k] * value[k] # Product value
        )

        for
            i in 1:num_source_pts, # Iterate over source points
            j in 1: num_sink_pts, # Iterate over sink points
            k in 1: num_species#, # Iterate over number of species
    ))


    # Solving the optimization problem
    optimize!(eco2RR)


    # Set objective value
    obj_value = JuMP.objective_value(eco2RR)

    #OPTIMIZATION SOLUTION WRANGLING
    #--------------------
    println("Optimizer termination status:")
    println(termination_status(eco2RR))
    println(primal_status(eco2RR))
    println("Objective function value:")
    println(obj_value)

    #Output filepath and variable name
    #--------------------

    C2H4_result_filepath = joinpath(@__DIR__, "results/flows/C2H4_flows.csv")
    HCOO_result_filepath = joinpath(@__DIR__, "results/flows/HCOO_flows.csv")
    CO_result_filepath = joinpath(@__DIR__, "results/flows/CO_flows.csv")


    C2H4_NPV_filepath = joinpath(@__DIR__, "results/npv/C2H4_NPV.csv")
    HCOO_NPV_filepath = joinpath(@__DIR__, "results/npv/HCOO_NPV.csv")
    CO_NPV_filepath = joinpath(@__DIR__, "results/npv/CO_NPV.csv")


    C2H4_trans_filepath = joinpath(@__DIR__, "results/econ/trans/C2H4_trans.csv")
    HCOO_trans_filepath = joinpath(@__DIR__, "results/econ/trans/HCOO_trans.csv")
    CO_trans_filepath = joinpath(@__DIR__, "results/econ/trans/CO_trans.csv")


    C2H4_tea_filepath = joinpath(@__DIR__, "results/econ/tea/C2H4_tea.csv")
    HCOO_tea_filepath = joinpath(@__DIR__, "results/econ/tea/HCOO_tea.csv")
    CO_tea_filepath = joinpath(@__DIR__, "results/econ/tea/CO_tea.csv")


    C2H4_ccu_filepath = joinpath(@__DIR__, "results/econ/ccu/C2H4_ccu.csv")
    HCOO_ccu_filepath = joinpath(@__DIR__, "results/econ/ccu/HCOO_ccu.csv")
    CO_ccu_filepath = joinpath(@__DIR__, "results/econ/ccu/CO_ccu.csv")


    C2H4_elec_filepath = joinpath(@__DIR__, "results/econ/electricity/C2H4_elec.csv")
    HCOO_elec_filepath = joinpath(@__DIR__, "results/econ/electricity/HCOO_elec.csv")
    CO_elec_filepath = joinpath(@__DIR__, "results/econ/electricity/CO_elec.csv")


    C2H4_prf_filepath = joinpath(@__DIR__, "results/econ/prf/C2H4_prf.csv")
    HCOO_prf_filepath = joinpath(@__DIR__, "results/econ/prf/HCOO_prf.csv")
    CO_prf_filepath = joinpath(@__DIR__, "results/econ/prf/CO_prf.csv")

    C2H4_water_rsrc_filepath = joinpath(@__DIR__, "results/resource/water/C2H4_water_rsrc.csv")
    HCOO_water_rsrc_filepath = joinpath(@__DIR__, "results/resource/water/HCOO_water_rsrc.csv")
    CO_water_rsrc_filepath = joinpath(@__DIR__, "results/resource/water/CO_water_rsrc.csv")

    C2H4_power_filepath = joinpath(@__DIR__, "results/resource/power/C2H4_power.csv")
    HCOO_power_filepath = joinpath(@__DIR__, "results/resource/power/HCOO_power.csv")
    CO_power_filepath = joinpath(@__DIR__, "results/resource/power/CO_power.csv")

    C2H4_electrolyzer_filepath = joinpath(@__DIR__, "results/econ/electrolyzer/C2H4_electrolyzer.csv")
    HCOO_electrolyzer_filepath = joinpath(@__DIR__, "results/econ/electrolyzer/HCOO_electrolyzer.csv")
    CO_electrolyzer_filepath = joinpath(@__DIR__, "results/econ/electrolyzer/CO_electrolyzer.csv")

    C2H4_separations_filepath = joinpath(@__DIR__, "results/econ/separations/C2H4_separations.csv")
    HCOO_separations_filepath = joinpath(@__DIR__, "results/econ/separations/HCOO_separations.csv")
    CO_separations_filepath = joinpath(@__DIR__, "results/econ/separations/CO_separations.csv")

    C2H4_battery_filepath = joinpath(@__DIR__, "results/econ/battery/C2H4_battery.csv")
    HCOO_battery_filepath = joinpath(@__DIR__, "results/econ/battery/HCOO_battery.csv")
    CO_battery_filepath = joinpath(@__DIR__, "results/econ/battery/CO_battery.csv")

    C2H4_maintenance_filepath = joinpath(@__DIR__, "results/econ/maintenance/C2H4_maintenance.csv")
    HCOO_maintenance_filepath = joinpath(@__DIR__, "results/econ/maintenance/HCOO_maintenance.csv")
    CO_maintenance_filepath = joinpath(@__DIR__, "results/econ/maintenance/CO_maintenance.csv")

    C2H4_water_filepath = joinpath(@__DIR__, "results/econ/water/C2H4_water.csv")
    HCOO_water_filepath = joinpath(@__DIR__, "results/econ/water/HCOO_water.csv")
    CO_water_filepath = joinpath(@__DIR__, "results/econ/water/CO_water.csv")

    # Initialize flow output array
    flow_results = zeros(Float64, num_source_pts, num_sink_pts, num_species)

    # Populate each output array with JuMP varaible value found in optimizer
    for i = 1:num_source_pts
        for j = 1:num_sink_pts
            for k in 1:num_species
                flow_results[i, j, k] = JuMP.value(flow[i, j, k])
            end
        end
    end

    # Initialize NPV output arrays for each molecule
    results_NPV = zeros(Float64, num_source_pts, num_sink_pts, num_species)

    # Initialize electricity cost output arrays for each molecule
    results_elec = zeros(Float64, num_source_pts, num_sink_pts, num_species)

    # Initialize electricity cost output arrays for each molecule
    results_trans = zeros(Float64, num_source_pts, num_sink_pts, num_species)

    # Initialize electricity cost output arrays for each molecule
    results_ccu = zeros(Float64, num_source_pts, num_sink_pts, num_species)

    # Initialize electricity cost output arrays for each molecule
    results_tea = zeros(Float64, num_source_pts, num_sink_pts, num_species)

    # Initialize electricity cost output arrays for each molecule
    results_prf = zeros(Float64, num_source_pts, num_sink_pts, num_species)

    #Initialize water consumption output arrays for each molecule
    results_water = zeros(Float64, num_source_pts, num_sink_pts, num_species)

    #Initialize power consumption output arrays for each molecule
    results_power = zeros(Float64, num_source_pts, num_sink_pts, num_species)

    #Initialize electrolyzer cost output arrays for each molecule
    results_electrolyzer = zeros(Float64, num_source_pts, num_sink_pts, num_species)

    #Initialize separations cost output arrays for each molecule
    results_separations = zeros(Float64, num_source_pts, num_sink_pts, num_species)

    #Initialize battery cost output arrays for each molecule
    results_battery = zeros(Float64, num_source_pts, num_sink_pts, num_species)

    #Initialize electrolyzer cost output arrays for each molecule
    results_maintenance = zeros(Float64, num_source_pts, num_sink_pts, num_species)

    #Initialize water cost output arrays for each molecule
    results_water_rsrc = zeros(Float64, num_source_pts, num_sink_pts, num_species)

    # Populate each output array with JuMP varaible value found in optimizer
    for i = 1:num_source_pts
        for j = 1:num_sink_pts
            for k = 1:num_species
                results_NPV[i, j, k] = sum(
                (
                - (D[i, j, k] * transport_A + transport_B) * flow_results[i, j, k] # Transportation cost
                -  electrolyzer_TEA(
                    flow_results[i, j, k],
                    product_phase[k],
                    product_mw[k],
                    product_density[k],
                    electrons_product[k],
                    electrons_co2[k],
                    system_potential[k],
                    selectivity[k],
                    hours_year,
                    current_density[k],
                    single_pass_conversion[k],
                    electrolyzer_cost,
                    electricity_price[i],
                    water_price,
                    water_rate[k],
                    crf)  # TEA system cost
                - flow_results[i, j, k] * carbons[k] * C_mw / product_mw[k] * CO2_cap_cost[i] # Carbon capture cost
                - flow_results[i, j, k] / hours_year / 3600 * 1e3 / product_mw[k] * electrons_product[k] * F * system_potential[k] / selectivity[k] / 1e3 * hours_year * electricity_price[i] # Electrochemical electricity price
                + flow_results[i, j, k] * value[k] # Product value
                )
                )

                if flow_results[i, j, k] != 0
                    results_electrolyzer[i, j, k],
                    results_separations[i, j, k],
                    results_battery[i, j, k],
                    results_maintenance[i, j, k],
                    results_water[i, j, k] =
                    electrolyzer_TEA_array(
                        flow_results[i, j, k],
                        product_phase[k],
                        product_mw[k],
                        product_density[k],
                        electrons_product[k],
                        electrons_co2[k],
                        system_potential[k],
                        selectivity[k],
                        hours_year,
                        current_density[k],
                        single_pass_conversion[k],
                        electrolyzer_cost,
                        electricity_price[i],
                        water_price,
                        water_rate[k],
                        crf)  # TEA system cost
                    results_trans[i, j, k] = (D[i, j, k] * transport_A + transport_B) * flow_results[i, j, k] # Transportation cost
                    results_tea[i, j, k] =  electrolyzer_TEA(
                        flow_results[i, j, k],
                        product_phase[k],
                        product_mw[k],
                        product_density[k],
                        electrons_product[k],
                        electrons_co2[k],
                        system_potential[k],
                        selectivity[k],
                        hours_year,
                        current_density[k],
                        single_pass_conversion[k],
                        electrolyzer_cost,
                        electricity_price[i],
                        water_price,
                        water_rate[k],
                        crf)  # TEA system cost
                    results_ccu[i, j, k] =  flow_results[i, j, k] * carbons[k] * C_mw / product_mw[k] * CO2_cap_cost[i] # Carbon capture cost
                    results_elec[i, j, k] =  flow_results[i, j, k] / hours_year / 3600 * 1e3 / product_mw[k] * electrons_product[k] * F * system_potential[k] / selectivity[k] / 1e3 * hours_year * electricity_price[i] # Electrochemical electricity price
                    results_prf[i, j, k] = flow_results[i, j, k] * value[k] # Product value
                    results_water_rsrc[i, j, k] = flow_results[i, j, k] / hours_year / 3600 * 1e3 / product_mw[k] * electrons_product[k] * F / selectivity[k] / H2O_electrons / F * H2O_mw / 1e3 * hours_year * 3600
                    results_power[i, j, k] = flow_results[i, j, k] / hours_year / 3600 * 1e3 / product_mw[k] * electrons_product[k] * F * system_potential[k] / selectivity[k] / 1e3 #kW
                end
            end
        end
    end

    # Write sink/source data to CSV files for data visualization
    # Useful for runs where points are excluded from optimization routine
    source_index_path = joinpath(@__DIR__, "results/index/source_index.csv")
    sink_index_path = joinpath(@__DIR__, "results/index/sink_index.csv")

    CSV.write(source_index_path, DataFrame(source_pts))
    CSV.write(sink_index_path, DataFrame(sink_pts))

    # Write output flow files to CSV
    CSV.write(C2H4_result_filepath, DataFrame(flow_results[:, :, 1]))
    CSV.write(HCOO_result_filepath, DataFrame(flow_results[:, :, 2]))
    CSV.write(CO_result_filepath, DataFrame(flow_results[:, :, 3]))

    # write output NPV dataframes to CSV
    CSV.write(C2H4_NPV_filepath, DataFrame(results_NPV[:, :, 1]))
    CSV.write(HCOO_NPV_filepath, DataFrame(results_NPV[:, :, 2]))
    CSV.write(CO_NPV_filepath, DataFrame(results_NPV[:, :, 3]))

    # write output electricity cost dataframes to CSV
    CSV.write(C2H4_trans_filepath, DataFrame(results_trans[:, :, 1]))
    CSV.write(HCOO_trans_filepath, DataFrame(results_trans[:, :, 2]))
    CSV.write(CO_trans_filepath, DataFrame(results_trans[:, :, 3]))

    # write output electricity cost dataframes to CSV
    CSV.write(C2H4_tea_filepath, DataFrame(results_tea[:, :, 1]))
    CSV.write(HCOO_tea_filepath, DataFrame(results_tea[:, :, 2]))
    CSV.write(CO_tea_filepath, DataFrame(results_tea[:, :, 3]))

    # write output electricity cost dataframes to CSV
    CSV.write(C2H4_ccu_filepath, DataFrame(results_ccu[:, :, 1]))
    CSV.write(HCOO_ccu_filepath, DataFrame(results_ccu[:, :, 2]))
    CSV.write(CO_ccu_filepath, DataFrame(results_ccu[:, :, 3]))

    # write output electricity cost dataframes to CSV
    CSV.write(C2H4_elec_filepath, DataFrame(results_elec[:, :, 1]))
    CSV.write(HCOO_elec_filepath, DataFrame(results_elec[:, :, 2]))
    CSV.write(CO_elec_filepath, DataFrame(results_elec[:, :, 3]))

    # write output electricity cost dataframes to CSV
    CSV.write(C2H4_prf_filepath, DataFrame(results_prf[:, :, 1]))
    CSV.write(HCOO_prf_filepath, DataFrame(results_prf[:, :, 2]))
    CSV.write(CO_prf_filepath, DataFrame(results_prf[:, :, 3]))

    # write output water consumption dataframes to CSV
    CSV.write(C2H4_water_rsrc_filepath, DataFrame(results_water_rsrc[:, :, 1]))
    CSV.write(HCOO_water_rsrc_filepath, DataFrame(results_water_rsrc[:, :, 2]))
    CSV.write(CO_water_rsrc_filepath, DataFrame(results_water_rsrc[:, :, 3]))

    # write output power consumption dataframes to CSV
    CSV.write(C2H4_power_filepath, DataFrame(results_power[:, :, 1]))
    CSV.write(HCOO_power_filepath, DataFrame(results_power[:, :, 2]))
    CSV.write(CO_power_filepath, DataFrame(results_power[:, :, 3]))

    # write output electrolyzer cost dataframes to CSV
    CSV.write(C2H4_electrolyzer_filepath, DataFrame(results_electrolyzer[:, :, 1]))
    CSV.write(HCOO_electrolyzer_filepath, DataFrame(results_electrolyzer[:, :, 2]))
    CSV.write(CO_electrolyzer_filepath, DataFrame(results_electrolyzer[:, :, 3]))

    # write output separations cost dataframes to CSV
    CSV.write(C2H4_separations_filepath, DataFrame(results_separations[:, :, 1]))
    CSV.write(HCOO_separations_filepath, DataFrame(results_separations[:, :, 2]))
    CSV.write(CO_separations_filepath, DataFrame(results_separations[:, :, 3]))

    # write output battery cost dataframes to CSV
    CSV.write(C2H4_battery_filepath, DataFrame(results_battery[:, :, 1]))
    CSV.write(HCOO_battery_filepath, DataFrame(results_battery[:, :, 2]))
    CSV.write(CO_battery_filepath, DataFrame(results_battery[:, :, 3]))

    # write output maintenance cost dataframes to CSV
    CSV.write(C2H4_maintenance_filepath, DataFrame(results_maintenance[:, :, 1]))
    CSV.write(HCOO_maintenance_filepath, DataFrame(results_maintenance[:, :, 2]))
    CSV.write(CO_maintenance_filepath, DataFrame(results_maintenance[:, :, 3]))

    # write output water cost dataframes to CSV
    CSV.write(C2H4_water_filepath, DataFrame(results_water[:, :, 1]))
    CSV.write(HCOO_water_filepath, DataFrame(results_water[:, :, 2]))
    CSV.write(CO_water_filepath, DataFrame(results_water[:, :, 3]))


    println("Optimization Routine Completed!")



end

@time co2_optimizer()
