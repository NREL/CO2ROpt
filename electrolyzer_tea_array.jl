


function electrolyzer_TEA_array(
    flow,
    product_phase,
    product_mw,
    product_density,
    electrons_product,
    electrons_co2,
    system_potential,
    selectivity,
    hours_year,
    current_density,
    single_pass_conversion,
    electrolyzer_cost,
    electricity_price,
    water_price,
    water_rate,
    crf)


    #_____________________________________________________________________________

    F = 96485 # [C/mol]
    g_kg = 1e3
    sec_hour = 3600
    W_MW = 1e6
    CO2_mw = 44 # [g/mol]
    CO2_density = 1.98 #[kg/m^3]
    l_m3 = 1e3 # [liters / m^3]
    sec_min = 60
    cm2_m2 = 100^2
    min_hour = 60
    h2_mw = 2.015 #[g/mol]
    h2_electrons = 2
    h2_density = 0.08375 #[kg/m^3]
    days_year = hours_year / 24
    hours_day = 24
    watergal_kg = 0.2642
    H2O_electrons = 4 # [mol e- / mol product]
    H2O_mw = 18 # [g/mol]

    # [total current]
    amps = flow * F / hours_year / sec_hour / (product_mw / g_kg) * electrons_product / selectivity

    #[m^2 active area]
    electrolyzer_area = amps / current_density / cm2_m2

    #[kg/day CO2]
    co2_flow_needed = amps / electrons_co2 / F * (CO2_mw / g_kg) * sec_hour * hours_day

    # [kg/hour CO2] inlet flow to electrolyzer
    co2_inlet_flow = co2_flow_needed / hours_day / single_pass_conversion

    # [m^3/hour CO2] outlet flow from electrolyzer
    co2_outlet_flow = (co2_inlet_flow - (co2_flow_needed / hours_day))/ CO2_density#co2_inlet_flow * (1 - single_pass_conversion) / CO2_density

    if product_phase == "Liquid"

        #[m^3/hour]
        gas_product_flow = 0.0
        #[l/min]
        liquid_product_flow = flow / product_density * l_m3 / hours_year / min_hour

    elseif product_phase == "Gas"

        #[m^3/hour]
        gas_product_flow = flow / product_density / hours_year
        #[l/min]
        liquid_product_flow = 0.0


    else
        println("Error in product phase classification!")

    end

    #[gal/day]
    water_flow = amps / H2O_electrons / F * (H2O_mw / g_kg) * watergal_kg * sec_hour * hours_day

    #[m^3/hour]
    h2_outlet_flow = amps * (1 - selectivity) / h2_electrons / F * h2_mw / g_kg / h2_density * sec_hour

    #[m^3/hour]
    gas_outlet_flow = co2_outlet_flow + gas_product_flow + h2_outlet_flow

    #[l/min]
    liquid_outlet_flow = liquid_product_flow

    # [$]
    electrolyzer_capex = electrolyzer_area * electrolyzer_cost

    # [$]
    electrolyzer_BOP = electrolyzer_capex * 0.35

    # [$] https://www.eia.gov/analysis/studies/electricity/batterystorage/pdf/battery_storage.pdf
    battery_capex = 1554 * (amps * system_potential / 1e3)

    #[$]
    PSA_capex = 1323.5 * gas_outlet_flow + 632838

    #[$]
    distil_capex = 2773.8 * liquid_outlet_flow + 1e6

    # [$]
    capex =
    electrolyzer_capex +
    electrolyzer_BOP +
    PSA_capex +
    distil_capex +
    battery_capex

    electrolyzer_capex *= crf
    electrolyzer_BOP *= crf
    PSA_capex *= crf
    distil_capex *= crf
    battery_capex *= crf

    # [$/year]
    maintenance_cost = capex * 0.1

    # [$/year]
    PSA_opex = 0.25 * gas_outlet_flow * hours_year * electricity_price

    # [$/year]
    distil_opex = 9895 * (gas_outlet_flow / 1000) * days_year

    # [$/year]
    water_opex = water_rate * flow * water_price

    # [$/year]
    opex =
    maintenance_cost +
    PSA_opex +
    distil_opex +
    water_opex

    # [$/year]
    annualized_capex = capex * crf

    # [$/year]
    annual_costs = annualized_capex + opex

    return electrolyzer_capex + electrolyzer_BOP, PSA_capex + PSA_opex + distil_capex + distil_opex, battery_capex, maintenance_cost, water_opex

end
