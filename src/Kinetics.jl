# -- INTERNAL METHODS ---------------------------------------------------------- #
function compute_transcription_kinetic_array(biophysical_dictionary::Dict{Symbol,Any}, model_dictionary::Dict{String,Any})

    # get stuff from the biophysical dictionary -
    RNAPII_copy_number = biophysical_dictionary[:copies_of_rnapII_per_cell]
    characteristic_length = biophysical_dictionary[:characteristic_transcript_length]
    transcription_elongation_rate = biophysical_dictionary[:transcription_elongation_rate]
    characteristic_initiation_time = biophysical_dictionary[:characteristic_initiation_time_transcription]
    KX = biophysical_dictionary[:transcription_saturation_constant]
    mass_of_single_cell = biophysical_dictionary[:mass_of_single_cell]
    fraction_of_water_per_cell = biophysical_dictionary[:fraction_of_water_per_cell]

    # fraction of dry weight -
    fraction_dry_cell = 1 - fraction_of_water_per_cell      # dimensionless

    # avagodros number -
    av_number = 6.02e23                                     # number/mol

    # what is the RNAP concentration -
    RNAPII_concentration = RNAPII_copy_number*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell) # nmol/gDW

    # initialize -
    transcription_kinetics_array = Float64[]

    # compute the kinetic limit -
    list_of_gene_nodes = model_dictionary["list_of_nodes"]
    number_of_nodes = length(list_of_gene_nodes)
    for node_index = 1:number_of_nodes

        # get the dictionary for this node -
        node_key = "N$(node_index)"
        node_dictionary = list_of_gene_nodes[node_key]

        # compute elongation constant -
        gene_length = parse(Float64,node_dictionary["parameters"]["gene_length_in_nt"])
        gene_copy_number = parse(Float64, node_dictionary["parameters"]["gene_copy_number_per_cell"])

        # compute kE -
        kE = transcription_elongation_rate*(1/gene_length)

        # compute kI -
        kI = (1/characteristic_initiation_time)

        # compute the tau factor -
        tau_factor = (kE/kI)

        # compute the gene concentration -
        G = gene_copy_number*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell)   # nmol/gDW

        # Compute the rate -
        value = kE*RNAPII_concentration*(G/(KX*tau_factor+(1+tau_factor)*G))*(3600)

        # push -
        push!(transcription_kinetics_array, value)
    end

    # return -
    return transcription_kinetics_array
end

function compute_translation_parameter_array(biophysical_dictionary::Dict{Symbol,Any}, model_dictionary::Dict{String,Any})

    # get some stuff from the biophysical_dictionary -
    mass_of_single_cell = biophysical_dictionary[:mass_of_single_cell]
    fraction_of_water_per_cell = biophysical_dictionary[:fraction_of_water_per_cell]
    ribosome_copy_number = biophysical_dictionary[:copies_of_ribosome_per_cell]
    translation_elongation_rate = biophysical_dictionary[:translation_elongation_rate]
    characteristic_initiation_time_translation = biophysical_dictionary[:characteristic_initiation_time_translation]
    KL = biophysical_dictionary[:translation_saturation_constant]

    # fraction of dry weight -
    fraction_dry_cell = 1 - fraction_of_water_per_cell      # dimensionless

    # avagodros number -
    av_number = 6.02e23                                     # number/mol

    # what is the ribosome concentration -
    ribosome_concentration = ribosome_copy_number*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell) # nmol/gDW

    # initialize -
    translation_parameters_array = Array{TranslationParameters,1}()

    # compute the translation parameters -
    list_of_gene_nodes = model_dictionary["list_of_nodes"]
    number_of_nodes = length(list_of_gene_nodes)
    for node_index = 1:number_of_nodes

        # get the dictionary for this node -
        node_key = "N$(node_index)"
        node_dictionary = list_of_gene_nodes[node_key]

        # compute elongation constant -
        protein_length = parse(Float64,node_dictionary["parameters"]["protein_length_in_aa"])

        # compute kE -
        kE = translation_elongation_rate*(1/protein_length)

        # compute kI -
        kI = (1/characteristic_initiation_time_translation)

        # compute the tau factor -
        tau_factor = (kE/kI)

        # compute the vmax -
        translation_vmax = kE*ribosome_concentration

        # build -
        translationParameters = TranslationParameters()
        translationParameters.vmax = translation_vmax*(3600)    # convert to hr
        translationParameters.tau_factor = tau_factor
        translationParameters.KL = KL
        push!(translation_parameters_array, translationParameters)
    end

    return translation_parameters_array
end
# -- INTERNAL METHODS ---------------------------------------------------------- #
