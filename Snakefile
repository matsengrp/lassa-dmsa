
# Initialize parameters
SEGMENTS = ["S"]
GENOME_SIZE_THRESHOLDS = [1000, 4000] # LASV S segment ~ 3kb and LASV L segment ~ 7kb
MAX_FRAC_N = 0 # Maximum fraction of ambiuous bases allowed
ACCESSIONS_TO_EXCLUDE = [
    "MK107912", # this sequences is an outlier with many more bases that do not align
    "MK281631", # this sequences has many 'X' amino acids
    "MK118006", # this sequences has many 'X' and more bases that do not align
    "MT193286", # this sequences has many 'X' and more bases that do not align
    "MT193287", # this sequences has many 'X' and more bases that do not align
    "K117948 ",# this sequences has many 'X' and more bases that do not align
    "MK117981", # this sequences has many 'X' and more bases that do not align
    "MK118036", # this sequences has many 'X' and more bases that do not align
    "MK118037", # this sequences has many 'X' and more bases that do not align
    "MH887756", # this sequences has many 'X' and more bases that do not align
    "MH157039", # this sequences has many 'X' and more bases that do not align
    "MH053490", # this sequences has many 'X' and more bases that do not align
    "MK117951", # this sequences has many 'X' and more bases that do not align
    "MK118005", # this sequences has many 'X' and more bases that do not align
    "MK118034", # this sequences has many 'X' and more bases that do not align
    "MK117850", # this sequences has many 'X' and more bases that do not align
    "MK117954", # this sequences has many 'X' and more bases that do not align
    "MK117977", # this sequences has many 'X' and more bases that do not align
    "MK118013", # more of a borderline one b/c of 'X'
    "MK118015", # this sequences has many 'X' and more bases that do not align
    "MK118017", # this sequences has many 'X' and more bases that do not align
    "MK118032", # this sequences has many 'X' and more bases that do not align
    "MH053523", # this sequences has many 'X' and more bases that do not align
    "MK118027", # more of a borderline one b/c has N bases instead of stop codon and adds trailing amino acids
]
DMS_WT_SEQ_ID = "Josiah_NC_004296_2018-08-13"


rule all:
    input:
        auspice = expand("auspice/lassa-{segment}.json", segment=SEGMENTS),
        auspice_root_sequence = expand("auspice/lassa-{segment}_root-sequence.json", segment=SEGMENTS)

rule files:
    params:
        accessions = "config/all_LASV_accessions_081023.txt", # List of accessions to download
        reference = "config/LASV_S_reference_NC_004296_RC.gb", # Josiah reference
        dropped_strains = "config/dropped_strains.txt",
        colors = "config/colors.tsv",
        auspice_config = "config/auspice_config.json",
        allow_missing_sites = "allowed_missing_sites.txt"

files = rules.files.params

rule download_data:
    message:
        """
        Processing the following accessions:
          - {input.accessions}
        """
    input:
        accessions = files.accessions
    params:
        genome_size_threshold_lower = GENOME_SIZE_THRESHOLDS[0],
        genome_size_threshold_upper = GENOME_SIZE_THRESHOLDS[1],
        desired_segment = SEGMENTS[0],
        max_frac_N = MAX_FRAC_N,
        accesstions_to_exclude = ACCESSIONS_TO_EXCLUDE,
    output:
        sequences = "results/unfiltered_{segment}.fasta",
        metadata = "results/{segment}_metadata.tsv",
    conda: 
        "my_profiles/dmsa-pred/dmsa_env.yaml"
    log:
        "logs/download_{segment}_data.txt"
    script:
        "scripts/download_NCBI_sequences.py"

rule filter:
    message:
        """
        Filtering to
          - excluding strains in {input.exclude}
        """
    input:
        sequences = "results/unfiltered_{segment}.fasta",
        metadata = "results/{segment}_metadata.tsv",
        exclude = files.dropped_strains
    output:
        sequences = "results/filtered_{segment}.fasta"
    conda: 
        "my_profiles/dmsa-pred/dmsa_env.yaml"
    log:
        "logs/filter_{segment}.txt"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output {output.sequences}
        """

rule align:
    message:
        """
        Aligning sequences to {params.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
    output:
        alignment = "results/aligned_{segment}.fasta"
    params:
        reference = files.reference
    conda: 
        "my_profiles/dmsa-pred/dmsa_env.yaml"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {params.reference} \
            --output {output.alignment} \
            --remove-reference \
            --fill-gaps
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree_raw_{segment}.nwk"
    params:
        method = "iqtree"
    conda: 
        "my_profiles/dmsa-pred/dmsa_env.yaml"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method {params.method}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - fix clock rate at {params.clock_rate}
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = "results/{segment}_metadata.tsv",
    output:
        tree = "results/tree_{segment}.nwk",
        node_data = "results/branch_lengths_{segment}.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_rate = 0.0006
    conda: 
        "my_profiles/dmsa-pred/dmsa_env.yaml"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --clock-rate {params.clock_rate} \
            --date-confidence \
            --date-inference {params.date_inference}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt_muts_{segment}.json"
    params:
        inference = "joint"
    conda: 
        "my_profiles/dmsa-pred/dmsa_env.yaml"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/aa_muts_{segment}.json",
        alignments = expand("results/translations/{{segment}}_{gene}.fasta", gene=["Glycoprotein", "Nucleoprotein"])
    conda: 
        "my_profiles/dmsa-pred/dmsa_env.yaml"
    log:
        "logs/translate_{segment}.txt"
        
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --genes Glycoprotein Nucleoprotein \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
            --alignment-output results/translations/{wildcards.segment}_%GENE.fasta
        """

# Note, this is hardcoded to only make predictions for the GPC gene in the S segment build
# additional 'segment' and 'gene' wildcards would need to be defined to make this more general
rule variant_escape_prediction:
    input:
        alignment = "results/translations/S_Glycoprotein.fasta"
    output:
        node_data = "results/dmsa-phenotype/{collection}/{experiment}_escape_prediction.json",
        pred_data = "results/dmsa-phenotype/{collection}/{experiment}_escape_prediction.csv"
    log:
        "logs/{collection}/{experiment}_escape_prediction.txt"
    params:
        basedir = lambda w: config["dmsa_phenotype_collections"].get(w.collection)['mut_effects_dir'],
        dms_wt_seq_id = DMS_WT_SEQ_ID,
        mut_effect_col = lambda w: config["dmsa_phenotype_collections"].get(w.collection)['mut_effect_col'],
        mutation_col = lambda w: config["dmsa_phenotype_collections"].get(w.collection)['mutation_col'],
        mut_effects_df = lambda w: os.path.join(
            config["dmsa_phenotype_collections"].get(w.collection)['mut_effects_dir'], 
            w.experiment
        ),
    conda:
        "my_profiles/dmsa-pred/dmsa_env.yaml"
    shell:
        """
        python my_profiles/dmsa-pred/dmsa_pred.py phenotype-prediction \
            --model-type additive \
            --alignment {input.alignment} \
            --dms-wt-seq-id {params.dms_wt_seq_id} \
            --mask-seqs-with-disallowed-aa-subs False \
            --min-pred-pheno 0.0 \
            --mut-effects-df {params.mut_effects_df} \
            --mut-effect-col {params.mut_effect_col} \
            --mutation-col {params.mutation_col} \
            --experiment-label {wildcards.experiment} \
            --output-json {output.node_data} \
            --output-df {output.pred_data} 2>&1 | tee {log}
        """

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = "results/{segment}_metadata.tsv",
    output:
        node_data = "results/traits_{segment}.json",
    params:
        columns = "country"
    conda: 
        "my_profiles/dmsa-pred/dmsa_env.yaml"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

def _get_variant_escape_node_data(wildcards):
    inputs=[]
    wildcards_dict = dict(wildcards)

    import glob
    for collection_name, collection_dict in config['dmsa_phenotype_collections'].items():

        # run the predictions using every csv in the glob
        requested_files = expand(
            rules.variant_escape_prediction.output.node_data,
            collection=collection_name,
            experiment=[
                os.path.basename(fp) 
                for fp in glob.glob(collection_dict['mut_effects_dir']+"/*.csv")
            ],
            **wildcards_dict
        )
        inputs.extend(requested_files)

    return inputs


rule auspice_config:
    message: "Getting auspice config for modification"
    input:
        default_auspice_config = files.auspice_config
    output:
        "results/dmsa_modified_auspice_config.json"
    params:
        path_config = workflow.overwrite_configfiles[0]
    conda:
        "my_profiles/dmsa-pred/dmsa_env.yaml"
    shell:
        """
        python my_profiles/dmsa-pred/modify_auspice_config.py \
            --auspice-config-path {input.default_auspice_config} \
            --snake-config-path {params.path_config} \
            --output-config-path {output}
        """

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = "results/{segment}_metadata.tsv",
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        escape_predictions = _get_variant_escape_node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        colors = files.colors,
        auspice_config = files.auspice_config
    output:
        auspice_json = "auspice/lassa-{segment}.json",
        root_sequence_json = "auspice/lassa-{segment}_root-sequence.json" 
    log:
        "logs/export_{segment}.txt"
    conda: 
        "my_profiles/dmsa-pred/dmsa_env.yaml"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.escape_predictions} \
            --include-root-sequence \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json} 2>&1 | tee {log}
        """
