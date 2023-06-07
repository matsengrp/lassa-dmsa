SEGMENTS = ["s"]

rule all:
    input:
        #auspice_tree = expand("auspice/lassa_{segment}_tree.json", segment=SEGMENTS),
        #auspice_meta = expand("auspice/lassa_{segment}_meta.json", segment=SEGMENTS)
        auspice = expand("auspice/lassa-{segment}.json", segment=SEGMENTS),
        auspice_root_sequence = expand("auspice/lassa-{segment}_root-sequence.json", segment=SEGMENTS)

rule files:
    params:
        dropped_strains = "config/dropped_strains.txt",
        reference = "config/lassa_{segment}.gb",
        colors = "config/colors.tsv",
        auspice_config = "config/auspice_config.json"

files = rules.files.params

#rule parse:
#    message: "Parsing fasta into sequences and metadata"
#    #input:
#    #    sequences = files.input_fasta
#    output:
#        # sequences = "results/sequences_{segment}.fasta",
#        sequences = lambda w: config["inputs"][f"{w.segment}"]["sequences"],
#        metadata = lambda w: config["inputs"][f"{w.segment}"]["metadata"]
#    #params:
#    #    fasta_fields = "strain accession segment date region country host authors title journal puburl"
#    #shell:
#    #    """
#    #    augur parse \
#    #        --sequences {input.sequences} \
#    #        --output-sequences {output.sequences} \
#    #        --output-metadata {output.metadata} \
#    #        --fields {params.fasta_fields}
#    #    """

rule filter:
    message:
        """
        Filtering to
          - excluding strains in {input.exclude}
        """
    input:
        sequences = lambda w: config["inputs"][f"{w.segment}"]["sequences"],
        metadata = lambda w: config["inputs"][f"{w.segment}"]["metadata"],
        #sequences = rules.parse.output.sequences,
        #metadata = rules.parse.output.metadata,
        exclude = files.dropped_strains
    output:
        sequences = "results/filtered_{segment}.fasta"
    #params:
    #    group_by = "country year",
    #    sequences_per_group = 2,
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output {output.sequences}
        """
# - {params.sequences_per_group} sequence(s) per {params.group_by!s}
#--group-by {params.group_by} \
#--sequences-per-group {params.sequences_per_group}

rule align:
    message:
        """
        Aligning sequences to {params.reference_name}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        # reference_name = lambda w: config["inputs"][f"{w.segment}"]["reference_name"]
    output:
        alignment = "results/aligned_{segment}.fasta"
    params:
        reference_name = lambda w: config["inputs"][f"{w.segment}"]["reference_name"]
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-name {params.reference_name} \
            --output {output.alignment} \
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
        # metadata = rules.parse.output.metadata
        # metadata = rules.filter.output.metadata
        metadata = lambda w: config["inputs"][f"{w.segment}"]["metadata"]
    output:
        tree = "results/tree_{segment}.nwk",
        node_data = "results/branch_lengths_{segment}.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_rate = 0.0006
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
        # reference_annotation = files.reference
        reference_annotation = lambda w: config["inputs"][f"{w.segment}"]["reference_annotation"]
    output:
        node_data = "results/aa_muts_{segment}.json",
        alignments = expand("results/translations/{{segment}}_{gene}.fasta", gene=["NP", "GPC"])
        
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --genes NP GPC \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference_annotation} \
            --output-node-data {output.node_data} \
            --alignment-output results/translations/{wildcards.segment}_%GENE.fasta
        """

rule polyclonal_escape_prediction:
    input:
        alignment = "results/translations/s_GPC.fasta"
    output:
        node_data = "results/{serum}_polyclonal_escape_prediction.json",
        pred_data = "results/{serum}_polyclonal_escape_prediction.csv"
    log:
        "logs/{serum}_polyclonal_escape_prediction.txt"
    params:
        dms_wt_seq_id = lambda w: config["polyclonal_serum_models"][f"{w.serum}"]["dms_wt_seq_id"],
        mut_effects_df = lambda w: config["polyclonal_serum_models"][f"{w.serum}"]["mut_effects_df"],
        mut_effect_col = lambda w: config["polyclonal_serum_models"][f"{w.serum}"]["mut_effect_col"],
        mutation_col = lambda w: config["polyclonal_serum_models"][f"{w.serum}"]["mutation_col"],
        activity_wt_df = lambda w: config["polyclonal_serum_models"][f"{w.serum}"]["activity_wt_df"],
        concentrations = lambda w: config["polyclonal_serum_models"][f"{w.serum}"]["concentrations"] if "concentrations" in config["polyclonal_serum_models"][f"{w.serum}"] else "0.0",
        icxx = lambda w: config["polyclonal_serum_models"][f"{w.serum}"]["icxx"] if "icxx" in config["polyclonal_serum_models"][f"{w.serum}"] else 0.0
    conda:
        "my_profiles/dmsa-pred/dmsa_env.yaml"
    resources:
        mem_mb=2000
    shell:
        """
        python my_profiles/dmsa-pred/dmsa_pred.py polyclonal-escape \
            --alignment {input.alignment} \
            --dms-wt-seq-id {params.dms_wt_seq_id} \
            --mut-effects-df {params.mut_effects_df} \
            --mut-effect-col {params.mut_effect_col} \
            --mutation-col {params.mutation_col} \
            --activity-wt-df {params.activity_wt_df} \
            --concentrations {params.concentrations} \
            --icxx {params.icxx} \
            --experiment-label {wildcards.serum} \
            --output-json {output.node_data} \
            --output-df {output.pred_data} 2>&1 | tee {log}
        """

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        # metadata = rules.parse.output.metadata
        metadata = lambda w: config["inputs"][f"{w.segment}"]["metadata"]
    output:
        node_data = "results/traits_{segment}.json",
    params:
        columns = "country"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

def _get_polyclonal_node_data(wildcards):
    inputs=[]
    if "polyclonal_serum_models" in config:
        inputs += list(expand(
            rules.polyclonal_escape_prediction.output.node_data, 
            serum=list(config["polyclonal_serum_models"])
        ))
    return inputs


rule auspice_config:
    message: "Getting auspice config for modification"
    input:
        default_auspice_config = files.auspice_config
    output:
        "results/dmsa_modified_auspice_config.json"
    params:
        path_config = workflow.overwrite_configfiles[0]
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
        # metadata = rules.parse.output.metadata,
        metadata = lambda w: config["inputs"][f"{w.segment}"]["metadata"],
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        escape_predictions = _get_polyclonal_node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        colors = files.colors,
        # auspice_config = files.auspice_config
        auspice_config = rules.auspice_config.output
    output:
        auspice_json = "auspice/lassa-{segment}.json",
        root_sequence_json = "auspice/lassa-{segment}_root-sequence.json" 
        #auspice_tree = "auspice/lassa_{segment}_tree.json",
        #auspice_meta = "auspice/lassa_{segment}_meta.json"
    log:
        "logs/export_{segment}.txt"
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
#--output-tree {output.auspice_tree} \
#--output-meta {output.auspice_meta}
