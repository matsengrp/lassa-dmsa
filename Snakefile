SEGMENTS = ["s"]

rule all:
    input:
        auspice_tree = expand("auspice/lassa_{segment}_tree.json", segment=SEGMENTS),
        auspice_meta = expand("auspice/lassa_{segment}_meta.json", segment=SEGMENTS)

rule files:
    params:
        input_fasta = "data/lassa_{segment}.fasta",
        dropped_strains = "config/dropped_strains.txt",
        reference = "config/lassa_{segment}.gb",
        colors = "config/colors.tsv",
        auspice_config = "config/auspice_config.json"

files = rules.files.params

rule parse:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = files.input_fasta
    output:
        sequences = "results/sequences_{segment}.fasta",
        metadata = "results/metadata_{segment}.tsv"
    params:
        fasta_fields = "strain accession segment date region country host authors title journal puburl"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields}
        """

rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - excluding strains in {input.exclude}
        """
    input:
        sequences = rules.parse.output.sequences,
        metadata = rules.parse.output.metadata,
        exclude = files.dropped_strains
    output:
        sequences = "results/filtered_{segment}.fasta"
    params:
        group_by = "country year",
        sequences_per_group = 2,
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned_{segment}.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
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
        metadata = rules.parse.output.metadata
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
        reference = files.reference
    output:
        node_data = "results/aa_muts_{segment}.json",
        alignments = expand("results/translations/{{segment}}_{gene}.fasta", gene=["NP", "GPC"])
        
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --genes NP GPC \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
            --alignment-output results/translations/{wildcards.segment}_%GENE.fasta
        """

rule polyclonal_escape_prediction:
    input:
        alignments = "results/translations/s_GPC.fasta"
    output:
        node_data = "results/polclonal_escape_prediction.json"
    log:
        "logs/polclonal_escape_prediction.txt"
    params:
        dms_wt_seq = "Josiah",
        mut_escape_df = "my_profiles/polyclonal-data/89F/89F_avg.csv",
        activity_wt_df = "my_profiles/polyclonal-data/89F/89F_epitope.csv",
        concentrations = "5.0,20.0,80.0",
    conda:
        "my_profiles/dmsa-pred/dmsa_env.yaml"
    resources:
        mem_mb=2000
    shell:
        """
        python my_profiles/dmsa-pred/dmsa_pred.py polyclonal-escape \
            --activity-wt-df {params.activity_wt_df} \
            --concentrations {params.concentrations} \
            --escape-column escape_mean \
            --alignment {input.alignments} \
            --mut-effects-df {params.mut_escape_df} \
            --dms-wt-seq-id {params.dms_wt_seq:q} \
            --experiment-label 89F \
            --output {output.node_data} 2>&1 | tee {log}
        """

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata
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

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        escape_predictions = rules.polyclonal_escape_prediction.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        colors = files.colors,
        auspice_config = files.auspice_config
    output:
        auspice_tree = "auspice/lassa_{segment}_tree.json",
        auspice_meta = "auspice/lassa_{segment}_meta.json"
    shell:
        """
        augur export v1 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.escape_predictions} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_meta}
        """
