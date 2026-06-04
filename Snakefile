# --- Config ---
# Load all parameters from the config file
configfile: "config.yaml"

# Import required modules
from itertools import combinations
import pandas as pd
import os

# Load .env if present (gitignored) so users can set OUTPUT_FOLDER to an absolute path
_env_path = os.path.join(os.path.dirname(workflow.snakefile), ".env")
if os.path.exists(_env_path):
    with open(_env_path) as _f:
        for _line in _f:
            _line = _line.strip()
            if _line and not _line.startswith("#") and "=" in _line:
                _k, _v = _line.split("=", 1)
                os.environ.setdefault(_k.strip(), _v.strip())

OUTPUT_FOLDER = os.environ.get("OUTPUT_FOLDER", config["output_folder"])
DATA_FOLDER = config["data_folder"]
SCRIPTS_FOLDER = config["scripts_folder"]

# Generate the list of seeds based on the config parameters. [count] seeds, starting at [start]
split_seeds = list(range(config["split_seed_start"], config["split_seed_start"] + config["split_seed_count"]))

wildcard_constraints:
    seed=r"\d+",
    n_datasets=r"\d+",
    test_study=r"[a-zA-Z0-9_]+",  # Study names like GSE37250_SA, USA, India
    mean=r"[\d\.\-]+", 
    variance=r"[\d\.]+",        # Matches numbers and dots
    classifier=r"[a-zA-Z0-9]+",
    adjuster=r"[a-zA-Z0-9_]+",    # Underscores allowed in adjuster names
    scenario_name=r"[a-zA-Z0-9_\-]+",  # Scenario names with underscores and dashes
    replicate=r"\d+"


# --- Final Target ---
# Define all available studies
ALL_STUDIES = ["GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M"]

# Helper function to get valid test studies for a given n_datasets
# n = number of TRAINING studies (test study is excluded from training)
# Any study can be test for any n < len(ALL_STUDIES)
def get_test_studies_for_n(n):
    return ALL_STUDIES  # All studies valid as test; pipeline excludes test from step-1 training

# Helper to generate all adjusted data combinations for specified adjusters
def get_all_adjusted_data_combinations():
    """Generate all adjusted data files for combat & combat_sup × n_datasets × test_study combinations"""
    # Only generate adjusted data for ComBat and ComBat supervised for investigation
    focus_adjusters = ["combat", "combat_sup"]
    combinations = []
    for adjuster in focus_adjusters:
        for n in config["num_datasets"]:
            for test_study in get_test_studies_for_n(n):
                combinations.append(
                    f"{OUTPUT_FOLDER}/adjusted_data/all_scenarios/{adjuster}_n{n}_test{test_study}_target.csv"
                )
    return combinations

# Generate all valid combinations
VIZ_TARGETS = []
for n in config["num_datasets"]:
    for test_study in get_test_studies_for_n(n):
        for method in ["pca", "lda", "umap"]:
            VIZ_TARGETS.append(
                f"{OUTPUT_FOLDER}/visualizations/grids/{method}_n{n}_test{test_study}.png"
            )

rule all:
    input:
        adjuster_main=OUTPUT_FOLDER + "/adjusters_on_classifiers.png",
        adjuster_relative_agg=OUTPUT_FOLDER + "/adjusters_on_classifiers_relative.png",
        avg_rank_study=OUTPUT_FOLDER + "/average_rank_by_classifier.png",
        mcc_rank_by_adjuster=OUTPUT_FOLDER + "/mcc_rank_by_adjuster.png",
        mcc_rank_by_adjuster_wide=OUTPUT_FOLDER + "/mcc_rank_by_adjuster_wide.png",
        mean_over_classifiers=OUTPUT_FOLDER + "/mean_over_classifiers.png",
        within_study_baseline=OUTPUT_FOLDER + "/within_study_baseline.csv",
        class_imbalance_results=OUTPUT_FOLDER + "/class_imbalanced_results.csv" if config.get("class_imbalance_analysis", {}).get("enabled", False) else [],
        class_imbalance_trend_plot=OUTPUT_FOLDER + "/class_imbalance_trend.png" if config.get("class_imbalance_analysis", {}).get("enabled", False) else [],
        sample_counts=OUTPUT_FOLDER + "/sample_counts.csv",
        gapr_heatmap=OUTPUT_FOLDER + "/gapr_heatmap.png",
        adjuster_similarity=OUTPUT_FOLDER + "/adjuster_similarity.png",
        adjuster_recommendation=OUTPUT_FOLDER + "/adjuster_recommendation.png",
        schematic_matrix=OUTPUT_FOLDER + "/schematic_matrix.png",
        pca_demonstration=OUTPUT_FOLDER + "/pca_demonstration.png",
        sa_uk_adjusted=expand(
            OUTPUT_FOLDER + "/adjusted_data/SA_UK/{adjuster}_target.csv",
            adjuster=config["adjusters"]
        ),
        all_adjusted=get_all_adjusted_data_combinations()





## How does each adjuster help mitigate batch effects, for each classifier?

# Step 1.3a: Create main adjuster plot
rule make_fig2_main:
    input:
        data=OUTPUT_FOLDER + "/adjusters_on_classifiers.csv",
        script=SCRIPTS_FOLDER + "/plot_adjusters_main_only.R",
        utils=SCRIPTS_FOLDER + "/adjuster_plot_utils.R",
        generator=SCRIPTS_FOLDER + "/generate_main_plot.R"
    output:
        OUTPUT_FOLDER + "/adjusters_on_classifiers.png"
    log:
        OUTPUT_FOLDER + "/logs/plot_adjusters_main.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10
    params:
        adjusters=",".join(config.get("viz_adjusters", config["adjusters"]))
    shell:
        """
        Rscript {input.script} -i {input.data} -o {output} --adjusters {params.adjusters} &> {log}
        """

# Step 1.3b: Create aggregated relative performance plot
rule make_fig2_relative_aggregated:
    input:
        data=OUTPUT_FOLDER + "/adjusters_on_classifiers.csv",
        script=SCRIPTS_FOLDER + "/plot_adjusters_relative.R",
        utils=SCRIPTS_FOLDER + "/adjuster_plot_utils.R",
        generator=SCRIPTS_FOLDER + "/generate_relative_plot_aggregated.R"
    output:
        plot=OUTPUT_FOLDER + "/adjusters_on_classifiers_relative.png",
        sig=OUTPUT_FOLDER + "/adjusters_on_classifiers_relative_significance.csv"
    log:
        OUTPUT_FOLDER + "/logs/plot_adjusters_relative.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10
    params:
        adjusters=",".join(config.get("viz_adjusters", config["adjusters"]))
    shell:
        """
        Rscript {input.script} -i {input.data} -o {output.plot} --adjusters {params.adjusters} &> {log}
        """

# Step 1.3c: Create average rank by classifier plot for 4-study case
rule make_avg_rank_study:
    input:
        data=OUTPUT_FOLDER + "/adjusters_on_classifiers.csv",
        script=SCRIPTS_FOLDER + "/plot_average_rank_by_classifier.R",
        utils=SCRIPTS_FOLDER + "/adjuster_plot_utils.R"
    output:
        OUTPUT_FOLDER + "/average_rank_by_classifier.png"
    log:
        OUTPUT_FOLDER + "/logs/plot_average_rank.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10
    params:
        adjusters=",".join(config.get("viz_adjusters", config["adjusters"])),
        n_datasets="all"
    shell:
        """
        Rscript {input.script} -i {input.data} -o {output} --adjusters {params.adjusters} --n-datasets {params.n_datasets} &> {log}
        """

# Step 1.3d: Create MCC rank by adjuster plot (adjuster on x, rank on y, colored by classifier)
rule make_mcc_rank_by_adjuster:
    input:
        data=OUTPUT_FOLDER + "/adjusters_on_classifiers.csv",
        script=SCRIPTS_FOLDER + "/plot_mcc_rank_by_adjuster.R",
        utils=SCRIPTS_FOLDER + "/adjuster_plot_utils.R"
    output:
        skinny=OUTPUT_FOLDER + "/mcc_rank_by_adjuster.png",
        wide=OUTPUT_FOLDER + "/mcc_rank_by_adjuster_wide.png"
    log:
        OUTPUT_FOLDER + "/logs/plot_mcc_rank_by_adjuster.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10
    params:
        adjusters=",".join(config.get("viz_adjusters", config["adjusters"])),
        n_datasets="all"
    shell:
        """
        Rscript {input.script} -i {input.data} -o {output.skinny} --adjusters {params.adjusters} --n-datasets {params.n_datasets} &> {log}
        Rscript {input.script} -i {input.data} -o {output.wide} --adjusters {params.adjusters} --n-datasets {params.n_datasets} --width 8 --height 7 &>> {log}
        """

# Step 1.3e: Mean performance gap averaged over classifiers
rule make_mean_over_classifiers:
    input:
        data=OUTPUT_FOLDER + "/adjusters_on_classifiers.csv",
        script=SCRIPTS_FOLDER + "/plot_mean_over_classifiers.R",
        utils=SCRIPTS_FOLDER + "/adjuster_plot_utils.R"
    output:
        OUTPUT_FOLDER + "/mean_over_classifiers.png"
    log:
        OUTPUT_FOLDER + "/logs/plot_mean_over_classifiers.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10
    params:
        adjusters=",".join(config.get("viz_adjusters", config["adjusters"]))
    shell:
        """
        Rscript {input.script} -i {input.data} -o {output} --adjusters {params.adjusters} &> {log}
        """

# Step 1.3f: GAPR Heatmap with R2E seriation
rule make_gapr_heatmap:
    input:
        data=OUTPUT_FOLDER + "/adjusters_on_classifiers.csv",
        script=SCRIPTS_FOLDER + "/plot_gapr_heatmap.R",
        utils=SCRIPTS_FOLDER + "/adjuster_plot_utils.R"
    output:
        OUTPUT_FOLDER + "/gapr_heatmap.png"
    log:
        OUTPUT_FOLDER + "/logs/plot_gapr_heatmap.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10
    params:
        adjusters=",".join(config.get("viz_adjusters", config["adjusters"]))
    shell:
        """
        Rscript {input.script} -i {input.data} -o {output} --adjusters {params.adjusters} &> {log}
        """
# Step 1.3g: Adjuster Similarity Heatmap
rule make_adjuster_similarity:
    input:
        data=OUTPUT_FOLDER + "/adjusters_on_classifiers.csv",
        script=SCRIPTS_FOLDER + "/plot_adjuster_similarity.R",
        utils=SCRIPTS_FOLDER + "/adjuster_plot_utils.R"
    output:
        OUTPUT_FOLDER + "/adjuster_similarity.png"
    log:
        OUTPUT_FOLDER + "/logs/plot_adjuster_similarity.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10
    params:
        adjusters=",".join(config.get("viz_adjusters", config["adjusters"]))
    shell:
        """
        Rscript {input.script} -i {input.data} -o {output} --adjusters {params.adjusters} &> {log}
        """


# Step 1.3h: Adjuster Recommendation Plot (Faceted by Classifier)
rule make_adjuster_recommendation:
    input:
        data=OUTPUT_FOLDER + "/adjusters_on_classifiers.csv",
        script=SCRIPTS_FOLDER + "/plot_adjuster_recommendation.R",
        utils=SCRIPTS_FOLDER + "/adjuster_plot_utils.R"
    output:
        OUTPUT_FOLDER + "/adjuster_recommendation.png"
    log:
        OUTPUT_FOLDER + "/logs/plot_adjuster_recommendation.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10
    params:
        adjusters=",".join(config.get("viz_adjusters", config["adjusters"])),
        top_n=10
    shell:
        """
        Rscript {input.script} -i {input.data} -o {output} --adjusters {params.adjusters} --top-n {params.top_n} &> {log}
        """

# Step 1.3i: Schematic Data Matrix (Presentation Diagram)
rule make_schematic_matrix:
    input:
        data=DATA_FOLDER + "/TB_real_data.RData",
        script=SCRIPTS_FOLDER + "/plot_schematic_matrix.R"
    output:
        OUTPUT_FOLDER + "/schematic_matrix.png"
    log:
        OUTPUT_FOLDER + "/logs/plot_schematic_matrix.log"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=5
    shell:
        """
        Rscript {input.script} -i {input.data} -o {output} &> {log}
        """


# Step 1.3j: PCA Demonstration (Before vs After)
rule make_pca_demonstration:
    input:
        data=DATA_FOLDER + "/TB_real_data.RData",
        script=SCRIPTS_FOLDER + "/plot_pca_demonstration.R"
    output:
        OUTPUT_FOLDER + "/pca_demonstration.png"
    log:
        OUTPUT_FOLDER + "/logs/plot_pca_demonstration.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10
    shell:
        """
        Rscript {input.script} -i {input.data} -o {output} &> {log}
        """


# Helper to generate all valid adjuster combinations
def get_adjuster_combinations():
    """Generate all valid adjuster × classifier × n_datasets × test_study combinations"""
    combinations = []
    for adjuster in config["adjusters"]:
        for classifier in config["classifiers"]:
            for n in config["num_datasets"]:
                for test_study in get_test_studies_for_n(n):
                    combinations.append(
                        f"{OUTPUT_FOLDER}/results/adjusters/individual/{adjuster}_{classifier}_{n}_{test_study}.csv"
                    )
    return combinations

# Helper to generate within-study CV combinations
def get_within_study_cv_combinations():
    """Generate all valid classifier × n_datasets × test_study combinations for within-study CV"""
    combinations = []
    for classifier in config["classifiers"]:
        for n in config["num_datasets"]:
            for test_study in get_test_studies_for_n(n):
                combinations.append(
                    f"{OUTPUT_FOLDER}/results/within_study_cv/individual/{classifier}_{n}_{test_study}.csv"
                )
    return combinations

# Helper to generate class imbalance scenario names statically from config.
# All C(6,2)=15 training pairs × all imbalance levels × remaining test studies × replicates.
# The imbalance_levels in config MUST be chosen so every combination is feasible;
# the R script will fail hard if any is not.
def get_class_imbalance_scenario_names():
    """Compute all scenario names from config — no file I/O at DAG-build time."""
    ci_config = config.get("class_imbalance_analysis", {})
    if not ci_config.get("enabled", False):
        return []

    imbalance_levels = ci_config.get("imbalance_levels", [])
    n_replicates = ci_config.get("n_replicates", 1)
    scenarios = []

    # All ordered pairs from combinations(ALL_STUDIES, 2)
    for i, s1 in enumerate(ALL_STUDIES):
        for s2 in ALL_STUDIES[i+1:]:
            pair_name = f"{s1}-{s2}"
            test_datasets = [s for s in ALL_STUDIES if s not in (s1, s2)]
            for imbal in imbalance_levels:
                imbal_int = int(imbal * 100)
                for test_ds in test_datasets:
                    for rep in range(1, n_replicates + 1):
                        scenarios.append(f"{pair_name}-imbal{imbal_int}-test{test_ds}-rep{rep}")

    return scenarios

# Pre-compute scenario names once (purely from config, no file reads)
CLASS_IMBALANCE_SCENARIOS = get_class_imbalance_scenario_names()

# Step 1.2: Aggregate all individual adjuster results (including within-study CV baseline)
rule aggregate_adjusters:
    input:
        adjuster_data=get_adjuster_combinations(),
        within_study_data=get_within_study_cv_combinations(),
        script=SCRIPTS_FOLDER + "/aggregate_results.R"
    output:
        OUTPUT_FOLDER + "/adjusters_on_classifiers.csv"
    log:
        OUTPUT_FOLDER + "/logs/aggregate_adjusters.log"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=20
    params:
        adjuster_dir=OUTPUT_FOLDER + "/results/adjusters/individual",
        within_study_dir=OUTPUT_FOLDER + "/results/within_study_cv/individual"
    shell:
        """
        Rscript {input.script} \
            --input-dir {params.adjuster_dir} \
            --input-dir {params.within_study_dir} \
            -o {output} \
            &> {log}
        """


# Step 1.1a: Run classification on real data with batch correction
rule classify_adjusters:
    input:
        data=DATA_FOLDER + "/TB_real_data.RData",
        # script=SCRIPTS_FOLDER + "/classify_adjusters.R",
        helper=SCRIPTS_FOLDER + "/helper.R"
    output:
        OUTPUT_FOLDER + "/results/adjusters/individual/{adjuster}_{classifier}_{n_datasets}_{test_study}.csv"
    log:
        OUTPUT_FOLDER + "/logs/classify_adjusters/{adjuster}_{classifier}_{n_datasets}_{test_study}.log"
    # group: "batch_real_group"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=300
    shell:
        """
        Rscript {SCRIPTS_FOLDER}/classify_adjusters.R -o {output} \
            --adjuster {wildcards.adjuster} \
            --classifier {wildcards.classifier} \
            --num-datasets {wildcards.n_datasets} \
            --test-study {wildcards.test_study} \
            |& tee {log}
        """


# Step 1.1b: Within-study cross-validation baseline (no batch effects)
rule classify_within_study_cv:
    input:
        data=DATA_FOLDER + "/TB_real_data.RData",
        script=SCRIPTS_FOLDER + "/classify_within_study_cv.R",
        helper=SCRIPTS_FOLDER + "/helper.R"
    output:
        OUTPUT_FOLDER + "/results/within_study_cv/individual/{classifier}_{n_datasets}_{test_study}.csv"
    log:
        OUTPUT_FOLDER + "/logs/within_study_cv/{classifier}_{n_datasets}_{test_study}.log"
    group: "within_study_cv_group"
    threads: 1
    resources:
        mem_mb=lambda wildcards: 32000 if wildcards.classifier in ["nnet", "logistic"] else 8000,
        runtime=180
    params:
        n_folds=3,
        n_features=2000  # 0 = use all genes (no feature selection)
    shell:
        """
        Rscript {input.script} -o {output} \
            --classifier {wildcards.classifier} \
            --num-datasets {wildcards.n_datasets} \
            --test-study {wildcards.test_study} \
            --n-folds {params.n_folds} \
            --n-features {params.n_features} \
            &> {log}
        """


# Step 1.1c: Aggregate within-study CV results for summary
rule aggregate_within_study_baseline:
    input:
        data=get_within_study_cv_combinations(),
        script=SCRIPTS_FOLDER + "/aggregate_results.R"
    output:
        OUTPUT_FOLDER + "/within_study_baseline.csv"
    log:
        OUTPUT_FOLDER + "/logs/aggregate_within_study_baseline.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10
    params:
        input_dir=OUTPUT_FOLDER + "/results/within_study_cv/individual"
    shell:
        """
        Rscript {input.script} \
            --input-dir {params.input_dir} -o {output} \
            &> {log}
        """


# Step 1.1d: Report sample counts per study
rule report_sample_counts:
    input:
        data=DATA_FOLDER + "/TB_real_data.RData",
        script=SCRIPTS_FOLDER + "/report_sample_counts.R"
    output:
        OUTPUT_FOLDER + "/sample_counts.csv"
    log:
        OUTPUT_FOLDER + "/logs/report_sample_counts.log"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=5
    shell:
        """
        Rscript {input.script} -i {input.data} -o {output} &> {log}
        """


## Visualization Pipeline: Visualize batch adjustment effects

# Step 2.1: Generate individual visualizations for each adjuster
rule visualize_batch_adjustment:
    input:
        data=DATA_FOLDER + "/TB_real_data.RData",
        script=SCRIPTS_FOLDER + "/visualize_batch_adjustment.R",
        helper=SCRIPTS_FOLDER + "/helper.R"
    output:
        pca=OUTPUT_FOLDER + "/visualizations/individual/pca/{adjuster}_n{n_datasets}_test{test_study}.png",
        lda=OUTPUT_FOLDER + "/visualizations/individual/lda/{adjuster}_n{n_datasets}_test{test_study}.png",
        umap=OUTPUT_FOLDER + "/visualizations/individual/umap/{adjuster}_n{n_datasets}_test{test_study}.png"
    log:
        OUTPUT_FOLDER + "/logs/visualize_batch_adjustment/{adjuster}_n{n_datasets}_test{test_study}.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=17
    params:
        output_dir=OUTPUT_FOLDER + "/visualizations/individual"
    shell:
        """
        Rscript {input.script} \
            --adjuster {wildcards.adjuster} \
            --num-datasets {wildcards.n_datasets} \
            --test-study {wildcards.test_study} \
            --output-dir {params.output_dir} \
            --reduce 0 \
            &> {log}
        """


# Custom pipeline for South Africa vs UK batch adjustment and target datasets
rule adjust_target_data_SA_UK:
    input:
        data=DATA_FOLDER + "/TB_real_data.RData",
        script=SCRIPTS_FOLDER + "/adjust_target_data_SA_UK.R",
        helper=SCRIPTS_FOLDER + "/helper.R"
    output:
        adjusted_target=OUTPUT_FOLDER + "/adjusted_data/SA_UK/{adjuster}_target.csv"
    log:
        OUTPUT_FOLDER + "/logs/adjust_target_data_SA_UK/{adjuster}.log"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10
    shell:
        """
        Rscript {input.script} \
            --adjuster {wildcards.adjuster} \
            --output-data {output.adjusted_target} \
            &> {log}
        """


# Save adjusted data for all scenarios (all n_datasets and test_study combinations)
rule adjust_target_data_all_scenarios:
    input:
        data=DATA_FOLDER + "/TB_real_data.RData",
        script=SCRIPTS_FOLDER + "/adjust_target_data_SA_UK.R",
        helper=SCRIPTS_FOLDER + "/helper.R"
    output:
        adjusted_target=OUTPUT_FOLDER + "/adjusted_data/all_scenarios/{adjuster}_n{n_datasets}_test{test_study}_target.csv"
    log:
        OUTPUT_FOLDER + "/logs/adjust_target_data_all_scenarios/{adjuster}_n{n_datasets}_test{test_study}.log"
    priority: 100
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10
    shell:
        """
        Rscript {input.script} \
            --adjuster {wildcards.adjuster} \
            --num-datasets {wildcards.n_datasets} \
            --test-study {wildcards.test_study} \
            --output-data {output.adjusted_target} \
            &> {log}
        """


# Step 2.2: Create comparison grids
rule create_visualization_grid:
    input:
        plots=expand(
            OUTPUT_FOLDER + "/visualizations/individual/{{method}}/{adjuster}_n{{n_datasets}}_test{{test_study}}.png",
            adjuster=config.get("viz_adjusters", config["adjusters"])
        ),
        script=SCRIPTS_FOLDER + "/create_visualization_grid.R"
    output:
        OUTPUT_FOLDER + "/visualizations/grids/{method}_n{n_datasets}_test{test_study}.png"
    log:
        OUTPUT_FOLDER + "/logs/visualization_grids/{method}_n{n_datasets}_test{test_study}.log"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=8
    params:
        input_dir=OUTPUT_FOLDER + "/visualizations/individual",
        adjusters=",".join(config.get("viz_adjusters", config["adjusters"]))
    shell:
        """
        Rscript {input.script} \
            --input-dir {params.input_dir} \
            --output-file {output} \
            --method {wildcards.method} \
            --num-datasets {wildcards.n_datasets} \
            --test-study {wildcards.test_study} \
            --adjusters {params.adjusters} \
            &> {log}
        """

## Class Imbalance Analysis Pipeline: Systematic class imbalance testing across all dataset combinations

# Helper functions for class imbalance analysis
# Helper to generate all valid class imbalance combinations
def get_class_imbalance_combinations():
    """Generate all valid class imbalance adjuster × classifier × scenario combinations"""
    ci_config = config.get("class_imbalance_analysis", {})
    if not ci_config.get("enabled", False):
        return []

    class_imbalance_adjusters = ci_config.get("adjusters", ["unadjusted", "combat", "combat_sup"])
    class_imbalance_classifiers = ci_config.get("classifiers", config["classifiers"])

    combinations_list = []
    for scenario in CLASS_IMBALANCE_SCENARIOS:
        for adjuster in class_imbalance_adjusters:
            for classifier in class_imbalance_classifiers:
                combinations_list.append(
                    f"{OUTPUT_FOLDER}/results/class_imbalance/individual/{adjuster}-{classifier}-{scenario}.csv"
                )

    return combinations_list

# Step 3.1: Run classification on class imbalanced data (subsetting done inline, no intermediate files)
if config.get("class_imbalance_analysis", {}).get("enabled", False):
    rule classify_class_imbalanced:
        input:
            data=DATA_FOLDER + "/TB_real_data.RData",
            script=SCRIPTS_FOLDER + "/classify_class_imbalanced.R",
            helper=SCRIPTS_FOLDER + "/helper.R"
        output:
            OUTPUT_FOLDER + "/results/class_imbalance/individual/{adjuster}-{classifier}-{scenario_name}.csv"
        log:
            OUTPUT_FOLDER + "/logs/classify_class_imbalanced/{adjuster}-{classifier}-{scenario_name}.log"
        group: "class_imbalance_group"
        priority: 50
        threads: 1
        resources:
            mem_mb=4000,
            runtime=20
        params:
            seed=config.get("class_imbalance_analysis", {}).get("seed", 123),
            imbalance_levels=",".join(str(x) for x in config.get("class_imbalance_analysis", {}).get("imbalance_levels", []))
        shell:
            """
            Rscript {input.script} \
                --adjuster {wildcards.adjuster} \
                --classifier {wildcards.classifier} \
                --scenario-name {wildcards.scenario_name} \
                --data {input.data} \
                --seed {params.seed} \
                --imbalance-levels {params.imbalance_levels} \
                -o {output} \
                &> {log}
            """

    # Step 3.3: Aggregate class imbalanced results
    rule aggregate_class_imbalanced_results:
        input:
            data=get_class_imbalance_combinations(),
            script=SCRIPTS_FOLDER + "/aggregate_class_imbalanced_results.R"
        output:
            OUTPUT_FOLDER + "/class_imbalanced_results.csv"
        log:
            OUTPUT_FOLDER + "/logs/aggregate_class_imbalanced_results.log"
        threads: 1
        resources:
            mem_mb=4000,
            runtime=10
        params:
            input_dir=OUTPUT_FOLDER + "/results/class_imbalance/individual"
        shell:
            """
            Rscript {input.script} \
                --input-dir {params.input_dir} \
                -o {output} \
                &> {log}
            """

    # Step 3.4: Create class imbalance trend plot
    rule plot_class_imbalance_trend:
        input:
            data=OUTPUT_FOLDER + "/class_imbalanced_results.csv",
            script=SCRIPTS_FOLDER + "/plot_class_imbalance_trend.R"
        output:
            OUTPUT_FOLDER + "/class_imbalance_trend.png"
        log:
            OUTPUT_FOLDER + "/logs/plot_class_imbalance_trend.log"
        threads: 1
        resources:
            mem_mb=4000,
            runtime=15
        shell:
            """
            Rscript {input.script} \
                --input-data {input.data} \
                -o {output} \
                &> {log}
            """