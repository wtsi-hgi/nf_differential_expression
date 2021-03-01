#!/usr/bin/env nextflow


def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}


process get_cell_label_list {
    // Get all of the cell labels in an anndata file
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    input:
        path(anndata)
        val(anndata_cell_label)

    output:
        path("cell_labels.csv", emit: cell_labels)

    script:
        runid = random_hex(16)
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "get_cell_label_list: ${process_info}"
        010-get_cell_label_list.py \
            --h5_anndata ${anndata} \
            --cell_label ${anndata_cell_label}
        """
}


process run_differential_expression {
    // Run differential expression
    // ------------------------------------------------------------------------
    label 'long_job'
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(anndata)
        val(cell_label_column)
        val(experiment_id)
        val(mean_cp10k_filter)
        each cell_label
        each model

    output:
        val(outdir, emit: outdir)
        tuple(
            val(runid), //need random hex to control grouping
            val(variable_target),
            val(cell_label),
            val(formula_clean),
            val(model.method),
            path("${outfile}_unfiltered-de_results.tsv.gz"),
            path("${outfile}_filtered-de_results.tsv.gz"),
            val(outdir),
            emit: results
        ) optional true
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        // cell_label comes in array-format.
        cell_label = cell_label[0] // Get first element.
        formula = "${model.formula}"
        formula_clean = "${model.formula}".replaceAll("_", "")
        formula_clean = "${formula_clean}".replaceAll(" ", "_")
        formula_clean = "${formula_clean}".replaceAll("~", "")
        formula_clean = "${formula_clean}".replaceAll("\\+", "_plus_")
        formula_clean = "${formula_clean}".replaceAll("_", "")
        formula_clean = "${formula_clean}".replaceAll("\\)", "")
        formula_clean = "${formula_clean}".replaceAll("\\(.*\\|", "ra_")
        formula_clean = "${formula_clean}".replaceAll("I\\(", "")
        formula_clean = "${formula_clean}".replaceAll("\\^", "power")
        formula_clean = "${formula_clean}".replaceAll("\\.", "pt")
        formula_clean = "${formula_clean}".replaceAll("\\/", "_div_")
        variable_target = "${model.variable_target}"
        variable_target_clean = "${model.variable_target}".replaceAll(
            "\\)", ""
        )
        variable_target_clean = "${variable_target_clean}".replaceAll(
            "I\\(", ""
        )
        variable_target_clean = "${variable_target_clean}".replaceAll(
            "\\^", "power"
        )
        variable_target_clean = "${variable_target_clean}".replaceAll(
            "\\.", "pt"
        )
        outdir = "${outdir_prev}/differential_expression/${variable_target_clean}"
        outdir = "${outdir}/cell_label=${cell_label}"
        outdir = "${outdir}/method=${model.method}___formula=${formula_clean}"
        // Sort out any variables that need to be cast
        cmd__varcast = ""
        if (model.variable_discrete != "") {  // add disc cov call
            cmd__varcast = "${cmd__varcast} --discrete_variables \"${model.variable_discrete}\""
        }
        if (model.variable_continuous != "") {  // add contin cov call
            cmd__varcast = "${cmd__varcast} --continuous_variables \"${model.variable_continuous}\""
        }
        // Make discrete levels command
        cmd__levels = ""
        if (model.variable_discrete_level != "") {  // add disc cov call
            cmd__levels = "--discrete_levels \"${model.variable_discrete_level}\""
        }
        outfile = "cell_label__${cell_label}"
        // Finally get the correct script
        base_method = model.method.split("::")[0]
        if (base_method == "edger") {
            method_script = "011-run_edger.R"
        } else if (base_method == "deseq") {
            method_script = "011-run_deseq.R"
        } else if (base_method == "mast") {
            method_script = "011-run_mast.R"
        }
        // Details on process
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "run_differential_expression: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        convert_h5ad_R.py \
            --h5ad_file ${anndata} \
            --variable_target "${variable_target}" \
            --cell_label_column "${cell_label_column}" \
            --cell_label "${cell_label}" \
            --experiment_key "${experiment_id}" \
            --output_dir de_input \
            ${cmd__varcast}
        011-run_differential_expression.R \
            --input_dir de_input \
            --cell_label_column "${cell_label_column}" \
            --cell_label "${cell_label}" \
            --experiment_key "${experiment_id}" \
            --formula "${model.formula}" \
            --variable_target "${variable_target}" \
            --method "${model.method}" \
            --method_script $baseDir/bin/${method_script} \
            --mean_cp10k_filter ${mean_cp10k_filter} \
            --out_file "${outfile}" \
            --cores_available ${task.cpus} \
            ${cmd__varcast} \
            ${cmd__levels}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process plot_dge_results {
    // Plot DGE results
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        tuple(
            val(runid),
            val(variable_target),
            val(cell_label),
            val(formula_clean),
            val(method),
            file(de_results_unfiltered),
            file(de_results_filtered),
            val(outdir_prev)
        )

    output:
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = outdir_prev
        outfile = "cell_label__${cell_label}"
        // Details on process
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_dge_results: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        012-plot_dge_results.R \
            --input_file ${de_results_unfiltered} \
            --target_var 'coef_value' \
            --sig_var 'qvalue_bh_percelltype' \
            --sig_threshold 0.05 \
            --sig_label 'FDR <= 0.05' \
            --out_file "${outfile}_unfiltered" \
            --verbose
        012-plot_dge_results.R \
            --input_file ${de_results_filtered} \
            --target_var 'coef_value' \
            --sig_var 'qvalue_bh_percelltype' \
            --sig_threshold 0.05 \
            --sig_label 'FDR <= 0.05' \
            --out_file "${outfile}_filtered" \
            --verbose
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

process serialize_de_files {
    // Serializes known markers for analysis
    // ------------------------------------------------------------------------
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    input:
        tuple(
            val(runid),
            val(variable_target),
            val(cell_label),
            val(formula),
            val(method),
            file(de_results_unfiltered),
            file(de_results_filtered),
            val(outdir_prev)
        )

    output:
        tuple(
            val(runid),
            val(variable_target),
            val(cell_label),
            val(formula),
            val(method),
            path("${runid}-${de_results_unfiltered}"),
            path("${runid}-${de_results_filtered}"),
            val(outdir),
            emit: results
        )

    script:
        //runid = random_hex(16)
        outdir = outdir_prev
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "serialize_de_files: ${process_info}"
        ln --physical ${de_results_unfiltered} ${runid}-${de_results_unfiltered}
        ln --physical ${de_results_filtered} ${runid}-${de_results_filtered}
        """
}

process merge_de_dataframes {
    // Merge resulting dataframes
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        val(run_correction)
        val(correction_config)
        tuple(
            val(condition),
            val(result_keys),
            file(result_paths)
        )

    output:
        val(outdir, emit: outdir)
        tuple(
            val(condition),
            path("${outfile}-de_results.tsv.gz"),
            emit: merged_results
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/differential_expression/${condition}/"
        outfile = "${condition}_merged"
        result_keys = result_keys.join(",")
        result_paths = result_paths.join(",")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "merge_de_dataframes: ${process_info}"
        echo "publish_directory: ${outdir}"
        merge_dataframes.py \
            --dataframe_keys ${result_keys} \
            --dataframe_paths ${result_paths} \
            --output_file  '${outfile}-de_results.tsv.gz'
        012-correct_pvals.R \
            --de_results '${outfile}-de_results.tsv.gz' \
            --grouping_cols 'de_method,formula_passed,coef_value' \
            --output_file '${outfile}-de_results.tsv.gz'
        if ${run_correction}; then
            012-ihw_correction.R \
                --de_results "${outfile}-de_results.tsv.gz" \
                --covariates ${correction_config.covariates} \
                --alpha ${correction_config.alpha} \
                --grouping_cols 'de_method,formula_passed' \
                --output_file "${outfile}-de_results.tsv.gz"
        fi
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process plot_merged_dge {
    // Generate plots from the merged data frames to evaluate
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        tuple(
            val(condition),
            path(merged_df)
        )
        each mean_expression_filter

    output:
        val(outdir, emit: outdir)
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/differential_expression/${condition}/"
        outfile = "${condition}_merged_de"
        // script automatically adds expression filter
        // outfile = "${condition}-mean_expr_filt__${mean_expression_filter}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_merged_dge: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        013-compare_de_results.py \
            --dataframe ${merged_df} \
            --columns_to_compare de_method,formula_passed \
            --mean_expression_filter ${mean_expression_filter} \
            --output_file ${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process run_idea {
    // Run iDEA for each DE result
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script
    label 'long_job'

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        tuple(
            val(runid),
            val(condition_column),
            val(cell_label),
            val(covariate_columns),
            val(method),
            file(de_results_unfiltered),
            file(de_results_filtered),
            val(outdir_prev)
        )
        val(max_beta_var)
        val(min_cr)
        each model

    output:
        tuple(
            val(runid),
            val(condition_column),
            val(cell_label),
            val(covariate_columns),
            val(method),
            path("${outfile}-gsea_results.tsv.gz"),
            val(idea_key),
            val(outdir),
            emit: results
        )
        path("${outfile}-gsea_results.Rds.gz")
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true
        path("plots/*.jpg") optional true

    script:
        runid = random_hex(16)
        // Generate key to represent iDEA run to group results
        idea_key = "max_bv=${max_beta_var}"
        idea_key = "${idea_key}::min_cr=${min_cr}"
        idea_key = "${idea_key}::tau_init=${model.tau}"
        idea_key = "${idea_key}::em=${model.em_iteration}"
        idea_key = "${idea_key}::mcmc=${model.mcmc_iteration}"
        idea_key = "${idea_key}::fit_tol=${model.fit_tolerance}"
        idea_key = "${idea_key}::fdr_n=${model.fdr_n_permutations}"
        idea_key = "${idea_key}::db=${model.database}"
        cmd_null = "" // Don't pass null args in
        if (model.beta != "") {
            cmd_null = "${cmd_null} --beta '${model.beta}'"
            idea_key = "${idea_key}::beta_init=${model.beta}"
        }
        idea_key = idea_key.replaceAll(",", "_")
        outdir = outdir_prev
        outdir = "${outdir}/idea"
        outdir = "${outdir}-${idea_key}"
        outfile = "cell_label__${cell_label}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "run_idea: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        015-run_idea.R \
            --de_results ${de_results_filtered} \
            --group_var 'coef_value' \
            --max_beta_var ${max_beta_var} \
            --min_cr ${min_cr} \
            --tau '${model.tau}' \
            --em_iter ${model.em_iteration} \
            --mcmc_iter ${model.mcmc_iteration} \
            --fit_toler '${model.fit_tolerance}' \
            --fdr_n_permutations ${model.fdr_n_permutations} \
            --database '${model.database}' \
            --n_cores ${task.cpus} \
            --output_file ${outfile} \
            --verbose \
            ${cmd_null}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        mv *jpg plots/ 2>/dev/null || true
        """
}

process run_fgsea {
    // Run fGSEA for each DE result
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script
    label 'long_job'

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        tuple(
            val(runid),
            val(condition_column),
            val(cell_label),
            val(covariate_columns),
            val(method),
            file(de_results_unfiltered),
            file(de_results_filtered),
            val(outdir_prev)
        )
        val(sample_size)
        val(score_type)
        each model
        each signed

    output:
        tuple(
            val(runid),
            val(condition_column),
            val(cell_label),
            val(covariate_columns),
            val(method),
            path("${outfile}-gsea_results.tsv.gz"),
            val(fgsea_key),
            val(outdir),
            emit: results
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true
        path("plots/*.jpg") optional true

    script:
        runid = random_hex(16)
        // Generate key to represent iDEA run to group results
        fgsea_key = "sample_size=${sample_size}"
        fgsea_key = "${fgsea_key}::score_type=${score_type}"
        fgsea_key = "${fgsea_key}::min_set_size=${model.min_set_size}"
        fgsea_key = "${fgsea_key}::max_set_size=${model.max_set_size}"
        fgsea_key = "${fgsea_key}::eps=${model.eps}"
        fgsea_key = "${fgsea_key}::db=${model.database}"
        fgsea_key = "${fgsea_key}::${signed}"
        fgsea_key = fgsea_key.replaceAll(",", "_")
        outdir = outdir_prev
        outdir = "${outdir}/fgsea"
        outdir = "${outdir}-${fgsea_key}"
        outfile = "cell_label__${cell_label}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        cmd_null = "" // Don't pass null args in
        if (signed == "unsigned") {
            cmd_null = "${cmd_null} --unsigned_ranking"
        }
        """
        echo "run_fGSEA: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        015-run_fgsea.R \
            --de_results ${de_results_filtered} \
            --group_var 'coef_value' \
            --ranking_var 'test_statistic' \
            --sample_size ${sample_size} \
            --score_type ${score_type} \
            --min_set_size ${model.min_set_size} \
            --max_set_size ${model.max_set_size} \
            --eps '${model.eps}' \
            --gsets_gene_matrix '$baseDir/data/gene_set_gene_matrix.tsv.gz' \
            --gsets_info_file '$baseDir/data/gene_set_info.tsv.gz' \
            --database '${model.database}' \
            --n_cores ${task.cpus} \
            --output_file ${outfile} \
            --verbose \
            ${cmd_null}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        mv *jpg plots/ 2>/dev/null || true
        """
}

process serialize_gsea_files {
    // Serializes known markers for analysis
    // ------------------------------------------------------------------------
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    input:
        tuple(
            val(runid),
            val(condition_column),
            val(cell_label),
            val(covariate_columns),
            val(method),
            path(results_file),
            val(idea_key),
            val(outdir_prev)
        )

    output:
        tuple(
            val(runid),
            val(condition_column),
            val(cell_label),
            val(covariate_columns),
            val(method),
            path("${runid}-${results_file}"),
            val(idea_key),
            val(outdir),
            emit: results
        )

    script:
        //runid = random_hex(16)
        outdir = outdir_prev
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "serialize_gsea_results: ${process_info}"
        echo "publish_directory: ${outdir}"
        ln --physical ${results_file} ${runid}-${results_file}
        """
}

process merge_gsea_dataframes {
    // Merge resulting dataframes from GSEA
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        tuple(
            val(condition),
            val(result_keys),
            file(result_paths)
        )

    output:
        val(outdir, emit: outdir)
        tuple(
            val(condition),
            path("${outfile}-gsea_results.tsv.gz"),
            emit: merged_results
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/differential_expression/${condition}/"
        outfile = "${condition}_merged"
        result_keys = result_keys.join(",")
        result_paths = result_paths.join(",")
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "merge_gsea_dataframes: ${process_info}"
        echo "publish_directory: ${outdir}"
        merge_dataframes.py \
            --dataframe_keys '${result_keys}' \
            --dataframe_paths '${result_paths}' \
            --output_file '${outfile}-gsea_results.tsv.gz'
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

process plot_gsea_results {
    // Generate plots from the merged data frames to evaluate
    // ------------------------------------------------------------------------
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        tuple(
            val(condition),
            path(merged_df)
        )

    output:
        val(outdir, emit: outdir)
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/differential_expression/${condition}/"
        outfile = "${condition}_merged_gsea"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_gsea_results: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        016-compare_gsea_results.py \
            --dataframe ${merged_df} \
            --columns_to_compare de_method,gsea_method,formula_passed,signed_ranking \
            --output_file ${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

workflow wf__differential_expression {
    take:
        outdir
        anndata
        anndata_cell_label
        experiment_key
        model
        de_merge_config
        de_plot_config
        fgsea_config
        idea_config
    main:
        // Get a list of all of the cell types
        get_cell_label_list(
            anndata,
            anndata_cell_label
        )
        // For each cell type compute differential expression for that cell
        // type
        cell_labels = get_cell_label_list.out.cell_labels
            .splitCsv(header: false, sep: ',')
            // .map{row -> tuple(
            //     row[0]
            // )}

        // Run DGE
        run_differential_expression(
            outdir,
            anndata,
            anndata_cell_label,
            experiment_key,
            model.mean_cp10k_filter,
            // '1',  // just run on first cluster for development
            cell_labels,  // run for all clusters for run time
            model.value
        )

        de_results = run_differential_expression.out.results

        // Plot results from each DGE run
        plot_dge_results(
            de_results
        )

        // Serialize input files to prep for merge
        serialize_de_files(
            de_results
        )

        // Group results by formula for merge
        de_results_merged = serialize_de_files.out.results
            .reduce([:]) { map, tuple ->
                def dataframe_key = "cell_label=" + tuple[2]
                dataframe_key += "::formula=" + tuple[3].replaceAll(
                    ",",
                    "-"
                )
                dataframe_key += "::method=" + tuple[4]

                def map_key = tuple[1] // structure map by condition
                def key_list = map[map_key]
                // Right now, we only want to carry down the filtered data
                if (!key_list) {
                    key_list = [[dataframe_key, tuple[6]]]
                } else {
                    key_list.add([dataframe_key, tuple[6]])
                }
                map[map_key] = key_list
                return(map)
            }
            .flatMap()
            .map {  entry ->
                combined_data = [entry.key, [], []]
                entry.value.each {
                    combined_data[1].add(it[0])
                    combined_data[2].add(it[1])
                }
                return(combined_data)
            }
        merge_de_dataframes(
            outdir,
            de_merge_config.ihw_correction.run_process,
            de_merge_config.ihw_correction.value,
            de_results_merged
        )
        // Basic plots of the differential expression results across all models
        plot_merged_dge(
            outdir,
            merge_de_dataframes.out.merged_results,
            de_plot_config.mean_expression_filter.value
        )

        // Run fGSEA on DE results
        if (fgsea_config.run_process) {
            run_fgsea(
                de_results,
                fgsea_config.sample_size,
                fgsea_config.score_type,
                fgsea_config.value,
                ["unsigned", "signed"]
            )
        }

        // Run iDEA on DE results
        if (idea_config.run_process) {
            run_idea(
                de_results,
                idea_config.max_beta_variance,
                idea_config.min_cover_rate,
                idea_config.value
            )
        }

        // Combine results of different GSEA methods
        gsea_results = null
        if (fgsea_config.run_process) {
            gsea_results = run_fgsea.out.results
        }
        if (idea_config.run_process) {
            if (gsea_results) {
                gsea_results = gsea_results.concat(run_idea.out.results)
            } else {
                gsea_results = run_idea.out.results
            }
        }

        // Serialize input files to prep for merge
        serialize_gsea_files(
            gsea_results
        )

        // Combine and compare all of the enrichment analysis results
        gsea_results_merged = serialize_gsea_files.out.results
            .reduce([:]) { map, tuple ->
                def dataframe_key = "cell_label=" + tuple[2]
                dataframe_key += "::covariates=" + tuple[3].replaceAll(
                    ",",
                    "-"
                )
                dataframe_key += "::method=" + tuple[4]
                dataframe_key += "__gsea_params=" + tuple[6]

                def map_key = tuple[1] // structure map by condition
                def key_list = map[map_key]
                if (!key_list) {
                    key_list = [[dataframe_key, tuple[5]]]
                } else {
                    key_list.add([dataframe_key, tuple[5]])
                }
                map[map_key] = key_list
                return(map)
            }
            .flatMap()
            .map {  entry ->
                combined_data = [entry.key, [], []]
                entry.value.each {
                    combined_data[1].add(it[0])
                    combined_data[2].add(it[1])
                }
                return(combined_data)
            }

        merge_gsea_dataframes(
            outdir,
            gsea_results_merged
        )

        // Basic plots of the GSEA results across all models
        plot_gsea_results(
            outdir,
            merge_gsea_dataframes.out.merged_results
        )

    emit:
        cell_labels = get_cell_label_list.out.cell_labels
}
