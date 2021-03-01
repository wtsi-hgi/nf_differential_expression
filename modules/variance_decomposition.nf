#!/usr/bin/env nextflow


def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}


process run_variance_decomposition {
    // Run variance decomposition
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    //maxForks 2         // hard to control memory usage. limit to 3 concurrent
    //label 'gpu'          // use GPU
    // label 'long_job'
    scratch false        // use tmp directory
    echo echo_mode       // echo output from script

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${task.publish_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(anndata)
        val(covariates)

    output:
        val(outdir, emit: outdir)
        path("variance_decomposed.tsv.gz", emit: decomposed_variance)
        path("plots/*.png") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/variance_decomposition/"
        outfile = "variance_decomposed"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "run_variance_decomposition: ${process_info}"
        rm -fr plots
        010-run_variance_decomposition.py \
            --h5_anndata ${anndata} \
            --covariates ${covariates} \
            --output_file ${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

workflow wf__variance_decomposition {
    take:
        outdir
        anndata
        covariates
        vardecomp_method_config

    main:

        if (vardecomp_method_config.run_process) {
            run_variance_decomposition(
                outdir,
                anndata,
                covariates
            )
        }

}
