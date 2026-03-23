process pipeline_provenance {

    tag { pipeline_name + " / " + pipeline_version }

    executor 'local'

    input:
    tuple val(nextflow_version), val(session_id), val(run_name), val(pipeline_name), val(pipeline_version), 
    val(user), val(timestamp_analysis_start), val(command_line), val(launch_dir), val(project_dir), val(work_dir),
    val(repo), val(commit_id), val(branch)

    script:
    """
    printf -- "- pipeline_name: ${pipeline_name}"
    printf -- "  pipeline_version: ${pipeline_version}\\n"
    printf -- "- nextflow_version: ${nextflow_version}\\n"
    printf -- "  nextflow_session_id: ${session_id}\\n"
    printf -- "  nextflow_run_name: ${run_name}\\n"
    printf -- "- timestamp_analysis_start: ${timestamp_analysis_start}\\n"
    printf -- "  user: ${user}\\n"
    printf -- "  command_executed: ${command_line}\\n"
    printf -- "  launch_directory: ${launch_dir}\\n"
    printf -- "  project_directory: ${project_dir}\\n"
    printf -- "  working_directory: ${work_dir}\\n"
    printf -- "- git_repo: ${repo}\\n"
    printf -- "  git_commit: ${commit_id}\\n"
    printf -- "  git_branch: ${branch}\\n"
    """
}

process collect_provenance {

    tag { sample_id }
    
    executor 'local'

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_*_provenance.yml", mode: 'copy'

    input:
    tuple val(sample_id), path(provenance_files)

    output:
    tuple val(sample_id), file("${sample_id}_*_provenance.yml")

    script:
    """
    cat ${provenance_files} > ${sample_id}_\$(date +%Y%m%d%H%M%S)_provenance.yml
    """
}
