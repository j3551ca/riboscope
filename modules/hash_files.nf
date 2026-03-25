process hash_files {

    tag { sample_id + " / " + file_type }

    input:
    tuple  val(sample_id), path(files_to_hash), val(file_type)

    output:
    tuple val(sample_id), path(files_to_hash), path("${sample_id}_${file_type}.sha256.csv"), emit: hashes

    script:
    """
    shasum -a 256 ${files_to_hash} | tr -s ' ' ',' > ${sample_id}_${file_type}.sha256.csv
    """
}

process print_hashed_records {

    input:
    tuple val(sample_id), val(file_type), val(sha256), val(filename)

    script:
    """
    echo "for provenance: ${sample_id},${file_type},${filename},${sha256}"
    """
}