params.targetsBed = "targets.bed"
params.samplesheet = "samples.analysis.tsv"
params.numChunks = 5
params.numlines = 100
params.disable = "false"
params.outputDir = "output-${params.numlines}lines"
params.insert_size = 500 // 500bp reported by wet lab for sequencing

Channel.from("${params.numlines}").map { nlines ->
    def lineChunkLabel = "${params.numlines}lineChunk"

    return([ lineChunkLabel, nlines ])
}.into { lineChunkLabel_ch; lineChunkLabel_ch2 }


Channel.from("${params.numChunks}").map { nchunks ->
    def nChunkLabel = "${nchunks}Chunk"

    return([ nChunkLabel, nchunks ])
    }.into { nChunkLabel_ch; nChunkLabel_ch2 }

Channel.fromPath("${params.samplesheet}")
    .splitCsv(header: true, sep: '\t')
    .map { row ->
        def tumorID = "${row.Tumor}"
        def normalID = "${row.Normal}"
        def comparisonID = "${tumorID}_${normalID}"
        def tumorBam = file("${row.Tumor_Bam}")
        def tumorBai = file("${row.Tumor_Bai}")
        def normalBam = file("${row.Normal_Bam}")
        def normalBai = file("${row.Normal_Bai}")

        return([ comparisonID, tumorID, normalID, tumorBam, tumorBai, normalBam, normalBai ])
    }
    .into { input_bams; input_bams2; input_bams3; input_bams4 }

Channel.fromPath( file(params.targetsBed) ).into{ targets_bed; targets_bed2; targets_bed3; targets_bed4 }
Channel.fromPath( file(params.ref_fa) ).into { ref_fasta; ref_fasta2; ref_fasta3; ref_fasta4; ref_fasta5; ref_fasta6; ref_fasta7; ref_fasta8 }
Channel.fromPath( file(params.ref_fai) ).into { ref_fai; ref_fai2; ref_fai3; ref_fai4; ref_fai5; ref_fai6; ref_fai7; ref_fai8 }
Channel.fromPath( file(params.ref_dict) ).into { ref_dict; ref_dict2; ref_dict3; ref_dict4; ref_dict5; ref_dict6; ref_dict7; ref_dict8 }
Channel.fromPath( file(params.ANNOVAR_DB_DIR) ).into { annovar_db_dir; annovar_db_dir2; annovar_db_dir3 }

// get the unique chromosomes in the targets bed file
Channel.fromPath( params.targetsBed )
            .splitCsv(sep: '\t')
            .map{row ->
                row[0]
            }
            .unique()
            .set{ chroms }
// chroms.subscribe { println "${it}" }



// ~~~~~~~ START PIPELINE ~~~~~~~~ //




// ~~~~~~~ NO CHUNK ~~~~~~~~ //
input_bams.combine(targets_bed)
    .combine(ref_fasta)
    .combine(ref_fai)
    .combine(ref_dict)
    .tap { input_noChunk_ch }
process pindel_noChunk {
    tag "${prefix}"
    publishDir "${params.outputDir}/variants", overwrite: true, mode: 'copy'

    input:
    set val(comparisonID), val(tumorID), val(normalID), file(tumorBam), file(tumorBai), file(normalBam), file(normalBai), file(targets_bed), file(ref_fasta), file(ref_fai), file(ref_dict) from input_noChunk_ch

    // output:
    // set val("${label}"), val(comparisonID), val(tumorID), val(normalID), file("${output_norm_vcf}") into variants_noChunk
    // file("${output_vcf}")
    // file("${multiallelics_stats}")
    // file("${realign_stats}")
    // when: params.disable != "true"

    script:
    label = "noChunk"
    prefix = "${comparisonID}.${label}"
    config_file = "pindel_config.txt"
    output_dir = "${prefix}.pindel_output"
    output_vcf = "${prefix}.vcf"
    """
    # make config file for Pindel
    printf "${tumorBam}\t${params.insert_size}\tTUMOR\n" > "${config_file}"
    printf "${normalBam}\t${params.insert_size}\tNORMAL\n" >> "${config_file}"

    mkdir "${output_dir}"

    pindel \
    --fasta "${ref_fasta}" \
    --config-file "${config_file}" \
    --output-prefix "${output_dir}/" \
    --number_of_threads \${NSLOTS:-\${NTHREADS:-1}} \
    --include "${targets_bed}"

    pindel2vcf \
    --pindel_output_root "${output_dir}/" \
    --reference "${ref_fasta}" \
    --reference_name hg19 \
    --reference_date 2012_03_15 \
    --gatk_compatible \
    --vcf "${output_vcf}"
    """
}
//
// process filter_vcf_noChunk {
//     // filter the .vcf to only include 'PASS' entries
//     tag "${prefix}"
//     publishDir "${params.outputDir}/variants", mode: 'copy', overwrite: true
//
//     input:
//     set val(label), val(comparisonID), val(tumorID), val(normalID), file(vcf) from variants_noChunk
//
//     output:
//     set val(label), val(comparisonID), val(tumorID), val(normalID), file("${output_vcf}") into variants_noChunk_filtered
//
//     script:
//     prefix = "${comparisonID}.${label}"
//     output_vcf = "${prefix}.filter.vcf"
//     """
//     # get the header
//     grep '^#' "${vcf}" > "${output_vcf}"
//     # get the 'PASS' entries
//     grep -v '^#' "${vcf}" | grep 'PASS' >> "${output_vcf}" || :
//     """
// }
//
// process vcf_to_tsv_noChunk {
//     tag "${prefix}"
//     publishDir "${params.outputDir}/variants", overwrite: true, mode: 'copy'
//
//     input:
//     set val(label), val(comparisonID), val(tumorID), val(normalID), file(vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from variants_noChunk_filtered.combine(ref_fasta6).combine(ref_fai6).combine(ref_dict6)
//
//     output:
//     set val(label), val(comparisonID), val(tumorID), val(normalID), file(vcf), file("${reformat_tsv}") into vcf_tsv_noChunk // to annotation
//
//     script:
//     caller = "MuTect2"
//     prefix = "${comparisonID}.${label}"
//     tsv_file = "${prefix}.tsv"
//     reformat_tsv = "${prefix}.reformat.tsv"
//     """
//     # convert VCF to TSV
//     # NOTE: automatically filters for only PASS entries
//     gatk.sh -T VariantsToTable \
//     -R "${ref_fasta}" \
//     -V "${vcf}" \
//     -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN -F NLOD -F TLOD \
//     -GF AD -GF DP -GF AF \
//     -o "${tsv_file}"
//
//     # reformat and adjust the TSV table for consistency downstream
//     # add extra columns to the VCF TSV file for downstream
//     reformat-vcf-table.py -c MuTect2 -s "${tumorID}" -i "${tsv_file}" | \
//     paste-col.py --header "Sample" -v "${tumorID}"  | \
//     paste-col.py --header "Tumor" -v "${tumorID}"  | \
//     paste-col.py --header "Normal" -v "${normalID}"  | \
//     paste-col.py --header "VariantCaller" -v "${caller}" > \
//     "${reformat_tsv}"
//     """
// }
//
//
//
//
// ~~~~~~~ CHROM CHUNK ~~~~~~~~ //
targets_bed2.combine(chroms).set { targets_chroms }
process bed_chromChunk {
    // split the targets.bed by chromosome
    executor "local"
    tag "${chrom}"

    input:
    set file(targets_bed), val(chrom) from targets_chroms

    output:
    set val("${label}"), val(chrom), file("${output_bed}") into chromChunk_targets

    script:
    label = "chromChunk"
    output_bed = "targets.${chrom}.bed"
    """
    subset_bed.py "${chrom}" "${targets_bed}" > "${output_bed}"
    """
}

input_bams2.combine(chromChunk_targets).map { comparisonID, tumorID, normalID, tumorBam, tumorBai, normalBam, normalBai, label, chrom, targets_bed ->
    return([ label, chrom, comparisonID, tumorID, normalID, tumorBam, tumorBai, normalBam, normalBai, targets_bed ])
}
.combine(ref_fasta2)
.combine(ref_fai2)
.combine(ref_dict2)
.tap { input_chromChunk_ch }
// .subscribe { println "${it}" }

process pindel_chromChunk {
    tag "${prefix}"
    publishDir "${params.outputDir}/variants", overwrite: true, mode: 'copy'

    input:
    set val(label), val(chrom), val(comparisonID), val(tumorID), val(normalID), file(tumorBam), file(tumorBai), file(normalBam), file(normalBai), file("targets.bed"), file(ref_fasta), file(ref_fai), file(ref_dict) from input_chromChunk_ch

    // output:
    // set val("${label}"), val(chrom), val(comparisonID), val(tumorID), val(normalID), file("${output_norm_vcf}") into variants_chromChunk
    // file("${output_vcf}")
    // file("${multiallelics_stats}")
    // file("${realign_stats}")
    //
    // when: params.disable != "true"

    script:
    prefix = "${comparisonID}.${label}.${chrom}"
    config_file = "pindel_config.txt"
    output_dir = "${prefix}.pindel_output"
    output_vcf = "${prefix}.vcf"
    """
    # make config file for Pindel
    printf "${tumorBam}\t${params.insert_size}\tTUMOR\n" > "${config_file}"
    printf "${normalBam}\t${params.insert_size}\tNORMAL\n" >> "${config_file}"

    mkdir "${output_dir}"

    pindel \
    --fasta "${ref_fasta}" \
    --config-file "${config_file}" \
    --output-prefix "${output_dir}/" \
    --number_of_threads \${NSLOTS:-\${NTHREADS:-1}} \
    --include "targets.bed"

    pindel2vcf \
    --pindel_output_root "${output_dir}/" \
    --reference "${ref_fasta}" \
    --reference_name hg19 \
    --reference_date 2012_03_15 \
    --gatk_compatible \
    --vcf "${output_vcf}"
    """
}

// process filter_vcf_chromChunk {
//     // filter the .vcf to only include 'PASS' entries
//     tag "${prefix}"
//     publishDir "${params.outputDir}/variants", mode: 'copy', overwrite: true
//
//     input:
//     set val(label), val(chrom), val(comparisonID), val(tumorID), val(normalID), file(vcf) from variants_chromChunk
//
//     output:
//     set val(label), val(chrom), val(comparisonID), val(tumorID), val(normalID), file("${output_vcf}") into variants_chromChunk_filtered
//
//     script:
//     prefix = "${comparisonID}.${label}.${chrom}"
//     output_vcf = "${prefix}.filter.vcf"
//     """
//     # get the header
//     grep '^#' "${vcf}" > "${output_vcf}"
//     # get the 'PASS' entries
//     grep -v '^#' "${vcf}" | grep 'PASS' >> "${output_vcf}" || :
//     """
// }
//
// process vcf_to_tsv_chromChunk {
//     tag "${prefix}"
//     publishDir "${params.outputDir}/variants", overwrite: true, mode: 'copy'
//
//     input:
//     set val(label), val(chrom), val(comparisonID), val(tumorID), val(normalID), file(vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from variants_chromChunk_filtered.combine(ref_fasta5).combine(ref_fai5).combine(ref_dict5)
//
//     output:
//     set val(label), val(chrom), val(comparisonID), val(tumorID), val(normalID), file(vcf), file("${reformat_tsv}") into vcf_tsv_chromChunk // to annotation
//
//     script:
//     caller = "MuTect2"
//     prefix = "${comparisonID}.${label}.${chrom}"
//     tsv_file = "${prefix}.tsv"
//     reformat_tsv = "${prefix}.reformat.tsv"
//     """
//     # convert VCF to TSV
//     # NOTE: automatically filters for only PASS entries
//     gatk.sh -T VariantsToTable \
//     -R "${ref_fasta}" \
//     -V "${vcf}" \
//     -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN -F NLOD -F TLOD \
//     -GF AD -GF DP -GF AF \
//     -o "${tsv_file}"
//
//     # reformat and adjust the TSV table for consistency downstream
//     # add extra columns to the VCF TSV file for downstream
//     reformat-vcf-table.py -c MuTect2 -s "${tumorID}" -i "${tsv_file}" | \
//     paste-col.py --header "Sample" -v "${tumorID}"  | \
//     paste-col.py --header "TUMOR" -v "${tumorID}"  | \
//     paste-col.py --header "Normal" -v "${normalID}"  | \
//     paste-col.py --header "VariantCaller" -v "${caller}" > \
//     "${reformat_tsv}"
//     """
// }
//
//
//
//
// ~~~~~~~ N CHUNK ~~~~~~~~ //
nChunkLabel_ch.combine(targets_bed3).set { nChunk_target }//.subscribe { println "${it}" }
process bed_nChunk {
    executor "local"

    input:
    set val(nChunkLabel), val(nchunks), file("targets") from nChunk_target

    output:
    file("targets.*") into chunked_targets

    script:
    """
    chunk-lines.py "targets" "${nchunks}"
    """
}
chunked_targets.flatten()
    .combine(nChunkLabel_ch2)
    .map { targets, chunkLabel, numChunks ->
        def targetChunkNum = "${targets.name}".findAll(/\d*$/)[0] // number at the end of the file basename
        return([ chunkLabel, targets, targetChunkNum ])
    }
    .combine(input_bams3)
    .map { chunkLabel, targets, targetChunkNum, comparisonID, tumorID, normalID, tumorBam, tumorBai, normalBam, normalBai ->
        return([ chunkLabel, targetChunkNum, comparisonID, tumorID, normalID, tumorBam, tumorBai, normalBam, normalBai, targets])
    }
    .combine(ref_fasta3)
    .combine(ref_fai3)
    .combine(ref_dict3)
    .tap { input_nChunk_ch }
    // .subscribe { println "${it}" }

process pindel_nChunk {
    tag "${prefix}"
    publishDir "${params.outputDir}/variants", overwrite: true, mode: 'copy'

    input:
    set val(chunkLabel), val(targetChunkNum), val(comparisonID), val(tumorID), val(normalID), file(tumorBam), file(tumorBai), file(normalBam), file(normalBai), file("targets.bed"), file(ref_fasta), file(ref_fai), file(ref_dict) from input_nChunk_ch

    // output:
    // set val(chunkLabel), val(targetChunkNum), val(comparisonID), val(tumorID), val(normalID), file("${output_norm_vcf}") into variants_nChunk
    // file("${output_vcf}")
    // file("${multiallelics_stats}")
    // file("${realign_stats}")
    //
    // when: params.disable != "true"

    script:
    prefix = "${comparisonID}.${chunkLabel}.${targetChunkNum}"
    config_file = "pindel_config.txt"
    output_dir = "${prefix}.pindel_output"
    output_vcf = "${prefix}.vcf"
    """
    # make config file for Pindel
    printf "${tumorBam}\t${params.insert_size}\tTUMOR\n" > "${config_file}"
    printf "${normalBam}\t${params.insert_size}\tNORMAL\n" >> "${config_file}"

    mkdir "${output_dir}"

    pindel \
    --fasta "${ref_fasta}" \
    --config-file "${config_file}" \
    --output-prefix "${output_dir}/" \
    --number_of_threads \${NSLOTS:-\${NTHREADS:-1}} \
    --include "targets.bed"

    pindel2vcf \
    --pindel_output_root "${output_dir}/" \
    --reference "${ref_fasta}" \
    --reference_name hg19 \
    --reference_date 2012_03_15 \
    --gatk_compatible \
    --vcf "${output_vcf}"
    """
}

//
// process filter_vcf_nChunk {
//     // filter the .vcf to only include 'PASS' entries
//     tag "${prefix}"
//     publishDir "${params.outputDir}/variants", mode: 'copy', overwrite: true
//
//     input:
//     set val(chunkLabel), val(targetChunkNum), val(comparisonID), val(tumorID), val(normalID), file(vcf) from variants_nChunk
//
//     output:
//     set val(chunkLabel), val(targetChunkNum), val(comparisonID), val(tumorID), val(normalID), file("${output_vcf}") into variants_nChunk_filtered
//
//     script:
//     prefix = "${comparisonID}.${chunkLabel}.${targetChunkNum}"
//     output_vcf = "${prefix}.filter.vcf"
//     """
//     # get the header
//     grep '^#' "${vcf}" > "${output_vcf}"
//     # get the 'PASS' entries
//     grep -v '^#' "${vcf}" | grep 'PASS' >> "${output_vcf}" || :
//     """
// }
//
//
// process vcf_to_tsv_nChunk {
//     tag "${prefix}"
//     publishDir "${params.outputDir}/variants", overwrite: true, mode: 'copy'
//
//     input:
//     set val(chunkLabel), val(targetChunkNum), val(comparisonID), val(tumorID), val(normalID), file(vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from variants_nChunk_filtered.combine(ref_fasta4).combine(ref_fai4).combine(ref_dict4)
//
//     output:
//     set val(chunkLabel), val(targetChunkNum), val(comparisonID), val(tumorID), val(normalID), file(vcf), file("${reformat_tsv}") into vcf_tsv_nChunk // to annotation
//
//     script:
//     caller = "MuTect2"
//     prefix = "${comparisonID}.${chunkLabel}.${targetChunkNum}"
//     tsv_file = "${prefix}.tsv"
//     reformat_tsv = "${prefix}.reformat.tsv"
//     """
//     # convert VCF to TSV
//     # NOTE: automatically filters for only PASS entries
//     gatk.sh -T VariantsToTable \
//     -R "${ref_fasta}" \
//     -V "${vcf}" \
//     -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN -F NLOD -F TLOD \
//     -GF AD -GF DP -GF AF \
//     -o "${tsv_file}"
//
//     # reformat and adjust the TSV table for consistency downstream
//     # add extra columns to the VCF TSV file for downstream
//     reformat-vcf-table.py -c MuTect2 -s "${tumorID}" -i "${tsv_file}" | \
//     paste-col.py --header "Sample" -v "${tumorID}"  | \
//     paste-col.py --header "Tumor" -v "${tumorID}"  | \
//     paste-col.py --header "Normal" -v "${normalID}"  | \
//     paste-col.py --header "VariantCaller" -v "${caller}" > \
//     "${reformat_tsv}"
//     """
// }
//
//
//
//
//
//
// // ~~~~~~~ LINE CHUNK ~~~~~~~~ //
lineChunkLabel_ch.combine(targets_bed4).set { lineChunk_target }
process targetChunkLines {
    executor "local"

    input:
    set val(lineChunkLabel), val(lineChunks), file("targets") from lineChunk_target

    output:
    file("targets.*") into line_chunked_targets

    script:
    """
    chunk-by-lines.py "targets" "${lineChunks}"
    """
}
line_chunked_targets.flatten()
    .combine(lineChunkLabel_ch2)
    .map { targets, chunkLabel, numLineChunks ->
        def targetChunkNum = "${targets.name}".findAll(/\d*$/)[0] // number at the end of the file basename
        return([ chunkLabel, targets, targetChunkNum ])
    }
    .combine(input_bams4)
    .map { chunkLabel, targets, targetChunkNum, comparisonID, tumorID, normalID, tumorBam, tumorBai, normalBam, normalBai ->
        return([ chunkLabel, targetChunkNum, comparisonID, tumorID, normalID, tumorBam, tumorBai, normalBam, normalBai, targets])
    }
    .combine(ref_fasta7)
    .combine(ref_fai7)
    .combine(ref_dict7)
    .tap { input_lineChunk_ch }
    // .subscribe { println "[line_chunked_targets] ${it}" }

process pindel_lineChunk {
    tag "${prefix}"
    publishDir "${params.outputDir}/variants", overwrite: true, mode: 'copy'

    input:
    set val(chunkLabel), val(targetChunkNum), val(comparisonID), val(tumorID), val(normalID), file(tumorBam), file(tumorBai), file(normalBam), file(normalBai), file("targets.bed"), file(ref_fasta), file(ref_fai), file(ref_dict) from input_lineChunk_ch

    // output:
    // set val(chunkLabel), val(targetChunkNum), val(comparisonID), val(tumorID), val(normalID), file("${output_norm_vcf}") into variants_lineChunk
    // file("${output_vcf}")
    // file("${multiallelics_stats}")
    // file("${realign_stats}")
    //
    // when: params.disable != "true"

    script:
    prefix = "${comparisonID}.${chunkLabel}.${targetChunkNum}"
    config_file = "pindel_config.txt"
    output_dir = "${prefix}.pindel_output"
    output_vcf = "${prefix}.vcf"
    """
    # make config file for Pindel
    printf "${tumorBam}\t${params.insert_size}\tTUMOR\n" > "${config_file}"
    printf "${normalBam}\t${params.insert_size}\tNORMAL\n" >> "${config_file}"

    mkdir "${output_dir}"

    pindel \
    --fasta "${ref_fasta}" \
    --config-file "${config_file}" \
    --output-prefix "${output_dir}/" \
    --number_of_threads \${NSLOTS:-\${NTHREADS:-1}} \
    --include "targets.bed"

    pindel2vcf \
    --pindel_output_root "${output_dir}/" \
    --reference "${ref_fasta}" \
    --reference_name hg19 \
    --reference_date 2012_03_15 \
    --gatk_compatible \
    --vcf "${output_vcf}"
    """
}

// process filter_vcf_lineChunk {
//     // filter the .vcf to only include 'PASS' entries
//     tag "${prefix}"
//     publishDir "${params.outputDir}/variants", mode: 'copy', overwrite: true
//
//     input:
//     set val(chunkLabel), val(targetChunkNum), val(comparisonID), val(tumorID), val(normalID), file(vcf) from variants_lineChunk
//
//     output:
//     set val(chunkLabel), val(targetChunkNum), val(comparisonID), val(tumorID), val(normalID), file("${output_vcf}") into variants_lineChunk_filtered
//
//     script:
//     prefix = "${comparisonID}.${chunkLabel}.${targetChunkNum}"
//     output_vcf = "${prefix}.filter.vcf"
//     """
//     # get the header
//     grep '^#' "${vcf}" > "${output_vcf}"
//     # get the 'PASS' entries
//     grep -v '^#' "${vcf}" | grep 'PASS' >> "${output_vcf}" || :
//     """
// }
//
//
// process vcf_to_tsv_lineChunk {
//     tag "${prefix}"
//     publishDir "${params.outputDir}/variants", overwrite: true, mode: 'copy'
//
//     input:
//     set val(chunkLabel), val(targetChunkNum), val(comparisonID), val(tumorID), val(normalID), file(vcf), file(ref_fasta), file(ref_fai), file(ref_dict) from variants_lineChunk_filtered.combine(ref_fasta8).combine(ref_fai8).combine(ref_dict8)
//
//     output:
//     set val(chunkLabel), val(targetChunkNum), val(comparisonID), val(tumorID), val(normalID), file(vcf), file("${reformat_tsv}") into vcf_tsv_lineChunk // to annotation
//
//     script:
//     caller = "MuTect2"
//     prefix = "${comparisonID}.${chunkLabel}.${targetChunkNum}"
//     tsv_file = "${prefix}.tsv"
//     reformat_tsv = "${prefix}.reformat.tsv"
//     """
//     # convert VCF to TSV
//     # NOTE: automatically filters for only PASS entries
//     gatk.sh -T VariantsToTable \
//     -R "${ref_fasta}" \
//     -V "${vcf}" \
//     -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F QUAL -F AC -F AN -F NLOD -F TLOD \
//     -GF AD -GF DP -GF AF \
//     -o "${tsv_file}"
//
//     # reformat and adjust the TSV table for consistency downstream
//     # add extra columns to the VCF TSV file for downstream
//     reformat-vcf-table.py -c MuTect2 -s "${tumorID}" -i "${tsv_file}" | \
//     paste-col.py --header "Sample" -v "${tumorID}"  | \
//     paste-col.py --header "Tumor" -v "${tumorID}"  | \
//     paste-col.py --header "Normal" -v "${normalID}"  | \
//     paste-col.py --header "VariantCaller" -v "${caller}" > \
//     "${reformat_tsv}"
//     """
// }
//
//
// // ~~~~~~~ ANNOTATE ~~~~~~~~ //
// // refactor the channels for consistency
// vcf_tsv_lineChunk.map { chunkLabel, targetChunkNum, comparisonID, tumorID, normalID, vcf, tsv ->
//     def label = "${comparisonID}.${chunkLabel}.${targetChunkNum}"
//     def lineChunkLabel = "${params.numlines}lineChunk"
//     return([ lineChunkLabel, label, vcf, tsv ])
// }.set { vcf_tsv_lineChunk_refactor }
//
// vcf_tsv_nChunk.map { chunkLabel, targetChunkNum, comparisonID, tumorID, normalID, vcf, tsv ->
//     def label = "${comparisonID}.${chunkLabel}.${targetChunkNum}"
//     def nChunkLabel = "${params.numChunks}Chunk"
//     return([ nChunkLabel, label, vcf, tsv ])
// }.set { vcf_tsv_nChunk_refactor }
//
// vcf_tsv_chromChunk.map { label, chrom, comparisonID, tumorID, normalID, vcf, tsv ->
//     def newlabel = "${comparisonID}.${label}.${chrom}"
//     def type = "chromChunk"
//     return([ type, newlabel, vcf, tsv ])
// }.set { vcf_tsv_chromChunk_refactor }
//
// vcf_tsv_noChunk.map { label, comparisonID, tumorID, normalID, vcf, tsv ->
//     def newlabel = "${comparisonID}.${label}"
//     def type = "noChunk"
//     return([ type, newlabel, vcf, tsv ])
// }.set { vcf_tsv_noChunk_refactor }
//
// // make sure there are variants in the TSV
// import java.nio.file.Files;
// vcf_tsv_nChunk_refactor.concat(vcf_tsv_chromChunk_refactor, vcf_tsv_noChunk_refactor, vcf_tsv_lineChunk_refactor)
//         .filter { type, label, vcf, tsv ->
//             long count = Files.lines(tsv).count()
//             if (count <= 1) println ">>> WARNING: file ${tsv} (${label}) does not have enough lines and will not be included"
//             count > 1
//         }
//         .set { vcfs_tsvs }
// // vcfs_tsvs.subscribe { println "[vcfs_tsvs] ${it}" }
//
// process annotate {
//     // annotate the VCF file
//     tag "${prefix}"
//     publishDir "${params.outputDir}/annotations", overwrite: true, mode: 'copy'
//
//     input:
//     set val(type), val(label), file(vcf), file(tsv), file(annovar_db_dir) from vcfs_tsvs.combine(annovar_db_dir)
//
//     output:
//     set val(type), val(label), file(vcf), file(tsv), file("${annovar_output_txt}"), file("${avinput_tsv}") into vcfs_tsvs_annotations
//
//     script:
//     prefix = "${label}"
//     avinput_file = "${prefix}.avinput"
//     avinput_tsv = "${avinput_file}.tsv"
//     annovar_output_txt = "${prefix}.${params.ANNOVAR_BUILD_VERSION}_multianno.txt"
//     """
//     table_annovar.pl "${vcf}" "${annovar_db_dir}" \
//     --buildver "${params.ANNOVAR_BUILD_VERSION}" \
//     --remove \
//     --protocol "${params.ANNOVAR_PROTOCOL}" \
//     --operation "${params.ANNOVAR_OPERATION}" \
//     --nastring . \
//     --vcfinput \
//     --onetranscript \
//     --outfile "${prefix}"
//     printf "Chr\tStart\tEnd\tRef\tAlt\tCHROM\tPOS\tID\tREF\tALT\n" > "${prefix}.avinput.tsv"
//     cut -f1-5,9-13 ${avinput_file} >>  "${avinput_tsv}"
//     """
// }
//
//
// process merge_tables {
//     // merge the annotation and vcf tables
//     tag "${prefix}"
//     // publishDir "${params.outputDir}", mode: 'copy', overwrite: true
//
//     input:
//     set val(type), val(label), file(vcf), file(tsv), file(annovar_txt), file(avinput_tsv) from vcfs_tsvs_annotations
//
//     output:
//     file("${output_annotations}") into merged_tables
//
//     script:
//     prefix = "${label}"
//     output_annotations = "${prefix}.annotations.tsv"
//     """
//     merge-vcf-tables.R "${tsv}" "${annovar_txt}" "${avinput_tsv}" tmp.tsv
//
//     cat tmp.tsv | \
//     paste-col.py --header "Type" -v "${type}"  | \
//     paste-col.py --header "Label" -v "${label}"  > \
//     "${output_annotations}"
//     """
//
// }
// // concatenate all the variant tables
// merged_tables.collectFile(name: "MuTect2.chunking_validation.annotations.tsv", storeDir: "${params.outputDir}", keepHeader: true, sort: false)
