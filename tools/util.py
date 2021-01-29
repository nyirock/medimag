def remove_suffix(text, suffix=""):
    if suffix and text.endswith(suffix):
        return text[:-len(suffix)]
    return text

def remove_common_substring(text, common_suffixes = [".gz", ".gzip", ".fastq", ".fasta", ".fq", ".fa", ".fna"]):

    for suffix in common_suffixes:
        text = text.replace(suffix, '')

    return text