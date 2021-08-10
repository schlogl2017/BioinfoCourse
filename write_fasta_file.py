#!usr/bin/en python


def write_fasta(data, file_name_to_save, dir_name, sud_dir=False, ssub=False):
    to_save = full_name(dir_name, sud_dir, ssub)
    if os.path.exists(to_save):
        pass
    else:
        make_me_a_folder(to_save)
    with gzip.open(f"{to_save}/{file_name_to_save}.fna.gz", "wt") as f:
        f.write(f">{name}_genome_shuffled\n")
        f.write(textwrap.fill(data, 80))
