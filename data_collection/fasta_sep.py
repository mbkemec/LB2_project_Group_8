def seperation(fasta_file, output_dir):  

    with open(fasta_file, "r") as f:
        seq_id = None
        seq = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):  # new sequence
                if seq_id is not None:
                    # Save previous sequence
                    with open(f"{output_dir}/{seq_id}.fasta", "w") as out:
                        out.write(f">{seq_id}\n{''.join(seq)}\n")

                # get new sequence ID
                seq_id = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)

        # Save the last sequence
        if seq_id is not None:
            with open(f"{output_dir}/{seq_id}.fasta", "w") as out:
                out.write(f">{seq_id}\n{''.join(seq)}\n")

if __name__=='__main__':    

    fasta_file = input('give me the name of the file you want to seperate one by one: ')
    output_dir = input('give me the name of the output direction: ')
    seperation(fasta_file, output_dir)