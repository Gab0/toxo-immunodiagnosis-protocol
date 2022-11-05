
import biotite.sequence as seq
import biotite.sequence.io.fasta as fasta
import biotite.sequence.graphics as graphics
import biotite.sequence.align as align
import matplotlib.pyplot as plt
import biotite.application.clustalo  as clustalo
from Bio import AlignIO, SeqIO
# Visualize the first 200 columns of the alignment
# Reorder alignments to reflect sequence distance


def plot_seq(filename, segments, output):
    # aln = fasta.get_alignment(filename)
    # seq = AlignIO.read(filename, 'clustal')

    fasta_file = fasta.FastaFile.read(filename)
    seqs = list(fasta.get_sequences(fasta_file).values())
    print(seqs)
    #
    dim = 10.0
    fig = plt.figure(figsize=(dim * 1.41, dim), dpi=600)

    app = clustalo.ClustalOmegaApp(seqs)
    app.start()
    app.join()

    labels = [
        k.replace("\\", "")
        for k in fasta.get_sequences(fasta_file).keys()
    ]

    alignment = app.get_alignment()

    print(alignment)
    if segments is None:
        segments = [(0, len(alignment))]

    for i, (f, t) in enumerate(segments):
        ax = fig.add_subplot(len(segments), 1, i + 1)

        graphics.plot_alignment_type_based(
            ax,
            alignment[f:t],
            labels=labels,
            show_numbers=False,
            color_scheme="clustalx",
            show_line_position=True,
        )

        ticks = ax.xaxis.get_ticklabels()
        xticks = ax.get_xticks()
        print(xticks)

        def proc_label(l):
            try:
                return str(round(l + f))
            except ValueError:
                return ""

        ax.set_xticklabels([proc_label(t) for t in xticks])
    fig.tight_layout()
    plt.savefig(output)


if __name__ == "__main__":
    plot_seq(
        "../Patentes/Abbott_p30/sequences.fasta",
        [(50, 90), (185, 235), (236, 285), (286, -1)],
        "seq-abbott-p30.png"
    )
    plot_seq("SRS29B_INTRON.fasta", None, "seq-srs29b-intron.png")
