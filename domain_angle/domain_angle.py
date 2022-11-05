import re
import time
import os
import argparse
import pickle

import pymol
from pymol import cmd

from autogromacs.Analysis import mdplots, mdanalysis

import orientation

# FIXME: This group selection must be loaded from a project-wide
# configuration file in the future.
groups = ["Natural", "Artificial", "Derivada"]


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--visual-mode",
                        action="store_true")
    parser.add_argument("-d", "--base-directory",
                        help="Path to the directory containing .pdb structures.")
    parser.add_argument("-s", "--session-path",
                        help="Path to the autoGromacs analysis session file.")
    parser.add_argument("-n", "--run-name",
                        help="Path to the autoGromacs analysis session file.")

    return parser.parse_args()


def subchain_selector(cutoff, direction):

    w = f"{cutoff}-1000"
    if direction < 0:
        w = f"1-{cutoff}"

    return f"resi {w}"


def detect_angle_for_frame(domain_selectors, objects):

    A = objects[0]
    B = objects[1]

    # align on chain A
    cmd.align(f"{A} & {domain_selectors[0]}", f"{B} & {domain_selectors[0]}")

    # measure rotation and displacement of chain B
    w = orientation.angle_between_domains(
        f"{A} & {domain_selectors[1]}",
        f"{B} & {domain_selectors[1]}"
    )

    return w


def analyze(movie_filepath: str):
    # two conformations of a two-chain structure
    cmd.delete('all')
    cmd.load(movie_filepath)
    cmd.load("snapshot_T_400_SALT_ALL_total_all.pdb")

    k = sorted(cmd.get_object_list('all'))

    domain_cutoff = 129
    domain_selectors = [
        subchain_selector(domain_cutoff, -1),
        subchain_selector(domain_cutoff, 1)
    ]

    nbstates = cmd.count_states(selection="(all)")
    angles = []
    assert nbstates, "Invalid nbstates!"
    for frame in range(1, nbstates + 1):
        if not frame % 10:
            print(f"Frame {frame} of {nbstates}")
        cmd.frame(frame)
        angle = detect_angle_for_frame(domain_selectors, k)
        angles.append(angle)

    return angles


def extract_name_from_file(filepath):
    return re.findall("([^_]+).pdb$", filepath)[0]


def rank_structure_name(name):
    final_digit = re.findall(r"\d+", name)
    if final_digit:
        return int(final_digit[0])

    return 0


def analyze_angles(arguments, base_path):
    if arguments.visual_mode:
        pymol.finish_launching(['pymol', '-q'])

    all_angles = []
    all_labels = []

    pdb_structures = sorted(os.listdir(base_path), key=rank_structure_name)
    for f in pdb_structures:
        print(f)
        if not f.startswith("movie"):
            continue
        print("Checking file...")

        try:
            angles = analyze(os.path.join(base_path, f))
            all_angles.append(angles)
            all_labels.append(extract_name_from_file(f))
        except pymol.CmdException as e:
            print(e)

        if arguments.visual_mode:
            time.sleep(5)

    return {
        "angles": all_angles,
        "labels": all_labels
    }


def main(arguments):
    for directory in os.listdir(arguments.base_directory):
        base_path = os.path.join(arguments.base_directory, directory)
        print(f"Checking angles for group at {base_path}.")
        analyze_group(arguments, base_path, directory)


def analyze_group(arguments, base_path, identifier):
    angle_file = f"angle_datat_{identifier}.obj"

    if os.path.isfile(angle_file):
        with open(angle_file, 'rb') as f:
            angle_data = pickle.load(f)
    else:
        angle_data = analyze_angles(arguments, base_path)
        with open(angle_file, 'wb') as f:
            pickle.dump(angle_data, f)

    dt = angle_data["angles"]
    assert dt, "Empty angle data!"
    print(dt)
    print(len(angle_data["angles"]))

    # Load session <optional>.
    if arguments.session_path:
        sessions = mdanalysis.load_session(arguments.session_path)
        for session in sessions:
            if "total" in session.plot_suffix:
                print(session.labels)
                session.select_simulation_indexes([16, 17, 18, 19])
                rmsd_series = session.rmsd_series

        print(":>")
        print(len(rmsd_series))

        #A = [x * 1.5 for x in angle_data["angles"]]
        dt = []
        for t, x in zip(rmsd_series, angle_data["angles"]):
            dt.append([t, x])

    run_length = determine_run_length(identifier)
    # Plot angle series.
    mdplots.show_rms_series_stacked(
        dt,
        angle_data["labels"],
        [run_length for _ in angle_data["labels"]],
        f"ts_angle_stacked_total_{arguments.run_name}_{identifier}",
        "ANGLES"
    )

    mdplots.show_rms_series_monolithic(
        dt,
        angle_data["labels"],
        [run_length for _ in angle_data["labels"]],
        f"ts_angle_mono_total_{arguments.run_name}_{identifier}",
        "ANGLES"
    )


def determine_run_length(identifier):
    """Very weak, ad hoc method to extract simulation length from the simulation identifier."""
    if "4K" in identifier:
        return 4000
    return 400


if __name__ == "__main__":
    arguments = parse_arguments()
    main(arguments)
