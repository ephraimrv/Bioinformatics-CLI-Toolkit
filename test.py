from Bio.SeqUtils.ProtParam import ProteinAnalysis


def calculate_mature_core(full_protein: str) -> str:
    """
    Biochemical algorithm to find the Double-Glycine cut site and
    trim off highly charged C-terminal tails based on hydrophobicity.
    """
    if "GG" not in full_protein:
        return full_protein  # If no cut site, use the whole sequence

    # Split at the FIRST 'GG' and keep everything after it
    parts = full_protein.split("GG", 1)
    mature_peptide = parts[1]

    MIN_LENGTH = 25
    if len(mature_peptide) < MIN_LENGTH:
        return mature_peptide

    # Scan the sequence for structural boundaries
    for i in range(MIN_LENGTH, len(mature_peptide)):
        current_residue = mature_peptide[i]

        # Rule A: Proline Helix-Breaker
        if current_residue == "P":
            return mature_peptide[: i + 1]

        # Rule B: Kyte-Doolittle Hydrophobicity Drop
        if i + 5 <= len(mature_peptide):
            window = mature_peptide[i : i + 5]
            analyzer = ProteinAnalysis(window)

            avg_hydro = analyzer.gravy()

            # If the window becomes highly charged/hydrophilic, snip it!
            if avg_hydro < -0.5:
                return mature_peptide[:i]

    return mature_peptide


if __name__ == "__main__":
    print("=========================================")
    print("  Bacteriocin Core Auto-Trimmer CLI")
    print("=========================================")
    print("Paste a sequence and press Enter.")
    print("Type 'quit' or press Ctrl+C to exit.\n")

    while True:
        try:
            # 1. Ask the user for input and clean off any accidental invisible spaces
            user_input = input("Paste protein sequence: ").strip()

            # 2. Check if the user wants to close the program
            if user_input.lower() in ["quit", "exit", "q"]:
                print("Exiting tool. Goodbye!")
                break

            # 3. If they accidentally hit Enter on an empty line, just ask again
            if not user_input:
                continue

            # 4. Run the sequence through your algorithm!
            result = calculate_mature_core(user_input.upper())

            # 5. Print the formatted result back to the terminal
            print(f"\n  [+] Calculated Core : {result}")
            print(f"  [i] Core Length     : {len(result)} amino acids\n")
            print("-" * 50)

        except KeyboardInterrupt:
            # Gracefully handle the user pressing Ctrl+C to force-quit
            print("\nForce quitting tool. Goodbye!")
            break
