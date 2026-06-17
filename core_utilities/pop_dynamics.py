"""
Age-Structured Population Dynamics Simulator

A simplified, discrete-time population simulator for species with defined
maturation periods, fixed fecundity, and strict maximum lifespans.
While derived from Leslie Matrix principles, this engine assumes uniform
survival/fecundity across all mature cohorts (ideal for insects, rodents,
or stationary-phase microbes via the Carrying Capacity bottleneck).
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.0.0"

import sys
import argparse
from collections import deque
from utils import base_parser


def simulate_population(
    generations: int,
    max_age: int,
    fecundity: int = 1,
    carrying_capacity: int | None = None,
    initial_pop: int = 1,
) -> int:
    """
    Calculates the total population using a sliding age window (O(n) time).

    Args:
        generations: Total duration of the simulation in time steps.
        max_age: The maximum lifespan of the organism in time steps.
        fecundity: Number of offspring produced by mature units per step.
        carrying_capacity: Optional environmental limit (K) for logistic growth.

    Returns:
        The total population at the end of the simulation.

    Raises:
        ValueError: If generations < 1, max_age < 1, or fecundity < 1.
    """

    if generations < 1:
        raise ValueError("Simulation must run for at least 1 generation.")
    if max_age < 1:
        raise ValueError("Maximum age (lifespan) must be at least 1 time step.")
    if fecundity < 1:
        raise ValueError("Fecundity (offspring rate) cannot be less than 1.")

    population_cohorts = deque([0] * max_age)
    population_cohorts[0] = initial_pop
    total_population = initial_pop

    for _ in range(generations - 1):
        mature_units = total_population - population_cohorts[0]
        theoretical_newborns = mature_units * fecundity

        if carrying_capacity is not None:
            newborns = min(
                theoretical_newborns, max(0, carrying_capacity - total_population)
            )
        else:
            newborns = theoretical_newborns

        dying_units = population_cohorts.pop()
        population_cohorts.appendleft(newborns)

        total_population = total_population - dying_units + newborns

    return total_population


def setup_cli() -> argparse.ArgumentParser:
    """Configures the command-line interface arguments."""
    parser = base_parser(
        "Population Dynamics Simulator", include_input=False, include_output=True
    )

    parser.add_argument(
        "-i",
        "--inoculum",
        type=int,
        default=1,
        dest="initial_pop",
        metavar="Initial population or inoculum",
        help="The starting population size (N_0). Default is 1.",
    )

    parser.add_argument(
        "-t",
        "-n",
        "--generations",
        type=int,
        required=True,
        dest="generations",
        metavar="Generations or Time Steps",
        help="Total duration of the simulation (generation/time steps).",
    )
    parser.add_argument(
        "-m",
        "--max-age",
        type=int,
        required=True,
        dest="max_age",
        metavar="Viability limit or lifespan",
        help="The viability limit or lifespan of the organism (in generations/time steps).",
    )
    parser.add_argument(
        "-f",
        "--fecundity",
        type=int,
        default=1,
        dest="fecundity",
        metavar="RATE",
        help="Offspring produced per mature unit (default: 1. Use 2 for binary fision).",
    )
    parser.add_argument(
        "-c",
        "--capacity",
        type=int,
        default=None,
        dest="carrying_capacity",
        metavar="Carrying Capacity",
        help="Maximum environmental carrying capacity (K) for logistics growth.",
    )

    return parser


def main() -> None:
    "Pipelin manager handling CLI routing and error suppression."

    args = setup_cli().parse_args()

    try:
        final_pop = simulate_population(
            generations=args.generations,
            max_age=args.max_age,
            fecundity=args.fecundity,
            carrying_capacity=args.carrying_capacity,
            initial_pop=args.initial_pop,
        )

    except ValueError as e:
        sys.exit(f"\nSimulation Halted: {e}")
    except KeyboardInterrupt:
        sys.exit(f"\nSimulation interrupted by user.")

    print("\n" + "=" * 45)
    print("      POPULATION DYNAMICS REPORT")
    print("=" * 45)
    print(f"Generations / Time Steps: {args.generations}")
    print(f"Max Age / Lifespan:       {args.max_age}")
    print(f"Fecundity Rate:           {args.fecundity}")
    if args.carrying_capacity:
        print(f"Carrying Capacity (K):    {args.carrying_capacity}")
    print("-" * 45)
    print(f"Final Population:         {final_pop:,}")
    print("=" * 45 + "\n")

    if args.output:
        try:
            args.output.write_text(f"{final_pop}\n", encoding="utf-8")
            print(f"Result securely logged to: {args.output.name}\n")
        except OSError as e:
            sys.exit(f"Error writing to {args.output.name} as {e}")


if __name__ == "__main__":
    main()
