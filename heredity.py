import csv
import itertools
import sys
import copy

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    # Create dict to represent network of people and traits
    if "target_network" in locals():
        target_network.clear()
    target_network = copy.deepcopy(people)

    # Add target number of genes to dict
    for person in target_network:
        if person in one_gene:
            target_network[person]["num_genes"] = 1
        elif person in two_genes:
            target_network[person]["num_genes"] = 2
        else:
            target_network[person]["num_genes"] = 0
    
    # Calculate probability that each person has target number of genes
    for person in target_network:
        calc_num_genes_probability(target_network, target_network[person]["name"])

    # Add trait values to dict
    for person in target_network:
        if person in have_trait:
            target_network[person]["trait"] = True
        else:
            target_network[person]["trait"] = False
    
    # Calculate conditional probability of each person's trait value given their number of genes
    for person in target_network:
        calc_conditional_trait_probability(target_network, person)
    
    # Calculate joint probability
    joint_prob = 1
    for person in target_network:
        joint_prob *= target_network[person]["p_trait"]

    return joint_prob
        

def calc_num_genes_probability(target_network, target_person):
    mother = target_network[target_person]["mother"]
    father = target_network[target_person]["father"]

    # Recursively call function to calculate parent probabilities
    if (mother != None) or (father != None):
        if target_network[mother] != None:
            calc_num_genes_probability(target_network, mother)
        if target_network[father] != None:
            calc_num_genes_probability(target_network, father)

    inheritance_probs = {
        # Probability that parent passes 0 genes based on parent's number of genes 
        0: {
            0: 1 - PROBS["mutation"],
            1: 0.5,
            2: PROBS["mutation"]
        },
        # Probability that parent passes 1 gene based on parent's number of genes
        1: {
            0: PROBS["mutation"],
            1: 0.5,
            2: 1 - PROBS["mutation"]
        }
    }

    target_person_genes = target_network[target_person]["num_genes"]

    # If parents not in database, assume general population probability
    if (mother == None) and (father == None):
        target_network[target_person]["p_genes"] = PROBS["gene"][target_person_genes]
    
    # If parents in database, calculate probability based on parent genes
    else:
        # If target_person_genes = 0, probability = p(0 from mother and 0 from father)
        if target_person_genes == 0:
            p_zero_from_mother = inheritance_probs[0][target_network[mother]["num_genes"]]
            p_zero_from_father = inheritance_probs[0][target_network[father]["num_genes"]]

            target_network[target_person]["p_genes"] = p_zero_from_mother * p_zero_from_father
        # If target_person_genes = 1, probability = p(1 from mother and 0 from father) or p(0 from mother and 1 from father)
        elif target_person_genes == 1:
            # 1 from mother, 0 from father
            p_one_from_mother = inheritance_probs[1][target_network[mother]["num_genes"]]
            p_zero_from_father = inheritance_probs[0][target_network[father]["num_genes"]]

            # 0 from mother, 1 from father
            p_one_from_father = inheritance_probs[1][target_network[father]["num_genes"]]
            p_zero_from_mother = inheritance_probs[0][target_network[mother]["num_genes"]]

            target_network[target_person]["p_genes"] = (p_one_from_mother * p_zero_from_father) + (p_zero_from_mother * p_one_from_father)
        # If target_person_genes = 2, probability = p(1 from mother and 1 from father)
        else:
            p_one_from_mother = inheritance_probs[1][target_network[mother]["num_genes"]]
            p_one_from_father = inheritance_probs[1][target_network[father]["num_genes"]]

            target_network[target_person]["p_genes"] = p_one_from_mother * p_one_from_father


def calc_conditional_trait_probability(target_network, target_person):
    num_genes = target_network[target_person]["num_genes"]
    have_trait = target_network[target_person]["trait"]

    # Lookup trait probability in `PROBS` dict and multiply by genes probability
    target_network[target_person]["p_trait"] = target_network[target_person]["p_genes"] * PROBS["trait"][num_genes][have_trait]


def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    for person in probabilities:
        # Update gene probabilities
        if person in one_gene:
            probabilities[person]["gene"][1] += p
        elif person in two_genes:
            probabilities[person]["gene"][2] += p
        else:
            probabilities[person]["gene"][0] += p
        # Update trait probabilities
        if person in have_trait:
            probabilities[person]["trait"][True] += p
        else:
            probabilities[person]["trait"][False] += p


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    # Normalize probability distributions by multiplying them by the inverse of their sum
    for person in probabilities:
        gene_normalization_factor = 0
        trait_normalization_factor = 0

        # Calculate sum of distributions
        for prob in probabilities[person]["gene"].values():
            gene_normalization_factor += prob
        for prob in probabilities[person]["trait"].values():
            trait_normalization_factor += prob

        # Take inverse
        gene_normalization_factor = 1 / gene_normalization_factor
        trait_normalization_factor = 1 / trait_normalization_factor

        # Multiply each probability by the appropriate normalizing factor
        for prob in probabilities[person]["gene"]:
            probabilities[person]["gene"][prob] *= gene_normalization_factor
        for prob in probabilities[person]["trait"]:
            probabilities[person]["trait"][prob] *= trait_normalization_factor
        

if __name__ == "__main__":
    main()
