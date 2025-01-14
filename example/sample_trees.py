import tskit
import numpy as np


def sample_trees(tree_sequence, n_trees=100, random_seed=None):
    """
    Randomly sample trees from a tree sequence and create a new tree sequence.

    Parameters:
    -----------
    tree_sequence : tskit.TreeSequence
        The input tree sequence to sample from
    n_trees : int, optional
        Number of trees to sample (default: 100)
    random_seed : int, optional
        Random seed for reproducibility

    Returns:
    --------
    tskit.TreeSequence
        A new tree sequence containing only the sampled trees
    """
    # Set random seed if provided
    if random_seed is not None:
        np.random.seed(random_seed)

    # Get total number of trees
    total_trees = tree_sequence.num_trees

    # Ensure we don't try to sample more trees than available
    n_trees = min(n_trees, total_trees)

    # Randomly select tree indices
    selected_indices = np.random.choice(total_trees, size=n_trees, replace=False)
    selected_indices.sort()  # Sort to maintain order

    # Get the breakpoints for these trees
    breakpoints = []
    trees = tree_sequence.trees()

    for tree_idx, tree in enumerate(trees):
        if tree_idx in selected_indices:
            if not breakpoints:
                breakpoints.append(tree.interval.left)
            breakpoints.append(tree.interval.right)

    # Create tables for the new tree sequence
    tables = tree_sequence.dump_tables()

    # Keep only the edges that fall within our selected intervals
    tables.edges.clear()
    for edge in tree_sequence.edges():
        for i in range(len(breakpoints) - 1):
            left = breakpoints[i]
            right = breakpoints[i + 1]

            # If edge overlaps with this interval, add it
            if edge.left < right and edge.right > left:
                new_left = max(edge.left, left)
                new_right = min(edge.right, right)
                tables.edges.add_row(
                    left=new_left,
                    right=new_right,
                    parent=edge.parent,
                    child=edge.child
                )

    # Return new tree sequence
    return tables.tree_sequence()


trees = tskit.load("chr22.part-02.relate.trees")
sampled_trees = sample_trees(trees, n_trees=100, random_seed=42)
sampled_trees.dump("chr22.part-02.relate.100.trees")
