import sys
import ete3
import os
os.environ['QT_QPA_PLATFORM']='offscreen'

tree = sys.argv[1]
outputdir = sys.argv[2]

# Load a tree structure from a newick file.
t = ete3.Tree(tree, format=1)

t.unroot()

# Traverse tree and name internal nodes
node_counter = 1
for node in t.traverse("preorder"):
    if not node.is_leaf():
        node.name = "Node" + str(node_counter)
        node_counter += 1

# Render the original tree
os.makedirs(os.path.join(outputdir, "tree_vis"), exist_ok=True)
t.render(
    os.path.join(outputdir, "tree_vis", f"fulltree.png"),
    w=1500,
)

def get_splits(tree):
    # List to hold the resulting splits
    splits = []

    # Iterate over all the tree edges
    for node in tree.traverse():
        if not node.is_root():
            # Detach the node to create a split
            parent = node.up
            parent.remove_child(node)

            # Add the tree post removal and subtree to the list
            splits.append((tree.copy(method="newick"), node.copy(method="newick")))

            # Reattach the node to continue the iteration
            parent.add_child(node)

    return splits


splits = get_splits(t)
for i, (part1, part2) in enumerate(splits):
    part1.render(
        os.path.join(outputdir, "tree_vis", f"split{i}_part1.png"),
        w=1500,
    )
    part2.render(
        os.path.join(outputdir, "tree_vis", f"split{i}_part2.png"),
        w=1500,
    )
    with open(os.path.join(outputdir, f"split{i}_part1.txt"), "w") as file1, open(
        os.path.join(outputdir, f"split{i}_part2.txt"), "w"
    ) as file2:
        for l in part1.get_leaves():
            file1.write(l.name + "\n")
        for l in part2.get_leaves():
            file2.write(l.name + "\n")
