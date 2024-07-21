## Python 3
# Filter matching topologies
# Created by:  Julija Pecerska <pece@zhaw.ch>
# ACGT ZHAW
# Created date: July 2024

from pathlib import Path
from ete3 import Tree

if __name__ == '__main__':
    path = Path("./new_2024.07.04.01_Prank_phyml_redo")
    all_tree_files = list(path.glob('*.phy_phyml_tree.txt'))
    total_datasets = len(all_tree_files)
    print("Number of datasets: {}".format(total_datasets))
    true_topo = Tree("((Mus:0.0, Rattus:0.0):0.0,(((Homo:0.0, Pan:0.0), Gorilla:0.0), Macaca:0.0));")
    true_topo.unroot()
    print("True topology:{}\n".format(true_topo.write(format=9)))
    true_topo_names = []
    topologies = {}
    for file in all_tree_files:
        tree = Tree(str(file))
        new = True
        for topo in topologies:
            diff = topo.compare(tree, unrooted=True)
            if diff["rf"] == 0:
                new = False
                break
        if new:
            topologies[tree] = 1
        else:
            topologies[topo] += 1
        if tree.compare(true_topo, unrooted=True)["rf"] == 0:
            true_topo_names.append(file.name.split('.')[0].split('_')[1])
    print("Number of different topologies: {}\n".format(len(topologies)))

    print("Inferred more than 1% of the time:")
    for topology, topology_count in sorted(topologies.items(), key=lambda x:x[1]):
        if topology_count > 0.01 * total_datasets:
            print("Inferred {} times:".format(topology_count))
            print("{}\n".format(topology.write(format=9)))

    print("Number of topologies matching the species: {}\n".format(len(true_topo_names)))
    true_topo_names.sort()
    with open("agreeing_topology_IDs.txt", "w") as output:
        output.write("\n".join(true_topo_names))

        
    

