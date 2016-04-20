def fitch(node):
    """Infers the state of characters for all ancestral nodes in a given tree 
    by majority rule.  Also tries to determine single-step events that could 
    lead to the gain and loss of characters.
    This function expects a bifurcating tree as an ete3-object, in which the 
    leafs contain the feature 'characters' (a set of all characters present
    in the leaf).
    """
    leaves_to_root(node)
    root_to_leaves(node)
    ancestral_states(node)

def leaves_to_root(node):
    """For one node, collect all the characters that are present in the 
    children.  If a single character is present in both children, set the 
    state for this character to present (1).  If an character is only 
    present in one of the children, set the state to unknown (0).  Because we 
    are only comparing characters that are present in at least one of the 
    children, getting an character which is not present in both children is 
    impossible.  
    """
    if node.is_leaf():
        return dict((arr, set([1])) for arr in node.characters.keys())

    state = dict()
    lstates = list()
    for child in node.get_children():
        lstates.append(leaves_to_root(child))

    lcharacters = [set(child_state.keys()) for child_state in lstates]
    characters = set()
    if lcharacters != []:
        characters = set(lcharacters[0])
    for child_characters in lcharacters[1]:
        characters = characters.union(child_characters)

    for character in characters:
        # if state defined for character in every child, take intersection
        every_child = {False if character not in child_state else True for child_state in lstates}
        character_states = [child_state.get(character, set([0])) for child_state in lstates]
        inter_state = character_states[0]
        for char_state in character_states[1:]:
            inter_state = inter_state.intersection(char_state)
        if inter_state != set():
            state[character] = inter_state
        else:
            union_state = character_states[0]
            for char_state in character_states[1:]:
                union_state = union_state.union(char_state)
            state[character] = union_state

    node.add_features(characters=state)
    return node.characters

def root_to_leaves(node):
    """Determine the state for uncertain characters in a node according to 
    the state it has in the parents.  If a domain status is uncertain at the 
    root, set it to present (1).  In other nodes, determine the status of 
    uncertain characters by looking that particular character up in the
    parent node.  If it is present in the parent, set it to present (1),
    otherwise delete it from that nodes repository.
    """
    if node.is_leaf():
        return
    elif node.is_root():
        # delete all characters which are definitely not there and set the
        # status for unsure ones to present
        node.characters = set(a for (a, c) in node.characters.items() 
                if 1 in c)
    else:
        # delete all characters that are definitely not present at this node
        # and all character with state 0, if the domain does not exist in
        # the parent.
        node_item = list(node.characters.items()) # copy
        for a, c in node_item:
            # remove characters that are not present
            if len(c) == 1 and c == set([0]):
                del node.characters[a]
            elif len(c) == 2 and c == set([0,1]):
                del node.characters[a]
        # count numbers are now irrelevant
        node.characters= set(node.characters)

    for child in node.get_children():
        root_to_leaves(child)

def ancestral_states(node):
    """Calculate which characters have been gained or lost on all nodes in 
    the tree.
    """
    d = set(node.characters)
    if node.is_root():
        p = set()
    else:
        p = set(node.up.characters)
    gained = d.difference(p)
    lost = p.difference(d)
    node.add_features(gained_characters=gained, lost_characters=lost)

    if node.is_leaf():
        return

    for child in node.get_children():
        ancestral_states(child)

if __name__ == '__main__':
    import ete3
    tree = ete3.Tree('(((A,B),B2),(C2,(C,D)));')
    (tree&'A').add_features(characters = {'A':1})
    (tree&'B').add_features(characters = {}) 
    (tree&'B2').add_features(characters = {}) 
    (tree&'C').add_features(characters = {})
    (tree&'C2').add_features(characters = {"A":1})
    (tree&'D').add_features(characters = {"A":1})
    fitch(tree)

    print(tree)
    print("\ntraverse root to leaves and right to left\n")
    print("name\tdomains\tgained domains\tlostdomains")
    for node in tree.traverse( ) :
        name = node.name
        if node.is_root() : 
            name = "root"
        print("{}\t{}\t{}\t{}".format(name, "/".join(node.characters), "/".join(node.gained_characters), "/".join(node.lost_characters)))
    
