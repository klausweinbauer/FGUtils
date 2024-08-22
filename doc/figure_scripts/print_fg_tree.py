from fgutils.fgconfig import print_tree, FGConfigProvider, FGTreeNode

provider = FGConfigProvider()
tree = provider.get_tree()
s = print_tree(tree)
