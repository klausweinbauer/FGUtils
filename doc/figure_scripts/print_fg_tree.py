from fgutils.fgconfig import print_tree, FGConfigProvider

provider = FGConfigProvider()
tree = provider.get_tree()
s = print_tree(tree)
