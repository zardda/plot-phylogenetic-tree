from matplotlib.gridspec import GridSpec
from parse_newick import *


tree_str = "((((A:0.48420347,B:0.57779907):0.20194337,C:0.58826710):0.25546475,(D:0.92001531,E:0.90135602):0.00743693,F:0.91475001):0.06125499,(G:0.70355473,(H:0.46362493,I:0.43955339):0.18145267):0.43415809)root;"
root = parse_newick(tree_str)

fig = plt.figure()

gs = GridSpec(2, 2, figure=fig)
ax = fig.add_subplot(gs[0,1],projection='polar')
values = np.random.randn(2, 9)
add_cicleheatmap(root,values,ax)
plot_polar_tree(root,ax=ax,color={"node_1":"red","node_5":"blue","node_2":"green"})

ax1 = fig.add_subplot(gs[0:,0])
plot_tree(root,ax1,color={"node_1":"red","node_5":"blue","node_2":"green"},align_name=True)

ax2 = fig.add_subplot(gs[1,1],projection='polar')
plot_polar_tree(root,ax2,color={"node_1":"red","node_5":"blue","node_2":"green"},align_name=True,add_label=True)

plt.show()
# tree_str = "((((FgRID1:0.48420347,NcRID1:0.57779907):0.20194337,Masc1:0.58826710):0.25546475,(DNMT2:0.92001531,Masc2:0.90135602):0.00743693,DNMT1:0.91475001):0.06125499,(Masc3:0.70355473,(NcDIM2:0.46362493,FgDIM2:0.43955339):0.18145267):0.43415809)root;"
# root = parse_newick(tree_str)
# #plot_tree(root,color={"node_1":"red","node_5":"blue","node_2":"green"},align_name=False)
# # print(get_parent_node(root))
# fig,ax = plt.subplots(subplot_kw={'projection': 'polar'})
# plot_polar_tree(root,ax,color={"node_1":"red","node_5":"blue","node_2":"green"})
# plt.show()