import re
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d
import numpy as np
import math
import matplotlib
import matplotlib.cm as cm
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# https://matplotlib.org/stable/gallery/axes_grid1/simple_colorbar.html#sphx-glr-gallery-axes-grid1-simple-colorbar-py
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

class Node:
    def __init__(self, name, dist=None):
        self.name = name
        self.dist = dist
        self.children = []

def parse_newick(nw_str):
    tokens = re.findall(r"""[\w.:-]+|[,();]""", nw_str)
    stack = [[]]
    root = None
    children = []
    i = 0
    for t in tokens:
        if t == '(':
            stack.append([])
        elif t == ',':
            pass
        elif t == ')':
            children = stack.pop()
        else:
            if ':' in t:
                name, dist = t.split(':')
                if name == "":
                    name = f'node_{i}'
                    i += 1
                node = Node(name, float(dist))
            else:
                node = Node(t,0)
            if children:
                node.children = children
                children = []
            stack[-1].append(node)
        # print(stack)    
    root = stack[-1][0]
    return root

def get_all_node(node):
    mylist = []
    for child in node.children:
        mylist.extend(get_all_node(child))
    mylist.append(node)
    return mylist

def get_ypos(node):
    ymax = get_leaf_num(node)
    outsild_name = get_outsild_name(node)
    mydict = {}
    for name,y in zip(outsild_name,range(ymax)):
        mydict[name] = y
    for each_node in get_all_node(node):
        if each_node.name in outsild_name:
            pass
        else:
            cur_node_ypos = [mydict.get(name) for name in get_name_leaf(each_node)]
            cur_ymax = max(cur_node_ypos)
            cur_ymin = min(cur_node_ypos)
            mydict[each_node.name] = 0.5 * (cur_ymax + cur_ymin)
    return mydict

def get_parent_node(node,parent_name=None):
    if parent_name is None:
        parent_name = {"root":[node.dist]}
    for child in node.children:
        parent_name[child.name] = [dist for dist in parent_name[node.name]]
        parent_name[child.name].append(child.dist)
        if child.children:
            parent_name.update(get_parent_node(child,parent_name = parent_name))
    return parent_name

def get_xpos(node):
    node_relate = get_parent_node(node)
    x_posdict = {}
    for k in node_relate.keys():
        x_pos = 0
        for v in node_relate[k]:
            x_pos += v
        x_posdict[k] = x_pos
    return x_posdict

def get_xmax(node):
    x_pos = [v for v in get_xpos(node).values()]
    return max(x_pos)

def get_num_of_subnodes(node):
    return len(node.children)

def get_leaf_num(node):
    # return the number of total outside leaf 
    c = 0
    for child in node.children:
        if child.children:
            c += get_leaf_num(child)
        else:
            c += 1
    return c

def get_name_leaf(node):
    leaf_name = []
    for child in node.children:
        leaf_name.append(child.name)
    return leaf_name

def get_outsild_name(node):
    name = []
    for child in node.children:
        if child.children:
            name.extend(get_outsild_name(child))
        else:
            name.append(child.name)
    return name

def get_xypos(node,polar=False):
    # 返回三个值，x坐标，y坐标，到父节点的距离
    tol_name = [node.name for node in get_all_node(node)]
    dist_dict = {}
    for node in get_all_node(node):
        dist_dict[node.name] = node.dist
    x_pos = get_xpos(node)
    if polar:
        y_pos = polar_ypos(node)
    else:
        y_pos = get_ypos(node)
    xy_pos = {}
    for name in tol_name:
        xy_pos[name] = [x_pos[name],y_pos[name],dist_dict[name]]
    return xy_pos

def plot_hlines(node,color={}):
    mylines = [] 
    xy_pos = get_xypos(node)
    if color:
        color_node = color_nodes(node,color)
    else:
        color_node = {}
    for k in xy_pos.keys():
        x = xy_pos[k][0]
        y = xy_pos[k][1]
        d = xy_pos[k][2]
        x2 = x - d
        if k in color_node.keys():
            mylines.append(Line2D([x2,x],[y,y],color=color_node[k]))
        else:
            mylines.append(Line2D([x2,x],[y,y],color="black"))
    return mylines

def add_name(node,ax,align_name=True,polar=False):
    outside_name = get_outsild_name(node)
    if polar:
        xy_pos = get_xypos(node,polar=polar)
        xmax = get_xmax(node) * 1.02
        for n in outside_name:
            y = xy_pos[n][1]
            y_rad = y/360 * 2 * math.pi
            if align_name:
                ax.text(x=y_rad,y=xmax,s=n,rotation=y,
                        verticalalignment='center',rotation_mode="anchor")
            else:
                x = xy_pos[n][0]
                print(x,y_rad)
                ax.text(x=y_rad,y=x,s=n,rotation=y,
                        verticalalignment='center',rotation_mode="anchor")
    else: 
        xy_pos = get_xypos(node) 
        xmax = get_xmax(node)
        for n in outside_name:
            y = xy_pos[n][1]
            if align_name:
                ax.text(xmax,y,n,verticalalignment='center')
            else:
                x = xy_pos[n][0]
                ax.text(x,y,n,verticalalignment='center')  

def get_xmax(node):
    return max([v for v in get_xpos(node).values()])

def get_max_min_ynode(node,polar=False):
    outsild_name = get_outsild_name(node)
    if polar:
        y_pos = polar_ypos(node)
    else:
        y_pos = get_ypos(node)
    max_min_ypos = {}
    for each_node in get_all_node(node):
        if each_node.name in outsild_name:
            pass
        else:
            cur_node_ypos = [y_pos.get(name) for name in get_name_leaf(each_node)]
            cur_ymax = max(cur_node_ypos)
            cur_ymin = min(cur_node_ypos)
            max_min_ypos[each_node.name] = [cur_ymax,cur_ymin]
    return max_min_ypos

def add_dash(node,polar=False):
    outside_name = get_outsild_name(node)
    xy_pos = get_xypos(node,polar=polar)
    dash_lines = []
    xmax = get_xmax(node)
    for n in outside_name:
        x = xy_pos[n][0]
        y = xy_pos[n][1]
        if polar:
            y_rad = y/360 * 2 * math.pi
            dash_lines.append(Line2D([y_rad,y_rad],[x,xmax+0.001],linestyle="--",color="gray"))
        else:
            dash_lines.append(Line2D([x,xmax],[y,y],linestyle="--",color="gray"))
    return dash_lines

def plot_vlines(node,color={}):
    my_vlines = []
    tar_ypos = get_max_min_ynode(node)
    x_pos = get_xpos(node)
    if color:
        color_node = color_nodes(node,color)
    else:
        color_node = {}
    for k in tar_ypos.keys():
        x = x_pos[k]
        y = tar_ypos[k]
        if k in color_node.keys():
            my_vlines.append(Line2D([x,x],y,color=color_node[k]))
        else:
            my_vlines.append(Line2D([x,x],y,color="black"))
    return my_vlines

def plot_tree(node,ax,color={},align_name=True):
    for l in plot_hlines(node,color=color):
        ax.add_line(l)
    for l in plot_vlines(node,color=color):
        ax.add_line(l)
    if align_name:
        for l in add_dash(node,polar=False):
            ax.add_line(l)
    add_name(node,ax=ax,align_name=align_name)
    xmax = get_xmax(node)
    ax.set_xlim(0,xmax*1.2)
    ymax = get_leaf_num(node)
    ax.set_ylim(0,ymax)
    # 设置ax的详细信息
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    ax.invert_yaxis()

# 绘制环形的图

def polar_ypos(node,start=0,end=350):
    normal_ypos = get_ypos(node)
    ymax = get_leaf_num(node)
    polar_ypos = {}
    for k in normal_ypos.keys():
        polar_ypos[k] = normal_ypos[k] / ymax * (end - start) + start
    return polar_ypos

def add_polar_avr(r,ax,start,end,color):
    for curve in [[[start, end], [r,r]]]:
        curve[0] = np.deg2rad(curve[0])
        x = np.linspace(curve[0][0], curve[0][1], 500)
        y = interp1d(curve[0], curve[1])(x)
        ax.plot(x,y,color=color)

def add_polar_curve(node,ax,color={}):
    tar_ypos = get_max_min_ynode(node,polar=True)
    x_pos = get_xpos(node)
    if color:
        color_node = color_nodes(node,color)
    else:
        color_node = {}
    for k in tar_ypos.keys():
        x = x_pos[k]
        y = tar_ypos[k]
        if k in color_node.keys():
            add_polar_avr(x,ax,y[1],y[0],color=color_node[k])
        else:
            add_polar_avr(x,ax,y[1],y[0],color="black")

def add_polar_lines(node,ax,color={}):
    xy_pos = get_xypos(node,polar=True)
    if color:
        color_node = color_nodes(node,color)
    else:
        color_node = {}
    for k in xy_pos.keys():
        x = xy_pos[k][0]
        y = np.deg2rad(xy_pos[k][1])
        d = xy_pos[k][2]
        x2 = x - d
        if k in color_node.keys():
            ax.plot([y,y],[x2,x],color=color_node[k])
        else:
            ax.plot([y,y],[x2,x],color="black")
    
def plot_polar_tree(node,ax,color=[],align_name=False,polar=True,add_label=False):
    ax = ax 
    add_polar_curve(node,ax,color=color)
    add_polar_lines(node,ax,color=color)
    if add_label:
        add_name(node,ax,align_name=align_name,polar=polar)
    if align_name:
        for l in add_dash(node,polar=polar):
            ax.add_line(l)    
    ax.yaxis.grid(False)
    ax.spines['polar'].set_visible(False)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])

def color_nodes(node,color:dict):
    all_node = get_all_node(node)
    color_dict = {}
    for node in all_node:
        tar_dict = {}
        if node.name in color.keys():
            tar_node = [n.name for n in get_all_node(node)]
            tar_color = [color[node.name]] * len(tar_node)
            tar_dict = {key:value for key,value in zip(tar_node,tar_color)}
        if tar_dict:
            color_dict.update(tar_dict)
    return color_dict

# x为极坐标中的r
# y为极坐标中的角度

def add_polar_avr_grid(values,width,ax,start_rad,end_rad):
    y_bin = get_ybin(values,start_rad,end_rad)
    ypos = get_ypos_polar(values,start_rad,end_rad)[0] + 0.5 * y_bin
    ymin = min(ypos) - y_bin
    ymax = max(ypos)
    x_bin = get_xbin(values,width)
    x_pos = list(get_xpos_polar(values,width)[0] - 0.5 * x_bin)
    x_pos.append(x_pos[-1] + x_bin)
    for r in x_pos:
        for curve in [[[ymin,ymax], [r,r]]]:
            x = np.linspace(curve[0][0], curve[0][1], 500)
            y = interp1d(curve[0], curve[1])(x)
            ax.plot(x,y,color="gray",lw=.5)

def add_vline_polar(values,start_rad,end_rad,width,ax):
    y_bin = get_ybin(values,start_rad,end_rad)
    ypos = list(get_ypos_polar(values,start_rad,end_rad)[0] + 0.5 * y_bin)
    ypos.append(min(ypos)-y_bin)
    r_min = get_rmin_max(values,width)[1]
    r_max = get_rmin_max(values,width)[2]
    for y in ypos:
        ax.plot([y,y],[r_min,r_max],color="gray",lw=.5)

def get_color(values):
    # 根据数值的大小将数值映射为对应的颜色
    data = values.reshape(-1)
    cmap = matplotlib.colormaps.get_cmap('bwr')
    norm = plt.Normalize(data.min(), data.max())
    color = cmap(norm(data))
    return color

def get_rmin_max(values,width):
    r = values.shape[0]
    r_min = width * 1.01
    r_max = r_min + 0.1 * r_min * r
    return r,r_min,r_max

def get_xpos_polar(values,width):
    r  = get_rmin_max(values,width)[0]
    r_min = get_rmin_max(values,width)[1]
    r_max = get_rmin_max(values,width)[2]
    x_bin = (r_max - r_min) / r
    x_pos = np.linspace(r_min, r_max-x_bin,r) + 0.5 * x_bin
    return x_pos,x_bin

def get_xbin(values,width):
    return get_xpos_polar(values,width)[1]

def get_ypos_polar(values,start_rad,end_rad):
    degree = values.shape[1]
    y_bin = (end_rad - start_rad) / degree
    y_pos = np.linspace(start_rad,end_rad-y_bin,degree)
    return y_pos,y_bin

def get_ybin(values,start_rad,end_rad):
    return get_ypos_polar(values,start_rad,end_rad)[1]

def get_value(values,start_rad,end_rad,width):
    x_pos = get_xpos_polar(values,width)[0]
    y_pos = get_ypos_polar(values,start_rad,end_rad)[0]
    xy_pos = [[x,y] for x in x_pos for y in y_pos]
    colors = get_color(values)
    for i,l in enumerate(xy_pos):
        l.append(tuple(colors[i]))
    return xy_pos

def plot_each_area(ax,pos:list,xbin,ybin):
    # 根据中心坐标得到每一小块四个端点的坐标，+- 1/2的间隔
    x_pos = pos[1]
    y_pos = pos[0]
    color = pos[2]
    x = np.linspace(x_pos-0.5*xbin, x_pos+0.5*xbin, 50)
    y = interp1d([x_pos-0.5*xbin, x_pos+0.5*xbin], [y_pos,y_pos])(x)
    y1 = y - 0.5*ybin
    y2 = y + 0.5*ybin
    ax.fill_between(x,y1,y2,where=y1<y2,facecolor=color)

def add_color_bar(values,ax):
    cax = inset_axes(ax, width="4%", height="20%", loc='upper right')
    cmap = matplotlib.colormaps.get_cmap('bwr')
    norm = plt.Normalize(values.min(), values.max())
    plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,ticks=np.linspace(-1, 1, 3))

def add_cicleheatmap(node,values,ax,start_rad=0,end_rad= 350/360*2*math.pi,):
    width = get_xmax(node)
    xy_pos = get_value(values,start_rad,end_rad,width)
    y_bin = get_ybin(values,start_rad,end_rad)
    x_bin = get_xbin(values,width)
    for pos in xy_pos:
        # print(pos,y_bin,x_bin)
        plot_each_area(ax,pos,y_bin,x_bin)
    add_polar_avr_grid(values,width,ax,start_rad,end_rad)
    add_vline_polar(values,start_rad,end_rad,width,ax)
    add_color_bar(values,ax)



# 柱状图的绘制类似
# 对于平面的矩阵图可以使用gridspec https://matplotlib.org/stable/gallery/subplots_axes_and_figures/gridspec_multicolumn.html
# 箱线图也类似

# 计算层级就是统计“）”的个数并加一
# def print_names_2(node):
#     for child in node.children:
#         print_names(child)
#     print(node.name)
# # 在最后手动加入root的标识
